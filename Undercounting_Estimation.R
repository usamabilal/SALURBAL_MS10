rm(list=ls())
source("MS10_SALURBAL_Helper.R")
library(readxl)
library(tidyverse)
library(data.table)
library(DDM)
library(lspline)
countries<-c("AR", "BR", "CL", "CR", "CO", "MX", "PE", "PA", "SV")
# get population and mortality data
country<-"CL"
age_limit<-75
mortality_dir<-"../SALURBAL_DATA/Mortality Data/DIN/"
pop_dir<-"../SALURBAL_DATA/Population Data/L1AD/Age Category/"
all_data<-map_dfr(countries, function(country){
  print(country)
  files<-list.files(mortality_dir)
  file<-files[grepl(country, files)]
  mortality<-as.data.frame(fread(paste0(mortality_dir, file)))
  mortality<-mortality %>% 
    filter(substr(SALID1, 4, 6)!="888") %>% 
    group_by(SALID1, DINAGE5C, DINMALE, YEAR) %>% 
    summarise(deaths=n()) %>% 
    mutate(iso2=country)
  files<-list.files(pop_dir)
  file<-files[grepl(country, files)]
  pop<-as.data.frame(fread(paste0(pop_dir, file)))
  pop<-pop %>% 
    filter(substr(SALID1, 4, 6)!="888") %>% 
    group_by(SALID1, PRJAGE5C, PRJMALE, YEAR) %>% 
    summarise(pop=sum(PRJL1ADPOP))
  pop<-pop %>% mutate(age=PRJAGE5C,
                      age=floor(age/5)*5,
                      age=ifelse(age>=age_limit, age_limit, age),
                      sex=ifelse(PRJMALE==1, "M", "F")) %>% 
    group_by(SALID1, YEAR, sex, age) %>% 
    summarise(pop=sum(pop))
  mortality<-mortality %>% 
    mutate(age=DINAGE5C,
           age=floor(age/5)*5,
           age=ifelse(age>=age_limit, age_limit, age),
           sex=ifelse(DINMALE==1, "M", "F")) %>% 
    group_by(SALID1, YEAR, sex, age) %>% 
    summarise(deaths=sum(deaths))
  both<-full_join(pop, mortality) 
  if (country=="SV"){
    both<-both %>% filter(YEAR%in%c(2010:2014))
  } else {
    both<-both %>% filter(YEAR%in%c(2012:2016))
  }
  both<-both %>% 
    mutate(deaths=replace_na(deaths, 0),
           iso2=country)
  both
})
l1s<-data.frame(iso2=c("AR", "BR", "CL", "CR", "CO", "MX", "PE", "PA", "SV"),
                country_name=c("Argentina", "Brazil", "Chile", "Costa Rica", "Colombia",
                               "Mexico", "Peru", "Panama", "El Salvador"), 
                stringsAsFactors = F)
# getting ex on the open ended group for each country from UNDP WPP 2019
dta1<-(read_excel("../../SALURBAL Shared/SALURBAL_MS10/Other_Data/UNDP Life Tables/WPP2019_MORT_F17_2_ABRIDGED_LIFE_TABLE_MALE.xlsx", skip=16)) %>%
  mutate(sex="M")
dta2<-(read_excel("../../SALURBAL Shared/SALURBAL_MS10/Other_Data/UNDP Life Tables/WPP2019_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE.xlsx", skip=16)) %>%
  mutate(sex="F")
UNDP<-map_dfr(list(dta1, dta2),function(temp){
  temp %>% 
    rename(country_name=`Region, subregion, country or area *`,
           age=`Age (x)`,
           le_country=`Expectation of life e(x)`) %>% 
    mutate(le_country=as.numeric(le_country),
           y1=as.numeric(substr(Period, 1, 4)),
           y2=as.numeric(substr(Period, 6, 9))) %>% 
    rowwise() %>% 
    mutate(year=mean(c(y1, y2))) %>% 
    select(country_name, year, age, le_country) %>% 
    right_join(l1s) %>% 
    group_by(iso2) %>% 
    # getting yearly LE using a linear interpolation
    group_modify(~{
      .x<-.x %>% filter(age==75)
      fit<-lm(le_country ~ lspline(year,knots = seq(1952.5, 2017.5, by=5), marginal=T),data = .x )
      emax_list<-expand.grid(year=1950:2020)
      emax_list$emax<-predict(fit, newdata=emax_list)
      emax_list
    }) %>% mutate(sex=unique(temp$sex))
})
# get pop in 2010/2012 and 2014/2016
pop<-all_data %>% select(SALID1, YEAR, sex, age, pop, iso2) %>% 
  filter((iso2=="SV"&YEAR%in%c(2010, 2014))|
           (iso2!="SV"&YEAR%in%c(2012, 2016))) %>% 
  ungroup() %>% 
  mutate(YEAR=ifelse(YEAR%in%c(2010, 2012), "y1", "y2")) %>% 
  spread(YEAR, pop) %>% 
  rename(pop1=y1, pop2=y2)
# get average yearly deaths from 2010-2014 or 2012-2016
mortality<-all_data %>% select(SALID1, YEAR, sex, age, deaths, iso2) %>% 
  filter((iso2=="SV"&YEAR%in%c(2010:2014))|
           (iso2!="SV"&YEAR%in%c(2012:2016))) %>% 
  group_by(sex, age, SALID1, iso2) %>% 
  summarise(deaths=mean(deaths))
# combine. Set date to mid-point through the year
both<-full_join(mortality, pop) %>% 
  mutate(date1=as.Date(ifelse(iso2=="SV", "2010-06-30", "2012-06-30")),
         date2=as.Date(ifelse(iso2=="SV", "2014-06-30", "2016-06-30"))) %>% 
  # midpoint year for eMax
  mutate(year=ifelse(iso2=="SV", 2012, 2014)) %>% 
  left_join(UNDP) %>% 
  rename(cod=SALID1) %>% 
  select(cod, sex, age, deaths, pop1, pop2, date1, date2, emax, iso2) %>% 
  arrange(cod, age)
# specific age bands for hill (30-65) and murray (50-70)
ages_hill<-seq(30, 65, by=5)
ages_murray<-seq(50, 70, by=5)
#.x<-both %>% filter(iso2=="AR", sex=="F") %>% ungroup() %>% select(-sex, -iso2)
# calculate correction factors for each country
correction<-both %>% group_by(iso2, sex) %>% 
  group_modify(~{
    emax<-.x %>% pull(emax) %>% unique
    .x<-.x %>% select(-emax)
    ddm_murray<-ddm(.x, eOpen=emax, exact.ages = ages_murray, deaths.summed = F)
    ddm_hill<-ddm(.x, eOpen=emax, exact.ages = ages_hill)
    ddm_auto<-ddm(.x, eOpen=emax)
    bind_rows(ddm_murray %>% mutate(ages="murray"), 
              ddm_hill %>% mutate(ages="hill"),
              ddm_auto %>% mutate(ages="auto")) %>% 
      select(cod, ggb, seg, ggbseg, ages)
  }) %>% 
  rowwise() %>% 
  rename(SALID1=cod) %>% 
  gather(method, ucnt, -iso2, -ages, -SALID1, -sex) %>% 
  # cap at 1
  mutate(ucnt=ifelse(ucnt>1, 1, ucnt))
# final datasets: iso2/sex/city_id/3 methods/3 age bands/coverage
save(correction, file="analytic files/undercounting_correction_bysex.rdata")

