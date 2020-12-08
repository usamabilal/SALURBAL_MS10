# load harmonized SALURBAL data and save analytic file
rm(list=ls())
library(data.table)
library(tidyverse)
library(DemoTools)
library(furrr)
source("MS10_SALURBAL_Helper.R")
# set age limit (open ended age group)
age_limit<-85
# get population
country_list<-c("AR", "BR", "CL", "CO", "CR", "MX", "PA", "PE", "SV")
pop_dir<-"../SALURBAL_DATA/Population Data/L1AD/Age Category/"
pop_files<-list.files(pop_dir, pattern="L1")
plan(multiprocess, .cleanup=T)
# country<-"SV"
population<-future_map_dfr(country_list, function(country){
  print(country)
  file<-pop_files[grepl(country, pop_files)]
  population<-fread(paste0(pop_dir, file),stringsAsFactors = F)
  # filter to 2010:2016
  population<-population %>% filter(YEAR%in%(2010:2016))
  agemax<-max(population$PRJAGE)
  population$year<-population$YEAR
  population$age<-population$PRJAGE
  population$male<-population$PRJMALE
  population<-population %>% group_by(SALID1, year, age, male) %>% 
    summarise(pop=sum(PRJL1ADPOP)) %>% 
    ungroup()
  # for population, if agemax<85, extend age to 85 (age_limit)
  if (agemax<age_limit) {
    population<-population %>% group_by(SALID1, year, male) %>% 
      group_modify(~{
        #.x<-population %>% filter(SALID1==population$SALID1[[1]], year==2012, male==1)
        .x<-.x %>% arrange(age)
        t<-graduate(Value = .x$pop, Age=.x$age, method="pclm",
                    OAnew =  115, keep0=T, constrain = T)
        data.frame(pop=t, age=as.numeric(names(t))) %>% 
          mutate(age=case_when(
            age==0 ~ 0,
            age%in%(1:4) ~ 1,
            age>=age_limit ~ age_limit,
            T ~ floor(age/5)*5
          )) %>% 
          group_by(age) %>% 
          summarise(pop=sum(pop))
      })
  } else {
    population<-population %>% 
      mutate(age=case_when(
        age==0 ~ 0,
        age%in%(1:4) ~ 1,
        age>=age_limit ~ age_limit,
        T ~ floor(age/5)*5
      )) %>% group_by(SALID1, year, age, male) %>% 
      summarise(pop=sum(pop))
  }
  # fill in any potential gaps
  template<-expand.grid(SALID1=unique(population$SALID1),
                        year=unique(population$year),
                        age=unique(population$age),
                        male=unique(population$male))
  population<-full_join(population, template) %>% 
    mutate(pop=replace_na(pop, 0),
           iso2=country)
  population
})

# get all-cause mortality
mortality_dir<-"../SALURBAL_DATA/Mortality Data/DIN/"
mortality_files<-list.files(mortality_dir, pattern="csv")
all_cause_mortality<-map_dfr(country_list, function(country){
  print(country)
  file<-mortality_files[grepl(country, mortality_files)]
  mortality<-fread(paste0(mortality_dir, file),stringsAsFactors = F)
  # filter to 2010:2016
  mortality<-mortality %>% filter(YEAR%in%(2010:2016)) %>% 
    mutate(year=YEAR,
           n=1,
           male=DINMALE,
           age=case_when(
             DINAGE5C==0 ~ 0,
             DINAGE5C%in%(1:4)~ 1,
             DINAGE5C>=age_limit ~ age_limit,
             T ~ floor(DINAGE5C/5)*5
           ))%>% 
    group_by(SALID1, year, age, male) %>% 
    summarise(deaths=sum(n))
  # fill in any potential gaps
  template<-expand.grid(SALID1=unique(mortality$SALID1),
                        year=unique(mortality$year),
                        age=unique(mortality$age),
                        male=unique(mortality$male))
  mortality<-full_join(mortality, template) %>% 
    mutate(deaths=replace_na(deaths, 0),
           iso2=country)
})

# cause-specific mortality (redistributed)
mortality_dir<-"../SALURBAL_DATA/Mortality Data/DIN/"
mortality_files<-list.files(mortality_dir, pattern="csv")
cause_specific_red<-map_dfr(country_list, function(country){
  print(country)
  file<-mortality_files[grepl(country, mortality_files)]
  mortality<-fread(paste0(mortality_dir, file),stringsAsFactors = F)
  mortality<-mortality %>% filter(YEAR%in%(2010:2016)) %>% 
    mutate(year=YEAR,
           n=1,
           male=DINMALE,
           age=case_when(
             DINAGE5C==0 ~ 0,
             DINAGE5C%in%(1:4)~ 1,
             DINAGE5C>=age_limit ~ age_limit,
             T ~ floor(DINAGE5C/5)*5
           ))%>% 
    group_by(SALID1, DINCOD_FINAL2, year, age, male) %>% 
    summarise(deaths=sum(n)) %>% 
    rename(DTHCOD_FINAL2=DINCOD_FINAL2)
  # fill in any potential gaps
  template<-expand.grid(SALID1=unique(mortality$SALID1),
                        DTHCOD_FINAL2=unique(mortality$DTHCOD_FINAL2),
                        year=unique(mortality$year),
                        age=unique(mortality$age),
                        male=unique(mortality$male))
  mortality<-full_join(mortality, template) %>% 
    mutate(deaths=replace_na(deaths, 0),
           iso2=country)
  mortality
})
# cause specific, not redistributed
mortality_dir<-"../SALURBAL_DATA/Mortality Data/DIN/"
mortality_files<-list.files(mortality_dir, pattern="csv")
cause_specific_nored<-map_dfr(country_list, function(country){
  print(country)
  file<-mortality_files[grepl(country, mortality_files)]
  mortality<-fread(paste0(mortality_dir, file),stringsAsFactors = F)
  mortality<-mortality %>% filter(YEAR%in%(2010:2016)) %>% 
    mutate(year=YEAR,
           n=1,
           male=DINMALE,
           age=case_when(
             DINAGE5C==0 ~ 0,
             DINAGE5C%in%(1:4)~ 1,
             DINAGE5C>=age_limit ~ age_limit,
             T ~ floor(DINAGE5C/5)*5
           ))%>% 
    group_by(SALID1, DINCOD_GHE2, year, age, male) %>% 
    summarise(deaths=sum(n)) %>% 
    rename(DTHCOD_GHE2=DINCOD_GHE2)
  # fill in any potential gaps
  template<-expand.grid(SALID1=unique(mortality$SALID1),
                        DTHCOD_GHE2=unique(mortality$DTHCOD_GHE2),
                        year=unique(mortality$year),
                        age=unique(mortality$age),
                        male=unique(mortality$male))
  
  mortality<-full_join(mortality, template) %>% 
    mutate(deaths=replace_na(deaths, 0),
           iso2=country)
  mortality
})
# get total population (unmatched with mortality)
total_population<-map_dfr(country_list, function(country){
  print(country)
  file<-pop_files[grepl(country, pop_files)]
  population<-fread(paste0(pop_dir, file),stringsAsFactors = F)
  population<-population %>% 
    rename(year=YEAR) %>% 
    group_by(SALID1, year) %>% 
    summarise(total_pop=sum(PRJL1ADPOP)) %>% 
    mutate(iso2=country)
  population
})

# exclude non-SALURBAL areas
population<-population %>% filter(substr(SALID1, 4, 6)!="888")
total_population<-total_population %>% filter(substr(SALID1, 4, 6)!="888")
all_cause_mortality<-all_cause_mortality %>% filter(substr(SALID1, 4, 6)!="888")
cause_specific_red<-cause_specific_red %>% filter(substr(SALID1, 4, 6)!="888")
cause_specific_nored<-cause_specific_nored %>% filter(substr(SALID1, 4, 6)!="888")

save(population, total_population,
     all_cause_mortality,
     cause_specific_red,
     cause_specific_nored,
     file="analytic files/all_data_mortality_population_corrected_level1_age14.RData")

# estimation of undercounting factors
rm(list=ls())
source("MS10_SALURBAL_Helper.R")
library(readxl)
library(tidyverse)
library(data.table)
library(DDM)
library(lspline)
age_limit<-85
# getting life expectancy on the open ended age group at the country level from UNDP's WPP 2019
l1s<-data.frame(iso2=c("AR", "BR", "CL", "CR", "CO", "MX", "PE", "PA", "SV"),
                country_name=c("Argentina", "Brazil", "Chile", "Costa Rica", "Colombia",
                               "Mexico", "Peru", "Panama", "El Salvador"), 
                stringsAsFactors = F)
dta1<-(read_excel("Other_Data/UNDP Life Tables/WPP2019_MORT_F17_2_ABRIDGED_LIFE_TABLE_MALE.xlsx", skip=16)) %>%
  mutate(sex="M")
dta2<-(read_excel("Other_Data/UNDP Life Tables/WPP2019_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE.xlsx", skip=16)) %>%
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
      .x<-.x %>% filter(age==age_limit)
      fit<-lm(le_country ~ lspline(year,knots = seq(1952.5, 2017.5, by=5), marginal=T),data = .x )
      emax_list<-expand.grid(year=1950:2020)
      emax_list$emax<-predict(fit, newdata=emax_list)
      emax_list
    }) %>% mutate(sex=unique(temp$sex))
})
# get pop in 2010/2012 and 2014/2016
load("analytic files/all_data_mortality_population_corrected_level1_age14.RData")
pop<-population %>% select(SALID1, year, male, age, pop, iso2) %>% 
  mutate(age=floor(age/5)*5) %>% 
  group_by(iso2, SALID1, year, male, age) %>% 
  summarise(pop=sum(pop)) %>% 
  filter((iso2=="SV"&year%in%c(2010, 2014))|
           (iso2!="SV"&year%in%c(2012, 2016))) %>% 
  ungroup() %>% 
  mutate(year=ifelse(year%in%c(2010, 2012), "y1", "y2")) %>% 
  spread(year, pop) %>% 
  rename(pop1=y1, pop2=y2)
# get average yearly deaths from 2010-2014 or 2012-2016
mortality<-all_cause_mortality %>% ungroup() %>% select(SALID1, year, male, age, deaths, iso2) %>% 
  mutate(age=floor(age/5)*5) %>% 
  group_by(iso2, SALID1, year, male, age) %>% 
  summarise(deaths=sum(deaths)) %>% 
  filter((iso2=="SV"&year%in%c(2010:2014))|
           (iso2!="SV"&year%in%c(2012:2016))) %>% 
  group_by(male, age, SALID1, iso2) %>% 
  summarise(deaths=mean(deaths))
# combine. Set date to mid-point through the year
both<-full_join(mortality, pop) %>% 
  mutate(date1=as.Date(ifelse(iso2=="SV", "2010-06-30", "2012-06-30")),
         date2=as.Date(ifelse(iso2=="SV", "2014-06-30", "2016-06-30"))) %>% 
  # midpoint year for eMax
  mutate(year=ifelse(iso2=="SV", 2012, 2014),
         sex=ifelse(male==1, "M", "F")) %>% 
  left_join(UNDP) %>% 
  rename(cod=SALID1) %>% 
  ungroup() %>% 
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
    ddm_murray<-ddm(.x, eOpen=emax, exact.ages = ages_murray)
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


## calculate undercounting irrespective of sex
# getting ex on the open ended group for each country from UNDP WPP 2019
dta_both<-(read_excel("../../SALURBAL Shared/SALURBAL_MS10/Other_Data/UNDP Life Tables/WPP2019_MORT_F17_1_ABRIDGED_LIFE_TABLE_BOTH_SEXES.xlsx", skip=16))
UNDP<-dta_both %>% 
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
    .x<-.x %>% filter(age==age_limit)
    fit<-lm(le_country ~ lspline(year,knots = seq(1952.5, 2017.5, by=5), marginal=T),data = .x )
    emax_list<-expand.grid(year=1950:2020)
    emax_list$emax<-predict(fit, newdata=emax_list)
    emax_list
  })

# combine. Set date to mid-point through the year
both<-both %>% group_by(cod, age, date1, date2, iso2) %>% 
  summarise(deaths=sum(deaths),
            pop1=sum(pop1),
            pop2=sum(pop2)) %>% 
  mutate(year=ifelse(iso2=="SV", 2012, 2014)) %>% 
  left_join(UNDP) %>% 
  select(-year)

# specific age bands for hill (30-65) and murray (50-70)
ages_hill<-seq(30, 65, by=5)
ages_murray<-seq(50, 70, by=5)
#.x<-both %>% filter(iso2=="AR") %>% ungroup() %>% select(-sex, -iso2)
# calculate correction factors for each country
correction<-both %>% group_by(iso2) %>% 
  group_modify(~{
    emax<-.x %>% pull(emax) %>% unique
    .x<-.x %>% select(-emax)
    ddm_murray<-ddm(.x, eOpen=emax, exact.ages = ages_murray)
    ddm_hill<-ddm(.x, eOpen=emax, exact.ages = ages_hill)
    ddm_auto<-ddm(.x, eOpen=emax)
    bind_rows(ddm_murray %>% mutate(ages="murray"), 
              ddm_hill %>% mutate(ages="hill"),
              ddm_auto %>% mutate(ages="auto")) %>% 
      select(cod, ggb, seg, ggbseg, ages)
  }) %>% 
  rowwise() %>% 
  rename(SALID1=cod) %>% 
  gather(method, ucnt, -iso2, -ages, -SALID1) %>% 
  # cap at 1
  mutate(ucnt=ifelse(ucnt>1, 1, ucnt))

# final datasets: iso2/city_id/3 methods/3 age bands/coverage
save(correction, file="analytic files/undercounting_correction.rdata")


## save database of exposures
rm(list=ls())
source("MS10_SALURBAL_Helper.R")
library(tidyverse)
library(data.table)
library(readxl)
age_limit<-85
country_list<-c("argentina", "brazil", "chile", "mexico", 
                "peru", "colombia", "salvador", "guatemala", "costarica", "panama", "nicaragua")
country_list2<-c("Argentina", "Brazil", "Chile", "Mexico", 
                 "Peru", "Colombia", "El Salvador", "Guatemala", "Costa Rica", "Panama", "Nicaragua")
short_labels<-data.frame(country=country_list,
                         country2=country_list2,
                         short=c("AR", "BR", "CL", "MX", "PE", "CO", "SV", "GT", "CR", "PA", "NI"), 
                         salid0=c(101, 102, 103, 204, 105, 104, 202, 203, 201, 206, 205),
                         stringsAsFactors = F)
load("analytic files/l1s.rdata")
l1s$salid0<-as.numeric(substr(l1s$SALID1, 1, 3))
l1s<-merge(l1s, short_labels, by="salid0", all=T)
l1s$city<-as.numeric(!is.na(l1s$SALID1))

# all countries but el salvador start in 2012
years<-data.frame(iso2=c("AR", "BR", "CL", "CO", "CR", "MX", "PA", "PE", "SV"))
years<-years %>% mutate(y1=ifelse(iso2=="SV", 2010, 2012))
years<-years %>% 
  mutate(
    # study period (y1-y1+4)
    y2=y1+4,
    # same period, but right before
    y1b=y1-5,
    y2b=y1-1)
years
# load population data
load("analytic files/all_data_mortality_population_corrected_level1_age14.RData")
pop_l1<-total_population %>% right_join(years) %>% 
  filter(year==y1 | year==y2) %>% 
  mutate(year=ifelse(year==y1,"y1", "y2")) %>% 
  spread(year, total_pop) %>% 
  mutate(growth=y2-y1,
         growth_pct=growth/y1,
         pop_baseline=y1) %>% 
  select(SALID1, growth, growth_pct, pop_baseline)
pop_growth_prev<-total_population %>% right_join(years) %>% 
  filter(year==y1b | year==y2b) %>% 
  mutate(year=ifelse(year==y1b,"y1", "y2")) %>% 
  select(SALID1, year, total_pop) %>% 
  spread(year, total_pop) %>% 
  mutate(growth_pct_prev=(y2-y1)/y1) %>% 
  select(SALID1, growth_pct_prev) %>% 
  # remove those without preivous years data
  filter(!is.na(growth_pct_prev))
pop_l1<-full_join(pop_l1, pop_growth_prev)
pop_l1<-left_join(pop_l1, l1s %>% select(SALID1, iso2))
summary(pop_l1)

# BEC variables
bec<-read.csv("../SALURBAL_DATA/BEC Data/BEC_L1AD_20200506.csv", stringsAsFactors = F)
bec$total_transport<-rowSums(bec[,colnames(bec)[grepl("BECPR", colnames(bec))]])
bec$total_transport_yn<-as.numeric(bec$total_transport>=1)
bec<-bec[,c("SALID1", colnames(bec)[grepl("BEC", colnames(bec))],"total_transport", "total_transport_yn")]
becux<-read.csv("../SALURBAL_DATA/BEC Data/BEC_L1UX_20200506.csv", stringsAsFactors = F)
becux<-becux[,c("SALID1", colnames(becux)[grepl("BEC", colnames(becux))])]
becux<-becux[,!grepl("GASPRICE", colnames(becux))]
bec<-full_join(bec, becux)
head(bec)


## SEC variables
sec<-fread("../SALURBAL_DATA/SEC Data/SEC_Census_L1AD_07162020.csv") %>% 
  filter(YEAR!=2017) %>% select(-ISO2)
# calculate SEI
sec<-sec %>% ungroup() %>% 
  mutate(educationsd=scale(CNSMINPR_L1AD, scale=T, center=T),
         watersd=scale(CNSWATINL1AD, scale=T, center=T),
         sewagesd=scale(CNSSEWNETL1AD, scale=T, center=T),
         overcrowdingsd=scale(CNSCROWD3RML1AD, scale=T, center=T)*(-1)) %>% 
  rowwise() %>% 
  mutate(sei=mean(c(educationsd, watersd, sewagesd, overcrowdingsd)))

all_exposure<-full_join(pop_l1, sec) %>% 
  full_join(bec) %>% 
  filter(!iso2%in%c("GT", "NI"),
         !is.na(iso2))
summary(all_exposure)
head(all_exposure)
table(all_exposure$iso2)
summary(all_exposure[!complete.cases(all_exposure),colSums(is.na(all_exposure))>0])
# some floor (brazil) crowd (peru) popconc (not used)
# some brt, tram, bikelane, etc.
# Peru gas prices
# and 1 missing for street design (Zarate Campana)
# plus 40 (AR, CR, PA, SV) missing growth previous years

# get age-adjusted mortality for the same period
load("analytic files/all_data_mortality_population_corrected_level1_age14.RData")
# get who std population from SEER
stds<-read.fwf(file="Other_Data/standard populations/stdpop.19ages.txt", widths=c(3, 3, 8),
               col.names = c("std", "age", "pop")) %>% 
  filter(std==10) %>% 
  mutate(age=case_when(
    age==0 ~ 0,
    age==1 ~ 1,
    T ~ (age-1)*5),
    std_pop=pop/sum(pop)) %>% 
  group_by(age) %>% 
  summarise(std_pop=sum(std_pop))
sum(stds$std_pop)
# get undercounting correction factors
load("analytic files/undercounting_correction_bysex.rdata")
correction<-correction %>% 
  mutate(male=ifelse(sex=="M", 1, 0)) %>% 
  group_by(SALID1, male) %>% 
  summarise(correction=hmean(ucnt)) 
# calculating aamr
numerator<-all_cause_mortality %>% inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(SALID1,male, age) %>% 
  summarise(deaths=sum(deaths)) %>% 
  left_join(correction) %>% 
  mutate(deaths=deaths/correction) %>% 
  group_by(SALID1, age) %>% 
  summarise(deaths=sum(deaths))
denominator<-population %>% inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(SALID1, age) %>% 
  summarise(pop=sum(pop))
# get age adjusted mortality
aamr<-full_join(numerator, denominator) %>% 
  full_join(stds) %>% 
  ungroup() %>% 
  mutate(rate=deaths/pop,
         adjusted=rate*std_pop) %>% 
  group_by(SALID1) %>% 
  summarise(adj_rate=sum(adjusted)*100000)
hist(aamr$adj_rate, breaks=50)
aamr %>% left_join(l1s) %>% arrange(desc(adj_rate))
aamr %>% left_join(l1s) %>% arrange((adj_rate))

#pop for adjustment covariates (>65, <15, 15_64)
popadj<-population %>% inner_join(years) %>% filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(age2=ifelse(age>=65, "a65", 
                     ifelse(age>=15, "a15", 
                            "a0"))) %>% 
  group_by(SALID1, age2) %>% 
  summarise(pop=sum(pop)) %>% 
  ungroup() %>% 
  spread(age2, pop) %>% 
  rowwise() %>% 
  mutate(total=sum(c(a0, a15, a65)),
         under15=a0/total,
         between15and64=a15/total,
         above65=a65/total) %>% 
  select(SALID1, under15, between15and64, above65)
all_exposure<-full_join(all_exposure, aamr) %>% 
  full_join(popadj) %>% select(-YEAR)
all_exposure<-all_exposure %>% filter(!iso2%in%c("NI", "GT"), !is.na(SALID1))
head(all_exposure)
summary(all_exposure)

save(all_exposure,years, file="analytic files/MS10_exposure_data.RData")

# Save SHP for the App
rm(list=ls())
library(sf)
library(rmapshaper)
shp = st_read('../SALURBAL_DATA/SHPs/SALURBAL_L1AD_11_30_18/SALURBAL_L1AD_11_30_18.shp')
shp <- st_transform(shp, "+init=epsg:4326") %>% ms_simplify()
save(shp, file="MS10/l1ad_shp.rdata")


# prepare outcome data for cause-specific part
rm(list=ls())
age_limit<-85
library(tidyverse)
load("analytic files/l1s.RData")
load("analytic files/all_data_mortality_population_corrected_level1_age14.RData")
load("analytic files/MS10_exposure_data.RData")
source("MS10_SALURBAL_Helper.R")
# get undercounting correction factors
load("analytic files/undercounting_correction_bysex.rdata")
correction<-correction %>% 
  mutate(male=ifelse(sex=="M", 1, 0)) %>% 
  group_by(SALID1, male) %>% 
  summarise(correction=hmean(ucnt)) 
# calculate two versions: 
# one for descriptives with upweighted mortality by undercounting factor
pop<-population %>% inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(SALID1, age) %>% 
  summarise(pop=sum(pop))
mortality<-cause_specific_red %>% inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(cause_cat=case_when(
    DTHCOD_FINAL2%in%c(20,380, 
                       420, 490, 530)~"cmnn",
    DTHCOD_FINAL2%in%c(600, 780)~"cancer",
    DTHCOD_FINAL2%in%c(1040,
                       790, 800, 810, 980, 1110, 1150, 1200, 1240, 1250, 1310, 1430, 1475)~"ncd",
    DTHCOD_FINAL2%in%c(1490)~"accident",
    DTHCOD_FINAL2%in%c(1560)~"violent"
  ),
  age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(SALID1, male, cause_cat, age) %>% 
  summarise(deaths=sum(deaths)) %>% 
  left_join(correction) %>% 
  mutate(deaths=deaths/correction) %>% 
  group_by(SALID1, cause_cat, age) %>% 
  summarise(deaths=sum(deaths)) %>% 
  full_join(pop)

mortality_cause<-mortality %>% 
  group_by(SALID1, cause_cat) %>% 
  summarise(deaths=sum(deaths),
            pop=sum(pop)) %>% 
  as_tibble() %>% 
  ungroup() %>% 
  spread(cause_cat, deaths) %>% 
  rowwise() %>% 
  mutate(total=sum(c(accident, cancer, cmnn, ncd, violent)),
         p_cmnn=cmnn/total,
         p_cancer=cancer/total,
         p_ncd=ncd/total,
         p_accident=accident/total,
         p_violent=violent/total)

# and one for models with downweighted population by undercounting factor
pop<-population %>% inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(SALID1, male, age) %>% 
  summarise(pop=sum(pop)) %>% 
  left_join(correction) %>% 
  mutate(pop=pop*correction) %>% 
  group_by(SALID1, age) %>% 
  summarise(pop=sum(pop))
mortality<-cause_specific_red %>% inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(cause_cat=case_when(
    DTHCOD_FINAL2%in%c(20,380, 
                       420, 490, 530)~"cmnn",
    DTHCOD_FINAL2%in%c(600, 780)~"cancer",
    DTHCOD_FINAL2%in%c(1040,
                       790, 800, 810, 980, 1110, 1150, 1200, 1240, 1250, 1310, 1430, 1475)~"ncd",
    DTHCOD_FINAL2%in%c(1490)~"accident",
    DTHCOD_FINAL2%in%c(1560)~"violent"
  ),
  age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(SALID1, cause_cat, age) %>% 
  summarise(deaths=sum(deaths)) %>% 
  full_join(pop)

mortality_cause_models<-mortality %>% 
  group_by(SALID1, cause_cat) %>% 
  summarise(deaths=sum(deaths),
            pop=sum(pop)) %>% 
  as_tibble() %>% 
  ungroup() %>% 
  spread(cause_cat, deaths) %>% 
  rowwise() %>% 
  mutate(total=sum(c(accident, cancer, cmnn, ncd, violent)),
         p_cmnn=cmnn/total,
         p_cancer=cancer/total,
         p_ncd=ncd/total,
         p_accident=accident/total,
         p_violent=violent/total)
cause_titles_coll<-c("cmnn", "cancer", "ncd", "accident", "violent")

# calculate AAPM
# get who std population from SEER
stds<-read.fwf(file="Other_Data/standard populations/stdpop.19ages.txt", widths=c(3, 3, 8),
               col.names = c("std", "age", "pop")) %>% 
  filter(std==10) %>% 
  mutate(age=case_when(
    age==0 ~ 0,
    age==1 ~ 1,
    T ~ (age-1)*5),
    std_pop=pop/sum(pop)) %>% 
  group_by(age) %>% 
  summarise(std_pop=sum(std_pop))
sum(stds$std_pop)
pop<-population %>% inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(SALID1, age) %>% 
  summarise(pop=sum(pop))
mortality_age<-cause_specific_red %>% inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(cause_cat=case_when(
    DTHCOD_FINAL2%in%c(20,380, 
                       420, 490, 530)~"cmnn",
    DTHCOD_FINAL2%in%c(600, 780)~"cancer",
    DTHCOD_FINAL2%in%c(1040,
                       790, 800, 810, 980, 1110, 1150, 1200, 1240, 1250, 1310, 1430, 1475)~"ncd",
    DTHCOD_FINAL2%in%c(1490)~"accident",
    DTHCOD_FINAL2%in%c(1560)~"violent"
  ),
  age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(SALID1, male, cause_cat, age) %>% 
  summarise(deaths=sum(deaths)) %>% 
  left_join(correction) %>% 
  mutate(deaths=deaths/correction) %>% 
  group_by(SALID1, cause_cat, age) %>% 
  summarise(deaths=sum(deaths)) %>% 
  left_join(pop) %>% 
  left_join(stds) %>% 
  mutate(rate=deaths/pop,
         adjrate=rate*std_pop) %>% 
  group_by(SALID1, cause_cat) %>% 
  summarise(adjrate=sum(adjrate))
mortality_age<-mortality_age %>% 
  full_join(mortality_age %>% 
              group_by(SALID1) %>% 
              summarise(adjrate_total=sum(adjrate)))
aapm<-mortality_age %>% 
  mutate(aapm=adjrate/adjrate_total) %>% 
  select(SALID1, cause_cat, aapm) %>% 
  spread(cause_cat, aapm) %>% 
  rename_at(-1, ~paste0("p_", ., "_ageadj"))
# calculate aapm for the entire country
pop<-population %>% inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(iso2, age) %>% 
  summarise(pop=sum(pop))
mortality_age<-cause_specific_red %>% inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(cause_cat=case_when(
    DTHCOD_FINAL2%in%c(20,380, 
                       420, 490, 530)~"cmnn",
    DTHCOD_FINAL2%in%c(600, 780)~"cancer",
    DTHCOD_FINAL2%in%c(1040,
                       790, 800, 810, 980, 1110, 1150, 1200, 1240, 1250, 1310, 1430, 1475)~"ncd",
    DTHCOD_FINAL2%in%c(1490)~"accident",
    DTHCOD_FINAL2%in%c(1560)~"violent"
  ),
  age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(iso2, SALID1, male, cause_cat, age) %>% 
  summarise(deaths=sum(deaths)) %>% 
  left_join(correction) %>% 
  mutate(deaths=deaths/correction) %>% 
  group_by(iso2, cause_cat, age) %>% 
  summarise(deaths=sum(deaths)) %>% 
  left_join(pop) %>% 
  left_join(stds) %>% 
  mutate(rate=deaths/pop,
         adjrate=rate*std_pop) %>% 
  group_by(iso2, cause_cat) %>% 
  summarise(adjrate=sum(adjrate))
mortality_age<-mortality_age %>% 
  full_join(mortality_age %>% 
              group_by(iso2) %>% 
              summarise(adjrate_total=sum(adjrate)))
aapm_country<-mortality_age %>% 
  mutate(aapm=adjrate/adjrate_total) %>% 
  select(iso2, cause_cat, aapm) %>% 
  spread(cause_cat, aapm) %>% 
  rename_at(-1, ~paste0("p_", ., "_ageadj"))
# calculate aapm for the the entire sample
pop<-population %>% inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(age) %>% 
  summarise(pop=sum(pop))
mortality_age<-cause_specific_red %>% inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(cause_cat=case_when(
    DTHCOD_FINAL2%in%c(20,380, 
                       420, 490, 530)~"cmnn",
    DTHCOD_FINAL2%in%c(600, 780)~"cancer",
    DTHCOD_FINAL2%in%c(1040,
                       790, 800, 810, 980, 1110, 1150, 1200, 1240, 1250, 1310, 1430, 1475)~"ncd",
    DTHCOD_FINAL2%in%c(1490)~"accident",
    DTHCOD_FINAL2%in%c(1560)~"violent"
  ),
  age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(iso2, SALID1, male, cause_cat, age) %>% 
  summarise(deaths=sum(deaths)) %>% 
  left_join(correction) %>% 
  mutate(deaths=deaths/correction) %>% 
  group_by(cause_cat, age) %>% 
  summarise(deaths=sum(deaths)) %>% 
  left_join(pop) %>% 
  left_join(stds) %>% 
  mutate(rate=deaths/pop,
         adjrate=rate*std_pop) %>% 
  group_by(cause_cat) %>% 
  summarise(adjrate=sum(adjrate))
mortality_age<-mortality_age %>%
  mutate(adjrate_total=sum(adjrate))
aapm_all<-mortality_age %>% 
  mutate(aapm=adjrate/adjrate_total,
         iso2="ZZALL") %>% 
  select(iso2, cause_cat, aapm) %>% 
  spread(cause_cat, aapm) %>% 
  rename_at(-1, ~paste0("p_", ., "_ageadj"))
# calculate % ill-defined
pop<-population %>% inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(SALID1, age) %>% 
  summarise(pop=sum(pop))
mortality<-cause_specific_nored %>% 
  inner_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(cause_cat=case_when(
    DTHCOD_GHE2%in%c(1610)~"illdefined_disease",
    DTHCOD_GHE2%in%c(1620)~"illdefined_injury",
    T~"non_illdefined"),
    age=ifelse(age>=age_limit, age_limit, age)) %>% 
  group_by(SALID1, male, cause_cat, age) %>% 
  summarise(deaths=sum(deaths)) %>% 
  left_join(correction) %>% 
  mutate(deaths=deaths/correction) %>% 
  group_by(SALID1, cause_cat, age) %>% 
  summarise(deaths=sum(deaths))
mortality_illdefined<-mortality %>% 
  group_by(SALID1, cause_cat) %>% 
  summarise(deaths=sum(deaths)) %>% 
  as_tibble() %>% 
  ungroup() %>% 
  spread(cause_cat, deaths) %>% 
  rowwise() %>% 
  mutate(total=sum(c(illdefined_disease, illdefined_injury, non_illdefined)),
         illdefined_total=sum(c(illdefined_disease, illdefined_injury)),
         p_illdefined_disease=illdefined_disease/total,
         p_illdefined_injury=illdefined_injury/total,
         p_illdefined_total=illdefined_total/total,
         p_non_illdefined=non_illdefined/total)

save(mortality_cause, mortality_cause_models,
     aapm,aapm_country,aapm_all,
     mortality_illdefined,
     cause_titles_coll, file="analytic files/MS9_outcome_data.rdata")
