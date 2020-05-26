# load harmonized SALRUBAL data and save analytic file
rm(list=ls())
library(data.table)
library(tidyverse)
# get undercounting factors, and create parameters of beta distribution
source("MS10_SALURBAL_Helper.R")
# get population and mortality
country_list<-c("AR", "BR", "CL", "CO", "CR", "GT", "MX", "PA", "PE", "SV")
mortality_dir<-"../SALURBAL_DATA/Mortality Data/DTH/All Cause/"
mortality_files<-list.files(mortality_dir, pattern="L1")
pop_dir<-"../SALURBAL_DATA/Population Data/L1AD/Age Category/"
pop_files<-list.files(pop_dir, pattern="L1")
#country<-"CR"
# load country by country
dta_corrected<-map_dfr(country_list, function(country){
  print(country)
  file<-mortality_files[grepl(country, mortality_files)]
  mortality<-fread(paste0(mortality_dir, file),stringsAsFactors = F)
  file<-pop_files[grepl(country, pop_files)]
  population<-fread(paste0(pop_dir, file),stringsAsFactors = F)
  # Get single-age population for Costa Rica to get 0-1 and 1-4
  if (country=="CR"){
    population<-fread("../SALURBAL_DATA/Population Data/L1AD/Single Age/PRJSCR_L1AD_201905014.csv")
  }
  avail_years<-unique(mortality$YEAR)
  avail_years<-avail_years[avail_years%in%unique(population$YEAR)]
  population<-population[population$PRJAGE!=999,]
  agemax<-max(population$PRJAGE)
  mortality[mortality$DTHAGE5C>=agemax, "DTHAGE5C"]<-agemax
  mortality$year<-mortality$YEAR
  mortality$age<-mortality$DTHAGE5C
  # make sure we have 5-year age categories except for 1-4
  mortality$age<-ifelse(mortality$age==1, mortality$age, floor(mortality$age/5)*5)
  mortality$male<-mortality$DTHMALE
  population$year<-population$YEAR
  population$age<-population$PRJAGE
  # make sure we have 5-year age categories except for 1-4
  population$age<-ifelse(population$age==1, population$age, floor(population$age/5)*5)
  population$male<-population$PRJMALE
  mortality<-mortality %>% group_by(SALID1, year, age, male) %>% 
    summarise(deaths=sum(DTHDEATHS))
  population<-population %>% group_by(SALID1, year, age, male) %>% 
    summarise(pop=sum(PRJL1ADPOP))
  both<-
    # missing deaths or pop are structural 0s
    both<-full_join(mortality, population) %>% 
    filter(year%in%avail_years) %>% 
    mutate(deaths=replace_na(deaths, 0),
           pop=replace_na(pop, 0),
           iso2=country)
  both
})
# exclude non-SALURBAL areas
dta_corrected<-dta_corrected %>% filter(substr(SALID1, 4, 6)!="888")
summary(dta_corrected)

# get total population (unmatched with mortality)
population<-map_dfr(country_list, function(country){
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
summary(population)

# load correciton factors
load("analytic files/undercounting_correction_bysex.rdata")
correction<-correction %>% 
  group_by(SALID1, ages, sex) %>% 
  summarise(phi2=hmean(ucnt),
            var2=hmean_var(ucnt)) %>% 
  rowwise() %>% 
  mutate(K2=(phi2*(1-phi2))/var2-1) %>%
  # for those that have all 1s, make K2 high (high certainty)
  mutate(K2=ifelse(is.nan(K2), 10000, K2)) %>% 
  mutate(a=K2*phi2,
         b=K2*(1-phi2)) %>% 
  mutate(b=ifelse(b==0, 1, b)) %>% 
  mutate(a=ifelse(b<1, a/b, a),
         b=ifelse(b<1, b/b, b)) %>% 
  select(SALID1, ages, sex, a, b)
save(dta_corrected, correction, population, file="analytic files/all_data_mortality_population_corrected_level1_age14.RData")

## save database of exposures
rm(list=ls())
library(tidyverse)
library(data.table)
library(readxl)
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
pop_l1<-population %>% right_join(years) %>% 
  filter(year==y1 | year==y2) %>% 
  mutate(year=ifelse(year==y1,"y1", "y2")) %>% 
  spread(year, total_pop) %>% 
  mutate(growth=y2-y1,
         growth_pct=growth/y1,
         pop_baseline=y1) %>% 
  select(SALID1, growth, growth_pct, pop_baseline)
pop_growth_prev<-population %>% right_join(years) %>% 
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
bec<-read.csv("../SALURBAL_DATA/BEC Data/BEC_L1AD_20191031.csv", stringsAsFactors = F)
bec$total_transport<-rowSums(bec[,colnames(bec)[grepl("BECPR", colnames(bec))]])
bec$total_transport_yn<-as.numeric(bec$total_transport>=1)
bec<-bec[,c("SALID1", colnames(bec)[grepl("BEC", colnames(bec))],"total_transport", "total_transport_yn")]
becux<-read.csv("../SALURBAL_DATA/BEC Data/BEC_L1UX_20191031.csv", stringsAsFactors = F)
becux<-becux[,c("SALID1", colnames(becux)[grepl("BEC", colnames(becux))])]
becux<-becux[,!grepl("GASPRICE", colnames(becux))]
bec<-full_join(bec, becux)
head(bec)


## SEC variables
sec<-fread("../SALURBAL_DATA/SEC Data/SEC_Census_L1AD_03162020.csv") %>% 
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

save(all_exposure,years, file="analytic files/MS10_exposure_data.RData")

# Save SHP for the App
rm(list=ls())
library(sf)
library(rmapshaper)
shp = st_read('../SALURBAL_DATA/SHPs/SALURBAL_L1AD_11_30_18/SALURBAL_L1AD_11_30_18.shp')
shp <- st_transform(shp, "+init=epsg:4326") %>% ms_simplify()
save(shp, file="MS10/l1ad_shp.rdata")
