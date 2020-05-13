rm(list=ls())
library(data.table)
library(tidyverse)
library(furrr)
library(rjags)
library(broom)
source("MS10_SALURBAL_Helper.R")
select<-dplyr::select
load("analytic files/all_data_mortality_population_corrected_level1_age14.RData")
load("analytic files/MS10_exposure_data.RData")
load("analytic files/l1s.rdata")
l1s<-l1s %>% mutate(country_name=as.character(country_name),
                    iso2=as.character(iso2))

dta<-dta_corrected %>% full_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(sex=ifelse(male==1, "M", "F"),
         age=ifelse(age>=75, 75, age)) %>% 
  #mutate(iso2=ifelse(iso2=="CR", "PA", iso2)) %>% 
  group_by(iso2, SALID1, age, sex) %>% 
  summarise(deaths_unc=sum(deaths),
            total=sum(pop)) %>% 
  mutate(rate=deaths_unc/total) %>% 
  left_join(correction %>% select(SALID1,sex,  a, b)) %>% 
  mutate(avg_corr=a/(a+b)) 
summary(dta)  


# JAGS model, by sex
# sex_var<-"M"
plan(multiprocess)
both_output_inmodel<-future_map(c("M","F"), function(sex_var){
  # JAGS
  bugs<-dta %>% arrange(SALID1) %>% 
    filter(sex==sex_var) %>% 
    select(iso2, SALID1, age, deaths_unc, total, a, b) %>% 
    rename(deaths=deaths_unc) %>% 
    as.data.frame
  # vector of ages
  ages<-as.character(unique(bugs$age))
  # matrix of deaths n by K
  deaths<-bugs %>% ungroup() %>% select(SALID1, age, deaths) %>% 
    spread(age, deaths)
  Y<-as.matrix(deaths[,ages])
  # -	the number of age groups (K)
  K<-ncol(Y)
  # matrix of population n by K
  pop<-bugs %>% ungroup() %>% select(SALID1, age, total) %>% 
    spread(age, total)
  n<-pop[,ages]
  # ids of cities
  city_ids<-pop$SALID1
  # -	the rate of undercounting (correction; vector for shapes of beta distribution)
  correction_a<-bugs %>% filter(!duplicated(SALID1)) %>% pull(a)
  correction_b<-bugs %>% filter(!duplicated(SALID1)) %>% pull(b)
  # -	the number of cities (Ncity)
  Ncity<-nrow(Y)
  # -	the number of countries (Ncountry) and ids
  country<-bugs[!duplicated(bugs$SALID1),"iso2"]
  country_name<-country[!duplicated(country)]
  country<-as.numeric(as.factor(country))
  Ncountry<-length(unique(country))
  # get starting values for age and country coefficients
  m<-glm(deaths~-1+as.factor(age)+as.factor(iso2)+offset(log(total)), 
         data=bugs, family="poisson")
  summary(m)
  betas<-summary(m)$coefficients
  beta.age<-as.numeric(betas[grepl("age", rownames(betas)),1])
  beta.country<-as.numeric(betas[grepl("iso2", rownames(betas)),1])
  beta.country<-c(0, beta.country)
  #beta.age
  #beta.country
  jags_model<-"model
      {
        for(i in 1:Ncity) { #loop through city
         for(k in 1:K) {   #loop through all ages (K=16)
            Y[i,k] ~ dpois(theta[i,k])
            theta[i,k] <- n[i,k] * lambda[i,k] * correction[i]
      	    # beta is the overall mean for each age
      	    # country is a vector of length Ncity with the id of the country to which the city belongs
      	    # betacountry is the deviation from the overall mean for each age in each country
            lambda[i,k] <- exp(beta[k] + z[i,k] + betacountry[country[i]])
          }
          correction[i] ~ dbeta(correction_a[i],correction_b[i])
          z[i,1] ~ dnorm(0, tau[1])
          for (k in 2:K){
      	    z[i,k]~dnorm(mu[i,k], tau[k])
      	    mu[i,k]<-rho*z[i,k-1]
          }
        }
        
        betacountry[1] ~ dnorm(0, 0.001)
        for (j in 2:Ncountry){
      	  betacountry[j]~dnorm(0, tau2)
        }
        
        rho~dbeta(18, 2)
        tau2~dgamma(0.1, 0.1)
        for(k in 1:K) {
      	  beta[k]~dnorm(0, 0.001)
      	  tau[k]~dgamma(0.1, 0.1)
        }
      }
      "
  jags = jags.model(textConnection(jags_model), 
                    data = list('Y' = as.matrix(Y), 
                                'n' = n,
                                'correction_a'=correction_a,
                                'correction_b'=correction_b,
                                'country'=country,
                                'Ncity'=Ncity,
                                'K'=K,
                                'Ncountry'=Ncountry),
                    inits=list(tau2=10,
                               rho=.9,
                               tau=rep(10, times=K),
                               correction=correction_a/(correction_a+correction_b),
                               betacountry=beta.country,
                               beta=beta.age,
                               z=array(0, dim=c(Ncity,K)),
                               .RNG.name="base::Wichmann-Hill",
                               .RNG.seed=333),
                    n.chains = 1, 
                    n.adapt =5000) 
  update(jags, 50000) 
  nsims=100000;nthin=nsims/1000
  output=coda.samples(jags,
                      c('lambda', 'correction'),
                      n.iter = nsims,thin=nthin)
  coefficients<-as.matrix(output[[1]])
  coefficients
})
# diagnosis
# get corrections for diagnosis
corrections<-map(both_output_inmodel, function(both_output){
  #set.seed(333)
  coefs<-as.data.frame(both_output)
  coefs %>% as_tibble %>%  mutate(id=1:n()) %>% 
    gather(varid, rate, -id) %>% 
    filter(grepl("correction", varid)) %>% 
    ungroup() %>% 
    # lambda[XXX,YY]
    # YY is the age group: 1 is 0-1, 2 is 1-4, 3 is 5-9...
    # XXX is the city
    mutate(city=as.numeric(substr(varid, 
                                  regexpr("\\[", varid)+1, regexpr("\\]", varid)-1))) %>% 
    select(-varid) 
})
# extract mortality rates
model_output<-map(both_output_inmodel, function(both_output){
  #set.seed(333)
  coefs<-as.data.frame(both_output)
  coefs %>% as_tibble %>%  mutate(id=1:n()) %>% 
    gather(varid, rate, -id) %>% 
    filter(grepl("lambda", varid)) %>% 
    ungroup() %>% 
    # lambda[XXX,YY]
    # YY is the age group: 1 is 0-1, 2 is 1-4, 3 is 5-9...
    # XXX is the city
    mutate(age=as.numeric(substr(varid, 
                                 regexpr("\\,", varid)+1, nchar(varid)-1)),
           age=ifelse(age==1, 0,
                      ifelse(age==2, 1,
                             (age-2)*5)),
           city=as.numeric(substr(varid, 
                                  regexpr("\\[", varid)+1, regexpr("\\,", varid)-1))) %>% 
    select(-varid) 
})
model_output[[1]]$sex<-"M"
model_output[[2]]$sex<-"F"

# getting original data back (to get totals mostly)
bugs<-dta %>% 
  select(iso2, SALID1, age, sex, total) %>% as.data.frame
ages<-as.character(unique(bugs$age))
pop<-bugs %>% ungroup() %>% select(SALID1, sex, age, total) %>% 
  spread(age, total)
poplong<-pop %>% gather(age, pop, -SALID1, -sex) %>% mutate(age=as.numeric(age))
city_ids<-data.frame(SALID1=pop$SALID1) %>% 
  filter(!duplicated(SALID1)) %>% 
  mutate(city=row_number())
other_data<-city_ids %>% full_join(poplong) %>% left_join(l1s)

# get number of deaths from the rate, 
# then get just the data we need, 
# and group to do paralel processing
model_output<-do.call(rbind, model_output) %>% 
  left_join(other_data) %>% 
  mutate(deaths=rate*pop) %>% 
  left_join(l1s) %>% 
  #filter(id%in%idstoget) %>% 
  select(SALID1, id, sex, rate, deaths, age, pop) %>% 
  group_by(SALID1) 


#diagnosing posteriors
# first correction factors
check<-corrections[[1]] %>% group_by(city) %>% 
  summarise(median=median(rate),
            iqr=diff(quantile(rate, probs=c(.25, .75))))
low_coverage<-check %>% arrange(median) %>% slice(1:20) %>% pull(city)
uncertain_coverage<-check %>% arrange(desc(iqr)) %>% slice(1:20) %>% pull(city)
ggplot(corrections[[1]] %>% 
         filter(city%in%low_coverage), aes(x=rate)) + 
  geom_histogram(bins=30)+facet_wrap(~city)
ggsave("diagnosis/correction_lowM.pdf")
ggplot(corrections[[1]] %>% 
         filter(city%in%uncertain_coverage), aes(x=rate)) + 
  geom_histogram(bins=30)+facet_wrap(~city)
ggsave("diagnosis/correction_wideM.pdf")
check<-corrections[[2]] %>% group_by(city) %>% 
  summarise(median=median(rate),
            iqr=diff(quantile(rate, probs=c(.25, .75))))
low_coverage<-check %>% arrange(median) %>% slice(1:20) %>% pull(city)
uncertain_coverage<-check %>% arrange(desc(iqr)) %>% slice(1:20) %>% pull(city)
ggplot(corrections[[2]] %>% 
         filter(city%in%low_coverage), aes(x=rate)) + 
  geom_histogram(bins=30)+facet_wrap(~city)
ggsave("diagnosis/correction_lowF.pdf")
ggplot(corrections[[2]] %>% 
         filter(city%in%uncertain_coverage), aes(x=rate)) + 
  geom_histogram(bins=30)+facet_wrap(~city)
ggsave("diagnosis/correction_wideF.pdf")
# then mortality rate
#.x<-model_output %>% filter(SALID1==101101);.y<-data.frame(SALID1==101101, stringsAsFactors=F)
options(scipen=999)
diagnosis<-model_output %>% 
  left_join(l1s) %>% 
  group_by(SALID1) %>% 
  group_map(~{
    p<-ggplot(.x,aes(x=id, y=rate*100000)) +
      geom_line(alpha=0.3)+
      scale_y_log10()+
      scale_x_continuous(breaks=c(0, 500, 1000))+
      annotation_logticks(sides="l") +
      facet_grid(sex~age) +
      labs(title=.y$SALID1)+
      theme_bw()
    p$SALID1<-.y$SALID1
    p
  })
pdf("diagnosis/all_cities.pdf", width=10, height=5)
diagnosis
dev.off()


groups<-group_split(model_output, keep=F)
keys<-group_keys(model_output)
# clean things up to make things easier to process
rm(list=setdiff(ls(), c("groups","keys", "l1s", "ax")))
select<-dplyr::select
gc()
source("MS10_SALURBAL_Helper.R")
plan(multiprocess)
#xx1<-groups[[1]];xx2<-1;.x<-xx1 %>% filter(id==1, sex=="M");.y<-data.frame(id=1, sex="M", stringAsFactors=F)
ale<-future_map2_dfr(groups,1:nrow(keys) , function(xx1, xx2){
  xx1 %>% group_by(id, sex) %>% 
    group_modify(~{
      library(broom)
      library(DemoTools)
      xx2<-keys[xx2,] %>% left_join(l1s, by=c("SALID1"))
      salid1<-xx2$SALID1
      if (.y$id==1) print(salid1)
      country_n<-unique(xx2$country_name)
      gomp<-gompertz_approximation(temp=.x %>% arrange(age), 
                                   max_original=75, min=45, max=85, just5age=T)
      ages<-as.numeric(gomp$age)
      rates<-as.numeric(gomp$final_rate)
      # abridged lifetable function from DemoTools package
      lt<-lt_abridged(nMx = rates, Age = ages, sex=tolower(.y$sex))
      lt<-as.data.frame(lt[which(lt$Age%in%c(0, 20, 40, 60)),"ex"])
      colnames(lt)<-"ale"
      lt$age<-c(0, 20, 40, 60)
      lt$SALID1<-salid1
      lt$r2<-unique(gomp$r2)
      lt
    })
})

save(ale, file="analytic files/Life_Expectancy_All_Iterations.rdata")

# estimate LE for each undercounting strategy (Murray, Hill)
rm(list=ls())
library(data.table)
library(tidyverse)
library(furrr)
library(rjags)
library(broom)
source("MS10_SALURBAL_Helper.R")
select<-dplyr::select
hmean<-function(a){
  1/mean(1/a)  
}
load("analytic files/all_data_mortality_population_corrected_level1_age14.RData")
load("analytic files/MS10_exposure_data.RData")
load("analytic files/l1s.rdata")
load("analytic files/undercounting_correction_bysex.rdata")
correction<-correction %>% 
  mutate(ages=case_when(
    grepl("ages_hill", type) ~ "ages_hill",
    grepl("ages_murray", type) ~ "ages_murray",
    grepl("ages_auto", type) ~ "ages_auto"
  )) %>% 
  # remove the main method (already used)
  filter(!grepl("ages_auto", type)) %>% 
  group_by(SALID1, ages, sex) %>% 
  summarise(phi2=hmean(ucnt),
            var2=var(ucnt)) %>% 
  rowwise() %>% 
  mutate(K2=(phi2*(1-phi2))/var2-1) %>%
  # for those that have all 1s, make K2 high (high certainty)
  mutate(K2=ifelse(is.nan(K2), 100000, K2)) %>% 
  mutate(a=K2*phi2,
         b=K2*(1-phi2)) %>% 
  select(-var2) %>% 
  mutate(b=ifelse(b==0, 1, b)) %>% 
  mutate(a=ifelse(b<1, a/b, a),
         b=ifelse(b<1, b/b, b)) %>% 
  rename(type=ages)
methods<-correction %>% pull(type) %>% unique
l1s<-l1s %>% mutate(country_name=as.character(country_name),
                    iso2=as.character(iso2))

dta_original<-dta_corrected %>% full_join(years) %>% 
  filter(year>=y1, year<=y2) %>% 
  ungroup() %>% 
  mutate(sex=ifelse(male==1, "M", "F"),
         age=ifelse(age>=75, 75, age)) %>% 
  #mutate(iso2=ifelse(iso2=="CR", "PA", iso2)) %>% 
  group_by(iso2, SALID1, age, sex) %>% 
  summarise(deaths_unc=sum(deaths),
            total=sum(pop)) %>% 
  mutate(rate=deaths_unc/total)
summary(dta_original)  
#ucnt_var<-methods[[1]]
ale_by_ucnt<-map(methods, function(ucnt_var){
  dta<-left_join(dta_original, correction %>% filter(type==ucnt_var) %>% select(-type))
  # JAGS model, by sex, same as above
  # sex_var<-"M"
  both_output_inmodel<-lapply(c("M","F"), function(sex_var){
    # JAGS
    bugs<-dta %>% arrange(SALID1) %>% 
      filter(sex==sex_var) %>% 
      select(iso2, SALID1, age, deaths_unc, total, a, b) %>% 
      rename(deaths=deaths_unc) %>% 
      as.data.frame
    # vector of ages
    ages<-as.character(unique(bugs$age))
    # matrix of deaths n by K
    deaths<-bugs %>% ungroup() %>% select(SALID1, age, deaths) %>% 
      spread(age, deaths)
    Y<-as.matrix(deaths[,ages])
    # -	the number of age groups (K)
    K<-ncol(Y)
    # matrix of population n by K
    pop<-bugs %>% ungroup() %>% select(SALID1, age, total) %>% 
      spread(age, total)
    n<-pop[,ages]
    # ids of cities
    city_ids<-pop$SALID1
    # -	the rate of undercounting (correction; vector for shapes of beta distribution)
    correction_a<-bugs %>% filter(!duplicated(SALID1)) %>% pull(a)
    correction_b<-bugs %>% filter(!duplicated(SALID1)) %>% pull(b)
    # -	the number of cities (Ncity)
    Ncity<-nrow(Y)
    # -	the number of countries (Ncountry) and ids
    country<-bugs[!duplicated(bugs$SALID1),"iso2"]
    country_name<-country[!duplicated(country)]
    country<-as.numeric(as.factor(country))
    Ncountry<-length(unique(country))
    # get starting values for age and country coefficients
    m<-glm(deaths~-1+as.factor(age)+as.factor(iso2)+offset(log(total)), 
           data=bugs, family="poisson")
    summary(m)
    betas<-summary(m)$coefficients
    beta.age<-as.numeric(betas[grepl("age", rownames(betas)),1])
    beta.country<-as.numeric(betas[grepl("iso2", rownames(betas)),1])
    beta.country<-c(0, beta.country)
    #beta.age
    #beta.country
    jags_model<-"model
      {
        for(i in 1:Ncity) { #loop through city
         for(k in 1:K) {   #loop through all ages (K=16)
            Y[i,k] ~ dpois(theta[i,k])
            theta[i,k] <- n[i,k] * lambda[i,k] * correction[i]
      	    # beta is the overall mean for each age
      	    # country is a vector of length Ncity with the id of the country to which the city belongs
      	    # betacountry is the deviation from the overall mean for each age in each country
            lambda[i,k] <- exp(beta[k] + z[i,k] + betacountry[country[i]])
          }
          correction[i] ~ dbeta(correction_a[i],correction_b[i])
          z[i,1] ~ dnorm(0, tau[1])
          for (k in 2:K){
      	    z[i,k]~dnorm(mu[i,k], tau[k])
      	    mu[i,k]<-rho*z[i,k-1]
          }
        }
        
        betacountry[1] ~ dnorm(0, 0.001)
        for (j in 2:Ncountry){
      	  betacountry[j]~dnorm(0, tau2)
        }
        
        rho~dbeta(18, 2)
        tau2~dgamma(0.1, 0.1)
        for(k in 1:K) {
      	  beta[k]~dnorm(0, 0.001)
      	  tau[k]~dgamma(0.1, 0.1)
        }
      }
      "
    jags = jags.model(textConnection(jags_model), 
                      data = list('Y' = as.matrix(Y), 
                                  'n' = n,
                                  'correction_a'=correction_a,
                                  'correction_b'=correction_b,
                                  'country'=country,
                                  'Ncity'=Ncity,
                                  'K'=K,
                                  'Ncountry'=Ncountry),
                      inits=list(tau2=10,
                                 rho=.9,
                                 tau=rep(10, times=K),
                                 correction=correction_a/(correction_a+correction_b),
                                 betacountry=beta.country,
                                 beta=beta.age,
                                 z=array(0, dim=c(Ncity,K)),
                                 .RNG.name="base::Wichmann-Hill",
                                 .RNG.seed=333),
                      n.chains = 1, #number of chains
                      n.adapt =5000) 
    update(jags, 50000) 
    nsims=100000;nthin=nsims/1000
    output=coda.samples(jags,
                        c('lambda', 'correction'),
                        n.iter = nsims,thin=nthin)
    coefficients<-as.matrix(output[[1]])
    coefficients
  })
  # extract mortality rates
  model_output<-map(both_output_inmodel, function(both_output){
    #set.seed(333)
    coefs<-as.data.frame(both_output)
    coefs %>% as_tibble %>%  mutate(id=1:n()) %>% 
      gather(varid, rate, -id) %>% 
      filter(grepl("lambda", varid)) %>% 
      ungroup() %>% 
      # lambda[XXX,YY]
      # YY is the age group: 1 is 0-1, 2 is 1-4, 3 is 5-9...
      # XXX is the city
      mutate(age=as.numeric(substr(varid, 
                                   regexpr("\\,", varid)+1, nchar(varid)-1)),
             age=ifelse(age==1, 0,
                        ifelse(age==2, 1,
                               (age-2)*5)),
             city=as.numeric(substr(varid, 
                                    regexpr("\\[", varid)+1, regexpr("\\,", varid)-1))) %>% 
      select(-varid) 
  })
  model_output[[1]]$sex<-"M"
  model_output[[2]]$sex<-"F"
  
  # getting original data back (to get totals mostly)
  bugs<-dta %>% 
    select(iso2, SALID1, age, sex, total) %>% as.data.frame
  ages<-as.character(unique(bugs$age))
  pop<-bugs %>% ungroup() %>% select(SALID1, sex, age, total) %>% 
    spread(age, total)
  poplong<-pop %>% gather(age, pop, -SALID1, -sex) %>% mutate(age=as.numeric(age))
  city_ids<-data.frame(SALID1=pop$SALID1) %>% 
    filter(!duplicated(SALID1)) %>% 
    mutate(city=row_number())
  other_data<-city_ids %>% full_join(poplong) %>% left_join(l1s)
  
  # get number of deaths from the rate, 
  # then get just the data we need, 
  # and group to do paralel processing
  model_output<-do.call(rbind, model_output) %>% 
    left_join(other_data) %>% 
    mutate(deaths=rate*pop) %>% 
    left_join(l1s) %>% 
    #filter(id%in%idstoget) %>% 
    select(SALID1, id, sex, deaths, age, pop) %>% 
    group_by(SALID1) 
  groups<-group_split(model_output, keep=F)
  keys<-group_keys(model_output)
  plan(multiprocess)
  #xx1<-groups[[1]];xx2<-1;.x<-xx1 %>% filter(id==1, sex=="M");.y<-data.frame(id=1, sex="M", stringAsFactors=F)
  ale<-future_map2_dfr(groups,1:nrow(keys) , function(xx1, xx2){
    xx1 %>% group_by(id, sex) %>% 
      group_modify(~{
        library(broom)
        library(DemoTools)
        xx2<-keys[xx2,] %>% left_join(l1s, by=c("SALID1"))
        salid1<-xx2$SALID1
        if (.y$id==1) print(salid1)
        country_n<-unique(xx2$country_name)
        gomp<-gompertz_approximation(temp=.x %>% arrange(age), 
                                     max_original=75, min=45, max=85, just5age=T)
        ages<-as.numeric(gomp$age)
        rates<-as.numeric(gomp$final_rate)
        # abridged lifetable function from DemoTools package
        lt<-lt_abridged(nMx = rates, Age = ages, sex=tolower(.y$sex))
        lt<-as.data.frame(lt[which(lt$Age%in%c(0, 20, 40, 60)),"ex"])
        colnames(lt)<-"ale"
        lt$age<-c(0, 20, 40, 60)
        lt$SALID1<-salid1
        lt$r2<-unique(gomp$r2)
        lt
      })
  })
  ale$method<-ucnt_var
  ale
})

save(ale_by_ucnt, methods, file="analytic files/Life_Expectancy_by_UCNT.rdata")
