# First, life expectancy part of paper
rm(list=ls())
library(data.table)
library(tidyverse)
library(broom)
library(readxl)
library(furrr)
library(broom.mixed)
library(scales)
library(lme4)
library(lspline)
library(grid)
library(gridExtra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(rmapshaper)
library(classInt)
run_models=F
source("MS10_SALURBAL_Helper.R")
select<-dplyr::select
load("analytic files/all_data_mortality_population_corrected_level1_age14.RData")
load("analytic files/l1s.rdata")
load("analytic files/MS9_outcome_data.rdata")
load("analytic files/Life_Expectancy_All_Iterations.rdata")
load("analytic files/MS10_exposure_data.RData")
ale<-ale %>% filter(age%in%c(0, 20, 40, 60))
# save app data from LE part
# need: SALID1, age, sex, le, lci, uci, iso2, city_link, country_name
dta_le<-ale %>% 
  group_by(SALID1, sex, age) %>% 
  summarise(le=median(ale),
            lci=quantile(ale, probs=0.025),
            uci=quantile(ale, probs=0.975)) %>% 
  left_join(l1s)
save(dta_le, file="MS10/data_LE.rdata")
file.copy(from="analytic files/l1s.rdata", 
          to="MS10/l1s.rdata", overwrite = T)
file.copy(from="analytic files/MS9_outcome_data.rdata", 
          to="MS10/MS9_outcome_data.rdata", overwrite = T)
file.copy(from="analytic files/MS10_exposure_data.RData", 
          to="MS10/MS10_exposure_data.RData", overwrite = T)


# median LE in long format for some of the descriptives
ale_median<-ale %>% group_by(SALID1, sex, age) %>% 
  summarise(le=median(ale),
            lci=quantile(ale, probs=0.025),
            uci=quantile(ale, probs=0.975),
            se=sqrt(var(ale)),
            dif_ci=uci-lci) %>% 
  mutate(rse=se/le*100) %>% 
  left_join(l1s)

# FIGURE 1
p1<-ale_median %>% 
  filter(age==0) %>% 
  group_by(age, sex) %>% 
  group_map(~{
    age<-.y$age
    sex<-.y$sex
    ylim<-c(min(floor(ale_median[ale_median$age==age, "le"]/5)*5),
            max(ceiling(ale_median[ale_median$age==age, "le"]/5)*5))
    ylab<-ifelse(age==0, "Life Expectancy at Birth", paste0("Life Expectancy at Age ", age))
    title<-ifelse(sex=="M", "Men", "Women")
    ylab2<-ifelse(sex=='M', ylab, "")
    ylab<-ifelse(sex=='M', "", ylab)
    x<-.x %>% left_join(l1s)
    ggplot(x, aes(x=iso2, y=le, group=iso2)) +
      geom_boxplot(aes(group=as.factor(iso2)), fill=NA, outlier.color = NA, width=0.5)+
      geom_jitter(aes(fill=as.factor(iso2)), width=0.1, height=0, alpha=1, size=2, 
                  color="black", pch=21) +
      guides(color=F, fill=F, size=F)+
      labs(x="",
           y=ylab,
           title=title)+
      scale_y_continuous(sec.axis=dup_axis(name = ylab2), limits=ylim)+
      theme_bw() +
      theme(legend.position = "bottom",
            legend.key.width = unit(50, "points"),
            panel.grid.major.x = element_blank(),
            axis.text.x=element_text(size=20, color="black"),
            axis.text.y=element_text(size=16, color="black"),
            axis.title.y=element_text(face="bold", size=20),
            plot.title=element_text(face="bold", size=25))
  })
pLEB<-arrangeGrob(grobs=p1[1:2], ncol=2)
plot(pLEB)
ggsave("results/Figure1_color.pdf", pLEB, width=15, height=5)
# re-doing figure in BnW
p1<-ale_median %>% 
  filter(age==0) %>% 
  group_by(age, sex) %>% 
  group_map(~{
    age<-.y$age
    sex<-.y$sex
    ylim<-c(min(floor(ale_median[ale_median$age==age, "le"]/5)*5),
            max(ceiling(ale_median[ale_median$age==age, "le"]/5)*5))
    ylab<-ifelse(age==0, "Life Expectancy at Birth", paste0("Life Expectancy at Age ", age))
    title<-ifelse(sex=="M", "Men", "Women")
    ylab2<-ifelse(sex=='M', ylab, "")
    ylab<-ifelse(sex=='M', "", ylab)
    x<-.x %>% left_join(l1s)
    ggplot(x, aes(x=iso2, y=le, group=iso2)) +
      geom_boxplot(aes(group=as.factor(iso2)), fill=NA, outlier.color = NA, width=0.5)+
      geom_jitter(fill="gray", width=0.1, height=0, alpha=1, size=2, 
                  color="black", pch=21) +
      guides(color=F, fill=F, size=F)+
      labs(x="",
           y=ylab,
           title=title)+
      scale_y_continuous(sec.axis=dup_axis(name = ylab2), limits=ylim)+
      theme_bw() +
      theme(legend.position = "bottom",
            legend.key.width = unit(50, "points"),
            panel.grid.major.x = element_blank(),
            axis.text.x=element_text(size=20, color="black"),
            axis.text.y=element_text(size=16, color="black"),
            axis.title.y=element_text(face="bold", size=20),
            plot.title=element_text(face="bold", size=25))
  })
pLEB<-arrangeGrob(grobs=p1[1:2], ncol=2)
plot(pLEB)
ggsave("results/Figure1.pdf", pLEB, width=15, height=5)
# save source data
ale_median %>% 
  filter(age==0) %>% ungroup() %>% select(sex, le, iso2) %>% 
  fwrite("Source/Source_Data_Figure1.csv")

# other LE descriptives for Supplement, Appendix, discussion, results, etc.
# first, descriptives on ranges, highest, lowest ciites, etc.
# range
# get ranges and medians by country
ranges<-ale_median %>% 
  left_join(l1s) %>% 
  group_by(iso2, age, sex) %>%  
  summarise(min=min(le),
            max=max(le),
            range=max(le)-min(le),
            median=median(le)) %>% 
  filter(age==0)
ranges_overall<-ale_median %>% 
  left_join(l1s) %>% 
  group_by(age, sex) %>%  
  summarise(min=min(le),
            max=max(le),
            range=max(le)-min(le),
            median=median(le)) %>% 
  filter(age==0)
ranges[order(ranges$sex, ranges$range),]
ranges[order(ranges$sex, ranges$median),]
ranges_overall
# top cities for discussion
ale_median %>% filter(age==0) %>% left_join(l1s) %>%
  group_by(sex) %>% arrange(desc(le)) %>% slice(1:6) %>% select(sex, le, lci, uci, iso2, city_link)
ale_median %>% filter(age==0) %>% left_join(l1s) %>%
  group_by(sex) %>% arrange((le)) %>% slice(1:6) %>% select(sex, le, lci, uci, iso2, city_link)

# comparison with other countries
country_list<-l1s %>% filter(!duplicated(iso2)) %>% select(country_name, iso2) %>%
  rename(country=country_name)

dta1<-(read_excel("Other_Data/UNDP Life Tables/WPP2019_MORT_F17_2_ABRIDGED_LIFE_TABLE_MALE.xlsx", skip=16)) %>%
  mutate(sex="M")
dta2<-(read_excel("Other_Data/UNDP Life Tables/WPP2019_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE.xlsx", skip=16)) %>%
  mutate(sex="F")

country_les<-map_dfr(list(dta1, dta2), function(temp){
  sexvar<-unique(temp$sex)
  temp<-temp %>%
    filter(Type=="Country/Area") %>% 
    rename(country=`Region, subregion, country or area *`,
           age=`Age (x)`,
           le_country=`Expectation of life e(x)`) %>%
    mutate(le_country=as.numeric(le_country),
           sex=sexvar) %>%
    mutate(y1=as.numeric(substr(Period, 1, 4)),
           y2=as.numeric(substr(Period, 6, 9))) %>% 
    rowwise() %>% 
    mutate(year=mean(c(y1, y2))) %>% 
    #filter(Period=="2010-2015") %>%
    filter(age%in%c(0, 20, 40, 60)) %>%
    select(country, year, age,sex,  le_country) %>% 
    group_by(country, age) %>% 
    group_modify(~{
      fit<-lm(le_country ~ lspline(year,knots = seq(1952.5, 2017.5, by=5), marginal=T),data = .x )
      le_list<-expand.grid(year=1950:2020)
      le_list$le_country<-predict(fit, newdata=le_list)
      le_list
    }) %>% 
    mutate(sex=sexvar)
  temp
})  %>% filter(year%in%c(2012:2016)) %>% 
  group_by(country, age, sex) %>% 
  summarise(le_country=mean(le_country))

country_les %>% 
  filter(age==0, sex=="F", le_country>(ale_median %>% 
           filter(age==0, sex=="F") %>% 
           pull(le) %>% max)) %>% 
  arrange(le_country) %>% print(n=20)
country_les %>% 
  filter(age==0, sex=="M", le_country>(ale_median %>% 
                                         filter(age==0, sex=="M") %>% 
                                         pull(le) %>% max)) %>% 
  arrange(le_country) %>% print(n=20)

country_les %>% 
  filter(age==0, sex=="F", le_country<(ale_median %>% 
                                         filter(age==0, sex=="F") %>% 
                                         pull(le) %>% min)) %>% 
  arrange(desc(le_country)) %>% print(n=20)
country_les %>% 
  filter(age==0, sex=="M", le_country<(ale_median %>% 
                                         filter(age==0, sex=="M") %>% 
                                         pull(le) %>% min)) %>% 
  arrange(desc(le_country)) %>% print(n=20)

p1<-ale_median %>% 
  left_join(country_les %>% rename(country_name=country)) %>% 
  group_by(age, sex) %>% 
  group_map(~{
    age<-.y$age
    sex<-.y$sex
    ylim<-c(min(floor(ale_median[ale_median$age==age, "le"]/5)*5),
            max(ceiling(ale_median[ale_median$age==age, "le"]/5)*5))
    ylab<-ifelse(age==0, "Life Expectancy at Birth", paste0("Life Expectancy at Age ", age))
    title<-ifelse(sex=="M", "Men", "Women")
    ylab2<-ifelse(sex=='M', ylab, "")
    ylab<-ifelse(sex=='M', "", ylab)
    x<-.x %>% left_join(l1s)
    ggplot(x, aes(x=iso2, y=le, group=iso2)) +
      geom_boxplot(aes(group=(iso2)), fill=NA, outlier.color = NA, width=0.5)+
      geom_jitter(aes(fill=(iso2)), width=0.075, height=0, alpha=1, size=2, 
                  color="black", pch=21) +
      geom_point(data=.x %>% filter(!duplicated(iso2)),
                 aes(x=iso2, y=le_country, fill=iso2),
                 position=position_nudge(x=+0.3),
                 pch=24, color="black",size=4)+
      guides(color=F, fill=F, size=F)+
      labs(x="",
           y=ylab,
           title=title)+
      scale_y_continuous(sec.axis=dup_axis(name = ylab2), limits=ylim)+
      theme_bw() +
      theme(legend.position = "bottom",
            legend.key.width = unit(50, "points"),
            panel.grid.major.x = element_blank(),
            axis.text.x=element_text(size=20, color="black"),
            axis.text.y=element_text(size=16, color="black"),
            axis.title.y=element_text(face="bold", size=20),
            plot.title=element_text(face="bold", size=25))
  })
p<-arrangeGrob(grobs=p1, ncol=2)
ggsave("results/ExtendedData1.pdf", p, width=15, height=5*4)
ggsave("results/ExtendedData1.eps", p, width=15, height=5*4)
# save source data
ale_median %>% 
  left_join(country_les %>% rename(country_name=country)) %>%
  ungroup() %>% 
  select(sex, age, le, iso2, le_country) %>% 
  fwrite("Source/Source_Data_ExtendedData1.csv")

# Extended Data 2
# comparison with other income groups
dta1<-(read_excel("Other_Data/UNDP Life Tables/WPP2019_MORT_F17_2_ABRIDGED_LIFE_TABLE_MALE.xlsx", skip=16)) %>%
  mutate(sex="M")
dta2<-(read_excel("Other_Data/UNDP Life Tables/WPP2019_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE.xlsx", skip=16)) %>%
  mutate(sex="F")
income_les<-map_dfr(list(dta1, dta2), function(temp){
  sexvar<-unique(temp$sex)
  temp<-temp %>%
    filter(Type=="Income Group") %>% 
    rename(country=`Region, subregion, country or area *`,
           age=`Age (x)`,
           le_country=`Expectation of life e(x)`) %>%
    mutate(le_country=as.numeric(le_country),
           sex=sexvar) %>%
    mutate(y1=as.numeric(substr(Period, 1, 4)),
           y2=as.numeric(substr(Period, 6, 9))) %>% 
    rowwise() %>% 
    mutate(year=mean(c(y1, y2))) %>% 
    #filter(Period=="2010-2015") %>%
    filter(age%in%c(0, 20, 40, 60)) %>%
    select(country, year, age,sex,  le_country) %>% 
    group_by(country, age) %>% 
    group_modify(~{
      fit<-lm(le_country ~ lspline(year,knots = seq(1952.5, 2017.5, by=5), marginal=T),data = .x )
      le_list<-expand.grid(year=1950:2020)
      le_list$le_country<-predict(fit, newdata=le_list)
      le_list
    }) %>% 
    mutate(sex=sexvar)
  temp
})  %>% filter(year%in%c(2012:2016)) %>% 
  group_by(country, age, sex) %>% 
  summarise(le_country=mean(le_country))

income_labels<-c("High-income countries", "Upper-middle-income countries",
                 "Middle-income countries","Lower-middle-income countries",
                 "Low-income countries")
#.x<-ale_median %>% filter(sex=="M", age==0);.y<-data.frame(sex="M", age=0, stringsAsFactors=F)
figure2<-ale_median %>% 
  filter(age==0) %>% 
  group_by(age, sex) %>% 
  group_map(~{
    incometemp<-income_les %>% filter(sex==.y$sex, age==.y$age) %>% 
      filter(!grepl("No income", country)) %>% arrange(le_country) %>% 
      ungroup() %>% 
      mutate(country=factor(country, levels=income_labels))
    .x<-.x %>% arrange(le) %>% 
      ungroup() %>% 
      mutate(id=row_number())
    title<-paste0(ifelse(.y$age==0, "Life Expectancy at Birth in ", 
                         paste0("Life Expectancy at Age ", .y$age, " in ")),
                  ifelse(.y$sex=="M", "Men", "Women"))
    title<-ifelse(.y$sex=="M", "Men", "Women")
    ylim1<-pmin(ale_median %>% filter(age==0) %>% pull(lci) %>% min,
                income_les %>% filter(age==0) %>% pull(le_country) %>% min)
    ylim2<-pmax(ale_median %>% filter(age==0) %>% pull(uci) %>% max,
                income_les %>% filter(age==0) %>% pull(le_country) %>% max)
    ylim<-c(ylim1, ylim2)
    ggplot(.x, aes(x=id, y=le)) +
      geom_errorbar(aes(ymin=lci, ymax=uci))+
      geom_point() +
      geom_hline(data=incometemp, lty=2, size=2,
                 aes(yintercept = le_country, color=country))+
      scale_x_continuous(breaks=.x$id, labels=.x$city_link, expand=c(0.01,0.01))+
      scale_y_continuous(limits=ylim)+
      scale_color_manual(values=brewer_pal(type="qual", palette=2)(5), name="")+
      guides(color=guide_legend(override.aes = list(size=2, linetype=1)))+
      labs(title=title, y="Life Expectancy (years)", x="")+
      theme_bw()+
      theme(legend.position = "bottom",
            axis.text.y=element_text(face="bold", size=14, color="black"),
            axis.title=element_text(face="bold", size=16, color="black"),
            panel.grid = element_blank(),
            axis.ticks.x=element_blank(),
            plot.title=element_text(face="bold", size=24, color="black"),
            axis.text.x=element_text(size=3, angle=90, 
                                     hjust=1, color="black"),
            legend.text=element_text(size=20, color="black"),
            legend.title=element_text(face="bold", size=20, color="black"))+
      theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
  })
legend<-get_legend(figure2[[1]])
figure2<-lapply(figure2, function(xx) xx+guides(color=F, fill=F))
figure2<-arrangeGrob(grobs=list(figure2[[1]], figure2[[2]]), ncol=2)
figure2<-arrangeGrob(grobs=list(figure2, legend), nrows=2,
                     heights=c(10, 1))
plot(figure2)
ggsave("results/ExtendedData2.pdf", figure2, width=20, height=11.4)
ggsave("results/ExtendedData2.eps", figure2, width=20, height=11.4)
# source data
ale_median %>% 
  filter(age==0) %>% 
  ungroup() %>% 
  select(iso2, sex, age, le, lci, uci) %>% 
  left_join(income_les %>% filter(age==0) %>% spread(country, le_country)) %>% 
  fwrite("Source/Source_Data_ExtendedData2.csv")
fread("Source/Source_Data_ExtendedData2.csv") %>% head

# now explore the uncertainty
# 95 CI ranges
ale_median %>% 
  filter(age==0) %>% 
  pull(dif_ci) %>% 
  quantile(probs=c(0.01, 0.025, 0.05, 0.25, .5, .75, .95, .975, .99))
ale_median %>% group_by(age, sex) %>% 
  summarise(b5=sum(dif_ci<5),
            total=n()) %>% 
  mutate(prop_b5=b5/total)

# RSE distribution and Figure
table(ale_median %>% filter(age==0, sex=="F") %>% pull(rse)<2)

sex_label<-c("Women", "Men")
names(sex_label)<-c("F", "M")
age_label<-paste0("LE at ", c("Birth", paste0("age ", c(20, 40, 60))))
names(age_label)<-c(0, 20, 40, 60)
figure_rse<-ggplot(ale_median, aes(x=iso2, y=rse, group=iso2)) +
  #geom_boxplot(aes(group=as.factor(iso2)), fill=NA, outlier.color = NA, width=0.5)+
  geom_jitter(aes(fill=as.factor(iso2)), width=0.1, height=0, alpha=1, size=2, 
              color="black", pch=21) +
  guides(color=F, fill=F, size=F)+
  labs(x="",
       y="Relative Standard Error (%)",
       title="")+
  #scale_y_continuous(sec.axis=dup_axis(name = ylab2), limits=ylim)+
  facet_grid(sex~age, labeller=labeller(age=age_label, sex=sex_label))+
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(50, "points"),
        panel.grid.major.x = element_blank(),
        axis.text.x=element_text(size=20, color="black"),
        axis.text.y=element_text(size=16, color="black"),
        axis.title.y=element_text(face="bold", size=20),
        strip.background = element_blank(),
        strip.text=element_text(face="bold", size=20),
        plot.title=element_text(face="bold", size=25))
plot(figure_rse)
ggsave("results/SupplementaryFigure52.pdf", figure_rse, width=20, height=10)

# plot undercounting
load("analytic files/undercounting_correction_bysex.rdata")
correction<-correction %>% 
  group_by(SALID1, sex) %>% 
  summarise(phi2=hmean(ucnt)) %>% 
  left_join(l1s)
ylim<-c(min(correction$phi2), 1)
sex_label<-c("Women", "Men")
names(sex_label)<-c("F", "M")
ggplot(correction, aes(x=iso2, y=phi2, group=iso2)) +
  geom_hline(yintercept = 1, lty=2)+
  #geom_boxplot(aes(group=as.factor(iso2)), fill=NA, outlier.color = NA, width=0.5)+
  geom_jitter(aes(fill=as.factor(iso2)), 
              width=0.1, height=0, alpha=1, pch=21, color="black", size=3) +
  guides(fill=F, size=F)+
  labs(x="", y="Completeness")+
  scale_y_continuous(trans = log_trans(), 
                     breaks = base_breaks(5),
                     labels=percent,
                     limits = ylim)+
  facet_wrap(~sex, labeller = labeller(sex=sex_label))+
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(50, "points"),
        axis.text.x=element_text(size=20, color="black"),
        axis.text.y=element_text(size=14, color="black"),
        axis.title.y=element_text(size=14, color="black", face="bold"),
        plot.title=element_text(face="bold", size=25),
        strip.text =element_text(size=20, color="black", face="bold"),
        strip.background = element_blank())
ggsave("results/SupplementaryFigure55.pdf", width=10, height=5)

# total variability
# ICCs taking into consideration variability between iterations
iccs<-ale %>% 
  group_by(age, sex) %>% 
  group_modify(~{
    .x<-.x %>% 
      left_join(all_exposure) 
    m<-lmer(ale~1+(1|iso2/SALID1), data=.x, weights=pop_baseline/1000000)
    vars<-as.data.frame(VarCorr(m)) %>% select(vcov) %>% 
      mutate(level=c("city", "country", "iteration")) %>% 
      spread(level, vcov) %>% 
      mutate(ICCiteration=iteration/(country+city+iteration),
             ICCcity=city/(country+city+iteration),
             ICCcountry=country/(country+city+iteration))
  })
iccs<-iccs %>% select(age, sex, ICCiteration,ICCcity, ICCcountry) %>% 
  mutate(ICCiteration=paste0(format(ICCiteration*100, nsmall=1, digits=1), "%"),
         ICCcity=paste0(format(ICCcity*100, nsmall=1, digits=1), "%"),
         ICCcountry=paste0(format(ICCcountry*100, nsmall=1, digits=1), "%"))
iccs<-full_join(iccs %>%  select(age, sex, ICCiteration) %>% 
                  spread(sex, ICCiteration) %>% 
                  rename(IterationF=F, IterationM=M),
                iccs %>%  select(age, sex, ICCcity) %>% 
                  spread(sex, ICCcity) %>% 
                  rename(CityF=F, CityM=M)) %>% 
  full_join(iccs %>%  select(age, sex, ICCcountry) %>%  
              spread(sex, ICCcountry) %>% 
              rename(CountryF=F, CountryM=M)) %>% 
  select(age, IterationF,CityF, CountryF, IterationM, CityM, CountryM)
iccs
fwrite(iccs, file="results/SupplementTable51.csv")


# Outcome models
vars_cont<-vars<-c("pop_baseline", "growth_pct", "BECPOPDENSL1AD","BECPTCHDENSL1AD", "BECADINTDENSL1AD",
                   "CNSMINPR_L1AD", "CNSWATINL1AD", "CNSSEWNETL1AD", "CNSCROWD3RML1AD", 
                   "sei")
# get standard deviations for table
vars_sds<-vars_cont[-c(1, length(vars_cont))]
sds<-all_exposure %>% ungroup() %>% 
  mutate(growth_pct=growth_pct*100) %>% 
  select(all_of(vars_sds)) %>% summarise_all(sd, na.rm=T) %>% 
  as.numeric()
sds_char<-as.character(sds)
sds_char[which(vars_sds=="growth_pct")]<-paste0(round(sds[which(vars_sds=="growth_pct")], digits=1), "%")
sds_char[grep("CNS", vars_sds)]<-paste0(round(sds[grep("CNS", vars_sds)], digits=1), "%")
sds_char[grep("BEC", vars_sds)]<-round(sds[grep("BEC", vars_sds)], digits=2)
sds_char[grep("BECPOPDENS", vars_sds)]<-round(sds[grep("BECPOPDENS", vars_sds)], digits=0)
# Units of population change? 50%
pop_coefficient<-1.5
sds_char<-c(paste0((pop_coefficient-1)*100, "%"), sds_char, "1 SD")
sds_char

ale_long<-full_join(ale, all_exposure %>% select(SALID1, iso2, BECPCTURBANL1AD, all_of(vars_cont)))
vars_multiv<-vars_cont[!grepl("CNS", vars_cont)]
#.x<-ale_long %>% filter(age==0, sex=="M", id==1)
# this runs LE models
if (run_models==T) {
  ale_long_grouped<-ale_long %>%group_by(sex, age, id) %>% 
    group_split(.keep=T )
  plan(multiprocess, .cleanup=T)
  models<-future_map_dfr(ale_long_grouped, function(xx){
    library(broom.mixed)
    library(lme4)
    temp<-xx
    temp<-temp %>% mutate(pop_baseline=log(pop_baseline*1000000))
    temp<-temp %>% mutate_at(vars_cont[-1], scale, center=T, scale=T) %>% 
      mutate_at(vars_cont[-1], as.numeric)
    # empty
    f<-as.formula(paste0("ale~", "+ BECPCTURBANL1AD + (1|iso2)"))
    model<-lmer(f,  data=temp)
    variances<-VarCorr(model)
    tau00_empty=variances %>% as.data.frame %>% filter(grp=="iso2") %>% pull(vcov)
    sigma2_empty=variances %>% as.data.frame %>% filter(grp=="Residual") %>% pull(vcov)
    # model 1
    model1<-map_dfr(vars_cont, function(var){
      f<-as.formula(paste0("ale~", var, "+ BECPCTURBANL1AD+(1|iso2)"))
      model<-lmer(f, data=temp, na.action = "na.omit")
      coefs<-model %>% tidy %>% filter(term==var)
      variances<-VarCorr(model)
      data.frame(var=var, 
                 coef=coefs %>% pull(estimate), 
                 se=coefs %>% pull(std.error), 
                 tau00=variances %>% as.data.frame %>% filter(grp=="iso2") %>% pull(vcov), 
                 sigma2=variances %>% as.data.frame %>% filter(grp=="Residual") %>% pull(vcov))   
    })
    # model 2
    f<-as.formula(paste0("ale~", paste(vars_multiv, collapse="+"), 
                         "+ BECPCTURBANL1AD+(1|iso2)"))
    model<-lmer(f, data=temp, na.action = "na.omit")
    coefs<-model %>% tidy %>% filter(term!="(Intercept)", 
                                     term!="BECPCTURBANL1AD",
                                     effect=="fixed")
    variances<-VarCorr(model)
    model2<-data.frame(var=vars_multiv,
                       coef_multiv=coefs %>% pull(estimate), 
                       se_multiv=coefs %>% pull(std.error), 
                       tau00_multiv=variances %>% as.data.frame %>% filter(grp=="iso2") %>% pull(vcov), 
                       sigma2_multiv=variances %>% as.data.frame %>% filter(grp=="Residual") %>% pull(vcov))   
    model1<-full_join(model1, model2, by="var")
    model1$tau00_empty<-tau00_empty
    model1$sigma2_empty<-sigma2_empty
    model1$id<-unique(xx$id)
    model1$age<-unique(xx$age)
    model1$sex<-unique(xx$sex)
    model1
  }, .progress=T)
  plan(multiprocess, .cleanup=T)
  # Sensitivity analysis: is there assoc with pop growth previous to study period?
  models_growth<-future_map_dfr(ale_long_grouped, function(xx){
    library(broom.mixed)
    library(lme4)
    temp<-xx
    temp<-temp %>% left_join(all_exposure %>% select(SALID1,growth_pct_prev ))
    # exclude cities for which we dont have previous period data
    temp<-temp %>% filter(!is.na(growth_pct_prev))
    temp<-temp %>% mutate(pop_baseline=log(pop_baseline*1000000))
    temp<-temp %>% mutate_at(vars_cont[-1], scale, center=T, scale=T) %>% 
      mutate_at(vars_cont[-1], as.numeric)
    temp$growth_pct_prev<-scale(temp$growth_pct_prev, center=T, scale=T)
    f<-as.formula(paste0("ale~", "growth_pct", "+ BECPCTURBANL1AD+(1|iso2)"))
    model1<-lmer(f, data=temp, na.action = "na.omit")
    coefs1<-model1 %>% tidy %>% filter(term=="growth_pct")
    f<-as.formula(paste0("ale~", "growth_pct_prev", "+ BECPCTURBANL1AD+(1|iso2)"))
    model2<-lmer(f, data=temp, na.action = "na.omit")
    coefs2<-model2 %>% tidy %>% filter(term=="growth_pct_prev")
    temp<-bind_rows(coefs1, coefs2) %>% select(term, estimate, std.error) %>% 
      rename(coef=estimate, se=std.error)
    temp$id<-unique(xx$id)
    temp$age<-unique(xx$age)
    temp$sex<-unique(xx$sex)
    temp
  }, .progress=T)
  save(models, models_growth, file="results/model_results_le.rdata")
} else {
  load("results/model_results_le.rdata")  
}


results_multiv<-models %>% 
  filter(!grepl("CNS", var)) %>% 
  group_by(age, sex, var) %>% 
  summarise(beta=mean(coef_multiv),
            se=sqrt(mean(se_multiv^2)+var(coef_multiv))) %>% 
  # re-scale population coefficients
  mutate(beta=ifelse(var=="pop_baseline", beta*log(pop_coefficient), beta),
         se=ifelse(var=="pop_baseline", se*log(pop_coefficient), se)) %>% 
  ungroup() %>% 
  mutate(lci=beta-1.96*se,
         uci=beta+1.96*se,
         final=paste0(round(beta, digits=2),
                      " [",
                      round(lci, digits=2),
                      ";",
                      round(uci, digits=2),
                      "]"))

results_multiv_table<-results_multiv %>% 
  filter(age==0) %>% 
  ungroup() %>% 
  select(sex, var, final) %>% 
  spread(sex, final) %>% select(var, M, F) %>% 
  rename(multiM=M, multiF=F) %>% 
  arrange(factor(var, levels=vars_cont), desc(var))%>% 
  mutate(multiM=replace(multiM, is.na(multiM), ""),
         multiF=replace(multiF, is.na(multiF), "")) %>% 
  select(var,multiM, multiF)

table1<-results_multiv_table
table1
write.csv(table1, "results/Table1.csv")
sds_char[which(vars%in%vars_multiv)]
# extract coefficients for figure with all coefficients
results_univ<-models %>%
  group_by(age, sex, var) %>% 
  summarise(beta=mean(coef),
            se=sqrt(mean(se^2)+var(coef))) %>% 
  # re-scale population coefficients
  mutate(beta=ifelse(var=="pop_baseline", beta*log(pop_coefficient), beta),
         se=ifelse(var=="pop_baseline", se*log(pop_coefficient), se)) %>% 
  ungroup() %>% 
  mutate(lci=beta-1.96*se,
         uci=beta+1.96*se,
         final=paste0(round(beta, digits=2),
                      " [",
                      round(lci, digits=2),
                      ";",
                      round(uci, digits=2),
                      "]")) 
#univ all ages
#temp<-results_univ %>% filter(sex=="M")
age_univ<-results_univ %>% 
  ungroup() %>% group_by(sex) %>% 
  group_map(~{
    temp<-.x 
    temp<-temp %>% 
      ungroup() %>% 
      arrange(factor(var, levels=c(vars_cont)), desc(var)) %>% 
      mutate(id1=as.numeric(factor(var, levels=c(vars_cont)), desc(var)),
             id2=as.numeric(as.factor(age)),
             id=id1*6+id2-6)
    
    labels<-temp %>% filter(age==20) %>% select(id, age, var) %>% 
      mutate(id=id+0.5) %>% 
      arrange(factor(var, levels=c(vars_cont)), desc(var)) %>% 
      mutate(label=c("Population", "Growth",
                     "Pop. Density", "Fragmentation",
                     "Connectivity",
                     "Education", "Water Access", "Sewage",
                     "Overcrowding", "Social Index"))
    temp<-temp %>% mutate(label="") %>% bind_rows(labels)
    title<-ifelse(.y$sex=="M", "Men", "Women")
    ggplot(temp, aes(x=id, y=beta, group=id)) +
      geom_hline(yintercept = 0, lty=2, color="black")+
      geom_errorbar(aes(ymin=lci, ymax=uci, color=as.factor(age))) +
      geom_point(aes(fill=as.factor(age)), size=2, pch=21, color="black")+
      scale_x_continuous(breaks=temp$id, labels=temp$label)+
      ylim(c(-1.2, +1.2))+
      scale_color_discrete(name="Life Expectancy at",
                           labels=c("Birth", paste0("Age ", c(20, 40, 60))))+
      scale_fill_discrete(name="Life Expectancy at",
                          labels=c("Birth", paste0("Age ", c(20, 40, 60))))+
      xlab("") + ylab("Change in Life Expectancy (95% CI)")+
      theme_classic() +
      ggtitle(title)+
      theme(legend.position = "bottom",
            axis.text.y=element_text(face="bold", size=14, color="black"),
            axis.title=element_text(face="bold", size=16, color="black"),
            axis.ticks.x=element_blank(),
            plot.title=element_text(face="bold", size=24, color="black"),
            axis.text.x=element_text(face="bold", size=14, angle=90, hjust=1, color="black"),
            legend.text=element_text(size=14, color="black"),
            legend.title=element_text(face="bold", size=16, color="black"))
  })
legend<-get_legend(age_univ[[1]])
age_univ<-lapply(age_univ, function(xx) xx+guides(color=F, fill=F))
age_univ<-arrangeGrob(grobs=list(age_univ[[2]], age_univ[[1]]), ncol=2,
                      top=textGrob("Univariable", gp=gpar(fontsize=30, face="bold")))
#multiv all ages
#temp<-results_multiv %>% filter(sex=="M")
age_multiv<-results_multiv %>% 
  ungroup() %>% group_by(sex) %>% 
  group_map(~{
    temp<-.x 
    temp<-temp %>% 
      ungroup() %>% 
      arrange(factor(var, levels=c(vars_multiv)), desc(var)) %>% 
      mutate(id1=as.numeric(factor(var, levels=c(vars_multiv)), desc(var)),
             id2=as.numeric(as.factor(age)),
             id=id1*6+id2-6)
    
    labels<-temp %>% filter(age==20) %>% select(id, age, var) %>% 
      mutate(id=id+0.5) %>% 
      arrange(factor(var, levels=c(vars_multiv)), desc(var)) %>% 
      mutate(label=c("Population", "Growth",
                     "Pop. Density", "Fragmentation",
                     "Connectivity","Social Index"))
    temp<-temp %>% mutate(label="") %>% bind_rows(labels)
    title<-ifelse(.y$sex=="M", "Men", "Women")
    ggplot(temp, aes(x=id, y=beta, group=id)) +
      geom_hline(yintercept = 0, lty=2, color="black")+
      geom_errorbar(aes(ymin=lci, ymax=uci, color=as.factor(age))) +
      geom_point(aes(fill=as.factor(age)), size=2, pch=21, color="black")+
      scale_x_continuous(breaks=temp$id, labels=temp$label)+
      ylim(c(-1.2, +1.2))+
      scale_color_discrete(name="Life Expectancy at",
                           labels=c("Birth", paste0("Age ", c(20, 40, 60))))+
      scale_fill_discrete(name="Life Expectancy at",
                          labels=c("Birth", paste0("Age ", c(20, 40, 60))))+
      xlab("") + ylab("Change in Life Expectancy (95% CI)")+
      theme_classic() +
      ggtitle(title)+
      theme(legend.position = "bottom",
            axis.text.y=element_text(face="bold", size=14, color="black"),
            axis.title=element_text(face="bold", size=16, color="black"),
            axis.ticks.x=element_blank(),
            plot.title=element_text(face="bold", size=24, color="black"),
            axis.text.x=element_text(face="bold", size=16, angle=90, hjust=1, color="black"),
            legend.text=element_text(size=16, color="black"),
            legend.title=element_text(face="bold", size=18, color="black"))
  })
age_multiv<-lapply(age_multiv, function(xx) xx+guides(color=F, fill=F))
age_multiv<-arrangeGrob(grobs=list(age_multiv[[2]], age_multiv[[1]]), ncol=2,
                        top=textGrob("Multivariable", gp=gpar(fontsize=30, face="bold")))

age_both<-arrangeGrob(grobs=list(age_univ, age_multiv), ncol=1)
age_both<-arrangeGrob(grobs=list(age_both, legend), nrows=2,
                      heights=c(20, 1))
ggsave("results/ExtendedData3.pdf", age_both, width=20, height=20/1.75)
ggsave("results/ExtendedData3.eps", age_both, width=20, height=20/1.75)
# source data
results_multiv %>% full_join(data.frame(var=vars_multiv, label=c("Population", "Growth",
                                              "Pop. Density", "Fragmentation",
                                              "Connectivity","Social Index"))) %>% 
  select(label, age, sex, beta, lci, uci) %>% 
  fwrite("Source/Source_Data_ExtendedData3.csv")


# MORTALITY PART OF PAPER ANALYSIS STARTS HERE
mortality_cause<-mortality_cause %>% 
  left_join(l1s) %>% 
  mutate(city_link=ifelse(grepl("Tucuman", city_link), "Tucuman", city_link))

# ranges, overall PM, etc. for reuslts/discussion
# calculate overall PM
mortality_cause %>% ungroup %>% summarise_if(is.numeric, sum) %>% 
  mutate( p_cmnn=cmnn/total,
          p_cancer=cancer/total,
          p_ncd=ncd/total,
          p_accident=accident/total,
          p_violent=violent/total) %>% 
  select(p_cmnn, p_cancer, p_ncd, p_accident, p_violent)
# ranges of PM
mortality_cause %>% 
  mutate(p_ncds=p_ncd+p_cancer,
         p_injuries=p_accident+p_violent) %>% 
  select(SALID1, p_cmnn, p_ncds, p_injuries, p_cancer, p_ncd, p_accident, p_violent) %>% 
  gather(type, prop, -SALID1) %>% 
  group_by(type) %>% 
  summarise(min=min(prop*100),
            max=max(prop*100))


# ICCs for appendix
# get iccs for descriptives
iccs<-map_dfr(cause_titles_coll, function(name){
  f<-as.formula(paste0("p_", name, " ~ 1 | iso2"))
  m1<-lmer(f, data=mortality_cause )
  data.frame(icc=ICC(m1), cause=name, stringsAsFactors = F)
}) %>%  
  mutate(icc2=1-icc) %>% 
  select(cause, icc, icc2) %>% 
  mutate(icc=paste0(format(icc*100, digits=0), "%"),
         icc2=paste0(format(icc2*100, digits=0), "%"))
# repeating icc with age-adjusted 
colnames(aapm)<-sub("_ageadj", "", colnames(aapm))

aapm<-aapm %>% left_join(l1s)
iccs_aapm<-map_dfr(cause_titles_coll, function(name){
  f<-as.formula(paste0("p_", name, " ~ 1 | iso2"))
  m1<-lmer(f, data=aapm )
  data.frame(icc_aapm=ICC(m1), cause=name, stringsAsFactors = F)
}) %>%  
  mutate(icc2_aapm=1-icc_aapm) %>% 
  select(cause, icc_aapm, icc2_aapm) %>% 
  mutate(icc_aapm=paste0(format(icc_aapm*100, digits=0), "%"),
         icc2_aapm=paste0(format(icc2_aapm*100, digits=0), "%"))
iccs<-full_join(iccs, iccs_aapm)
iccs
fwrite(iccs, file="results/SupplementTable54.csv")


# table of ranges and overall by country
t1<-mortality_cause %>% mutate(n=1) %>% 
  group_by(iso2) %>% 
  summarise_if(is.numeric, sum) %>% 
  mutate( p_cmnn=cmnn/total,
          p_cancer=cancer/total,
          p_ncd=ncd/total,
          p_accident=accident/total,
          p_violent=violent/total,
          n=n) %>% 
  select(iso2, p_cmnn, p_cancer, p_ncd, p_accident, p_violent) %>% 
  gather(type, prop, -iso2) %>% 
  mutate(prop=prop*100)
t2<-mortality_cause %>% 
  select(iso2, SALID1, p_cmnn, p_cancer, p_ncd, p_accident, p_violent) %>% 
  gather(type, prop, -iso2, -SALID1) %>% 
  group_by(iso2, type) %>% 
  summarise(min=min(prop*100),
            max=max(prop*100))
# for all countries also
t1all<-mortality_cause %>% 
  mutate(n=1, 
         iso2="ZZALL") %>% 
  group_by(iso2) %>% 
  summarise_if(is.numeric, sum) %>% 
  mutate( p_cmnn=cmnn/total,
          p_cancer=cancer/total,
          p_ncd=ncd/total,
          p_accident=accident/total,
          p_violent=violent/total,
          n=n) %>% 
  select(iso2, p_cmnn, p_cancer, p_ncd, p_accident, p_violent) %>% 
  gather(type, prop, -iso2) %>% 
  mutate(prop=prop*100)
t2all<-mortality_cause %>%
  mutate(iso2="ZZALL") %>% 
  select(iso2, SALID1, p_cmnn, p_cancer, p_ncd, p_accident, p_violent) %>% 
  gather(type, prop, -iso2, -SALID1) %>% 
  group_by(iso2, type) %>% 
  summarise(min=min(prop*100),
            max=max(prop*100))
# n cities
ncities<-bind_rows(mortality_cause %>% 
                     group_by(iso2) %>% 
                     summarise(n=n()),
                   mortality_cause %>% 
                     mutate(iso2="ZZALL") %>% 
                     group_by(iso2) %>% 
                     summarise(n=n()))
# put together
minmax<-full_join(t1, t2) %>% 
  bind_rows(full_join(t1all, t2all)) %>% 
  mutate(value=ifelse(iso2!="CR",paste0(paste0(round(prop, digits=1), "%"), 
                                        " [",round(min, digits=1),
                                        ";",round(max, digits=1),
                                        "]"),
                      paste0(paste0(round(prop, digits=1), "%")))) %>%
  select(iso2, type, value) %>% 
  spread(type, value) %>% 
  full_join(ncities) %>% 
  select(iso2, n, p_cmnn, p_cancer, p_ncd, p_accident, p_violent) 
## add ill-defined diseases and injuries
mortality_illdefined<-mortality_illdefined %>% left_join(l1s)
t1<-mortality_illdefined %>% mutate(n=1) %>% 
  group_by(iso2) %>% 
  summarise_if(is.numeric, sum) %>% 
  mutate( p_illdefined_disease=illdefined_disease/total,
          p_illdefined_injury=p_illdefined_injury/total) %>% 
  select(iso2, p_illdefined_disease,p_illdefined_injury ) %>% 
  gather(type, prop, -iso2) %>% 
  mutate(prop=prop*100)
t2<-mortality_illdefined %>% 
  select(iso2, SALID1, p_illdefined_disease, p_illdefined_injury) %>% 
  gather(type, prop, -iso2, -SALID1) %>% 
  group_by(iso2, type) %>% 
  summarise(min=min(prop*100),
            max=max(prop*100))
# for all countries also
t1all<-mortality_illdefined %>% 
  mutate(iso2="ZZALL") %>% 
  group_by(iso2) %>% 
  summarise_if(is.numeric, sum) %>% 
  mutate( p_illdefined_disease=illdefined_disease/total,
          p_illdefined_injury=p_illdefined_injury/total) %>% 
  select(iso2, p_illdefined_disease,p_illdefined_injury ) %>% 
  gather(type, prop, -iso2) %>% 
  mutate(prop=prop*100)
t2all<-mortality_illdefined %>% 
  mutate(iso2="ZZALL") %>% 
  select(iso2, SALID1, p_illdefined_disease, p_illdefined_injury) %>% 
  gather(type, prop, -iso2, -SALID1) %>% 
  group_by(iso2, type) %>% 
  summarise(min=min(prop*100),
            max=max(prop*100))
# put together
minmax_ill<-full_join(t1, t2) %>% 
  bind_rows(full_join(t1all, t2all))
t1<-minmax_ill %>% 
  filter(grepl("disease", type)) %>% 
  mutate(value=ifelse(iso2!="CR",paste0(paste0(round(prop, digits=1), "%"), 
                                        " [",round(min, digits=1),
                                        ";",round(max, digits=1),
                                        "]"),
                      paste0(paste0(round(prop, digits=1), "%")))) %>%
  select(iso2, type, value) %>% 
  spread(type, value) %>% 
  select(iso2, p_illdefined_disease) 
t2<-minmax_ill %>% 
  filter(grepl("injury", type)) %>% 
  mutate(value=ifelse(iso2!="CR",paste0(paste0(round(prop, digits=2), "%"), 
                                        " [",round(min, digits=2),
                                        ";",round(max, digits=2),
                                        "]"),
                      paste0(paste0(round(prop, digits=2), "%")))) %>%
  select(iso2, type, value) %>% 
  spread(type, value) %>% 
  select(iso2, p_illdefined_injury) 

minmax<-full_join(minmax, t1) %>% full_join(t2)
minmax
fwrite(minmax, file="results/SupplementTable53.csv")


# Ill-defined deaths
mortality_illdefined<-mortality_illdefined %>% left_join(l1s)
#i<-1
cols<-gg_color_hue(9)
cause_titles_ill<-c("illdefined_disease", "illdefined_injury")
p<-map(seq_along(cause_titles_ill), function(i){
  titles<-c("Ill-defined Diseases","Injuries of Ill-defined Intent")
  cause<-cause_titles_ill[[i]]
  var<-paste0("p_", cause)
  temp<-mortality_illdefined
  temp$label<-temp$country_name
  #title<-paste0(titles[[i]], " (ICC=", iccs[[i]], ")")
  title<-titles[[i]]
  temp2<-temp
  temp2$iso2<-"ZZ"
  temp2$label<-"All"
  temp<-rbind(temp, temp2)
  colnames(temp)[colnames(temp)==var]<-"outcome"
  ggplot(temp, ggplot2::aes(x=iso2, y=outcome, group=label)) +
    #geom_boxplot(aes(group=as.factor(label)), fill=NA, outlier.color = NA, width=0.5)+
    geom_jitter(aes(fill=as.factor(iso2)),
                color="black", pch=21,
                width=0.1, height=0) +
    guides(fill=F, size=F, alpha=F)+
    xlab("") + ylab("Proportionate Mortality") +
    scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,.35),
                       breaks=seq(0, 0.7, by=0.1))+
    scale_x_discrete(labels=c(sort(unique(mortality_illdefined$iso2)), "All"))+
    scale_fill_manual(values=c(cols, "gray"))+
    scale_alpha_manual(values=c(0.5, 0.2))+
    ggtitle(title)+
    ggplot2::theme_bw() +
    theme(legend.position = "bottom",
          legend.key.width = unit(50, "points"),
          axis.text.x=element_text(size=13, face="bold",color="black"),
          axis.text.y=element_text(size=13, face="bold",color="black"),
          axis.title.y=element_text(size=16, face="bold",color="black"),
          plot.title=element_text(face="bold", size=16))
})
pmortality<-arrangeGrob(grobs=p, ncol=2)
plot(pmortality)
ggsave("results/SupplementaryFigure56.pdf", pmortality, width=11*2/2, height=7.5*1/2)


# FIGURE3
colors<-c("#55ACEE", "#FF5757", "#F4BC33", "#8E13A2", "#03989E")
# get ordering based on violent deaths
country_order<-mortality_cause %>% group_by(iso2) %>% 
  summarise(violent=sum(violent),
            total=sum(total),
            p_violent=violent/total) %>% 
  arrange(p_violent) %>% pull(iso2) %>% 
  data.frame(country_id=1:length(unique(mortality_cause$iso2)), stringsAsFactors = F) %>% 
  rename(iso2='.')
mortality_cause<-full_join(mortality_cause, country_order)
figure_data<-mortality_cause %>% 
  ungroup() %>% 
  arrange(country_id, p_violent) %>% 
  mutate(new_id=row_number()) %>% 
  select(SALID1, city_link, new_id, country_id, iso2,
         p_cmnn, p_cancer, p_ncd, p_accident, p_violent) %>% 
  gather(id, prop, -SALID1,-city_link, -new_id, -country_id, -iso2) %>% 
  mutate(id=case_when(
    id=="p_cmnn"~"1cmnn",
    id=="p_cancer"~"2cancer",
    id=="p_ncd"~"3ncd",
    id=="p_accident"~"4accident",
    id=="p_violent"~"5violent"
  ),
  city_label=paste0(city_link, ", ", iso2))
# separating countries
limits<-figure_data %>% group_by(iso2) %>% 
  summarise(limit=max(new_id)+0.5) %>% 
  arrange(limit) %>% slice(-n()) %>% pull(limit)

# middle: where to place name of country
f2<-ggplot(figure_data, aes(x=as.factor(new_id), y=prop, fill=as.factor(id)))+
  geom_bar(width = 1, stat = "identity", color=NA, size=0)+
  geom_vline(xintercept = limits, lty=2, color="black")+
  scale_fill_manual(labels=c("Communicable/Maternal/Neonatal/Nutritional", 
                             "Cancer", "CVD and other NCD",
                             "Unintentional Injuries", "Violent Injuries"),
                    name="Cause",
                    values=colors)+
  scale_alpha_manual(values=c(0,1))+
  coord_cartesian(clip = "off", ylim=c(0, 1))+
  scale_x_discrete(labels=unique(figure_data$city_label))+
  scale_y_continuous(expand=c(0,0),
                     labels=percent,
                     limits=c(0, 1.03),
                     breaks=seq(0, 1, by=.2),
                     sec.axis = sec_axis(~.*1, name = "Proportionate Mortality",
                                         labels=percent,
                                         breaks=seq(0, 1, by=.2)))+
  labs(x="", y="Proportionate Mortality", title="")+
  guides(fill=guide_legend(nrow=2))+
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=26, color="black"),
        legend.title=element_text(size=30, color="black", face="bold"),
        plot.title=element_text(face="bold", size=20),
        #plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.x=element_text(face="bold", color="black"),
        axis.text.x=element_text(color="black", size=4,angle=90, hjust=1, vjust=.5),
        axis.title.y=element_text(face="bold", color="black", size=30),
        axis.text.y=element_text(color="black", size=26))
#annotate country
{
  middle<-figure_data %>% group_by(iso2) %>% 
    summarise(min=min(new_id),
              max=max(new_id)) %>% 
    rowwise() %>% 
    mutate(middle=round(mean(c(min, max)))) %>% 
    select(iso2, middle) %>% 
    arrange(middle) %>% left_join(l1s %>% filter(!duplicated(iso2))) %>% 
    ungroup()
  
  annotate_size<-5
  for (i in which(!middle %>% pull(iso2) %in%c("CR", "PA", "SV"))){
    f2<-f2+annotate("text", 
                    label=middle %>% slice(i) %>% pull(country_name), 
                    x=middle %>% slice(i) %>% pull(middle), 
                    y=.4,
                    fontface="bold", size=annotate_size)
  }
  f2<-f2+annotate("text", 
                  label="Costa Rica", 
                  x=middle %>% filter(iso2=="CR") %>% pull(middle) %>% `-`(10), 
                  y=.55,
                  fontface="bold", size=annotate_size)
  f2<-f2+annotate("segment", 
                  x=middle %>% filter(iso2=="CR") %>% pull(middle) %>% `-`(10), 
                  xend = middle %>% filter(iso2=="CR") %>% pull(middle),
                  y=.55,
                  yend=.4,arrow=arrow())
  f2<-f2+annotate("text", 
                  label="Panama", 
                  x=middle %>% filter(iso2=="PA") %>% pull(middle) %>% `+`(10), 
                  y=.5,
                  fontface="bold", size=annotate_size)
  f2<-f2+annotate("segment", 
                  x=middle %>% filter(iso2=="PA") %>% pull(middle) %>% `+`(10), 
                  xend = middle %>% filter(iso2=="PA") %>% pull(middle),
                  y=.5,
                  yend=.45,arrow=arrow())
  f2<-f2+annotate("text", 
                  label="El Salvador", 
                  x=middle %>% filter(iso2=="SV") %>% pull(middle) %>% `-`(10), 
                  y=.55,
                  fontface="bold", size=annotate_size)
  f2<-f2+annotate("segment", 
                  x=middle %>% filter(iso2=="SV") %>% pull(middle) %>% `-`(10), 
                  xend = middle %>% filter(iso2=="SV") %>% pull(middle),
                  y=.55,
                  yend=.4,arrow=arrow())
  
  
  }
#FINAL
f2+theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
ggsave("results/Figure3.pdf", f2+theme(axis.text.x=element_blank(), axis.ticks.x = element_blank()), 
       width=20, height=20/(10/6))

figure_data %>% select(city_link, iso2, id, prop) %>% spread(id, prop) %>% 
  rename(city=city_link, CMNN=`1cmnn`, Cancer=`2cancer`, NCDs=`3ncd`, Unintentional=`4accident`, Violence=`5violent`) %>% 
  fwrite("Source/Source_Data_Figure3.csv")
# same data for ED7 and ED5
figure_data %>% select(city_link, iso2, id, prop) %>% spread(id, prop) %>% 
  rename(city=city_link, CMNN=`1cmnn`, Cancer=`2cancer`, NCDs=`3ncd`, Unintentional=`4accident`, Violence=`5violent`) %>% 
  fwrite("Source/Source_Data_ExtendedData7.csv")
figure_data %>% select(city_link, iso2, id, prop) %>% spread(id, prop) %>% 
  rename(city=city_link, CMNN=`1cmnn`, Cancer=`2cancer`, NCDs=`3ncd`, Unintentional=`4accident`, Violence=`5violent`) %>% 
  fwrite("Source/Source_Data_ExtendedData5.csv")


#age-adjusted version of figure 3
aapm<-aapm %>% 
  mutate(city_link=ifelse(grepl("Tucuman", city_link), "Tucuman", city_link)) 
# get ordering based on violent deaths
country_order<-aapm %>% group_by(iso2) %>% 
  # imperfect
  summarise(p_violent=mean(p_violent)) %>% 
  arrange(p_violent) %>% pull(iso2) %>% 
  data.frame(country_id=1:length(unique(aapm$iso2)), stringsAsFactors = F) %>% 
  rename(iso2='.')
aapm<-full_join(aapm, country_order)
figure_data<-aapm %>% 
  ungroup() %>% 
  arrange(country_id, p_violent) %>% 
  mutate(new_id=row_number()) %>% 
  select(SALID1, city_link, new_id, country_id, iso2,
         p_cmnn, p_cancer, p_ncd, p_accident, p_violent) %>% 
  gather(id, prop, -SALID1,-city_link, -new_id, -country_id, -iso2) %>% 
  mutate(id=case_when(
    id=="p_cmnn"~"1cmnn",
    id=="p_cancer"~"2cancer",
    id=="p_ncd"~"3ncd",
    id=="p_accident"~"4accident",
    id=="p_violent"~"5violent"
  ),
  city_label=paste0(city_link, ", ", iso2))
# separating countries
limits<-figure_data %>% group_by(iso2) %>% 
  summarise(limit=max(new_id)+0.5) %>% 
  arrange(limit) %>% slice(-n()) %>% pull(limit)

# middle: where to place name of country
f2<-ggplot(figure_data, aes(x=as.factor(new_id), y=prop, fill=as.factor(id)))+
  geom_bar(width = 1, stat = "identity", color=NA, size=0)+
  geom_vline(xintercept = limits, lty=2, color="gray")+
  scale_fill_manual(labels=c("Communicable/Maternal/Neonatal/Nutritional", 
                             "Cancer", "CVD and other NCD",
                             "Unintentional Injuries", "Violent Injuries"),
                    name="Cause",
                    values=colors)+
  scale_alpha_manual(values=c(0,1))+
  coord_cartesian(clip = "off", ylim=c(0, 1))+
  scale_x_discrete(labels=unique(figure_data$city_label))+
  scale_y_continuous(expand=c(0,0),
                     labels=percent,
                     limits=c(0, 1.03),
                     breaks=seq(0, 1, by=.2),
                     sec.axis = sec_axis(~.*1, name = "Proportionate Mortality",
                                         labels=percent,
                                         breaks=seq(0, 1, by=.2)))+
  labs(x="", y="Proportionate Mortality", title="")+
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=16, color="black"),
        legend.title=element_text(size=20, color="black", face="bold"),
        plot.title=element_text(face="bold", size=20),
        #plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.x=element_text(face="bold", color="black"),
        axis.text.x=element_text(color="black", size=4,angle=90, hjust=1, vjust=.5),
        axis.title.y=element_text(face="bold", color="black"),
        axis.text.y=element_text(color="black"))
#annotate country
{
  middle<-figure_data %>% group_by(iso2) %>% 
    summarise(min=min(new_id),
              max=max(new_id)) %>% 
    rowwise() %>% 
    mutate(middle=round(mean(c(min, max)))) %>% 
    select(iso2, middle) %>% 
    arrange(middle) %>% left_join(l1s %>% filter(!duplicated(iso2)))
  
  annotate_size<-5
  for (i in which(!middle %>% pull(iso2) %in%c("CR", "PA", "SV"))){
    f2<-f2+annotate("text", 
                    label=middle %>% ungroup() %>% slice(i) %>% pull(country_name), 
                    x=middle %>%  ungroup() %>% slice(i) %>% pull(middle), 
                    y=.4,
                    fontface="bold", size=annotate_size)
  }
  f2<-f2+annotate("text", 
                  label="Costa Rica", 
                  x=middle %>% filter(iso2=="CR") %>% pull(middle) %>% `-`(10), 
                  y=.55,
                  fontface="bold", size=annotate_size)
  f2<-f2+annotate("segment", 
                  x=middle %>% filter(iso2=="CR") %>% pull(middle) %>% `-`(10), 
                  xend = middle %>% filter(iso2=="CR") %>% pull(middle),
                  y=.55,
                  yend=.4,arrow=arrow())
  f2<-f2+annotate("text", 
                  label="Panama", 
                  x=middle %>% filter(iso2=="PA") %>% pull(middle) %>% `+`(10), 
                  y=.5,
                  fontface="bold", size=annotate_size)
  f2<-f2+annotate("segment", 
                  x=middle %>% filter(iso2=="PA") %>% pull(middle) %>% `+`(10), 
                  xend = middle %>% filter(iso2=="PA") %>% pull(middle),
                  y=.5,
                  yend=.45,arrow=arrow())
  f2<-f2+annotate("text", 
                  label="El Salvador", 
                  x=middle %>% filter(iso2=="SV") %>% pull(middle) %>% `-`(10), 
                  y=.55,
                  fontface="bold", size=annotate_size)
  f2<-f2+annotate("segment", 
                  x=middle %>% filter(iso2=="SV") %>% pull(middle) %>% `-`(10), 
                  xend = middle %>% filter(iso2=="SV") %>% pull(middle),
                  y=.55,
                  yend=.4,arrow=arrow())
  
  
  }
f2<-f2+theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
f2
ggsave("results/ExtendedData6.pdf", f2, width=20, height=20/(10/6))
ggsave("results/ExtendedData6.eps", f2, width=20, height=20/(10/6))
figure_data %>% select(city_link, iso2, id, prop) %>% spread(id, prop) %>% 
  rename(city=city_link, CMNN=`1cmnn`, Cancer=`2cancer`, NCDs=`3ncd`, Unintentional=`4accident`, Violence=`5violent`) %>% 
  fwrite("Source/Source_Data_ExtendedData6.csv")


# FIGURE 4
# melt, run line of exposure in x axis and prop mortality on y, showing 5 lines for 5 outcomes
mortality<-mortality_cause %>% select(SALID1, p_cmnn, p_cancer, p_ncd, p_accident, p_violent) %>% 
  gather(type, prop, -SALID1) %>% 
  left_join(all_exposure) %>% 
  mutate(type=case_when(
    type=="p_cmnn"~"1cmnn",
    type=="p_cancer"~"2cancer",
    type=="p_ncd"~"3ncd",
    type=="p_accident"~"4accident",
    type=="p_violent"~"5violent"
  ))

color_labels<-c("Communicable/Maternal/Neonatal/Nutritional","Cancer", "CVD and other NCDs", "Unintentional Injuries", "Violent Injuries")
mortality<-mortality %>% 
  mutate(exp=sei)
lowess_est<-mortality %>% group_by(type) %>% 
  group_modify(~{
    t<-loess(prop~exp, data=.x)
    min<-min(.x$exp, na.rm=T)
    max<-max(.x$exp, na.rm=T)
    newdata<-data.frame(exp=seq(min, max, length.out = 1000))
    newdata$prop<-predict(t, newdata=newdata)
    newdata$type<-unique(.x$type)
    newdata
  })
xvar<-"Social Environment Index"
ytitle1<-ytitle2<-"Proportionate Mortality"
# get data to plot a "rug" at the bottom with the actual data
exposure<-mortality %>% filter(!duplicated(SALID1)) %>% select(SALID1, exp)
# include 5 estimates of %
buckets<-5
q<-quantile(exposure$exp, probs = seq(0, 1, by=1/(buckets-1)), na.rm=T)
# alternative: take range... and divide it over 5 
range<-c(min(exposure$exp, na.rm=T),max(exposure$exp, na.rm=T))
rangediff<-abs(diff(range))/(buckets-1)
q<-sapply(0:(buckets-1), function(xx) range[1]+rangediff*xx)
# move extremes inwards a tiny bit so you can see the full number
q[1]<-q[1]+(max(exposure$exp, na.rm=T)-min(exposure$exp, na.rm=T))*0.075
q[buckets]<-q[buckets]-(max(exposure$exp, na.rm=T)-min(exposure$exp, na.rm=T))*0.075
q<-data.frame(value=q, quant=1:buckets)
# find closest to each quantile
q<-map_dfr(q %>% pull(quant), function(quant_temp){
  xx<-q %>% filter(quant==quant_temp)
  t<-lowess_est
  t$dif<-t$exp-xx$value
  t<-t %>% arrange(abs(t$dif)) %>% filter(!duplicated(type)) %>% 
    ungroup() %>% 
    slice(1:5) %>% 
    select(type, prop, exp) %>% 
    mutate(id=6-row_number())
  t$y[t$id==1]<-t$prop[t$id==1]*0.5
  for (i in 2:5){
    t$y[t$id==i]<-sum(t$prop[t$id%in%c(1:(i-1))])+t$prop[t$id==i]*0.5
  }
  t
})
# figure ys out: half of the range
# so for the last one: half of itself, for next one, previous+half of itself
q$label<-paste0(format(q$prop*100, digits=1, nsmall=1), "%")
p1<-ggplot(lowess_est)+
  geom_area(aes(x=exp, y=prop, fill=as.factor(type)))+
  geom_linerange(data=exposure,
                 aes(x=exp, ymin=0, ymax=0.02), alpha=0.4, size=1.5)+
  geom_text(data=q, aes(x=exp, y=y, label=label), color="white", size=8)+
  scale_fill_manual(labels=color_labels,
                    name="Cause",
                    values=colors)+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_y_continuous(expand=c(0,0),
                     breaks=seq(0, 1, by=.25),
                     limits=c(0, 1.001),
                     sec.axis = sec_axis(~.*1, name = "",
                                         labels=percent,
                                         breaks=seq(0, 1, by=.25)),
                     labels=percent)+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Social Environment Index", y="Proportionate Mortality",
       title="Not age-adjusted")+
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title.x=element_text(face="bold", color="black", size=20),
        axis.title.y=element_text(face="bold", color="black", size=20),
        axis.text.y=element_text(color="black", size=16),
        axis.text.x=element_text(color="black", size=16),
        plot.title=element_text(face="bold", size=20),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black", face="bold"))
p1

# age-adjusted version
mortality<-aapm %>% select(SALID1, p_cmnn, p_cancer, p_ncd, p_accident, p_violent) %>% 
  gather(type, prop, -SALID1) %>% 
  left_join(all_exposure) %>% 
  mutate(type=case_when(
    type=="p_cmnn"~"1cmnn",
    type=="p_cancer"~"2cancer",
    type=="p_ncd"~"3ncd",
    type=="p_accident"~"4accident",
    type=="p_violent"~"5violent"
  ))

color_labels<-c("Communicable/Maternal/Neonatal/Nutritional","Cancer", "CVD and other NCDs", "Unintentional Injuries", "Violent Injuries")
mortality<-mortality %>% 
  mutate(exp=sei)
lowess_est<-mortality %>% group_by(type) %>% 
  group_modify(~{
    t<-loess(prop~exp, data=.x)
    min<-min(.x$exp, na.rm=T)
    max<-max(.x$exp, na.rm=T)
    newdata<-data.frame(exp=seq(min, max, length.out = 1000))
    newdata$prop<-predict(t, newdata=newdata)
    newdata$type<-unique(.x$type)
    newdata
  })
xvar<-"Social Environment Index"
ytitle1<-ytitle2<-"Proportionate Mortality"
# get data to plot a "rug" at the bottom with the actual data
exposure<-mortality %>% filter(!duplicated(SALID1)) %>% select(SALID1, exp)
# include 5 estimates of %
buckets<-5
q<-quantile(exposure$exp, probs = seq(0, 1, by=1/(buckets-1)), na.rm=T)
# alternative: take range... and divide it over 5 
range<-c(min(exposure$exp, na.rm=T),max(exposure$exp, na.rm=T))
rangediff<-abs(diff(range))/(buckets-1)
q<-sapply(0:(buckets-1), function(xx) range[1]+rangediff*xx)
# move extremes inwards a tiny bit so you can see the full number
q[1]<-q[1]+(max(exposure$exp, na.rm=T)-min(exposure$exp, na.rm=T))*0.075
q[buckets]<-q[buckets]-(max(exposure$exp, na.rm=T)-min(exposure$exp, na.rm=T))*0.075
q<-data.frame(value=q, quant=1:buckets)
# find closest to each quantile
q<-map_dfr(q %>% pull(quant), function(quant_temp){
  xx<-q %>% filter(quant==quant_temp)
  t<-lowess_est
  t$dif<-t$exp-xx$value
  t<-t %>% arrange(abs(t$dif)) %>% filter(!duplicated(type)) %>% 
    ungroup() %>% 
    slice(1:5) %>% 
    select(type, prop, exp) %>% 
    mutate(id=6-row_number())
  t$y[t$id==1]<-t$prop[t$id==1]*0.5
  for (i in 2:5){
    t$y[t$id==i]<-sum(t$prop[t$id%in%c(1:(i-1))])+t$prop[t$id==i]*0.5
  }
  t
})
# figure ys out: half of the range
# so for the last one: half of itself, for next one, previous+half of itself
q$label<-paste0(format(q$prop*100, digits=1, nsmall=1), "%")
p2<-ggplot(lowess_est)+
  geom_area(aes(x=exp, y=prop, fill=as.factor(type)))+
  geom_linerange(data=exposure,
                 aes(x=exp, ymin=0, ymax=0.02), alpha=0.4, size=1.5)+
  geom_text(data=q, aes(x=exp, y=y, label=label), color="white", size=8)+
  scale_fill_manual(labels=color_labels,
                    name="Cause",
                    values=colors)+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  scale_y_continuous(expand=c(0,0),
                     breaks=seq(0, 1, by=.25),
                     limits=c(0, 1.001),
                     sec.axis = sec_axis(~.*1, name = "Proportionate Mortality",
                                         labels=percent,
                                         breaks=seq(0, 1, by=.25)),
                     labels=percent)+
  scale_x_continuous(expand=c(0,0))+
  labs(x="Social Environment Index", y="",
       title="Age-adjusted")+
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title.x=element_text(face="bold", color="black", size=20),
        axis.title.y=element_text(face="bold", color="black", size=20),
        axis.text.y=element_text(color="black", size=16),
        axis.text.x=element_text(color="black", size=16),
        plot.title=element_text(face="bold", size=20),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black", face="bold"))
p2
legend<-get_legend(p1)
p1<-p1+guides(fill=F)
p2<-p2+guides(fill=F)
pall<-arrangeGrob(grobs=list(p1, p2), ncol=2)
pall<-arrangeGrob(grobs=list(pall, legend), ncol=1, heights=c(10, 1))
ggsave("results/Figure4.pdf",pall, width=20, height=7)

# all of the other variables
vars<-c("adj_rate","pop_baseline", "growth_pct", 
        "BECPOPDENSL1AD","BECPTCHDENSL1AD", "BECADINTDENSL1AD",
        "CNSMINPR_L1AD", "CNSWATINL1AD", "CNSSEWNETL1AD", "CNSCROWD3RML1AD")
xtitles<-c("Age-Adjusted All-Cause Mortality", "Population at Baseline (log scale)", "% Growth in Population over 5 years",
           "Population Density (Population/km2)","Patch Density (Patches/km2)","Intersection Density (intersections/km2)",
           "% with Completed primary Education or Above", "% HHs with Piped Water in the Household","% HHs with Connection to the Sewage Network","% Overcrowded HHs (>3 people per room)")
plottitles<-c("Mortality Levels","City Size", "Population Growth",
              "Population Density", "Fragmentation", "Street Connectivity",
              "Educational Attainment", "Water Access","Sewage Access", "Overcrowding")
iso2list<-c("AR", "BR", "CL", "CO", "CR", "MX", "PA", "PE", "SV")
# correct fragmentation and street connectivity and density for PCT URBAN (residuals)
all_exposure<-all_exposure %>% 
  ungroup() %>% 
  mutate(BECPTCHDENSL1AD_adj=all_exposure %>% 
           lm(formula=BECPTCHDENSL1AD~BECPCTURBANL1AD) %>% 
           augment %>% pull(.resid),
         BECPOPDENSL1AD_adj=all_exposure %>% 
           lm(formula=BECPOPDENSL1AD~BECPCTURBANL1AD) %>% 
           augment %>% pull(.resid))
all_exposure$BECADINTDENSL1AD_adj<-all_exposure %>% 
  lm(formula=BECADINTDENSL1AD~BECPCTURBANL1AD, na.action = na.exclude) %>% 
  augment(data=all_exposure) %>% pull(.resid)
# melt, run line of exposure in x axis and prop mortality on y, showing 5 lines for 5 outcomes
mortality<-mortality_cause %>% select(SALID1, p_cmnn, p_cancer, p_ncd, p_accident, p_violent) %>% 
  gather(type, prop, -SALID1) %>% 
  left_join(all_exposure) %>% 
  mutate(type=case_when(
    type=="p_cmnn"~"1cmnn",
    type=="p_cancer"~"2cancer",
    type=="p_ncd"~"3ncd",
    type=="p_accident"~"4accident",
    type=="p_violent"~"5violent"
  ))

color_labels<-c("Communicable/Maternal/Neonatal/Nutritional","Cancer", "CVD and other NCDs", "Unintentional Injuries", "Violent Injuries")
mortality$pop_baseline<-log10(mortality$pop_baseline)
mortality$growth_pct<-mortality$growth_pct*100
# lowess plots with equal intervals 
plots<-map(seq_along(vars), function(i){
  print(i)
  var<-vars[[i]]
  temp<-mortality %>% 
    rename(exp:=!!var)
  loess_mortality<-function(x, y){
    # x<-temp %>% filter(type=="p_cmnn")
    t<-loess(prop~exp, data=x)
    min<-min(x$exp, na.rm=T)
    max<-max(x$exp, na.rm=T)
    newdata<-data.frame(exp=seq(min, max, length.out = 1000))
    newdata$prop<-predict(t, newdata=newdata)
    newdata$type<-unique(x$type)
    newdata
  }
  lowess_est<-temp %>% group_by(type) %>% 
    group_modify(loess_mortality)
  
  xvar<-xtitles[i]
  ytitle1<-ifelse(i%in%c(1, 6), "Proportionate Mortality", "")
  ytitle2<-ifelse(i%in%c(5, 10), "Proportionate Mortality", "")
  exposure<-temp[!duplicated(temp$SALID1),c("SALID1", "exp")]
  # nudge min and max a little bit to the center so they show up (just a tiny bit)
  exposure<-exposure[order(exposure$exp),]
  exposure[1,"exp"]<-exposure[1,"exp"]+(max(exposure$exp, na.rm=T)-min(exposure$exp, na.rm=T))*0.01
  exposure<-exposure[order(exposure$exp, decreasing = T),]
  exposure[1,"exp"]<-exposure[1,"exp"]-(max(exposure$exp, na.rm=T)-min(exposure$exp, na.rm=T))*0.01
  # include around 5 estimates of %
  buckets<-5
  q<-quantile(exposure$exp, probs = seq(0, 1, by=1/(buckets-1)), na.rm=T)
  # alternative: take range... and divide it over 5 
  range<-c(min(exposure$exp, na.rm=T),max(exposure$exp, na.rm=T))
  rangediff<-abs(diff(range))/(buckets-1)
  q<-sapply(0:(buckets-1), function(xx) range[1]+rangediff*xx)
  # move extremes inwards a tiny bit so you can see the full number
  q[1]<-q[1]+(max(exposure$exp, na.rm=T)-min(exposure$exp, na.rm=T))*0.075
  q[buckets]<-q[buckets]-(max(exposure$exp, na.rm=T)-min(exposure$exp, na.rm=T))*0.075
  q<-data.frame(value=q, quant=1:buckets)
  # find closest to each quantile
  q<-map_dfr(q %>% pull(quant), function(quant_temp){
    xx<-q %>% filter(quant==quant_temp)
    t<-lowess_est
    t$dif<-t$exp-xx$value
    t<-t[order(abs(t$dif)),]
    t<-t[1:5,c("type", "prop", "exp")]
    t$id<-5:1
    t$y[t$id==1]<-t$prop[t$id==1]*0.5
    for (i in 2:5){
      t$y[t$id==i]<-sum(t$prop[t$id%in%c(1:(i-1))])+t$prop[t$id==i]*0.5
    }
    t
  })
  # figure ys out: half of the range
  # so for the last one: half of itself, for next one, previous+half of itself
  q$label<-paste0(format(q$prop*100, digits=1, nsmall=1), "%")
  p<-ggplot(lowess_est)+
    geom_area(aes(x=exp, y=prop, fill=as.factor(type)))+
    geom_linerange(data=exposure,
                   aes(x=exp, ymin=0, ymax=0.02), alpha=0.4, size=1.5)+
    geom_text(data=q, aes(x=exp, y=y, label=label), color="white", size=8)+
    scale_fill_manual(labels=color_labels,
                      name="Cause",
                      values=colors)+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))+
    scale_y_continuous(expand=c(0,0),
                       breaks=seq(0, 1, by=.25),
                       sec.axis = sec_axis(~.*1, name = ytitle2,
                                           labels=percent,
                                           breaks=seq(0, 1, by=.25)),
                       labels=percent)+
    scale_x_continuous(expand=c(0,0))+
    labs(x=xvar, y=ytitle1, title=plottitles[[i]])+
    theme_classic() +
    theme(legend.position = "bottom",
          legend.text=element_text(size=30, color="black"),
          legend.title=element_text(size=40, color="black", face="bold"),
          plot.title=element_text(face="bold", size=40),
          axis.title.x=element_text(face="bold", color="black", size=30),
          axis.title.y=element_text(face="bold", color="black", size=40),
          axis.text.y=element_text(color="black", size=30),
          axis.text.x=element_text(color="black", size=30))
  if(var=="pop_baseline"){
    p<-p+scale_x_continuous(expand=c(0,0),
                            breaks=log10(c(2*10^5, 5*10^5, 10^6, 2*10^6, 5*10^6, 10^7)),
                            labels=c("200K", "500K", "1M", "2M", "5M", "10M"))
  } 
  p
})
legend<-get_legend(plots[[1]])
plots<-lapply(plots, function(p) p+guides(fill=F))
pall<-arrangeGrob(grobs=plots[1:10], ncol=5)
pall<-arrangeGrob(grobs=list(pall, legend), nrow=2, heights=c(15, 1))
ggsave("results/ExtendedData8.pdf",pall, width=25, height=12.5, 
       units="in", scale=2, limitsize = F)
ggsave("results/ExtendedData8.eps",pall, width=25, height=12.5, 
       units="in", scale=2, limitsize = F)

# age-adjusted version
mortality<-aapm %>% select(SALID1, p_cmnn, p_cancer, p_ncd, p_accident, p_violent) %>% 
  gather(type, prop, -SALID1) %>% 
  left_join(all_exposure) %>% 
  mutate(type=case_when(
    type=="p_cmnn"~"1cmnn",
    type=="p_cancer"~"2cancer",
    type=="p_ncd"~"3ncd",
    type=="p_accident"~"4accident",
    type=="p_violent"~"5violent"
  ))

color_labels<-c("Communicable/Maternal/Neonatal/Nutritional","Cancer", "CVD and other NCDs", "Unintentional Injuries", "Violent Injuries")
mortality$pop_baseline<-log10(mortality$pop_baseline)
mortality$growth_pct<-mortality$growth_pct*100
# lowess plots with equal intervals 
plots<-map(seq_along(vars), function(i){
  print(i)
  var<-vars[[i]]
  temp<-mortality %>% 
    rename(exp:=!!var)
  loess_mortality<-function(x, y){
    # x<-temp %>% filter(type=="p_cmnn")
    t<-loess(prop~exp, data=x)
    min<-min(x$exp, na.rm=T)
    max<-max(x$exp, na.rm=T)
    newdata<-data.frame(exp=seq(min, max, length.out = 1000))
    newdata$prop<-predict(t, newdata=newdata)
    newdata$type<-unique(x$type)
    newdata
  }
  lowess_est<-temp %>% group_by(type) %>% 
    group_modify(loess_mortality)
  
  xvar<-xtitles[i]
  ytitle1<-ifelse(i%in%c(1, 6), "Proportionate Mortality", "")
  ytitle2<-ifelse(i%in%c(5, 10), "Proportionate Mortality", "")
  exposure<-temp[!duplicated(temp$SALID1),c("SALID1", "exp")]
  # nudge min and max a little bit to the center so they show up (just a tiny bit)
  exposure<-exposure[order(exposure$exp),]
  exposure[1,"exp"]<-exposure[1,"exp"]+(max(exposure$exp, na.rm=T)-min(exposure$exp, na.rm=T))*0.01
  exposure<-exposure[order(exposure$exp, decreasing = T),]
  exposure[1,"exp"]<-exposure[1,"exp"]-(max(exposure$exp, na.rm=T)-min(exposure$exp, na.rm=T))*0.01
  # include around 5 estimates of %
  buckets<-5
  q<-quantile(exposure$exp, probs = seq(0, 1, by=1/(buckets-1)), na.rm=T)
  # alternative: take range... and divide it over 5 
  range<-c(min(exposure$exp, na.rm=T),max(exposure$exp, na.rm=T))
  rangediff<-abs(diff(range))/(buckets-1)
  q<-sapply(0:(buckets-1), function(xx) range[1]+rangediff*xx)
  # move extremes inwards a tiny bit so you can see the full number
  q[1]<-q[1]+(max(exposure$exp, na.rm=T)-min(exposure$exp, na.rm=T))*0.075
  q[buckets]<-q[buckets]-(max(exposure$exp, na.rm=T)-min(exposure$exp, na.rm=T))*0.075
  q<-data.frame(value=q, quant=1:buckets)
  # find closest to each quantile
  q<-map_dfr(q %>% pull(quant), function(quant_temp){
    xx<-q %>% filter(quant==quant_temp)
    t<-lowess_est
    t$dif<-t$exp-xx$value
    t<-t[order(abs(t$dif)),]
    t<-t[1:5,c("type", "prop", "exp")]
    t$id<-5:1
    t$y[t$id==1]<-t$prop[t$id==1]*0.5
    for (i in 2:5){
      t$y[t$id==i]<-sum(t$prop[t$id%in%c(1:(i-1))])+t$prop[t$id==i]*0.5
    }
    t
  })
  # figure ys out: half of the range
  # so for the last one: half of itself, for next one, previous+half of itself
  q$label<-paste0(format(q$prop*100, digits=1, nsmall=1), "%")
  p<-ggplot(lowess_est)+
    geom_area(aes(x=exp, y=prop, fill=as.factor(type)))+
    geom_linerange(data=exposure,
                   aes(x=exp, ymin=0, ymax=0.02), alpha=0.4, size=1.5)+
    geom_text(data=q, aes(x=exp, y=y, label=label), color="white", size=8)+
    scale_fill_manual(labels=color_labels,
                      name="Cause",
                      values=colors)+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))+
    scale_y_continuous(expand=c(0,0),
                       breaks=seq(0, 1, by=.25),
                       sec.axis = sec_axis(~.*1, name = ytitle2,
                                           labels=percent,
                                           breaks=seq(0, 1, by=.25)),
                       labels=percent)+
    scale_x_continuous(expand=c(0,0))+
    labs(x=xvar, y=ytitle1, title=plottitles[[i]])+
    theme_classic() +
    theme(legend.position = "bottom",
          legend.text=element_text(size=30, color="black"),
          legend.title=element_text(size=40, color="black", face="bold"),
          plot.title=element_text(face="bold", size=40),
          axis.title.x=element_text(face="bold", color="black", size=30),
          axis.title.y=element_text(face="bold", color="black", size=40),
          axis.text.y=element_text(color="black", size=30),
          axis.text.x=element_text(color="black", size=30))
  if(var=="pop_baseline"){
    p<-p+scale_x_continuous(expand=c(0,0),
                            breaks=log10(c(2*10^5, 5*10^5, 10^6, 2*10^6, 5*10^6, 10^7)),
                            labels=c("200K", "500K", "1M", "2M", "5M", "10M"))
  } 
  p
})
legend<-get_legend(plots[[1]])
plots<-lapply(plots, function(p) p+guides(fill=F))
pall<-arrangeGrob(grobs=plots[1:10], ncol=5)
pall<-arrangeGrob(grobs=list(pall, legend), nrow=2, heights=c(15, 1))
ggsave("results/ExtendedData9.pdf",pall, width=25, height=12.5, 
       units="in", scale=2, limitsize = F)
ggsave("results/ExtendedData9.eps",pall, width=25, height=12.5, 
       units="in", scale=2, limitsize = F)

# Regression models
vars<-c("adj_rate","pop_baseline", "growth_pct", 
        "BECPOPDENSL1AD","BECPTCHDENSL1AD", "BECADINTDENSL1AD",
        "CNSMINPR_L1AD", "CNSWATINL1AD", "CNSSEWNETL1AD", "CNSCROWD3RML1AD", 
        "sei")

source("MS10_SALURBAL_Helper.R")
# prepare data for NB models
mortality<-mortality_cause_models %>% 
  select(SALID1, cmnn, cancer, ncd, accident, violent, pop) %>% 
  gather(type, deaths, -SALID1, -pop) %>%  
  mutate(lnpop=log(pop)) %>% 
  left_join(all_exposure) %>% 
  mutate(pop_baseline=log(pop_baseline)) %>% 
  mutate_at(vars[-2],
            scale_new) %>% 
  mutate(pop_baseline=pop_baseline-mean(pop_baseline),
         BECPCTURBANL1AD=scale_new(BECPCTURBANL1AD),
         under15=scale_new(under15),
         above65=scale_new(above65)) %>% 
  mutate(type=case_when(
    type=="cmnn"~"1cmnn",
    type=="cancer"~"2cancer",
    type=="ncd"~"0ncd",
    type=="accident"~"4accident",
    type=="violent"~"5violent"
  )) %>% 
  mutate(growth_pct_prev=scale_new(growth_pct_prev))
mortality %>% select(all_of(vars), BECPCTURBANL1AD) %>% summary
# THIS RUNS ALL mortality MODELS
if (run_models==T){
  # single variable model
  plan(multiprocess, .cleanup = T)
  models_univ <- future_map(vars, function(var){
    print(var)
    if (grepl("BEC", var)){
      f<-as.formula(paste0("deaths~offset(lnpop)+type*(", var, "+BECPCTURBANL1AD)+(1|iso2/SALID1)"))
    } else {
      f<-as.formula(paste0("deaths~offset(lnpop)+type*(", var, ")+(1|iso2/SALID1)"))
    }
    m<-glmer.nb(f, data=mortality, verbose=T,
                control=glmerControl(optimizer = "bobyqa",
                                     optCtrl=list(maxfun=20000)))
    m
  }, .progress = T)
  # add age as an interaction
  plan(multiprocess, .cleanup = T)
  models_univ_ageadj <- future_map(vars, function(var){
    print(var)
    if (grepl("BEC", var)){
      f<-as.formula(paste0("deaths~offset(lnpop)+type*(under15+above65+", var, "+BECPCTURBANL1AD)+(1|iso2/SALID1)"))
    } else {
      f<-as.formula(paste0("deaths~offset(lnpop)+type*(under15+above65+", var, ")+(1|iso2/SALID1)"))
    }
    m<-glmer.nb(f, data=mortality, verbose=T,
                control=glmerControl(optimizer = "bobyqa",
                                     optCtrl=list(maxfun=20000)))
    m
  }, .progress = T)
  # add age and mortality as an interaction
  plan(multiprocess)
  models_univ_mortality_ageadj <- future_map(vars[-1], function(var){
    print(var)
    if (grepl("BEC", var)){
      f<-as.formula(paste0("deaths~offset(lnpop)+type*(under15+above65+adj_rate+", var, "+BECPCTURBANL1AD)+(1|iso2/SALID1)"))
    } else {
      f<-as.formula(paste0("deaths~offset(lnpop)+type*(under15+above65+adj_rate+", var, ")+(1|iso2/SALID1)"))
    }
    m<-glmer.nb(f, data=mortality, verbose=T,
                control=glmerControl(optimizer = "bobyqa",
                                     optCtrl=list(maxfun=20000)))
    m
  }, .progress = T)
  # multiv models
  plan(multiprocess, .cleanup = T)
  models_multiv<-future_map(1:4, function(i){
    multiv<-paste(vars[!grepl("CNS", vars)],collapse="+")
    if (i==1){
      # adjusted by % under age 15 and above 65
      f<-as.formula(paste0("deaths~+offset(lnpop)+type*(under15+above65+", multiv, "+BECPCTURBANL1AD)+(1|iso2/SALID1)"))
      m<-glmer.nb(f, data=mortality, verbose=T,
                                    control=glmerControl(optimizer = "bobyqa",
                                                         optCtrl=list(maxfun=20000)))
    } else if (i==2){
      # model excluding cities above pct 90
      limit_pct90<-quantile(mortality_illdefined$p_illdefined_total, probs=.9)
      include<-mortality_illdefined %>% filter(p_illdefined_total<limit_pct90) %>% pull(SALID1)
      f<-as.formula(paste0("deaths~+offset(lnpop)+type*(under15+above65+", multiv, "+BECPCTURBANL1AD)+(1|iso2/SALID1)"))
      m<-glmer.nb(f, data=mortality %>% filter(SALID1%in%include), verbose=T,
                                               control=glmerControl(optimizer = "bobyqa",
                                                                    optCtrl=list(maxfun=20000)))
    } else if (i==3){
      # model using pop growth_prev (comparison model)
      # figure out restricted set of cities
      include<-all_exposure %>% filter(!is.na(growth_pct_prev)) %>% pull(SALID1)
      f<-as.formula(paste0("deaths~+offset(lnpop)+type*(under15+above65+", 
                           sub("growth_pct", "growth_pct", multiv), "+BECPCTURBANL1AD)+(1|iso2/SALID1)"))
      m<-glmer.nb(f, data=mortality %>% filter(SALID1%in%include), verbose=T,
                             control=glmerControl(optimizer = "bobyqa",
                                                  optCtrl=list(maxfun=20000)))
    } else if (i==4){
      # model using pop growth_prev (prev model)
      # figure out restricted set of cities
      include<-all_exposure %>% filter(!is.na(growth_pct_prev)) %>% pull(SALID1)
      f<-as.formula(paste0("deaths~+offset(lnpop)+type*(under15+above65+", 
                           sub("growth_pct", "growth_pct_prev", multiv), "+BECPCTURBANL1AD)+(1|iso2/SALID1)"))
      m<-glmer.nb(f, data=mortality %>% filter(SALID1%in%include), verbose=T,
                                  control=glmerControl(optimizer = "bobyqa",
                                                       optCtrl=list(maxfun=20000)))
    }
    m
  }, .progress = T)
  model_multiv_ageadj<-models_multiv[[1]]
  model_multiv_ageadj_restricted<-models_multiv[[2]]
  models_growth_mortality<-list(models_multiv[[3]], models_multiv[[4]])
  # save
  save(models_univ, models_univ_ageadj,
       models_univ_mortality_ageadj,
       model_multiv_ageadj,
       model_multiv_ageadj_restricted,
       models_growth_mortality,
       file="results/model_results_mortality.rdata")
} else {
  load("results/model_results_mortality.rdata")  
}



# extract coefficients for age-adjusted multivariable model
# population will be re-scaled by log(1+increase)
# in this case: 50% (log(1.5))
pop_coefficient<-1.5
coefs_multiv<-as.data.frame(summary(model_multiv_ageadj)$coefficients)
coefs_multiv$var<-rownames(coefs_multiv)
n<-model_multiv_ageadj@frame %>% pull(SALID1) %>% unique %>% length
coefs_multiv<-coefs_multiv[grepl("\\:", coefs_multiv$var),]
coefs_multiv<-coefs_multiv[!grepl("under15|above65|BECPCT", coefs_multiv$var),]
coefs_multiv$exp<-substr(coefs_multiv$var, regexpr("\\:", coefs_multiv$var)+1, nchar(coefs_multiv$var))
coefs_multiv$var<-substr(coefs_multiv$var, 5, nchar(coefs_multiv$var))
coefs_multiv$var<-substr(coefs_multiv$var, 1, regexpr("\\:", coefs_multiv$var)-1)
colnames(coefs_multiv)[1:2]<-c("est", "se")
coefs_multiv<-coefs_multiv %>% 
  select(var, exp, est, se) %>% 
  # re-scale population
  mutate(beta=ifelse(var=="pop_baseline", est*log(pop_coefficient), est),
         se=ifelse(var=="pop_baseline", se*log(pop_coefficient), se))

sds<-all_exposure %>% 
  ungroup() %>% 
  mutate(growth_pct=growth_pct*100) %>% 
  select(vars) %>% summarise_all(sd, na.rm=T) %>% 
  gather(exp, sd) %>% 
  mutate(sd=round(sd, digits=1)) %>% 
  mutate(sd=replace(sd, exp=="pop_baseline", paste0((pop_coefficient-1)*100, "%"))) %>% 
  mutate(sd=replace(sd, exp=="sei", "1 SD")) 

table_multiv<-coefs_multiv %>% 
  group_by(exp) %>% 
  group_modify(t3_function) %>% 
  left_join(sds) %>% 
  select(exp, sd, cmnn, cancer, cvd, accident, violent) %>% 
  arrange(factor(exp, levels = c(vars)), desc(exp)) 
table_multiv
fwrite(table_multiv, "results/Table2.csv")

# plot all model results
models_univ_list<-list(models_univ, models_univ_ageadj, models_univ_mortality_ageadj)
exclude_vars<-c("BECPCT", "BECPCT|under15|above65", "BECPCT|under15|above65|adj_rate")
models_multiv_list<-list(model_multiv_ageadj, 
                         model_multiv_ageadj_restricted)
multiv_names<-c("4multiv_ageadj",
                "5multi_ageadj_restricted")
coefs_figure_univ<-map2_dfr(models_univ_list, exclude_vars, function(model, exclude){
  model_name<-case_when(
    exclude=="BECPCT" ~ "1univ",
    exclude=="BECPCT|under15|above65" ~ "2univ_ageadj",
    exclude=="BECPCT|under15|above65|adj_rate" ~ "3univ_mortality_ageadj",
  )
  temp<-map_dfr(model, function(model_temp){
    tidy(model_temp) %>% filter(grepl("\\:", term), effect=="fixed") %>% 
      mutate(model=model_name,
             var=substr(term, regexpr("\\:", term)+1, nchar(term)),
             outcome=substr(term, 1, regexpr("\\:", term)-1)) %>% 
      select(var, model, outcome, estimate, std.error) %>% 
      filter(!grepl(exclude, var))
  })
  temp
})

coefs_figure_multiv<-map2_dfr(models_multiv_list, multiv_names, function(model, model_name){
  tidy(model) %>% filter(grepl("\\:", term), effect=="fixed") %>% 
    mutate(model=model_name,
           var=substr(term, regexpr("\\:", term)+1, nchar(term)),
           outcome=substr(term, 1, regexpr("\\:", term)-1)) %>% 
    select(var, model, outcome, estimate, std.error) 
})

all<-bind_rows(coefs_figure_univ,coefs_figure_multiv)
# reference for CVD
all<-bind_rows(all, 
               expand.grid(var=unique(all$var), 
                           model=unique(all$model), 
                           outcome="type3cvd", 
                           estimate=0, std.error=NA,
                           stringsAsFactors = F))
all<-all %>% filter(!var%in%c("BECPCTURBANL1AD","under15", "above65")) %>% 
  arrange(model, var, outcome)

all<-all %>% 
  mutate(var=factor(var, levels = c(vars)),
         outcome=factor(outcome),
         model=factor(model)) %>% 
  arrange(var, outcome, model) %>% 
  mutate(id=(as.numeric(outcome)-1)*(length(unique(all$model))+1)+
           as.numeric(model)) %>% 
  mutate(beta=exp(estimate),
         lci=exp(estimate-1.96*std.error),
         uci=exp(estimate+1.96*std.error),
         label=case_when(
           outcome=="type1cmnn" ~ "CMNN",
           outcome=="type2cancer" ~ "Cancer",
           outcome=="type3cvd" ~ "CVD/NCDs\n(Ref.)",
           outcome=="type4accident" ~ "Unintentional\nInjuries",
           outcome=="type5violent" ~ "Violence\nInjuries",
         ),
         label=replace(label, model!="3univ_mortality_ageadj", ""))
vars_labels<-c("Mortality Levels","City Size", "Population Growth",
               "Population Density", "Fragmentation", "Street Connectivity",
               "Educational Attainment", "Water Access","Sewage Access", "Overcrowding",
               "Social Environment Index")

# plot: univariable, univariable age-adjusted, univariable-mortality and age, multivariable, sensitivity

names(vars_labels)<-vars
p<-ggplot(all, aes(x=id, y=beta)) +
  geom_hline(yintercept = 1, lty=2, color="black")+
  geom_linerange(aes(color=model,ymin=lci, ymax=uci)) +
  geom_point(aes(fill=model, shape=model), color="black", size=2) +
  scale_x_continuous(breaks=all %>% filter(var=="pop_baseline") %>% pull(id), 
                     labels=all %>% filter(var=="pop_baseline") %>% pull(label))+
  scale_y_continuous(trans = log_trans(), breaks = base_breaks(3))+
  scale_fill_brewer(name="Model", labels=c("Univariable\n", 
                                           "Univariable\n(adjusted for age)",
                                           "Univariable\n(adjusted for age and mortality)",
                                           "Multivariable\n(adjusted for age)", 
                                           "Multivariable\n(adjusted for age,\nrestricted by % ill-defined deaths)"),
                    palette = 2, type="qual")+
  scale_color_brewer(name="Model", labels=c("Univariable\n", 
                                            "Univariable\n(adjusted for age)",
                                            "Univariable\n(adjusted for age and mortality)",
                                            "Multivariable\n(adjusted for age)", 
                                            "Multivariable\n(adjusted for age,\nrestricted by % ill-defined deaths)"),
                     palette = 2, type="qual")+
  scale_shape_manual(name="Model", labels=c("Univariable\n", 
                                            "Univariable\n(adjusted for age)",
                                            "Univariable\n(adjusted for age and mortality)",
                                            "Multivariable\n(adjusted for age)", 
                                            "Multivariable\n(adjusted for age,\nrestricted by % ill-defined deaths)"),
                     values=c(ifelse(grepl("univ", unique(all$model)), 21, 22)))+
  facet_wrap(~var, labeller = labeller(var=vars_labels), ncol=3) +
  guides(fill=guide_legend(override.aes = list(size=5)), color=F)+
  theme_bw() +
  labs(x="",
       y="RR (95% CI)",
       title="Comparison of Models")+
  theme(axis.ticks.x=element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        axis.text.x=element_text(size=10, face="bold",color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title.y=element_text(size=16, face="bold",color="black"),
        legend.title =element_text(size=16, face="bold",color="black"),
        legend.text  =element_text(size=14, color="black"),
        plot.title=element_text(face="bold", size=16),
        strip.text = element_text(size=12, face="bold",color="black")
  )
p
ggsave("results/ExtendedData10.pdf",p,  width=13.5, height=10)
ggsave("results/ExtendedData10.eps",p,  width=13.5, height=10)
all %>% select(-label) %>% 
  full_join(data.frame(var=unique(all$var), label=vars_labels)) %>% 
  full_join(data.frame(model=unique(all$model), model_label=c("Univariable\n", 
                                                               "Univariable\n(adjusted for age)",
                                                               "Univariable\n(adjusted for age and mortality)",
                                                               "Multivariable\n(adjusted for age)", 
                                                               "Multivariable\n(adjusted for age,\nrestricted by % ill-defined deaths)"))) %>% 
  mutate(model_label=sub("\\\n", " ", model_label),
         outcome_label=case_when(
           outcome=="type1cmnn" ~ "CMNN",
           outcome=="type2cancer" ~ "Cancer",
           outcome=="type3cvd" ~ "CVD/NCDs\n(Ref.)",
           outcome=="type4accident" ~ "Unintentional\nInjuries",
           outcome=="type5violent" ~ "Violence\nInjuries",
         )) %>% 
  select(label, model_label, outcome_label, beta, lci, uci) %>% 
  fwrite("Source/Source_Data_ExtendedData10.csv")

# plot results of growth sens analysis: first LE
models_growth_results<-models_growth %>% 
  group_by(age, sex, term) %>% 
  summarise(beta=mean(coef),
            se=sqrt(mean(se^2)+var(coef))) %>% 
  mutate(lci=beta-1.96*se,
         uci=beta+1.96*se,
         final=paste0(round(beta, digits=2),
                      " [",
                      round(lci, digits=2),
                      ";",
                      round(uci, digits=2),
                      "]")) %>% 
  arrange(sex, age, term) %>% 
  mutate(id=as.numeric(as.factor(sex))*13+
           as.numeric(as.factor(age))*3+
           as.numeric(as.factor(term))*-16)
models_growth_results$id<-
  as.numeric(as.factor(models_growth_results$sex))*13+
  as.numeric(as.factor(models_growth_results$age))*3+
  as.numeric(as.factor(models_growth_results$term))-16
models_growth_results<-bind_rows(models_growth_results,
                                 models_growth_results %>% 
                                   group_by(sex) %>% 
                                   summarise(id=mean(id)) %>% 
                                   mutate(label=ifelse(sex=="F", "Women", "Men"))) %>% 
  arrange(id) %>% 
  mutate(label=ifelse(is.na(label), "", label))

growth_le<-ggplot(models_growth_results %>% filter(!is.na(beta)), aes(x=id, y=beta))+
  geom_hline(yintercept = 0, lty=2)+
  geom_errorbar(aes(ymin=lci, ymax=uci, color=as.factor(age)))+
  geom_point(aes(fill=as.factor(age), shape=term), size=4, color="black") +
  scale_x_continuous(breaks=models_growth_results$id, 
                     labels=paste0(models_growth_results$label, "\n"))+
  scale_y_continuous(limits=c(-0.65, 0.65), breaks=c(-0.6, -0.3, 0, 0.3, 0.6))+
  scale_color_discrete(name="Life Expectancy at",
                       labels=c("Birth", paste0("Age ", c(20, 40, 60))))+
  scale_fill_discrete(name="Life Expectancy at",
                      labels=c("Birth", paste0("Age ", c(20, 40, 60))))+
  scale_linetype_discrete(name="Variable",
                          labels=c("Concurrent Growth", "Growth 5 years prior"))+
  scale_shape_manual(name="Variable",values=c(21, 22),
                     labels=c("Concurrent Growth", "Growth 5 years prior"))+
  labs(x="", y="Change in Life Expectancy (95% CI)",
       title="Life Expectancy")+
  guides(color = guide_legend(override.aes = list(shape = 16)),
         fill=F)+
  theme_classic() +
  theme(legend.position = "bottom",
        legend.box="vertical",
        axis.text.y=element_text(face="bold", size=14, color="black"),
        axis.title=element_text(face="bold", size=16, color="black"),
        axis.ticks.x=element_blank(),
        plot.title=element_text(face="bold", size=24, color="black"),
        axis.text.x=element_text(face="bold", size=14,color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(face="bold", size=16, color="black"))
growth_le

# growth sensitivity analysis in mortality
coefs_growth<-map_dfr(models_growth_mortality, function(model){
  coefs<-as.data.frame(tidy(model)) %>% 
    filter(grepl("growth_pct", term)) %>% 
    mutate(lci=exp(estimate-1.96*std.error),
           uci=exp(estimate+1.96*std.error),
           beta=exp(estimate),
           outcome=substr(term, 1, regexpr("\\:", term)-1),
           type=ifelse(grepl("prev", term), "previous", "concurrent")) %>% 
    select(outcome, type, beta, lci, uci) %>% 
    filter(outcome!="")
  coefs<-coefs %>% bind_rows(coefs %>% 
                               slice(1) %>% 
                               mutate(outcome="type3cvd", beta=1, lci=NA, uci=NA))
  coefs
}) %>% arrange(outcome, type) %>% 
  mutate(id=as.numeric(as.factor(outcome))*3+as.numeric(as.factor(type))-3)

growth_mortality<-ggplot(coefs_growth, aes(x=id, y=beta)) +
  geom_hline(yintercept = 1, lty=2, color="black")+
  geom_errorbar(aes(color=outcome,ymin=lci, ymax=uci)) +
  geom_point(aes(fill=outcome, shape=type), color="black", size=4) +
  scale_x_continuous(breaks=coefs_growth %>% filter(type=="concurrent") %>% pull(id)+.5, 
                     labels=c("CMNN","Cancer","CVD/NCDs\n(Ref.)","Unintentional\nInjuries","Violence\nInjuries"))+
  scale_y_continuous(trans = log_trans(), breaks = base_breaks(3),
                     limits=c(0.83, 1/0.83))+
  scale_fill_manual(values=colors, name="Cause")+
  scale_color_manual(values=colors, name="Cause",
                     labels=c("CMNN","Cancer","CVD/NCDs","Unintentional","Violence"))+
  scale_shape_manual(name="Variable",values=c(21, 22),
                     labels=c("Concurrent Growth", "Growth 5 years prior"))+
  labs(x="", y="RR of Population Growth (95% CI)",
       title="Proportionate Mortality")+
  guides(color = guide_legend(override.aes = list(shape = 16)),
         fill=F)+
  theme_classic() +
  theme(legend.position = "bottom",
        legend.box="vertical",
        axis.text.y=element_text(face="bold", size=14, color="black"),
        axis.title=element_text(face="bold", size=16, color="black"),
        axis.ticks.x=element_blank(),
        plot.title=element_text(face="bold", size=24, color="black"),
        axis.text.x=element_text(face="bold", size=14,color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(face="bold", size=16, color="black"))
growth_mortality
p<-arrangeGrob(grobs=list(growth_le, growth_mortality), ncol=2)
ggsave("results/ExtendedData4.pdf", p, width=20, height=7.5)
ggsave("results/ExtendedData4.eps", p, width=20, height=7.5)

bind_rows(models_growth_results %>% 
  mutate(outcome_label="LE", type=ifelse(term=="growth_pct", "concurrent", "previous")) %>% 
  select(outcome_label, age, sex, type, beta, lci, uci) %>% 
  filter(!is.na(beta)),
  coefs_growth %>% 
  mutate(outcome_label=case_when(
    outcome=="type1cmnn" ~ "CMNN",
    outcome=="type2cancer" ~ "Cancer",
    outcome=="type3cvd" ~ "CVD/NCDs\n(Ref.)",
    outcome=="type4accident" ~ "Unintentional\nInjuries",
    outcome=="type5violent" ~ "Violence\nInjuries",
  )) %>% 
  select(outcome_label, type, beta, lci, uci)) %>% 
  fwrite("Source/Source_Data_ExtendedData4.csv")

# MAPS
world <- ne_countries(scale = "medium", returnclass = "sf")
shp = st_read('../SALURBAL_DATA/SHPs/SALURBAL_L1AD_11_30_18/SALURBAL_L1AD_11_30_18.shp') %>% 
  ms_simplify() %>% 
  st_transform(crs=st_crs(world))
bbox<-st_bbox(shp)
cats=5
colors_map<-colorRampPalette(c("red", "yellow", "darkgreen"))(cats)

jenks_le<-ale_median %>% filter(age==0) %>% 
  select(SALID1, sex, le) %>% 
  group_by(sex) %>% 
  mutate(jenks=cut(le, breaks=classIntervals(le,n=cats, style="jenks")$brks, include.lowest = T)) %>% 
  select(-le) %>% spread(sex, jenks) %>% 
  mutate(across(where(is.factor), droplevels)) %>% 
  rename(men=`M`, women=`F`)
jenks_mortality<-mortality_cause %>% select(SALID1, p_cmnn, p_cancer, p_ncd, p_accident, p_violent) %>% 
  gather(type, prop, -SALID1) %>% 
  mutate(prop=round(prop*100)) %>% 
  group_by(type) %>% 
  mutate(jenks=cut(prop, breaks=classIntervals(prop,n=cats, style="jenks")$brks, include.lowest = T)) %>% 
  select(-prop) %>% spread(type, jenks) %>% 
  mutate(across(where(is.factor), droplevels))
shp2<-right_join(shp, full_join(jenks_le, jenks_mortality))

maps_le<-map(c("women", "men"), function(sex){
  temp<-shp2
  title<-case_when(
    sex=="men" ~ "Men",
    sex=="women" ~ "Women")
  colnames(temp)[colnames(temp)==sex]<-"jenks"
  ggplot(data = world) +
    geom_sf(fill="lightgray") +
    geom_sf(data=temp, 
            aes(geometry=geometry, fill=jenks, color=jenks))+
    # annotation_scale(location = "bl", width_hint = 0.5) +
    # annotation_north_arrow(location = "bl", which_north = "true",
    #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
    #                        style = north_arrow_fancy_orienteering)+
    scale_fill_manual(values=colors_map, name="Life Expectancy")+
    scale_color_manual(values=colors_map, name="Life Expectancy")+
    # scale_fill_gradient2(low="red", mid="yellow", high="darkgreen",
    #                      midpoint=median(shp2$men), name="LE in Men", n.breaks=7)+
    # scale_color_gradient2(low="red", mid="yellow", high="darkgreen",
    #                       midpoint=median(shp2$men), name="LE in Men", n.breaks=7)+
    guides(color=guide_legend(), fill=guide_legend())+
    coord_sf(xlim = c(bbox$xmin, bbox$xmax), 
             ylim = c(bbox$ymin, bbox$ymax), expand = expansion(mult=0.02))+
    labs(title=title)+
    theme_void() +
    theme(panel.background = element_rect(fill = "aliceblue"),
          plot.title=element_text(size=20, face="bold"),
          legend.position=c(0.1, 0.5))
})
maps_le[[2]]<-maps_le[[2]]+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)
pall<-arrangeGrob(grobs=maps_le, ncol=2)
ggsave("results/Figure2.pdf", pall, width=10, height=5)

ale_median %>% filter(age==0) %>% 
  ungroup() %>% 
  select(iso2, city_link, sex, le) %>% 
  rename(city=city_link) %>% 
  fwrite("Source/Source_Data_Figure2.csv")

maps_mortality<-map(cause_titles_coll, function(cause){
  temp<-shp2
  cause_name<-paste0("p_", cause)
  title<-case_when(
    cause=="cmnn" ~ "CMNN",
    cause=="cancer" ~ "Cancer",
    cause=="ncd" ~ "CVD/NCDs",
    cause=="accident" ~ "Unintentional Injuries",
    cause=="violent" ~ "Violent Injuries",
  )
  colnames(temp)[colnames(temp)==cause_name]<-"jenks"
  #colors_map<-colorRampPalette(c("white", colors[which(cause_titles_coll==cause)]))(cats)
  ggplot(data = world) +
    geom_sf(fill="lightgray") +
    geom_sf(data=temp, 
            aes(geometry=geometry, fill=jenks, color=jenks))+
    # annotation_scale(location = "bl", width_hint = 0.5) +
    # annotation_north_arrow(location = "bl", which_north = "true",
    #                        pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
    #                        style = north_arrow_fancy_orienteering)+
    scale_fill_manual(values=rev(colors_map), name="Proportion")+
    scale_color_manual(values=rev(colors_map), name="Proportion")+
    # scale_fill_gradient2(low="red", mid="yellow", high="darkgreen",
    #                      midpoint=median(shp2$men), name="LE in Men", n.breaks=7)+
    # scale_color_gradient2(low="red", mid="yellow", high="darkgreen",
    #                       midpoint=median(shp2$men), name="LE in Men", n.breaks=7)+
    guides(color=guide_legend(), fill=guide_legend())+
    coord_sf(xlim = c(bbox$xmin, bbox$xmax), 
             ylim = c(bbox$ymin, bbox$ymax), expand = expansion(mult=0))+
    labs(title=title)+
    theme_void() +
    theme(panel.background = element_rect(fill = "aliceblue"),
          plot.title=element_text(size=20, face="bold"),
          legend.text = element_text(size=14),
          legend.title = element_text(size=14, face="bold"),
          legend.position=c(0.1, 0.5))
})
maps_mortality[[5]]<-maps_mortality[[5]]+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)
pall<-arrangeGrob(grobs=maps_mortality, ncol=3)
ggsave("results/ExtendedData7.pdf", pall, width=24/5*3, height=10)
ggsave("results/ExtendedData7.eps", pall, width=24/5*3, height=10)



## supplement descriptive table
# Descriptive table
# first categorize
all_exposure$pop_group<-as.numeric(cut(all_exposure$pop_baseline, 
                                       breaks = c(0.5*10^5, 2.5*10^5, 5*10^5, 1*10^6, 5*10^6, 25*10^6), include.lowest = T))
table(all_exposure$pop_group)
# n cities, % each category, 
# AAMR, Size, Growth, Pop Density, Fragmentation, Street Connectivity, 
# Primary Ed, Piped Water, Sewage Connection, Overcrowding, SEI
cause_titles_coll<-c("cmnn", "cancer", "ncd", "accident", "violent")
dta<-full_join(all_exposure, mortality_cause) %>% 
  full_join(mortality_illdefined %>% select(SALID1, p_illdefined_disease, p_illdefined_injury))
# calculate LE and add
letable<-ale %>% filter(age==0) %>% 
  group_by(SALID1, sex) %>% 
  summarise(le=median(ale)) %>% 
  spread(sex, le) %>% 
  rename(leF=`F`, leM=`M`)
dta<-full_join(dta, letable)

vars<-c("leM", "leF", paste0("p_", cause_titles_coll),"p_illdefined_disease", "p_illdefined_injury",
        "adj_rate","pop_baseline", "growth_pct", "BECPOPDENSL1AD","BECPTCHDENSL1AD", "BECADINTDENSL1AD",
        "CNSMINPR_L1AD", "CNSWATINL1AD", "CNSSEWNETL1AD", "CNSCROWD3RML1AD", 
        "sei")
all_exposure$pop_group<-as.numeric(cut(all_exposure$pop_baseline, breaks = c(0*10^5, 2.5*10^5, 5*10^5, 1*10^6, 5*10^6, 25*10^6), include.lowest = T))
table(all_exposure$pop_group)
all_exposure %>% group_by(pop_group) %>% 
  summarise(min(pop_baseline), max(pop_baseline))
# rescale stuff for talbe
dta$pop_baseline<-dta$pop_baseline/1000000
dta$growth_pct<-dta$growth_pct*100
dta$BECPOPDENSL1AD<-dta$BECPOPDENSL1AD/1000
var<-vars[[1]]
#x<-dta %>% filter(pop_group==1)
options(scipen=999)
times100<-function(vector) vector*100
percentit<-function(vector) paste0(vector, "%")

table1<-dta %>% group_by(pop_group) %>% 
  group_modify(t1_function) %>% 
  spread(pop_group, V1) %>% 
  arrange(factor(var, levels = c("n", vars)), desc(var)) %>% 
  mutate('0'=dta %>% ungroup %>% do(t1_function(.)) %>% pull(V1)) %>% 
  select(1, '0', 2:6)
table1
fwrite(table1, "results/SupplementTable57.csv")
# total pop
all_exposure %>% pull(pop_baseline) %>% sum %>% `/` (1000000)


# Extended Data 1 
## ggtern for some reason conflicts with the themes of ggplot2, so if run, this section requires R restarting if you are running any other plot
library(ggtern)
mortality_cause<-mortality_cause %>% rowwise() %>% 
  mutate(p_ncds=sum(c(p_ncd, p_cancer)),
         p_injuries=sum(c(p_accident, p_violent)))
p1<-ggtern(data = mortality_cause, aes(z = p_cmnn, x = p_ncds, y = p_injuries)) +
  geom_point(aes(fill = country_name),size = 3,shape = 21,color = "black") +
  scale_fill_discrete(name="")+
  scale_T_continuous(limits=c(0, .6), name="+ Injuries")+
  scale_L_continuous(limits=c(0.4, 1), name="+ NCDs")+
  scale_R_continuous(limits=c(0, 0.6), name="+ CMNNs")+
  Larrowlab("% NCDs") + Rarrowlab("% CMNNs") + Tarrowlab("% Injuries")+
  theme_rgbw()+
  guides(fill=guide_legend(override.aes = list(size=6)))+
  theme(legend.position="bottom",
        legend.text = element_text(face="bold", size=14),
        legend.background = element_blank())
p1
ggsave("results/ExtendedData5.pdf", p1, width=7.5, height=7.5/1.33)
ggsave("results/ExtendedData5.eps", p1, width=7.5, height=7.5/1.33)


# finally, just making sure I am not disclosing anything I cant disclose in the source data
files<-list.files("Source/", pattern="csv", full.names =T)
map(files, function(file) fread(file) %>% colnames)
map(files, function(file) fread(file) %>% head)

