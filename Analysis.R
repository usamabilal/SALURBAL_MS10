# Main Tables and Figures
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
source("MS10_SALURBAL_Helper.R")
select<-dplyr::select
load("analytic files/all_data_mortality_population_corrected_level1_age14.RData")
load("analytic files/l1s.rdata")
load("analytic files/Life_Expectancy_All_Iterations.rdata")
load("analytic files/MS10_exposure_data.RData")
# get median LE for basic descriptives
ale_median<-ale %>% group_by(SALID1, sex, age) %>% 
  summarise(le=median(ale)) %>% 
  spread(sex, le) %>% 
  rename(leM='M',
         leF='F')
vars<-c("leM", "leF", "pop_baseline", "growth_pct", "BECPOPDENSL1AD","BECPTCHDENSL1AD", "BECADINTDENSL1AD",
        "CNSMINPR_L1AD", "CNSWATINL1AD", "CNSSEWNETL1AD", "CNSCROWD3RML1AD", 
        "sei")
all_exposure$pop_group<-as.numeric(cut(all_exposure$pop_baseline, breaks = c(0*10^5, 2.5*10^5, 5*10^5, 1*10^6, 5*10^6, 25*10^6), include.lowest = T))
table(all_exposure$pop_group)
all_exposure %>% group_by(pop_group) %>% 
  summarise(min(pop_baseline), max(pop_baseline))
# Number of Cities and  Number of Large Cities (>2M?)
ncities<-as.numeric(table(all_exposure$pop_group))
# rescale stuff
all_exposure$pop_baseline<-all_exposure$pop_baseline/1000000
all_exposure$growth_pct<-all_exposure$growth_pct*100
all_exposure$BECPOPDENSL1AD<-all_exposure$BECPOPDENSL1AD/1000
all_exposure<-full_join(all_exposure, ale_median %>% filter(age==0) %>% select(-age))
var<-vars[[1]]
#temp<-all_exposure %>% filter(pop_group==1)

t1_function<-function(x, y){
  t<-lapply(vars, function(var){
    if (is.na(var)){
      paste0("")
    } else {
      summary<-x %>% pull(var) %>% quantile(probs=c(.5, .25, .75), na.rm=T) %>% format(digits=1, nsmall=1)
      paste0(summary[2], " [",
             summary[1], ";",
             summary[3], "]")
    }
  })
  t<-append(t, nrow(x), after=0)
  t<-as.data.frame(do.call(rbind, t), stringsAsFactors = F)
  t$var<-c("n", vars)
  t
}
table1<-all_exposure %>% group_by(pop_group) %>% 
  group_modify(t1_function) %>% 
  spread(pop_group, V1) %>% 
  arrange(factor(var, levels = c("n", vars)), desc(var)) %>% 
  mutate('0'=all_exposure %>% ungroup %>% do(t1_function(.)) %>% pull(V1)) %>% 
  select(1, '0', 2:6)
table1
fwrite(table1, "results/ExtendedData9.csv")

# median LE in long format ofr some of the descriptives
ale_median<-ale %>% group_by(SALID1, sex, age) %>% 
  summarise(le=median(ale),
            lci=quantile(ale, probs=0.025),
            uci=quantile(ale, probs=0.975),
            se=sqrt(var(ale)),
            dif_ci=uci-lci) %>% 
  mutate(rse=se/le*100) %>% 
  left_join(l1s)

# 95 CI ranges
ale_median %>% 
  filter(age==0) %>% 
  pull(dif_ci) %>% 
  quantile(probs=c(0.01, 0.025, 0.05, 0.25, .5, .75, .95, .975, .99))
ale_median %>% group_by(age, sex) %>% 
  summarise(b5=sum(dif_ci<5),
            total=n()) %>% 
  mutate(prop_b5=b5/total)

# total variability
# ICCs taking into consideration variability between iterations
iccs<-ale %>% 
  group_by(age, sex) %>% 
  group_modify(~{
    .x<-.x %>% 
      left_join(all_exposure) 
    m<-lmer(ale~1+(1|iso2/SALID1), data=.x, weights=pop_baseline)
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
  
full_join(iccs %>%  select(age, sex, ICCiteration) %>% 
            spread(sex, ICCiteration) %>% 
            rename(IterationF=F, IterationM=M),
          iccs %>%  select(age, sex, ICCcity) %>% 
            spread(sex, ICCcity) %>% 
            rename(CityF=F, CityM=M)) %>% 
  full_join(iccs %>%  select(age, sex, ICCcountry) %>%  
            spread(sex, ICCcountry) %>% 
            rename(CountryF=F, CountryM=M)) %>% 
  select(age, IterationF,CityF, CountryF, IterationM, CityM, CountryM) %>% 
  fwrite(file="results/ExtendedData2.csv")

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
ranges[order(ranges$sex, ranges$range),]
ranges[order(ranges$sex, ranges$median),]

# top cities for discussion
ale_median %>% filter(age==0) %>% left_join(l1s) %>%
  group_by(sex) %>% arrange(desc(le)) %>% slice(1:5)
ale_median %>% filter(age==0) %>% left_join(l1s) %>%
  group_by(sex) %>% arrange((le)) %>% slice(1:5)

# figure 1
p1<-ale_median %>% 
  filter(age==0) %>% 
  group_by(age, sex) %>% 
  group_map(~{
    age<-.y$age
    sex<-.y$sex
    ylim<-c(min(floor(ale_median[ale_median$age==age, "le"]/5)*5),
            max(ceiling(ale_median[ale_median$age==age, "le"]/5)*5))
    ylab<-ifelse(age==0, "Life Expectancy at Birth", paste0("Life Expectancy at Age ", age))
    # icc<-iccs[iccs$sex==sex&iccs$age==age,"icc"]
    # icc<-paste0(format(icc*100, nsmall=1, digits=1), "%")
    #title<-paste0("Life Expectancy in ", ifelse(sex=="M", "Men", "Women"))
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
ggsave("results/Figure1.pdf", pLEB, width=15, height=5)


# comparison with other countries.
# # figure 1 with country-level references
country_list<-l1s %>% filter(!duplicated(iso2)) %>% select(country_name, iso2) %>%
  rename(country=country_name)
dta1<-(read_excel("Other_Data/UNDP Life Tables/WPP2019_MORT_F17_2_ABRIDGED_LIFE_TABLE_MALE.xlsx", skip=16)) %>%
  mutate(sex="M")
dta2<-(read_excel("Other_Data/UNDP Life Tables/WPP2019_MORT_F17_3_ABRIDGED_LIFE_TABLE_FEMALE.xlsx", skip=16)) %>%
  mutate(sex="F")
#temp<-dta1
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

country_les %>% filter(age==0, sex=="F", le_country>83.7) %>% arrange(le_country) %>% print(n=40)
country_les %>% filter(age==0, sex=="M", le_country>78) %>% arrange(le_country) %>% print(n=40)

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

p1<-ale_median %>% 
  left_join(country_les %>% rename(country_name=country)) %>% 
  group_by(age, sex) %>% 
  group_map(~{
    age<-.y$age
    sex<-.y$sex
    ylim<-c(min(floor(ale_median[ale_median$age==age, "le"]/5)*5),
            max(ceiling(ale_median[ale_median$age==age, "le"]/5)*5))
    ylab<-ifelse(age==0, "Life Expectancy at Birth", paste0("Life Expectancy at Age ", age))
    # icc<-iccs[iccs$sex==sex&iccs$age==age,"icc"]
    # icc<-paste0(format(icc*100, nsmall=1, digits=1), "%")
    #title<-paste0("Life Expectancy in ", ifelse(sex=="M", "Men", "Women"))
    title<-ifelse(sex=="M", "Men", "Women")
    ylab2<-ifelse(sex=='M', ylab, "")
    ylab<-ifelse(sex=='M', "", ylab)
    x<-.x %>% left_join(l1s)
    ggplot(x, aes(x=iso2, y=le, group=iso2)) +
      geom_boxplot(aes(group=as.factor(iso2)), fill=NA, outlier.color = NA, width=0.5)+
      geom_jitter(aes(fill=as.factor(iso2)), width=0.1, height=0, alpha=1, size=2, 
                  color="black", pch=21) +
      geom_point(data=.x %>% filter(!duplicated(iso2)),
                 aes(x=iso2, y=le_country),
                 pch=17, color="black",size=4)+
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
ggsave("results/ExtendedData1.pdf", p, width=15, height=5*3)

# RSE Figure
sex_label<-c("Women", "Men")
names(sex_label)<-c("F", "M")
age_label<-paste0("LE at ", c("Birth", paste0("age ", c(20, 40, 60))))
names(age_label)<-c(0, 20, 40, 60)
figure_rse<-ggplot(ale_median, aes(x=iso2, y=rse, group=iso2)) +
  geom_boxplot(aes(group=as.factor(iso2)), fill=NA, outlier.color = NA, width=0.5)+
  geom_jitter(aes(fill=as.factor(iso2)), width=0.1, height=0, alpha=1, size=2, 
              color="black", pch=21) +
  guides(color=F, fill=F, size=F)+
  labs(x="",
       y="Relative Standard Error",
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
ggsave("results/ExtendedData10.pdf", figure_rse, width=20, height=10)

income_labels<-c("High-income countries", "Upper-middle-income countries",
                 "Middle-income countries","Lower-middle-income countries",
                 "Low-income countries")
#.x<-ale_median %>% filter(sex=="M", age==0);.y<-data.frame(sex="M", age=0, stringsAsFactors=F)
figureS12<-ale_median %>% 
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
            legend.text=element_text(size=16, color="black"),
            legend.title=element_text(face="bold", size=16, color="black"))+
      theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
  })
legend<-get_legend(figureS12[[1]])
figureS12<-lapply(figureS12, function(xx) xx+guides(color=F, fill=F))
figureS12<-arrangeGrob(grobs=list(figureS12[[1]], figureS12[[2]]), ncol=2)
figureS12<-arrangeGrob(grobs=list(figureS12, legend), nrows=2,
                       heights=c(10, 1))
plot(figureS12)
ggsave("results/Figure2.pdf", figureS12, width=20, height=7.5)



# exposures: distributions, correlations
vars_cont<-vars[-c(1, 2, 12)]
cols<-gg_color_hue(9)
iso2list<-c("AR", "BR", "CL", "CO", "CR", "MX", "PA", "PE", "SV")

xtitles<-c("Population (log scale)", "% Growth in Population over 5 years",
           "Population Density (Population/km2)","Patch Density (Patches/km2)",
           "Intersection Density (Intersections/km2)",
           "% with Completed primary Education or Above", "% Households with Piped Water in the Household","% of Households with Connection to the Sewage Network","% Overcrowded Households (>3 people per room)",
           "Social Environment Index (SD)")
plottitles<-c("City Size", "Population Growth",
              "Population Density", "Patch Fragmentation", 
              "Street Connectivity",
              "Educational Attainment", "Water Access","Sewage Connection", "Overcrowding",
              "Social Environment Index")
p<-lapply(seq_along(vars_cont), function(i){
  print(i)
  var<-vars_cont[[i]]
  temp<-all_exposure
  if (var=="pop_baseline"){
    temp$pop_baseline<-temp$pop_baseline*1000000
  }
  temp<-temp[temp$iso2%in%iso2list,]
  temp$label<-temp$iso2
  title<-paste0(plottitles[[i]])
  temp2<-temp
  temp2$iso2<-"ZZ"
  temp2$label<-"ZZ"
  temp<-rbind(temp, temp2)
  
  p<-ggplot(temp, aes_string(x="iso2", y=var, group="label")) +
    geom_boxplot(aes(group=as.factor(label)), fill=NA, outlier.color = NA, width=0.5)+
    geom_jitter(aes(fill=as.factor(label)), width=0.1, height=0, pch=21, color="black") +
    guides(fill=F, size=F, alpha=F)+
    xlab("") + ylab(xtitles[[i]]) +
    # scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,0.71),
    #                    breaks=seq(0, 0.7, by=0.1))+
    scale_x_discrete(labels=c(iso2list, "Overall"))+
    scale_fill_manual(values=c(cols, "black"))+
    scale_alpha_manual(values=c(0.5, 0.2))+
    ggtitle(title)+
    theme_classic() +
    theme(legend.position = "bottom",
          legend.key.width = unit(50, "points"),
          axis.text.x=element_text(size=13, face="bold",color="black"),
          axis.text.y=element_text(size=13, face="bold",color="black"),
          axis.title.y=element_text(size=16, face="bold",color="black"),
          plot.title=element_text(face="bold", size=16))
  if (var=="pop_baseline"){
    p<-p+
      scale_y_log10(breaks=c(10^5, 2.5*10^5, 5*10^5, 
                             10^6, 2.5*10^6, 5*10^6,
                             10^7, 2.5*10^7),
                    labels=c("100K", "250K", "500K",
                             "1M", "2.5M", "5M",
                             "10M", "25M")
      )+
      annotation_logticks(sides="l")
  }
  p
})
pmortality<-arrangeGrob(grobs=p[1:9], ncol=3)
ggsave("results/ExtendedData8.pdf",pmortality, width=20, height=20/1.6, units="in", scale=1.5, limitsize = F)


# outcome models
# run a thousand models of LE ~ predictor
# getting var SDs
vars_cont<-vars<-c("pop_baseline", "growth_pct", "BECPOPDENSL1AD","BECPTCHDENSL1AD", "BECADINTDENSL1AD",
                   "CNSMINPR_L1AD", "CNSWATINL1AD", "CNSSEWNETL1AD", "CNSCROWD3RML1AD", 
                   "sei")

select<-dplyr::select
vars_sds<-vars_cont[-c(1, length(vars_cont))]
sds<-all_exposure %>% ungroup() %>% 
  mutate(BECPOPDENSL1AD=BECPOPDENSL1AD*1000) %>% 
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
models<-ale_long %>% 
  #filter(id%in%1:10) %>% 
  group_by(sex, age, id) %>% 
  group_modify(~{
    temp<-.x
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
    model1
  })

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
results_univ_table<-results_univ %>% 
  filter(age==0) %>% 
  ungroup() %>% 
  select(sex, var, final) %>% 
  spread(sex, final) %>% 
  select(var, M, F) %>% 
  arrange(factor(var, levels=vars_cont), desc(var)) %>% 
  mutate(sds=sds_char) %>% 
  select(var, sds, M, F)

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
  mutate(sds=sds_char[which(vars_cont%in%vars_multiv)]) %>% 
  select(var, sds,multiM, multiF)

results_univ_table
write.csv(results_univ_table, "results/Table1_univariate.csv")
results_multiv_table
write.csv(results_multiv_table, "results/Table2_multiv.csv")


results_var_univ<-models %>% 
  group_by(age, sex, var) %>% 
  summarise(sigma2=mean(sigma2),
            sigma2_empty=mean(sigma2_empty),
            tau00=mean(tau00),
            tau00_empty=mean(tau00_empty)) %>% 
  rowwise() %>% 
  mutate(city_change=mean((sigma2-sigma2_empty)/sigma2_empty)*100,
         country_change=mean((tau00-tau00_empty)/tau00_empty)*100,
         total_change=mean(((tau00+sigma2)-(tau00_empty+sigma2_empty))/
                                (tau00_empty+sigma2_empty))*100) %>% 
  filter(age==0) %>% 
  ungroup() %>% 
  select(var, sex, city_change, country_change, total_change) %>% 
  gather(var2, value, city_change, country_change, total_change) %>% 
  mutate(id=paste0(var2, sex)) %>% 
  select(var, id, value) %>% 
  spread(id, value) %>% 
  select(var, 
         city_changeM, city_changeF,
         country_changeM, country_changeF,
         total_changeM, total_changeF)
results_var_multiv<-models %>% 
  filter(var=="pop_baseline") %>% 
  group_by(age, sex) %>% 
  summarise(sigma2_multiv=mean(sigma2_multiv),
            sigma2_empty=mean(sigma2_empty),
            tau00_multiv=mean(tau00_multiv),
            tau00_empty=mean(tau00_empty)) %>% 
  rowwise() %>% 
  filter(age==0) %>% 
  mutate(city_change=((sigma2_multiv-sigma2_empty)/sigma2_multiv)*100,
            country_change=((tau00_multiv-tau00_empty)/tau00_empty)*100,
            total_change=(((tau00_multiv+sigma2_multiv)-(tau00_empty+sigma2_empty))/
                                (tau00_empty+sigma2_empty))*100) %>% 
  mutate(var="multiv") %>% 
  gather(var2, value, city_change, country_change, total_change) %>% 
  mutate(id=paste0(var2, sex)) %>% 
  select(var, id, value) %>% 
  spread(id, value) %>% 
  select(var, 
         city_changeM, city_changeF,
         country_changeM, country_changeF,
         total_changeM, total_changeF)

tableS5<-bind_rows(results_var_univ, results_var_multiv) %>% 
  arrange(factor(var, levels=c(vars_cont, "multiv")), desc(var)) %>% 
  mutate_at(-1, ~paste0(round(.x, digits=1),"%"))
tableS5

#univ all ages
#temp<-results_univ %>% filter(sex=="M")
figureS7<-results_univ %>% 
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
legend<-get_legend(figureS7[[1]])
figureS7<-lapply(figureS7, function(xx) xx+guides(color=F, fill=F))
figureS7<-arrangeGrob(grobs=list(figureS7[[2]], figureS7[[1]]), ncol=2)
figureS7<-arrangeGrob(grobs=list(figureS7, legend), nrows=2,
                      heights=c(10, 1))
plot(figureS7)
ggsave("results/ExtendedData3.pdf", figureS7, width=20, height=7.5)

#multiv all ages
#temp<-results_multiv %>% filter(sex=="M")
figureS8<-results_multiv %>% 
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
            axis.text.x=element_text(face="bold", size=14, angle=90, hjust=1, color="black"),
            legend.text=element_text(size=14, color="black"),
            legend.title=element_text(face="bold", size=16, color="black"))
  })

legend<-get_legend(figureS8[[1]])
figureS8<-lapply(figureS8, function(xx) xx+guides(color=F, fill=F))
figureS8<-arrangeGrob(grobs=list(figureS8[[2]], figureS8[[1]]), ncol=2)
figureS8<-arrangeGrob(grobs=list(figureS8, legend), nrows=2,
                      heights=c(10, 1))
plot(figureS8)
ggsave("results/ExtendedData3.pdf", figureS8, width=20, height=7.5)

save(models, file="results/models_main_analysis.rdata")
# does pop growth previous to this affect?
#.x<-ale_long %>% filter(sex=="F", age==0)
models_growth<-ale_long %>% #filter(id%in%1:10) %>% 
  group_by(sex, age, id) %>% 
  group_modify(~{
    temp<-.x
    temp<-temp %>% left_join(all_exposure %>% select(SALID1,growth_pct_prev ))
    temp<-temp %>% filter(!is.na(growth_pct_prev))
    temp$growth_pct<-scale(temp$growth_pct, center=T, scale=T)
    temp$growth_pct_prev<-scale(temp$growth_pct_prev, center=T, scale=T)
    f<-as.formula(paste0("ale~", "growth_pct", "+ BECPCTURBANL1AD+(1|iso2)"))
    model1<-lmer(f, data=temp, na.action = "na.omit")
    coefs1<-model1 %>% tidy %>% filter(term=="growth_pct")
    f<-as.formula(paste0("ale~", "growth_pct_prev", "+ BECPCTURBANL1AD+(1|iso2)"))
    model2<-lmer(f, data=temp, na.action = "na.omit")
    coefs2<-model2 %>% tidy %>% filter(term=="growth_pct_prev")
    bind_rows(coefs1, coefs2) %>% select(term, estimate, std.error) %>% 
      rename(coef=estimate, se=std.error)
  })

# sensitivity analysis (same correction and 3 other methods)
## same correction method, but restricting to those with >=90% [same analysis as above]
include<-correction %>% 
  filter(a/(a+b) > 0.9) %>% select(SALID1, sex) 
models_sens10<-ale_long %>% #filter(id%in%1:10) %>% 
  right_join(include) %>% 
  group_by(sex, age, id) %>% 
  group_modify(~{
    temp<-.x
    temp<-temp %>% mutate(pop_baseline=log(pop_baseline*1000000))
    temp<-temp %>% mutate_at(vars_cont[-1], scale, center=T, scale=T) %>% 
      mutate_at(vars_cont[-1], as.numeric)
    library(broom)
    library(lme4)
    temp<-.x
    temp<-temp %>% mutate(pop_baseline=log(pop_baseline*1000000))
    temp<-temp %>% mutate_at(vars_cont[-1], scale, center=T, scale=T) %>% 
      mutate_at(vars_cont[-1], as.numeric)
    f<-as.formula(paste0("ale~", paste(vars_multiv, collapse="+"), 
                         "+ BECPCTURBANL1AD+(1|iso2)"))
    model<-lmer(f, data=temp, na.action = "na.omit")
    coefs<-model %>% tidy %>% filter(term!="(Intercept)", 
                                     term!="BECPCTURBANL1AD",
                                     effect=="fixed")
    model2<-data.frame(var=vars_multiv,
                       coef_multiv=coefs %>% pull(estimate), 
                       se_multiv=coefs %>% pull(std.error))   
    model2 %>% mutate(sens="sens10")
  })
#new_ale<-ale_by_ucnt[[1]]
load("analytic files/Life_Expectancy_by_UCNT.rdata")
plan(multiprocess)
models_sensother<-future_map_dfr(ale_by_ucnt, function(new_ale){
  ale_long_ucnt<-full_join(new_ale, all_exposure %>% select(SALID1, iso2, BECPCTURBANL1AD, vars_cont))
  ale_long_ucnt %>% 
    #filter(id%in%1:10) %>% 
    group_by(sex, age, id) %>% 
    group_modify(~{
      library(broom.mixed)
      library(lme4)
      temp<-.x
      temp<-temp %>% mutate(pop_baseline=log(pop_baseline*1000000))
      temp<-temp %>% mutate_at(vars_cont[-1], scale, center=T, scale=T) %>% 
        mutate_at(vars_cont[-1], as.numeric)
      f<-as.formula(paste0("ale~", paste(vars_multiv, collapse="+"), 
                           "+ BECPCTURBANL1AD+(1|iso2)"))
      model<-lmer(f, data=temp, na.action = "na.omit")
      coefs<-model %>% tidy %>% filter(term!="(Intercept)", 
                                       term!="BECPCTURBANL1AD",
                                       effect=="fixed")
      model2<-data.frame(var=vars_multiv,
                         coef_multiv=coefs %>% pull(estimate), 
                         se_multiv=coefs %>% pull(std.error))   
      model2 %>% mutate(sens=unique(ale_long_ucnt$method))
    })
})

save(models_growth, models_sens10,models_sensother, file="results/models_sens_analysis.rdata")
#load("results/models_main_analysis.rdata")
#load("results/models_sens_analysis.rdata")
#growth sens analysis
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

ggplot(models_growth_results %>% filter(!is.na(beta)), aes(x=id, y=beta))+
  geom_hline(yintercept = 0, lty=2)+
  geom_errorbar(aes(ymin=lci, ymax=uci, color=as.factor(age)))+
  geom_point(aes(fill=as.factor(age), shape=term), size=4, color="black") +
  scale_x_continuous(breaks=models_growth_results$id, labels=models_growth_results$label)+
  ylim(c(-0.2, 1))+
  scale_color_discrete(name="Life Expectancy at",
                       labels=c("Birth", paste0("Age ", c(20, 40, 60))))+
  scale_fill_discrete(name="Life Expectancy at",
                      labels=c("Birth", paste0("Age ", c(20, 40, 60))))+
  scale_linetype_discrete(name="Variable",
                          labels=c("Concurrent Growth", "Growth 5 years prior"))+
  scale_shape_manual(name="Variable",values=c(21, 22),
                       labels=c("Concurrent Growth", "Growth 5 years prior"))+
  labs(x="", y="Change in Life Expectancy (95% CI)")+
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
ggsave("results/ExtendedData6.pdf", width=10, height=7.5)


models_sens<-bind_rows(models_sens10, models_sensother)

results_sens<-models_sens %>%
  group_by(age, sens, sex, var) %>% 
  summarise(beta=mean(coef_multiv),
            se=sqrt(mean(se_multiv^2)+var(coef_multiv))) %>% 
  ungroup() %>% 
  # re-scale population coefficients
  mutate(beta=ifelse(var=="pop_baseline", beta*log(pop_coefficient), beta),
         se=ifelse(var=="pop_baseline", se*log(pop_coefficient), se)) %>% 
  mutate(lci=beta-1.96*se,
         uci=beta+1.96*se,
         final=paste0(round(beta, digits=2),
                      " [",
                      round(lci, digits=2),
                      ";",
                      round(uci, digits=2),
                      "]")) %>% 
  bind_rows(results_multiv %>% 
              mutate(sens="main")) %>% arrange(sex, age, var, sens) 
# .x<-results_sens %>% filter(age==0)%>% filter(sex=="M");.y<-data.frame(sex="M", stringsAsFactors=F)
ylim<-c(min(results_sens$lci), max(results_sens$uci))
methods_figure<-c("main", "sens10", methods)
figureS9<-results_sens %>% filter(age==0) %>% group_by(sex) %>% 
  group_map(~{
    xx<-.x
    xx$id1<-as.numeric(factor(xx$var, levels = c(vars_multiv)))
    xx$sens<-factor(xx$sens, levels = methods_figure)
    xx$id2<-as.numeric(xx$sens)
    xx<-xx[order(xx$id1, xx$id2),]
    xx$id<-xx$id1*(length(methods_figure)+1)+xx$id2-(length(methods_figure)+1)
    labels<-xx %>% filter(!duplicated(var)) %>% 
      ungroup() %>% 
      select(var) %>% 
      mutate(label=c("Population", "Growth",
                     "Pop. Density", "Fragmentation",
                     "Connectivity","Social Index")) %>% 
      mutate(id2=2)
    xx<-xx %>% left_join(labels) %>% 
      mutate(label=ifelse(is.na(label), "", label))
    #title<-paste0("Life Expectancy at Birth in ", ifelse(.y$sex=="M", "Men", "Women"))
    xx$main<-grepl("main", xx$sens)
    model_label<-c("Main model\n(best fitting ages)", 
                   "Restriction model\n(>90% coverage by main model)",
                   "Hill age-bands (30-65)",
                   "Murray age-bands (50-70)")
    title<-ifelse(unique(.y$sex)=="M", "Men", "Women")
    ggplot(xx, aes(x=id, y=beta, group=id)) +
      geom_hline(yintercept = 0, lty=2, color="black")+
      geom_errorbar(aes(ymin=lci, ymax=uci, color=sens)) +
      geom_point(aes(fill=sens), size=2, color="black", pch=21) +
      scale_x_continuous(breaks=xx$id, labels=xx$label)+
      ylim(ylim)+
      # scale_color_discrete(name="Model",
      #                      labels=model_label)+
      # scale_fill_discrete(name="Model",
      #                      labels=model_label)+
      scale_color_discrete(name="Type", labels=model_label)+
      scale_fill_discrete(name="Type", labels=model_label)+
      scale_shape_manual(name="Type",values=c(21:24))+
      xlab("") + ylab("Change in Life Expectancy (95% CI)")+
      theme_classic() +
      ggtitle(title)+
      guides(fill=guide_legend(nrow=1,byrow=TRUE),
             color=guide_legend(nrow=1,byrow=TRUE))+
      theme(legend.position = "bottom",
            legend.box = "vertical",
            axis.text.y=element_text(face="bold", size=14, color="black"),
            axis.title=element_text(face="bold", size=16, color="black"),
            axis.ticks.x=element_blank(),
            plot.title=element_text(face="bold", size=24, color="black"),
            axis.text.x=element_text(face="bold", size=14, angle=90, hjust=1, color="black"),
            legend.text=element_text(size=14, color="black"),
            legend.title=element_text(face="bold", size=16, color="black"))
  })
legend<-get_legend(figureS9[[1]])
figureS9<-lapply(figureS9, function(xx) xx+guides(color=F, fill=F))
figureS9<-arrangeGrob(grobs=list(figureS9[[2]], figureS9[[1]]), ncol=2)
figureS9<-arrangeGrob(grobs=list(figureS9, legend), nrows=2,
                      heights=c(10, 1))
plot(figureS9)
ggsave("results/ExtendedData5.pdf", figureS9, width=20, height=7.5)



# plot undercounting
correction_plot<-correction %>% 
  mutate(phi2=a/(a+b)) %>% 
  left_join(l1s)
ylim<-c(min(correction_plot$phi2), 1)
sex_label<-c("Women", "Men")
names(sex_label)<-c("F", "M")
methods_label<-c("Best-fitting age bands", "Hill bands (30-65)", "Murray bands (50-70)")
names(methods_label)<-correction_plot$ages %>% unique
ggplot(correction_plot, aes(x=iso2, y=phi2, group=iso2)) +
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
  facet_grid(sex~ages, labeller = labeller(sex=sex_label,
                                           ages=methods_label))+
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.width = unit(50, "points"),
        axis.text.x=element_text(size=20, color="black"),
        axis.text.y=element_text(size=14, color="black"),
        axis.title.y=element_text(size=14, color="black", face="bold"),
        plot.title=element_text(face="bold", size=25),
        strip.text =element_text(size=20, color="black", face="bold"),
        strip.background = element_blank())
ggsave("results/ExtendedData7.pdf", width=15, height=10)

# save app data
# need: SALID1, age, sex, le, lci, uci, iso2, city_link, country_name
dta<-ale %>% 
  group_by(SALID1, sex, age) %>% 
  summarise(le=median(ale),
            lci=quantile(ale, probs=0.025),
            uci=quantile(ale, probs=0.975)) %>% 
  left_join(l1s)
save(dta, file="MS10/data.rdata")



