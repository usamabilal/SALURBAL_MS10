# helper functions
# gompertz approximation
gompertz_approximation<-function(temp, max_original=80, min=45, max=90, just5age=F,
                                 returnage14=F){
  # From Chetty et al, JAMA 2016
  # We estimate the Gompertz parameters α and β for each income percentile, gender, and
  # year using maximum likelihood estimation (MLE). Under the Gompertz model, the number of
  # deaths at each age follows a Binomial distribution with a probability of death given by eα + βage
  # Hence, the Gompertz model is equivalent to a generalized linear model (GLM) with a binomial
  # probability distribution and a log link function with a single predictor (age). We estimate this
  # GLM separately for each income percentile, gender, and year using Stata’s glm function
  temp2<-temp %>% filter(age%in%c(min:(max_original-1))) %>% 
    mutate(deaths=ifelse(deaths==0, 1, deaths),
           lnrate=log(deaths/pop))
  r2<-temp2 %>% lm(formula=lnrate~age) %>% glance %>% pull(r.squared)
  model<-temp2 %>% glm(formula=cbind(deaths, pop-deaths)~age,
                family=binomial(link="log"))
  if (just5age==T){
    predict<-expand.grid(age=seq(min, max, by=5))  
  } else {
    predict<-expand.grid(age=min:max)
  }
  temp2<-predict %>% mutate(predicted_rate=model %>% augment(newdata=predict) %>% pull(.fitted) %>% exp) %>% 
    full_join(temp, by="age") %>% arrange(age) %>% 
    mutate(original_rate=deaths/pop,
           final_rate=ifelse(is.na(original_rate), predicted_rate, original_rate),
           final_rate=ifelse(age==max_original, predicted_rate, final_rate),
           r2=r2) %>% 
    select(age, final_rate, r2)
}

# ggplot functions (thanks stackoverflow...)
# https://stackoverflow.com/questions/14255533/pretty-ticks-for-log-normal-scale-using-ggplot2-dynamic-not-manual
base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}
# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n, alpha=1) {
  hues = seq(15, 375, length = n + 1)
  t<-hcl(h = hues, l = 65, c = 100)[1:n]
  adjustcolor(t, alpha)
}
# unclear where i got this from
get_legend<-function(plot){
  grobs<-ggplotGrob(plot)$grobs
  legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]  
  return(legend)
}
# get ICC from an 2-level lme4 model
ICC_lmer <- function(out) {
  varests <- as.data.frame(VarCorr(out))
  varests<-varests[,"vcov"]
  return(varests[1]/sum(varests))
}
# harmonic mean
hmean<-function(a){
  1/mean(1/a)  
}
# sample variance
hmean_var<-function(x){
  var<-sum((x-hmean(x))^2)/(length(x)-1)
  return(var)
}
scale_new<-function(x){
  as.numeric(scale(x, center=T, scale=T))
}
# function to take coefficients and make a table with CIs
t3_function<-function(x, y){
  #xx<-coefs[coefs$exp=="mortality",]
  xx<-x
  xx$exp<-exp(xx$est)
  xx$lci<-exp(xx$est-1.96*xx$se)
  xx$uci<-exp(xx$est+1.96*xx$se)
  xx$coefs<-paste0(format(xx$exp, digits=2, nsmall=2)," (",
                   format(xx$lci, digits=2, nsmall=2), ";",
                   format(xx$uci, digits=2, nsmall=2), ")")
  data.frame(cmnn=xx$coefs[1], 
             cancer=xx$coefs[2], 
             cvd="1 (Ref.)", 
             accident=xx$coefs[3], 
             violent=xx$coefs[4], stringsAsFactors = F)
}
# table 1 function
t1_function<-function(x, y){
  t<-lapply(vars, function(var){
    if (is.na(var)){
      paste0("")
    } else if (grepl("^p_", var)){
      summary<-x %>% pull(var) %>% quantile(probs=c(.5, .25, .75), na.rm=T) %>% times100 %>% format(digits=1, nsmall=1) 
      paste0(summary[1], " [",
             summary[2], ";",
             summary[3], "]")
    } else if (grepl("^CNS", var)){
      summary<-x %>% pull(var) %>% quantile(probs=c(.5, .25, .75), na.rm=T)  %>% format(digits=1, nsmall=1) 
      paste0(summary[1], " [",
             summary[2], ";",
             summary[3], "]")
    } else {
      summary<-x %>% pull(var) %>% quantile(probs=c(.5, .25, .75), na.rm=T) %>% format(digits=1, nsmall=1)
      paste0(summary[1], " [",
             summary[2], ";",
             summary[3], "]")
    }
  })
  t<-append(t, nrow(x), after=0)
  t<-as.data.frame(do.call(rbind, t), stringsAsFactors = F)
  t$var<-c("n", vars)
  t
}
ICC <- function(out) {
  varests <- as.data.frame(VarCorr(out))
  varests<-varests[,"vcov"]
  return(varests[1]/sum(varests))
}
