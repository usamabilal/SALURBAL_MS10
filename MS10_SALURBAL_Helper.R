# helper functions
# harmonic mean
hmean<-function(a){
  1/mean(1/a)  
}
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
# ggplot functions
base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}
gg_color_hue <- function(n, alpha=1) {
  hues = seq(15, 375, length = n + 1)
  t<-hcl(h = hues, l = 65, c = 100)[1:n]
  adjustcolor(t, alpha)
}
get_legend<-function(plot){
  grobs<-ggplotGrob(plot)$grobs
  legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]  
  return(legend)
}

# get ICC from an lme4 model
ICC_lmer <- function(out) {
  varests <- as.data.frame(VarCorr(out))
  varests<-varests[,"vcov"]
  return(varests[1]/sum(varests))
}