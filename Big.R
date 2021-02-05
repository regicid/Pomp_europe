library(panelPomp)
names = c("all","just_b","just_f","no_int","just_diff","gdp_only")
job = list()
job2 = list()
mifs_pomp = list()
name = "all"
PARAM = c("a","b","c","d","e","f","z","sigma","sigma_obs","N_0")
rwsd = rw.sd(a=.2,z = .2,sigma=.1,sigma_obs = .1,N_0=ivp(.1),d = .1,c=.1,b=.1,e=.1,f=.1)
Csnippet("double eps = rnorm(0,pow(sigma,2));
         N = z*N + pow(e,2)*cum + pow(f,2)*gdp*cum + a*gdp  + pow(d,2)*diff + pow(b,2)*gdp*diff + c + eps;
         ") -> evol_diff
mif3 <- function(a,sigma,N_0,sigma_obs,z,d,c,b,e,f){
  k = panelPomp::mif2(mf1,Nmif = 2500,shared = c(a = a,sigma = sigma,N_0 = N_0,sigma_obs = sigma_obs,z = z,d = d,b=b,c=c,e=e,f=f), specific = Model_diff@specific)
  
  return(k)
}

library(reshape)
library(tidyr)
library(dplyr)
library(doParallel)
library(rslurm)

Countries<-c("Italy","France","Poland","Germany","United Kingdom","Iberia","Russia","Belgium","Netherlands","Scandinavia")
Results = read.csv("Results.csv",row.names = 1)
Distances = read.csv("Distances.csv",row.names = 1)

R = Results[c("Countries","Date","Nobs")]
R  = spread(R,"Countries","Nobs")
rownames(R) = R$Date
R = R[-1]
R = R[Countries]

Csnippet("
         N = N_0;
         ") -> rinit
Csnippet("
         double eps_obs = rnorm(0,pow(sigma_obs,2));
         Nobs = N + eps_obs;
         ") -> rmeas


Csnippet("
         lik = dnorm(Nobs,N,pow(sigma_obs,2),give_log);
         ") -> dmeas
eval = function(i){
  mf = mifs_pomp[[i]]
  replicate(1000, mf %>% panelPomp::pfilter(Np = 10000) %>% panelPomp::logLik()) %>% logmeanexp(se=TRUE) -> ll
  return(c(coef(mf),loglik=ll[1],loglik.se=ll[2]))
}

Pomps = list()

for(country in Countries){
  data = dplyr::filter(Results,Countries==country)
  z = which(!is.na(data$gdp))
  Gdp = data[,c("Date","gdp")]
  Gdp$diff = rowSums(R/untable(Distances[country,],12)**2)
  data = data[c("Date","Nobs")]
  Gdp$cum = vector("numeric",length = nrow(Gdp))
  Gdp$cum[2:nrow(Gdp)] = cumsum(data$Nobs)[1:nrow(Gdp)-1]
  
  data %>%
    pomp(
      times="Date", t0=1300,
      rinit=rinit,
      covar = covariate_table(Gdp, times= "Date"),
      rprocess=discrete_time(evol_diff,delta.t = 50),
      dmeasure = dmeas,
      obsnames = c("Nobs"),
      statenames=c("N"),
      rmeasure=rmeas,
      paramnames=PARAM,
      covarnames = c("gdp","diff","cum")
    ) -> Pomps[[country]]
  
}

p = rep(0,length(PARAM))
names(p) = PARAM

Model_diff = panelPomp(Pomps,shared = p)
lower = c(a = 0.2,sigma=0.5,N_0 = 0,sigma_obs=0.2,z = 1.1,d = -0.6,b = -0.5,c =-1,e=-0.3,f=-0.3)
upper = c(a = 0.6,sigma=1, N_0 = 0,sigma_obs=0.6,z = 1.3,d = 1,b = 0.5,c = 0,e=-0.3,f=-0.3)

sobolDesign(lower = lower[PARAM], upper = upper[PARAM], nseq = 48) -> guesses


Model_diff %>%
  panelPomp::mif2(
    shared.start=unlist(guesses[6,]),
    specific.start = Model_diff@specific,
    Np=5000,
    Nmif=3,
    cooling.fraction.50=0.9,
    cooling.type="geometric",
    rw.sd= rwsd,
    pars = PARAM
  ) -> mf1
param = coef(mf1)

pkg = c("panelPomp")
options =  list("firstgen","benoit2c@gmail.com","ALL")
names(options) = c("partition","mail-user","mail-type")

job[[name]] = slurm_apply(f = mif3,params = guesses,jobname = paste(name,"0",sep="_"),nodes = 15,cpus_per_node = 12,pkgs = pkg,slurm_options = options,add_objects = c("mf1","Model_diff"))

mifs_pomp = get_slurm_out(job[[name]],wait = TRUE)















name = "just_b"
PARAM = c("a","b","c","z","sigma","sigma_obs","N_0")
rwsd = rw.sd(a=.2,z = .2,sigma=.1,sigma_obs = .1,N_0=ivp(.1),c=.1,b=.1)
Csnippet("double eps = rnorm(0,pow(sigma,2));
         N = z*N + a*gdp + pow(b,2)*gdp*diff + c + eps;
         ") -> evol_diff
mif3 <- function(a,sigma,N_0,sigma_obs,z,c,b){
  k = panelPomp::mif2(mf1,Nmif = 2500,shared = c(a = a,sigma = sigma,N_0 = N_0,sigma_obs = sigma_obs,z = z,b=b,c=c), specific = Model_diff@specific)
  
  return(k)
}
Pomps = list()

for(country in Countries){
  data = dplyr::filter(Results,Countries==country)
  z = which(!is.na(data$gdp))
  Gdp = data[,c("Date","gdp")]
  Gdp$diff = rowSums(R/untable(Distances[country,],12)**2)
  data = data[c("Date","Nobs")]
  Gdp$cum = vector("numeric",length = nrow(Gdp))
  Gdp$cum[2:nrow(Gdp)] = cumsum(data$Nobs)[1:nrow(Gdp)-1]
  
  data %>%
    pomp(
      times="Date", t0=1300,
      rinit=rinit,
      covar = covariate_table(Gdp, times= "Date"),
      rprocess=discrete_time(evol_diff,delta.t = 50),
      dmeasure = dmeas,
      obsnames = c("Nobs"),
      statenames=c("N"),
      rmeasure=rmeas,
      paramnames=PARAM,
      covarnames = c("gdp","diff","cum")
    ) -> Pomps[[country]]
  
}

p = rep(0,length(PARAM))
names(p) = PARAM

Model_diff = panelPomp(Pomps,shared = p)
lower = c(a = 0.2,sigma=0.5,N_0 = 0,sigma_obs=0.2,z = 1.1,d = -0.6,b = -0.5,c =-1,e=-0.3,f=-0.3)
upper = c(a = 0.6,sigma=1, N_0 = 0,sigma_obs=0.6,z = 1.3,d = 1,b = 0.5,c = 0,e=-0.3,f=-0.3)

sobolDesign(lower = lower[PARAM], upper = upper[PARAM], nseq = 48) -> guesses


Model_diff %>%
  panelPomp::mif2(
    shared.start=unlist(guesses[6,]),
    specific.start = Model_diff@specific,
    Np=5000,
    Nmif=3,
    cooling.fraction.50=0.9,
    cooling.type="geometric",
    rw.sd= rwsd,
    pars = PARAM
  ) -> mf1
param = coef(mf1)


options =  list("firstgen","benoit2c@gmail.com","ALL")
names(options) = c("partition","mail-user","mail-type")

job[[name]] =  slurm_apply(f = mif3,params = guesses,jobname = paste(name,"0",sep="_"),nodes = 15,cpus_per_node = 12,pkgs = pkg,slurm_options = options,add_objects = c("mf1","Model_diff"))

name = "just_f"
PARAM = c("a","c","f","z","sigma","sigma_obs","N_0")
rwsd = rw.sd(a=.2,z = .2,sigma=.1,sigma_obs = .1,N_0=ivp(.1),c=.1,f=.1)
Csnippet("double eps = rnorm(0,pow(sigma,2));
         N = z*N  + pow(f,2)*gdp*cum + a*gdp + c + eps;
         ") -> evol_diff
mif3 <- function(a,sigma,N_0,sigma_obs,z,c,f){
  k = panelPomp::mif2(mf1,Nmif = 2500,shared = c(a = a,sigma = sigma,N_0 = N_0,sigma_obs = sigma_obs,z = z,c=c,f=f), specific = Model_diff@specific)
  
  return(k)
}
Pomps = list()

for(country in Countries){
  data = dplyr::filter(Results,Countries==country)
  z = which(!is.na(data$gdp))
  Gdp = data[,c("Date","gdp")]
  Gdp$diff = rowSums(R/untable(Distances[country,],12)**2)
  data = data[c("Date","Nobs")]
  Gdp$cum = vector("numeric",length = nrow(Gdp))
  Gdp$cum[2:nrow(Gdp)] = cumsum(data$Nobs)[1:nrow(Gdp)-1]
  
  data %>%
    pomp(
      times="Date", t0=1300,
      rinit=rinit,
      covar = covariate_table(Gdp, times= "Date"),
      rprocess=discrete_time(evol_diff,delta.t = 50),
      dmeasure = dmeas,
      obsnames = c("Nobs"),
      statenames=c("N"),
      rmeasure=rmeas,
      paramnames=PARAM,
      covarnames = c("gdp","diff","cum")
    ) -> Pomps[[country]]
  
}

p = rep(0,length(PARAM))
names(p) = PARAM

Model_diff = panelPomp(Pomps,shared = p)
lower = c(a = 0.2,sigma=0.5,N_0 = 0,sigma_obs=0.2,z = 1.1,d = -0.6,b = -0.5,c =-1,e=-0.3,f=-0.3)
upper = c(a = 0.6,sigma=1, N_0 = 0,sigma_obs=0.6,z = 1.3,d = 1,b = 0.5,c = 0,e=-0.3,f=-0.3)

sobolDesign(lower = lower[PARAM], upper = upper[PARAM], nseq = 48) -> guesses


Model_diff %>%
  panelPomp::mif2(
    shared.start=unlist(guesses[6,]),
    specific.start = Model_diff@specific,
    Np=5000,
    Nmif=3,
    cooling.fraction.50=0.9,
    cooling.type="geometric",
    rw.sd= rwsd,
    pars = PARAM
  ) -> mf1
param = coef(mf1)

options =  list("firstgen","benoit2c@gmail.com","ALL")
names(options) = c("partition","mail-user","mail-type")

job[[name]] =  slurm_apply(f = mif3,params = guesses,jobname = paste(name,"0",sep="_"),nodes = 15,cpus_per_node = 12,pkgs = pkg,slurm_options = options,add_objects = c("mf1","Model_diff"))


name = "no_int"
PARAM = c("a","c","d","e","z","sigma","sigma_obs","N_0")
rwsd = rw.sd(a=.2,z = .2,sigma=.1,sigma_obs = .1,N_0=ivp(.1),d = .1,c=.1,e=.1)
Csnippet("double eps = rnorm(0,pow(sigma,2));
         N = z*N + pow(e,2)*cum + a*gdp  + pow(d,2)*diff + c + eps;
         ") -> evol_diff
mif3 <- function(a,sigma,N_0,sigma_obs,z,d,c,e){
  k = panelPomp::mif2(mf1,Nmif = 2500,shared = c(a = a,sigma = sigma,N_0 = N_0,sigma_obs = sigma_obs,z = z,d = d,c=c,e=e), specific = Model_diff@specific)
  
  return(k)
}
Pomps = list()

for(country in Countries){
  data = dplyr::filter(Results,Countries==country)
  z = which(!is.na(data$gdp))
  Gdp = data[,c("Date","gdp")]
  Gdp$diff = rowSums(R/untable(Distances[country,],12)**2)
  data = data[c("Date","Nobs")]
  Gdp$cum = vector("numeric",length = nrow(Gdp))
  Gdp$cum[2:nrow(Gdp)] = cumsum(data$Nobs)[1:nrow(Gdp)-1]
  
  data %>%
    pomp(
      times="Date", t0=1300,
      rinit=rinit,
      covar = covariate_table(Gdp, times= "Date"),
      rprocess=discrete_time(evol_diff,delta.t = 50),
      dmeasure = dmeas,
      obsnames = c("Nobs"),
      statenames=c("N"),
      rmeasure=rmeas,
      paramnames=PARAM,
      covarnames = c("gdp","diff","cum")
    ) -> Pomps[[country]]
  
}

p = rep(0,length(PARAM))
names(p) = PARAM

Model_diff = panelPomp(Pomps,shared = p)
lower = c(a = 0.2,sigma=0.5,N_0 = 0,sigma_obs=0.2,z = 1.1,d = -0.6,b = -0.5,c =-1,e=-0.3,f=-0.3)
upper = c(a = 0.6,sigma=1, N_0 = 0,sigma_obs=0.6,z = 1.3,d = 1,b = 0.5,c = 0,e=-0.3,f=-0.3)

sobolDesign(lower = lower[PARAM], upper = upper[PARAM], nseq = 48) -> guesses


Model_diff %>%
  panelPomp::mif2(
    shared.start=unlist(guesses[6,]),
    specific.start = Model_diff@specific,
    Np=5000,
    Nmif=3,
    cooling.fraction.50=0.9,
    cooling.type="geometric",
    rw.sd= rwsd,
    pars = PARAM
  ) -> mf1
param = coef(mf1)


options =  list("firstgen","benoit2c@gmail.com","ALL")
names(options) = c("partition","mail-user","mail-type")

job[[name]] =  slurm_apply(f = mif3,params = guesses,jobname = paste(name,"0",sep="_"),nodes = 15,cpus_per_node = 12,pkgs = pkg,slurm_options = options,add_objects = c("mf1","Model_diff"))



name = "just_diff"
PARAM = c("c","d","e","z","sigma","sigma_obs","N_0")
rwsd = rw.sd(z = .2,sigma=.1,sigma_obs = .1,N_0=ivp(.1),d = .1,c=.1,e=.1)
Csnippet("double eps = rnorm(0,pow(sigma,2));
         N = z*N + pow(e,2)*cum  + pow(d,2)*diff + c + eps;
         ") -> evol_diff
mif3 <- function(sigma,N_0,sigma_obs,z,d,c,e){
  k = panelPomp::mif2(mf1,Nmif = 2500,shared = c(sigma = sigma,N_0 = N_0,sigma_obs = sigma_obs,z = z,d = d,c=c,e=e), specific = Model_diff@specific)
  
  return(k)
}
Pomps = list()

for(country in Countries){
  data = dplyr::filter(Results,Countries==country)
  z = which(!is.na(data$gdp))
  Gdp = data[,c("Date","gdp")]
  Gdp$diff = rowSums(R/untable(Distances[country,],12)**2)
  data = data[c("Date","Nobs")]
  Gdp$cum = vector("numeric",length = nrow(Gdp))
  Gdp$cum[2:nrow(Gdp)] = cumsum(data$Nobs)[1:nrow(Gdp)-1]
  
  data %>%
    pomp(
      times="Date", t0=1300,
      rinit=rinit,
      covar = covariate_table(Gdp, times= "Date"),
      rprocess=discrete_time(evol_diff,delta.t = 50),
      dmeasure = dmeas,
      obsnames = c("Nobs"),
      statenames=c("N"),
      rmeasure=rmeas,
      paramnames=PARAM,
      covarnames = c("gdp","diff","cum")
    ) -> Pomps[[country]]
  
}

p = rep(0,length(PARAM))
names(p) = PARAM

Model_diff = panelPomp(Pomps,shared = p)
lower = c(a = 0.2,sigma=0.5,N_0 = 0,sigma_obs=0.2,z = 1.1,d = -0.6,b = -0.5,c =-1,e=-0.3,f=-0.3)
upper = c(a = 0.6,sigma=1, N_0 = 0,sigma_obs=0.6,z = 1.3,d = 1,b = 0.5,c = 0,e=-0.3,f=-0.3)

sobolDesign(lower = lower[PARAM], upper = upper[PARAM], nseq = 48) -> guesses


Model_diff %>%
  panelPomp::mif2(
    shared.start=unlist(guesses[6,]),
    specific.start = Model_diff@specific,
    Np=5000,
    Nmif=3,
    cooling.fraction.50=0.9,
    cooling.type="geometric",
    rw.sd= rwsd,
    pars = PARAM
  ) -> mf1
param = coef(mf1)


options =  list("firstgen","benoit2c@gmail.com","ALL")
names(options) = c("partition","mail-user","mail-type")

job[[name]] =  slurm_apply(f = mif3,params = guesses,jobname = paste(name,"0",sep="_"),nodes = 15,cpus_per_node = 12,pkgs = pkg,slurm_options = options,add_objects = c("mf1","Model_diff"))














name = "gdp_only"
PARAM = c("a","c","z","sigma","sigma_obs","N_0")
rwsd = rw.sd(a=.2,z = .2,sigma=.1,sigma_obs = .1,N_0=ivp(.1),c=.1)
Csnippet("double eps = rnorm(0,pow(sigma,2));
         N = z*N + a*gdp + c + eps;
         ") -> evol_diff
mif3 <- function(a,sigma,N_0,sigma_obs,z,c){
  k = panelPomp::mif2(mf1,Nmif = 2500,shared = c(a = a,sigma = sigma,N_0 = N_0,sigma_obs = sigma_obs,z = z,c=c), specific = Model_diff@specific)
  
  return(k)
}
Pomps = list()

for(country in Countries){
  data = dplyr::filter(Results,Countries==country)
  z = which(!is.na(data$gdp))
  Gdp = data[,c("Date","gdp")]
  Gdp$diff = rowSums(R/untable(Distances[country,],12)**2)
  data = data[c("Date","Nobs")]
  Gdp$cum = vector("numeric",length = nrow(Gdp))
  Gdp$cum[2:nrow(Gdp)] = cumsum(data$Nobs)[1:nrow(Gdp)-1]
  
  data %>%
    pomp(
      times="Date", t0=1300,
      rinit=rinit,
      covar = covariate_table(Gdp, times= "Date"),
      rprocess=discrete_time(evol_diff,delta.t = 50),
      dmeasure = dmeas,
      obsnames = c("Nobs"),
      statenames=c("N"),
      rmeasure=rmeas,
      paramnames=PARAM,
      covarnames = c("gdp","diff","cum")
    ) -> Pomps[[country]]
  
}

p = rep(0,length(PARAM))
names(p) = PARAM

Model_diff = panelPomp(Pomps,shared = p)
lower = c(a = 0.2,sigma=0.5,N_0 = 0,sigma_obs=0.2,z = 1.1,d = -0.6,b = -0.5,c =-1,e=-0.3,f=-0.3)
upper = c(a = 0.6,sigma=1, N_0 = 0,sigma_obs=0.6,z = 1.3,d = 1,b = 0.5,c = 0,e=-0.3,f=-0.3)

sobolDesign(lower = lower[PARAM], upper = upper[PARAM], nseq = 48) -> guesses


Model_diff %>%
  panelPomp::mif2(
    shared.start=unlist(guesses[6,]),
    specific.start = Model_diff@specific,
    Np=1000,
    Nmif=100,
    cooling.fraction.50=0.9,
    cooling.type="geometric",
    rw.sd= rwsd,
    pars = PARAM
  ) -> mf1
param = coef(mf1)

data = data.frame(mf1@pconv.rec)
data$iteration = as.numeric(row.names(mf1@pconv.rec))
data = gather(data,"key","value",-iteration)
ggplot(data,aes(iteration,value)) + geom_line()+ facet_wrap(vars(key),scales = "free")


options =  list("secondgen","benoit2c@gmail.com","ALL")
names(options) = c("partition","mail-user","mail-type")

job[[name]] =  slurm_apply(f = mif3,params = guesses,jobname = paste(name,"0",sep="_"),nodes = 15,cpus_per_node = 12,pkgs = pkg,slurm_options = options,add_objects = c("mf1","Model_diff"))

for(name in names){
  mifs_pomp[[name]] = get_slurm_out(job,wait = TRUE)
  file1 = paste("~/mifs_pomp_",name,sep="")
  file1 = paste(file1, ".RDS",sep="")
  saveRDS(mifs_pomp[[name]],file1)
}

for(name in names[1:5]){
  mifs = mifs_pomp[[name]]
  d = data.frame(i = 1:length(mifs))
  job2[[name]] = slurm_apply(f = eval,params = d,jobname = paste("evaluation",name,sep="_"),nodes = 15,cpus_per_node = 12,pkgs = pkg,slurm_options = options,add_objects = c("mifs"))
}

name = "gdp_only"
mifs = mifs_pomp[[name]]
d = data.frame(i = 1:length(mifs))
job2[[name]] = slurm_apply(f = eval,params = d,jobname = paste("evaluation",name,sep="_"),nodes = 15,cpus_per_node = 12,pkgs = pkg,slurm_options = options,add_objects = c("mifs"))
options =  list("secondgen","benoit2c@gmail.com","ALL")
names(options) = c("partition","mail-user","mail-type")

for(name in names){
  estimates = get_slurm_out(job2[[name]], outtype = "table",wait=TRUE)
  file2 = paste("~/estimates_",name,sep="")
  file2 = paste(file2, ".csv",sep="")
  write.csv(estimates,file2)
}
