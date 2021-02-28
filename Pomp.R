library(tidyr)
library(panelPomp)
library(reshape)
library(rslurm)
library(doParallel)
library(dplyr)
library(foreach)
analysis = "MLE2"
options =  list("dellgen","benoit2c@gmail.com","ALL")
names(options) = c("partition","mail-user","mail-type")
pkg = c("panelPomp")
cpus = c(12,12,16,16,40,48)
names(cpus) = c("firstgen","fastgen","secondgen","lastgen","dellgen","gpu")
#partitions = c(rep("dellgen",7),"gnt")
##Importations
Countries<-c("Italy","France","Poland","Germany","United Kingdom","Iberia","Russia","Belgium","Netherlands","Scandinavia")
Results = read.csv("Results.csv",row.names = 1)
Results$Nobs = log(Results$Nobs+1)
Results$Nobs = Results$Nobs/sd(Results$Nobs)
Results$gdp = Results$gdp - min(Results$gdp)
Distances = read.csv("Distances.csv",row.names = 1)

#Create the diffusion matrix
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
         Nobs = N+eps_obs;
         ") -> rmeas
Csnippet("
         lik = dnorm(Nobs,N,pow(sigma_obs,2),give_log);
         ") -> dmeas
Csnippet("double eps = fmax(rnorm(1,pow(sigma,2)),0);
         double eps2 = rnorm(0,pow(sigma2,2));
         N = z*N*eps + pow(e,2)*cum + pow(f,2)*gdp*cum + a*gdp  + pow(d,2)*diff + pow(b,2)*gdp*diff + c + eps2;
         ") -> evol_diff



names = c("all","just_b","just_f","no_int","just_diff","gdp_only","expon","linear","constant")
PARAM = c("bla","a","b","c","d","e","f","z","sigma","sigma_obs","N_0","sigma2")
job = list()
job2 = list()
mifs_pomp = list()
unused_parameters = list()
unused_parameters[[1]] = c(1,12)
unused_parameters[[2]] = c(1,5,6,7,12)
unused_parameters[[3]] = c(1,3,5,6,12)
unused_parameters[[4]] = c(1,3,7,12)
unused_parameters[[5]] = c(1,2,3,7,12)
unused_parameters[[6]] = c(1,3,5,6,7,12)
unused_parameters[[7]] = c(1,2,3,5,6,7,12)
unused_parameters[[8]] = c(1,3,5,6,7,8,12)
unused_parameters[[9]] = c(1,2,3,5,6,7,8,12)


names(unused_parameters) = names
names = names[-3][-3]
submit_job <- function(nmif=10000,np=15000,
                       cooling_fraction=.95,n=40){
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
        times="Date", t0=1300,rinit=rinit,
        covar = covariate_table(Gdp, times= "Date"),
        rprocess=discrete_time(evol_diff,delta.t = 50),
        dmeasure = dmeas,obsnames = c("Nobs"),
        statenames=c("N"),rmeasure=rmeas,
        paramnames=PARAM[-1],covarnames = c("gdp","diff","cum")
      ) -> Pomps[[country]]}
  p = rep(0,length(PARAM[-1]))
  names(p) = PARAM[-1]
  Model_diff = panelPomp(Pomps,shared = p)
  lower = c(a=.0,b=-.5,c=0,d=-.5,e=-.2,f=-.1,z=.6,sigma=.4,sigma_obs=.1,N_0 = -.1,sigma2=.1)
  upper = c(a=.3,b=.5,c=.2,d=.5,e=.2,f=.1,z=1,sigma=.6,sigma_obs=.3,N_0 = .1,sigma2=.3)
  lower[unused_parameters[[name]]-1] = 0
  upper[unused_parameters[[name]]-1] = 0
  sobolDesign(lower = lower[PARAM[-1]], upper = upper[PARAM[-1]], nseq = n) -> guesses
  rwsd = rw.sd(a=.1,b=.1,c=.1,d=.1,e=.1,f=.1,z=.2,sigma=.1,sigma_obs=.1,N_0 = ivp(.1),sigma2=.1)
  rwsd@call[PARAM[unused_parameters[[name]]][-1]] = 0
  Model_diff %>%
    panelPomp::mif2(
      shared.start=unlist(guesses[6,]),
      specific.start = Model_diff@specific,
      Np=np,
      Nmif=3,
      cooling.fraction.50=cooling_fraction,
      cooling.type="geometric",
      rw.sd= rwsd,
      pars = PARAM
    ) -> mf1
  mif3 <- function(a,sigma,N_0,sigma_obs,z,d,c,b,e,f,sigma2){
    k = panelPomp::mif2(mf1,Nmif = nmif,shared = c(a = a,sigma = sigma,N_0 = N_0,sigma_obs = sigma_obs,z = z,d = d,b=b,c=c,e=e,f=f,sigma2=sigma2), specific = Model_diff@specific)
    return(k)
  }
  
  k = slurm_apply(f = mif3,params = guesses,jobname = paste(name,analysis,sep="_"),nodes = 15,cpus_per_node = cpus[options[[1]]],pkgs = pkg,slurm_options = options,add_objects = c("mf1","Model_diff"))
  return(k)
}

for(name in names){
  job[[name]] = submit_job()
}

for(name in names){
  mifs_pomp[[name]] = get_slurm_out(job[[name]],wait = TRUE)
  #cleanup_files(job[[name]])
  file1 = paste(paste("mifs_pomp_",analysis,sep="_"),name,sep="")
  file1 = paste(file1, ".RDS",sep="")
  saveRDS(mifs_pomp[[name]],file1)
}


for(name in names){
  mifs = mifs_pomp[[name]]
  d = data.frame(i = 1:length(mifs))
  job2[[name]] = slurm_apply(f = eval,params = d,jobname = paste(paste("ev",name,sep="_"),analysis,sep="_"),nodes = 15,cpus_per_node = cpus[options[[1]]],pkgs = pkg,slurm_options = options,add_objects = c("mifs"))
}

for(name in names){
  estimates = get_slurm_out(job2[[name]], outtype = "table",wait=TRUE)
  #cleanup_files(job2[[name]])
  file2 = paste(paste("estimates_",analysis,sep="_"),name,sep="")
  file2 = paste(file2, ".csv",sep="")
  write.csv(estimates,file2)
}

