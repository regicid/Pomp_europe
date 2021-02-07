library(panelPomp)
library(reshape)
library(rslurm)
library(doParallel)
library(dplyr)
library(foreach)

options =  list("dellgen","benoit2c@gmail.com","ALL")
names(options) = c("partition","mail-user","mail-type")
cpus = c(12,12,16,16,40)
names(cpus) = c("firstgen","fastgen","secondgen","lastgen","dellgen")
cores = cpus[options[[1]]]

##Importations
Countries<-c("Italy","France","Poland","Germany","United Kingdom","Iberia","Russia","Belgium","Netherlands","Scandinavia")
Results = read.csv("Results.csv",row.names = 1)
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
         double eps_obs = rnorm(1,pow(sigma_obs,2));
         Nobs = N*eps_obs;
         ") -> rmeas
Csnippet("
         lik = dnorm(Nobs,N,pow(N,2)*pow(sigma_obs,2),give_log);
         ") -> dmeas
Csnippet("double eps = rnorm(0,pow(sigma,2));
         N = z*N + pow(e,2)*cum + pow(f,2)*gdp*cum + a*gdp  + pow(d,2)*diff + pow(b,2)*gdp*diff + c + eps;
         ") -> evol_diff


#mif3 <- function(a,sigma,N_0,sigma_obs,z,d,c,b,e,f){

names = c("all","just_b","just_f","no_int","just_diff","gdp_only")
PARAM = c("bla","a","b","c","d","e","f","z","sigma","sigma_obs","N_0")
rwsd = rw.sd(a=.1,b=.1,c=.1,d=.1,e=.1,f=.1,z=.1,sigma=.1,sigma_obs=.1,N_0 = ivp(.1))
job = list()
job2 = list()
mifs_pomp = list()
unused_parameters = list()
unused_parameters[[1]] = c(1)
unused_parameters[[2]] = c(1,5,6,7)
unused_parameters[[3]] = c(1,3,5,6)
unused_parameters[[4]] = c(1,3,7)
unused_parameters[[5]] = c(1,2,3,7)
unused_parameters[[6]] = c(3,5,6,7)
names(unused_parameters) = names

submit_job <- function(unused_parameters,nmif=2,np=2,
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
  lower = c(a=.0,b=-.5,c=0,d=-.5,e=-.2,f=-.1,z=.6,sigma=.4,sigma_obs=.1,N_0 = -.1)
  upper = c(a=.3,b=.5,c=.2,d=.5,e=.2,f=.1,z=1,sigma=.6,sigma_obs=.3,N_0 = .1)
  lower[unused_parameters[[name]]-1] = 0
  upper[unused_parameters[[name]]-1] = 0
  sobolDesign(lower = lower[PARAM[-1]], upper = upper[PARAM[-1]], nseq = n) -> guesses
  rwsd@call[unused_parameters[[name]]-1] = 0
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
  mif3 <- function(a,sigma,N_0,sigma_obs,z,d,c,b,e,f){
    k = panelPomp::mif2(mf1,Nmif = nmif,shared = c(a = a,sigma = sigma,N_0 = N_0,sigma_obs = sigma_obs,z = z,d = d,b=b,c=c,e=e,f=f), specific = Model_diff@specific)
    return(k)
  }
  pkg = c("panelPomp")
  job[[name]] = slurm_apply(f = mif3,params = guesses,jobname = paste(name,"0",sep="_"),nodes = 15,cpus_per_node = cores,pkgs = pkg,slurm_options = options,add_objects = c("mf1","Model_diff"))
}

for(name in names){
  submit_job(unused_parameters = unused_parameters[[name]])
}

for(name in names){
  mifs_pomp[[name]] = get_slurm_out(job[[name]],wait = TRUE)
  file1 = paste("mifs_pomp_",name,sep="")
  file1 = paste(file1, ".RDS",sep="")
  saveRDS(mifs_pomp2[[name]],file1)
}


for(name in names){
  mifs = mifs_pomp[[name]]
  d = data.frame(i = 1:length(mifs))
  job2[[name]] = slurm_apply(f = eval,params = d,jobname = paste("evaluation",name,sep="_"),nodes = 15,cpus_per_node = 16,pkgs = pkg,slurm_options = options,add_objects = c("mifs"))
}

for(name in names){
  estimates = get_slurm_out(job2[[name]], outtype = "table",wait=TRUE)
  file2 = paste("estimates_",name,sep="")
  file2 = paste(file2, ".csv",sep="")
  write.csv(estimates,file2)
}
