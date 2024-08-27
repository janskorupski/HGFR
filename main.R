library(devtools)
load_all()

set.seed('1234')


### DEMO ----
# just a short demo to show how the implemented classes work
if(T){

  # let's generate some example stimuli input
  u = sample(0:1, 20, replace=T, prob=c(.5, .5))
  u = c(u, sample(0:1, 20, replace=T, prob=c(.2, .8)) )

  # now we can create a class with that input
  mathys = HGF_binary(u=u)
  # the init function automatically calls mathys = fit(mathys, method="VB") if u is supplied
  plot(mathys)

  # the simulations have to be run explicitly, as they are time-consuming
  mathys = fit(mathys, method="MCMC")
  # this plots the empirical kernel estimation of distribution vs VB estimation
  plot.distributions(mathys, timestamps = c(2,4,17,19), levels = c(2,3))


  # the same goes for the continuous version
  u = c( rnorm(100, sd = 1) + 1:100/10, rnorm(100, sd=0.2)+10 )

  aaa = HGF_continuous(u=u)
  plot(aaa)

  aaa = fit(aaa, method="MCMC")
  plot.distributions(aaa, timestamps = c(1,10,40,80), levels = c(1,2,3), logscale=T)

}


### SIMULATION SETUP BINARY STIMULI ----

if( file.exists("simTable") ){
  readRDS("simTable")
}else{
  kappaSpace = (1:20)/10
  omegaSpace = (-20:20)/10
  thetaSpace = (1:10)/10

  KLlvl2 = NA
  KLlvl3 = NA
  DecisionLogLikelyhood = NA
  DecisionAccuracy = NA

  completed = FALSE

  simTable = expand.grid(kappa=kappaSpace,
                         omega=omegaSpace,
                         theta=thetaSpace,
                         KLlvl2=KLlvl2,
                         KLlvl3=KLlvl3,
                         DecisionLogLikelyhood=DecisionLogLikelyhood,
                         DecisionAccuracy=DecisionAccuracy,
                         completed = completed)
}

if(file.exists("simArchives")){
  readRDS("simArchives")
}else{
  simArchives = list()
}


compress_to_means = function(table){ apply(table, 2, mean) }
compress_to_vars = function(table){ apply(table, 2, var) }
compressHGFInstance = function(instance){

  cnames = colnames(instance@simulations[[1]])

  means = sapply(instance@simulations, compress_to_means)
  means = t(means)
  colnames(means) = cnames

  vars = sapply(instance@simulations, compress_to_vars)
  vars = t(vars)
  colnames(vars) = cnames

  instance@simulations = list()
  instance@simulations[[1]] = means
  instance@simulations[[2]] = vars

  return(instance)
}

KL_div = function(distribution, sample){return(NA)}
get_dec_acc = function(instance){
  vb_dec = s(instance@moments[[1]][,2]) > 0.5
  mc_dec = instance@sim_preds_ > 0.5
  acc = mean(vb_dec == mc_dec)
  return(acc)
}


### SIMULATIONS PART I - BINARY STIMULI ----

# set computation parameters (timers)
{
time_to_pass = 60
save_interval = 10

break_time = 60*15
work_time = 60*60*4
}

# initialise helping variables
{
iter_no = which.min(simTable$completed)

start_time = Sys.time()
end_time = start_time + time_to_pass
last_save_time = Sys.time()
last_break_time = Sys.time()

pb = txtProgressBar(min=0, max= nrow(simTable), style=3)
while( (iter_no <= nrow(simTable)) & (Sys.time() < end_time) ){

  pars = as.list(simTable[iter_no, c("kappa","omega","theta")])
  instance = HGF_binary(u=u,
                        parameters = pars)
  instance = fit(instance, method="MCMC")

  # store data in the table
  simTable[iter_no,"DecisionAccuracy"] = get_dec_acc(instance)

  # store the simulations in archive list
  simArchives[[iter_no]] = compressHGFInstance(instance)

  # save the calculations
  if(Sys.time() - last_save_time > save_interval){
    saveRDS(simTable,"simTable")
    saveRDS(simArchives,"simArchives")
  }

  # update iteration number
  simTable$completed[iter_no] = TRUE
  iter_no = which.min(simTable$completed)

  setTxtProgressBar(pb, sum(simTable$completed))

  if(Sys.time() - last_break_time > work_time){
    last_break_time = Sys.time()
    Sys.sleep(break_time)
  }
}
close(pb)
}

### PLOTTING RESULTS PART I - BINARY STIMULI ----

ggplot(simTable[simTable$completed,]) +
  geom_tile(aes(x=kappa,y=omega,fill=DecisionAccuracy))

agree_stat = function(instance){
  vb_dec = s(instance@moments[[1]][,2]) - 0.5
  mc_dec = instance@sim_preds_ - 0.5
  return(data.frame(kappa=instance@parameters$kappa,
                    time=1:length(instance@u),
                    value=abs(vb_dec-mc_dec)*sign(vb_dec*mc_dec)))
}

lapply(simArchives[1:5], agree_stat)
plot_data_1 = Reduce(rbind, lapply(simArchives[1:5], agree_stat))
ggplot(plot_data_1) +
  geom_tile(aes(y=kappa, x=time, fill=value))

### end ----


