library(devtools)
load_all()

set.seed('1234')


### DEMO ----
# just a short demo to show how the implemented classes work
if(FALSE){

  # let's generate some example stimuli input
  u = sample(0:1, 20, replace=T, prob=c(.5, .5))
  u = c(u, sample(0:1, 10, replace=T, prob=c(.2, .8)) )

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

### SIMULATIONS SETUP ----



### SIMULATIONS PART I - BINARY STIMULI ----



### end ----

