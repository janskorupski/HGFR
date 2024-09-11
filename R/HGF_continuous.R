
# Define the HGF_continuous subclass inheriting from HGF
setClass("HGF_continuous",
         contains = "HGF",
         prototype = list(

           no_of_levels    = 3,
           no_of_moments   = 2,
           parameters      = list("kappa"=rep(1.4,3), "omega"=rep(-2.2, 3), "theta"=0.5, "alpha"=1),
           priors          = list("mu" = c(0, 0, 0), "sigma" = c(1, 1, 1)),
           update_formulas = list(),
           u               = NA_real_,
           moments         = list(),
           simulations     = list()

         )
)

# Constructor function for HGF_continuous objects
HGF_continuous <- function(
    u = integer(),
    no_of_levels = 3,
    no_of_moments = 2,
    parameters = list("kappa"=rep(1.4,3), "omega"=rep(-2.2, 3), "theta"=0.5, "alpha"=1),
    priors = list("mu" = c(0, 0, 0), "sigma" = c(1, 1, 1)),
    update_formulas = list()){
  u = c(NaN, u)

  # creating the moments' matrix
  moments = list()
  for( i in 1:no_of_moments ){
    moments[[i]] = matrix(ncol=no_of_levels)
  }
  names(moments) = c("mu", "sigma", "skewness", "kurtosis")[1:min(4,no_of_moments)]

  # assigning the priors as the first row of moments' matrix
  if( length(priors) != 0 ){
    for( i in 1:no_of_moments ){

      if( length(priors[[i]]) != no_of_levels ){
        throw("Wrong number of priors")
      }

      moments[[i]][1,] = priors[[i]]
    }
  }

  object = new("HGF_continuous",
               no_of_levels = no_of_levels,
               no_of_moments = no_of_moments,
               parameters = parameters,
               priors = priors,
               update_formulas = update_formulas,
               u = u,
               moments = moments,
               simulations = list())

  if( length(object@u) > 1 ){
    object = fit(object, method="VB")
  }

  return(object)

}

HGF_continuous_VB = function(learner){


  if( ! is.nan(learner@u[1]) ){
    learner@u = c(NaN, learner@u)
  }

  for( moment_no in 1:learner@no_of_moments){
    learner@moments[[moment_no]] =
      matrix(nrow = (length(learner@u)),
             ncol = learner@no_of_levels)
    learner@moments[[moment_no]][1,] = learner@priors[[moment_no]]
  }

  v = matrix(nrow = (length(learner@u)),
             ncol = learner@no_of_levels)
  pi_hat = matrix(nrow = (length(learner@u)),
                  ncol = learner@no_of_levels)
  delta = matrix(nrow = (length(learner@u)),
                 ncol = learner@no_of_levels)

  kappa = c(learner@parameters$kappa)
  omega = c(learner@parameters$omega)
  theta = learner@parameters$theta
  alpha = learner@parameters$alpha

  mu = learner@moments$mu
  sigma = learner@moments$sigma

  for( t in 2:(length(learner@u)) ){

    # LEVEL 1
    if( is.na(learner@u[t]) ){ # TWEAK THIS LATER -------------------------------<<<
      print("dupa")
      learner@moments[[1]][t, 1] = learner@moments[[1]][t-1, 1]
    }

    v[t, 1] = exp( kappa[1] * mu[t-1, 2] + omega[1] )


    pi_hat[t, 1] = 1/( sigma[t-1, 1] + v[t, 1] )

    delta_u = learner@u[t] - mu[t-1, 1]  ### ------------------------------------<<< there's probably a mistake in the Mathys paper

    sigma[t, 1] = 1/(pi_hat[t, 1] + 1/alpha)
    mu[t, 1] = mu[t-1, 1] + (sigma[t, 1]/ alpha)*delta_u

    delta[t, 1] = ( (sigma[t, 1] + (mu[t, 1] - mu[t-1, 1])**2 ) /
                      (sigma[t-1, 1] + v[t, 1]) ) - 1

    # OTHER LEVELS (assumed normal) # based on Mathys et al. 2014 equations 9-14
    for( level in 2:learner@no_of_levels){

      # eq (11)
      if( level != learner@no_of_levels ){
        v[t, level] = exp( kappa[level] * mu[t-1, level+1] + omega[level] )
      }else{
        v[t, level] = theta
      }

      # eq (12)
      mu_hat = mu[t-1, level]

      # eq (13)
      pi_hat[t, level] = 1/( sigma[t-1, level] + v[t, level] )

      # eq (10)
      pi = pi_hat[t, level] +
        (1/2)*( kappa[level-1] *
                  v[t, level-1] *
                  pi_hat[t, level-1] )**2 *
        (1 + (1 - ((sigma[t-1, level-1])/( v[t, level-1] )))*delta[t, level-1] )

      sigma[t, level] = 1/pi

      # eq (9)
      mu[t, level] = mu_hat +
        (1/2)*kappa[level-1]*v[t, level-1]*
        (pi_hat[t, level-1]*sigma[t, level])*
        delta[t, level-1]

      # eq (14)
      delta[t, level] = ( (sigma[t, level] + (mu[t, level] - mu_hat)**2 ) /
                            (sigma[t-1, level] + v[t, level]) ) - 1


    }

    learner@moments$mu = mu
    learner@moments$sigma = sigma

  }
  return( learner )
}

HGF_continuous_BF = function(learner){
  N = 1e4

  if( ! is.nan(learner@u[1]) ){
    learner@u = c(NaN, learner@u)
  }

  kappa = learner@parameters$kappa
  omega = learner@parameters$omega
  theta = learner@parameters$theta
  alpha = learner@parameters$alpha

  ## First time based on priors
  learner@simulations[[1]] = list()

  for( level in 1:learner@no_of_levels){

    learner@simulations[[1]][[paste0("x", level)]] = rnorm(N,
                                                           mean = learner@priors$mu[level],
                                                           sd = sqrt(learner@priors$sigma[level]))

  }
  learner@simulations[[1]][["w"]] = rep(1, N)
  learner@simulations[[1]] = as.data.frame(learner@simulations[[1]])


  for( t in 2:(length(learner@u)) ){

    prob = learner@simulations[[t-1]][["w"]]

    prob[is.na(prob)] = 0
    prob = prob/sum(prob)
    sample_idx = sample(1:N, size=N, replace = T, prob = prob)

    if(any(learner@simulations[[t-1]]$w[sample_idx] == 0)){
      print("dupa")
    }

    res = list()

    x_prev = learner@simulations[[t-1]][[paste0('x', learner@no_of_levels)]][sample_idx]
    res[[paste0("x", learner@no_of_levels)]] =
      rnorm(N,
            mean = x_prev,
            sd = sqrt(theta))

    for( level in (learner@no_of_levels-1):1){

      x_prev = learner@simulations[[t-1]][[paste0('x', level)]][sample_idx]

      x_above = res[[paste0("x", level+1)]]
      vars = exp( kappa[level]*x_above + omega[level])

      res[[paste0("x", level)]] = rnorm(N,
                                        mean = x_prev,
                                        sd = sqrt(vars))
    }


    res[["w"]] = dnorm(learner@u[t],
                       mean = res[["x1"]],
                       sd = sqrt(alpha))

    learner@simulations[[t]] = as.data.frame(res)

  }
  return( learner )
}


# Define a method for fit specific to HGF_continuous
setMethod("fit",
          signature(object = "HGF_continuous", method = "character"),
          function(object, method="VB", ...) {
            # Implementation of fit method specific to HGF_continuous
            if( method == "VB"){
              return( HGF_continuous_VB(object))
            }
            if( method == "BF" ){
              return( HGF_continuous_BF(object))
            }

          }
)

# Define a method for plot specific to HGF_continuous
setMethod("plot",
          signature(object = "HGF_continuous"),
          function(object, ...) {
            # Implementation of plot method specific to HGF_continuous
            learner = object
            N = length(learner@u)

            plots = list()

            plots[[1]] = ggplot() +
              geom_line(aes(x=1:N,
                            y=learner@u))

            for( level in 1:learner@no_of_levels ){

              plot_df = data.frame(time=1:N,
                                   y=learner@moments$mu[,level],
                                   ymin=learner@moments$mu[,level]-learner@moments$sigma[,level]*1,
                                   ymax=learner@moments$mu[,level]+learner@moments$sigma[,level]*1)

              plots[[level+1]] = ggplot(plot_df) +
                geom_line(aes(x=time,
                              y=y)) +
                geom_ribbon(aes(x=1:N,
                                ymin=ymin,
                                ymax=ymax),
                            alpha=0.2)

            }

            plot_grid(plotlist=plots, ncol=1)
          }
)



setMethod("plot.distributions",
          signature(object = "HGF_continuous", timestamps="numeric", levels="numeric"),
          function(object, timestamps, levels=c(2,3), ...) {

            plots = list()
            datalist = list()

            for( t in timestamps ){
              for( level in levels ){

                sample = object@simulations[[t]][,paste0("x",level)]
                weights = object@simulations[[t]][,"w"]

                sample = sample[ weights > 0 ]
                weights = weights[ weights > 0]/sum(weights, na.rm = T)

                sample_limits = find_optimal_plot_limits_3(sample, weights)

                sample_min = sample_limits[1]
                sample_max = sample_limits[2]

                mu = object@moments$mu[t, level]
                sd = object@moments$sigma[t, level]
                vb_min = mu - 3*sqrt(sd)
                vb_max = mu + 3*sqrt(sd)

                xlim_min = min(c( sample_min,
                                  vb_min))
                xlim_max = max(c(sample_max,
                                 vb_max))

                x = seq(xlim_min, xlim_max, length.out=100)
                px = dnorm(x, mean = mu, sd = sqrt(sd))

                df1 = data.frame(x=sample, w=weights)
                df2 = data.frame(x=x,px=px)

                df1 = df1[apply(!is.na(df1), 1, any),]
                df1 = df1[df1$w!=0,]



                plots[[length(plots) + 1]] = ggplot() +
                  geom_density(aes(x=x, weight=w),
                               data = df1 ,
                               fill = "lightblue",
                               alpha=0.3) +
                  geom_line(aes(x=x, y=px),
                            data=df2) +
                  xlim(c(xlim_min,xlim_max)) +
                  ggtitle(paste( "level " , as.character(level),
                                 " time " , as.character(t)))

              }

            }

            final = plot_grid(plotlist=plots, ncol = length(levels))
            return(final)

          }
)



find_optimal_plot_limits = function(sample, threshold=0.005, breaks=100){

  sample_hist = hist(sample, plot = F, breaks = breaks)

  # 'density' sometimes fails for some reason
  densities = sample_hist$counts/sum(sample_hist$counts)

  max_dens = max(sample_hist$density)

  to_keep = which(sample_hist$density > max_dens*threshold)

  x_min = sample_hist$breaks[min(to_keep)]
  x_max = sample_hist$breaks[max(to_keep) + 1]

  return(c(x_min, x_max))

}


find_optimal_plot_limits_3 = function(obs_sample, obs_weights,
                                      cutoff=0.01,
                                      breaks=100,
                                      threshold=0.005){

  to_keep = !is.na(obs_sample) & !is.na(obs_weights)
  obs_sample = obs_sample[to_keep]
  obs_weights = obs_weights[to_keep]

  obs_sample = obs_sample[obs_weights != 0 ]
  obs_weights = obs_weights[obs_weights != 0 ]

  increasing = order(obs_sample)
  obs_sample = obs_sample[increasing]

  min_pos = round(length(obs_sample)*cutoff)
  max_pos = length(obs_sample) - min_pos + 1
  obs_sample = obs_sample[min_pos:max_pos]
  obs_weights = obs_weights[min_pos:max_pos]

  sample_hist = hist(obs_sample, plot = F, breaks = breaks)

  # 'density' sometimes fails for some reason
  densities = sample_hist$counts/sum(sample_hist$counts)

  max_dens = max(sample_hist$density)

  to_keep = which(sample_hist$density > max_dens*threshold)

  x_min = sample_hist$breaks[min(to_keep)]
  x_max = sample_hist$breaks[max(to_keep)]

  return(c(x_min, x_max))

}


find_optimal_plot_limits_2 = function(obs_sample, obs_weights=NA,
                                      threshold=0.005,
                                      kernel_width = NA,
                                      breaks = 100,
                                      order_obs_sample = T,
                                      cutoff_rate = 0.01){ # this argument has to caluclated by function


  if( !all(is.na(obs_weights)) ){

    if( length(obs_weights) != length(obs_sample)){
      errorCondition(paste0("length(obs_weights) == ",
                            length(obs_weights),
                            " != length(obs_sample) == ",
                            length(obs_sample)))
    }

  }else{
    obs_weights = rep(1, length(obs_sample)) # make more efficient later, this wastes time <---------------<<<<<
  }

  to_keep = (!is.na(obs_sample)) & (!is.na(obs_weights))
  obs_sample = obs_sample[to_keep  ]
  obs_weights = obs_weights[to_keep]


  if(order_obs_sample ){

    s_order = order(obs_sample)
    obs_sample = obs_sample[s_order]
    obs_weights = obs_weights[s_order]
  }

  if( is.na(kernel_width) ){

    cutoff_min = round(length(obs_sample)*cutoff_rate)
    cutoff_max = length(obs_sample) - cutoff_min

    obs_sample = obs_sample[cutoff_min:cutoff_max]
    obs_weights = obs_weights[cutoff_min:cutoff_max]

    kernel_width = ( obs_sample[length(obs_sample)] - obs_sample[1] )/breaks

  }

  wnon = matrix(0, nrow=length(obs_sample), ncol=1) # weighted number of neighbours
  for( i in 1:( length(obs_sample) - 1 ) ){
    j=i+1
    while( j <= length(obs_sample) && obs_sample[j] - obs_sample[i] < kernel_width ){

      wnon[i] = wnon[i] + obs_weights[j]
      wnon[j] = wnon[j] + obs_weights[i]

      j = j + 1
    }

  }

  max_wnon = max(wnon)

  i=1
  while( wnon[i] < max_wnon*threshold ){i=i+1}
  x_min = obs_sample[i]

  i = length(obs_sample)
  while( wnon[i] < max_wnon*threshold ){ i = i-1 }
  x_max = obs_sample[i]

  return(c(x_min, x_max))

}

set.seed(1234)

u = c( rnorm(100, sd = 1) + 1:100/10, rnorm(100, sd=0.2)+10 )

aaa = HGF_continuous(u=u)
plot(aaa)

aaa = fit(aaa, method="BF")
plot.distributions(aaa, timestamps = c(1,10,40,80), levels = c(1,2,3))



