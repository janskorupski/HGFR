
# Define the HGF_binary subclass inheriting from HGF
setClass("HGF_binary",
         contains = "HGF",
         slots = c(

           no_of_levels    = "numeric",
           no_of_moments   = "numeric",
           parameters      = "list",
           priors          = "list",
           update_formulas = "list",
           u               = "numeric",
           moments         = "list",
           simulations     = "list",
           sim_preds_      = "numeric"

         ),
         prototype = list(

           no_of_levels    = 3,
           no_of_moments   = 2,
           parameters      = list("kappa"=1.4, "omega"=-2.2, "theta"=0.5),
           priors          = list("mu" = c(0.5, 0, 0), "sigma" = c(NaN, 1, 1)),
           update_formulas = list(),
           u               = NA_real_,
           moments         = list(),
           simulations     = list(),
           sim_preds_      = numeric()

         )
)

# Constructor function for HGF_binary objects
HGF_binary <- function(
                    u = integer(),
                    no_of_levels = 3,
                    no_of_moments = 2,
                    parameters = list("kappa"=1.4, "omega"=-2.2, "theta"=0.5),
                    priors = list("mu" = c(0.5, 0, 0), "sigma" = c(NaN, 1, 1)),
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

  object = new("HGF_binary",
      no_of_levels = no_of_levels,
      no_of_moments = no_of_moments,
      parameters = parameters,
      priors = priors,
      update_formulas = update_formulas,
      u = u,
      moments = moments,
      simulations = list(),
      sim_preds_ = numeric())

  if( length(object@u) > 1 ){
    object = fit(object, method="VB")
  }

  return(object)

}

mathys_VB = function(learner){


  if( ! is.nan(learner@u[1]) ){
    learner@u = c(NaN, learner@u)
  }

  for( moment_no in 1:learner@no_of_moments){
    learner@moments[[moment_no]] =
      matrix(nrow = (length(learner@u)),
             ncol = learner@no_of_levels)
    learner@moments[[moment_no]][1,] = learner@priors[[moment_no]]
  }



  kappa = learner@parameters$kappa
  omega = learner@parameters$omega
  theta = learner@parameters$theta


  for( t in 2:(length(learner@u)) ){

    # LEVEL 1
    learner@moments$mu[t, 1] = learner@u[t]
    if( is.na(learner@u[t]) ){ # TWEAK THIS LATER -------------------------------<<<
      print("dupa")
      learner@moments[[1]][t, 1] = learner@moments[[1]][t-1]
    }


    # LEVEL 2
    mu_hat_1 = s( learner@moments$mu[t-1, 2] )
    delta_hat_1 = learner@moments$mu[t, 1] - mu_hat_1

    ### CHECK THAT LATER, IT'S DIFFERENT IN THE PAPER !!!
    sigma_hat_1 = mu_hat_1*(1 - mu_hat_1)

    sigma_hat_2 = learner@moments$sigma[t-1, 2] +
      exp( kappa * learner@moments$mu[t-1, 3] + omega )

    learner@moments$sigma[t, 2] = 1/( 1/sigma_hat_2 + sigma_hat_1)
    learner@moments$mu[t, 2] = learner@moments$mu[t-1, 2] + learner@moments$sigma[t, 2] * delta_hat_1

    # LEVEL 3
    exp_kmo = exp(kappa * learner@moments$mu[t-1, 3] + omega)
    w_2 = (exp_kmo)/
      ( learner@moments$sigma[t-1, 2] + exp_kmo )
    r_2 = (exp_kmo - learner@moments$sigma[t-1, 2]) /
      (exp_kmo + learner@moments$sigma[t-1, 2])
    delta_2 = (learner@moments$sigma[t, 2] + (learner@moments$mu[t, 2] - learner@moments$mu[t-1, 2])**2)/
      (learner@moments$sigma[t-1, 2] + exp_kmo) - 1

    pi_hat_3 = 1/( learner@moments$sigma[t-1, 3] + learner@parameters$theta )
    pi_3 = pi_hat_3 + (learner@parameters$kappa**2)/2 * w_2 * (w_2 + r_2 * delta_2)

    learner@moments$sigma[t, 3] = 1/pi_3

    learner@moments$mu[t, 3] = learner@moments$mu[t-1, 3] + learner@moments$sigma[t, 3] * learner@parameters$kappa/2 *w_2 * delta_2
  }
  return( learner )
}

mathys_RS = function(learner){
  N = 1e4

  if( ! is.nan(learner@u[1]) ){
    learner@u = c(NaN, learner@u)
  }


  kappa = learner@parameters$kappa
  omega = learner@parameters$omega
  theta = learner@parameters$theta

  ## First time based on priors
  x3 = rnorm(N,
             mean = learner@priors$mu[3],
             sd = sqrt(learner@priors$sigma[3]))
  x2 = rnorm(N,
             mean = learner@priors$mu[2],
             sd = sqrt(learner@priors$sigma[2]))
  x1 = rbinom(N,
              size = 1,
              prob = learner@priors$mu[1])
  learner@simulations[[1]] = data.frame(x1=x1, x2=x2, x3=x3)
  learner@sim_preds_[1] = learner@priors$mu[1]


  for( t in 2:(length(learner@u)) ){

    res = matrix(nrow = 0, ncol = learner@no_of_levels)

    prob_of_ones=c()

    while(dim(res)[1] < N){
      m = N - dim(res)[1]

      prev_id = sample(1:N, m, replace = T)

      x3_prev = learner@simulations[[t-1]]$x3[prev_id]
      x3 = rnorm(m,
                 mean = x3_prev,
                 sd = sqrt(theta))
      sd_2 = exp(kappa * x3 + omega)
      x2_prev = learner@simulations[[t-1]]$x2[prev_id]
      x2 = rnorm(m,
                 mean = x2_prev,
                 sd = sqrt(sd_2))
      x1 = rbinom(m,
                  size = 1,
                  prob = s(x2))

      prob_of_ones = c(prob_of_ones, s(x2))

      add = matrix(c(x1,x2,x3), ncol = learner@no_of_levels)
      add = add[ add[,1] == learner@u[t] , ]
      res = rbind(res, add)

    }

    learner@simulations[[t]] = as.data.frame(res)
    learner@sim_preds_[t] = mean(prob_of_ones)
    names(learner@simulations[[t]]) = c("x1", "x2", "x3")

  }
  return( learner )
}

# Define a method for fit specific to HGF_binary
setMethod("fit",
          signature(object = "HGF_binary", method = "character"),
          function(object, method="VB", ...) {
            # Implementation of fit method specific to HGF_binary
            if( method == "VB"){
              return( mathys_VB(object))
            }
            if( method == "RS" ){
              return( mathys_RS(object))
            }

          }
)

# Define a method for plot specific to HGF_binary
setMethod("plot",
          signature(object = "HGF_binary"),
          function(object, ...) {
            # Implementation of plot method specific to HGF_binary
            learner = object
            N = length(learner@u)

            lvl1 = ggplot() +
              geom_point(aes(x=1:N,
                             y=learner@moments$mu[,1])) +
              geom_line(aes(x=1:N,
                            y=s( learner@moments$mu[, 2] )))


            lvl2 = ggplot() +
              geom_line(aes(x=1:N,
                            y=learner@moments$mu[,2])) +
              geom_ribbon(aes(x=1:N,
                              ymin=learner@moments$mu[,2]-learner@moments$sigma[,2]*1,
                              ymax=learner@moments$mu[,2]+learner@moments$sigma[,2]*1),
                          alpha=0.2)

            lvl3 = ggplot() +
              geom_line(aes(x=1:N,
                            y=learner@moments$mu[,3])) +
              geom_ribbon(aes(x=1:N,
                              ymin=learner@moments$mu[,3]-learner@moments$sigma[,3]*1,
                              ymax=learner@moments$mu[,3]+learner@moments$sigma[,3]*1),
                          alpha=0.2)



            plot_grid(lvl3, lvl2, lvl1, ncol=1)
          }
)


# Define a method for fit specific to HGF_binary
setMethod("plot.distributions",
          signature(object = "HGF_binary", timestamps="numeric", levels="numeric"),
          function(object, timestamps=NaN, levels=c(2,3), ...) {

            plots = list()
            datalist = list()

            for( t in timestamps ){
              for( level in levels ){


                sample = object@simulations[[t]][,level]
                sample = sample[ abs(sample - mean(sample)) < 2.5*sd(sample) ] # < -----------------------------<<<<< to fix
                sample_min = min(sample)
                sample_max = max(sample)

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

                df1 = data.frame(x=sample)
                df2 = data.frame(x=x,px=px)

                plots[[length(plots) + 1]] = ggplot() +
                  geom_density(aes(x=x),
                               data = df1 ,
                               fill = "lightblue",
                               alpha=0.3) +
                  geom_line(aes(x=x, y=px),
                            data=df2) +
                  ggtitle(paste( "level " , as.character(level),
                                 " time " , as.character(t)))

              }

            }

            final = plot_grid(plotlist=plots, ncol = length(levels))
            return(final)

          }
)

