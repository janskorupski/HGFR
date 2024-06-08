library(R.oo)
library(methods)
library(ggplot2)
library(cowplot)

s = function(x){
  return( 1/(1 + exp(-x)) )
}

setClass("HGF",
         slots = c(

           no_of_levels    = "numeric",
           no_of_moments   = "numeric",
           parameters      = "list",
           priors          = "list",
           update_formulas = "list",
           u               = "numeric",
           moments         = "list",
           simulations     = "list"

         ),
         prototype = list(

           no_of_levels    = NA_integer_,
           no_of_moments   = NA_integer_,
           parameters      = list(),
           priors          = list(),
           update_formulas = list(),
           u               = NA_real_,
           moments         = list(),
           simulations     = list()

         )
)


HGF <- function(no_of_levels = integer(),
                        no_of_moments = integer(),
                        parameters = list(),
                        priors = list(),
                        update_formulas = list(),
                        u = integer()){

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

  new("HGF",
            no_of_levels = no_of_levels,
            no_of_moments = no_of_moments,
            parameters = parameters,
            priors = priors,
            update_formulas = update_formulas,
            u = u,
            moments = moments,
            simulations = list())

}

# Define the fit generic function
setGeneric("fit",
           function(object, method, ...)
             standardGeneric("fit")
)

# Define a method for fit
setMethod("fit",
          signature(object = "HGF", method = "character"),
          function(object, method, ...) {
            # Implementation of fit method
            # This will be implemented differently based on the method chosen (e.g., "variational", "mcmc")
            stop("fit method not implemented")
          }
)

# Define a plot method
setGeneric("plot",
           function(object, ...)
             standardGeneric("plot")
)

# Define a method for plot
setMethod("plot",
          signature(object = "HGF"),
          function(object, ...) {
            # Implementation of plot method
            stop("plot method not implemented")
          }
)


setGeneric("plot.distributions",
           function(object, timestamps=NaN, levels=NaN,...) standardGeneric("plot.distributions"),
)



