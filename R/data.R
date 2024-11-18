#' Example Data From a Simulation Study
#'
#' A simulated dataset used to demonstrate mediator screening and analysis functions.
#'
#' @format ## `ExampleData`
#' A list containing the following elements:
#' \describe{
#'   \item{M}{A 200x2000 matrix of mediators. Each row corresponds to an observation, and each column represents a potential mediator.}
#'   \item{y}{A vector of length 200 representing the outcome for each observation.}
#'   \item{x}{A vector of length 200 representing the exposure for each observation.}
#'   \item{alp_vec}{A parameter vector of length 200 that relates the exposure variable to the mediators}
#'   \item{beta_vec}{A parameter vector of length 200 that relates the mediators to the outcome variable}
#' }
#' @source Simulated data for illustrative purposes.
"ExampleData"
