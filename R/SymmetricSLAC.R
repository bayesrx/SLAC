#' Implementation of SLAC model used in simulation and application of data in
#' "SLAC: A Data Augmentation Approach to Modeling Multiplex Imaging Data."
#' @param mltplx_object object of class MltplxObject (from DIMPLE package).
#' @param r radius of interaction between different points.
#' @param width width of observation window of observed point process.
#' @param height height of observation window of observed point process.
#' @param n_aux_each number of auxiliary points per cell type. Can be either a
#' single scalar number, or a named vector with length equal to the number of
#' cell types in `mltplx_object`, and names corresponding to the types of cells in
#' `mltplx_object`.
#' @param fit_method method by which glm can be fit. "glm" will use the `glm`
#' function, "glmnet" will use `cv.glmnet`, and "spikeslab" will use a spike-and-
#' slab formulation implemented by the `BoomSpikeSlab` package.
#' @param n_iter number of iterations to run the spike-and-slab model for. Must
#' be provided when `fit_method` is set to "spikeslab."
#' @param saturation optional. This term can be used to specify a maximum number
#' of neighbors that a point can have, in the manner of Rajala et al. 2019. Its
#' usage is generally not recommended unless you know what you are doing.
#' @export
SymmetricSLAC = function(mltplx_object,
                           r, width, height,
                           n_aux_each,
                           fit_method = c("glm", "glmnet", "spikeslab"),
                           n_iter = NULL,
                           saturation = Inf){
  #------------------------------------------------------------
  # Preliminary computations
  n_cell_types = mltplx_object$mltplx_image$ppp$marks %>% unique %>% length
  print("Generating auxiliary points...")
  neighbor_result = GetNeighborTibble(mltplx_object, r, width, height,
                                saturation, n_aux_each)
  neighbors = neighbor_result$neighbor_tib

  #------------------------------------------------------------
  # Model fitting
  if(length(fit_method) > 1){
    stop(paste0("fit method must be `glm`, `glmnet`, or `spikeslab`"))
  }

  if(fit_method == "glm"){
    print("Fitting glm...")
    model_fit = glm(y ~ 0+.,
                    data = neighbors %>%
                      NeighborModelMatrix() %>%
                      as_tibble(),
                    family = binomial())
  } else if(fit_method == "glmnet"){
    print("Fitting glmnet...")
    x = neighbors %>% select(-y) %>% NeighborModelMatrix()
    model_fit = glmnet::cv.glmnet(x = x,
                                  y = neighbors$y,
                                  family = "binomial",
                                  penalty.factor = c(rep(0, n_cell_types),
                                                     rep(1, ncol(x) - n_cell_types))
    )
  } else if(fit_method == "spikeslab"){
    method_settings = list()
    neighbor_model_matrix = neighbors %>%
          select(-y) %>%
          NeighborModelMatrix()

    boom_prior = BoomSpikeSlab::LogitZellnerPrior(
        predictors = neighbor_model_matrix,
        successes = neighbors$y,
        prior.inclusion.probabilities = c(
          rep(1, n_cell_types),
          rep(0.5, ncol(neighbor_model_matrix) - n_cell_types)
        )
      )

    method_settings$formula = y ~ 0+.

    method_settings$data = neighbors %>%
                            NeighborModelMatrix() %>%
                            as_tibble()
    method_settings$prior = boom_prior
    if(is.null(n_iter)){
      stop("`n_iter` must be provided when fitting the spike-and-slab model")
    }
    method_settings$niter = n_iter

    print("Fitting spike and slab...")
    model_fit = do.call(BoomSpikeSlab::logit.spike, method_settings)
  } else {
    stop(paste0("fit method must be `glm`, `glmnet`, or `spikeslab`"))
  }

  results_object = list(
    fit = model_fit,
    raw_points = neighbor_result$points,
    raw_data = neighbors %>%
                  NeighborModelMatrix() %>%
                  as_tibble()
  )

  return(results_object)

}
