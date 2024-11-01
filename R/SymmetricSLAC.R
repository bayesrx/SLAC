#' Implementation of SLAC model used in simulation and application of
#' @export
SymmetricSLAC = function(mltplx_object,
                           r, width, height,
                           n_dummy_each,
                           fit_method = c("glm", "glmnet", "spikeslab"),
                           n_iter = NULL,
                           saturation = Inf){
  #------------------------------------------------------------
  # Preliminary computations
  n_cell_types = mltplx_object$mltplx_image$ppp$marks %>% unique %>% length
  print("Generating dummy points...")
  neighbor_result = GetNeighborTibble(mltplx_object, r, width, height,
                                saturation, n_dummy_each)
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
