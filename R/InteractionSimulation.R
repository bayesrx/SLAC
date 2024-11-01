#' Function for simulating data with positive or negative interaction
#' @param n_each numeric vector where each entry is the number of cells to be
#' simulated for a given type of point
#' @param width number indicating the width of the desired window of simulation
#' @param height number indicating the height of the desired window of simulation
#' @param r_interaction number indicating the radius at which two points can be
#' said to be interacting
#' @param pos_m symmetric T x T matrix where T is the number of types of points.
#' Each entry should be between 0 and 1, and rows should sum to at most 1. Row i,
#' column j represents the probability of a point of type i *positively*
#' interacting with a point of type j, i.e. being generated within `r_interaction`
#' microns of a point of type j.
#' @param neg_m Defined analogously to `pos_m`, but for negative interactions.
#' @param existing TODO: take existing configuration as starting point
#' @param always_random_interaction If `TRUE`, a point is always chosen to
#' randomly interact positively or negatively with another point. If `FALSE`,
#' the type of interaction will only be selected randomly if the point interacts
#' both positively and negatively with points of other types.
#' @export
#' @import spatstat.geom
#' @import purrr
#' @import tibble
#' @import dplyr
InteractionSimulation = function(n_each, width, height, r_interaction,
                                 pos_m, neg_m, existing = NULL,
                                 always_random_interaction = FALSE){
  #------------------------------------------------------------
  # Preamble
  n_types = length(n_each)
  cur_type = 0
  pos_m = cbind(pos_m, 1 - rowSums(pos_m))
  neg_m = cbind(neg_m, 1 - rowSums(neg_m))
  full_window = spatstat.geom::owin(c(0, width), c(0, height))
  type_regions = map(1:n_types, ~ spatstat.geom::owin(c(0, 0), c(0, 0)))
  n_rejects = numeric(n_types)
  sim_data = tibble(
    x = numeric(),
    y = numeric(),
    type = numeric()
  )

  while(any(n_each > 0)){
    # if(n_each[cur_type] == 0){
    cur_type = `if`(cur_type + 1 > n_types, 1, cur_type + 1)
    while(n_each[cur_type] == 0){
      cur_type = `if`(cur_type + 1 > n_types, 1, cur_type + 1)
    }
    # }
    n_each[cur_type] = n_each[cur_type] - 1

    if(!always_random_interaction && pos_m[cur_type, n_types + 1] == 1){
      # No positive interactions for this type
      is_pos = FALSE
    }else if(!always_random_interaction && neg_m[cur_type, n_types + 1] == 1){
      # No negative interactions for this type
      is_pos = TRUE
    }else{
      # If it interacts positively and negatively, choose randomly
      is_pos = rbinom(1, 1, 0.5) %>% as.logical()
    }

    if(is_pos){
      cur_probs = pos_m[cur_type, ] %>% as.numeric()
    } else {
      cur_probs = neg_m[cur_type, ] %>% as.numeric()
    }

    which_interaction = sample(1:(n_types + 1), 1, prob = cur_probs)

    any_of_that_type = length(which(sim_data$type == which_interaction)) > 0
    if(which_interaction > n_types || !any_of_that_type){
      sim_data = bind_rows(sim_data, c(x = runif(1, 0, width),
                                       y = runif(1, 0, height),
                                       type = cur_type))
    }else if(is_pos){
      anchor_cell_index = sample(which(sim_data$type == which_interaction), 1)
      anchor_x = sim_data$x[anchor_cell_index]
      anchor_y = sim_data$y[anchor_cell_index]
      r = runif(1, 0, r_interaction)
      theta = runif(1, 0, 2*pi)
      sim_data = bind_rows(sim_data, c(x = r * cos(theta) + anchor_x,
                                       y = r * sin(theta) + anchor_y,
                                       type = cur_type))
    } else {
      neg_sim = tryCatch({
        rpoint(1, 1,
               win = spatstat.geom::intersect.owin(type_regions[[which_interaction]],
                                    full_window) %>%
                 spatstat.geom::complement.owin(frame = full_window))
      },
      error = function(e){
        if(e$message == "Gave up after 1000 trials, 0 points accepted"){
          n_rejects[cur_type] <<- n_rejects[cur_type] + 1
          return(rpoint(1, 1, win = full_window))
        } else {
          stop(e$message)
        }
      }
      )
      sim_data = bind_rows(sim_data, c(x = neg_sim$x[1],
                                       y = neg_sim$y[1],
                                       type = cur_type))
    }

    type_regions[[cur_type]] = spatstat.geom::union.owin(type_regions[[cur_type]],
                                          disc(radius = r_interaction,
                                               centre = c(
                                                 sim_data$x[nrow(sim_data)],
                                                 sim_data$y[nrow(sim_data)]
                                               )))

  }

  for(i in 1:length(n_rejects)){
    rej_count = n_rejects[i]
    if(rej_count > 0){
      warning(paste(rej_count, "points of type", i, "could not be",
                    "simulated with negative interaction"))
    }
  }

  return(sim_data)
}
