#' Given a ppp object, creates a tibble of points
#' and dummy points with neighbor counts. For internal use by `SymmetricSlac`
#' function.
#' @import dplyr
#' @import magrittr
#' @import tibble
#' @import purrr
GetNeighborTibble = function(obj, r, saturation, n_aux_each){
  if(class(obj)[1] != "ppp"){
    stop(paste0("`obj` argument to `GetNeighborTibble` should be `ppp` object"))
  }
  n = obj$n

  #------------------------------------------------------------
  # Marks
  actual_marks = obj$marks %>% as.character()
  unique_marks = actual_marks %>% unique() %>% sort()
  mark_pairs = tidyr::crossing(m1 = unique_marks, m2 = unique_marks) %>%
    filter(
      m1 <= m2
    )

  #------------------------------------------------------------
  # Checking dummy labels
  if(length(n_aux_each) == 1){
    n_aux_each = rep(n_aux_each, length(unique_marks))
    names(n_aux_each) = unique_marks
  }

  dummy_marks = names(n_aux_each)
  dummy_counts = unname(n_aux_each)

  if(length(dummy_marks) != length(unique_marks)){
    stop("length of dummy counts does not match length of unique marks")
  }

  if(! all( dummy_marks %in% unique_marks )){
    stop("not all dummy marks in `n_aux_each` are present in data")
  }

  if(! all( unique_marks %in% dummy_marks )){
    stop("not all marks in data are present in `n_aux_each`")
  }

  #------------------------------------------------------------
  # Dummy points
  dummy_sim = runifpoint(sum(n_aux_each), win = obj$window)
  dummy_x = dummy_sim$x
  dummy_y = dummy_sim$y
  dummy_coords = matrix(c(dummy_x, dummy_y), ncol = 2)

  #------------------------------------------------------------
  # Regular points
  coords_mat = matrix(c(obj$x, obj$y), ncol = 2)

  #------------------------------------------------------------
  # Full matrix
  full_coords = rbind(coords_mat, dummy_coords)
  full_dist = dist(full_coords) %>% as.matrix()
  diag(full_dist) = Inf
  full_dist = full_dist[,1:n]

  #------------------------------------------------------------
  # Tibble
  neighbor_tib = tibble(
    y = c(rep(1, n), rep(0, sum(n_aux_each))),
    marks = c(actual_marks, rep(dummy_marks, dummy_counts))
  )

  for(m in unique_marks){
    reduced_dists = full_dist[,which(actual_marks == m)]
    if(reduced_dists %>% dim %>% is.null){
      reduced_dists = reduced_dists %>% as.matrix() %>% t()
    }
    neighbor_tib[m] = (reduced_dists < r) %>%
      rowSums %>%
      map_dbl(~ min(.x, saturation))
  }

  for(i in 1:length(unique_marks)){
    m1 = unique_marks[i]
    for(j in i:length(unique_marks)){
      m2 = unique_marks[j]
      mark_str = paste0(m1, ":", m2)
      m1_selection = neighbor_tib$marks == m1
      m2_selection = neighbor_tib$marks == m2
      if(m1 != m2){
        neighbor_tib[mark_str] = (neighbor_tib[m2]*m1_selection) +
          (neighbor_tib[m1]*m2_selection)
      } else {
        neighbor_tib[mark_str] = (neighbor_tib[m2]*m1_selection)
      }
    }
  }

  results_obj = list(
    neighbor_tib = neighbor_tib %>%
                    select(y, marks, contains(":")),
    points = full_coords %>%
      magrittr::set_colnames(c("x", "y")) %>%
      as_tibble() %>%
      mutate(
        marks = neighbor_tib$marks,
        observed = neighbor_tib$y
      )
  )

  return(results_obj)
}
