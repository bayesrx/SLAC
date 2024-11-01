#' Creates model matrix for neighbors. Designed to handle edge case where there
#' is only one type of point without breaking.
#' @param neighbor_tibble `neighbor_tib` from `GetNeighborTibble` function.
#' @import dplyr
#' @import tibble
NeighborModelMatrix = function(neighbor_tibble){
  mark_levels = unique(neighbor_tibble$marks)
  model_matrix = neighbor_tibble
  for(mark in mark_levels){
    cur_mark = paste0("mark_", mark)
    model_matrix[cur_mark] = (neighbor_tibble$marks == mark) %>% as.numeric()
  }
  model_matrix = model_matrix %>%
    select(-marks) %>%
    select(starts_with("mark_"), everything()) %>%
    as.matrix()
  return(model_matrix)
}
