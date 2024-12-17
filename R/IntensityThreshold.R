#' @import spatstat.explore
IntensityThreshold = function(mltplx_obj, ps, bw, threshold) {
  point_intensity = spatstat.explore::density.ppp(
    mltplx_obj$mltplx_image$ppp,
    sigma = bw,
    eps = ps,
    diggle = TRUE
  )

  points_x = mltplx_obj$mltplx_image$ppp$x
  points_y = mltplx_obj$mltplx_image$ppp$y

  x_min = point_intensity$xrange[1]
  x_step = point_intensity$xstep
  pixel_x = point_intensity$xcol
  y_min = point_intensity$yrange[1]
  y_step = point_intensity$ystep
  pixel_y = point_intensity$yrow

  final_window = owin(c(0, 0), c(0, 0))

  for (c in 1:ncol(point_intensity$v)){
    for (r in 1:nrow(point_intensity$v)){
      cur_window_x_coord = (c - 1) * x_step + x_min
      cur_window_y_coord = (r - 1) * y_step + y_min

      # add/subtract 1 to create overlap and remove border artefacts
      cur_window = owin(c(cur_window_x_coord - 1, cur_window_x_coord + x_step + 1),
                        c(cur_window_y_coord - 1, cur_window_y_coord + y_step + 1))

      points_x_within = points_x > cur_window_x_coord &
                        points_x <= cur_window_x_coord + x_step
      points_y_within = points_y > cur_window_y_coord &
                        points_y <= cur_window_y_coord + y_step
      points_within = any(points_x_within & points_y_within)


      if(point_intensity$v[r, c] > threshold || points_within) {
        final_window = union.owin(final_window, cur_window)
      }

    }
  }

  final_ppp = ppp(
    x = mltplx_obj$mltplx_image$ppp$x,
    y = mltplx_obj$mltplx_image$ppp$y,
    marks = mltplx_obj$mltplx_image$ppp$marks,
    window = final_window
  )

  return(final_ppp)
}
