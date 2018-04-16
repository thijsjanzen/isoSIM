joyplot_frequencies <- function(frequencies,
                                time_points,
                                picked_ancestor = "ALL"
                                )
{
  vz <- subset(frequencies,
               frequencies$time %in% time_points)
  vz$ancestor <- as.factor(vz$ancestor)

  if(picked_ancestor == "ALL") {
    p1 <- ggplot2::ggplot(vz, ggplot2::aes(x = vz$location,
                         y = as.factor(vz$time),
                         height = vz$frequency,
                         fill = vz$ancestor)) +
            ggridges::geom_ridgeline(scale = 1.3) +
            ggplot2::ylab("Time")
    return(p1)
  } else {
    vy <- subset(vz, vz$ancestor == picked_ancestor)
    p1 <- ggplot2::ggplot(vy, ggplot2::aes(x = vy$location,
                         y = as.factor(vy$time),
                         height = vy$frequency)) +
            ggridges::geom_ridgeline(scale = 1.3,
                                     fill = "lightblue") +
            ggplot2::ylab("Time")
    return(p1)
  }
}

plot_start_end <- function(results,
                           picked_ancestor = "ALL") {

  a1 <- results$initial_frequency
  a2 <- results$final_frequency

  a1_m <- dplyr::mutate(a1, timepoint = "start")
  a2_m <- dplyr::mutate(a2, timepoint = "end")

  to_plot_m <- rbind(a1_m, a2_m)

  if(picked_ancestor == "ALL") {
    to_plot <- to_plot_m

    p1 <- ggplot2::ggplot(to_plot, ggplot2::aes(x = to_plot$location,
                              y = to_plot$frequency,
                              colour = to_plot$ancestor,
                              group = interaction(to_plot$ancestor,
                                                  to_plot$timepoint))) +
      ggplot2::geom_line(ggplot2::aes(lty = to_plot$timepoint))
  } else {

    to_plot <- dplyr::filter(to_plot_m,
                      ancestor == picked_ancestor)

    p1 <- ggplot2::ggplot(to_plot, ggplot2::aes(x = location,
                              y = to_plot$frequency,
                              colour = to_plot$ancestor,
                              group = interaction(to_plot$ancestor,
                                                  to_plot$timepoint))) +
      ggplot2::geom_line(ggplot2::aes(lty = to_plot$timepoint))
  }
  return(p1)
}