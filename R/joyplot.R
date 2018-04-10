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