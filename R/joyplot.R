joyplot_frequencies <- function(frequencies,
                                time_points,
                                picked_ancestor = "ALL"
                                )
{
  vz <- subset(frequencies,
               frequencies$time %in% time_points)
  vz$ancestor <- as.factor(vz$ancestor)

  if(picked_ancestor == "ALL") {
    p1 <- ggplot2::ggplot(vz, aes(x = location,
                         y = as.factor(time),
                         height = frequency,
                         fill = ancestor)) +
            ggridges::geom_ridgeline(scale = 1.3) +
            ggplot2::ylab("Time")
    return(p1)
  } else {
    vy <- subset(vz, vz$ancestor == picked_ancestor)
    p1 <- ggplot2::ggplot(vy, aes(x = location,
                         y = as.factor(time),
                         height = frequency),) +
            ggridges::geom_ridgeline(scale = 1.3,
                                     fill = "lightblue") +
            ggplot2::ylab("Time")
    return(p1)
  }
}