extractTrajectories <- function(variable, nSamples)
{
  output <- NULL
  samples <- sort(sample.int(999, nSamples))
  
  for(i in samples)
  {
    tmp <- rbi::extract_sample(posterior, i)
    tmp2 <- tmp[variable]
    tmp3 <- unlist(tmp2, recursive = FALSE, use.names=FALSE)
    output <- cbind(output, tmp3[[2]])
  }
  return(output)
}

extractTrajectoriesTibble <- function(variable, nSamples)
{

  outputDF <- data.frame(extractTrajectories(variable, nSamples)) %>%
    dplyr::tibble() %>%
    tibble::add_column(time = 1:32, .before = "X1") %>%
    tidyr::pivot_longer(cols = starts_with("X"), 
                        names_to = "realisation") %>% 
    dplyr::mutate(realisation = stringr::str_remove(realisation, "X")) %>%
    dplyr::group_by(realisation) %>%
    dplyr::arrange(realisation)
  
}
 
plottingFun <- function(variable, colorIndex, plotTitle, nSamples, withData = FALSE)
{
  
  trajectoryTibble <- extractTrajectoriesTibble(variable, nSamples)
  dataMedian <- dplyr::tibble(time = 1:32,
                           value = subset(trajectories, var == variable)$Median)
  
  plot <- ggplot2::ggplot(data = trajectoryTibble,
                          ggplot2::aes(x = time, 
                                       y = value)) +
  ggplot2::stat_smooth(ggplot2::aes(group = realisation),
                       geom ='line',
                       alpha = 0.015, 
                       se = FALSE,
                       color = viridis::viridis(200)[colorIndex]) + 
  ggplot2::stat_smooth(ggplot2::aes(x = time, y = value),
                       geom ='line',
                       data = dataMedian,
                       alpha = 0.8, 
                       se = TRUE,
                       color = viridis::viridis(200)[colorIndex]) +
  # ggplot2::geom_ribbon(data = trajectoryTibble, 
  #                      ggplot2::aes(ymin=`2.5%`, ymax= `97.5%`), 
  #                      fill = viridis::viridis(200)[colorIndex],
  #                      alpha = 0.1) +
  ggplot2::geom_vline(xintercept = 14, linetype = 4, color = "#ff6961", size = 0.5) + 
  ggplot2::ggtitle(plotTitle) + 
  ggplot2::theme(axis.title.x = element_blank(),
                 axis.title.y = element_blank())  +
  ggplot2::scale_x_continuous(breaks = breaksPlot, 
                              labels = datesPlot) +
  ggplot2::scale_y_continuous(limits = c(-1, NA)) +
  ggplot2::theme(plot.title = element_text(family = "Source Sans Pro",
                                            angle = 0, 
                                            hjust = 1, 
                                            size = 9)) +
  ggplot2::theme(axis.text.x = element_text(family = "Source Sans Pro",
                                            angle = 90, 
                                            hjust = 1, 
                                            size = 7)) +
  ggplot2::theme(axis.text.y = element_text(family = "Source Sans Pro",
                                            angle = 0, 
                                            hjust = 1, 
                                            size = 7))
  
  if(withData == TRUE & variable == "Z_sk")
  {
    plot <- plot + ggplot2::geom_point(data = dataTogether, 
                                       ggplot2::aes(x = time, y = symp),
                                       color = "black")
  }
  
  if(withData == TRUE & variable == "Z_n")
  {
    plot <- plot + ggplot2::geom_point(data = dataTogether, 
                                       ggplot2::aes(x = time, y = non_symp),
                                       color = "black")
  }
  
  return(plot)
}
