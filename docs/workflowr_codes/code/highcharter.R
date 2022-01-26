# Line chart


plot(y~x ,
     lwd=4 , type="l" , bty="n" , ylab="value of y (decreasing)" , col=rgb(0.2,0.4,0.6,0.8) )

ggplot( aes(x=x, y=y)) +
geom_line() +
geom_point()library(tidyverse)

library(hrbrthemes)
library(plotly)
library(patchwork)
library(babynames)
library(viridis)

data <- data.frame(x,y)

data %>%
  ggplot( aes(x=x, y=y, label=y)) +
  geom_path()+
  geom_point(shape=21, color="black", fill="#69b3a2", size=6) +
  theme_gray()+
  ggtitle("Evolution of bitcoin price")

library(highcharter)
x <- c("v1.0.0 (2019-11-13)", "v1.2.2 (2020-12-10)", "v1.3.0 (2021-01-04)", "v1.4.0 (2021-03-09)", "v1.5.0 (2021-06-22)")
y <- c(15, 242, 278, 448, 457)
x_annot <- c("Initial Release", "Integrate LISA2", "Support multi-sample scATAC-seq", "Add Signic", "Add Chromap and new documentations" )

highchart() %>%
  hc_chart(type = "line") %>%
  hc_title(text = "MAESTRO Anaconda User Statistics") %>%
  hc_subtitle(text = "Source: https://anaconda.org/liulab-dfci/maestro/files ") %>%
  hc_xAxis(categories = x, title= list(text = "Released versions")) %>%
  hc_yAxis(title = list(text = "Cumulative number of downloads")) %>%
  hc_plotOptions(line = list(
    dataLabels = list(enabled = TRUE),
    enableMouseTracking = FALSE)
  ) %>%
  hc_series(
    list(
      name = "MAESTRO User Growth",
      data = y
    )
  ) %>% hc_annotations(
    list(
      labels =
        list(
          list(
            point = list(xAxis = "v1.0.0 (2019-11-13)", y= 15, yAxis = 0),
            text ="Initial Release"
            ),
          list(
            point = list(xAxis = "v1.2.2 (2020-12-10)", y= 242, yAxis = 0),
            text = "Integrate LISA2"
            ),
          list(
            point = list(xAxis = "v1.3.0 (2021-01-04)", y= 278, yAxis = 0),
            text = "Support multi-sample scATAC-seq"
            ),
          list(
            point = list(xAxis = "v1.4.0 (2021-03-09)", y= 448, yAxis = 0),
            text = "Add Signic"
            ),
          list(
            point = list(xAxis = "v1.5.0 (2021-06-22)", y= 457, yAxis = 0),
            text = "Add Chromap and new documentations"
            )
          )
      )
  )



highchart() %>%
  hc_add_series(
    data = c(29.9, 71.5, 106.4, 129.2, 144.0, 176.0, 135.6, 148.5, 216.4, 194.1, 95.6, 54.4)
  ) %>%
  hc_xAxis(
    tickInterval = 0.5,
    gridLineWidth = 1
  ) %>%
  hc_annotations(
    list(
      labels =
        list(
          list(
            point = list(x = 3, y = 129.2, xAxis = 0, yAxis = 0),
            text = "x: {x}<br/>y: {y}"
          ),
          list(
            point = list(x = 9, y = 194.1, xAxis = 0, yAxis = 0),
            text = "x: {x}<br/>y: {y}"
          ),
          list(
            point = list(x = 5, y = 100, xAxis = 0),
            text = "x: {x}<br/>y: {point.plotY} px"
          ),
          list(
            point = list(x = 0, y = 0),
            text = "x: {point.plotX} px<br/>y: {point.plotY} px"
          )
        )
    )
  )
