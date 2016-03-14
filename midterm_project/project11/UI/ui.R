setwd("E:/Spencer/forecast")
library("shiny")

shinyUI(fluidPage(
  
  titlePanel(title=h2("Forecast Studio")),
  
    sidebarPanel(
      uiOutput("region"),
      uiOutput("plant"),
      uiOutput("material")
      ),
  
  
    mainPanel(
      tabsetPanel(type="tab",
                  tabPanel("Home",img(src="umich.png", height = 300, width = 900),img(src="supermrkt.jpg", height = 300, width = 900)),
                  tabPanel("Data",plotOutput("data")),
                  tabPanel("Holt winter's forecast",tableOutput("accuracy")),
                  tabPanel("Forecast",verbatimTextOutput("forecasted")),
                  tabPanel("Decomposition",plotOutput("decomposition")),
                  tabPanel("Strutural Time Series",tableOutput("stsm"),plotOutput("stplot"))
                  #tabPanel("Str Plot",plotOutput("stplot"))
                  
                  )
              )
)
)


