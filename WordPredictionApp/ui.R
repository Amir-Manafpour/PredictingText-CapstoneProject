#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("3-gram Word Prediction App"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(position = "left",
        sidebarPanel(
            textInput("sentInput",
                      "Enter text below:")
        ),

        # Show a plot of the generated distribution
        mainPanel(
            h2("Prediction Results"),
            DT::dataTableOutput("predResults")
        )
    )
))
