#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(BiocManager)
options(repos = BiocManager::repositories())
library(RNASeqPower)
library(shinythemes)

# Define validate function
# returns a message if condition is true
fn_validate <- function(input,message) if(input) print(message)

# Define UI for application
ui <- fluidPage(

    # Add theme
    theme = shinytheme("slate") , 
    
    # Application title
    titlePanel("RNA-Seq | Power analysis"),
  
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("in_pa_est","Variable to estimate",choices=c("n","cv","effect","alpha","power"),selected=1,multiple=FALSE),
            uiOutput("ui_pa")),

        # Show a plot of the generated distribution
        mainPanel(
            verbatimTextOutput("out_pa")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$ui_pa <- renderUI({
        div(
            textInput("in_pa_depth","Sequencing depth",value=100),
            if(input$in_pa_est != "n")  textInput("in_pa_n","Sample size",value=12),
            if(input$in_pa_est != "cv")  textInput("in_pa_cv","Coefficient of variation",value=0.4),
            if(input$in_pa_est != "effect")  textInput("in_pa_effect","Effect",value=2),
            if(input$in_pa_est != "alpha")  textInput("in_pa_alpha","Alpha",value=0.05),
            if(input$in_pa_est != "power")  textInput("in_pa_power","Power",value=0.8)
        )
    })
    
    output$out_pa <- renderPrint({
        depth <- as.numeric(unlist(strsplit(gsub(" ","",input$in_pa_depth),",")))
        validate(fn_validate(any(is.na(depth)),"Sequencing depth must be a numeric."))
        
        if(input$in_pa_est != "n") {
            n <- as.numeric(unlist(strsplit(gsub(" ","",input$in_pa_n),",")))       
            validate(fn_validate(any(is.na(n)),"Sample size must be a numeric."))
        }
        
        if(input$in_pa_est != "cv") {
            cv <- as.numeric(unlist(strsplit(gsub(" ","",input$in_pa_cv),",")))
            validate(fn_validate(any(is.na(cv)),"Coefficient of variation must be a numeric."))
        }
        
        if(input$in_pa_est != "effect") {
            effect <- as.numeric(unlist(strsplit(gsub(" ","",input$in_pa_effect),",")))
            validate(fn_validate(any(is.na(effect)),"Effect must be a numeric."))
        }
        
        if(input$in_pa_est != "alpha")  {
            alpha <- as.numeric(unlist(strsplit(gsub(" ","",input$in_pa_alpha),",")))
            validate(fn_validate(any(is.na(alpha)),"Alpha must be a numeric."))
            validate(fn_validate(any(alpha>=1|alpha<=0),"Alpha must be a numeric between 0 and 1."))
        }
        
        if(input$in_pa_est != "power")  {
            power <- as.numeric(unlist(strsplit(gsub(" ","",input$in_pa_power),",")))
            validate(fn_validate(any(is.na(power)),"Power must be a numeric."))
            validate(fn_validate(any(power>=1|power<=0),"Power must be a numeric between 0 and 1."))
        }
        
        switch(input$in_pa_est,
               "n"=rnapower(depth=depth, cv=cv, effect=effect, alpha=alpha, power=power),
               "cv"=rnapower(depth=depth, n=n, effect=effect, alpha=alpha, power=power),
               "effect"=rnapower(depth=depth, cv=cv, n=n, alpha=alpha, power=power),
               "alpha"=rnapower(depth=depth, cv=cv, effect=effect, n=n, power=power),
               "power"=rnapower(depth=depth, cv=cv, effect=effect, alpha=alpha, n=n)
        )
        
    })
}



# Run the application 
shinyApp(ui = ui, server = server)
