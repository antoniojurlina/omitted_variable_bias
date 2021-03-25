library(shiny)
library(shinyWidgets)
library(tidyverse)
library(MASS)
library(ggthemes)

set.seed(5678)

betas <- matrix(c(-3, 2, 0.5), 
                nrow = 3, 
                ncol = 1)

alphas <-  matrix(c(1.5, -3), 
                  nrow = 2, 
                  ncol = 1)

sigmas <- matrix(c(1, 0.75), 
                 nrow = 2, 
                 ncol = 1)

myOLS <- function(data, outcomes) {
    Y <- outcomes
    X <- data
    
    degrees_of_freedom <- nrow(X)-ncol(X)
    
    beta_estimates <- ginv(t(X)%*%X)%*%t(X)%*%Y
    
    error_estimates   <- Y - X%*%beta_estimates
    sigma_estimate    <- (t(error_estimates)%*%error_estimates) / degrees_of_freedom
    variance_estimate <- ginv(t(X)%*%X)
    standard_error    <- as.matrix(sqrt(diag(variance_estimate)))
    
    RSS       <- t(error_estimates)%*%error_estimates         
    y_bar     <- mean(Y)
    yybar     <- Y-y_bar
    TSS       <- t(yybar)%*%yybar 
    R_squared <- 1 - RSS/TSS 
    
    t_statistic <- beta_estimates / standard_error
    p_values     <- 2*(1-pt(abs(t_statistic), degrees_of_freedom)) %>%
        round(5)
    
    estimates <- tibble("estimate"    = c(beta_estimates), 
                        "std_error"   = c(standard_error), 
                        "t_statistic" = c(t_statistic), 
                        "p_value"     = c(p_values))
    
    list("estimates" = estimates,
         "r_squared" = R_squared)
}
clean_results <- function(output, omit = TRUE) {
    results <- c()
    
    for(i in 1:dim(output)[3]) {
        results <- rbind(results, output[,,i])
    }
    
    if(omit == FALSE) {
        beta_0 <- results[seq(1,(dim(output)[3]*2), 3),]
        beta_1 <- results[seq(2,(dim(output)[3]*2), 3),]
    }
    
    if(omit == TRUE) {
        beta_0 <- results[seq(1,(dim(output)[3]*2), 2),]
        beta_1 <- results[seq(2,(dim(output)[3]*2), 2),]
    }
    
    return(list("beta_0" = as_tibble(beta_0), 
                "beta_1" = as_tibble(beta_1)))
}

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    theme = shinythemes::shinytheme("journal"),

    # Application title
    titlePanel("Omitted Variable Bias"),

    sidebarLayout(
        sidebarPanel(
            progressBar(id = "pb1", 
                        value = 50, 
                        display_pct = FALSE, 
                        status = "warning"),
            
            sliderInput("observations", label = "Select number of observations", min = 0, 
                        max = 1000, value = 750),
            
            progressBar(id = "pb2", 
                        value = 50, 
                        display_pct = FALSE, 
                        status = "warning"),
            
            sliderInput("repetitions", label = "Select number of repetitions", min = 0, 
                        max = 1000, value = 750),
            
            awesomeCheckbox(
                inputId = "omit",
                label = "Omit variable", 
                value = FALSE
            ),
            
            uiOutput('equation1'),
            uiOutput('equation2'),
            
            br(),
            
            actionBttn(
                inputId = "calculate",
                label = "Calculate",
                style = "float", 
                color = "danger"
            )
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel("Introduction", value = "intro",
                         "Consider the linear model",
                         "$$y_i = \\beta_0 + \\beta_1 X_i + \\beta_2 D_i + 
                         \\epsilon_{1i}, \\; \\epsilon_{1} \\sim N(0,\\sigma_{1}^2)$$",
                         "where",
                         "$$D_i = \\alpha_0 + \\alpha_1 X_i + \\epsilon_{2i},
                         \\; \\epsilon_{2} \\sim N(0,\\sigma_{2}^2)$$",
                         "and",
                         "$$X_i \\sim N(10,5).$$",
                         "Effects on the dependent variable which are not captured by 
                         the model are collected in the error term, which we so far 
                         assumed to be uncorrelated with the regressor. This is one of
                         the key assumptions regarding OLS, namely that",
                         br(),
                         "$$E(e_i|X_i)=0.$$",
                         br(),
                         "However, omitting an independent variable that does have an
                         effect on the dependent one, while also being correlated 
                         with another regressor, would introduce bias. Specifically,
                         if two variables are correlated, and they both affect the outcome,
                         leaving one of them out would cause it to end up in the error 
                         term. Furthermore, the remaining variable would now be correlated
                         with the error term and the sampling distribution of its estimate
                         would not be dispersed around the true mean anymore. This effect 
                         would not get fixed with a larger sample size, either.",
                         br(), br(),
                         "In order to demonstrate this effect, I will generate psuedo-data,
                         in such a way that one of the two variables is correlated with both 
                         the outcome and the remaining variable. The pseudo data will be 
                         generated given the equations above, with",
                         br(),
                         "$$\\beta_0 = -3, \\beta_1 = 2, \\beta_2 = 0.5, \\sigma_1 = 1$$",
                         br(),
                         "$$\\alpha_0 = 1.5, \\alpha_1 = -3, \\sigma_2 = 0.75.$$",
                         br(),
                         "Omitted variable bias is the bias in the OLS estimator that arises 
                         when the regressor, ", tags$b("X"), "is correlated with an omitted variable. 
                         For omitted variable bias to occur, two conditions must be fulfilled:",
                         br(), br(),
                         strong("1. X is correlated with the omitted variable."),
                         br(),
                         strong("2. The omitted variable is a determinant of the dependent 
                                variable Y."),
                         br(), br(),
                         "Omitting the variable D will satisfy these
                         conditions. Resulting plots and data serve to verify our expectations. 
                         The user selects the number of observations for the generated data, as well
                         as the number of Monte Carlo repetitions for obtaining the sampling 
                         distribution of the parameter estimates.",
                         br(), br(),
                         "So, ", tags$b("CALCULATE"), " and head on over to Plots tab.",
                         HTML("<br><br><br>")
                         ),
                tabPanel("Plots", value = "plots",
                         plotly::plotlyOutput("histogram"),
                         HTML("<br><br>"),
                         strong("Estimates compared to the actual values:"),
                         plotly::plotlyOutput("histogram2"),
                         plotly::plotlyOutput("histogram3")),
                tabPanel("Table", 
                         "Average of estimates across all repetitions:",
                         tableOutput("estimates_table"),
                         "Individual estimates:",
                         tableOutput("complete_data")),
                tabPanel("Code", 
                         "Custom code for Ordinary Least Squares in R:",
                         verbatimTextOutput("codeR"),
                         "Custom code for Ordinary Least Squares in MATLAB:",
                         verbatimTextOutput("codeMATLAB"))
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    color1 <- reactive({
        if (input$observations < 250) {
            return("danger")
        } else if (input$observations >= 250 & input$observations < 750) {
            return("warning")
        } else {
            return("success")
        }
    })
    
    color2 <- reactive({
        if (input$repetitions < 250) {
            return("danger")
        } else if (input$repetitions >= 250 & input$repetitions < 750) {
            return("warning")
        } else {
            return("success")
        }
    })
    
    observeEvent(input$observations, {
        updateProgressBar(session = session, 
                          id = "pb1", 
                          value = input$observations/10, 
                          status = color1())
    })
    
    observeEvent(input$repetitions, {
        updateProgressBar(session = session, 
                          id = "pb2", 
                          value = input$repetitions/10, 
                          status = color2())
    })
    
    output$codeR <- renderPrint({
        cat("myOLS <- function(data, outcomes) {",
            "\n",
            "    Y <- outcomes",
            "    X <- data",
            "\n",
            "    degrees_of_freedom <- nrow(X)-ncol(X)",
            "\n",
            "    beta_estimates <- ginv(t(X)%*%X)%*%t(X)%*%Y",
            "\n",
            "    error_estimates   <- Y - X%*%beta_estimates",
            "    sigma_estimate    <- (t(error_estimates)%*%error_estimates) / degrees_of_freedom",
            "    variance_estimate <- ginv(t(X)%*%X)",
            "    standard_error    <- as.matrix(sqrt(diag(variance_estimate)))",
            "\n",
            "    RSS       <- t(error_estimates)%*%error_estimates"    ,     
            "    y_bar     <- mean(Y)",
            "    yybar     <- Y-y_bar",
            "    TSS       <- t(yybar)%*%yybar ",
            "    R_squared <- 1 - RSS/TSS ",
            "\n",
            "    t_statistic <- beta_estimates / standard_error",
            "    p_values     <- 2*(1-pt(abs(t_statistic), degrees_of_freedom)) %>%",
            "                    round(5)",
            "\n",
            "    estimates <- tibble(\"estimate\"    = c(beta_estimates), ",
            "                        \"std_error\"   = c(standard_error), ",
            "                        \"t_statistic\" = c(t_statistic)," ,
            "                        \"p_value\"     = c(p_values))",
            "\n",
            "    list(\"estimates\" = estimates,",
            "         \"r_squared\" = R_squared)",
        "}",
        sep = "\n")
    })
    
    output$codeMATLAB <- renderPrint({
        cat("function [bhat, se, r2] = myols(X,y)",
            "% MYOLS Perform OLS regression",
            "%   [bhat, se, r2] = MYOLS(x,y). x is the data input. y is the output",
            "\n",
            "%% Determine dimensions of data",
            "[n, k] = size(X);   % sample size and number of parameters estimated",
            "df = n-k;           % degrees of freedom",
            "\n",
            "%% OLS parameter estimates",
            "bhat = inv(X'*X)*X'*y % use OLS -- the matrix version is simplest.",
            "\n",
            "%% Standard errors",
            "ehat = y - X*bhat;          % residuals",
            "sighat = (ehat'*ehat) / df;    % calculate sigma hat",
            "Vhat = inv(X'*X)*sighat;    % Var-Cov matrix on parameter estimates",
            "se = sqrt(diag(Vhat)) % standard errors ",
            "\n",
            "%% Coefficient of determination",
            "RSS = ehat'*ehat;           % residual sum of squares",
            "ybar = mean(y);",
            "yybar = y-ybar;",
            "TSS = yybar'*yybar;         % total sum of squares",
            "r2 = 1 - RSS/TSS           % R-squared",
            sep = "\n")
    })
    
    output$equation1 <- renderUI({
        withMathJax(helpText('Data generating process: 
                             $$y_i = -3 + 2X_i + 0.5D_i + \\epsilon_i$$'))
    })
    
    output$equation2 <- renderUI({
        if(input$omit == TRUE) {
            withMathJax(helpText('Your model: 
                                 $$y_i = -3 + 2X_i + \\epsilon_i$$'))
        } else {
            withMathJax(helpText('Your model: 
                                 $$y_i = -3 + 2X_i + 0.5D_i + \\epsilon_i$$'))
        }
        
    })
    
    pseudo_data <- eventReactive(input$calculate, {
        
        if(input$omit == FALSE) {
            results1 <- array(numeric(), c(3, 4, input$repetitions))
        } 
        if(input$omit == TRUE) {
            results1 <- array(numeric(), c(2, 4, input$repetitions))
        } 
        
        for(i in 1:input$repetitions) {
            X_i <- as.matrix(rnorm(1:input$observations, mean = 10, sd = 5))
            
            error_1 <- as.matrix(rnorm(1:input$observations), mean = 0, 
                                 sd = sigmas[1, 1])
            error_2 <- as.matrix(rnorm(1:input$observations), mean = 0, 
                                 sd = sigmas[2, 1])
            
            D <- alphas[1,1] + alphas[2,1]*X_i + error_2
            
            y <- betas[1,1] + betas[2,1]*X_i + betas[3,1]*D + error_1
            
            pseudo_data <- list("X_i" = X_i, 
                                "D"   = D, 
                                "Y"   = y)
            
            outcomes <- pseudo_data$Y
            
            if(input$omit == FALSE) {
                data <- cbind(rep(1, input$observations), 
                              pseudo_data$X_i, 
                              pseudo_data$D)
            } 
            
            if(input$omit == TRUE) {
                data <- cbind(rep(1, input$observations), 
                              pseudo_data$X_i)
            }
            
            results <- myOLS(data, outcomes)[[1]]
            
            results1[,,i] <- (as.matrix(results))
        }
        
        return(results1)
        
    })
    
    results <- reactive({
        results <- clean_results(pseudo_data(), isolate(input$omit))
        results <- bind_rows(results, .id = "parameter")
        
        names(results) <- c("parameter", "estimate",
                            "standard_error", "t_statistic",
                            "p_value")
        
        return(results)
    })
    
    output$estimates_table <- renderTable({
        results() %>% 
            dplyr::group_by(parameter) %>%
            summarize(estimate = mean(estimate),
                      standard_error = mean(standard_error),
                      t_statistic = mean(t_statistic),
                      p_value = mean(p_value)) %>%
            bind_cols(betas[1:2,]) %>%
            rename("actual" = "...6") %>%
            dplyr::select(parameter, actual, 
                          estimate, standard_error, 
                          t_statistic, p_value)
    })
    
    output$complete_data <- renderTable({
        results()
    })
    
    output$histogram <- plotly::renderPlotly({
        results() %>%
            ggplot() +
            geom_histogram(aes(estimate), 
                           fill = "royalblue", 
                           color = "navyblue", 
                           bins = 30) +
            #geom_vline(xintercept = betas[1], color = "red") +
            facet_wrap(.~parameter, scales = "free") +
            theme_minimal() +
            theme(axis.title.y = element_blank())
    })
    
    output$histogram2 <- plotly::renderPlotly({
        results() %>%
            filter(parameter == "beta_0") %>%
            ggplot() +
            geom_histogram(aes(estimate), 
                           fill = "royalblue", 
                           color = "navyblue", 
                           bins = 30) +
            geom_vline(xintercept = betas[1], 
                       color = "red",
                       linetype = "dashed") +
            theme_minimal() +
            xlab("Estimate of Beta 0 (blue) vs. actual (red)") +
            theme(axis.title.y = element_blank())
    })
    
    output$histogram3 <- plotly::renderPlotly({
        results() %>%
            filter(parameter == "beta_1") %>%
            ggplot() +
            geom_histogram(aes(estimate), 
                           fill = "royalblue", 
                           color = "navyblue", 
                           bins = 30) +
            geom_vline(xintercept = betas[2], 
                       color = "red",
                       linetype = "dashed") +
            theme_minimal() +
            xlab("Estimate of Beta 1 (blue) vs. actual (red)") +
            theme(axis.title.y = element_blank())
    })
    

}

# Run the application 
shinyApp(ui = ui, server = server)
