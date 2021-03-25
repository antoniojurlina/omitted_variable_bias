# Omitted Variable Bias Shiny App

In statistics, omitted-variable bias (OVB) occurs when a statistical model leaves out one or more relevant variables. The bias results in the model attributing the effect of the missing variables to those that were included.

More specifically, OVB is the bias that appears in the estimates of parameters in a regression analysis, when the assumed specification is incorrect in that it omits an independent variable that is a determinant of the dependent variable and correlated with one or more of the included independent variables.

This is a repository for a Shiny app I made that explores the effect of estimating a model with an omitted variable. In the app I provide code for a matrix-based estimation of OLS parameters in R and in MATLAB. Underneath it all is a Monte Carlo simulation that repeats the estimation process in order to find a mean around which the estimates converge.  

The app itself can be found here: https://ant-onio.shinyapps.io/omitted_variable_bias/
