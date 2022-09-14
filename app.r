library(shiny); library(popbio);library(truncnorm)

#function for simulation populations based on known lambdas
pop_knownLambda <- function(size,time,ext_thr,lambdas){
  r.lamb <- sample(lambdas,time,T)
  
  N <- c(size,size*cumprod(r.lamb)) 
  N[which(N<ext_thr)]=0
  zero.check <- cumprod(N)
  zero.check[which(zero.check>0)]=1
  N*zero.check
}

#function for simulating populations based on mean & variance of lambdas
##variance is entered as percent, i.e. coefficient of variation
mean.lamb=1.01;covar.lamb=10
pop_LambdaDist <- function(size,time,ext_thr,mean_lamb,covar_lamb){
  sd_lamb=((covar_lamb/100)*mean_lamb)
  r.lamb <- rtruncnorm(n=time,a=0,mean=mean_lamb,sd=sd_lamb)
  
  N <- c(size,size*cumprod(r.lamb)) 
  N[which(N<ext_thr)]=0
  zero.check <- cumprod(N)
  zero.check[which(zero.check>0)]=1
  N*zero.check
}

#function for calculating mean and variance log(lambda) from user inputs (mean & CV)
lambda_calc <- function(mean_lamb,covar_lamb){
  mean_loglamb <- log(mean_lamb)
  sd_lamb=((covar_lamb/100)*mean_lamb)
  #can't use log(sd_lamb)
  var_loglamb <- var(log(rtruncnorm(n=10000,a=0,mean=mean_lamb,sd=sd_lamb)))
  out <- c(mean_loglamb,sd_lamb,var_loglamb)  
}

ui <-fluidPage(
  titlePanel("Stochasticity, extinction, and population viability analyses"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Population Parameters"),
      numericInput("N", "Starting population size:", 10, min = 1, max = 100),
      numericInput("et", "Extinction threshold", 3, min = 0, max = 500),
      # radioButtons(inputId = "fixed_values",label = "Fixed values?",
      #              c("Yes!"="T",
      #                "No, means & variances!"="F")),
      # textAreaInput(inputId = "lambdas",label="Enter lambdas here (or leave blank)", value=1,width="100%",height="50%"),
      numericInput("mu_lambda", "Average lambda", 0.95, min = 0.05, max = 2),
      numericInput("CV_lambda", "% variation in lambda *", 15, min = 0, max = 100),
      h6("*as coefficient of variation in lambda"),
      
      h3("Simulation Parameters"),
      numericInput("time","Generations (max: 150)", 20, min = 1, max = 150),
      numericInput("sims", "Number of populations to simulate (max: 25):", 10, min = 1, max = 25),
      actionButton("action",label="Go!")
    ),                
    
    
    mainPanel(
      h4(" "),
      h3("Population Output"),
      h3(" "),
      h4("PVA estimates"),
      textOutput("ext_yr"),
      textOutput("end_size"),
      textOutput("ext_ct"),
      textOutput("ext_prp"),
      h3(" "),
      h4("Simulated Populations",width="100%",height="400px"),
      plotOutput("simPops",width="100%",height="400px"),
      h4("Extinction Probability by Generation*"),
      plotOutput("CDF_plot",width="100%",height="400px"),
      h3(" "),
      h3(" "),
      h3(" "),
      p("* based on Equation 4 of Morris & Doak (2002) Quantitative Conservation Biology: Theory and practice of population viability analysis. Sinauer Associates Inc."),
      h3(" "),
      h3(" "),
      h3(" "),
      p("R and Shiny code available from Ned Dochtermann @ https://github.com/DochtermannLab"),
      p("Code modifiable according to CC BY-NC-ND 4.0")
    )
  )  )

server <- function(input, output) {
  observeEvent(input$action,{
    sims <- reactiveValues()
    
    sims$reps <- replicate(input$sims,pop_LambdaDist(size=input$N,time=input$time,ext_thr=input$et,
                                                     mean_lamb=input$mu_lambda,covar_lamb=input$CV_lambda))
    
    sims$loglam <-   lambda_calc(input$mu_lambda,covar_lamb = input$CV_lambda)[3]
    
    sims$CDF_out <- extCDF(mu=log(input$mu_lambda),sig2=sims$loglam,Nc=input$N,Ne=input$et,tmax=100000)
    
    output$simPops <- renderPlot({
      matplot((sims$reps),ylim=c(0,max((sims$reps))),type='l',lty=1,lwd=2,
              xlab="Generation (t)",ylab="Population size (N)")
      #      axis(2,at = 0:max(log(sims$reps)),labels=c(0,10^(1:max(log(sims$reps)))))
      abline(h=log10(input$et))
    })
    
    
    output$CDF_plot <- renderPlot({
      plot(sims$CDF_out[1:input$time],ylim=c(0,1),xlim=c(1,input$time),
           ylab="Probability of extinction",xlab="Generation (t)",type='l')
    }
    )
    output$ext_yr <- renderText({
      paste("Average time to extinction: ",min(which(sims$CDF_out>0.5)))
    })
    output$end_size <- renderText({
      paste("Average population size at t_max:",round(mean(sims$reps[dim(sims$reps)[1],]),1))
    })
    output$ext_ct <- renderText({
      paste("Number of extinct simulated populations:",sum(sims$reps[dim(sims$reps)[1],]<input$et))
    })
    output$ext_prp <- renderText({
      paste("Proportion of simulated populations extinct:",sum(sims$reps[dim(sims$reps)[1],]<input$et)/dim(sims$reps)[2])
    })
    
  }
  )
}

# Run the app ----
shinyApp(ui = ui, server = server)
