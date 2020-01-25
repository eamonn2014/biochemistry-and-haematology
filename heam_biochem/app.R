#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rshiny ideas from on https://gallery.shinyapps.io/multi_regression/
# https://r.789695.n4.nabble.com/Whiskers-on-the-default-boxplot-graphics-td2195503.html
# https://www.r-bloggers.com/whisker-of-boxplot/
# https://journals.sagepub.com/doi/pdf/10.1177/1536867X0900900309
# https://journals.sagepub.com/doi/pdf/10.1177/1536867X1301300214
# https://www.stata.com/support/faqs/graphics/box-plots-and-logarithmic-scales/
# https://stats.stackexchange.com/questions/112705/boxplots-with-lognormally-distributed-data
# https://www.r-graph-gallery.com/96-boxplot-with-jitter.html
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
library(shiny) 
library(nlme)
library(VCA)
library(MASS)
require(tidyverse)
require(ggplot2)
library(shinyWidgets)

options(max.print=1000000)
fig.width <- 1200
fig.height <- 550
fig.height2 <- 450
library(shinythemes)        # more funky looking apps
p1 <- function(x) {formatC(x, format="f", digits=1)}
p2 <- function(x) {formatC(x, format="f", digits=2)}
options(width=100)
#set.seed(12345) #reproducible
 

# function to create longitudinal data

is.even <- function(x){ x %% 2 == 0 }

mdata <- function( n,beta0, beta1, beta2, ar.val, sigma, tau0, tau1, tau01, m ) {
    
    p <- round(runif(n,4,m))
    obs <- unlist(sapply(p, function(x) c(1, sort(sample(2:m, x-1, replace=FALSE)))))
    dat <- data.frame(id=rep(1:n, times=p), obs=obs)
    dat$trt = 1*is.even(as.numeric(as.character(dat$id)))
    dat$trt1 <- ifelse(dat$obs==1,  0, dat$trt) # new manipulate so trt effect after first visit 
    mu  <- c(0,0)
    S   <- matrix(c(1, tau01, tau01, 1), nrow=2)
    tau <- c(tau0, tau1)
    S   <- diag(tau) %*% S %*% diag(tau)
    U   <- mvrnorm(n, mu=mu, Sigma=S)  # rows = patients
    dat$eij <- unlist(sapply(p, function(x) arima.sim(model=list(ar=ar.val), n=x) * sqrt(1-ar.val^2) * sigma))
    dat$yij <- (beta0 + rep(U[,1], times=p)) + (beta1 + rep(U[,2], times=p)) * (dat$obs) + dat$eij + beta2*dat$trt1 # beta2*dat$trt NEW
    dat$trt1 <- NULL # new i want treatmetn effect after first visit
    z <-  groupedData(yij ~ obs + trt | id, data=dat)   #TRT NEW
    
    return(z)
    
}

biochemistry <- c(
    "Fasting blood glucose (mg/dL)"  ,
    "HbA1c glycosylated hemoglobin (%)",  
    "Non fasting blood Glucose (mg/dL)" , 
    "Calcium (mmol/L)"  ,
    "Chloride (mmol/L)" ,
    "Creatine phosphokinase (U/L)",  
    "Inorganic phophorous (mmol/L)",
    "Magnesium (mmol/L)"   ,
    "Potassium (mmol/L)"  ,
    "Sodium (mmol/L)" ,
    "Albumin (g/L)"   ,
    "Blood urea nitrogen (mmol/L)"  ,
    "Creatinine (umol/L)" ,
    "Total protein (g/L)"  ,
    "Uric acid (umol/L)"  ,
    "Direct bilirubin (umol/L)",  
    "Gamma glutamyltransferase (U/L)",  
    "Indirect bilirubin (umol/L)" ,
    "Total bilirubin (umol/L)" ,
    "Alkaline phosphatase (U/L)", 
    "High density lipoprotein (mmol/L)" , 
    "Low density lipoprotein (mmol/L)"  ,
    "Lipase (U/L)",
    "Total cholesterol (mmol/L)",  
    "Triglycerides (mmol/L)"  
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ui <- fluidPage(theme = shinytheme("paper"), #https://www.rdocumentation.org/packages/shinythemes/versions/1.1.2
                
                
                 # tags$style(HTML(" 
                 #          .skin-blue .sidebar-menu > li.active > a { border-left-color: #e6d4d9;}
                 #   ")), 
                
                
                 setBackgroundColor(
              #      color = c("#F7FBFF", "#2171B5"),
                    color = c("#ecf0ba", "#bad1f0"), ##  e6e1d4
                   gradient = "radial",
                    direction = c("bottom", "left")
                 ),
                # .skin-blue .sidebar-menu > li.active > a {
                #   border-left-color: #ff0000;
                # }
                
                shinyUI(pageWithSidebar(
                    
                    #ui <-shinyUI(pageWithSidebar(
                  
                  
                    
                    headerPanel("Presenting the results of diagnostic tests routinely ordered to determine general 
                    health status"),
                    
             
               
                    
                    #sidebarLayout(  #new
                    # Sidebar with a slider and selection inputs
                    
                
                      sidebarPanel( #tags$hr()  ,
        #               
        #                 tags$style(HTML(" 
        #                   .skin-blue .sidebar-menu > li.active > a { border-left-color: #ff3300;}
        #            ")), 
        # #                 tags$head(tags$style(
        #                   HTML('
        #  #sidebar {
        #     background-color: #dec4de;
        # }
        # 
        # body, label, input, button, select { 
        #   font-family: "Arial";
        # }')
        #                 )),
                        
                        
                   #      tags$style(HTML(" 
                   #          .sidebar-default .sidebar-brand {color: cyan;}
                   #  .sidebar { background-color: lightblue;}
                   #          
                   # ")),
       #                  tags$head(tags$style(
       #                    HTML('
       #   sidebarpanel {
       #      background-color: #ffffff;
       #  }
       # 
       # ')
       #                  )),
                        
            
                    div(p("Typically for clinical trials and non-interventional studies reams of
                    output tables are generated presenting biochemistry and haemotology test results, perhaps summarised over time.
                   
                          How useful it is presenting a multitude of output tables of test result summary statistics 
                          is open to debate. Firstly, the presentation is extremely dry. 
                          There is the issue too, that it is
                          quite possible a different data set can give rise to similar summary statistics [1]. 
                          Spotting errors in the data is also not easy and
                          tracking individual patient profiles not possible.
                          I would argue that these parameters should always be 
                          plotted. R Shiny is an ideal medium to communicate this information
                          in a more useful and exciting way.  We will focus on biochemistry tests that are routinely ordered to determine a patient's general 
                    health status. 
                          I will make an attempt using simulated data.")),
                        
                        div(
                            
                            
                            actionButton(inputId='ab1', label="R code",   icon = icon("th"), 
                                         onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/Boxplots/master/app.R', '_blank')"),   
                            actionButton("resample", "Simulate a new sample"),
                            br(), br(),
                            
                            
                            div(strong("Select the parameters using the sliders below"),p(" ")),
                            
                            
                            div(("  
                           
                            We can select an overall plot, showing 
                            all patient profiles and boxplots across visits. Boxplots are generated using the ggplot2 package [2].
                            Each biochemistry test can be inspected one by one. 
                            The '2 Select plot' = 'Individual' allows
                            the inspection of patient profiles of choice to be displayed, by typing in the subject identifier 
                            into the third input option. 
                            Here the patient IDs are 
                 just numeric starting at 1. Selecting patient ID 999 will 
                 show all patients. There is a third plot option, 'Individual all tests' which displays a matrix of scatterplots
                                 all the test results for a single patient. The patient in question is the first entry typed in the the third option.  ")),
                            br(),
                            selectInput("Plot",
                                        strong("1. Select which biochemistry test to present"),
                                        choices=biochemistry),
                            
                            selectInput("Plot1",
                                        strong("2. Select plot"),
                                        choices=c("Overall","Individual","Individual all tests")),
                            
                     
                            textInput('vec1', "3. Select patient. If '2 select plot' 'Individual' is selected, enter sample ID(s) (comma delimited); 
                                      enter 999 to show all profiles; If 'Individual all tests' is selected, all test results for the first ID only are presented", "1,2,3,4"),
                            
                            
                            sliderInput("V",
                                        "4. Maximum visit number in data simulation including baseline",
                                        min=3, max=10, step=1, value=8, ticks=FALSE),
                            
                            sliderInput("VV",
                                        "5. Estimate treatment effect at this visit",
                                        min=1, max=10, step=1, value=4, ticks=FALSE),
                            
                            
                            # sliderInput("N",
                            #             "Select the total number of data points",
                            #             min=3, max=500, step=1, value=100, ticks=FALSE),
                            # 
                            # sliderInput("Whisker",
                            #             "Select the length of whiskers (multiples of the IQR)",
                            #             min=0, max=3, step=.5, value=1.5, ticks=TRUE),
                            # 
                            # sliderInput("outliers",
                            #             "Highlight 'outliers'",
                            #             min=0, max=1, step=1, value=0, ticks=FALSE),
                            # 
                            # sliderInput("dp",
                            #             "Show me the data!",
                            #             min=0, max=1, step=1, value=1, ticks=FALSE),
                            
                            div(p("References:")),  
                            
                            # tags$a(href = "https://en.wikipedia.org/wiki/Exploratory_data_analysis", "[1] Tukey, J.W. 'Exploratory Data Analysis'"),
                            # div(p(" ")),
                            # tags$a(href = "https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/boxplot", "[2] R boxplot"),
                            # div(p(" ")),
                            tags$a(href = "https://ggplot2.tidyverse.org/reference/geom_boxplot.html", "[1] Boxplots using ggplot2"),
                            div(p(" ")),
                            tags$a(href = "https://en.wikipedia.org/wiki/Anscombe%27s_quartet", "[2] Anscombe's quartet"),
                            div(p(" ")),
                         
                            # tags$a(href = "https://r.789695.n4.nabble.com/Whiskers-on-the-default-boxplot-graphics-td2195503.html", "[5] Whiskers on the default boxplot"),
                            # div(p(" ")),
                            
                        )
                        
                    ),
                    
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tab panels
                    mainPanel(
                      
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        #    tabsetPanel(type = "tabs", 
                        navbarPage(       
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
                            tags$style(HTML(" 
                            .navbar-default .navbar-brand {color: cyan;}
                            .navbar-default .navbar-brand:hover {color: blue;}
                            .navbar { background-color: #bad1f0;}
                            .navbar-default .navbar-nav > li > a {color:black;}
                            .navbar-default .navbar-nav > .active > a,
                            .navbar-default .navbar-nav > .active > a:focus,
                            .navbar-default .navbar-nav > .active > a:hover {color: pink;background-color: purple;}
                            .navbar-default .navbar-nav > li > a:hover {color: black;background-color:yellow;text-decoration:underline;}
                            .navbar-default .navbar-nav > li > a[data-value='t1'] {color: red;background-color: pink;}
                            .navbar-default .navbar-nav > li > a[data-value='t2'] {color: blue;background-color: lightblue;}
                            .navbar-default .navbar-nav > li > a[data-value='t3'] {color: green;background-color: lightgreen;}
     
                   ")), 
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end of section to add colour     
                            tabPanel("Plotting the data", 
                                 #    h2("Plotting the data"),
                                     div(plotOutput("reg.plot3", width=fig.width, height=fig.height)),  
                                     
                                     
                                    h3(" "),
                                     
                                     p(strong(
                                      "On this first tab, first select the test that one wants to examine. 
The plot selection 'Overall' allows a visual comparison between the Active and Placebo responses for all patients over all visits.
The plot selection 'Individual' allows one to select any particular patient and visualise their profile. More than one patient can be examined simutaneously. 
The third 'Select plot' option allows one to look at only patient at a time but all their test results simultaneously.
Moving on to the Summary statistics, we see the said statistics for the selected test. These are typically what is presented.
These are often supplemented with change from baseline and percentage change frome baseline. However the purpose of a parallel-group randomized trial is 
to compare the parallel groups, not to look at change from baseline. 
                                      Baseline should always be an adjustment covariate (only). That is precisely what we do on the next tab, where the GLS model output is presented. 
                                      You can imagine this could be adjusted for other important covariates possibly age and sex.")) ,
                                     
                 #                     
                  p(strong("The next tab presents the model treatment effect estimates graphically and in a table for each visit with 95CIs. We do not make any adjustments for modelling over 20 diagnostic tests. ")),
                  
                 #  Often you may hear it is said, 'the limits of the boxes represent Q1 and Q3 and the 
                 #  whiskers extend +/- 1.5 x Q3-Q1'. 
                 #  This is not necessarily true. There are many variations on the box plot, it is therefore better to explicitly state what is being presented.")) ,
                 #                     div(""),
                 #                     p(strong("I hope you agree the plot on the right gives a better understanding of the data distributions. 
                 # In this programming exercise select 'Show me the data!' and deselect 'Highlight 'outliers' '.")) ,
                                     
                            ) ,
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            tabPanel("Summary statistics", value=3, 
                                     #  div( verbatimTextOutput("table2")),
                                     h3("Summary statistics, typically generated as outputs for clinicial scrutiny"),
                                     DT::dataTableOutput("table2"),
                                     #
                                    # tags$head(tags$style("#dummy table {background-color: red; }", media="screen", type="text/css")),
                                     
                            ) ,
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            tabPanel("Statistical modelling", value=3, 
                                     h3("Modelling"),
                                     # div(plotOutput("reg.plot2", width=fig.width, height=fig.height)),  
                                     # h3("Figure 2 Top panel untransformed data, bottom panel using a natural log transformation"),
                                     # 
                                     # h4("Lets first look at the relationship between the hinges vrs. quartiles.") ,
                                     # p(strong("Boxplot stats for the untransformed data, as presented in the top panel:")) ,
                                     # div( verbatimTextOutput("table2")),
                                     # p(strong("Now summarise the same data. Q1 and Q3 may not match the hinges. Paraphrasing the 'box.plot.stats' help file...'The ‘hinges’ are versions of Q1 and Q3'. See what happens when n is odd. The hinges and the quartiles match.")) ,
                                     # div( verbatimTextOutput("table3")),
                                     # 
                                     # p(strong("Now check out the boxplot stats based on the log transformed data, after exponentiating (bottom panel), check with the raw median:")) ,
                                     # div( verbatimTextOutput("table4")),
                                     # p(strong("In the same way summarise the log transformed data then exponentiating: ")) ,
                                     # div( verbatimTextOutput("table5")),
                                     
                                     # p(strong("Use Harrell's rms function 'contrast' to estimate the treatment effect at desired visit by selecting 'What timepoint to estimate trt B- trt A?:'")),
                                     # div(class="span7", verbatimTextOutput("reg.summary4")),
                                     p(strong("Be patient...A generalised least squares model fit with unstructured covariance structure will appear, 
                                     the reference level for the visit variable is selected using the slider '5. Estimate treatment effect at this visit'.
                                     This means we can simply read off the treatment effect ' trt=Placebo ' directly from the model output, 
                                              for the treatment effect estimate comparing Placebo - Active at that particular visit. Note often a log transformation of laboratory test data will be fruitful as the data is often skewed and negative values are not expected.")),
                                     div(class="span7", verbatimTextOutput("table4")),
                                      div(class="span7", verbatimTextOutput("reg.summary2")),
                                     
                 #                     
                 #                     h4("Look again at how to calculate the boxplot statistics when transforming") ,
                 #                     
                 #                     p(strong("The length of the whiskers is typically set at 1.5xIQR, that is, the 
                 # whiskers extend to the most extreme data point which is no more than 1.5 times the interquartile 
                 # range from the box. A value of zero, using the R function, causes the whiskers to extend to the data extremes. 
                 #          So it may not be as simple as saying the whiskers extend 1.5xIQR.")),
                 #                     
                 #                     
                 #                     div(""),
                 #                     
                 #                     p(strong("To summarise, apply the whisker calculation rule on the scale used to draw the boxplot. Therefore when using a transformation 
                 # the calculations must 
                 # be on the transformed scale. Do not mix scales. So, do not calculate the whiskers based on the 1.5xIQR rule
                 # on the raw scale and then log transform 
                 #          to present the boxplot.  
                 #          The median and hinges will be the logs of the original median and hinges, they will match with odd sample sizes, but not generally with even. 
                 #          The step which determines the length of the whiskers will change. Hence the fences generally will differ. 
                 #          (The location of the end of the whiskers is referred to as the fences). ")) ,
                 #                     
                 #                     
                 #                     p(strong("")),
                 #                     h3("Take home messages; when the sample size is even, the hinges may not necessarily equal Q1 and Q3 (raw and transformed) and the transformed median may be slightly different to the raw median.
                 #      The fences may not necessarily equal the fences on the raw scale for any sample size.
                 #    Explicitly state how you constructed your boxplots.")
                            ) ,
                            
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                             tabPanel("Plot of the treatment effect estimates", 
                                      h3("Plot of the treatment effect estimates"),
                                      div(plotOutput("reg.plote", width=fig.width, height=fig.height2)),  
                                      h3("Table of the treatment effect estimates"),
                                     # p(strong("Use Harrell's rms function 'contrast' to estimate the treatment effect at each visit:")),
                                      DT::dataTableOutput("reg.summary4"),
                                     # div(class="span7", verbatimTextOutput("reg.summary4")),
                                      
                                      
                             ) ,
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                          
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                          
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            tabPanel("Diagnostics",
                                     h3("Three residual plots to check for absence of trends in central tendency and in variability"),
                                     div(plotOutput("res.plot", width=fig.width, height=fig.height)),       
                                     p(strong("Upper left panel shows the baseline score on the x-axis. Upper right panel shows shows time on the x-axis. Bottom left panel is the QQ plot for checking normality of residuals from the GLS fit.")),
                            ),
              
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                             tabPanel("Data listing", value=3, 
                                      h3("Data listing"),
                                      DT::dataTableOutput("table1"),
                                      #div( verbatimTextOutput("table1")),
                             ) 
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        )
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    )
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end tab panels 
                    
                )
                )
)

server <- shinyServer(function(input, output   ) {
    
    
    # --------------------------------------------------------------------------
    # This is where a new sample is instigated 
    random.sample <- reactive({
        
        # Dummy line to trigger off button-press
        foo <- input$resample
        Whisker<-   isolate(input$Whisker )
        outliers <- isolate(input$outliers )
        n=(input$N )
        dp=isolate(input$dp )
        
        return(list( n=n ,  Whisker=Whisker , outliers=outliers, dp=dp )) 
        
    })
    
    
    make.data <- reactive({
        
        sample <- random.sample()
        
        is.even <- function(x){ x %% 2 == 0 }
        
        mdata <- function( n,beta0, beta1, beta2, ar.val, sigma, tau0, tau1, tau01, m ) {
            
            p <- round(runif(n,4,m))
            obs <- unlist(sapply(p, function(x) c(1, sort(sample(2:m, x-1, replace=FALSE)))))
            dat <- data.frame(id=rep(1:n, times=p), obs=obs)
            dat$trt = 1*is.even(as.numeric(as.character(dat$id)))
            dat$trt1 <- ifelse(dat$obs==1,  0, dat$trt) # new manipulate so trt effect after first visit 
            mu  <- c(0,0)
            S   <- matrix(c(1, tau01, tau01, 1), nrow=2)
            tau <- c(tau0, tau1)
            S   <- diag(tau) %*% S %*% diag(tau)
            U   <- mvrnorm(n, mu=mu, Sigma=S)  # rows = patients
            dat$eij <- unlist(sapply(p, function(x) arima.sim(model=list(ar=ar.val), n=x) * sqrt(1-ar.val^2) * sigma))
            dat$yij <- (beta0 + rep(U[,1], times=p)) + (beta1 + rep(U[,2], times=p)) * (dat$obs) + dat$eij + beta2*dat$trt1 # beta2*dat$trt NEW
            dat$trt1 <- NULL # new i want treatmetn effect after first visit
            z <-  groupedData(yij ~ obs + trt | id, data=dat)   #TRT NEW
            
            return(z)
            
        }
        
        i=input$V
        
       # beta1 = slope
      #  beta2 = tr effect
        
        
        # trt <- c(-5,.4,-5,0,2,2,.1,.1,.1,.1,.1,.05,5,.5,.5,.1,1)
        
        
        z <- mdata(n=99,beta0=100, beta1=0, beta2=-5, ar.val=.9, sigma=5, tau0=7, tau1=1, tau01=0.01, m=i)
        df1 <- cbind(test=biochemistry[1], z)

        z <- mdata(n=85,beta0=5.6, beta1=0, beta2=-.04, ar.val=.9, sigma=.4, tau0=.2, tau1=.04, tau01=0.1, m=i)
        df2 <- cbind(test=biochemistry[2], z)

        z <- mdata(n=102,beta0=110, beta1=0, beta2=-5, ar.val=.9, sigma=5, tau0=.5, tau1=.5, tau01=0.01, m=i)
        df3 <- cbind(test=biochemistry[3], z)

        z <- mdata(n=102,beta0=2.4, beta1=.1, beta2=0, ar.val=.9, sigma=.2, tau0=.5, tau1=.01, tau01=0.2, m=i)
        df4 <- cbind(test=biochemistry[4], z)

        z <- mdata(n=90,beta0=104, beta1=1, beta2=2, ar.val=.8, sigma=4, tau0=5, tau1=1, tau01=.6, m=i)
        df5 <- cbind(test=biochemistry[5], z)

        z <- mdata(n=90,beta0=95, beta1=0, beta2=2, ar.val=.8, sigma=4, tau0=15, tau1=1, tau01=.6, m=i)
        df6 <- cbind(test=biochemistry[6], z)

        z <- mdata(n=90,beta0=1.23, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df7 <- cbind(test=biochemistry[7], z)


        z <- mdata(n=90,beta0=1, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df8 <- cbind(test=biochemistry[8], z)


        z <- mdata(n=99,beta0=4, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df9 <- cbind(test=biochemistry[9], z)

        z <- mdata(n=99,beta0=141, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df10 <- cbind(test=biochemistry[10], z)

        z <- mdata(n=111,beta0=44, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df11 <- cbind(test=biochemistry[11], z)

        z <- mdata(n=88,beta0=5, beta1=0, beta2=.05, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df12 <- cbind(test=biochemistry[12], z)

        z <- mdata(n=112,beta0=70, beta1=0, beta2=.5, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df13 <- cbind(test=biochemistry[13], z)

        z <- mdata(n=112,beta0=72, beta1=0, beta2=.5, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df14 <- cbind(test=biochemistry[14], z)

        z <- mdata(n=97,beta0=3.4, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df15 <- cbind(test=biochemistry[15], z)

        z <- mdata(n=103,beta0=24, beta1=0, beta2=1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df16 <- cbind(test=biochemistry[16], z)

        z <- mdata(n=105,beta0=6, beta1=0, beta2=.5, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df18 <- cbind(test=biochemistry[18], z)

        z <- mdata(n=100,beta0=9, beta1=0, beta2=.0, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df19<- cbind(test=biochemistry[19], z)

        z <- mdata(n=103,beta0=64, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df20<- cbind(test=biochemistry[20], z)

        z <- mdata(n=97,beta0=1.3, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df21<- cbind(test=biochemistry[21], z)

        z <- mdata(n=98,beta0=3.3, beta1=0, beta2=-.2, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df22<- cbind(test=biochemistry[22], z)

        z <- mdata(n=104,beta0=34, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df23<- cbind(test=biochemistry[23], z)

        z <- mdata(n=99,beta0=5.11, beta1=0, beta2=.0, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df24<- cbind(test=biochemistry[24], z)

        z <- mdata(n=100,beta0=1.23, beta1=0, beta2=.0, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df17<- cbind(test=biochemistry[17], z)
        
        z <- mdata(n=100,beta0=1.23, beta1=.1, beta2=.0, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=i)
        df25<- cbind(test=biochemistry[25], z)


        d1 <- do.call(rbind, list(df1 = df1, df2 = df2, df3=df3, df4=df4, df5=df5, df6=df6, df7=df7,
                                  df8=df8,df9=df9,df10=df10,df11=df11,df12=df12,df13=df13,df14=df14,
                                  df15=df15, df16=df16,df17=df17, df18=df18, df19=df19, df20=df20,
                                  df21=df21,df22=df22, df23=df23, df24=df24, df25=df25))

        colnames(d1)[colnames(d1)=="yij"] <- "hillest"
        colnames(d1)[colnames(d1)=="trt"] <- "tailindex"
        colnames(d1)[colnames(d1)=="obs"] <- "memorypar"
        colnames(d1)[colnames(d1)=="id"] <- "rep"
        d1$memorypar <- factor(d1$memorypar)
        d1$tailindex <- factor(d1$tailindex)
        d1$tailindex <- ifelse(d1$tailindex %in% 1, "Active","Placebo" )
           
        return(list(  d1=d1)) 
        
        
    })
    
    make.data2 <- reactive({
        
        sample <- random.sample()
       
        
    })
    
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    fit.regression <- reactive({
        
       d <- make.data()$d1
        
       target <- input$Plot
            
       d <- d[d$test %in% target,]
        
       d$time <- d$memorypar
           
       require(rms)
       
      # d$time <- factor(d$time) 
       
       d$time <- relevel(d$time, ref=input$VV)
        
       d$trt <- factor(d$tailindex) 
       
       d$yij <- d$hillest
       
       # new~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       #convert time=1 to baseline var
       d <- d[, c("test", "rep", "time", "trt", "yij")]
       d$id=1:nrow(d)  # add index
       baseline <- d[d$time %in% 1,] # create baseline removing baseline
       names(baseline) <-  c("test", "rep", "time", "trt", "baseline", "id")
       baseline$id <- NULL
       baseline$time <- NULL
       d2 <- d[!d$time ==1,]         # create follow up
       both <- merge (baseline , d2   , all=TRUE)
       both$rep <- as.numeric(as.character(both$rep))
       both <- plyr::arrange(both, rep, time)
       both$time <- as.numeric(as.character(both$time))
       d <- both
       d$time= d$time-1
       d$time<-factor(d$time)
       d$time <- relevel(d$time, ref=input$VV)
       z<-d
       #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       
       
       ddz <<- datadist(d)  # need the double in this environ <<
        
        options(datadist='ddz')
        
            
            fit.res <-  
                tryCatch(Gls(yij  ~ baseline+ time * trt ,
                             correlation=corSymm(form=~ as.numeric(time)|rep),
                             weights=varIdent(form=~1|time),
                             d, x=TRUE,
                             na.action=na.exclude ), 
                         error=function(e) e)
 
            return(list(fit.res=fit.res , target = target, z=z ))
    })     
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # treatment effect estimate
    output$reg.summary4 = DT::renderDataTable({
   
 
        f <- fit.regression()
        
        fit <- f$fit.res
        
        time. <- rep(1:(input$V-1))
        
        k1 <- contrast(fit, list(time=time.,  trt = 'Placebo'),
                            list(time=time.,  trt = 'Active'))
        
        # k1 <- contrast(fit , list(time= input$VV, trt="Active"), 
        #                    list(time= input$VV, trt="Placebo"), type="average") 
        
        
        x <- as.data.frame(k1[c('time', 'Contrast', 'Lower', 'Upper')]) 
        
        #return(list( est2=x ))
        
        library(DT)
        # 
         # ff <-   x %>%
         #   datatable(  ) %>%
         #  formatRound(
         #    columns=c('time', 'Contrast', 'Lower', 'Upper'), digits=2)  
        # datatable(x,   options = list(dom = 't' ) )
         
         #https://datatables.net/reference/option/
         datatable(x,   
            
                   rownames = FALSE,
         
        options = list(
           searching = FALSE,
           pageLength = input$V-1,
           paging=FALSE,
           lengthMenu = FALSE ,
           lengthChange = FALSE,
           autoWidth = TRUE
           # colReorder = TRUE,
           # deferRender = TRUE,
           # scrollY = 200,
           # scroller = T
         ))  %>%
             formatRound(
                columns=c('time', 'Contrast', 'Lower', 'Upper'), digits=c(0,2,2,2)  )
         
         # output$myTable <- renderDataTable(df, 
         #                                   options = list(pageLength = 15, lengthChange = FALSE),
         #                                   rownames= FALSE)
         # 
        # datatable(m) %>% formatCurrency(1:2, '\U20AC', digits = 0) %>% formatRound(  3)
       #  datatable(x) %>%
       
        #   formatRound('time', 'Contrast', 'Lower', 'Upper',digits=2) #%>%
          
 
         
         
         
    })     
    
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    # Plot the estimated trt effect  
    
    output$reg.plote <- renderPlot({         
     
        xl <- xlab( 'Visit')
         
        f <- fit.regression()
        
        fit <- f$fit.res
        
        time. <-  rep(1:(input$V-1))
        
        k1 <- contrast(fit, list(time=time.,  trt = 'Placebo'),
                            list(time=time.,  trt = 'Active'))
        
        k1 <- as.data.frame(k1[c('time', 'Contrast', 'Lower', 'Upper')]) 
        
        mi <- floor(min(k1$Lower))
        ma <- ceiling(max(k1$Upper))
        
        ggplot (k1, aes(x=time. , y=Contrast, group=1)) + geom_point () + geom_line () +
            ylim(mi,ma) +
            xlim(1, input$V-1) +
            scale_x_continuous(breaks=c(time.)) +
            ylab( 'Placebo - Active')+ xl +
            geom_errorbar(aes(ymin=Lower, ymax=Upper ), width =0) +
            ggtitle(paste0("Outcome measure ", input$Plot ,"; treatment effect estimate at each visit with 95% CI")) +
                geom_hline(aes(yintercept = 0, colour = 'red'), linetype="dashed") +
            theme_bw() +
            theme(legend.position="none") +
            theme(panel.background=element_blank(),
                  # axis.text.y=element_blank(),
                  # axis.ticks.y=element_blank(),
                  # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
                  # stop axis being clipped
                  plot.title=element_text(size = 18), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
                  legend.text=element_text(size=14),
                  legend.title=element_text(size=14),
                  legend.position="none",
                  axis.text.x  = element_text(size=15),
                  axis.text.y  = element_text(size=15),
                  axis.line.x = element_line(color="black"),
                  axis.line.y = element_line(color="black"),
                  plot.caption=element_text(hjust = 0, size = 7),
                  strip.text.x = element_text(size = 16, colour = "black", angle = 0),
                  axis.title.y = element_text(size = rel(1.5), angle = 90),
                  axis.title.x = element_text(size = rel(1.5), angle = 0),
                  strip.background = element_rect(colour = "black", fill = "white"))
           
      
        
        
    }) 
    
    
    
    
    # --------------------------------------------------------------------------
    # -----------------------------------------------OVERALL PLOT
    # ---------------------------------------------------------------------------
    
           output$reg.plot3 <- renderPlot({         
            
               d <- make.data()$d1
               
               if (input$Plot1 == "Overall") {
                   
            
            target <- input$Plot
            
            d <- d[d$test %in% target,]
            
            # lets get counts to put in ribbons
            d$tailindex <- factor(d$tailindex)
            dx <- unique(d[,c("rep","tailindex")])
            table(dx$tailindex)
            n <- as.vector(table(dx$tailindex))
            levels(d$tailindex) <- c(paste0("Active N=",n[1]), paste0("Placebo N=",n[2])) 
            
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            pd <- position_dodge(.4)
            pr1=NULL
            pr1 <- ggplot(d,aes(x=memorypar ,y=hillest,color=tailindex, fill=tailindex ))  + 
                stat_boxplot(geom = "errorbar", width = 0.3) +
                geom_boxplot( outlier.colour = NA  ) +  # removed fill=NA
                geom_line(aes(group=rep), position = pd,  alpha=0.6, linetype="dotted")   + 
                scale_size_manual( values = c( 1) ) +
                geom_point(aes(fill=tailindex, group=rep), pch=1, size=1, alpha=0.3, position = pd ) +
                stat_summary(fun.y=mean, geom="point", shape=3, size=2, colour="black", stroke=1.5,
                             position=pd, show.legend=FALSE) +
                scale_color_manual(name = "Treatment", values = c("blue", "darkgreen")) +
                scale_fill_manual(name = "Treatment", values = c("lightblue", "green")) +
                facet_wrap(~tailindex , ncol=2)    +
                labs(caption = "- The upper whisker is located at the smaller of the maximum y value and Q3 + 1.5xIQR, whereas the lower whisker is located at the larger of the smallest y value and Q1 – 1.5xIQR\n- The median is the horizontal line inside each box and the mean denoted by the cross\n -Individual patient profiles are denoted by dotted lines\n- A small amount of jitter is added to the data to aid visualisation.") +
                
                geom_text(data = d %>% group_by( memorypar, tailindex) %>%
                              summarise(Count = n()) %>%
                              ungroup %>%
                              mutate(hillest=min((d$hillest)) - 0.05 * diff(range((d$hillest)))),
                          aes(label = paste0("n = ", Count)),
                          position = pd, size=3, show.legend = FALSE) 
            
            #####################
            print(pr1 + labs(y=target, x = "Visit") + 
                      ggtitle(paste0("There are N=",
                                     length(unique(d$rep)),  
    " patients with data at baseline, presenting all patient profiles, with boxplots and the number of patient values at each visit") ) +
                      theme_bw() +
                      theme(legend.position="none") +
                      theme(#panel.background=element_blank(),
                            # axis.text.y=element_blank(),
                            # axis.ticks.y=element_blank(),
                            # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
                            # stop axis being clipped
                            plot.title=element_text(size = 18), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
                            legend.text=element_text(size=14),
                            legend.title=element_text(size=14),
                            legend.position="none",
                            axis.text.x  = element_text(size=15),
                            axis.text.y  = element_text(size=15),
                            axis.line.x = element_line(color="black"),
                            axis.line.y = element_line(color="black"),
                            plot.caption=element_text(hjust = 0, size = 11),
                            strip.text.x = element_text(size = 16, colour = "black", angle = 0),
                            axis.title.y = element_text(size = rel(1.5), angle = 90),
                            axis.title.x = element_text(size = rel(1.5), angle = 0),
                            strip.background = element_rect(colour = "black", fill = "gray96"),
                           panel.background = element_rect(fill = 'gray96', colour = 'white'),
                           plot.background = element_rect(fill = 'gray96', colour = 'white'),#
                           
                            
                            
                      ) 
            )
           
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Individual profiles
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
               }  else  if (input$Plot1 == "Individual") { 
                    
                   
                   i <- as.numeric(unlist(strsplit(input$vec1,",")))

                   target <- input$Plot
                   
                   d <- d[d$test %in% target,]
                   
                   d$tailindex <- factor(d$tailindex)
                   dx <- unique(d[,c("rep","tailindex")])
                   table(dx$tailindex)
                   n <- as.vector(table(dx$tailindex))
                   levels(d$tailindex) <- c(paste0("Active N=",n[1]), paste0("Placebo N=",n[2])) 
                   #  if 999 is entered all subjects are shown
                   
                   if("999" %in% i) {
                       
                       dd<-d
                       
                   } else {
                       
                       dd <- d[d$rep %in% i,]
                   }
                   
                   
                   sel <- unique(dd$tailindex) #if only one arm , dont show the empty arm
                   d <- d[d$tailindex %in% sel,]
                   
                    dx <- unique(dd[,c("rep","tailindex")])
                   
                     nn <- as.vector(table(dx$tailindex))
                     
                     levels(d$tailindex)  <- levels(dd$tailindex) <- 
                       c(paste0("Active N=",n[1], " with " ,nn[1]," patient profile(s) shown"), 
                        paste0("Placebo N=",n[2], " with " ,nn[2]," patient profile(s) shown"))
                   
                      
                     
                   pd <- position_dodge(.4)
                   pr1=NULL
                   pr1<-ggplot(d,aes(x=memorypar ,y=hillest,color=tailindex, fill=tailindex )) + 
                        stat_boxplot(geom = "errorbar", width = 0.3) +
                       geom_boxplot( outlier.colour = NA) +#,alpha=0.1, color="lightblue",)  +  
                       
                       geom_line(data = dd,
                                 aes(group=rep,x = memorypar, y = hillest),  size = .6, linetype="dashed") +
                           scale_size_manual( values = c( 1) ) +
                       
                       geom_point(data=dd, aes(fill=tailindex, group=rep), 
                                  pch=19, size=3, colour='black',alpha=0.7, position = pd ) +

                       geom_point(aes(fill=tailindex, group=rep), 
                                  pch=1, size=2, alpha=0.2, position = pd ) +
                  
                       stat_summary(fun.y=mean, geom="point", shape=3, size=2, colour="black", stroke=1.5,
                                    position=pd, show.legend=FALSE) +
                     
                     
                       scale_color_manual(name = "Treatment", values = c("blue", "darkgreen") ) +
                       scale_fill_manual(name = "Treatment", values = c("lightblue", "green") ) +
                       
                       
                       facet_wrap(~tailindex , ncol=2)    +
                     labs(caption = "- The upper whisker is located at the smaller of the maximum y value and Q3 + 1.5xIQR, whereas the lower whisker is located at the larger of the smallest y value and Q1 – 1.5xIQR\n- The median is the horizontal line inside each box and the mean denoted by the cross\n -Individual patient profiles are denoted by dotted lines\n- A small amount of jitter is added to the data to aid visualisation.")+
                   
                       geom_text(data = dd %>% group_by( memorypar, tailindex) %>%
                                     summarise(Count = n()) %>%
                                     ungroup %>%
                                     mutate(hillest=min((d$hillest)) - 0.05 * diff(range((d$hillest)))),
                                 aes(label = paste0("n = ", Count)),
                                 position = pd, size=5, show.legend = FALSE) 
                   
                   #####################
                   print(pr1 + labs(y=target, x = "Visit") + 
                            # ggtitle(paste0("N=",length(unique(d$rep)),
                             #               " is the total number of patients; profiles with the number of patient values at visit") ) +
                          # ggtitle(paste0("There are N=",length(unique(d$rep)),  " patients with data at baseline, boxplots at each visit with all patient profiles and the number of patient values at each visit") ) +
                           
                           
                           ggtitle(paste0("There are N=",
                                          length(unique(d$rep)),  
                                          " patients with data at baseline, presenting selected patient profiles") ) +
                           theme_bw() +
                             theme(legend.position="none") +
                             theme(panel.background=element_blank(),
                                   # axis.text.y=element_blank(),
                                   # axis.ticks.y=element_blank(),
                                   # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
                                   # stop axis being clipped
                                   plot.title=element_text(size = 18), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
                                   legend.text=element_text(size=14),
                                   legend.title=element_text(size=14),
                                   legend.position="none",
                                   axis.text.x  = element_text(size=15),
                                   axis.text.y  = element_text(size=15),
                                   axis.line.x = element_line(color="black"),
                                   axis.line.y = element_line(color="black"),
                                   plot.caption=element_text(hjust = 0, size = 11),
                                   strip.text.x = element_text(size = 16, colour = "black", angle = 0),
                                   axis.title.y = element_text(size = rel(1.5), angle = 90),
                                   axis.title.x = element_text(size = rel(1.5), angle = 0),
                                   strip.background = element_rect(colour = "black", fill = "white")
                             ) 
                   )
                   
                }   else  if (input$Plot1 ==  "Individual all tests") {
                   
                   
                                  i <- as.numeric(unlist(strsplit(input$vec1,",")))
                                
                                  d <-d1
                                  
                                  dd <- d[d$rep %in% i[1],]
                                  
                                  sel <- unique(dd$tailindex)  
                                  
                                  d <- dd[dd$tailindex %in% sel,]
               
                                  require(lattice)
              
                                  lattice.options(panel.error=NULL)
                                  
                                  colnames(d)[colnames(d)=="hillest"] <- "value"
                                   
                                  colnames(d)[colnames(d)=="memorypar"] <- "Visit"
                                   
                                  xyplot(value ~ Visit | test, main=paste0( input$Plot ,"; all observed results for patient ", i[1],""), par.settings=list(par.main.text=list(cex=2)),
                                         par.strip.text=list(cex=.7),
                                         group = test, data = d,
                                           type = c("p" ,"l"),  scales = "free") 

               }
               
               
               
 
        })
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # diagnostics
    
    output$res.plot  <- renderPlot({       
    
          f <- fit.regression()
          fit <- f$fit.res
          
          # d <- make.data()$d1
          d <- f$z
          
          target <- input$Plot
          d <- d[d$test %in% target,]
          d2 <- d
          
          d2$resid <- r <- resid(fit)
          
          d2$fitted <- fitted(fit)
          
          yl <- ylab('Residuals') 
          
          xl <- xlab("time")
          
          p1 <- ggplot(d2 , aes(x=fitted , y=resid)) + geom_point () + yl 
          
          p3 <- ggplot(d2 , aes(x=time , y=resid)) +  geom_point () + yl  + xl +
              stat_summary(fun.data ="mean_sdl", geom='smooth') 
          
          p4 <- ggplot(d2 , aes(sample=resid)) + stat_qq() +
              geom_abline(intercept=mean(r), slope=sd(r)) + yl 
          
          gridExtra::grid.arrange(p1,  p3, p4, ncol=2) #
    
    })
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # listing of simulated data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    output$table1 <- DT::renderDataTable({
        
        foo<- make.data()$d1
        
        target <- input$Plot
        
        foo <- foo[foo$test %in% target,]
        
        foo$eij <- NULL 
        
        names(foo) <- c("Biochemistry test", "ID", "Visit", "Treatment","Response")
        
        rownames(foo) <- NULL
        
        
        library(DT)
        
        foo <-   foo %>%
          datatable(  ) %>%
          formatRound(
            columns= c("Biochemistry test", "ID", "Visit", "Treatment","Response"), digits=4)  
        
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # summary stats
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   # output$table2 <- renderPrint({
      output$table2 = DT::renderDataTable({
      
        foo<- make.data()$d1
        
        target <- input$Plot
        
        foo <- foo[foo$test %in% target,]
      
      
      
      
      f<-plyr::ddply(foo, c("test", "memorypar","tailindex"), summarise,
                  min=min(hillest),mean = mean(hillest), sd = sd(hillest, na.rm=TRUE),
                  sem = sd(hillest)/sqrt(length(hillest)),  Q1=quantile(hillest, 0.25)    , 
                  median=median(hillest),   Q3=quantile(hillest, 0.75)  , max=max(hillest)  )
      
        names(f) <- c("Biochemistry test",  "Visit", "Treatment","Minimum", "Mean" , "SD", "SE", "Q1","Median","Q3", "Maximum")
      
        rownames(f) <- NULL
     
        library(DT)

        ff <-   f %>%
        datatable(  ) %>%
        formatRound(
          columns=c("Minimum", "Mean" , "SD", "SE", "Q1","Median","Q3", "Maximum"), digits=2)  
        
       
        #              
                    
        # datatable(f,  options = list(
        #   columnDefs = list(list(className = 'dt-center', targets = 5)),
        #   
        #   pageLength = 5,
        #   lengthMenu = c(2, 10, 15, 20)
        # ), )
        
                    
      
    #  return(ff)
      
    })
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # estimate at specified time
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # output$reg.summary4 <- renderPrint({
    # 
    #     summary4 <- fit.est()$ff
    # 
    #     return(list(summary4))
    # 
    # })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # model output
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    output$reg.summary2 <- renderPrint({
        
        summary <- fit.regression()$fit.res
        
        return(list(summary))
        
    })  
    
    
    output$table4 <- renderPrint({
      
      f <- fit.regression()$target
      
      
        
        # dd <-  make.data2()$d 
        # 
        # rangez <-    input$Whisker 
        # 
        # bs <- boxplot.stats(dd$logy, coef = rangez, do.conf = FALSE, do.out = FALSE)$stats
        # bs <- as.numeric((as.vector(exp(bs))))  #p2
        # names(bs ) <- c("Lower whisker", "Lower ‘hinge’", "Median", "Upper ‘hinge’" ,"Upper whisker")
        # 
        # return(  print(bs, row.names = FALSE)) 
        
    })
    
    output$table5 <- renderPrint({
        
        # dd <-  make.data2()$d 
        # 
        # f <- summary(dd$logy)[c(2,3,5)]
        # f <- as.matrix(f);
        # # f <-p2(f)
        # f <- exp(as.numeric(f))
        # #  f<-(p2(f))
        # f<-as.data.frame(t(f))
        # # names(f ) <- c("Minimum", "1st.Quartile", "Median", "Mean", "3rd.Quartile", "Maximum")
        # names(f ) <- c( "1st.Quartile", "Median",  "3rd.Quartile")
        # return( print(f, row.names = FALSE)) 
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
})

# Run the application 
shinyApp(ui = ui, server = server)