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
options(max.print=1000000)
fig.width <- 1200
fig.height <- 550
library(shinythemes)        # more funky looking apps
p1 <- function(x) {formatC(x, format="f", digits=1)}
p2 <- function(x) {formatC(x, format="f", digits=2)}
options(width=100)
set.seed(12345) #reproducible

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
    "Hba1c glycosylated hemoglobin (%)",  
    "non fasting blood Glucose (mg/dL)" , 
    "Calcium mmol/L"  ,
    "Chloride mmol/L" ,
    "Creatine phosphokinase (CPK) U/L",  
    "Inorganic phophorous mmol/L",
    "Magnesium mmol/L"   ,
    "Potassium mmol/L"  ,
    "Sodium mmol/L" ,
    "Albumin g/L"   ,
    "Blood urea nitrogen BUN mmol/L"  ,
    "Creatinine (umol/L)" ,
    "Total protein g/L"  ,
    "Uric acid umol/L"  ,
    "Direct bilirubin umol/L",  
    "Gamma glutamyltransferase",  
    "Indirect bilirubin umol/L" ,
    "Total bilirubin umol/L" ,
    "Alkaline phosphatase U/L", 
    "High density lipoprotein mmol/L" , 
    "Low density lipoprotein mmol/L"  ,
    "Lipase (U/L)",
    "Total cholesterol mmol/L",  
    "Triglycerides mmol/L"  
)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ui <- fluidPage(theme = shinytheme("journal"),
                
                shinyUI(pageWithSidebar(
                    
                    #ui <-shinyUI(pageWithSidebar(
                    
                    headerPanel("Presenting biochemistry and haematology data"),
                    
                    #sidebarLayout(  #new
                    # Sidebar with a slider and selection inputs
                    
                    sidebarPanel( 
                        
                    div(p("Typically reams and reams of outputs are generated (but thankfully not printed) consisting of summary statistics,
                          perhaps overtime, for many biochemistry and haemotology parameters. How useful this is is open to debate,
                          but in my experience these outputs are hardly looked at, perhaps one or two parameters are singled out and 
                        scrutinised. This is not suprising because the presentation is very dry. There is the issue also that it is
                          quite possible that a different data set can give rise to similar summary statistics [1]. Spotting errors 
                          too in the data is not easy with this kind of presentation. I would argue that these parameters should always be 
                          plotted. RShiny is an ideal medium to communicate this information in a more informative and exciting way. 
                          I will make an attempt using simulated data.")),
                        
                        div(
                            
                            
                            actionButton(inputId='ab1', label="R code",   icon = icon("th"), 
                                         onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/Boxplots/master/app.R', '_blank')"),   
                            actionButton("resample", "Simulate a new sample"),
                            br(), br(),
                            
                            
                            div(strong("Select the parameters using the sliders below"),p(" ")),
                            
                            
                            div(("Both the number of data points and the length of the whiskers can be varied. 
                 The default whisker length uses the common 1.5xIQR rule. Selecting zero will 
                 extend the whiskers to the extreme datapoints.
                 There is the option to show all the data on the plot 'Show me the data!'. 
                 Could you guess what the data would look like?
                 Unlikely, it is possible the same boxplot can be constructed from many different data sets [3].
                 Therefore if possible and feasible also plot all the data points. There is also the option
                 to 'Highlight 'outliers'', that is present the 'outliers'. Be careful as plotting the individual data points also will double up the 'outlying'
                 data points. Therefore if you are programming boxplots and showing the individual data, turn off the 
                 outlier option. This advice is also pertinent if generating boxplots using the ggplot2 package [4].
                 Another sample can be taken from the same data generating mechanism by clicking 'Simulate a new sample'.")),
                            br(),
                            
                            selectInput("Plot",
                                        strong("Select plot preference "),
                                        choices=biochemistry),
                            
                            textInput('vec1', 'Enter sample id (comma delimited)', "1,2,3,4"),
                            
                            sliderInput("N",
                                        "Select the total number of data points",
                                        min=3, max=500, step=1, value=100, ticks=FALSE),
                            
                            sliderInput("Whisker",
                                        "Select the length of whiskers (multiples of the IQR)",
                                        min=0, max=3, step=.5, value=1.5, ticks=TRUE),
                            
                            sliderInput("outliers",
                                        "Highlight 'outliers'",
                                        min=0, max=1, step=1, value=0, ticks=FALSE),
                            
                            sliderInput("dp",
                                        "Show me the data!",
                                        min=0, max=1, step=1, value=1, ticks=FALSE),
                            
                            div(p("References:")),  
                            
                            tags$a(href = "https://en.wikipedia.org/wiki/Exploratory_data_analysis", "[1] Tukey, J.W. 'Exploratory Data Analysis'"),
                            div(p(" ")),
                            tags$a(href = "https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/boxplot", "[2] R boxplot"),
                            div(p(" ")),
                            tags$a(href = "https://en.wikipedia.org/wiki/Anscombe%27s_quartet", "[3] Anscombe's quartet"),
                            div(p(" ")),
                            tags$a(href = "https://ggplot2.tidyverse.org/reference/geom_boxplot.html", "[4] Boxplots using ggplot2"),
                            div(p(" ")),
                            tags$a(href = "https://r.789695.n4.nabble.com/Whiskers-on-the-default-boxplot-graphics-td2195503.html", "[5] Whiskers on the default boxplot"),
                            div(p(" ")),
                            
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
                            .navbar { background-color: lightgrey;}
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
                                     
                                     div(plotOutput("reg.plot", width=fig.width, height=fig.height)),  
                                     
                                     
                                     h3("Figure 1 Left panel untransformed data, right panel natural log transformation labelled with antilogs"),
                                     
                                     p(strong("Boxplots are simple graphical characterisations of continuous variables. 
        Boxplots were proposed by Tukey, J.W. in his 1977 book 'Exploratory Data Analysis' (only in 1977!).
                  In this example, we re-express the data using a natural logarithmic transformation.
                  Next, prior to presentation, perform the boxplot calculations. So first, perform the natural
                  log transformation on the data. Secondly create a boxplot, 
                  that is calculate the median, the two hinges and the whiskers. 
                  Then we present the boxplot, finally replace the axis log values with the antilog values.")) ,
                                     
                                     
                                     p(strong("There are subtleties in the boxplot calculations.
                  Often you may hear it is said, 'the limits of the boxes represent Q1 and Q3 and the 
                  whiskers extend +/- 1.5 x Q3-Q1'. 
                  This is not necessarily true. There are many variations on the box plot, it is therefore better to explicitly state what is being presented.")) ,
                                     div(""),
                                     p(strong("I hope you agree the plot on the right gives a better understanding of the data distributions. 
                 In this programming exercise select 'Show me the data!' and deselect 'Highlight 'outliers''.")) ,
                                     
                            ) ,
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            tabPanel("Wait...what? Aren't the whiskers different?", value=3, 
                                     div(plotOutput("reg.plot2", width=fig.width, height=fig.height)),  
                                     h3("Figure 2 Top panel untransformed data, bottom panel using a natural log transformation"),
                                     
                                     h4("Lets first look at the relationship between the hinges vrs. quartiles.") ,
                                     p(strong("Boxplot stats for the untransformed data, as presented in the top panel:")) ,
                                     div( verbatimTextOutput("table2")),
                                     p(strong("Now summarise the same data. Q1 and Q3 may not match the hinges. Paraphrasing the 'box.plot.stats' help file...'The ‘hinges’ are versions of Q1 and Q3'. See what happens when n is odd. The hinges and the quartiles match.")) ,
                                     div( verbatimTextOutput("table3")),
                                     
                                     p(strong("Now check out the boxplot stats based on the log transformed data, after exponentiating (bottom panel), check with the raw median:")) ,
                                     div( verbatimTextOutput("table4")),
                                     p(strong("In the same way summarise the log transformed data then exponentiating: ")) ,
                                     div( verbatimTextOutput("table5")),
                                     
                                     h4("Look again at how to calculate the boxplot statistics when transforming") ,
                                     
                                     p(strong("The length of the whiskers is typically set at 1.5xIQR, that is, the 
                 whiskers extend to the most extreme data point which is no more than 1.5 times the interquartile 
                 range from the box. A value of zero, using the R function, causes the whiskers to extend to the data extremes. 
                          So it may not be as simple as saying the whiskers extend 1.5xIQR.")),
                                     
                                     
                                     div(""),
                                     
                                     p(strong("To summarise, apply the whisker calculation rule on the scale used to draw the boxplot. Therefore when using a transformation 
                 the calculations must 
                 be on the transformed scale. Do not mix scales. So, do not calculate the whiskers based on the 1.5xIQR rule
                 on the raw scale and then log transform 
                          to present the boxplot.  
                          The median and hinges will be the logs of the original median and hinges, they will match with odd sample sizes, but not generally with even. 
                          The step which determines the length of the whiskers will change. Hence the fences generally will differ. 
                          (The location of the end of the whiskers is referred to as the fences). ")) ,
                                     
                                     
                                     p(strong("")),
                                     h3("Take home messages; when the sample size is even, the hinges may not necessarily equal Q1 and Q3 (raw and transformed) and the transformed median may be slightly different to the raw median.
                      The fences may not necessarily equal the fences on the raw scale for any sample size.
                    Explicitly state how you constructed your boxplots.")
                            ) ,
                            
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            tabPanel("The data", value=3, 
                                     div( verbatimTextOutput("table1")),
                            ) ,
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            tabPanel(" ", 
                            ) ,
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            tabPanel(" ", 
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
        n<-  sample$n  
        
        y <- rlnorm(n, .7, 1.5) 
        x <- factor(sample(3, length(y), repl = TRUE))
        
        d <- data.frame(x=x, y=y)
        d$logy <- log(d$y) # log the data
        
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
        
        z <- mdata(n=99,beta0=100, beta1=0, beta2=-5, ar.val=.9, sigma=5, tau0=7, tau1=1, tau01=0.01, m=6)
        df1 <- cbind(test=biochemistry[1], z)

        z <- mdata(n=85,beta0=5.6, beta1=0, beta2=-.4, ar.val=.9, sigma=.02, tau0=.2, tau1=.02, tau01=0.01, m=6)
        df2 <- cbind(test=biochemistry[2], z)

        z <- mdata(n=102,beta0=110, beta1=0, beta2=-5, ar.val=.9, sigma=5, tau0=.5, tau1=.5, tau01=0.01, m=6)
        df3 <- cbind(test=biochemistry[3], z)

        z <- mdata(n=102,beta0=2.4, beta1=.1, beta2=0, ar.val=.9, sigma=.2, tau0=.5, tau1=.01, tau01=0.2, m=6)
        df4 <- cbind(test=biochemistry[4], z)

        z <- mdata(n=90,beta0=104, beta1=1, beta2=2, ar.val=.8, sigma=4, tau0=5, tau1=1, tau01=.6, m=6)
        df5 <- cbind(test=biochemistry[5], z)

        z <- mdata(n=90,beta0=95, beta1=0, beta2=2, ar.val=.8, sigma=4, tau0=15, tau1=1, tau01=.6, m=6)
        df6 <- cbind(test=biochemistry[6], z)

        z <- mdata(n=90,beta0=1.23, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df7 <- cbind(test=biochemistry[7], z)


        z <- mdata(n=90,beta0=1, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df8 <- cbind(test=biochemistry[8], z)


        z <- mdata(n=99,beta0=4, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df9 <- cbind(test=biochemistry[9], z)

        z <- mdata(n=99,beta0=141, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df10 <- cbind(test=biochemistry[10], z)

        z <- mdata(n=111,beta0=44, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df11 <- cbind(test=biochemistry[11], z)

        z <- mdata(n=88,beta0=5, beta1=0, beta2=.05, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df12 <- cbind(test=biochemistry[12], z)

        z <- mdata(n=112,beta0=70, beta1=0, beta2=.5, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df13 <- cbind(test=biochemistry[13], z)

        z <- mdata(n=112,beta0=72, beta1=0, beta2=.5, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df14 <- cbind(test=biochemistry[14], z)

        z <- mdata(n=97,beta0=3.4, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df15 <- cbind(test=biochemistry[15], z)

        z <- mdata(n=103,beta0=24, beta1=0, beta2=1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df16 <- cbind(test=biochemistry[16], z)

        z <- mdata(n=105,beta0=6, beta1=0, beta2=.5, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df18 <- cbind(test=biochemistry[18], z)

        z <- mdata(n=100,beta0=9, beta1=0, beta2=.0, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df19<- cbind(test=biochemistry[19], z)

        z <- mdata(n=103,beta0=64, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df20<- cbind(test=biochemistry[20], z)

        z <- mdata(n=97,beta0=1.3, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df21<- cbind(test=biochemistry[21], z)

        z <- mdata(n=98,beta0=3.3, beta1=0, beta2=-.2, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df22<- cbind(test=biochemistry[22], z)

        z <- mdata(n=104,beta0=34, beta1=0, beta2=.1, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df23<- cbind(test=biochemistry[23], z)

        z <- mdata(n=99,beta0=5.11, beta1=0, beta2=.0, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df24<- cbind(test=biochemistry[24], z)

        z <- mdata(n=100,beta0=1.23, beta1=0, beta2=.0, ar.val=.8, sigma=.1, tau0=.10, tau1=.01, tau01=.6, m=6)
        df17<- cbind(test=biochemistry[17], z)



        d1 <- do.call(rbind, list(df1 = df1, df2 = df2, df3=df3, df4=df4, df5=df5, df6=df6, df7=df7,
                                  df8=df8,df9=df9,df10=df10,df11=df11,df12=df12,df13=df13,df14=df14,
                                  df15=df15, df16=df16,df17=df17, df18=df18, df19=df19, df20=df20,
                                  df21=df21,df22=df22, df23=df23, df24=df24))



        colnames(d1)[colnames(d1)=="yij"] <- "hillest"
        colnames(d1)[colnames(d1)=="trt"] <- "tailindex"
        colnames(d1)[colnames(d1)=="obs"] <- "memorypar"
        colnames(d1)[colnames(d1)=="id"] <- "rep"
        d1$memorypar <- factor(d1$memorypar)
        d1$tailindex <- factor(d1$tailindex)
        d1$tailindex <- ifelse(d1$tailindex %in% 1, "Active","Placebo" )
        
        
          
        return(list(  d1=d1))# rangez=rangez, outliers=outliers, n=n, dp=dp))
        
        
    })
    
    make.data2 <- reactive({
        
        sample <- random.sample()
        rangez <-    sample$Whisker       # multiples of IQR length of whiskers, 0 means out to maximum
        n<-          sample$n  
        
        y<- c(rbeta(n-4, 2,6)*35,  sample(50:150,4, replace=FALSE))
        x <- factor(sample(3, length(y), repl = TRUE))
        
        d <- data.frame(x=x, y=y)
        d$logy <- log(d$y) # log the data
        
        return(list(d=d , y=d$y, rangez=rangez))# 
        
    })
    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # ---------------------------------------------------------------------------
    
           output$reg.plot3 <- renderPlot({         
            
            target <- input$Plot
            
            d <- make.data()$d1
             
            d <- d[d$test %in% target,]
            
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
                
                geom_text(data = d %>% group_by( memorypar, tailindex) %>%
                              summarise(Count = n()) %>%
                              ungroup %>%
                              mutate(hillest=min((d$hillest)) - 0.05 * diff(range((d$hillest)))),
                          aes(label = paste0("n = ", Count)),
                          position = pd, size=3, show.legend = FALSE) 
            
            #####################
            print(pr1 + labs(y=target, x = "Visit") + 
                      ggtitle(paste0("N=",length(unique(d$rep))," patient profiles with the number of patient values at visit") ) +
                      theme_bw() +
                      theme(legend.position="none") 
            )
           
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
        })
    
    # Plot a scatter of the data  
    
    output$reg.plot <- renderPlot({         
        
        
        
        d <- make.data()$d1
        
        
        ####################################SINGLE PLOT
       
        i <- as.numeric(unlist(strsplit(input$vec1,",")))
       
        target <- input$Plot
        
        d <- d[d$test %in% target,]
    
        dd <- d[d$rep %in% i,]
        
        sel <- unique(dd$tailindex)
        
        d <- d[d$tailindex %in% sel,]
        
        
        pd <- position_dodge(.4)
        pr1=NULL
        pr1<-ggplot(d,aes(x=memorypar ,y=hillest,color=tailindex, fill=tailindex )) + 
            #stat_boxplot(geom = "errorbar", width = 0.3) +
            geom_boxplot( outlier.colour = NA,alpha=0.1, color="lightblue",)  +  
            geom_boxplot(data = d,
                         aes(x = memorypar, y = hillest,  fill = tailindex ),outlier.shape = NA  , alpha=.2 ) + 
            
            # facet_wrap(~tailindex , ncol=2)    +
            geom_line(data = dd,
                      aes(group=rep,x = memorypar, y = hillest),  size = .6, linetype="dashed") +
            
            scale_size_manual( values = c( 1) ) +
            # geom_point(aes(fill=tailindex, group=rep), pch=1, size=1, alpha=0.3, position = pd ) +
            # stat_summary(fun.y=mean, geom="point", shape=3, size=2, colour="black", stroke=1.5,
            #             position=pd, show.legend=FALSE) +
            scale_color_manual(name = "Treatment", values = c("blue", "darkgreen") ) +
            scale_fill_manual(name = "Treatment", values = c("lightblue", "green") ) +
            
            
            facet_wrap(~tailindex , ncol=2)    +
            
            geom_text(data = dd %>% group_by( memorypar, tailindex) %>%
                          summarise(Count = n()) %>%
                          ungroup %>%
                          mutate(hillest=min((d$hillest)) - 0.05 * diff(range((d$hillest)))),
                      aes(label = paste0("n = ", Count)),
                      position = pd, size=3, show.legend = FALSE) 
        
        #####################
        print(pr1 + labs(y=target, x = "Visit") + 
                  ggtitle(paste0("N=",length(unique(dd$rep))," patient profiles with the number of patient values at visit") ) +
                  theme_bw() +
                  theme(legend.position="none") 
        )
        
        
        
        
    })  
    #---------------------------------------------------------------------------
    
    output$reg.plot2 <- renderPlot({         
        
       
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # listing of simulated data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    output$table1 <- renderPrint({
        
        return(make.data2()$d)
        
    })
    
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    output$table2 <- renderPrint({
        
        dd <-  make.data2()$d 
        
        rangez <-    input$Whisker 
        
        bs <- boxplot.stats(dd$y, coef = rangez, do.conf = FALSE, do.out = FALSE)$stats
        bs <- as.numeric((as.vector(bs))) #p2
        names(bs ) <- c("Lower whisker", "Lower ‘hinge’", "Median", "Upper ‘hinge’" ,"Upper whisker")
        
        return(  print(bs, row.names = FALSE)) 
        
    })
    
    output$table3 <- renderPrint({
        
        dd <-  make.data2()$d 
        
        f <- summary(dd$y)[c(2,3,5)]
        f <- as.matrix(f);
        #f <-p2(f)   # here
        f <- as.numeric(f)
        # f<-(p2(f))
        f<-as.data.frame(t(f))
        #    names(f ) <- c("Minimum", "1st.Quartile", "Median", "Mean", "3rd.Quartile", "Maximum")
        names(f ) <- c( "1st.Quartile", "Median",  "3rd.Quartile")
        
        return( print(f, row.names = FALSE)) 
        
    })
    
    
    output$table4 <- renderPrint({
        
        dd <-  make.data2()$d 
        
        rangez <-    input$Whisker 
        
        bs <- boxplot.stats(dd$logy, coef = rangez, do.conf = FALSE, do.out = FALSE)$stats
        bs <- as.numeric((as.vector(exp(bs))))  #p2
        names(bs ) <- c("Lower whisker", "Lower ‘hinge’", "Median", "Upper ‘hinge’" ,"Upper whisker")
        
        return(  print(bs, row.names = FALSE)) 
        
    })
    
    output$table5 <- renderPrint({
        
        dd <-  make.data2()$d 
        
        f <- summary(dd$logy)[c(2,3,5)]
        f <- as.matrix(f);
        # f <-p2(f)
        f <- exp(as.numeric(f))
        #  f<-(p2(f))
        f<-as.data.frame(t(f))
        # names(f ) <- c("Minimum", "1st.Quartile", "Median", "Mean", "3rd.Quartile", "Maximum")
        names(f ) <- c( "1st.Quartile", "Median",  "3rd.Quartile")
        return( print(f, row.names = FALSE)) 
        
    })
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
})

# Run the application 
shinyApp(ui = ui, server = server)