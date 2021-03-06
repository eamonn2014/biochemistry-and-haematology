#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rshiny ideas from on https://gallery.shinyapps.io/multi_regression/
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
library(nlme)
library(VCA)
library(MASS)
require(tidyverse)
require(ggplot2)


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
                     
                      
                      
                      
                      
                      
                      
                      
                      
                      
                           
                      
                      
            Plot=biochemistry[1] #"1. Select which biochemistry test to present"),
            
            
            
            target <-  Plot               
             Plot1= "Overall" # c("Overall","Individual","Individual all tests")),
                      
             # strong("3. Select patient. If '2 select plot' 'Individual' is selected, enter sample ID(s) (comma delimited); 
             #                          enter 999 to show all profiles; If 'Individual all tests' is selected, all test results for the first ID only are presented"), "1,2,3,4"),

              vec1= c(1,2,3,4)
              
                      
             V=8#. Maximum visit number in data simulation including baseline"
                                
                      
            VV=4#. Estimate treatment effect at this visit"),
                                
 
    
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
    
    i=V
    
    # beta1 = slope
    # beta2 = tr effect
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
    
    d1$memorypar <- d1$memorypar - 1
    d1$memorypar <- factor(d1$memorypar)
    
    d1$tailindex <- factor(d1$tailindex)
    d1$tailindex <- ifelse(d1$tailindex %in% 1, "Active","Placebo" )
    
 

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    d <-  d1
    
    target <-  Plot
    
    d <- d[d$test %in% target,]
    
    d$time <- d$memorypar
    
    require(rms)
    
    d$time <- relevel(d$time, ref=VV)
    
    d$trt <- factor(d$tailindex) 
    
    d$yij <- d$hillest
    
    # new~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #convert time=1 to baseline var
    d <- d[, c("test", "rep", "time", "trt", "yij")]
    d$id=1:nrow(d)  # add index
    baseline <- d[d$time %in% 0,] # create baseline removing baseline
    names(baseline) <-  c("test", "rep", "time", "trt", "baseline", "id")
    baseline$id <- NULL
    baseline$time <- NULL
    d2 <- d[!d$time ==0,]         # create follow up
    both <- merge (baseline , d2   , all=TRUE)
    both$rep <- as.numeric(as.character(both$rep))
    both <- plyr::arrange(both, rep, time)
    both$time <- as.numeric(as.character(both$time))
    d <- both
    d$time<-factor(d$time)
    d$time <- relevel(d$time, ref=VV)
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
    
   
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # treatment effect estimate
   
    fit <- fit.res
    
    time. <- rep(1:(V-1))
    
    k1 <- contrast(fit, list(time=time.,  trt = 'Placebo'),
                   list(time=time.,  trt = 'Active'))
    
    
    x <- as.data.frame(k1[c('time', 'Contrast', 'Lower', 'Upper')]) 
    
    namez <- c("Follow-up Visit", "Placebo - Active estimate", "Lower 95%CI","Upper 95%CI")
    
    names(x) <- namez
    
    namez2 <- c("Placebo - Active estimate", "Lower 95%CI","Upper 95%CI")
    
    library(DT)
    
    #https://datatables.net/reference/option/
    datatable(x,   
              
              rownames = FALSE,
              
              options = list(
                searching = FALSE,
                pageLength = V-1,
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
        columns=c(namez2), digits=c(2)  )
    
    
      
  
  #---------------------------------------------------------------------------
  #---------------------------------------------------------------------------
  #---------------------------------------------------------------------------
  # Plot the estimated trt effect  
 
    time. <-  rep(1:(V-1))
    
    xl <- xlab( 'Follow up visit')
    
    k1 <- contrast(fit, list(time=time.,  trt = 'Placebo'),
                   list(time=time.,  trt = 'Active'))
    
    k1 <- as.data.frame(k1[c('time', 'Contrast', 'Lower', 'Upper')]) 
    
    mi <- floor(min(k1$Lower))
    ma <- ceiling(max(k1$Upper))
    
    ggplot (k1, aes(x=time. , y=Contrast, group=1)) + geom_point () + geom_line () +
      ylim(mi,ma) +
      #   xlim(1, input$V-1) +
      xlim(1, V) +
      scale_x_continuous(breaks=c(time.)) +
      #   ggpubr::grids(linetype = "dashed") +
      
      ylab( 'Placebo - Active')+ xl +
      geom_errorbar(aes(ymin=Lower, ymax=Upper ), width =0) +
      ggtitle(paste0("Outcome measure ", Plot ,"; treatment effect estimate at each follow up visit with 95% CI")) +
      geom_hline(aes(yintercept = 0, colour = 'red'), linetype="dashed") +
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
        plot.caption=element_text(hjust = 0, size = 7),
        strip.text.x = element_text(size = 16, colour = "black", angle = 0),
        axis.title.y = element_text(size = rel(1.5), angle = 90),
        axis.title.x = element_text(size = rel(1.5), angle = 0),
        panel.grid.major.x = element_line(color = "grey80", linetype="dotted", size = 1),
        panel.grid.major.y = element_line(color = "grey80", linetype="dotted", size = 1),
        strip.background = element_rect(colour = "black", fill = "#ececf0"),
        panel.background = element_rect(fill = '#ececf0', colour = '#ececf0'),
        plot.background = element_rect(fill = '#ececf0', colour = '#ececf0')
      )
    
 
  
  # --------------------------------------------------------------------------
  # -----------------------------------------------OVERALL PLOT
  # ---------------------------------------------------------------------------
           
    
    d <-  d1
   Plot1 = "Overall" 
      
      target <- Plot
      
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
      
      
      print(pr1 + labs(y=target, x = 'Visit (0 denotes the baseline visit)') +    
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
                strip.background = element_rect(colour = "black", fill = "#ececf0"),
                panel.background = element_rect(fill = '#ececf0', colour = '#ececf0'),
                plot.background = element_rect(fill = '#ececf0', colour = '#ececf0'),#
              ) 
      )
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Individual profiles
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
    Plot1 = "Individual"
      
  
        
        
        i <-  vec1
        
        target <-  Plot
        
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
        
        print(pr1 + labs(y=target, x = 'Visit (0 denotes the baseline visit)') + 
                
                ggtitle(paste0("There are N=",
                               length(unique(d$rep)),  
                               " patients with data at baseline, presenting selected patient profiles") ) +
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
                  strip.background = element_rect(colour = "black", fill = "#ececf0"),
                  panel.background = element_rect(fill = '#ececf0', colour = '#ececf0'),
                  plot.background = element_rect(fill = '#ececf0', colour = '#ececf0'),#
                ) 
        )
        
      Plot1 ==  "Individual all tests"
      
      
      i <- vec1
      
      
      
      d <-  d1
        
        dd <- d[d$rep %in% i[1],]
        
        sel <- unique(dd$tailindex)  
        
        d <- dd[dd$tailindex %in% sel,]
        
        library(lattice)
        
        lattice.options(panel.error=NULL)
        
        colnames(d)[colnames(d) %in% "hillest"] <- "value"
        
        colnames(d)[colnames(d) %in% "memorypar"] <- "Visit"
        
        xy <<- xyplot(value ~ Visit | test,
                      main=paste0( Plot ,"; all observed results for patient ", i[1]," allocated to ",sel,""), 
                      par.settings=list(par.main.text=list(cex=2)),
                      par.strip.text=list(cex=.7),
                      group = test, data = d,
                      xlab="Visit (0 denotes the baseline visit)",
                      type = c("p" ,"l"),  scales = "free") 
        
        print(xy)
        
        
    
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # diagnostics
  
    

    
    d <- z
    
    target <- Plot
    d <- d[d$test %in% target,]
    d2 <- d
    
    d2$resid <- r <- resid(fit)
    
    d2$fitted <- fitted(fit)
    
    yl <- ylab('Residuals') 
    
    xl <- xlab("time")
    
    p1 <- ggplot(d2 , aes(x=fitted , y=resid)) + geom_point (   colour="#69b3a2") + yl 
    
    p3 <- ggplot(d2 , aes(x=time , y=resid )) +  geom_point ( colour="#69b3a2") + yl  + xl +
      stat_summary(fun.data ="mean_sdl", geom='smooth') 
    
    p4 <- ggplot(d2 , aes(sample=resid )) + stat_qq(colour="#69b3a2") +
      geom_abline(intercept=mean(r), slope=sd(r)  ,  colour="black") + 
      xlab('Residuals')   +
      ggtitle( " ")  
    
    # p5 <- d2 %>%
    #   ggplot( aes(x=r)) +
    #   geom_histogram( fill="#69b3a2", color="#e9ecef", alpha=0.9) + #binwidth=1, 
    #  theme(
    #     plot.title = element_text(size=15)
    #   ) 
    library(gridExtra)
    library(grid)
    df <- data.frame(Residuals = r)
    p5 <- ggplot(df, aes(x = Residuals)) + 
      geom_histogram(aes(y =..density..),
                     #breaks = seq(-50, 50, by = 2), 
                     colour = "black", 
                     fill = "#69b3a2") +
      stat_function(fun = dnorm, args = list(mean = 0, sd = sigma(fit)  ))
    
    grid.arrange(p1,  p3, p4,p5, ncol=2,
                 
                 top = textGrob(paste0(Plot, " GLS model fit diagnostics"),gp=gpar(fontsize=20,font=3)))
    #+
    #main=paste0(input$Plot, "GLS model fit diagnostics")  #
    
 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  # listing of simulated data
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
    foo<-  d1
    
    target <-  Plot
    
    foo <- foo[foo$test %in% target,]
    
    foo$eij <- NULL 
    
    names(foo) <- c("Biochemistry test", "ID", "Visit", "Treatment","Response")
    
    rownames(foo) <- NULL
    library(DT)
    
    datatable(foo,   
              
              rownames = TRUE,
              
              options = list(
                searching = TRUE,
                pageLength = V-1,
                paging=FALSE,
                lengthMenu = FALSE ,
                lengthChange = FALSE,
                autoWidth = FALSE
                # colReorder = TRUE,
                # deferRender = TRUE,
                # scrollY = 200,
                # scroller = T
              ))  %>%
      formatRound(
        columns= c("Response"), digits=c(4)  )
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  # summary stats
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
 
    
    foo<-  d1
    
    target <-  Plot
    
    foo <- foo[foo$test %in% target,]
    
    f<-plyr::ddply(foo, c("test", "memorypar","tailindex"), summarise,
                   min=min(hillest),mean = mean(hillest), sd = sd(hillest, na.rm=TRUE),
                   sem = sd(hillest)/sqrt(length(hillest)),  Q1=quantile(hillest, 0.25)    , 
                   median=median(hillest),   Q3=quantile(hillest, 0.75)  , max=max(hillest)  )
    
    names(f) <- c("Biochemistry test",  "Visit", "Treatment","Minimum", "Mean" , "SD", "SE", "Q1","Median","Q3", "Maximum")
    
    rownames(f) <- NULL
    
    
    library(DT)
    datatable(f,   
              rownames = TRUE,
              options = list(
                searching = TRUE,
                pageLength = V-1,
                paging=FALSE,
                lengthMenu = FALSE ,
                lengthChange = FALSE,
                autoWidth = FALSE
                # colReorder = TRUE,
                # deferRender = TRUE,
                # scrollY = 200,
                # scroller = T
              ))  %>%
      formatRound(
        columns= c("Minimum", "Mean" , "SD", "SE", "Q1","Median","Q3", "Maximum"), digits=c(2)  )
    
    
  