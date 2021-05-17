library(shiny)
library(dplyr)
library(readr)
library(DT)
library(pls)
#library(knitr)
#library(mdatools)
library(factoextra)
library(extrafont)
library(tidyverse)
library(ggplot2)


server <- function(input, output, session) {
  
  #_____________________________________________________________________________________
  #generic functions
  
  #increase file size uploads
  options(shiny.maxRequestSize=30*1024^2)
  
  #stop app when closing
  session$onSessionEnded(function() {
    stopApp()
  })
  
  #create SNV function
  SNV<-function(spectra){                                             
    spectra<-as.matrix(spectra)
    spectrat<-t(spectra)
    spectrat_snv<-scale(spectrat,center=TRUE,scale=TRUE)
    spectra_snv<-t(spectrat_snv)
    return(spectra_snv)
  }
  
  #get uploaded file name function
  file_name <- reactive({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    return (stringi::stri_extract_first(str = inFile$name, regex = ".*(?=\\.)"))
  })
  
  #read in data function
  getData <- reactive({
    df <- read.csv(input$file1$datapath)
    
    #tidy: remove unnamed columns, rename, lowercase values
    df <- dplyr::select(df, -X) #remove duplicate Image column without header
    df <- dplyr::select(df, -contains("Count")) #remove PAD Count
    df <- dplyr::select(df, -contains("Amox")) #remove Amox %
    df <- rename(df, DrugPct = contains("Drug")) #rename Drug % column
    df <- rename(df, PAD = contains("PAD")) #rename PAD S# column
    df <- rename(df, First = contains("First")) #rename First 20 column
    df <- rename(df, Mid = contains("Mid")) #rename Mid 20 column
    df <- rename(df, Last = contains("Last")) #rename Last 20 column
    df$Contains <- tolower(df$Contains) #lowercase all values in Contains column
    df$Contains <- ifelse(df$Contains=="ascorbicacid","ascorbic acid",df$Contains)
    df$Contains <- ifelse(df$Contains=="ferroussulfate","ferrous sulfate",df$Contains)
    df$Contains <- ifelse(df$Contains=="promethazinehydrochloride","promethazine hydrochloride",df$Contains)
    
    
    #create api variable; use Drug % column values for API of analysis; 0 all else
    df <- df %>% mutate(api = ifelse(df$Contains==input$selectAPI,df$DrugPct,0))
    
    
    set <- dplyr::select(df, !!input$images)
    
    df <- cbind(df, ImageUse=set[,1])
    
  })
  
  
  #_____________________________________________________________________________________  
  #Predict tab
  #getPredict() - reactive function for predicting test data APIs
  getPredict <- reactive({
    model <- getModel()[[1]] #retrieve PLSR model
    df <- getData() #retrieve original data to get testing dataset
    
    #lane selection
    df1 <- dplyr::select(df, -matches("[[:digit:]]\\.")) #exclude any color data (include ID columns and controls)
    if(input$controls == TRUE) {
      df1 <- df1 #include control columns (non-computer data only)
    } else {
      df1 <- dplyr::select(df1, -matches(".r$|.g$|.b$")) #remove control columns (non-computer data only)
    }
    
    #select color data
    df2 <- dplyr::select(df, matches("[[:digit:]]\\.")) #color data only
    lanesVar <- paste(unlist(input$laneBands),collapse=",") #capture lanes selected for analysis
    lanes <- paste0("^(",paste(unlist(input$laneBands),collapse = "|"),")")
    df2 <- dplyr::select(df2, matches(lanes)) #filter color data by selected lanes
    
    df <- cbind(df1,df2) #rejoin ID/control dataset with color dataset
    
    #divide data set into training and testing sets
    
    #ImageUse <- input$images
    
    training <- filter(df, ImageUse=="Training")
    
    test <- filter(df, ImageUse=="Test")
    
    
    #allow for manually selecting samples for training model
    if(input$rowselect == 1) {
      rows <- input$contents_rows_selected
      training <- training[rows,]
      test <- rbind(test, test[-rows,])
    }
    
    #samps, contains, PAD: get identifying data
    samps <- as.vector(dplyr::select(test, Image))
    
    contains <- as.vector(dplyr::select(test, Contains))
    
    PAD <- as.vector(dplyr::select(test, PAD))
    
    #x = independent variables for testing model
    x <- as.matrix(dplyr::select(test, matches(".r$|.g$|.b$|[[:digit:]]\\."))) #include RGB
    
    #scaling
    #x <- scale(x)
    
    #MSC
    x <- if(input$msc == TRUE) {
      msc(x)
    } else {
      x
    }
    
    #SNV
    x <- if(input$snv == TRUE) {
      SNV(x)
    } else {
      x
    }
    
    y <- as.vector(t(dplyr::select(test, matches("api"))))
    
    y[is.na(y)] <- 0
    
    x2 <- data.frame(yy=I(data.matrix(y)), xx=I(data.matrix(x)))
    
    preds <- predict(model, ncomp = input$components, newdata = x2)
    
    error <- RMSEP(model, newdata=x2, ncomp=input$components)
    
    error <- error$val[2]
    
    model_error <- RMSEP(model, ncomp=input$components)
    
    model_error <- model_error$val[3]
    
    params <- paste(input$algorithm, input$images,
                    ifelse(input$center==TRUE,"Centered", "No Centering"),
                    ifelse(input$scale==TRUE,"Scaled", "No Scaling"), 
                    ifelse(input$msc==TRUE, "MSC", "No MSC"),
                    ifelse(input$snv==TRUE, "SNV", "No SNV"),sep=", "
    )
    
    pred_table <- cbind(Image=samps, PAD=PAD, Contains=contains, 
                        API=as.numeric(y),
                        #API=as.numeric(exp(y)),
                        Predicted.API=as.numeric(as.vector(preds)), 
                        Error=as.numeric((abs(y-as.vector(preds)))), Test_RMSEP=error, 
                        CV_RMSEP=model_error, 
                        Delta=abs(model_error-error),
                        #RMSE=sqrt(mean((y-as.vector(preds))^2)),
                        Factors=input$components, Lanes=lanesVar, Parameters=params, file=file_name()
    )
    
    pred_table <- as.data.frame(pred_table, stringsAsFactors=FALSE)
    
    #RMSEP1 <- as.data.frame(RMSEP(model, newdata=x, ncomp=input$components)) 
    
    return(pred_table)
    #return(x)
    
  }) #end predict
  
  #Predict values output
  output$predictvalues <- renderDT({
    getPredict()
  })
  
  #data download handler
  output$data1Download <- downloadHandler(
    filename = function() {
      paste0("PLS_PREDICT_",toupper(input$selectAPI), ".csv")
    },
    content = function(file) {
      write.csv(getPredict(), file, row.names = FALSE)
    }
  ) #end data download handler
  
  #Find Min. RMSEP button
  observeEvent(input$lowRMSEP, {
    model <- getModel()[[2]]
    comps <- dim(model)[2]
    rmseps <- c()
    for(i in 1:comps){
      model <- getModel()[[1]] #retrieve PLSR model
      df <- getData() #retrieve original data to get testing dataset
      
      #lane selection
      df1 <- dplyr::select(df, -matches("[[:digit:]]\\.")) #exclude any color data (include ID columns and controls)
      if(input$controls == TRUE) {
        df1 <- df1 #include control columns (non-computer data only)
      } else {
        df1 <- dplyr::select(df1, -matches(".r$|.g$|.b$")) #remove control columns (non-computer data only)
      }
      
      #select color data
      df2 <- dplyr::select(df, matches("[[:digit:]]\\.")) #color data only
      lanesVar <- paste(unlist(input$laneBands),collapse=",") #capture lanes selected for analysis
      lanes <- paste0("^(",paste(unlist(input$laneBands),collapse = "|"),")")
      df2 <- dplyr::select(df2, matches(lanes)) #filter color data by selected lanes
      
      df <- cbind(df1,df2) #rejoin ID/control dataset with color dataset
      
      #divide data set into training and testing sets
      
      #ImageUse <- input$images
      
      training <- filter(df, ImageUse=="Training")
      
      test <- filter(df, ImageUse=="Test")
      
      
      #allow for manually selecting samples for training model
      if(input$rowselect == 1) {
        rows <- input$contents_rows_selected
        training <- training[rows,]
        test <- rbind(test, test[-rows,])
      }
      
      #samps, contains, PAD: get identifying data
      samps <- as.vector(dplyr::select(test, Image))
      
      contains <- as.vector(dplyr::select(test, Contains))
      
      PAD <- as.vector(dplyr::select(test, PAD))
      
      #x = independent variables for testing model
      x <- as.matrix(dplyr::select(test, matches(".r$|.g$|.b$|[[:digit:]]\\."))) #include RGB
      
      #MSC
      x <- if(input$msc == TRUE) {
        msc(x)
      } else {
        x
      }
      
      #SNV
      x <- if(input$snv == TRUE) {
        SNV(x)
      } else {
        x
      }
      
      y <- as.vector(t(dplyr::select(test, matches("api"))))
      
      y[is.na(y)] <- 0
      
      x2 <- data.frame(yy=I(data.matrix(y)), xx=I(data.matrix(x)))
      
      preds <- predict(model, ncomp = input$components, newdata = x2)
      
      error <- RMSEP(model, newdata=x2, ncomp=i)
      
      error <- error$val[2]
      
      rmseps[i] <- error
      
    }
    output$minRMSEP <- renderPrint(paste0(which.min(rmseps)," components | RMSEP=",min(rmseps)))
  })
  
  #_____________________________________________________________________________________  
  #Model tab
  #getModel() - reactive function for PLSR modeling
  getModel <- reactive({
    
    df <- getData()
    
    ctr <- if(input$center == TRUE) {
      TRUE
    } else {
      FALSE
    }
    
    scl <- if(input$scale == TRUE) {
      TRUE
    } else {
      FALSE
    }
    
    #leave one out CV
    val <- if(input$validation == TRUE) {
      "LOO"
    } else {
      "none"
    }
    
    #lane selection
    df1 <- dplyr::select(df, -matches("[[:digit:]]\\."))
    if(input$controls == TRUE) {
      df1 <- df1
    } else {
      df1 <- dplyr::select(df1, -matches(".r$|.g$|.b$"))
    }
    
    #select color data
    df2 <- dplyr::select(df, matches("[[:digit:]]\\."))
    lanesVar <- paste(unlist(input$laneBands),collapse=",")
    lanes <- paste0("^(",paste(unlist(input$laneBands),collapse = "|"),")")
    df2 <- dplyr::select(df2, matches(lanes))
    
    df <- cbind(df1,df2)
    
    #divide data set
    
    #ImageUse <- input$images
    
    training <- filter(df, ImageUse=="Training")
    
    test <- filter(df, ImageUse=="Test")
    
    if(input$rowselect == 1) {
      rows <- input$contents_rows_selected
      training <- training[rows,]
      test <- rbind(test, training[-rows,])
    }
    
    #get x and y
    x <- as.matrix(dplyr::select(training, matches(".r$|.g$|.b$|[[:digit:]]\\.")))
    
    y <- as.vector(t(dplyr::select(training, matches("api"))))
    
    #MSC
    x <- if(input$msc == TRUE) {
      msc(x)
    } else {
      x
    }
    
    #SNV
    x <- if(input$snv == TRUE) {
      SNV(x)
    } else {
      x
    }
    
    
    training2 <- data.frame(yy=I(y), xx=I(x))
    
    model <- plsr(yy~xx, method=input$algorithm, data=training2, scale=scl, center=ctr, validation=val)
    #model <- plsr(yy~xx, method=input$algorithm, data=training2, validation=val)
    
    return(list(model, x))
    
  }) #end getModel reactive
  
  #PLSR visuals
  
  #model output
  output$modelCoefs <- renderPrint({
    model <- getModel()[[1]]
    coefs <- coef(model,ncomp=input$components,intercept=TRUE)
    coefs_int <- if(input$scale==TRUE) {
      int <- coefs[1]
      names(int) <- "(Intercept)"
      coefs2 <- coefs[-1]/model$scale
      c(int,coefs2)
    } else {
      coefs
    }
    coefs_int
  }) #end output
  
  #data download handler
  output$downloadModel <- downloadHandler(
    filename = function() {
      paste0("PLS_COEFS_",toupper(input$selectAPI),"_",toupper(file_name()),"_",toupper(input$images),"_",input$components,"_FACTORS",".csv")
    },
    content = function(file) {
      model <- getModel()[[1]]
      #scale <- getModel()[[3]]
      coefs <- coef(model,ncomp=input$components,intercept=TRUE)
      coefs_int <- if(input$scale==TRUE) {
        int <- coefs[1]
        names(int) <- "(Intercept)"
        coefs2 <- coefs[-1]/model$scale
        c(int,coefs2)
      } else {
        coefs
      }
      write.csv(coefs_int, file, row.names = TRUE)
    }
  ) #end data download handler
  
  #model output
  output$plsSummary <- renderPrint({
    summary(getModel()[[1]]) #summary
    #coef(getModel()[[1]]) #coefficients
    #fitted(getModel()[[1]]) #fitted values
    #residuals(getModel()[[1]]) #error
  }) #end mdist output
  
  #model output
  output$plsScores <- renderDT({
    scs <- scores(getModel()[[1]])
    scs <- scs[,]
    as.data.frame(scs)
  }) #end output
  
  #get scores matrix
  getScores <- reactive({
    sc <- scores(getModel()[[1]])
    sc <- sc[,]
    sc
  })
  
  #data download handler
  output$downloadScores <- downloadHandler(
    filename = function() {
      paste0("PLS_SCORES", ".csv")
    },
    content = function(file) {
      write.csv(getScores(), file)
    }
  ) #end data download handler
  
  #model output
  output$plsLoadings <- renderDT({
    load <- loadings(getModel()[[1]])
    load <- load[,]
    load
  }) #end output
  
  #data download handler
  output$downloadLoadings <- downloadHandler(
    filename = function() {
      paste0("PLS_LOADINGS", ".csv")
    },
    content = function(file) {
      write.csv(loadings(getModel()[[1]]), file, row.names = TRUE)
    }
  ) #end data download handler
  
  #PLSR biplot output
  output$biplot <- renderPlot({
    #plot(RMSEP(getModel()[[1]]))
    corrplot(getModel()[[1]], labels="names")
  }) #end output
  
  #plot download handler
  output$plot1 <- downloadHandler(
    filename = function() {
      paste0("BIPLOT", ".png")
    },
    content = function(file) {
      png(file, width = 1200)
      corrplot(getModel()[[1]], labels="names")
      dev.off()
    }
  ) #end plot download handler
  
  #PLSR loadings plot output
  output$loadplot <- renderPlot({
    #plot(RMSEP(getModel()[[1]]))
    plot(getModel()[[1]], "loadings", comps = 1:2, legendpos = "topleft",
         labels = "names")
    abline(h = 0)
  }) #end output
  
  #plot download handler
  output$plot2 <- downloadHandler(
    filename = function() {
      paste0("LOADINGSPLOT", ".png")
    },
    content = function(file) {
      png(file, width = 1200)
      plot(getModel()[[1]], "loadings", comps = 1:2, legendpos = "topleft",
           labels = "names")
      dev.off()
    }
  ) #end plot download handler
  
  #PLSR regression coefficient plot output
  output$coefplot <- renderPlot({
    plot(getModel()[[1]], plottype = "coef", ncomp=1:2, legendpos="bottomleft", labels = "names")
  }) #end output
  
  #plot download handler
  output$plot3 <- downloadHandler(
    filename = function() {
      paste0("COEFPLOT", ".png")
    },
    content = function(file) {
      png(file, width = 1200)
      plot(getModel()[[1]], plottype = "coef", ncomp=1:2, legendpos="bottomleft", labels = "names")
      dev.off()
    }
  ) #end plot download handler
  
  #component suggest output
  output$compsuggest <- renderPlot({
    selectNcomp(getModel()[[1]], method="onesigma", plot=TRUE)
  }) #end output
  
  
  #PCA visuals (OPTIONAL)
  # #PCA summary
  # output$prcompSummary <- renderPrint({
  #   PCs <- prcomp(getModel()[[2]])
  #   summary(PCs)
  # })
  
  #PCA biplot output
  # output$biplot <- renderPlot({
  #   PCs <- prcomp(getModel()[[2]])
  #   fviz_pca_biplot(PCs, axes = c(1,2), label = "var", alpha = .5)
  #   fviz_pca_var(PCs, axes = c(1,2), col.var = "contrib", label = "var", alpha = .5)
  #   biplot(getModel()[[1]], which = "loadings")
  # }) #end output
  
  # #PCA loadings
  # output$plsSummary <- renderPrint({
  #   pca <- prcomp(getModel()[[2]])
  #   pca$rotation
  # }) 
  
  # #plot download handler
  # output$downloadPCAPlot <- downloadHandler(
  #   filename = function() {
  #     paste0("COEFPLOT", ".png")
  #   },
  #   content = function(file) {
  #     png(file, width = 1200)
  #     plot(getModel()[[1]], plottype = "coef", ncomp=1:3, legendpos="bottomleft", labels = "names")
  #     dev.off()
  #   }
  # ) #end PCA plot download handler
  # 
  # #data download handler (selected data)
  # output$downloadPCAData <- downloadHandler(
  #   filename = function() {
  #     paste0("PLS_SUMMARY", ".csv")
  #   },
  #   content = function(file) {
  #     pca <- prcomp(getModel()[[2]])
  #     content <- data.frame(pca$x)
  #     write.csv(content, file, row.names = TRUE)
  #   }
  # ) #end data download handler
  
  
  
  #_____________________________________________________________________________________
  #Upload tab
  
  #spectra output
  output$spectra <- renderPlot({
    
    df <- getData()
    
    #divide data set
    #ImageUse <- input$images
    training <- filter(df, ImageUse=="Training")
    test <- filter(df, ImageUse=="Test")
    
    #If manually selecting model data
    if(input$rowselect == 1) {
      rows <- input$contents_rows_selected
      training <- training[rows,]
      test <- rbind(test, test[-rows,])
    }
    
    #toggle between test and train data
    if(input$datatype == "test") {
      #x <- (test[,-c(1:7)])
      x <- dplyr::select(test, matches(".r$|.g$|.b$|[[:digit:]]\\."))
    } else if(input$datatype == "training") {
      x <- dplyr::select(training, matches(".r$|.g$|.b$|[[:digit:]]\\."))
    } else {
      x <- dplyr::select(df, matches(".r$|.g$|.b$|[[:digit:]]\\."))
    }
    
    x <- dplyr::select(df, matches(".r$|.g$|.b$|[[:digit:]]\\."))
    
    mdaplot(x, type = 'l')
    
  })
  
  
  #data table output
  output$contents <- renderDT({
    
    #req(input$file1)
    
    tryCatch(
      {
        df <- getData()
        
        #divide data set
        #ImageUse <- input$images
        training <- filter(df, ImageUse=="Training")
        test <- filter(df, ImageUse=="Test")
        
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    
    #toggle between test and train data
    if(input$datatype == "test") {
      return(test[,1:8])
    } else if(input$datatype == "training") {
      return(training[,1:8])
    } else {
      return(df[,1:8])
    }
    
    # cols <- dim(df)[2]    
    # 
    # return(df[,1:8])
    
    options = list(scrollX = TRUE)
    
  })
  
  output$info = renderPrint(
    input$contents_rows_selected
  )
  
  #boxplot
  output$boxplot <- renderPlot({
    df <- getPredict()
    
    #add concatenated contains_api variable
    df$API_Contains <- ifelse(df$Contains==input$selectAPI, paste0(df$Contains," ",df$API,"%"), paste0(df$Contains))
    
    #delete oscorespls from Params
    df$Parameters <- str_replace(df$Parameters, "oscorespls, ","")
    
    #add model var
    df$Model <- paste0(df$Parameters,"\n",df$Lanes,"\nFactors=", df$Factors,
                       "\nRMSEP=", round(df$Test_RMSEP,2))
    
    #loop
    for (i in 1:nlevels(factor(df$file)) ) {
      #factor df$file
      df$file <- as.factor(df$file)
      
      #summarize
      predict_summary <- df %>% filter(df$file==levels(df$file)[i]) %>% 
        group_by(API_Contains, Lanes, Parameters) %>% 
        summarise(mean=mean(Predicted.API), sd=sd(Predicted.API))
      
      # #write csv
      # filename <- paste0(gsub("\\-.*","", levels(df$file)[i]),"_summary.csv")
      # write.csv(predict_summary , 
      #           paste0("c:/users/ssortijas/downloads/",filename))
      
      df2 <- df %>% filter(df$file == levels(df$file)[i])
      
      for (j in 1:nlevels(factor(df2$Lanes))) {
        #graphs
        # graphname <- paste0("c:/users/ssortijas/downloads/", 
        #                     gsub("\\-.*","", levels(df$file)[i]),"_",
        #                     gsub(",","",levels(factor(df$Lanes))[j]),".png")
        
        graph <- df2 %>%
          mutate(
            API_Contains = if(input$selectAPI == "amoxicillin") {
              fct_relevel(API_Contains, "amoxicillin 20%", 
                          "amoxicillin 50%",
                          "amoxicillin 80%", "amoxicillin 100%"
                          # ,
                          # "albendazole","ampicillin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin","doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "ciprofloxacin") {
              fct_relevel(API_Contains, "ciprofloxacin 20%", 
                          "ciprofloxacin 50%",
                          "ciprofloxacin 80%", "ciprofloxacin 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "azithromycin") {
              fct_relevel(API_Contains, "azithromycin 20%", 
                          "azithromycin 50%",
                          "azithromycin 80%", "azithromycin 100%",
                          "albendazole", "amoxicillin", "ampicillin",
                          "benzyl penicillin",
                          "calcium carbonate", "ceftriaxone",
                          "ciprofloxacin", "doxycycline",
                          "hydroxychloroquine", "isoniazid","lactose",
                          "tetracycline")
            } else if(input$selectAPI == "chloroquine") {
              fct_relevel(API_Contains, "chloroquine 20%", 
                          "chloroquine 50%",
                          "chloroquine 80%", "chloroquine 100%",
                          "albendazole", "amoxicillin", "ampicillin",
                          "azithromycin", "benzyl penicillin",
                          "calcium carbonate", "ceftriaxone",
                          "ciprofloxacin", "doxycycline",
                          "hydroxychloroquine", "isoniazid","lactose",
                          "tetracycline")
            } else if(input$selectAPI == "hydroxychloroquine") {
              fct_relevel(API_Contains, "hydroxychloroquine 20%", 
                          "hydroxychloroquine 50%",
                          "hydroxychloroquine 80%", "hydroxychloroquine 100%"
                          ,
                          "albendazole", "amoxicillin", "ampicillin",
                          "azithromycin", "benzyl penicillin",
                          "calcium carbonate", "ceftriaxone",
                          "ciprofloxacin", "doxycycline",
                          "hydroxychloroquine", "isoniazid","lactose",
                          "tetracycline")
            } else if(input$selectAPI == "isoniazid") {
              fct_relevel(API_Contains, "isoniazid 20%", 
                          "isoniazid 50%",
                          "isoniazid 80%", "isoniazid 100%",
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          "tetracycline"
              )
            } else if(input$selectAPI == "ampicillin") {
              fct_relevel(API_Contains, "ampicillin 20%", 
                          "ampicillin 50%",
                          "ampicillin 80%", "ampicillin 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "benzyl penicillin") {
              fct_relevel(API_Contains, "benzyl penicillin 20%", 
                          "benzyl penicillin 50%",
                          "benzyl penicillin 80%", "benzyl penicillin 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "ceftriaxone") {
              fct_relevel(API_Contains, "ceftriaxone 20%", 
                          "ceftriaxone 50%",
                          "ceftriaxone 80%", "ceftriaxone 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "doxycycline") {
              fct_relevel(API_Contains, "doxycycline 20%", 
                          "doxycycline 50%",
                          "doxycycline 80%", "doxycycline 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "tetracycline") {
              fct_relevel(API_Contains, "tetracycline 20%", 
                          "tetracycline 50%",
                          "tetracycline 80%", "tetracycline 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "albendazole") {
              fct_relevel(API_Contains, "albendazole 20%", 
                          "albendazole 50%",
                          "albendazole 80%", "albendazole 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "epinephrine") {
              fct_relevel(API_Contains, "epinephrine 20%", 
                          "epinephrine 50%",
                          "epinephrine 80%", "epinephrine 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "ethambutol") {
              fct_relevel(API_Contains, "ethambutol 20%", 
                          "ethambutol 50%",
                          "ethambutol 80%", "ethambutol 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "ferrous sulfate") {
              fct_relevel(API_Contains, "ferrous sulfate 20%", 
                          "ferrous sulfate 50%",
                          "ferrous sulfate 80%", "ferrous sulfate 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "ripe") {
              fct_relevel(API_Contains, "ripe 20%", 
                          "ripe 50%",
                          "ripe 80%", "ripe 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "rifampicin") {
              fct_relevel(API_Contains, "rifampicin 20%", 
                          "rifampicin 50%",
                          "rifampicin 80%", "rifampicin 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "sulfamethoxazole") {
              fct_relevel(API_Contains, "sulfamethoxazole 20%", 
                          "sulfamethoxazole 50%",
                          "sulfamethoxazole 80%", "sulfamethoxazole 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "pyrazinamide") {
              fct_relevel(API_Contains, "pyrazinamide 20%", 
                          "pyrazinamide 50%",
                          "pyrazinamide 80%", "pyrazinamide 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            } else if(input$selectAPI == "promethazine hydrochloride") {
              fct_relevel(API_Contains, "promethazine hydrochloride 20%", 
                          "promethazine hydrochloride 50%",
                          "promethazine hydrochloride 80%", "promethazine hydrochloride 100%"
                          # ,
                          # "albendazole", "amoxicillin", "ampicillin",
                          # "azithromycin", "benzyl penicillin",
                          # "calcium carbonate", "ceftriaxone",
                          # "ciprofloxacin", "doxycycline",
                          # "hydroxychloroquine", "isoniazid","lactose",
                          # "tetracycline"
              )
            }
            
            , 
            Model = fct_relevel(Model,levels(factor(Model))[7],levels(factor(Model))[1],
                                levels(factor(Model))[9],levels(factor(Model))[5],
                                levels(factor(Model))[3],levels(factor(Model))[8],
                                levels(factor(Model))[2],levels(factor(Model))[10],
                                levels(factor(Model))[6],levels(factor(Model))[4])) %>%
          filter(Lanes == levels(factor(df$Lanes))[j]) %>%
          ggplot(aes(x=API_Contains, y=Predicted.API, fill=API_Contains)) +
          geom_boxplot() + facet_grid(.~Model) +
          geom_hline(yintercept = c(0,20,50,80,100), linetype="dashed", color="red") +
          theme_bw() + ggtitle(paste0(gsub("\\-.*","", levels(df$file)[i]),", ", input$selectAPI," PLSR Prediction with Lanes ", levels(factor(df$Lanes))[j]), input$images) +
          theme(text=element_text(family="Gill Sans MT",size=14), 
                legend.title = element_blank(),
                axis.ticks.x = element_blank(),axis.text.x = element_blank()) +
          xlab("") + ylab("Predicted API")
        
        #png(file=graphname, width = 1234, height = 484)
        # print(graph)
        # dev.off()
        
      }
    }
    return(graph)
  }) #end output
  
  
  
  
}