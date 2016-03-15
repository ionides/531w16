library("shiny")
shinyServer(function(input,output){
  
  data<-read.csv("Spencer_Sales.csv",header=TRUE,stringsAsFactors=FALSE)
  data$date<-as.Date(as.character(data$date),format="%d-%m-%Y")
  #data$quantity<-as.numeric(data$quantity)
  #data$cog_value<-as.numeric(data$cog_value)
  #data$net_value<-as.numeric(data$net_value)
  
  
  region_var<-reactive({
    as.character(unique(data$region))
  })
  
  output$region<-renderUI({
    selectInput("region_sel","Select Region",choices=region_var())
  })
    
  
  plant_var<-reactive({
    as.character(unique(data[which(data$region==input$region_sel),]$plant))
  })
  
  output$plant<-renderUI({
    selectInput("plant_sel","Select Plant",choices=plant_var())
  })
  
  material_var<-reactive({
    
    library("stringr")
    z<-unique(data[which(data$region==input$region_sel & data$plant==input$plant_sel),]$material)
    for(i in 1:length(z)){
      z[i]<-paste(unlist(str_split(str_trim(gsub("[[:punct:]]","",z[i]))," ")),collapse="_")
    }
    unique(sort(z))
  })
  
  output$material<-renderUI({
    selectInput("material_sel","Select SKU",choices=material_var())
  })
  
  Test<-reactive({
    Agg_data<-data[which(data$region==input$region_sel & data$plant==input$plant_sel),]
    Agg<-aggregate(Agg_data$quantity,by=list(Agg_data$date,Agg_data$material),sum)
    names(Agg)<-c("Date","Material","Quantity")
    library("reshape2")
    library("stringr")
    library("forecast")
    Dcast_Sales<-dcast(Agg,Date~Material,fill=0)
    Dcast_Sales_order<-Dcast_Sales[order(Dcast_Sales$Date),]
    index<-c(1:round(0.95*(nrow(Dcast_Sales_order))))
    Test<-Dcast_Sales_order[-index,]
    for(k in 2:ncol(Dcast_Sales_order)){
      Test[,k]<-ifelse(Test[,k]==0,mean(Test[,k]),Test[,k])
    }
    Test
  })
  
    Training<-reactive({
    
    Agg_data<-data[which(data$region==input$region_sel & data$plant==input$plant_sel),]
    Agg<-aggregate(Agg_data$quantity,by=list(Agg_data$date,Agg_data$material),sum)
    names(Agg)<-c("Date","Material","Quantity")
    library("reshape2")
    library("stringr")
    library("forecast")
    Dcast_Sales<-dcast(Agg,Date~Material,fill=0)
    Dcast_Sales_order<-Dcast_Sales[order(Dcast_Sales$Date),]
    index<-c(1:round(0.95*(nrow(Dcast_Sales_order))))
    Training<-Dcast_Sales_order[index,]
    Test<-Dcast_Sales_order[-index,]
    colnames<-names(Dcast_Sales_order)
    HW<-as.character(1:ncol(Dcast_Sales_order))
    SD<-as.character(1:ncol(Dcast_Sales_order))
    number_forecasts<-nrow(Test)
    ci<-90
    for(k in 2:ncol(Dcast_Sales_order)){
      Training[,k]<-ifelse(Training[,k]==0,mean(Training[,k]),Training[,k])
    }
    for(i in 1:ncol(Dcast_Sales_order)){
      HW[i]<-paste(unlist(str_split(str_trim(gsub("[[:punct:]]","",colnames[i]))," ")),collapse="_")
    }
    
    for(j in 2:ncol(Dcast_Sales_order)){
      
      assign(HW[j],hw((ts(Training[,j],frequency=7)),h=number_forecasts,seasonal="additive",level=ci))
    }
    
    for(j in 2:ncol(Training)){
      
      SD[j]<-sd(eval(parse(text=HW[j]))$residuals)
    }
    
    for(k in 2:ncol(Training)){
      for(l in 1:nrow(Training)){
        Training[l,k]<-ifelse(abs(eval(parse(text=HW[k]))$residuals[l])>(2*as.numeric(SD[k])),eval(parse(text=HW[k]))$fitted[l]+(2*as.numeric(SD[k])*sign(eval(parse(text=HW[k]))$residuals[l])),Training[l,k])
      }
    }
    as.data.frame(Training)
    })
    
    
  
  model<-reactive({
      HW<-as.character(1:ncol(Training()))
      colnames<-names(Training())
      for(i in 1:ncol(Training())){
        HW[i]<-paste(unlist(str_split(str_trim(gsub("[[:punct:]]","",colnames[i]))," ")),collapse="_")
      }
    for(m in 2:ncol(Training())){
      
      assign(HW[m],hw((ts(Training()[,m],frequency=7)),h=5,seasonal="additive",level=90))
    }
    
    
    Outlist<-list(HW[2],(eval(parse(text=HW[2]))))
    #Outlist_name<-as.list(1:ncol(Dcast_Sales_order))
    for(e in 2:ncol(Training())){
      Outlist[[e]]<-list(HW[e],(eval(parse(text=HW[e]))))
      names(Outlist)[e]<-HW[e]
      
    }
    
    Outlist
    
    #HW
    #Outlist
    #model_out<-as.data.frame(cbind(Outlist,Outlist_name))
    #names(model_out)<-c("material","model")
    #model_out
    #nrow(Test[which(Test$material==input$material_sel),])

})

  
  output$accuracy<-renderTable({
    Accuracy<-as.list(rep(NA,length=length(model())))
    
    for(n in 2:length(model())){
      Accuracy[[n]]<-accuracy(model()[[n]][[2]]$mean,Test()[,n],test=which(Test()[,n]!=0))
    }
    
    
    Accuracy_cb<-as.list(rep(NA,length=length(model())))
    
    for(n in 2:length(model())){
      Accuracy_cb[[n]]<-cbind.data.frame(names(Test())[n],as.data.frame(Accuracy[[n]],stringsAsFactors=FALSE,row.names=NULL))
    }
    
    Accuracy_df<-do.call(rbind.data.frame,Accuracy_cb)
    
    head(Accuracy_df[order(Accuracy_df$MAPE),],20)
  })

  
  output$forecasted<-renderPrint({
    #names(Test())
    forecast<-as.data.frame(model()[[input$material_sel]])[,c(2,3,4)]
    forecast$Date<-weekdays(seq(Test()[1,1],length.out=nrow(Test()),by=1))
    forecast
  })

  output$decomposition<-renderPlot({
    library(timeSeries)
    y <-ts(Training()[,paste(unlist(str_split(input$material_sel,"_")),collapse=" ")],freq=7)
    z<-plot(stl(y,s.window="periodic"))
    z
    
  })

  output$stsm<-renderTable({
    y <-ts(Training()[,paste(unlist(str_split(input$material_sel,"_")),collapse=" ")],freq=7)
    X<-StructTS(y, type = "BSM")
    D<- forecast(X,h=nrow(Test()))
    test_data<-Test()[,paste(unlist(str_split(input$material_sel,"_")),collapse=" ")]
    train<-Training()[,paste(unlist(str_split(input$material_sel,"_")),collapse=" ")]
    accuracy(D$fitted,train,test=which(train!=0))
  })

  output$stplot<-renderPlot({
    y <-ts(Training()[,paste(unlist(str_split(str_trim(input$material_sel),"_")),collapse=" ")],freq=7)
    X<-StructTS(y, type = "BSM")
    plot(forecast(X,h=5),
         ylab="Sales")
  })
  
  output$data<-renderPlot({
    plot(Training()[,paste(unlist(str_split(str_trim(input$material_sel),"_")),collapse=" ")],type="l",ylab="Sales",xlab="Time")
  })
  
})
  
