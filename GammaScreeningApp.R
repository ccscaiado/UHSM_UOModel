#Works. 

library(shinycssloaders)
library(shiny)
library(shinythemes)
library(ggplot2)
library(plotly)
library(mvtnorm)

#Combined RenalFailure
#RenalFailure1 = read.csv("/Volumes/Jordan/UHSM_2016/RenalFailure.csv")
#RenalFailure2 = read.csv("/Volumes/Jordan/UHSM_2017/RenalFailure.csv") 
#RenalFailure3 = read.csv("/Volumes/Jordan/UHSM_2017_vol2/RenalFailure.csv")
#RenalFailure1 = subset(RenalFailure1, select = -c(PseudoId))
#RenalFailure1 = RenalFailure1[, names(RenalFailure2)]
#RenalFailure3 = RenalFailure3[, names(RenalFailure2)]
#RenalFailure <-rbind(RenalFailure1, RenalFailure2, RenalFailure3)

RenalFailure = read.csv("/Volumes/Jordan/UHSM_2017/RenalFailure.csv")

#Combined PatientIndex
#PatientIndex1 = read.csv("/Volumes/Jordan/UHSM_2016/PatientIndex.csv")
#PatientIndex1 = subset(PatientIndex1, select = -c(FirstStage2c, AKI2Reason, FirstStage3c, AKI3Reason, CPBTime, ProcDetailsFinal, OtherActualCardiaProcsFinal, AddtionOtherNonCardiacSugeryFinal))
#PatientIndex2 = read.csv("/Volumes/Jordan/UHSM_2017/PatientIndex.csv")
#PatientIndex2 = subset(PatientIndex2, select = -c(FirstAKI1UO, FirstAKI1AbsoluteCr, FirstAKI1RelativeCr, CPB, ProcDetails, OtherActualCardProcs, AdditionalOtherNonCardiacSurgery))
#PatientIndex3 = read.csv("/Volumes/Jordan/UHSM_2017_vol2/PatientIndex.csv")
#PatientIndex3 = subset(PatientIndex3, select = -c(FirstAKI1UO, FirstAKI1AbsoluteCr, FirstAKI1RelativeCr, CPB, ProcDetails, OtherActualCardProcs, AdditionalOtherNonCardiacSurgery))
#PatientIndex = rbind(PatientIndex1, PatientIndex2, PatientIndex3)

PatientIndex = read.csv("/Volumes/Jordan/UHSM_2017/PatientIndex.csv")

#Combined FlowSheet
#FlowSheet1 = read.csv("/Volumes/Jordan/UHSM_2016/FlowSheet.csv", stringsAsFactors = FALSE)
#FlowSheet1 = subset(FlowSheet1, select = -c(PseudoId, I.E, VVECMO))
#FlowSheet2 = read.csv("/Volumes/Jordan/UHSM_2017/FlowSheet.csv", stringsAsFactors = FALSE)
#FlowSheet3 = read.csv("/Volumes/Jordan/UHSM_2017_vol2/FlowSheet.csv", stringsAsFactors = FALSE)
#FlowSheet = rbind(FlowSheet1, FlowSheet2, FlowSheet3)
FlowSheet = read.csv("/Volumes/Jordan/UHSM_2017/FlowSheet.csv", stringsAsFactors = FALSE)

#Combined Fluids
#Fluids1 = read.csv("/Volumes/Jordan/UHSM_2016/FinalAKIFluids.csv")
#Fluids2 = read.csv("/Volumes/Jordan/UHSM_2017/Fluids.csv")
#Fluids3 = read.csv("/Volumes/Jordan/UHSM_2017_vol2/Fluids 2.csv")
#Fluids1 = subset(Fluids1, select = -c(PseudoId, Quarter, Drain1..mL...Hrly.Dra.1, Furosemide.60mg.1, Albumin.5....1))
#Fluids2 = Fluids1[, names(Fluids1)]
#Fluids3 = Fluids2[, names(Fluids1)]
#Fluids <- rbind(Fluids1, Fluids2, Fluids3)

#Fluids$Furosemide.10mg.ml.Rate[which(Fluids$Furosemide.10mg.ml.Rate == 0)] = NA
#Fluids$Furosemide.10mg.ml.Volume[which(Fluids$Furosemide.10mg.ml.Volume == 0)] = NA
#Fluids$Furosemide.10mg[which(Fluids$Furosemide.10mg == 0)] = NA
#Fluids$Furosemide.20mg[which(Fluids$Furosemide.20mg == 0)] = NA
#Fluids$Furosemide.40mg[which(Fluids$Furosemide.40mg == 0)] = NA
#Fluids$Furosemide.60mg[which(Fluids$Furosemide.60mg == 0)] = NA
#Fluids$Furosemide.80mg[which(Fluids$Furosemide.80mg == 0)] = NA

#Combined Biochemistry
#Biochemistry1 = read.csv("/Volumes/Jordan/UHSM_2016/Biochemistry.csv", stringsAsFactors = FALSE)
#Biochemistry1 = subset(Biochemistry1, select = -c(PseudoId, OpStart.x, Creat, Alb, Bili))
#Biochemistry2 = read.csv("/Volumes/Jordan/UHSM_2017/Biochemistry.csv", stringsAsFactors = FALSE)
#colnames(Biochemistry2)[2] = "NewPseudoId"
#Biochemistry2 = subset(Biochemistry2, select = -c(TIME, Creatinine, Albumin, Bilirubin, OpStart))
#Biochemistry3 = read.csv("/Volumes/Jordan/UHSM_2017_vol2/Biochemistry.csv", stringsAsFactors = FALSE)
#Biochemistry3 = subset(Biochemistry3, select = -c(TIME, Creatinine, Albumin, Bilirubin, OpStart))
#Biochemistry <- rbind(Biochemistry1, Biochemistry2, Biochemistry3)

Biochemistry = read.csv("/Volumes/Jordan/UHSM_2017/Biochemistry.csv", stringsAsFactors = FALSE)
colnames(Biochemistry)[2] = "NewPseudoId"

inputV = 0.1
inputFreeParameters = 2
inputm = "0.55,-0.1"
inputC = "0.01,0.001"
inputG = "1,1,0,1"
inputFF = "1,0"
inputdelta= 0.85
inputdel = 0.8
inputd = 0.1
inputn = 1
inputS = 0.1

ui <- navbarPage(
  theme = shinytheme("cerulean"),
  "DLM",
  navbarMenu("DLM",
             tabPanel("DLM Test",
                      sidebarLayout(
                        sidebarPanel(
                          fluidRow(
                            
                            column(5, 
                                   numericInput(inputId = "k", label = "Enter k", value = 6, min = 1)),
                            
                            column(5, 
                                   numericInput(inputId = "Observation", label = "Enter Obs", value = 0, min = 0))
                            
                          ),
                          
                          fluidRow(
                            
                            column(7, 
                                   selectInput(inputId = 'Id', label = 'Select Patient', unique(RenalFailure$NewPseudoId)))
                          ),
                          fluidRow(
                            textOutput("FirstAKI1UO"),
                            textOutput("FirstFilter"),
                            textOutput("Gender"),
                            textOutput("Age"),
                            textOutput("Weight"),
                            textOutput("Height"),
                            textOutput("EuroSCORE"),
                            textOutput("ProcDetails"),
                            textOutput("Urgency")
                          )
                          
                        ),
                        mainPanel(
                          verbatimTextOutput("JointProbability"),
                          withSpinner(plotlyOutput("scat"), type = 5), 
                          withSpinner(tableOutput("df"), type = 5),
                          withSpinner(plotlyOutput("scat3"), type = 5),
                          withSpinner(plotlyOutput("scat4"), type = 5),
                          withSpinner(plotlyOutput("scat5"), type = 5),
                          withSpinner(plotlyOutput("scat6"), type = 5)
                        )
                      )
             )
  )
)

# Define server logic ----
server <- function(input, output) {
  
  output$df <-renderTable({
    
    SubsetRenalFailureId = subset(RenalFailure, NewPseudoId == input$Id) 
    SubsetRenalFailureId2 = SubsetRenalFailureId[-1, ]
    PatientIndexId = subset(PatientIndex, NewPseudoId == input$Id)
    SubsetRenalFailureId2$Urine60 = SubsetRenalFailureId2$Urine60/PatientIndexId$Weight
    SubsetRenalFailureId2$Urine60[which(SubsetRenalFailureId2$Urine60==0)] = 0.1
    SubsetRenalFailureId2$Urine60 = log(SubsetRenalFailureId2$Urine60)
    y = SubsetRenalFailureId2$Urine60
    
    Cdiag <- as.numeric(unlist(strsplit(inputC,",")))
    C <- array(diag(Cdiag, nrow = inputFreeParameters), dim = c(inputFreeParameters,inputFreeParameters,1,length(y) + 2))
    mvec <- as.numeric(unlist(strsplit(inputm,",")))
    m <- array(mvec, dim = c(inputFreeParameters,1,1,length(y) + 2))
    d <- array(inputd, dim = c(1,1,1,length(y) + 2))
    n <- array(inputn, dim = c(1,1,1,length(y) + 2))
    S <- array(inputS, dim = c(1,1,1,length(y) + 2))
    a <- array(mvec, dim = c(inputFreeParameters,1,input$k +1 ,length(y) + 2))
    R <- array(diag(Cdiag, nrow = inputFreeParameters), dim = c(inputFreeParameters,inputFreeParameters,input$k +1 ,length(y) + 2))
    f <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    Q <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    P <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    LowerInterval <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    UpperInterval <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    A <- array(rep(NA,inputFreeParameters), dim = c(inputFreeParameters,1,1,length(y) + 1))
    e <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    
    Gvec <- as.numeric(unlist(strsplit(inputG,",")))
    G = matrix(Gvec, nrow = inputFreeParameters, byrow =TRUE)
    FF = as.numeric(unlist(strsplit(inputFF,",")))
    
    V=inputV
    W <- array(diag(Cdiag, nrow = inputFreeParameters), dim = c(inputFreeParameters,inputFreeParameters,1,length(y) + 1))
    
    for(i in 1:(length(y) + 1)){
      W[,,,i] <- (1/inputdelta - 1)*G%*%C[,,1,i]%*%t(G) 
      for(j in 1:input$k){
        a[,,j+1,i] <- G%*%a[,,j,i]
        R[,,j+1,i] <- G%*%R[,,j,i]%*%t(G) + W[,,,i]
        f[,,j,i] <- FF%*%a[,,j+1,i]
        Q[,,j,i] <- FF%*%R[,,j+1,i]%*%FF + S[,,,i]
        P[,,j,i] <- pnorm(log(0.3),f[,,j,i], sqrt(Q[,,j,i]))
        LowerInterval[,,j,i] = qnorm(c(0.025, 0.975), mean = f[,,j,i], sd = sqrt(Q[,,j,i]))[1] 
        UpperInterval[,,j,i] = qnorm(c(0.025, 0.975), mean = f[,,j,i], sd = sqrt(Q[,,j,i]))[2]
      }
      A[,,1,i] <- R[,,2,i]%*%FF/Q[,,1,i]
      e[,,,i] <- y[i:(i+input$k -1)] - f[,,,i]
      n[,,,i+1] <- inputdel*n[,,,i] + 1
      if(is.na(e[,,1,i])){
        n[,,,i+1] <- n[,,,i]
        d[,,,i+1] <- d[,,,i]
        m[,,1,i+1] <- m[,,1,i]
      } else {
        n[,,,i+1] <- inputdel*n[,,,i] + 1
        d[,,,i+1] <- inputdel*d[,,,i] + S[,,,i]*e[,,1,i]^2/Q[,,1,i]
        m[,,1,i+1] <- a[,,2,i] + A[,,1,i]*e[,,1,i]
      }
      S[,,,i+1] = d[,,,i+1]/n[,,,i+1]
      C[,,1,i+1] <- (S[,,,i+1]/S[,,,i])*(R[,,2,i] - A[,,1,i]%*%t(A[,,1,i])*Q[,,1,i])
      a[,,1,i+1] <- m[,,1,i+1]
      R[,,1,i+1] <- C[,,1,i+1]
    }
    
    df = data.frame((input$Observation+1):(input$Observation+input$k),exp(f[,,,input$Observation+1]), exp(Q[,,,input$Observation+1]), exp(LowerInterval[,,,input$Observation+1]), exp(y[(input$Observation+1):(input$Observation+input$k)]),exp(y[(input$Observation+1):(input$Observation+input$k)]) - exp(f[,,,input$Observation+1]), format(P[,,,input$Observation+1], digits = 5))#, UnderPredict[,,,input$Observation+1]) 
    colnames(df) = c("X","f", "Q", "LowerInterval", "y", "e", "P")
    return(df)
    
  })
  
  output$scat <- renderPlotly({
    
    SubsetRenalFailureId = subset(RenalFailure, NewPseudoId == input$Id) 
    SubsetRenalFailureId2 = SubsetRenalFailureId[-1, ]
    PatientIndexId = subset(PatientIndex, NewPseudoId == input$Id)
    SubsetRenalFailureId2$Urine60 = SubsetRenalFailureId2$Urine60/PatientIndexId$Weight
    SubsetRenalFailureId2$Urine60[which(SubsetRenalFailureId2$Urine60==0)] = 0.1
    SubsetRenalFailureId2$Urine60 = log(SubsetRenalFailureId2$Urine60)
    y = SubsetRenalFailureId2$Urine60
    
    Cdiag <- as.numeric(unlist(strsplit(inputC,",")))
    C <- array(diag(Cdiag, nrow = inputFreeParameters), dim = c(inputFreeParameters,inputFreeParameters,1,length(y) + 2))
    mvec <- as.numeric(unlist(strsplit(inputm,",")))
    m <- array(mvec, dim = c(inputFreeParameters,1,1,length(y) + 2))
    d <- array(inputd, dim = c(1,1,1,length(y) + 2))
    n <- array(inputn, dim = c(1,1,1,length(y) + 2))
    S <- array(inputS, dim = c(1,1,1,length(y) + 2))
    a <- array(mvec, dim = c(inputFreeParameters,1,input$k +1 ,length(y) + 2))
    R <- array(diag(Cdiag, nrow = inputFreeParameters), dim = c(inputFreeParameters,inputFreeParameters,input$k +1 ,length(y) + 2))
    f <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    Q <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    P <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    LowerInterval <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    UpperInterval <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    A <- array(rep(NA,inputFreeParameters), dim = c(inputFreeParameters,1,1,length(y) + 1))
    e <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    
    Gvec <- as.numeric(unlist(strsplit(inputG,",")))
    G = matrix(Gvec, nrow = inputFreeParameters, byrow =TRUE)
    FF = as.numeric(unlist(strsplit(inputFF,",")))
    
    V=inputV
    W <- array(diag(Cdiag, nrow = inputFreeParameters), dim = c(inputFreeParameters,inputFreeParameters,1,length(y) + 1))
    
    for(i in 1:(length(y) + 1)){
      W[,,,i] <- (1/inputdelta - 1)*G%*%C[,,1,i]%*%t(G) 
      for(j in 1:input$k){
        a[,,j+1,i] <- G%*%a[,,j,i]
        R[,,j+1,i] <- G%*%R[,,j,i]%*%t(G) + W[,,,i]
        f[,,j,i] <- FF%*%a[,,j+1,i]
        Q[,,j,i] <- FF%*%R[,,j+1,i]%*%FF + S[,,,i]
        P[,,j,i] <- pnorm(log(0.3),f[,,j,i], sqrt(Q[,,j,i]))
        LowerInterval[,,j,i] = qnorm(c(0.025, 0.975), mean = f[,,j,i], sd = sqrt(Q[,,j,i]))[1] 
        UpperInterval[,,j,i] = qnorm(c(0.025, 0.975), mean = f[,,j,i], sd = sqrt(Q[,,j,i]))[2]
      }
      A[,,1,i] <- R[,,2,i]%*%FF/Q[,,1,i]
      e[,,,i] <- y[i:(i+input$k -1)] - f[,,,i]
      n[,,,i+1] <- inputdel*n[,,,i] + 1
      if(is.na(e[,,1,i])){
        n[,,,i+1] <- n[,,,i]
        d[,,,i+1] <- d[,,,i]
        m[,,1,i+1] <- m[,,1,i]
      } else {
        n[,,,i+1] <- inputdel*n[,,,i] + 1
        d[,,,i+1] <- inputdel*d[,,,i] + S[,,,i]*e[,,1,i]^2/Q[,,1,i]
        m[,,1,i+1] <- a[,,2,i] + A[,,1,i]*e[,,1,i]
      }
      S[,,,i+1] = d[,,,i+1]/n[,,,i+1]
      C[,,1,i+1] <- (S[,,,i+1]/S[,,,i])*(R[,,2,i] - A[,,1,i]%*%t(A[,,1,i])*Q[,,1,i])
      a[,,1,i+1] <- m[,,1,i+1]
      R[,,1,i+1] <- C[,,1,i+1]
    }
    
    df = data.frame((input$Observation+1):(input$Observation+input$k),exp(f[,,,input$Observation+1]), exp(Q[,,,input$Observation+1]), exp(LowerInterval[,,,input$Observation+1]), exp(UpperInterval[,,,input$Observation+1]), exp(y[(input$Observation+1):(input$Observation+input$k)]), exp(y[(input$Observation+1):(input$Observation+input$k)]) - exp(f[,,,input$Observation+1]))
    colnames(df) = c("X","f", "Q", "LowerInterval", "UpperInterval", "y", "e")
    dfActualUO = data.frame(1:(input$k+input$Observation), exp(y[1:(input$k+input$Observation)]))
    colnames(dfActualUO) = c("ActualUO", "y")
    
    ggplot()+
      geom_point(data=df, aes(x=X, y =f, colour = "f"))+
      geom_point(data=dfActualUO, aes(x=ActualUO, y = y, colour = "UO"))+
      geom_point(data=df, aes(x=X, y =LowerInterval, colour = "IntervalLower"), shape=95)+
      geom_segment(data=df, aes(x=X, y = LowerInterval, xend=X, yend = f))+
      geom_hline(aes(yintercept = 0.5, colour = "0.5line"))+
      geom_hline(aes(yintercept = 0.3, colour = "0.3line"))+
      scale_colour_manual("", 
                          breaks = c("f", "UO", "IntervalLower", "0.5line", "0.3line"),
                          values = c("purple", "red", "red", "black", "black"))+
      labs(title="UO")+
      labs(x = "Time(Hours)")+
      labs(y = "UO(ml/kg)")+
      theme_bw()
    
  })
  
  output$JointProbability <-renderPrint({
    
    SubsetRenalFailureId = subset(RenalFailure, NewPseudoId == input$Id) 
    SubsetRenalFailureId2 = SubsetRenalFailureId[-1, ]
    PatientIndexId = subset(PatientIndex, NewPseudoId == input$Id)
    SubsetRenalFailureId2$Urine60 = SubsetRenalFailureId2$Urine60/PatientIndexId$Weight
    SubsetRenalFailureId2$Urine60[which(SubsetRenalFailureId2$Urine60==0)] = 0.1
    SubsetRenalFailureId2$Urine60 = log(SubsetRenalFailureId2$Urine60)
    y = SubsetRenalFailureId2$Urine60
    
    Cdiag <- as.numeric(unlist(strsplit(inputC,",")))
    C <- array(diag(Cdiag, nrow = inputFreeParameters), dim = c(inputFreeParameters,inputFreeParameters,1,length(y) + 2))
    mvec <- as.numeric(unlist(strsplit(inputm,",")))
    m <- array(mvec, dim = c(inputFreeParameters,1,1,length(y) + 2))
    d <- array(inputd, dim = c(1,1,1,length(y) + 2))
    n <- array(inputn, dim = c(1,1,1,length(y) + 2))
    S <- array(inputS, dim = c(1,1,1,length(y) + 2))
    a <- array(mvec, dim = c(inputFreeParameters,1,input$k +1 ,length(y) + 2))
    R <- array(diag(Cdiag, nrow = inputFreeParameters), dim = c(inputFreeParameters,inputFreeParameters,input$k +1 ,length(y) + 2))
    f <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    Q <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    P <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    LowerInterval <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    UpperInterval <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    A <- array(rep(NA,inputFreeParameters), dim = c(inputFreeParameters,1,1,length(y) + 1))
    e <- array(rep(NA,1), dim = c(1,1,input$k,length(y) + 1))
    Qij <- array(NA, dim=c(1, 1, input$k,length(y)+1))
    
    Gvec <- as.numeric(unlist(strsplit(inputG,",")))
    G = matrix(Gvec, nrow = inputFreeParameters, byrow =TRUE)
    FF = as.numeric(unlist(strsplit(inputFF,",")))
    
    V=inputV
    W <- array(diag(Cdiag, nrow = inputFreeParameters), dim = c(inputFreeParameters,inputFreeParameters,1,length(y) + 1))
    
    for(i in 1:(length(y) + 1)){
      W[,,,i] <- (1/inputdelta - 1)*G%*%C[,,1,i]%*%t(G) 
      for(j in 1:input$k){
        a[,,j+1,i] <- G%*%a[,,j,i]
        R[,,j+1,i] <- G%*%R[,,j,i]%*%t(G) + W[,,,i]
        f[,,j,i] <- FF%*%a[,,j+1,i]
        Q[,,j,i] <- FF%*%R[,,j+1,i]%*%FF + S[,,,i]
        P[,,j,i] <- pnorm(log(0.3),f[,,j,i], sqrt(Q[,,j,i]))
        LowerInterval[,,j,i] = qnorm(c(0.025, 0.975), mean = f[,,j,i], sd = sqrt(Q[,,j,i]))[1] 
        UpperInterval[,,j,i] = qnorm(c(0.025, 0.975), mean = f[,,j,i], sd = sqrt(Q[,,j,i]))[2]
      }
      A[,,1,i] <- R[,,2,i]%*%FF/Q[,,1,i]
      e[,,,i] <- y[i:(i+input$k -1)] - f[,,,i]
      n[,,,i+1] <- inputdel*n[,,,i] + 1
      if(is.na(e[,,1,i])){
        n[,,,i+1] <- n[,,,i]
        d[,,,i+1] <- d[,,,i]
        m[,,1,i+1] <- m[,,1,i]
      } else {
        n[,,,i+1] <- inputdel*n[,,,i] + 1
        d[,,,i+1] <- inputdel*d[,,,i] + S[,,,i]*e[,,1,i]^2/Q[,,1,i]
        m[,,1,i+1] <- a[,,2,i] + A[,,1,i]*e[,,1,i]
      }
      S[,,,i+1] = d[,,,i+1]/n[,,,i+1]
      C[,,1,i+1] <- (S[,,,i+1]/S[,,,i])*(R[,,2,i] - A[,,1,i]%*%t(A[,,1,i])*Q[,,1,i])
      a[,,1,i+1] <- m[,,1,i+1]
      R[,,1,i+1] <- C[,,1,i+1]
    }
  
    for(i in 1:(length(y)+1)){
      Qij[,,,i] <- apply(R[,,-1,i], 3, function(x) FF%*%G%*%x%*%FF)
    }
    
    JointProbability = rep(0,(length(y) + 1))
    ProlongedHighRisk = rep(0,(length(y) + 1))
    
    for(i in 1:(length(y) + 1)){
      Joint = matrix(c(rep(Qij[,,,i], each = input$k)), byrow = TRUE, nrow = input$k)
      diag(Joint) = Q[,,,i]
      Joint[lower.tri(Joint)] <- t(Joint)[lower.tri(Joint)]
      JointProbability[i] = pmvnorm(lower = -Inf*rep(1,input$k), upper = rep(log(0.3),input$k), mean = f[,,,i], sigma = Joint, algorithm = Miwa())[1]
    }
    
    for(i in 2:(length(y) + 1)){
      if(JointProbability[i] > 0.8 & ProlongedHighRisk[i-1] == 0){
        ProlongedHighRisk[i] = 1
      } else if(JointProbability[i] > 0.8 & ProlongedHighRisk[i-1] != 0)
        ProlongedHighRisk[i] = ProlongedHighRisk[i-1] + 1
    }
    
    lst=list()
    if(JointProbability[input$Observation+1] >0.8){
      lst[[1]] = paste("The probability that the next", input$k, "UOs are below 0.3 is", JointProbability[input$Observation + 1])
      lst[[2]] = paste("This patient has been at high risk for", ProlongedHighRisk[input$Observation + 1], "hours")
      if(!is.na(PatientIndexId$FirstFilter)){
        lst[[3]] = paste("This patient was filtered at hour",which(SubsetRenalFailureId2$Time == as.character(PatientIndexId$FirstFilter)))
      }
      return(lst)
    } else {
      lst[[1]] = paste("The probability that the next", input$k, "UOs are below 0.3 is", JointProbability[input$Observation + 1])
      if(!is.na(PatientIndexId$FirstFilter)){
        lst[[2]] = paste("This patient was filtered at hour",which(SubsetRenalFailureId2$Time == as.character(PatientIndexId$FirstFilter)))
      }
      return(lst)
    }
    
  })
  
  output$FirstAKI1UO <- renderText({ 
    FirstAKI1UO = factor(PatientIndex$FirstAKI1UO[which(PatientIndex$NewPseudoId == input$Id)])
    paste("This patients FirstAKI1UO is at", FirstAKI1UO)
  })
  
  output$FirstFilter <- renderText({ 
    FirstFilter = factor(PatientIndex$FirstFilter[which(PatientIndex$NewPseudoId == input$Id)]) 
    paste("This patients FirstFilter is at",
          FirstFilter)
  })
  
  output$Gender = renderText({
    i = which(PatientIndex$NewPseudoId == input$Id)
    paste("Gender: ", PatientIndex$Gender[i])
  })
  
  output$Age = renderText({
    i = which(PatientIndex$NewPseudoId == input$Id)
    paste("Age: ", PatientIndex$Age[i])
  })
  
  output$Weight = renderText({
    i = which(PatientIndex$NewPseudoId == input$Id)
    paste("Weight(kg): ", PatientIndex$Weight[i])
  })
  
  output$Height = renderText({
    i = which(PatientIndex$NewPseudoId == input$Id)
    paste("Height(cm): ", PatientIndex$Height[i])
  })
  
  output$EuroSCORE = renderText({
    i = which(PatientIndex$NewPseudoId == input$Id)
    paste("EuroSCORE: ", PatientIndex$EuroSCORE[i])
  })
  
  output$ProcDetails = renderText({
    i = which(PatientIndex$NewPseudoId == input$Id)
    paste("ProcDetails: ", PatientIndex$ProcDetails[i])
  })
  
  output$Urgency = renderText({
    i = which(PatientIndex$NewPseudoId == input$Id)
    paste("Urgency: ", PatientIndex$Urgency[i])
  })
  
  output$scat3 <- renderPlotly({
    
    SubsetFlowSheeti = subset(FlowSheet, NewPseudoId == input$Id)
    Difftime = c(1, as.numeric(difftime(SubsetFlowSheeti$Result.DT[2:length(SubsetFlowSheeti$Result.DT)], SubsetFlowSheeti$Result.DT[1:(length(SubsetFlowSheeti$Result.DT)-1)], units="hours")))
    Times = cumsum(Difftime)
    RealTimes = Times[which(Times<=(input$Observation+1))]
    RealTimes = RealTimes - 1
    SubsetFlowSheetiRealTime = SubsetFlowSheeti[1:length(which(Times<=(input$Observation+1))),]
    
    ggplot(data=SubsetFlowSheetiRealTime, aes(x = RealTimes))+
      geom_point(aes(y = as.numeric(as.character(CVP)), colour = "CVP"))+
      geom_hline(aes(yintercept = 5, colour = "5line"))+
      scale_colour_manual("", 
                          breaks = c("CVP", "5line"),
                          values = c("red", "black"))+
      labs(title="CVP")+
      labs(x = "Time(hours)")+
      labs(y = "CVP")+
      theme_bw()
    
  })
  
  output$scat4 <- renderPlotly({
    
    Biochemistryi = subset(Biochemistry, NewPseudoId == input$Id)
    Biochemistryi = Biochemistryi[order(Biochemistryi$PostOpUsandEsTime),]
    
    FirstUreaCreatinineTime = as.character(Biochemistryi$PostOpUsandEsTime[1]) #This will need to be used to find the corresponding time and hence position in the RenalFailure table. So that the first Urea and Creatine are plotted at time, say 6, if the first measurements are made at the sixth hour of UOs for that person
    SubsetRenalFailurei = subset(RenalFailure, NewPseudoId == input$Id)
    FirstUreaCreatininePosition = which(SubsetRenalFailurei$Time == FirstUreaCreatinineTime)
    
    if(length(FirstUreaCreatininePosition)>0){
      FirstUreaCreatininePosition = FirstUreaCreatininePosition
    } else {
      FirstUreaCreatininePosition = as.numeric(difftime(FirstUreaCreatinineTime, SubsetRenalFailurei$Time[1], units="hours")) + 1
    }
    
    if(length(Biochemistryi$PostOpUsandEsTime)>1){
      Difftime = c(FirstUreaCreatininePosition, as.numeric(difftime(Biochemistryi$PostOpUsandEsTime[2:length(Biochemistryi$PostOpUsandEsTime)], Biochemistryi$PostOpUsandEsTime[1:(length(Biochemistryi$PostOpUsandEsTime)-1)], units="hours")))
    } else {
      Difftime = FirstUreaCreatininePosition
    }
    
    Times = cumsum(Difftime)
    RealTimes = Times[which(Times<=(input$Observation+1))]
    RealTimes = RealTimes - 1
    BiochemistryiRealTime = Biochemistryi[1:length(which(Times<=(input$Observation+1))),]
    
    if(length(RealTimes)>0){
      ggplot(data=BiochemistryiRealTime, aes(x = RealTimes))+
        geom_point(aes(y=Urea, colour = "Urea"))+
        geom_point(aes(y=Creatinine, colour = "Creatinine"))+
        scale_colour_manual("", 
                            breaks = c("Urea", "Creatinine"),
                            values = c("red", "black"))+
        labs(title="Urea and Creatinine")+
        labs(x = "Time(hours)")+
        labs(y = "Urea and Creatinine")+
        theme_bw()
    } else {
      ggplot()+theme_bw()
    }
    
  })
  
  output$scat5 <- renderPlotly({
    
    SubsetFlowSheeti = subset(FlowSheet, NewPseudoId == input$Id)
    
    difftime = as.numeric(difftime(SubsetFlowSheeti$Result.DT[2:length(SubsetFlowSheeti$Result.DT)], SubsetFlowSheeti$Result.DT[1:(length(SubsetFlowSheeti$Result.DT)-1)], units="hours"))
    FirstTime = as.character(SubsetFlowSheeti$Result.DT[1])
    SubsetRenalFailurei = subset(RenalFailure, NewPseudoId == input$Id)
    FirstPosition = which(SubsetRenalFailurei$Time == FirstTime)
    
    if(length(FirstPosition)>0){
      FirstPosition = FirstPosition
    } else {
      FirstPosition = as.numeric(difftime(FirstTime, SubsetRenalFailurei$Time[1], units="hours")) + 1
    }
    
    if(length(SubsetFlowSheeti$Result.DT)>1){
      Difftime = c(FirstPosition, difftime)
    } else {
      Difftime = FirstPosition
    }
    
    Times = cumsum(Difftime)
    RealTimes = Times[which(Times<=(input$Observation+1))]
    RealTimes = RealTimes - 1
    SubsetFlowSheetiRealTime = SubsetFlowSheeti[1:length(which(Times<=(input$Observation+1))),]
    
    if(length(RealTimes)>0){
      ggplot(data=SubsetFlowSheetiRealTime, aes(x = RealTimes))+
        geom_point(aes(y=pH, colour = "pH"))+
        geom_point(aes(y=as.numeric(as.character(BE.B.)), colour = "BE.B."))+
        geom_point(aes(y=Lac, colour = "Lac"))+
        geom_point(aes(y=as.numeric(as.character(K)), colour = "K"))+
        scale_colour_manual("", 
                            breaks = c("pH", "BE.B.", "Lac", "K"),
                            values = c("red", "black", "blue", "green4"))+
        labs(title="pH, BE.B., Lac, K")+
        labs(x = "Time(hours)") +
        labs(y = "pH, BE.B., Lac, K")+
        theme_bw()
    } else {
      ggplot()+theme_bw()
    }
    
  })
  
  output$scat6 <- renderPlotly({
    
    SubsetFlowSheeti = subset(FlowSheet, NewPseudoId == input$Id)
    Difftime = c(1, as.numeric(difftime(SubsetFlowSheeti$Result.DT[2:length(SubsetFlowSheeti$Result.DT)], SubsetFlowSheeti$Result.DT[1:(length(SubsetFlowSheeti$Result.DT)-1)], units="hours")))
    Times = cumsum(Difftime)
    Replace = which(is.na(SubsetFlowSheeti$ReliableART.M))
    SubsetFlowSheeti$ReliableART.M[Replace] = as.numeric(as.character(SubsetFlowSheeti$NBP.M[Replace]))
    RealTimes = Times[which(Times<=(input$Observation+1))]
    RealTimes = RealTimes - 1
    SubsetFlowSheetiRealTime = SubsetFlowSheeti[1:length(which(Times<=(input$Observation+1))),]
    
    ggplot(data=SubsetFlowSheetiRealTime, aes(x = RealTimes))+
      geom_point(aes(y=ReliableART.M, colour = "ReliableART.M"))+
      scale_colour_manual("", 
                          breaks = c("ReliableART.M"),
                          values = c("black"))+
      labs(title="ReliableART.M")+
      labs(x = "Time(hours)")+
      labs(y = "ReliableART.M")+
      theme_bw()
    
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)