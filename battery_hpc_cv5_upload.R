
## import the NSW electric load data in Feburary, 2018

DataSet1<-readRDS("data2017.rds")

## import the NSW electric load data in March, 2019

DataSet2<-readRDS("data2018.rds")

## training set
CyclVariabs<-function(a){
  t<-c(1:a)
  SC<-sin(2*pi/48*t)
  CC<-cos(2*pi/48*t)
  X_Cycle<-list(SC=SC,CC=CC)
}

## the training data: Feburary, 2019
TrainSet<-DataSet1$TOTALDEMAND
TrainSize<-length(TrainSet)
## the test data: March 2019
TestSet<-DataSet2$TOTALDEMAND
TestSize<-length(TestSet)

A<-CyclVariabs(TrainSize+TestSize)
XTrainCycle<-cbind(A$SC[49:TrainSize],A$CC[49:TrainSize])
XTestCycle<-cbind(A$SC[seq(TrainSize+49,TrainSize+TestSize,1)],A$CC[seq(TrainSize+49,TrainSize+TestSize,1)])
##Prepare for trainning

InputTrain<-array()
for (i in c(1:48)){
  InputTrain<-cbind(InputTrain, TrainSet[i:(TrainSize-49+i)])
}
InputTrain<-InputTrain[,2:49]
OutputTrain<-TrainSet[49:TrainSize]
k2Train <- DataSet1$RRP[49:TrainSize]
##Prepare for test

InputTest<-array()
for (i in c(1:48)){
  InputTest<-cbind(InputTest, TestSet[i:(TestSize-49+i)])
}
InputTest<-InputTest[,2:49]
OutputTest<-TestSet[49:TestSize]
k2Test <- DataSet2$RRP[49:TestSize]

############Data preprocessing#######

XYTrainSet<-as.data.frame(cbind(InputTrain,OutputTrain))
XYTestSet<-as.data.frame(cbind(InputTest,OutputTest))

c.TrainSet <- as.data.frame(cbind(XTrainCycle,XYTrainSet))
total.tr <- nrow(c.TrainSet)
c.TestSet <- as.data.frame(cbind(XTestCycle,XYTestSet))

###############################################
## Asymmetric LAV regression model training####
###############################################

library(L1pack)

Scale_LAV<-lad(formula = OutputTrain ~ ., data = c.TrainSet, method = "BR")


AsySVR_2<-function(beta, TrainX, TrainY, k1, k2, Const, cap, stp)
  
{
  eps1 <- cap*stp
  eps2 <- cap-eps1
  beta<-array(beta)
  TrainX_extra<-rep(1, times=dim(TrainX)[1])
  TrainX<-cbind(TrainX_extra,TrainX)
  Xmatrix<-as.matrix(TrainX)
  Err<-Xmatrix %*% beta-TrainY
  # under-ordered, error cost k2, eps1
  InsenErr1<- ifelse(Err< -eps1,-k1*(Err+eps1),0)
  # over-ordered, error cost k1, eps2,  
  InsenErr2<- ifelse(Err> eps2, k2*(Err-eps2), 0)
  
  SampleLLC<-InsenErr1+InsenErr2
  #rentC <- pcost*cap*nrow(TrainX)
  L2<-sum((beta[2:dim(TrainX)[2]])^2)
  TrainLLC<- L2+Const*(mean(SampleLLC)) # was sum here
}

AsySVR_2_ecost <- function(beta, TrainX, TrainY, k1, k2, Const, cap, stp){
  eps1 <- cap*stp
  eps2 <- cap-eps1
  beta<-array(beta)
  TrainX_extra<-rep(1, times=dim(TrainX)[1])
  TrainX<-cbind(TrainX_extra,TrainX)
  Xmatrix<-as.matrix(TrainX)
  Err<-Xmatrix %*% beta-TrainY
  # under-ordered, error cost k1, eps1
  InsenErr1<- ifelse(Err< -eps1,-k1*(Err+eps1),0)
  # over-ordered, error cost k2, eps2,  
  InsenErr2<- ifelse(Err> eps2, k2*(Err-eps2), 0)
  
  SampleLLC<-InsenErr1+InsenErr2
  #rentC <- pcost*cap
  TrainLLC<- (mean(SampleLLC))*48 #MDC daily cost
  
}
#For each set of k1, k2, Const, cap, stp, pcost -> beta

cv.id <- sample(1:5, nrow(c.TrainSet), replace=T)
cv.TrainSet<-split(c.TrainSet, cv.id)

cv.group <- c(1,2,3,4,5)


runid = 10#as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

KK2 <- c(70,80,90)#90

KK1 <- c(400,600,800,0)
Constset <-c(1,10) #c(0.001, 0.01, 0.1,1,10,100)

lnumk1 <- length(KK1)*length(Constset)

kk2 <- ceiling(runid/lnumk1)
kk1 <- ceiling(runid/length(Constset)) - (kk2-1)*length(KK1)

ccid <- runid - ceiling(runid/length(Constset)-1)*length(Constset)#runid - ((kk1-1)*4+kk2-1)*length(Constset)

Const <-Constset[ccid]


k2 <- KK2[kk2]/2
if(kk1==4){
  filename <- paste0("k1_",k2*2,"_k2_realtime_",Const,"_cv.rds") 
}else if(kk1<4){
  k1tr <- k1te <- c(400,600,800)[kk1]/2 #,1200
  filename <- paste0("k1_",k2*2,"_k2_",k1tr*2,"_",Const,"_cv.rds") 
}

Cset <- c(1,5,10,20,40,60,80,100,150,200,250,300,400,500,600,700,800,900,1000)

res <- data.frame()
for(cap in Cset){
  cat("\n")
  for(stp in seq(0,1,0.1)){
    cat(".")
    validate_Ecost  <- 0
    for(i in c(1:5)){
      c.TrainSet <- do.call(rbind,cv.TrainSet[cv.group[-i]])
      c.TestSet <- cv.TrainSet[[i]]
      
      if(kk1==4){
        idte <- which(cv.id==i)
        k1tr <-k1Train[-idte]/2 
        k1te <-k1Train[idte]/2
      }else if(kk2<4){
        k1tr <- k1te <- c(400,600,800)[kk1]/2  #1200
      }
      
      AsySVR_2_Coef<-optim(as.array(Scale_LAV$coefficients), function(beta) 
        AsySVR_2(beta, c.TrainSet[,1:50], c.TrainSet[,51], k1=k1tr,k2=k2, Const=Const, cap=cap, stp=stp), 
        method="L-BFGS-B")
      
      
      AsySVR_2_test_ecost <-   AsySVR_2_ecost(AsySVR_2_Coef$par, c.TestSet[,1:50], c.TestSet[,51], k1=k1te,k2=k2, Const=Const, cap=cap, stp=stp)
      #alpha_test <-mean(as.matrix(cbind(rep(1, dim(c.TestSet)[1]),c.TestSet[,1:50]))%*% as.array(AsySVR_2_Coef$par) - OutputTest>0)
      
      validate_Ecost <- validate_Ecost + AsySVR_2_test_ecost*nrow(c.TestSet)
    }
    
    
    res_tem <- as.data.frame(t(c(Const,cap,stp,validate_Ecost/total.tr)))
    
    res <- rbind(res,res_tem)
    
  }
}


saveRDS(res,filename)

