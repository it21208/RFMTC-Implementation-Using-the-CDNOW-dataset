library("Rcpp",character.only=TRUE)
library("optimization",character.only=TRUE)
library("optimx",character.only=TRUE)

# Application and Run Type constants
RUN_TYPE <- "TRAINING"  # TRAINING , TEST
RUN_TITLE <- paste("CDNOW RFMTC", RUN_TYPE, "DATASET")
# Data File constants
TransSampleDSFile <- "CDNOW_sample_trans.txt"
ChurnMasterFile <- "CDNOW_master_churn.txt"
CustSampleDSFile <- "CDNOW_sample_cust.txt"
TrainingDSFile <- "CDNOW_training.txt"
TestDSFile <- "CDNOW_test.txt"
# Optimization constants
TRAINING_OPTIMALS <- "TRAINING_OPTIMALS.csv"
# Data Frame Column Name constants
idColName <- "ID"
recencyColName = "Recency (months)"
frequencyColName = "Frequency (times)"
monetaryColName = "Monetary ($)"
timesColName = "Time (months)"
churnColName = "Churn (0/1)"
#
##############  2. FUNCTIONS #######################
#
### 2.1 CALCULATE (T-R)/(F-1) column ###
# ... with calculation =IF(F<>1,(T-R)/(F-1),0)
func_T_R_F_1 = function(r, f, t) {
  if (f != 1) 
    return((t - r) / (f - 1))
  else 
    return(0)
}
#
### 2.2 CALCULATE E[X[L]|L=1] ###
# ... formula 22 of the RFMTC paper
make.func_E_X_L_1 <- function(Q, g, b) {
  function(r, f, t) {
    return((g * f / (t - r + b)) * (1 - Q) ^ (r + 1))
  }
}
#
### 2.3 CALCULATE MOVING AVERAGE  ###
#
calcMovingAverage <- function(v, m) {
  len = length(v)
  ret_v <- rep(0, times = len)
  for (i in 1:len) {
    # find the position of the start and last element to summarize
    if (i <= m) {
      startPos <- 1
      endPos <- i + m
    } else if (i >= len - m) {
      endPos <- len
      startPos <- i - m
    } else {
      startPos <- i - m
      endPos <- i + m
    }
    # accumulate the elements in range [startPos..endPos]
    sum <- 0
    for (j in startPos:endPos)
      sum <- sum + v[j]
    # find average and store it in new vector
    ret_v[i] <- sum / (endPos - startPos + 1)
  }
  return(ret_v)
}
#
### 2.4  PARTITION CUSTOMER DATA SET TO TRAINING & TEST DATA SETS ###
#
partitionCustomerDSToTrainingTestDS <-
  function(CustSampleDSFile,
           TrainingDSFile,
           TestDSFile,
           RunType) {
    # RUN_TYPE is {"TRAINING" | "TEST"}
    # print function display header
    printf("partitionCustomerDSToTrainingTestDS()")
    # read custMasterFile
    custSDF <- read.csv(CustSampleDSFile, header = TRUE)
    # add autonumber column to create user-ids
    len <- nrow(custSDF)
    custSDF <- as.data.frame(cbind(1:len, custSDF))
    names(custSDF)[1] <- idColName
    # split training and test
    splitNo <- round(0.75*len)
    trainDF <- custSDF[custSDF[, idColName] <= splitNo,]
    testDF <- custSDF[custSDF[, idColName] > splitNo,]
    # store data frames to files
    write.csv(trainDF, TrainingDSFile, row.names = FALSE)
    write.csv(testDF, TestDSFile, row.names = FALSE)
    # set RunFile
    if (RunType == "TRAINING")
      RunFile <- TrainingDSFile
    else
      RunFile <- TestDSFile
    return(RunFile)
  }
# func_P_B_i <- function(churn_vector,m) {
#   new_churnVector <- 0
#   for (i in 1:length(churn_vector)){
#     if (i <= m) {
#       sum <- 0
#       for (j in i:(i+m)){sum = sum + churn_vector[j]}
#       k = i
#       while (k>1){
#         sum = sum + churn_vector[k-1]
#         k=k-1}
#       new_churnVector[i] = sum/(m+i)
#     }else if (i>=(length(churn_vector)-m)){
#       sum <- 0
#       for (j in i:(i-m)){sum = sum + churn_vector[j]}
#       k = i
#       while(k<length(churn_vector)){
#         sum = sum + churn_vector[k+1]
#         k=k+1}
#       new_churnVector[i]=sum/(length(churn_vector)-i+m+1)
#     }else{
#       sum <- 0
#       for (j in i:(i+m)) { sum = sum + churn_vector[j]}
#       for (k in (i-1):(i-m)) { sum = sum + churn_vector[k]}
#       new_churnVector[i]=sum/(2*m+1)}}
#   return(new_churnVector)}
#
### 2.5 CREATE OBJECTIVE FUNCTION TO MINIMIZE ###
# calculate sum of squares of residual errors
make.OBJECTIVE <- function(m, rfmtcReadyDF) {
  function(x) {
    Q <- x[1]
    g <- x[2]
    b <- x[3]
    # calculate E[X[L]|L=1]
    exl1 <- make.func_E_X_L_1(Q, g, b)
    rfmtcReadyDF$`E[X[L]|L=1]` <-
      mapply(function(r, f, t) exl1(r, f, t),
             rfmtcReadyDF$`Recency (months)`,
             rfmtcReadyDF$`Frequency (times)`,
             rfmtcReadyDF$`Time (months)`
      )
    # sort data on E[X[L]|L=1] in desc order
    rfmtcReadyDF <- rfmtcReadyDF[order(-rfmtcReadyDF$`E[X[L]|L=1]`), ]
    # calculate P[B]
    rfmtcReadyDF$`P[B]` <- calcMovingAverage(rfmtcReadyDF$`Churn (0/1)`, m)
    # calculater squares of residual errors
    rfmtcReadyDF$`(E[X|L=1]-P[B])^2` <-
      (rfmtcReadyDF$`E[X[L]|L=1]` - rfmtcReadyDF$`P[B]`) ^ 2
    return(sum(rfmtcReadyDF$`(E[X|L=1]-P[B])^2`))
  }
}
#
### 2.6 PRINTF ###
#
printf <- function(...) {
  invisible(print(sprintf(...)))
}
#
### 2.7 RFMTC DATASET PREPARATION ###
# RunFile schema : (CustID, Recency, Frequency, Monetary, Time, Churn)
prepareRFMTCdataset <- function(RunFile) {
  # print function display header
  printf("prepareRFMTCdataset()")
  # read RunFile with header
  rfmtcDF <- read.csv(RunFile, header = TRUE)
  len <- length(rfmtcDF[, 1])
  # add autonumber column to create user-ids and
  # create the rfmtcDF from the rfmtcDF and 4 blank score columns
  rfmtcDF <- as.data.frame(cbind(
    1:nrow(rfmtcDF),
    rfmtcDF[, 2:6],
    rep(0, times = len),
    rep(0, times = len),
    rep(0, times = len),
    rep(0, times = len)
  ))
  # define appropriate column labels
  names(rfmtcDF) <- c(
    idColName,
    recencyColName,
    frequencyColName,
    monetaryColName,
    timesColName,
    churnColName,
    "(T-R)/(F-1)", 
    "E[X[L]|L=1]", 
    "P[B]", 
    "(E[X|L=1]-P[B])^2"
  )
  return(rfmtcDF)
}
#
### 2.8 FIND RANGE CONSTRAINTS FOR PARAMETER (b) ###
#
calc_b_RangeConstraints <- function(rfmtcReadyDF) {
  # print function display header
  printf("calc_b_RangeConstraints()")
  rfmtcReadyDF$`(T-R)/(F-1)` <- 
    mapply(
      function(r, f, t)
        func_T_R_F_1(r, f, t),
      rfmtcReadyDF$`Recency (months)`,
      rfmtcReadyDF$`Frequency (times)`,
      rfmtcReadyDF$`Time (months)`
    )
  # calculate the average for non-zero values
  sum <- sum(rfmtcReadyDF$`(T-R)/(F-1)`)
  cnt <- sum(rfmtcReadyDF$`(T-R)/(F-1)` > 0)
  trf1_MEAN <- sum/cnt
  b_RangeDetailsList <- list(MEAN=trf1_MEAN, MIN=0.25 * trf1_MEAN, MAX=4 * trf1_MEAN)
  resList <- list(rfmtcReadyDF, b_RangeDetailsList)
  return(resList)
}
#
### 2.9 DISCOVER OPTIMAL SET OF INITIAL VALUES (m, Q, g, b) ###
###     TO USE FOR THE OPTIMIZATION ###
#
discoverOptimalInitialParams <- function(rfmtcReadyDF,
                                         m_range,
                                         Q_range,
                                         g_range,
                                         b_range,
                                         b_MIN,
                                         b_MAX) {
  # print function display header
  printf("discoverOptimalInitialParams()")
  attempts = length(m_range)*length(Q_range)*length(g_range)*length(b_range)
  resultData <- as.data.frame(cbind(c(1:attempts), rep(0, times=attempts),
                                    rep(0, times=attempts), rep(0, times=attempts),
                                    rep(0, times=attempts), rep(0, times=attempts), 
                                    rep(0, times=attempts), rep(0, times=attempts)))
  names(resultData) <- c("A/A", "m", "Q", "g", "b", "optim_sa", "optim_nm", "optimx")
  resultIdx <- 1
  for (m in m_range) {
    objfunc <- make.OBJECTIVE(m, rfmtcReadyDF)
    for (Q in Q_range) {
      for (g in g_range) {
        for (b in b_range) {
          resultData[resultIdx, "m"] <- m
          resultData[resultIdx, "Q"] <- Q
          resultData[resultIdx, "g"] <- g
          resultData[resultIdx, "b"] <- b
          printf("A/A = %d, m = %f, Q = %0.2f, g = %0.3f, b = %0.4f\n", resultIdx, m, Q, g, b)
          # Optimization optim_sa
          resultData[resultIdx, "optim_sa"]<- optim_sa(
            fun = objfunc,
            maximization = FALSE,
            start = c(Q, g, b),
            lower = c(0.01, 0.25, b_MIN),
            upper = c(0.2, 4, b_MAX),
            trace = FALSE, # TRUE
            control = list(
              t0 = 100,
              nlimit = 100,
              t_min = 0.1,
              dyn_rf = FALSE,
              rf = 1,
              r = 0.7
            )
          )[2]
          # Optimization optim_nm
          # resultData[resultIdx, "optim_nm"] <- optim_nm(
          #   fun = objfunc,
          #   k = 3,
          #   start = c(Q, g, b),
          #   maximum = FALSE,
          #   trace = TRUE, # TRUE
          #   alpha = 1,
          #   beta = 2,
          #   gamma = 1 / 2,
          #   delta = 1 / 2,
          #   tol = 0.00001,
          #   exit = 250,
          #   edge = 1
          # )[2]
          # Optimization optimx
          # resultData[resultIdx, "optimx"] <- optimx(
          #   par = c(Q, g, b),
          #   fn = objfunc,
          #   gr = NULL,
          #   hess = NULL,
          #   lower = c(0.01, 0.25, b_MIN),
          #   upper = c(0.2, 4, b_MAX),
          #   method = c("Nelder-Mead", "L-BFGS-B"),
          #   itnmax = NULL,
          #   hessian = FALSE
          # )[1,4]
          resultIdx <- resultIdx + 1
        }
      }
    }
  }
  # find method with minimal value
  indx <- 1 # which.min(c(min(resultData$`optim_sa`), min(resultData$`optim_nm`), min(resultData$`optimx`)))
  methodName <- c("optim_sa", "optim_nm", "optimx")[indx]
  # sort resultData in ascending order in column <methodName>
  resultData <- resultData[order(resultData[,methodName]),]
  # create result list <m, Q, g, b, "methodname">
  res <- c(resultData$m[1], resultData$Q[1], resultData$g[1], resultData$b[1], indx)
  return(res)
}
#
### 2.10 FIND OPTIMUM VALUES FOR (Q, g, b) THAT MINIMIZE ###
###      THE OBJECTIVE FUNCTION ###
#
findOptimumValues <- function(x, rfmtcReadyDF, b_MIN, b_MAX) {
  # print function display header
  printf("findOptimumValues()")
  m <- x[1]
  Q <- x[2]
  g <- x[3]
  b <- x[4]
  methodIndx <- x[5]
  # construct objective function
  objfunc <- make.OBJECTIVE(m, rfmtcReadyDF)
  # choose optimization method and run it
  if (methodIndx == 1) {
    optim_sa_res <- optim_sa(
      fun = objfunc,
      maximization = FALSE,
      start = c(Q, g, b),
      lower = c(0.01, 0.25, b_MIN),
      upper = c(0.2, 4, b_MAX),
      trace = FALSE, # TRUE
      control = list(
        t0 = 100,
        nlimit = 100,
        t_min = 0.1,
        dyn_rf = FALSE,
        rf = 1,
        r = 0.7
      )
    ) 
    par <- optim_sa_res$par
    vector_Q_g_b = c(par[1], par[2], par[3])
  } else if (methodIndx == 2) {
    optim_nm_res <- optim_nm(
      fun = objfunc,
      k = 3,
      start = c(Q, g, b),
      maximum = FALSE,
      trace = FALSE, # TRUE
      alpha = 1,
      beta = 2,
      gamma = 1 / 2,
      delta = 1 / 2,
      tol = 0.00001,
      exit = 200,
      edge = 1
    )
    par <- optim_nm_res$par
    vector_Q_g_b = c(par[1], par[2], par[3])
  } else {
    optimx_res <- optimx(
      par = c(Q, g, b),
      fn = objfunc,
      gr = NULL,
      hess = NULL,
      lower = c(0.01, 0.25, b_MIN),
      upper = c(0.2, 4, b_MAX),
      method = c("Nelder-Mead", "L-BFGS-B"),
      itnmax = NULL,
      hessian = FALSE
    )
    vector_Q_g_b = c(optimx_res$p1[1], optimx_res$p2[1], optimx_res$p3[1])
  }
  return(vector_Q_g_b)
}
#
### 2.11 CALCULATION OF E[X|L=1] WITH OPTIMAL (Q, g, b) ###
#
exl1_OptimalCalc <- function(x, rfmtcReadyDF) {
  # print function display header
  printf("exl1_OptimalCalc()")
  Q <- x[1]
  g <- x[2]
  b <- x[3]
  # create func_E_X_L_1() with optimal (Q, g, b)
  exl1 <- make.func_E_X_L_1(Q, g, b)
  # calculate column E[X[L]|L=1]
  rfmtcReadyDF$`E[X[L]|L=1]` <-
    mapply(function(r, f, t) exl1(r, f, t),
           rfmtcReadyDF$`Recency (months)`,
           rfmtcReadyDF$`Frequency (times)`,
           rfmtcReadyDF$`Time (months)`
    )
  # sort dataframe on E[X[L]|L=1] in desc order
  rfmtcReadyDF <- rfmtcReadyDF[order(-rfmtcReadyDF$`E[X[L]|L=1]`), ]
  # calculate P[B]
  rfmtcReadyDF$`P[B]` <- calcMovingAverage(rfmtcReadyDF$`Churn (0/1)`, m)
  # calculate squares of residual errors
  rfmtcReadyDF$`(E[X|L=1]-P[B])^2` <-
    (rfmtcReadyDF$`E[X[L]|L=1]` - rfmtcReadyDF$`P[B]`) ^ 2
  return(rfmtcReadyDF)
}
#
### 2.12 CONVERT TRANSACTION DATASETS TO RFM-RFMTC CUSTOMER DATASETS ###
#
convertTransactionToRFMTCreadyDataset <-
  function(transSampleFile,
           custSampleFile,
           churnMasterFile) {
    # print function display header
    printf("convertTransactionToRFMTCreadyDataset()")
    # read transSampleFile
    transSDF <- read.table(transSampleFile, header = FALSE)
    # remove ID column (col-2) and CDs column (col-4) from sample data frame
    transSDF <- as.data.frame(cbind(transSDF[, c(1, 3, 5)]))
    # set the column names for data frame
    names(transSDF) <- c("CustID", "Date", "Amount")
    # read churn master file
    churnDF <- read.csv(churnMasterFile, header = TRUE)
    # set the column names for churn data frame
    names(churnDF) <- c("CustID", "Churn")
    # tranform transaction sample data frame to customer sample data frame
    custSDF <- calcAppendRFMTC(transSDF, churnDF)
    # write custSampleFile
    write.csv(custSDF, custSampleFile, row.names = FALSE)
    return(NULL)
  }
#
### 2.13 CALCULATE & APPEND R,F,M,T,C COLUMNS TO CUSTOMER DATASETS ###
#  names(transDF) is c("CustID", "Date", "Amount")
#  names(churnDF) is c("CustID", "Churn")
calcAppendRFMTC <- function(transDF, churnDF) {
  # set the max(transaction date)+1
  endDate <- as.Date("19980701", "%Y%m%d")
  # a) convert Date col from numeric to date type
  transDF[, "Date"] <-
    as.Date(as.character(transDF[, "Date"]), "%Y%m%d")
  # b) find the max transaction date per CustID
  t_maxDate <-
    aggregate(transDF[, "Date"], list(transDF[, "CustID"]), max)
  names(t_maxDate) <- c("CustID", "MaxDate")
  t_maxDate <- t_maxDate[order(t_maxDate[, "CustID"]), ]
  # c) find the min transaction date per CustID
  t_minDate <-
    aggregate(transDF[, "Date"], list(transDF[, "CustID"]), min)
  names(t_minDate) <- c("CustID", "MinDate")
  t_minDate <- t_minDate[order(t_minDate[, "CustID"]), ]
  # d) find the average Amount per CustID
  t_AvgAmount <-
    aggregate(transDF[, "Amount"], list(transDF[, "CustID"]), mean)
  names(t_AvgAmount) <- c("CustID", "AvgAmount")
  t_AvgAmount <- t_AvgAmount[order(t_AvgAmount[, "CustID"]), ]
  # e) find the Frequency
  t_Freq <-
    aggregate(transDF[, "CustID"], by = list(transDF[, "CustID"]), length)
  names(t_Freq) <- c("CustID", "Frequency")
  t_Freq <- t_Freq[order(t_Freq[, "CustID"]), ]
  # f) find the Churn
  churnDF <- churnDF[churnDF[, "CustID"] %in% transDF[, "CustID"], ]
  churnDF <- churnDF[order(churnDF[, "CustID"]), ]
  # sort both transDF and churnDF with ID in ascending order
  transDF <- transDF[order(transDF[, "CustID"]), ]
  # g) merge columns and convert Recency and Time to months!
  custDF <-
    as.data.frame(cbind(
      (endDate - t_maxDate[, "MaxDate"]),
      t_Freq[, "Frequency"],
      t_AvgAmount[, "AvgAmount"],
      (endDate - t_minDate[, "MinDate"]),
      churnDF[, "Churn"]
    ))
  names(custDF) <-
    c("Recency", "Frequency", "Monetary", "Time", "Churn")
  custDF[, "Recency"] <- custDF[, "Recency"] / 30
  custDF[, "Time"] <- custDF[, "Time"] / 30
  return(custDF)
}
#
### 3. MAIN PROGRAM ###
printf("%s", RUN_TITLE)
printf("-------------------------------")
## PART A
# convert transaction datasets to RFM-RFMTC customer datasets
convertTransactionToRFMTCreadyDataset(
  TransSampleDSFile,
  CustSampleDSFile,
  ChurnMasterFile
)
## PART B
# partition customer dataset to training & test datasets
RunFileName <-
  partitionCustomerDSToTrainingTestDS(
    CustSampleDSFile,
    TrainingDSFile,
    TestDSFile,
    RUN_TYPE
  )
## PART C
# RFMTC dataset preparation
rfmtcReadyDF <- prepareRFMTCdataset(RunFileName)
if (RUN_TYPE == "TRAINING") {
  # find range constraints for parameter (b)
  resList <- calc_b_RangeConstraints(rfmtcReadyDF)
  rfmtcReadyDF <- resList[[1]]
  b_RangeDetailsList <- resList[[2]]
  # discover optimal set of initial values (m, Q, g, b)
  resultList <-
    discoverOptimalInitialParams(rfmtcReadyDF,
        6:6, #  4:4
        seq(0.01, 0.2, length.out = 3), # seq(0.01, 0.2, length.out = 3)
        seq(0.25, 4, length.out = 3),  # seq(0.25, 4, length.out = 3)
        seq(b_RangeDetailsList[["MIN"]], b_RangeDetailsList[["MAX"]], length.out = 3),
        b_RangeDetailsList[["MIN"]],
        b_RangeDetailsList[["MAX"]])
  # find optimum values for (Q, g, b) that minimize the objective function
  optimumParams <- findOptimumValues(resultList, 
                                     rfmtcReadyDF,
                                     b_RangeDetailsList[["MIN"]], 
                                     b_RangeDetailsList[["MAX"]])
  # save to external file the Optimal (m, Q, g, b)
  write.csv(
    c(resultList[1], optimumParams[1], optimumParams[2], optimumParams[3]), 
    TRAINING_OPTIMALS)
} 
# read from external file the Optimal (m, Q, g, b)
optimumParams <- read.csv(TRAINING_OPTIMALS)
m <- optimumParams$x[1]
Q <- optimumParams$x[2]
g <- optimumParams$x[3]
b <- optimumParams$x[4]
# calculation of E[X|L=1] with optimal (Q, g, b) & write data to file
write.csv(exl1_OptimalCalc(c(Q,g,b), rfmtcReadyDF), 
                          paste("OUR_RFMTC_", 
                          RUN_TYPE, ".csv", 
                          sep = ""), 
          row.names = FALSE)