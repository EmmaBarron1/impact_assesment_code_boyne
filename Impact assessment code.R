# /----------------------------------------------------------------------- #
#  This code is demonstration of conducting a climate change impacts assessment
#  by first calibrating/validating a conceptual rainfall runoff model with 
#  observed data and then forcing the model with the CORDEX ensemble of 11 
#  GCM/RCM runs with RCP8.5. Finally, code is presented for assessing changes 
#  in monthly flows and floods in the Boyne catchment. 

#Code written by Conor Murphy for GY655 Impacts and Policy
# -----------------------------------------------------------------------/ #

#install the GR4J model available in the airGR package and set WD
install.packages("airGR")
library(airGR)

setwd("C:/Users/cmmurphy/OneDrive - Maynooth University/Misc/Documents/Workshop")

# /----------------------------------------------------------------------- #
#  Prepapre input data for running GR4J model for Boyne at Slane Castle    #
# -----------------------------------------------------------------------/ #

# Read catchment data.
BasinObs      <- read.csv("BasinObs_7012.csv")
BasinObs$Date <- as.POSIXct(BasinObs$Date, format = "%d/%m/%Y")

# Define warm-up period.
IndPeriod_WarmUp <- seq(which(format(BasinObs$Date, format = "%Y-%m-%d") == "1987-01-01"),
                        which(format(BasinObs$Date, format = "%Y-%m-%d") == "1987-12-31"))

# Prepare input data.
InputsModel <- CreateInputsModel(
  FUN_MOD = RunModel_GR4J,
  DatesR  = BasinObs$Date,
  Precip  = BasinObs$P,
  PotEvap = BasinObs$E
)

# /----------------------------------------------------------------------- #
#  Calibration- set up run options and run calibration                                                     #
# -----------------------------------------------------------------------/ #

# Define calibration period.
IndPeriod_Cal <- seq(which(format(BasinObs$Date, format = "%Y-%m-%d") == "1988-01-01"),
                     which(format(BasinObs$Date, format = "%Y-%m-%d") == "2002-12-31"))

# Prepare calibration run options.
RunOptions_Cal <- CreateRunOptions(
  FUN_MOD          = RunModel_GR4J,
  InputsModel      = InputsModel,
  IndPeriod_WarmUp = IndPeriod_WarmUp,
  IndPeriod_Run    = IndPeriod_Cal
)

# Prepare objective function.
InputsCrit_Cal <- CreateInputsCrit(
  FUN_CRIT    = "ErrorCrit_NSE",
  InputsModel = InputsModel,
  RunOptions  = RunOptions_Cal,
  Obs         = BasinObs$Qmm[IndPeriod_Cal]
)

# Prepare calibration options.
CalibOptions <- CreateCalibOptions(
  FUN_MOD   = RunModel_GR4J,
  FUN_CALIB = "Calibration_Michel"
)

# Run calibration procedure.
OutputsCalib <- Calibration_Michel(
  FUN_MOD      = RunModel_GR4J,
  InputsModel  = InputsModel,
  RunOptions   = RunOptions_Cal,
  InputsCrit   = InputsCrit_Cal,
  CalibOptions = CalibOptions
)

# Get optimum parameters.
Param <- OutputsCalib$ParamFinalR

# Run the model for the calibration period (needed to get initial states for validation).
OutputsModel_Cal <- RunModel_GR4J(
  InputsModel = InputsModel,
  RunOptions  = RunOptions_Cal,
  Param       = Param
)

# /---------------------------------------------------------------------- #
#  Validation - set up run options and validate                           #
# ----------------------------------------------------------------------/ #

# Define validation period.
IndPeriod_Val <- seq(which(format(BasinObs$Date, format = "%Y-%m-%d") == "2003-01-01"),
                     which(format(BasinObs$Date, format = "%Y-%m-%d") == "2017-12-31"))

# Prepare validation run options.
RunOptions_Val <- CreateRunOptions(
  FUN_MOD          = RunModel_GR4J,
  InputsModel      = InputsModel,
  IndPeriod_WarmUp = 0L,
  IniStates        = OutputsModel_Cal$StateEnd,
  IndPeriod_Run    = IndPeriod_Val
)

# Prepare objective function.
InputsCrit_Val <- CreateInputsCrit(
  FUN_CRIT    = "ErrorCrit_NSE",
  InputsModel = InputsModel,
  RunOptions  = RunOptions_Val,
  Obs         = BasinObs$Qmm[IndPeriod_Val]
)

# Run the model for the validation period.
OutputsModel_Val <- RunModel_GR4J(
  InputsModel = InputsModel,
  RunOptions  = RunOptions_Val,
  Param       = Param
)

# Calculate efficiency.
OutputsCrit_Val <- ErrorCrit_NSE(
  InputsCrit   = InputsCrit_Val,
  OutputsModel = OutputsModel_Val
)

# /---------------------------------------------------------------------- #
#  Full simulation                                                        #
# ----------------------------------------------------------------------/ #

# Define run period.
IndPeriod_Run <- seq(which(format(BasinObs$Date, format = "%Y-%m-%d") == "1988-01-01"),
                     which(format(BasinObs$Date, format = "%Y-%m-%d") == "2017-12-31"))

# Prepare run options.
RunOptions <- CreateRunOptions(
  FUN_MOD          = RunModel_GR4J,
  InputsModel      = InputsModel,
  IndPeriod_WarmUp = IndPeriod_WarmUp,
  IndPeriod_Run    = IndPeriod_Run
)

# Prepare objective function.
InputsCrit <- CreateInputsCrit(
  FUN_CRIT    = "ErrorCrit_NSE",
  InputsModel = InputsModel,
  RunOptions  = RunOptions,
  Obs         = BasinObs$Qmm[IndPeriod_Run]
)

# Run the model.
OutputsModel <- RunModel_GR4J(
  InputsModel = InputsModel,
  RunOptions  = RunOptions,
  Param       = Param
)

# Calculate efficiency.
OutputsCrit <- ErrorCrit_NSE(
  InputsCrit   = InputsCrit,
  OutputsModel = OutputsModel
)

# Plot results.
plot.OutputsModel(OutputsModel, Qobs = BasinObs$Qmm[IndPeriod_Run])


###############################################################################


# Future Simulations using CORDEX RCP8.5
# Just load the RCP4.5 scenarios instead to run different ES

###############################################################################


# Dates (1976-2100).
DatesR <- seq(as.POSIXlt("1976-01-01"), as.POSIXlt("2100-12-31"), by = "days")

#DatesR <- DatesR[-which(format(DatesR, "%m-%d") == "02-29")] # Remove leap days.

#Import CORDEX ensemble data
Precip  <- read.csv("pr_07012_cordex_rcp85.csv")
Precip<-as.data.frame(Precip)
Temp    <- read.csv("tas_07012_cordex_rcp85.csv")

#Calculate PET Oudin

PotEvap <- NULL

for (i in colnames(Temp)) {

  PE <- PE_Oudin(JD = as.POSIXlt(DatesR)$yday + 1,
                          Temp = Temp[[i]],
                          Lat = 53.7, LatUnit = "deg")
  
  PotEvap<-cbind(PotEvap, PE)
}

colnames(PotEvap)<-c("m1", "m2","m3", "m4","m5", "m6","m7", "m8","m9", "m10","m11" )
PotEvap<-as.data.frame(PotEvap)

  
  # Run model for each ensemble member.
  Q <- data.frame(mapply(function(Precip, PotEvap) {
    
    # Create InputsModel.
    InputsModel <- CreateInputsModel(
      FUN_MOD = RunModel_GR4J,
      DatesR  = DatesR,
      Precip  = Precip,
      PotEvap = PotEvap
    )
    
    # Create RunOptions.
    RunOptions <- CreateRunOptions(
      FUN_MOD          = RunModel_GR4J,
      InputsModel      = InputsModel,
      IndPeriod_WarmUp = 1:365,
      IndPeriod_Run    = 1:length(DatesR)
    )
    
    # Run GR4J.
    OutputsModel <- RunModel_GR4J(
      InputsModel = InputsModel,
      RunOptions  = RunOptions,
      Param       = Param
    )
    
    # Return model output.
    return(OutputsModel$Qsim)
    
  }, Precip, PotEvap))
  
  # Add dates and parameter set.
  Q.sim           <- cbind(as.character(DatesR[1:length(DatesR)]), Q)
  colnames(Q.sim) <- c("Date", paste0("Q_Ens_", 1:11))
  
  write.table(Q.sim, file="Q.sim.csv")
  
  ###############################################################################
  
  # Calculate Percent change for future time periods and plot
  
  ###############################################################################
  
  #first add a year month day column
  data<-Q.sim
  # Create three 'empty' columns (to put seperated dates)
  #cbind() function combines vector, matrix or data frame by columns
  data <- cbind(NA, NA, NA, data)
  
  # Take a look!
  head(data)
  tail(data)
  
  # Convert 'Date' column to R date object using 'as.Date'
  data$Date <- as.Date(data$Date)
  
  # Check has this worked with 'str' (Displays structure of R objects) 
  str(data)
  
  # Split 'Date' column into three seperate columns
  data[ ,1] <- as.numeric(format(data$Date, format = "%d")) # Day
  data[ ,2] <- as.numeric(format(data$Date, format = "%m")) # Month 
  data[ ,3] <- as.numeric(format(data$Date, format = "%Y")) # Year 
  
  
  # Set column names
  colnames(data) <- c("Day", "Month", "Year", "Date", paste0("Q_Ens_", 1:11))
  
  # Take a look!
  head(data)
  
  #split into reference and future time periods
  ref<-data[data$Date >= "1976-01-01" & data$Date <= "2005-12-31", ]
  near<-data[data$Date >= "2010-01-01" & data$Date <= "2039-12-31", ]
  mid<-data[data$Date >= "2040-01-01" & data$Date <= "2069-12-31", ]
  far<-data[data$Date >= "2070-01-01" & data$Date <= "2099-12-31", ]
  
  # calculate monthly means for each period
  ref.month <- aggregate(ref[,5:15], by = list(ref$Month), FUN = mean, na.rm = TRUE)
  near.month <- aggregate(near[,5:15], by = list(near$Month), FUN = mean, na.rm = TRUE)
  mid.month <- aggregate(mid[,5:15], by = list(mid$Month), FUN = mean, na.rm = TRUE)
  far.month <- aggregate(far[,5:15], by = list(far$Month), FUN = mean, na.rm = TRUE)

  #calculate percent change
  near.diff<-(((near.month[,2:12]- ref.month[,2:12])/ref.month[,2:12])*100)
  near.med<-apply(near.diff,1,function(x) quantile(x, probs = .5, na.rm=TRUE))
  mid.diff<-(((mid.month[,2:12]- ref.month[,2:12])/ref.month[,2:12])*100)
  mid.med<-apply(mid.diff,1,function(x) quantile(x, probs = .5, na.rm=TRUE))
  far.diff<-(((far.month[,2:12]- ref.month[,2:12])/ref.month[,2:12])*100)
  far.med<-apply(far.diff,1,function(x) quantile(x, probs = .5, na.rm=TRUE))
  
  
  #Create Plot of Percent change for each future period
  
  #plot near
  plot(mid.month$Group.1, near.med,
       type="b",
       col="red",
       lwd=2,
       main="Near Future 2020s",
       ylim=c(-50, 80),
       xlab="Month",
       ylab= "Percent change",
       xlim=c(1, 12))
  
  abline(h=0)
  
  for (i in colnames(near.diff[,1:11])) {
  
  lines(ref.month$Group.1, near.diff[[i]], col="grey", type="b")
    
  }
  
  #plot mid
  plot(mid.month$Group.1, mid.med,
       type="b",
       col="red",
       lwd=2,
       main="Mid Future 2050s",
       ylim=c(-50, 80),
       xlab="Month",
       ylab= "Percent change",
       xlim=c(1, 12))
  
  abline(h=0)
  
  for (i in colnames(mid.diff[,1:11])) {
    
    lines(ref.month$Group.1, mid.diff[[i]], col="grey", type="b")
    
  }
  
  
  #plot far
  plot(mid.month$Group.1, far.med, 
       type="b",
       col="red",
       lwd=2,
       main="Far Future 2080s",
       ylim=c(-50, 80),
       xlab="Month",
       ylab= "Percent change",
       xlim=c(1, 12))
  
  abline(h=0)
  
  for (i in colnames(far.diff[,1:11])) {
    
    lines(ref.month$Group.1, far.diff[[i]], col="grey", type="b")
    
  }
  
#################################################################################
################################################################################
  
  #extreme value analysis to examine changes in floods of different return levels
  
  #--------------------------------------------------------------------------------
  
  # Load required libraries
  install.packages("dplyr")
  library(dplyr)
  install.packages("ggplot2")
  library(ggplot2)
  install.packages("extRemes")
  library(extRemes) # For GEV fitting
  
  # Load your data (adjust the path as needed)
  data <- Q
  head(data)
  
  # Insert a date column with a daily sequence from 1st January 1976 to 31st December 2100
  start_date <- as.Date("1976-01-01")
  end_date <- as.Date("2100-12-31")
  data$Date <- seq.Date(from = start_date, to = end_date, by = "day")
  
  # Calculate the average flow across all ensemble members
  data$Flow <- rowMeans(data[, 1:11], na.rm = TRUE)
  
  # Define periods
  ref <- data %>% filter(Date >= as.Date("1976-01-01") & Date <= as.Date("2005-12-31"))
  near <- data %>% filter(Date >= as.Date("2010-01-01") & Date <= as.Date("2039-12-31"))
  mid <- data %>% filter(Date >= as.Date("2040-01-01") & Date <= as.Date("2069-12-31"))
  far <- data %>% filter(Date >= as.Date("2070-01-01") & Date <= as.Date("2099-12-31"))
  
  
  # Find the annual maxima for the reference period
  annual_maxima <- ref %>%
    mutate(year = format(Date, "%Y")) %>%
    group_by(year) %>%
    summarise(max_flow = max(Flow, na.rm = TRUE))
  
  # Fit the GEV model
  gev_fit <- fevd(annual_maxima$max_flow, type = "GEV")
  
 ############################################################################# 
  
  # Diagnostic plots for GEV fit during reference period
  par(mfrow = c(2, 2))  # Arrange plots in a 2x2 grid
  plot(gev_fit)  # This will generate the probability plot, Q-Q plot, return level plot, and density plot
  par(mfrow = c(1, 1))  # Reset plotting layout
  
  summary(gev_fit)
  
  ###########################################################################
  
  # Function to fit GEV and calculate return levels
  fit_gev_and_return_levels <- function(flow_data, return_periods = c(2, 10, 20, 50, 100)) {
    # Find the annual maxima
    annual_maxima <- flow_data %>%
      mutate(year = format(Date, "%Y")) %>%
      group_by(year) %>%
      summarise(max_flow = max(Flow, na.rm = TRUE))
    
    # Fit the GEV model
    gev_fit <- fevd(annual_maxima$max_flow, type = "GEV")
    
    
    # Calculate return levels for given return periods
    return_levels <- return.level(gev_fit, return.period = return_periods)
    
    return(return_levels)
  }
  
  # Fit GEV and get return levels for each period
  return_periods <- c(2, 10, 20, 50, 100)
  ref_return_levels <- fit_gev_and_return_levels(ref, return_periods)
  near_return_levels <- fit_gev_and_return_levels(near, return_periods)
  mid_return_levels <- fit_gev_and_return_levels(mid, return_periods)
  far_return_levels <- fit_gev_and_return_levels(far, return_periods)
  
  # Combine return levels into a dataframe for plotting
  return_levels_df <- data.frame(
    Period = rep(c("1976-2005", "2010-2039", "2040-2069", "2070-2099"), each = length(return_periods)),
    Return_Period = rep(return_periods, 4),
    Return_Level = c(ref_return_levels, near_return_levels, mid_return_levels, far_return_levels)
  )
  
  # Plot the return levels
  ggplot(return_levels_df, aes(x = Return_Period, y = Return_Level, color = Period)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_x_log10() +
    labs(title = "Return Levels for Different Periods",
         x = "Return Period (years)",
         y = "Flood Return Level",
         color = "Period") +
    theme_minimal()

  ##############################################################################
  ##############################################################################