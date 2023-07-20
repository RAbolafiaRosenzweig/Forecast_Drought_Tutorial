##Tutorial for drought prediction for a single pixel at 1-, 2- and 3-month lead times


#load in libraries
library(scales)
library(nlme)
library(mgcv)

#define directory on your computer that the "Tutorial" folder is brought to:
global_dir = "/Users/abolafia/ForecastDrought/"

#directory where hydroclimate and drought data are located (1 file for each 4km pixel with time series of each variable corresponding to dates in datelist):
in_data_dir = paste(global_dir,"/Tutorial/", sep = "")

##================================================================================================================================================
##================================================================================================================================================
##PART 1: This first part of the script goes over how to find the best set of predictors for summer drought
##================================================================================================================================================
##================================================================================================================================================

#Load in time range:
datelist = read.table( paste(in_data_dir,"monthly_datelist.csv",sep=""),sep="," )

#define pre-summer data 
idx_pre_summer = which(datelist[,2] == 11  | datelist[,2] == 12  | datelist[,2] == 1  | datelist[,2] == 2) #3-month lead time
#idx_pre_summer = which(datelist[,2] == 11  | datelist[,2] == 12  | datelist[,2] == 1  | datelist[,2] == 2 |  datelist[,2] == 3) #2-month lead time
#idx_pre_summer = which(datelist[,2] == 11  | datelist[,2] == 12  | datelist[,2] == 1  | datelist[,2] == 2 |  datelist[,2] == 3 |  datelist[,2] == 4) #1-month lead time
datelist_pre_summer = datelist[idx_pre_summer,]

#define summer dates (June-Aug):
idx_summer = which(datelist[,2] == 6 | datelist[,2] == 7  | datelist[,2] == 8)
datelist_summer = datelist[idx_summer,]

#define years:
years = unique(datelist[,1])
nyears = length(years)

#define list of coordinates there is data for:
Coords = read.table( paste(in_data_dir,"/CONUS404_coords.csv",sep="") , sep=",")

#for the tutorial choose a random coordinate from the WUS:
c = 20000  #good example of a pixel with weak drought memory, where the modeling approach provides strong improvement vs. auto R
#c = 30000 #good example where model primarily relies on auto R to get accurate prediction of drought
#other data are saved locally on Ronnie's storage device here: /Volumes/Pruina_External_Elements/ForecastDrought/Data/Covariates_CONUS404_Grid/

lat = Coords[c,1]
lon = Coords[c,2]
in_data_filename = sprintf("Covariates_%.4f_%.4f.csv",lat,lon)
read_filename = paste(in_data_dir,in_data_filename,sep="")
current_data = read.table(read_filename,sep=",")
nvars = length(current_data[1,])

#calculate seasonally averaged climate conditions:
store_pre_summer_data = c()
store_summer_data = c()
for (y in 1:nyears){
  current_year = years[y]
  idx_current_pre_summer = idx_pre_summer[which(datelist_pre_summer[,1] == current_year)]
  idx_current_summer = idx_summer[which(datelist_summer[,1] == current_year)]
  
  current_pre_summer_data = current_data[idx_current_pre_summer,]
  current_Summer_data = current_data[idx_current_summer,]
  
  #define average pre_summer conditions:
  mean_pre_summer = colMeans(current_pre_summer_data)
  #define average summer conditions:
  mean_summer = colMeans(current_Summer_data)
  
  #store_data (monthly averages of pre-summer and summer data w/ 1 data point for each year):
  store_pre_summer_data = rbind(store_pre_summer_data,mean_pre_summer)
  store_summer_data = rbind(store_summer_data,mean_summer)
}


#aggregate predictor data into dataframes:
pre_summer_PDO = store_pre_summer_data[,1]
pre_summer_AMO = store_pre_summer_data[,2]
pre_summer_SMprctile = store_pre_summer_data[,3]
pre_summer_PDSI = store_pre_summer_data[,4]
pre_summer_Temp = store_pre_summer_data[,5]
pre_summer_VPD = store_pre_summer_data[,6]
pre_summer_prcp = store_pre_summer_data[,7]
pre_summer_ET = store_pre_summer_data[,8]
pre_summer_PET = store_pre_summer_data[,9]
pre_summer_SWE = store_pre_summer_data[,14]

#bring SM-percentile to percentile space (0-1) 
#using eq. 2 from: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2022WR033734
sorted_idx = order(pre_summer_SMprctile)
sorted = pre_summer_SMprctile[sorted_idx]
rank=1:length(sorted_idx)
rank[sorted_idx] = rank
pre_summer_SMprctile = (rank - min(rank))/( max(rank) - min(rank))
pre_summer_SMprctile = rescale(pre_summer_SMprctile,to=c(0,1))

#replace NaNs w/0 for SWE
idx_na<-which(is.na(pre_summer_SWE))
pre_summer_SWE[idx_na] = 0
idx_zero = which(pre_summer_SWE==0)

#this determines whether or not to use SWE as a predictor
if (max(pre_summer_SWE) == min(pre_summer_SWE) ){ #if there is no temporal variation in SWE, do not include it as a predictor
  Covariates_DF_all = data.frame(pre_summer_PDO=pre_summer_PDO,pre_summer_AMO=pre_summer_AMO,pre_summer_SMprctile=pre_summer_SMprctile,pre_summer_PDSI=pre_summer_PDSI,pre_summer_Temp=pre_summer_Temp,pre_summer_VPD=pre_summer_VPD,pre_summer_prcp=pre_summer_prcp,pre_summer_ET=pre_summer_ET,pre_summer_PET=pre_summer_PET)
}else if (max(pre_summer_SWE)<0.1) { #if SWE is very small in all years do not include it as a predictor
  Covariates_DF_all = data.frame(pre_summer_PDO=pre_summer_PDO,pre_summer_AMO=pre_summer_AMO,pre_summer_SMprctile=pre_summer_SMprctile,pre_summer_PDSI=pre_summer_PDSI,pre_summer_Temp=pre_summer_Temp,pre_summer_VPD=pre_summer_VPD,pre_summer_prcp=pre_summer_prcp,pre_summer_ET=pre_summer_ET,pre_summer_PET=pre_summer_PET)
}else if (length(idx_zero)>=10) { #if there is 0 SWE for at least 10 years, do not include it as a predictor
  Covariates_DF_all = data.frame(pre_summer_PDO=pre_summer_PDO,pre_summer_AMO=pre_summer_AMO,pre_summer_SMprctile=pre_summer_SMprctile,pre_summer_PDSI=pre_summer_PDSI,pre_summer_Temp=pre_summer_Temp,pre_summer_VPD=pre_summer_VPD,pre_summer_prcp=pre_summer_prcp,pre_summer_ET=pre_summer_ET,pre_summer_PET=pre_summer_PET)
}else{ #if the above conditions do not apply, use SWE as a predictor
  Covariates_DF_all = data.frame(pre_summer_PDO=pre_summer_PDO,pre_summer_AMO=pre_summer_AMO,pre_summer_SMprctile=pre_summer_SMprctile,pre_summer_PDSI=pre_summer_PDSI,pre_summer_Temp=pre_summer_Temp,pre_summer_VPD=pre_summer_VPD,pre_summer_prcp=pre_summer_prcp,pre_summer_ET=pre_summer_ET,pre_summer_PET=pre_summer_PET,pre_summer_SWE=pre_summer_SWE)
}

#make sure there are no all 0 predictors, which will cause error in PCA function:
DIM = dim(Covariates_DF_all)
npreds = DIM[2]
for (ipred in 1:npreds){
  current_pred = Covariates_DF_all[,ipred]
  IDX_nan<-which(is.na(current_pred))
  current_pred[IDX_nan] = 0
  if ( min(current_pred)==0 & max(current_pred)==0){
    current_pred=runif(length(current_pred), 0, 1) #reset to randomized vector
    Covariates_DF_all[,ipred]=current_pred
  } else{
    Covariates_DF_all[,ipred]=current_pred
  }
}

#summer drought conditions to predict:
summer_SMprctile = as.vector(store_summer_data[,3])
summer_PDSI = as.vector(store_summer_data[,4])


#bring SM-percentile to percentile space (0-1) 
#using eq. 2 from: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2022WR033734
sorted_idx = order(summer_SMprctile)
sorted = summer_SMprctile[sorted_idx]
rank=1:length(sorted_idx)
rank[sorted_idx] = rank
summer_SMprctile = (rank - min(rank))/( max(rank) - min(rank))
summer_SMprctile = rescale(summer_SMprctile,to=c(0,1))

##===============FIND BEST MODEL LOOPING THROUGH ALL POSSIBLE COMBINATIONS OF PREDICTORS#===============
#THIS USES GCV AS CRITERIA
#define best N-member ensemble
nvars = npreds
x=1:nvars
concurvity_threshold = 0.4 #this is applied to reduce overfitting issues (ie non-linear relationships between predictors), we only use models w/ low concurvity (more infor here: https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/concurvity.html)
store_GCV_PDSI = c()
store_GCV_SMprctile = c()
for (n_predictors in 1:nvars){
  print(c(c,n_predictors))
  ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
  combo_IDs<-combn(x,n_predictors)
  for (i in 1:ncombos){
    if (n_predictors == 1){
      predictor1 <-Covariates_DF_all[,combo_IDs[1,i]]
      
      #PCA combo model for predicting PDSI:
      Covariates_DF = data.frame(predictor1=predictor1)
      df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
      pcs = df.pca$x
      pc1 = pcs[,1]
      current_DF_pcs <- data.frame(Y=summer_PDSI,predictor1=pc1)
      mod_PDSI <- gam(data=current_DF_pcs,Y~s(predictor1),family="gaussian",control=list(mgcv.tol=1e-4))
      
      #PCA combo model for predicting SM-prctile:
      current_DF_pcs <- data.frame(Y=summer_SMprctile,predictor1=pc1)
      mod_SMprctile <- gam(data=current_DF_pcs,Y~s(predictor1),family="gaussian",control=list(mgcv.tol=1e-4))
    }
    
    if (n_predictors == 2){
      predictor1 <-Covariates_DF_all[,combo_IDs[1,i]]
      predictor2 <-Covariates_DF_all[,combo_IDs[2,i]]
      
      #PCA combo model for predicting PDSI:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2)
      df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
      pcs = df.pca$x
      pc1 = pcs[,1]
      pc2 = pcs[,2]
      current_DF_pcs <- data.frame(Y=summer_PDSI,predictor1=pc1,predictor2=pc2)
      mod_PDSI <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2),family="gaussian",control=list(mgcv.tol=1e-4))
      
      #PCA combo model for predicting SM-prctile:
      current_DF_pcs <- data.frame(Y=summer_SMprctile,predictor1=pc1,predictor2=pc2)
      mod_SMprctile <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2),family="gaussian",control=list(mgcv.tol=1e-4))
    }
    
    if (n_predictors == 3){
      predictor1 <-Covariates_DF_all[,combo_IDs[1,i]]
      predictor2 <-Covariates_DF_all[,combo_IDs[2,i]]
      predictor3 <-Covariates_DF_all[,combo_IDs[3,i]]
      
      #PCA combo model for predicting PDSI:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3)
      df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
      pcs = df.pca$x
      pc1 = pcs[,1]
      pc2 = pcs[,2]
      pc3 = pcs[,3]
      current_DF_pcs <- data.frame(Y=summer_PDSI,predictor1=pc1,predictor2=pc2,predictor3=pc3)
      mod_PDSI <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3),family="gaussian",control=list(mgcv.tol=1e-4))
      
      #PCA combo model for predicting SM-prctile:
      current_DF_pcs <- data.frame(Y=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3)
      mod_SMprctile <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3),family="gaussian",control=list(mgcv.tol=1e-4))
    }
    
    if (n_predictors == 4){
      predictor1 <-Covariates_DF_all[,combo_IDs[1,i]]
      predictor2 <-Covariates_DF_all[,combo_IDs[2,i]]
      predictor3 <-Covariates_DF_all[,combo_IDs[3,i]]
      predictor4 <-Covariates_DF_all[,combo_IDs[4,i]]
      
      #PCA combo model for predicting PDSI:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4)
      df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
      pcs = df.pca$x
      pc1 = pcs[,1]
      pc2 = pcs[,2]
      pc3 = pcs[,3]
      pc4 = pcs[,4]
      current_DF_pcs <- data.frame(Y=summer_PDSI,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4)
      mod_PDSI <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
      
      #PCA combo model for predicting SM-prctile:
      current_DF_pcs <- data.frame(Y=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4)
      mod_SMprctile <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
    }
    
    if (n_predictors == 5){
      predictor1 <-Covariates_DF_all[,combo_IDs[1,i]]
      predictor2 <-Covariates_DF_all[,combo_IDs[2,i]]
      predictor3 <-Covariates_DF_all[,combo_IDs[3,i]]
      predictor4 <-Covariates_DF_all[,combo_IDs[4,i]]
      predictor5 <-Covariates_DF_all[,combo_IDs[5,i]]
      
      #PCA combo model for predicting PDSI:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5)
      df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
      pcs = df.pca$x
      pc1 = pcs[,1]
      pc2 = pcs[,2]
      pc3 = pcs[,3]
      pc4 = pcs[,4]
      pc5 = pcs[,5]
      current_DF_pcs <- data.frame(Y=summer_PDSI,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5)
      mod_PDSI <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
      
      #PCA combo model for predicting SM-prctile:
      current_DF_pcs <- data.frame(Y=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5)
      mod_SMprctile <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
    }
    
    if (n_predictors == 6){
      predictor1 <-Covariates_DF_all[,combo_IDs[1,i]]
      predictor2 <-Covariates_DF_all[,combo_IDs[2,i]]
      predictor3 <-Covariates_DF_all[,combo_IDs[3,i]]
      predictor4 <-Covariates_DF_all[,combo_IDs[4,i]]
      predictor5 <-Covariates_DF_all[,combo_IDs[5,i]]
      predictor6 <-Covariates_DF_all[,combo_IDs[6,i]]
      
      #PCA combo model for predicting PDSI:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6)
      df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
      pcs = df.pca$x
      pc1 = pcs[,1]
      pc2 = pcs[,2]
      pc3 = pcs[,3]
      pc4 = pcs[,4]
      pc5 = pcs[,5]
      pc6 = pcs[,6]
      current_DF_pcs <- data.frame(Y=summer_PDSI,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6)
      mod_PDSI <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
      
      #PCA combo model for predicting SM-prctile:
      current_DF_pcs <- data.frame(Y=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6)
      mod_SMprctile <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
    }
    
    if (n_predictors == 7){
      predictor1 <-Covariates_DF_all[,combo_IDs[1,i]]
      predictor2 <-Covariates_DF_all[,combo_IDs[2,i]]
      predictor3 <-Covariates_DF_all[,combo_IDs[3,i]]
      predictor4 <-Covariates_DF_all[,combo_IDs[4,i]]
      predictor5 <-Covariates_DF_all[,combo_IDs[5,i]]
      predictor6 <-Covariates_DF_all[,combo_IDs[6,i]]
      predictor7 <-Covariates_DF_all[,combo_IDs[7,i]]
      
      #PCA combo model for predicting PDSI:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7)
      df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
      pcs = df.pca$x
      pc1 = pcs[,1]
      pc2 = pcs[,2]
      pc3 = pcs[,3]
      pc4 = pcs[,4]
      pc5 = pcs[,5]
      pc6 = pcs[,6]
      pc7 = pcs[,7]
      current_DF_pcs <- data.frame(Y=summer_PDSI,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7)
      mod_PDSI <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
      
      #PCA combo model for predicting SM-prctile:
      current_DF_pcs <- data.frame(Y=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7)
      mod_SMprctile <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
    }
    
    if (n_predictors == 8){
      predictor1 <-Covariates_DF_all[,combo_IDs[1,i]]
      predictor2 <-Covariates_DF_all[,combo_IDs[2,i]]
      predictor3 <-Covariates_DF_all[,combo_IDs[3,i]]
      predictor4 <-Covariates_DF_all[,combo_IDs[4,i]]
      predictor5 <-Covariates_DF_all[,combo_IDs[5,i]]
      predictor6 <-Covariates_DF_all[,combo_IDs[6,i]]
      predictor7 <-Covariates_DF_all[,combo_IDs[7,i]]
      predictor8 <-Covariates_DF_all[,combo_IDs[8,i]]
      
      #PCA combo model for predicting PDSI:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8)
      df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
      pcs = df.pca$x
      pc1 = pcs[,1]
      pc2 = pcs[,2]
      pc3 = pcs[,3]
      pc4 = pcs[,4]
      pc5 = pcs[,5]
      pc6 = pcs[,6]
      pc7 = pcs[,7]
      pc8 = pcs[,8]
      current_DF_pcs <- data.frame(Y=summer_PDSI,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8)
      mod_PDSI <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
      
      #PCA combo model for predicting SM-prctile:
      current_DF_pcs <- data.frame(Y=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8)
      mod_SMprctile <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
    }
    
    if (n_predictors == 9){
      predictor1 <-Covariates_DF_all[,combo_IDs[1,i]]
      predictor2 <-Covariates_DF_all[,combo_IDs[2,i]]
      predictor3 <-Covariates_DF_all[,combo_IDs[3,i]]
      predictor4 <-Covariates_DF_all[,combo_IDs[4,i]]
      predictor5 <-Covariates_DF_all[,combo_IDs[5,i]]
      predictor6 <-Covariates_DF_all[,combo_IDs[6,i]]
      predictor7 <-Covariates_DF_all[,combo_IDs[7,i]]
      predictor8 <-Covariates_DF_all[,combo_IDs[8,i]]
      predictor9 <-Covariates_DF_all[,combo_IDs[9,i]]
      
      #PCA combo model for predicting PDSI:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8,predictor9=predictor9)
      df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
      pcs = df.pca$x
      pc1 = pcs[,1]
      pc2 = pcs[,2]
      pc3 = pcs[,3]
      pc4 = pcs[,4]
      pc5 = pcs[,5]
      pc6 = pcs[,6]
      pc7 = pcs[,7]
      pc8 = pcs[,8]
      pc9 = pcs[,9]
      
      current_DF_pcs <- data.frame(Y=summer_PDSI,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8,predictor9=pc9)
      mod_PDSI <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
      
      #PCA combo model for predicting SM-prctile:
      current_DF_pcs <- data.frame(Y=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8,predictor9=pc9)
      mod_SMprctile <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
    }
    
    if (n_predictors == 10){
      predictor1 <-Covariates_DF_all[,combo_IDs[1,i]]
      predictor2 <-Covariates_DF_all[,combo_IDs[2,i]]
      predictor3 <-Covariates_DF_all[,combo_IDs[3,i]]
      predictor4 <-Covariates_DF_all[,combo_IDs[4,i]]
      predictor5 <-Covariates_DF_all[,combo_IDs[5,i]]
      predictor6 <-Covariates_DF_all[,combo_IDs[6,i]]
      predictor7 <-Covariates_DF_all[,combo_IDs[7,i]]
      predictor8 <-Covariates_DF_all[,combo_IDs[8,i]]
      predictor9 <-Covariates_DF_all[,combo_IDs[9,i]]
      predictor10 <-Covariates_DF_all[,combo_IDs[10,i]]
      #PCA combo model for predicting PDSI:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8,predictor9=predictor9,predictor10=predictor10)
      df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
      pcs = df.pca$x
      pc1 = pcs[,1]
      pc2 = pcs[,2]
      pc3 = pcs[,3]
      pc4 = pcs[,4]
      pc5 = pcs[,5]
      pc6 = pcs[,6]
      pc7 = pcs[,7]
      pc8 = pcs[,8]
      pc9 = pcs[,9]
      pc10 = pcs[,10]
      
      current_DF_pcs <- data.frame(Y=summer_PDSI,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8,predictor9=pc9,predictor10=pc10)
      mod_PDSI <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
      
      #PCA combo model for predicting SM-prctile:
      current_DF_pcs <- data.frame(Y=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8,predictor9=pc9,predictor10=pc10)
      mod_SMprctile <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
    }
    
    X = concurvity(mod_PDSI,full=FALSE)
    X=X$estimate
    idx <- which(X==1)
    X[idx] = 0
    Max_Concuvrity <- max(X)
    if (Max_Concuvrity <= concurvity_threshold){
      GCV = mod_PDSI$gcv.ubre
      store_data = c(n_predictors,i,GCV)
      store_GCV_PDSI = rbind(store_GCV_PDSI,store_data)
    }
    
    X = concurvity(mod_SMprctile,full=FALSE)
    X=X$estimate
    idx <- which(X==1)
    X[idx] = 0
    Max_Concuvrity <- max(X)
    if (Max_Concuvrity <= concurvity_threshold){
      GCV = mod_SMprctile$gcv.ubre
      store_data = c(n_predictors,i,GCV)
      store_GCV_SMprctile = rbind(store_GCV_SMprctile,store_data)
    }
  }
}

#choose 20 best ensemble members for PDSI predicting models
n_outmembers = 20 
PDSI_GCV = store_GCV_PDSI[,3]
sorted_pdsi_GCV = order(PDSI_GCV, decreasing = FALSE)
idx_best_pdsi <- sorted_pdsi_GCV[1:n_outmembers]
best_GCV_pdsi = PDSI_GCV[idx_best_pdsi]
best_npred_pdsi = store_GCV_PDSI[idx_best_pdsi,1]
best_comboID_pdsi = store_GCV_PDSI[idx_best_pdsi,2]

#choose 20 best ensemble members for SM-percentile predicting models
SMprctile_GCV = store_GCV_SMprctile[,3]
sorted_SM_GCV = order(SMprctile_GCV, decreasing = FALSE)
idx_best_SM <- sorted_SM_GCV[1:n_outmembers]
best_GCV_SM = SMprctile_GCV[idx_best_SM]
best_npred_SM = store_GCV_SMprctile[idx_best_SM,1]
best_comboID_SM = store_GCV_SMprctile[idx_best_SM,2]

#key outputs from part 1 - IDs defining top 20 models for PDSI and SM-percentile
Bestmods_pdsi <- data.frame(best_npred_pdsi=best_npred_pdsi,best_comboID_pdsi=best_comboID_pdsi,best_GCV_pdsi=best_GCV_pdsi)
Bestmods_sm <- data.frame(best_npred_SM=best_npred_SM,best_comboID_SM=best_comboID_SM,best_GCV_SM=best_GCV_SM)

##================================================================================================================================================
##================================================================================================================================================
##END PART 1
##================================================================================================================================================
##================================================================================================================================================

##================================================================================================================================================
##================================================================================================================================================
##PART 2 uses the best model selection from part 1 to make predictions on out-of-bag data (ie, cross-validate)
##================================================================================================================================================
##================================================================================================================================================


#loop through each of the best models, (1) define model formula first, then (2) make prediction with formula on out-of-bag data
nvars = npreds
x=1:nvars
#initialize outputs for predictions for summer drought
Store_drop1_SMprctile = c()
Store_drop1_PDSI = c()


for (ens in 1:n_outmembers){
  #(1) define model formula first
  print(c(c,ens))
  n_predictors_pdsi = Bestmods_pdsi[ens,1]
  comboID_pdsi = Bestmods_pdsi[ens,2]
  
  n_predictors_sm = Bestmods_sm[ens,1]
  comboID_sm = Bestmods_sm[ens,2]
  
  ncombos_pdsi = factorial(nvars)/(factorial(n_predictors_pdsi)*factorial(nvars-n_predictors_pdsi))
  combo_IDs_pdsi<-combn(x,n_predictors_pdsi)
  
  ncombos_sm = factorial(nvars)/(factorial(n_predictors_sm)*factorial(nvars-n_predictors_sm))
  combo_IDs_sm<-combn(x,n_predictors_sm)
  
  if (n_predictors_pdsi == 1){
    predictor1 <-Covariates_DF_all[,combo_IDs_pdsi[1,comboID_pdsi]]
    #PCA combo model for predicting PDSI:
    Covariates_DF = data.frame(predictor1=predictor1)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    current_DF_pcs_pdsi <- data.frame(Y1=summer_PDSI,Y2=summer_SMprctile,predictor1=pc1)
    mod_PDSI <- gam(data=current_DF_pcs_pdsi,Y1~s(predictor1),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_pdsi == 2){
    predictor1 <-Covariates_DF_all[,combo_IDs_pdsi[1,comboID_pdsi]]
    predictor2 <-Covariates_DF_all[,combo_IDs_pdsi[2,comboID_pdsi]]
    
    #PCA combo model for predicting PDSI:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    current_DF_pcs_pdsi <- data.frame(Y1=summer_PDSI,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2)
    mod_PDSI <- gam(data=current_DF_pcs_pdsi,Y1~s(predictor1)+s(predictor2),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_pdsi == 3){
    predictor1 <-Covariates_DF_all[,combo_IDs_pdsi[1,comboID_pdsi]]
    predictor2 <-Covariates_DF_all[,combo_IDs_pdsi[2,comboID_pdsi]]
    predictor3 <-Covariates_DF_all[,combo_IDs_pdsi[3,comboID_pdsi]]
    
    #PCA combo model for predicting PDSI:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    current_DF_pcs_pdsi <- data.frame(Y1=summer_PDSI,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3)
    mod_PDSI <- gam(data=current_DF_pcs_pdsi,Y1~s(predictor1)+s(predictor2)+s(predictor3),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_pdsi == 4){
    predictor1 <-Covariates_DF_all[,combo_IDs_pdsi[1,comboID_pdsi]]
    predictor2 <-Covariates_DF_all[,combo_IDs_pdsi[2,comboID_pdsi]]
    predictor3 <-Covariates_DF_all[,combo_IDs_pdsi[3,comboID_pdsi]]
    predictor4 <-Covariates_DF_all[,combo_IDs_pdsi[4,comboID_pdsi]]
    
    #PCA combo model for predicting PDSI:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    current_DF_pcs_pdsi <- data.frame(Y1=summer_PDSI,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4)
    mod_PDSI <- gam(data=current_DF_pcs_pdsi,Y1~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_pdsi == 5){
    predictor1 <-Covariates_DF_all[,combo_IDs_pdsi[1,comboID_pdsi]]
    predictor2 <-Covariates_DF_all[,combo_IDs_pdsi[2,comboID_pdsi]]
    predictor3 <-Covariates_DF_all[,combo_IDs_pdsi[3,comboID_pdsi]]
    predictor4 <-Covariates_DF_all[,combo_IDs_pdsi[4,comboID_pdsi]]
    predictor5 <-Covariates_DF_all[,combo_IDs_pdsi[5,comboID_pdsi]]
    
    #PCA combo model for predicting PDSI:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    pc5 = pcs[,5]
    current_DF_pcs_pdsi <- data.frame(Y1=summer_PDSI,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5)
    mod_PDSI <- gam(data=current_DF_pcs_pdsi,Y1~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_pdsi == 6){
    predictor1 <-Covariates_DF_all[,combo_IDs_pdsi[1,comboID_pdsi]]
    predictor2 <-Covariates_DF_all[,combo_IDs_pdsi[2,comboID_pdsi]]
    predictor3 <-Covariates_DF_all[,combo_IDs_pdsi[3,comboID_pdsi]]
    predictor4 <-Covariates_DF_all[,combo_IDs_pdsi[4,comboID_pdsi]]
    predictor5 <-Covariates_DF_all[,combo_IDs_pdsi[5,comboID_pdsi]]
    predictor6 <-Covariates_DF_all[,combo_IDs_pdsi[6,comboID_pdsi]]
    
    #PCA combo model for predicting PDSI:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    pc5 = pcs[,5]
    pc6 = pcs[,6]
    current_DF_pcs_pdsi <- data.frame(Y1=summer_PDSI,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6)
    mod_PDSI <- gam(data=current_DF_pcs_pdsi,Y1~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_pdsi == 7){
    predictor1 <-Covariates_DF_all[,combo_IDs_pdsi[1,comboID_pdsi]]
    predictor2 <-Covariates_DF_all[,combo_IDs_pdsi[2,comboID_pdsi]]
    predictor3 <-Covariates_DF_all[,combo_IDs_pdsi[3,comboID_pdsi]]
    predictor4 <-Covariates_DF_all[,combo_IDs_pdsi[4,comboID_pdsi]]
    predictor5 <-Covariates_DF_all[,combo_IDs_pdsi[5,comboID_pdsi]]
    predictor6 <-Covariates_DF_all[,combo_IDs_pdsi[6,comboID_pdsi]]
    predictor7 <-Covariates_DF_all[,combo_IDs_pdsi[7,comboID_pdsi]]
    
    #PCA combo model for predicting PDSI:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    pc5 = pcs[,5]
    pc6 = pcs[,6]
    pc7 = pcs[,7]
    current_DF_pcs_pdsi <- data.frame(Y1=summer_PDSI,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7)
    mod_PDSI <- gam(data=current_DF_pcs_pdsi,Y1~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_pdsi == 8){
    predictor1 <-Covariates_DF_all[,combo_IDs_pdsi[1,comboID_pdsi]]
    predictor2 <-Covariates_DF_all[,combo_IDs_pdsi[2,comboID_pdsi]]
    predictor3 <-Covariates_DF_all[,combo_IDs_pdsi[3,comboID_pdsi]]
    predictor4 <-Covariates_DF_all[,combo_IDs_pdsi[4,comboID_pdsi]]
    predictor5 <-Covariates_DF_all[,combo_IDs_pdsi[5,comboID_pdsi]]
    predictor6 <-Covariates_DF_all[,combo_IDs_pdsi[6,comboID_pdsi]]
    predictor7 <-Covariates_DF_all[,combo_IDs_pdsi[7,comboID_pdsi]]
    predictor8 <-Covariates_DF_all[,combo_IDs_pdsi[8,comboID_pdsi]]
    
    #PCA combo model for predicting PDSI:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    pc5 = pcs[,5]
    pc6 = pcs[,6]
    pc7 = pcs[,7]
    pc8 = pcs[,8]
    current_DF_pcs_pdsi <- data.frame(Y1=summer_PDSI,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8)
    mod_PDSI <- gam(data=current_DF_pcs_pdsi,Y1~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_pdsi == 9){
    predictor1 <-Covariates_DF_all[,combo_IDs_pdsi[1,comboID_pdsi]]
    predictor2 <-Covariates_DF_all[,combo_IDs_pdsi[2,comboID_pdsi]]
    predictor3 <-Covariates_DF_all[,combo_IDs_pdsi[3,comboID_pdsi]]
    predictor4 <-Covariates_DF_all[,combo_IDs_pdsi[4,comboID_pdsi]]
    predictor5 <-Covariates_DF_all[,combo_IDs_pdsi[5,comboID_pdsi]]
    predictor6 <-Covariates_DF_all[,combo_IDs_pdsi[6,comboID_pdsi]]
    predictor7 <-Covariates_DF_all[,combo_IDs_pdsi[7,comboID_pdsi]]
    predictor8 <-Covariates_DF_all[,combo_IDs_pdsi[8,comboID_pdsi]]
    predictor9 <-Covariates_DF_all[,combo_IDs_pdsi[9,comboID_pdsi]]
    
    #PCA combo model for predicting PDSI:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8,predictor9=predictor9)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    pc5 = pcs[,5]
    pc6 = pcs[,6]
    pc7 = pcs[,7]
    pc8 = pcs[,8]
    pc9 = pcs[,9]
    
    current_DF_pcs_pdsi <- data.frame(Y1=summer_PDSI,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8,predictor9=pc9)
    mod_PDSI <- gam(data=current_DF_pcs_pdsi,Y1~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_pdsi == 10){
    predictor1 <-Covariates_DF_all[,combo_IDs_pdsi[1,comboID_pdsi]]
    predictor2 <-Covariates_DF_all[,combo_IDs_pdsi[2,comboID_pdsi]]
    predictor3 <-Covariates_DF_all[,combo_IDs_pdsi[3,comboID_pdsi]]
    predictor4 <-Covariates_DF_all[,combo_IDs_pdsi[4,comboID_pdsi]]
    predictor5 <-Covariates_DF_all[,combo_IDs_pdsi[5,comboID_pdsi]]
    predictor6 <-Covariates_DF_all[,combo_IDs_pdsi[6,comboID_pdsi]]
    predictor7 <-Covariates_DF_all[,combo_IDs_pdsi[7,comboID_pdsi]]
    predictor8 <-Covariates_DF_all[,combo_IDs_pdsi[8,comboID_pdsi]]
    predictor9 <-Covariates_DF_all[,combo_IDs_pdsi[9,comboID_pdsi]]
    predictor10 <-Covariates_DF_all[,combo_IDs_pdsi[10,comboID_pdsi]]
    #PCA combo model for predicting PDSI:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8,predictor9=predictor9,predictor10=predictor10)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    pc5 = pcs[,5]
    pc6 = pcs[,6]
    pc7 = pcs[,7]
    pc8 = pcs[,8]
    pc9 = pcs[,9]
    pc10 = pcs[,10]
    current_DF_pcs_pdsi <- data.frame(Y1=summer_PDSI,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8,predictor9=pc9,predictor10=pc10)
    mod_PDSI <- gam(data=current_DF_pcs_pdsi,Y1~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_sm == 1){
    #PCA combo model for predicting SM-prctile:
    predictor1 <-Covariates_DF_all[,combo_IDs_sm[1,comboID_sm]]
    #PCA combo model for predicting SMprctile:
    Covariates_DF = data.frame(predictor1=predictor1)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    current_DF_pcs_sm <- data.frame(Y1=summer_SMprctile,Y2=summer_SMprctile,predictor1=pc1)
    mod_SMprctile <- gam(data=current_DF_pcs_sm,Y2~s(predictor1),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_sm == 2){
    #PCA combo model for predicting SM-prctile:
    predictor1 <-Covariates_DF_all[,combo_IDs_sm[1,comboID_sm]]
    predictor2 <-Covariates_DF_all[,combo_IDs_sm[2,comboID_sm]]
    
    #PCA combo model for predicting SMprctile:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    current_DF_pcs_sm <- data.frame(Y1=summer_SMprctile,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2)
    mod_SMprctile <- gam(data=current_DF_pcs_sm,Y2~s(predictor1)+s(predictor2),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_sm == 3){
    #PCA combo model for predicting SM-prctile:
    predictor1 <-Covariates_DF_all[,combo_IDs_sm[1,comboID_sm]]
    predictor2 <-Covariates_DF_all[,combo_IDs_sm[2,comboID_sm]]
    predictor3 <-Covariates_DF_all[,combo_IDs_sm[3,comboID_sm]]
    
    #PCA combo model for predicting SMprctile:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    current_DF_pcs_sm <- data.frame(Y1=summer_SMprctile,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3)
    mod_SMprctile <- gam(data=current_DF_pcs_sm,Y2~s(predictor1)+s(predictor2)+s(predictor3),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_sm == 4){
    #PCA combo model for predicting SM-prctile:
    predictor1 <-Covariates_DF_all[,combo_IDs_sm[1,comboID_sm]]
    predictor2 <-Covariates_DF_all[,combo_IDs_sm[2,comboID_sm]]
    predictor3 <-Covariates_DF_all[,combo_IDs_sm[3,comboID_sm]]
    predictor4 <-Covariates_DF_all[,combo_IDs_sm[4,comboID_sm]]
    
    #PCA combo model for predicting SMprctile:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    current_DF_pcs_sm <- data.frame(Y1=summer_SMprctile,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4)
    mod_SMprctile <- gam(data=current_DF_pcs_sm,Y2~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_sm == 5){
    #PCA combo model for predicting SM-prctile:
    predictor1 <-Covariates_DF_all[,combo_IDs_sm[1,comboID_sm]]
    predictor2 <-Covariates_DF_all[,combo_IDs_sm[2,comboID_sm]]
    predictor3 <-Covariates_DF_all[,combo_IDs_sm[3,comboID_sm]]
    predictor4 <-Covariates_DF_all[,combo_IDs_sm[4,comboID_sm]]
    predictor5 <-Covariates_DF_all[,combo_IDs_sm[5,comboID_sm]]
    
    #PCA combo model for predicting SMprctile:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    pc5 = pcs[,5]
    current_DF_pcs_sm <- data.frame(Y1=summer_SMprctile,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5)
    mod_SMprctile <- gam(data=current_DF_pcs_sm,Y2~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_sm == 6){
    #PCA combo model for predicting SM-prctile:
    predictor1 <-Covariates_DF_all[,combo_IDs_sm[1,comboID_sm]]
    predictor2 <-Covariates_DF_all[,combo_IDs_sm[2,comboID_sm]]
    predictor3 <-Covariates_DF_all[,combo_IDs_sm[3,comboID_sm]]
    predictor4 <-Covariates_DF_all[,combo_IDs_sm[4,comboID_sm]]
    predictor5 <-Covariates_DF_all[,combo_IDs_sm[5,comboID_sm]]
    predictor6 <-Covariates_DF_all[,combo_IDs_sm[6,comboID_sm]]
    
    #PCA combo model for predicting SMprctile:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    pc5 = pcs[,5]
    pc6 = pcs[,6]
    current_DF_pcs_sm <- data.frame(Y1=summer_SMprctile,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6)
    mod_SMprctile <- gam(data=current_DF_pcs_sm,Y2~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_sm == 7){
    #PCA combo model for predicting SM-prctile:
    predictor1 <-Covariates_DF_all[,combo_IDs_sm[1,comboID_sm]]
    predictor2 <-Covariates_DF_all[,combo_IDs_sm[2,comboID_sm]]
    predictor3 <-Covariates_DF_all[,combo_IDs_sm[3,comboID_sm]]
    predictor4 <-Covariates_DF_all[,combo_IDs_sm[4,comboID_sm]]
    predictor5 <-Covariates_DF_all[,combo_IDs_sm[5,comboID_sm]]
    predictor6 <-Covariates_DF_all[,combo_IDs_sm[6,comboID_sm]]
    predictor7 <-Covariates_DF_all[,combo_IDs_sm[7,comboID_sm]]
    
    #PCA combo model for predicting SMprctile:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    pc5 = pcs[,5]
    pc6 = pcs[,6]
    pc7 = pcs[,7]
    current_DF_pcs_sm <- data.frame(Y1=summer_SMprctile,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7)
    mod_SMprctile <- gam(data=current_DF_pcs_sm,Y2~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_sm == 8){
    #PCA combo model for predicting SM-prctile:
    predictor1 <-Covariates_DF_all[,combo_IDs_sm[1,comboID_sm]]
    predictor2 <-Covariates_DF_all[,combo_IDs_sm[2,comboID_sm]]
    predictor3 <-Covariates_DF_all[,combo_IDs_sm[3,comboID_sm]]
    predictor4 <-Covariates_DF_all[,combo_IDs_sm[4,comboID_sm]]
    predictor5 <-Covariates_DF_all[,combo_IDs_sm[5,comboID_sm]]
    predictor6 <-Covariates_DF_all[,combo_IDs_sm[6,comboID_sm]]
    predictor7 <-Covariates_DF_all[,combo_IDs_sm[7,comboID_sm]]
    predictor8 <-Covariates_DF_all[,combo_IDs_sm[8,comboID_sm]]
    
    #PCA combo model for predicting SMprctile:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    pc5 = pcs[,5]
    pc6 = pcs[,6]
    pc7 = pcs[,7]
    pc8 = pcs[,8]
    current_DF_pcs_sm <- data.frame(Y1=summer_SMprctile,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8)
    mod_SMprctile <- gam(data=current_DF_pcs_sm,Y2~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_sm == 9){
    #PCA combo model for predicting SM-prctile:
    predictor1 <-Covariates_DF_all[,combo_IDs_sm[1,comboID_sm]]
    predictor2 <-Covariates_DF_all[,combo_IDs_sm[2,comboID_sm]]
    predictor3 <-Covariates_DF_all[,combo_IDs_sm[3,comboID_sm]]
    predictor4 <-Covariates_DF_all[,combo_IDs_sm[4,comboID_sm]]
    predictor5 <-Covariates_DF_all[,combo_IDs_sm[5,comboID_sm]]
    predictor6 <-Covariates_DF_all[,combo_IDs_sm[6,comboID_sm]]
    predictor7 <-Covariates_DF_all[,combo_IDs_sm[7,comboID_sm]]
    predictor8 <-Covariates_DF_all[,combo_IDs_sm[8,comboID_sm]]
    predictor9 <-Covariates_DF_all[,combo_IDs_sm[9,comboID_sm]]
    
    #PCA combo model for predicting SMprctile:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8,predictor9=predictor9)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    pc5 = pcs[,5]
    pc6 = pcs[,6]
    pc7 = pcs[,7]
    pc8 = pcs[,8]
    pc9 = pcs[,9]
    
    current_DF_pcs_sm <- data.frame(Y1=summer_SMprctile,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8,predictor9=pc9)
    mod_SMprctile <- gam(data=current_DF_pcs_sm,Y2~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  if (n_predictors_sm == 10){
    #PCA combo model for predicting SM-prctile:
    predictor1 <-Covariates_DF_all[,combo_IDs_sm[1,comboID_sm]]
    predictor2 <-Covariates_DF_all[,combo_IDs_sm[2,comboID_sm]]
    predictor3 <-Covariates_DF_all[,combo_IDs_sm[3,comboID_sm]]
    predictor4 <-Covariates_DF_all[,combo_IDs_sm[4,comboID_sm]]
    predictor5 <-Covariates_DF_all[,combo_IDs_sm[5,comboID_sm]]
    predictor6 <-Covariates_DF_all[,combo_IDs_sm[6,comboID_sm]]
    predictor7 <-Covariates_DF_all[,combo_IDs_sm[7,comboID_sm]]
    predictor8 <-Covariates_DF_all[,combo_IDs_sm[8,comboID_sm]]
    predictor9 <-Covariates_DF_all[,combo_IDs_sm[9,comboID_sm]]
    predictor10 <-Covariates_DF_all[,combo_IDs_sm[10,comboID_sm]]
    
    #PCA combo model for predicting SMprctile:
    Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8,predictor9=predictor9,predictor10=predictor10)
    df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
    pcs = df.pca$x
    pc1 = pcs[,1]
    pc2 = pcs[,2]
    pc3 = pcs[,3]
    pc4 = pcs[,4]
    pc5 = pcs[,5]
    pc6 = pcs[,6]
    pc7 = pcs[,7]
    pc8 = pcs[,8]
    pc9 = pcs[,9]
    pc10 = pcs[,10]
    current_DF_pcs_sm <- data.frame(Y1=summer_SMprctile,Y2=summer_SMprctile,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8,predictor9=pc9,predictor10=pc10)
    mod_SMprctile <- gam(data=current_DF_pcs_sm,Y2~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4),family="gaussian",control=list(mgcv.tol=1e-4))
  }
  
  #(2) make prediction with formula on out-of-bag data
  WYs = years #1982-2020 (total analysis period)
  index = 1:length(WYs)
  yhat_SMprctile=1
  yhat_PDSI=1
  for(Y in WYs[1]:WYs[length(WYs)]){
    #drop 1 year
    IDX_drop <-which(WYs == Y)
    index1=index[-IDX_drop]
    
    #PCA combo model:
    dropped_DF_pcs_pdsi <- current_DF_pcs_pdsi[index1,]
    dropped_DF_pcs_sm <- current_DF_pcs_sm[index1,]
    
    #model PDSI (fit on all data minus the dropped year)
    fit_drop_gam_PDSI = gam(mod_PDSI$formula,data=dropped_DF_pcs_pdsi,family="gaussian",control=list(mgcv.tol=1e-4))
    #model SM 
    fit_drop_gam_SMprctile = gam(mod_SMprctile$formula,data=dropped_DF_pcs_sm,family="gaussian",control=list(mgcv.tol=1e-4))
    
    #define the point that was dropped
    newdata_pdsi = current_DF_pcs_pdsi[IDX_drop,] 
    newdata_sm = current_DF_pcs_sm[IDX_drop,] 
    
    #predict on the dropped point using the above fit-model
    
    #PDSI prediction
    yhat <- predict(fit_drop_gam_PDSI,newdata=newdata_pdsi,se.fit = TRUE)
    yhat_PDSI[IDX_drop] =yhat$fit
    

    #SM-percentile prediction
    yhat <- predict(fit_drop_gam_SMprctile,newdata=newdata_sm,se.fit = TRUE)
    yhat_SMprctile[IDX_drop] =yhat$fit
  }
  
  Store_drop1_SMprctile = cbind(Store_drop1_SMprctile,as.matrix(yhat_SMprctile))
  Store_drop1_PDSI = cbind(Store_drop1_PDSI,as.matrix(yhat_PDSI))
}

#get the ensemble mean predictions
PDSI_ens_mean = rowMeans(Store_drop1_PDSI)
SMprctile_ens_mean = rowMeans(Store_drop1_SMprctile)

#you can bring SM-percentile prediction back to percentile space (0-1)
#using eq. 2 from: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2022WR033734
sorted_idx = order(SMprctile_ens_mean)
sorted_SMprctile_ens_mean = SMprctile_ens_mean[sorted_idx]
rank=1:length(sorted_idx)
rank[sorted_idx] = rank
SMprctile_ens_mean_rescale = (rank - min(rank))/( max(rank) - min(rank))
SMprctile_ens_mean_rescale = rescale(SMprctile_ens_mean_rescale,to=c(0,1))

#similarly, PDSI can be redefined to be applied to new constraints (eg historically observed constraints)

#plot results and calculate summary stats:
Ysm = summer_SMprctile
Ypdsi = summer_PDSI

plot(pre_summer_PDSI,Ypdsi,xlab="pre-summer PDSI",ylab="monitored summer PDSI")
plot(PDSI_ens_mean,Ypdsi,xlab="predicted summer PDSI",ylab="monitored summer PDSI")

plot(pre_summer_SMprctile,Ysm,xlab="pre-summer soil moisture percentile",ylab="monitored summer soil moisture percentile")
plot(SMprctile_ens_mean_rescale,Ysm,xlab="predicted summer soil moisture percentile",ylab="monitored summer soil moisture percentile")

#report R between predicted and monitored conditions:
R_pdsi = cor(PDSI_ens_mean,Ypdsi)
R_sm = cor(SMprctile_ens_mean_rescale,Ysm)
#a great deal of the predictive information comes from the auto correlation of drought
#calculate the correlation between the pre-summer drought index and summer drought index to consider added value of the models (vs just using autocorrelation):
auotR_pdsi = cor(pre_summer_PDSI,Ypdsi)
autoR_sm = cor(pre_summer_SMprctile,Ysm)

sprintf("R for predicted PDSI = %.4f",R_pdsi)
sprintf("R for pre-summer PDSI = %.4f",auotR_pdsi)

sprintf("R for predicted SM = %.4f",R_sm)
sprintf("R for pre-summer SM = %.4f",autoR_sm)

#these results can be exported to be processed in a different platform (eg, python, excel, matlab, etc):
##e.g., with write.table()

