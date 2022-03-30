# POX-C calculations and Analysis

## Emma Link, last updated 3/30/2022
## A script for merging plate reads with metadata and calculating POX-C values. 
## You will need your raw data, your plate maps, and your weight/notes metadata. I reccommend adding treatment and blocks to your metadata.

## Before starting here, use ExcelPlateRead.R to get raw excel files into plater format and read into R. 
## Also use ExcelPlateRead.R to get plate key metadata into plater format and read into R. 

## See git repo for example data and metadata files that are ready to be run in this script.

## What we're going to do: 
### 1) Merge plate reads, PlateMetadata, and sample Metadata
### 2) Calculate standard curve slope and intercept for each run 
### 3) Calculate amount of POX-C
### 4) Visualizations and statistical analysis 

library(stringr)
library(tidyverse)
library(broom)

  
## 1) Merge plate reads, PlateMetadata, and sample Metadata (weights and notes from assay)

### First, we're going to save all of your datasets under a different name so you don't have to change names throughout this script.
### Currently is set up to run with the example datasets. You'll have to change the string to set your wd 
setwd("/.SEECRS-POX-C")
RunData <- read.csv("./Data/ExampleRawData")
#For some reason the example file is reading in a random first column, just delete it, but you shouldn't do this with your data probably 
RunData <- RunData[,-c(1)]
RunMap <- read.csv("./Data/ExamplePlateMap")
RunMap <- RunMap[,-c(1)]

#RunMedata <- *insert your sample metadata here*

### First we need to pull out the batch number from the plate name of both plate and metadata files
### depending on how you've named your files, you may need to pull a different string 

RunData$Batch <- str_sub(RunData$Plate, 17, 17)

RunMap_longer <- pivot_longer(RunMap, cols = starts_with("Batch"), names_to = "BatchName")
RunMap <- na.omit(RunMap_longer)
RunMap$Batch <- str_sub(RunMap$BatchName, 7,7)

### Now, we merge plate data and plate metadata. Here, merging by Batch (1 batch in each plate) and well id (ex. A5)

RunData_merged <- merge(RunData, RunMap, by = c("Batch", "Wells"))
RunData_merged <- subset(RunData_merged, select = -c(Plate.x, Plate.y, BatchName))
RunData_merged <- rename(RunData_merged, PlateRead = ...2 )
RunData_merged <- rename(RunData_merged, Plot = value)

### Now, we take depth out of the plot number string and add to it's own column 

RunData_merged$Depth <- str_sub(RunData_merged$Plot, 4,6)

for (i in 1:nrow(RunData_merged)){
  if(RunData_merged[i, "Depth"] == ""){RunData_merged[i, "Depth"] = "0-15"}
  else if(RunData_merged[i, "Depth"] == "-0" || RunData_merged[i, "Depth"] == ".0"){RunData_merged[i, "Depth"] = "0-15"}
  else if(RunData_merged[i, "Depth"] == "-15"){RunData_merged[i, "Depth"] = "15-30"}
  else if(RunData_merged[i, "Depth"] == "-30"){RunData_merged[i, "Depth"] = "15-30"}
}
### Now, we merge in the weight metadata. Recommend also putting block and treatment in with weight 
RunData_merged <- merge(RunData_merged, RunMetadata, by = c("Batch", "Plot"))


## 2) Calculate standard curve slope and intercept for each run 

### Calculate averages, sd, and cv for all plots/standard wells. Search for cv > 10
RunData_merged <- RunData_merged %>%
  group_by(Batch, Plot) %>%
  mutate(Plot_mean = mean(PlateRead), sd = sd(PlateRead), cv = sd(PlateRead)/mean(PlateRead)*100)

### Assign concentrations 

for (i in 1:nrow(RunData_merged)) {
  if(RunData_merged[i, "Plot"] == "SLTN_BLNK") {RunData_merged[i, "Concentration"] = "0"}
  else if (RunData_merged[i, "Plot"] == "0.05mM") {RunData_merged[i, "Concentration"] = "0.005"}
  else if (RunData_merged[i, "Plot"] == "0.15mM") {RunData_merged[i, "Concentration"] = "0.015"}
  else if (RunData_merged[i, "Plot"] == "0.2mM") {RunData_merged[i, "Concentration"] = "0.02"}
  else if (RunData_merged[i, "Plot"] == "0.1mM") {RunData_merged[i, "Concentration"] = "0.01"}
  else {RunData_merged[i, "Concentration"] = "0.02"}
}
RunData_merged$Concentration <- as.numeric(RunData_merged$Concentration)

### Calculate curve slope and intercept in another dataframe
####Y axis is molarity, x axis is absorbance values. See written protocol
StandardCurves <- RunData_merged %>%
  group_by(Batch) %>%
  filter(Plot %in% c("0.05mM", "SLTN_BLNK", "0.15mM", "0.1mM", "0.2mM")) %>%
  mutate(Intercept = tidy(lm(Concentration ~ Plot_mean))[1,2], Slope = tidy(lm(Concentration ~ Plot_mean))[2,2])

### Remove duplicate entries from standard curves and keep only batch, intercept, and slope columns as integer (not list)
StandardCurves <- StandardCurves[,c(1,10,11)]
StandardCurves$Slope <- StandardCurves$Slope[1]
StandardCurves$Intercept <- StandardCurves$Intercept[1]
StandardCurves <- distinct(StandardCurves)


### Remove from general dataset entries that are not sample reads 
### First, define function that is the opposite of %in% 
'%notin%' = Negate('%in%')
RunData_merged <- RunData_merged[RunData_merged$Plot %notin% c("0.05mM", "SLTN_BLNK", "0.15mM", "0.1mM", "0.2mM"),]

###Merge standard curves in with the general reads 
RunData_merged <- merge(RunData_merged, StandardCurves, by = "Batch")


## 3) Calculate amount of POX-C

### POXC (mg kg-1 soil) = [0.02 M/ L - (a + b × Abs)] × (9000 mg C/M) × (0.02 L/ 0.0025 kg soil), see https://docs.google.com/document/d/16QBG4-aZM86iIVg-GIjTXB-V69pyZ12Qzvw4IdvYATU/edit
### Note that currently weight of soil is fixed; we need to change it so that weight is using the recorded weight.                                   

RunData_merged <- RunData_merged %>%
  mutate(POXC = ((0.02-(Intercept + (Slope * Plot_mean)))*9000*(0.02/0.0025))) 

### Keep only the distinct entries (get rid of the triplicates since we calculate POX on the averages anyway)
RunData_final <- distinct(RunData_merged[,-c(2,3)])

view(RunData_final)

### Do any cleanup of the dataset that you need
## for examples, my Plot names need cleanup; I need to pull the depth indicators out of the names
###I'm going to do this by truncating the length of the plot name to 3 
RunData_final$Plot <- substr(RunData_final$Plot, 1,3)

### Save a csv copy of your dataset at this point - you're done with calculations! 
### make sure you update the location and name you want to write to
write.csv(RunData_final, "./Data/ExampleFinalDataset_2021")

## 4) Visualizations and statistical analysis
### This is Emma's pipeline, doesn't necessarily have to be your pipeline 
### Currently set up for PARCE which is a blocked experiment, hence block effect investigation 
