#load required library
library(xcms)
#set number of cores
noOfCores <- 4

## Use socket based parallel processing on Windows systems
if (.Platform$OS.type == "unix") {
  register(bpstart(MulticoreParam(noOfCores)))
} else {
  register(bpstart(SnowParam(noOfCores)))
}

### change the working directory to the samples directory

old_dir <- setwd("/home/sangeeta/Downloads/analysis_xcms")


library(plyr)


files <- list.dirs(".")
print(files)

for (i in files[-1]) {
   print(i)
}


for (i in files[-1]) {
     mzMLfiles <- list.files(i)
     s_groups <- rep(i, length(mzMLfiles))


    pheno <- data.frame(
    sample_name= mzMLfiles,
    sample_group= sub("./", "", s_groups),  stringsAsFactors = FALSE)

    full_df <- do.call(rbind, pheno)
    write.csv(full_df, "first.csv")
}



do_analysis <- function(folder_name){

  mzxml <- list.files(folder_name)
  usefiles <- list.files(folder_name, full.names= TRUE)
  s_groups <- rep(sub("./", "", folder_name), length(mzMLfiles))
  pd <- data.frame(sample_name= mzxml,
    sample_group = s_groups, stringsAsFactors = FALSE)
  rawData <- readMSData(usefiles, centroided. = TRUE, mode = "onDisk", pdata= new("NAnnotatedDataFrame", pd))

  bpis <- chromatogram(rawData, aggregationFun = "max")
  plot(bpis)
  cwp <- CentWaveParam(snthresh = 5, noise = 1000, peakwidth = c(3, 15), ppm = 10)
  processedData <- findChromPeaks(rawData, param = cwp)
  plotChromPeakImage(processedData, binSize = 10)
  processedData <- adjustRtime(processedData, param = ObiwarpParam())
  plotAdjustedRtime(processedData)
  pdp <- PeakDensityParam(sampleGroups = processedData$sample_group,
                        minFraction = 0.10)
  processedData <- groupChromPeaks(processedData, param = pdp)
  featuresDef <- featureDefinitions(processedData)
  featuresIntensities <- featureValues(processedData, value = "into")
  dataTable <- merge(featuresDef, featuresIntensities, by = 0, all = TRUE)
  dataTable <- dataTable[, !(colnames(dataTable) %in% c("peakidx"))]
  write.table(dataTable, paste(s_groups[[1]], ".csv"), sep = ",", quote = FALSE,
            row.names = FALSE)

  }


  main <- function(){
      for (i in files[-1]) {
              do_analysis(i)
}
}


## Running the main function to start the analysis

if (!interactive()) {
  main()
}
