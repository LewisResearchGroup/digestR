################################################################################
################################################################################
##                                                                            ##
##                                                                            ##
##  digestR version 1.0.0, Tools for viewing and anlyzing protein catabolism.    ##
##    Copyright (C) 2023, Dimitri Desmonts de Lamache, Raied Aburashed,       ##
##          Travis A. Bingemann, Ian A. Lewis under GPL-3                     ##
##                                                                            ##
##    This program is free software: you can redistribute it and/or modify    ##
##    it under the terms of the GNU General Public License as published by    ##
##    the Free Software Foundation, either version 3 of the License, or       ##
##    any later version.                                                      ##
##                                                                            ##
##    This program is distributed in the hope that it will be useful,         ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of          ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ## 
##    GNU General Public License for more details.                            ##
##                                                                            ##
##    A copy of the GNU General Public License can be found at:               ##
##    www.r-project.org/Licenses/GPL-3                                        ##
##                                                                            ##
##                                                                            ##
################################################################################
################################################################################

library(tcltk)
library(Rcpp)
library(ggplot2)
library(ggridges)
library(dplyr)


################################################################################
#
#Digestomics MK4
#
################################################################################

extractPathAndName <- function(sPath)
{
  
  if ( is.null(sPath) || nchar(sPath) == 0 || is.na(sPath) )
  {
    return('')
  }
  
  vIndices <- gregexpr('\\', sPath, fixed = TRUE)
  vEndices <- gregexpr('.', sPath, fixed = TRUE)
  if(vIndices[[1]][1] == -1)
  {
    vIndices <- gregexpr('/', sPath, fixed = TRUE)
  }
  
  if(vIndices[[1]][1] == -1)
  {
    indexStart <- 1
    path <- ''
  }else
  {
    indexStart <- vIndices[[1]][length(vIndices[[1]])] + 1
    path <- substr(sPath, 1, vIndices[[1]][length(vIndices[[1]])])
  }
  
  sName <- substr(sPath, indexStart, vEndices[[1]][1] - 1)
  vIndices <- gregexpr('_', sName, fixed = TRUE)
  
  if(vIndices[[1]][1] == -1)
  {
    indexStart <- 1
  }else
  {
    indexStart <- vIndices[[1]][1]
  }	
  sName <- substr(sName, indexStart, nchar(sName))
  
  vIndices <- gregexpr('[[:alnum:]]', sName, fixed = FALSE)
  
  if(vIndices[[1]][1] == -1)
  {
    indexStart <- 1
  }else
  {
    indexStart <- vIndices[[1]][1]
  }		
  
  result <- data.frame(path = path, name = substr(sName, indexStart, nchar(sName)), stringsAsFactors = FALSE)
  
  return(result)
}

extractFileNameNoExt <- function(sPath)
{
  
  if ( is.null(sPath) || nchar(sPath) == 0 || is.na(sPath) )
  {
    return('')
  }
  
  vIndices <- gregexpr('\\', sPath, fixed = TRUE)
  vEndices <- gregexpr('.', sPath, fixed = TRUE)
  if(vIndices[[1]][1] == -1)
  {
    vIndices <- gregexpr('/', sPath, fixed = TRUE)
  }
  
  if(vIndices[[1]][1] == -1)
  {
    indexStart <- 1
  }else
  {
    indexStart <- vIndices[[1]][length(vIndices[[1]])] + 1
  }
  
  return(substr(sPath, indexStart, vEndices[[1]][length(vEndices[[1]])] - 1))
}

extractFileName <- function(sPath)
{
  
  if ( is.null(sPath) || nchar(sPath) == 0 || is.na(sPath) )
  {
    return('')
  }
  
  vIndices <- gregexpr('\\', sPath, fixed = TRUE)
  if(vIndices[[1]][1] == -1)
  {
    vIndices <- gregexpr('/', sPath, fixed = TRUE)
  }
  
  if(vIndices[[1]][1] == -1)
  {
    indexStart <- 1
  }else
  {
    indexStart <- vIndices[[1]][length(vIndices[[1]])] + 1
  }
  
  return(substr(sPath, indexStart, nchar(sPath)))
}

copyFolderNamesToNewLocation<- function()
{
  srcDir <- tkchooseDirectory()
  destDir <- tkchooseDirectory()
  destDir <- paste(destDir, '/', sep='')
  
  folders <- list.dirs(path = file.path(srcDir), recursive = FALSE)
  
  for (i in 1:length(folders))
  {
    tmp <- gregexpr('/', folders[i])[[1]]
    sName <- substr(folders[i], tmp[length(tmp)]+1, nchar(folders[i]))
    dir.create(paste(destDir, sName, sep=''))
  }
}

fixFiles<- function(lSpecies)
{
  dir <- tkchooseDirectory()
  
  subDirs <- list.dirs(path = file.path(dir), recursive = FALSE)
  
  for (i in 1:length(subDirs))
  {
    
    sPath <- paste(subDirs[i], "/", "DIANA", sep="", collapse="")
    lFilesDCF <- list.files(path = sPath, pattern = '*.dcf', all.files=TRUE, full.names = FALSE, recursive = FALSE, ignore.case = TRUE)
    lDCFNames <- vector(mode="character", length = length(lFilesDCF))
    
    if (length(lFilesDCF) > 0)
    {
      for (j in 1:length(lFilesDCF))
      {
        lDCFNames[j] <- extractFileNameNoExt(lFilesDCF[j])
        result <- readDIANA(paste(sPath, '/', lFilesDCF[j], sep='', collapse=''), lSpecies)
        writeDIANA(paste(sPath, '/', lDCFNames[j], '_fixed.dcf', sep='', collapse=''), lSpecies, result)
      }			
    }
  }
}

renameCSVFiles<- function(myPath)
{
  
  subDirs <- list.dirs(path = file.path(myPath), full.names = TRUE, recursive = FALSE)
  
  for (i in 1:length(subDirs))
  {
    lFilesCSV <- list.files(path = subDirs[i], pattern = '*.csv', all.files=TRUE, full.names = FALSE, recursive = FALSE, ignore.case = TRUE)
    
    if (length(lFilesCSV) > 0)
    {
      for (j in 1:length(lFilesCSV))
      {
        lNewName <- strsplit(lFilesCSV[j], "___", fixed = TRUE)[[1]][2]
        lNewName <- strsplit(lNewName, ".", fixed = TRUE)[[1]][1]
        lNewName <- paste(lNewName, "csv", sep=".")
        lNewName <- paste(subDirs[i], '/', lNewName, sep="")
        file.rename(paste(subDirs[i], '/', lFilesCSV[j], sep=''), lNewName)
      }			
    }
  }
}

drawAndLabelAxisGenes <- function(in.folder, xVals, xDivs, labPos, xLabels, xlty, xAxis)
{
  xBounds <- vector(mode= "numeric", length = length(xDivs) + 1)
  xBounds <- c(1, xDivs)
  xAxisLabels <- unlist(strsplit(xAxis, "", fixed = TRUE))
  op <- par('mgp')
  par(mgp=c(1,0,0))
  
  axis(side=1, at= xVals, labels = xAxisLabels, line=0,
       lty=xlty, tck=-0.01, 
       cex.axis = in.folder$graphics.par$cex.axis * 0.5)
  
  axis(side=1, at= labPos, labels = xLabels, line=1, lty=0, tck=0)
  axis(side=1, at= xBounds, labels = FALSE, line=1, col=NA, col.ticks=in.folder$graphics.par$col.axis, lty=1, tck=in.folder$graphics.par$xtck)
  
  par(mgp=op)
}

drawAndLabelAxis <- function(in.folder, xVals, labPos, xLabels, xlty)
{
  xBounds <- vector(mode= "numeric", length = length(xVals) + 1)
  xBounds <- c(1, xVals)
  
  axis(side=1, at= xBounds, labels = FALSE, 
       lty=xlty, tck=in.folder$graphics.par$xtck, 
       cex.axis=in.folder$graphics.par$cex.axis)
  
  axis(side=1, at= labPos, labels = xLabels, 
       lty=xlty, tck=0)	
}

drawAndLabelGenes <- function(in.folder, axisRange, span, labels, xlty)
{
  first <- axisRange[1]
  last <- axisRange[2]
  print('span')
  print(span)
  
  xBounds <- c(span$x1[1], span$x2)
  
  if(xBounds[1] < first)
    xBounds[1] <- first
  
  if(xBounds[length(xBounds)] > last)
    xBounds[length(xBounds)] <- last
  
  labPos <- vector(mode= "numeric", length = length(span))
  
  i <- 1
  
  while (i < length(xBounds))
  {
    labPos[i] <- (xBounds[i] + xBounds[i+1]) / 2	
    i <- i+1
  }
  
  axis(side=1, at= xBounds, labels = FALSE, 
       lty=xlty, tck=in.folder$graphics.par$xtck, 
       cex.axis=in.folder$graphics.par$cex.axis)
  
  axis(side=1, at= labPos, labels = labels, 
       lty=xlty, tck=0)	
}

drawChromeBounds <- function(xValues)
{
  yBottomSegment <- par("usr")[3]
  iHeight <- abs((par("usr")[4] - par("usr")[3]) * 0.015)
  yTopSegment <- yBottomSegment + iHeight
  
  yHorizPos <- yTopSegment - (iHeight/2)
  
  xValues <- c(as.integer("0"), xValues)
  segments(xValues, yTopSegment, xValues, yBottomSegment, lty = "solid", col = "black", lwd=1)
  segments(min(xValues), yHorizPos, max(xValues), yHorizPos, lty = "solid", col = "black", lwd=1)	
}

getChromPosition <- function(lSpecies, lBounds = NULL, uBounds = NULL)
{
  if(is.null(lBounds))
    lBounds <- 1
  
  if(is.null(uBounds))
    uBounds <- length(lSpecies$genes$chrom)
  
  # in case bounds are passed in reverse order
  if (lBounds > uBounds)
  {
    tmp <- lBounds
    lBounds <- uBounds
    uBounds <- tmp
  }
  
  if (lSpecies$name == "human")
  {
    tar <- c(as.character(1:22), "x", "y")
  }else{
    chrom <- lSpecies$genes$chrom
    tar <- unique(lSpecies$genes$chrom)
  }
  
  idxL <- vector(mode = "integer", length = length(tar))
  idxU <- vector(mode = "integer", length = length(tar))
  idxChromLabelPos <- vector(mode = "integer", length = length(tar))
  names <- vector(mode = "character", length = length(tar))
  lTmp <- 0
  uTmp <- 0
  count <- 0
  # set up vectors for indices of upper and lower chromosome boundaries along x-axis
  for (i in 1:length(tar))
  {
    uTmp <- max(which(lSpecies$genes$chrom == tar[i]))
    if (uTmp <= lBounds)
      next
    
    if (uTmp > uBounds)
      uTmp <- uBounds
    
    lTmp <- min(which(lSpecies$genes$chrom == tar[i]))
    
    if (lTmp >= uBounds)
      break
    
    if (lTmp < lBounds)
      lTmp <- lBounds
    
    if ((lTmp >= lBounds) && (uTmp <= uBounds) && (lTmp < uTmp))
    {
      count <- count + 1
      idxL[count] <- geneToGenome(lTmp, lSpecies) 
      idxU[count] <- geneToGenome(uTmp, lSpecies)
      idxChromLabelPos[count] <- (idxL[count] + idxU[count]) / 2 	
      names[count] <- tar[i]
    }
  }	
  
  if (length(idxL) > count)
  {
    length(idxL) <- count
    length(idxU) <- count
    length(idxChromLabelPos) <- count
    length(names) <- count
  }
  
  return(data.frame(lower = idxL, upper = idxU, mid = idxChromLabelPos, name = names, stringsAsFactors = FALSE))
}

plotByChrom <- function(lSpecies, map, yLabel = 'yAxis')
{
  chromP <- getChromPosition(lSpecies)
  
  maxY <- max(map) * 1.05 # maximum y value + 5% to provide padding in plot
  yLim <- c(0, maxY)
  
  plot(map, xlab = "Chromosome", xaxt = "n", ylab = yLabel, type="l", ylim = yLim)
  drawChromeBounds(chromP$upper)
  axis(side = 1, at = chromP$mid, labels = tar, tick = FALSE)
}

## only pass range of data that we are interested in 
retrieveGenesOfInterest <- function(lSpecies, geneMap, nThreshold, lBound, uBound)
{
  
  if(uBound < lBound)
  {
    tmp <- lBound
    lBound <- uBound
    uBound <- tmp
  }
  
  if (nThreshold != -1)
  {
    indices <- which(geneMap >= nThreshold)
    
    if (length(indices) > 0)
    {
      # geneIndices are relative to the current subset of data being looked at 
      geneIndices <- indices
      # indices are relative to lBound, must add lBound to indices (-1) to get absolute indices
      indices <- indices + lBound - 1
      names <- lSpecies$genes$name[indices] 
      
      gOI <- data.frame(genomeAxis = indices, geneAxis = geneIndices, name = names, stringsAsFactors = FALSE)			
    }else
    {
      gOI <- data.frame(genomeAxis = -1, geneAxis = -1, name = "", stringsAsFactors = FALSE)
    }
  }else
  {
    gOI <- data.frame(genomeAxis = -1, geneAxis = -1, name = "", stringsAsFactors = FALSE)
  }
  return(gOI)
}

initDIANA <- function(sPath = '')
{
  if (sPath == '')
    sPath <- getwd()
  
  assign("DIANA_VERSION", 0.1, inherits=FALSE, envir=.GlobalEnv)
  assign("fileList", NULL, envir = .GlobalEnv)
  assign("fileDirectory", NULL, envir = .GlobalEnv)
  assign("resultList", NULL, envir = .GlobalEnv)
  assign("xMin", NULL, envir = .GlobalEnv)
  assign("yMin", NULL, envir = .GlobalEnv)
  assign("xMax", NULL, envir = .GlobalEnv)
  assign("yMax", NULL, envir = .GlobalEnv)
  assign("bZoom", FALSE, envir = .GlobalEnv)
  assign("lSpecies", NULL, envir = .GlobalEnv)
  assign("currImg", NULL, env = .GlobalEnv)
  
  #	source(paste(sPath, '/2016_08_11_digestR_2.0.0_plus_BRMB.r', collapse='', sep=''))
  #	source(paste(sPath, '/2016_03_16_mzPlot.r', collapse='', sep=''))
}

getGenotypeList <- function(resultList)
{
  sGenotypeList <- vector(mode="character", length=length(resultList))
  for(i in 1:length(resultList))
  {
    sGenotypeList[i] <- resultList[[i]]$Genotype
  }
  return(sGenotypeList)
}

generateDIANAFileName <- function(result)
{
  sFileName <- ""
  
  if (!is.null(result))
  {
    sName <- result$File
    ## find index of final period
    lastPeriod <- gregexpr(".", sName, fixed = TRUE)[[1]][length(gregexpr(".", sName, fixed = TRUE)[[1]])]
    ## find index of final slash
    lastSlash <- gregexpr("/", sName, fixed = TRUE)[[1]][length(gregexpr("/", sName, fixed = TRUE)[[1]])]
    ## extract filename and add genotype to beginning
    sFileName <- paste(substring(sName, lastSlash+1, lastPeriod), "dcf", sep = "", collapse="")	
  }
  return(sFileName)
}

testDCF <- function()
{
  result <- readDIANA()
  lSpecies <- get("lSpecies", envir=.GlobalEnv)
  plotByChrom(lSpecies, result$geneMap)
  windows()
  plotByGene(lSpecies, "HBB", result$aaMap)
}

loadSpecies <- function(fileName)
{
  lSpecies <- NULL
  switch(fileName,
         'chrHs.csv'=
           {
             lSpecies <- loadHumanGenes(fileName)
           },
         'sickle.csv'=
           {
             lSpecies <- loadHumanGenes(fileName)
           },
         'chr3D7.csv'=
           {
             lSpecies <- loadPfGenes(fileName)
           },
         'pa14.csv'=
           {
             lSpecies <- loadPA14Genes(fileName)
           },
         'pa01.csv'=
           {
             lSpecies <- loadPAO1Genes(fileName)
           }, 	
         'chrBtaurus1.csv'=
           {
             lSpecies <- loadTaurusGenes(fileName)
           }
  )
  return(lSpecies)
}

batchConvert <- function()
{
  
  ## load a species object based on the user choice made in processGUI
  lSpecies <- loadSpecies(globalSettings$speciesFiles[globalSettings$processSpeciesID])
  
  if(is.null(lSpecies))
  {
    return(0)
  }else
  {
    myAssign(in.name = "species", in.object = lSpecies)
  }
  
  ## get user input on the file(s) for processing
  dir <- choose.dir()
  
  if (dir.exists(dir))
  {
    if (!globalSettings$processSingleFile)
    {
      subDirs <- list.dirs(path = file.path(dir), full.names = TRUE, recursive = FALSE)
      if (length(subDirs) == 0)
        subDirs <- dir
    }else
    {
      subDirs <- ''
    }
    
    for (i in 1:length(subDirs))
    {
      if(!globalSettings$processSingleFile)
      {
        lFilesCSV <- list.files(path = subDirs[i], pattern = '*.csv', all.files=TRUE, full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
      }else
      {
        lFilesCSV <- myOpen(initialdir = dir, multiple = TRUE)
      }
      #			startTime <- Sys.time()
      if (length(lFilesCSV) > 0) ## only process if there are files selected
      {
        for (j in 1:length(lFilesCSV))
        {
          
          result <- processFile(lFilesCSV[j], lSpecies)
          saveFileName <- strsplit(lFilesCSV[j], '.', fixed = TRUE)
          saveFileName <- paste(unlist(saveFileName[[1]][1:length(saveFileName[[1]])-1]), collapse='')
          saveFileName <- paste(saveFileName, 'dcf', sep='.')
          
          if(is.numeric(result))
          {
            if(result == -1)
              print(paste0('Mascot file could not be read or found, ', saveFileName, " could not be generated."))
            else
              print(paste0('No matches found, ', saveFileName, " could not be generated."))	
          }else
          {
            sName <- writeDIANA(saveFileName, lSpecies, result) 
            print(paste0(sName, ' saved at ', format(Sys.time(), "%H:%M"), '.'))
          }
        }
      }
      #			endTime <- Sys.time()
      #			print(endTime - startTime)
    }
  }
}

renameFiles <- function()
{
  dir <- tkchooseDirectory()
  
  lFilesDCF <- list.files(path = file.path(dir), pattern = '*.DCF', all.files = FALSE, full.names = FALSE, recursive = TRUE, ignore.case = TRUE)
  
  for (i in 1:length(lFilesDCF))
  {
    innerDir <- extractDir(lFilesDCF[i])
    fileName <- extractFileName(lFilesDCF[i])
    
    vIndices <- gregexpr('_', fileName, fixed = TRUE)
    if(vIndices[[1]][1] > 0)
    {
      newFileName <- substr(fileName, vIndices[[1]][1] + 1, nchar(fileName))
      newFileName <- paste(dir, "//", innerDir, newFileName, sep="", collapse="")
      oldFileName <- paste(dir, "//", lFilesDCF[i], sep="", collapse="")
      file.rename(oldFileName, newFileName)
    }
    
  }
}

compareFileLists <- function()
{
  dir <- tkchooseDirectory()
  
  subDirs <- list.dirs(path = file.path(dir), recursive = FALSE)
  
  df <- data.frame(filename = character(), bProcessed = logical(), stringsAsFactors = FALSE)
  
  for (i in 1:length(subDirs))
  {
    lFilesCSV <- list.files(path = file.path(subDirs[i]), pattern = '*.CSV', all.files = TRUE, full.names = FALSE, recursive = FALSE, ignore.case = TRUE)
    df <- rbind(df, data.frame(filename = paste(subDirs[i], '/', sep="", collapse=""),	bProcessed = FALSE), stringsAsFactors = FALSE)
    
    if (length(lFilesCSV) > 0)
    {
      
      sPath <- paste(subDirs[i], "/", "DIANA", sep="", collapse="")
      lFilesDCF <- list.files(path = sPath, pattern = '*.DCF', all.files=TRUE, full.names = FALSE, recursive = FALSE, ignore.case = TRUE)
      lDCFNames <- vector(mode="character", length = length(lFilesDCF))
      
      if (length(lFilesDCF) > 0)
      {
        for (j in 1:length(lFilesDCF))
        {
          lDCFNames[j] <- extractFileNameNoExt(lFilesDCF[j])
        }
        
        for (j in 1:length(lFilesCSV))
        {
          df <- rbind(df, data.frame(filename = paste(subDirs[i], '/', lFilesCSV[j], sep="", collapse=""), 
                                     bProcessed = (extractFileNameNoExt(lFilesCSV[j]) %in% lDCFNames), stringsAsFactors = FALSE))	
        }				
      }
    }
  }
  return(df)	
  saveFileName <- mySave(defaultextension = '.csv', filetypes = list('csv' = 'Comma separated file'),	initialfile = '', title= 'Save')
  write.csv(df, saveFileName)
}

## assumes human species for now!!!
convertFile <- function()
{
  
  src <- myOpen("csv", list(csv = "Comma Separated Values File, xls = Excel File"), multiple = FALSE)
  dst <- mySave(defaultextension='.dcf', title= 'Save',
                filetypes=list('dcf'='DIANA Compressed File'))
  
  if ((src != '') & (file.exists(src)) & (dst != ''))
  {
    lHuman <- get("lSpecies", envir = .GlobalEnv)
    result <- alignAndCalculate(src, lHuman)
    if (!is.null(result))
      writeDIANA(dst, lHuman, result)
  }
}

convertMascotResult <- function()
{
  lSpecies <- get("lSpecies", envir = .GlobalEnv)
  in.file <- '' 
  result <- NULL
  
  in.file <- myOpen("csv", list(csv = "Comma Separated Values File, xls = Excel File"), multiple = FALSE)
  
  if (in.file != '')
    result <- alignAndCalculate(in.file, lSpecies)
  
  if (!is.null(result))
  {
    outPath <- myDir()
    if (outPath != '')
      writeDIANA(paste(outPath, generateDIANAFileName(result), sep = "/", collapse=""), lSpecies, result)
  }
}

buildFileFolder <- function()
{
  fileFolder <- list(
    file.par = buildFileFolderPar(),
    graphics.par = buildGraphicsPar(),
    data = NULL, w1 = NULL,	w2 = NULL)
  
}

buildFileFolderPar <- function()
{
  file.par <- list(	file.name = '', file.size = -1, date.modified = 0, 
                    ##axis = '', nucleus = '',
                    ##						matrix_size = as.integer(-1), block_size = as.integer(-1), 
                    upfield_ppm = -1, downfield_ppm = -1, 
                    center_ppm = -1, matrix_location = -1, 
                    species_location = -1, match_list_location = -1, endian = 'little', number_dimensions = as.integer(1), 
                    min_intensity = 0, max_intensity = 0, zero_offset = -1, user_title = '', genotype = '',
                    noise_override = -1, noise_multiplier = 3,
                    noise_est = -1, noise_max = -1, noise_mean = -1, noise_low_hinge = -1, noise_upper_hinge = -1,  
                    noise_min = -1, noise_sd = -1, aa_noise_est = -1, aa_noise_max = -1, aa_noise_mean = -1,
                    aa_noise_low_hinge = -1, aa_noise_upper_hinge = -1, aa_noise_min = -1, aa_noise_sd = -1,
                    aa_min_intensity = -1,	aa_max_intensity = -1, 	aa_zero_offset = -1, 
                    ##						aa_downfield_ppm = -1, aa_upfield_ppm = -1, aa_center_ppm = -1, 
                    files_merged = 0, genotype_override = '', title_override = '', force_reload = FALSE, 
                    appended_vectors = 0, vecInfo = data.frame(len = 0, loc = 0, label = '', stringsAsFactors = FALSE))
  return(file.par)
}

buildGraphicsPar <- function()
{
  ## Create defaultSettings
  graphics.par <- list(adj=0.5, ann=TRUE, ask=FALSE, bg="black", bty="7", cex=1, 
                       cex.axis=.95, cex.lab=1, cex.main=1, cex.sub=1, col="white",
                       col.axis="white", col.lab="white", col.main="white", col.sub="white", 
                       crt=0, err=0, family="", fg="white", fig=c(0, 1, 0, 1),
                       fin=c(10, 6.98958333333333), font=1, font.axis=1, font.lab=1, font.main=2, 
                       font.sub=1, lab=c(5, 5, 7), las=0, lend="round", lheight=1, ljoin="round", 
                       lmitre=10, lty="solid", lwd=1, mai=c(1.02, 0.82, 0.82, 0.42), 
                       mar=c(2.25, 2.25, 2.5, 1), mex=1, mfcol=c(1, 1), mfg=c(1, 1, 1, 1), 
                       mfrow=c(1, 1),	mgp=c(3, 1, 0), mkh=0.001, new=FALSE, oma=c(0, 0, 0, 0), 
                       omd=c(0, 1, 0, 1), omi=c(0, 0, 0, 0), pch=1, pin=c(8.76, 5.14958333333333), 
                       plt=c(0.082, 0.958, 0.145931445603577, 0.882682563338301), ps=12, pty="m", 
                       smo=1, srt=0, tck=NA, tcl=-0.5, usr=c(0, 10, 0,	10), xaxp=c(2, 8, 3), 
                       xaxs="r", xaxt="s", xlog=FALSE, xpd=FALSE, yaxp=c(20, 80, 6), yaxs="r", 
                       yaxt="s", ylog=FALSE, pos.color="blue", neg.color="green", conDisp=c(TRUE, TRUE), 
                       nlevels=20,	clevel=6, type="auto", theta=10, phi=10, asp=4, position.1D=.1, 
                       offset=0, proj.color='yellow', proj.type="l", proj.mode=FALSE, proj.direct=1, 
                       filter=function(x){range(x)[which.max(abs(range(x)))]}, peak.disp=FALSE, 
                       peak.color='white', peak.cex=.7, peak.pch='x', peak.labelPos='top', 
                       peak.noiseFilt=0, thresh.1D=6, roiMain=TRUE, roi.multi=TRUE, roiMax=TRUE, 
                       roi.bcolor=c('red', 'white'), roi.tcolor=c('red', 'white'),	
                       roi.lwd=c(2, 1), roi.lty=c('solid', 'dashed'), roi.labelPos='center', 
                       roi.cex=c(.95, .95), roi.noiseFilt=2, roi.w1=0, roi.w2=0, roi.pad=15, 
                       xtck=NA, ytck=NA, size.sub=c(6, 1.5), size.main=c(12, 7.25), #size.main=c(6, 4.25), 	
                       size.multi=c(6, 4.25), mar.sub=c(.1, .1, 1.6, .1), 
                       mar.multi=c(0, 0, 0, 0), cex.roi.multi=1, cex.files.multi=1, 
                       cex.roi.sub=1, overlay.text=TRUE, autoBackup=TRUE, sdi=TRUE, update=TRUE,
                       wd=path.expand('~'))	
  return(graphics.par)
}

################################################################################
##                                                                            ##
##               (Binary) file reading utilities                 			  ##
##                                                                            ##
################################################################################
# returns file.par
dianaHead <- function (file.name = NULL, print.info = TRUE)
{
  endFormat <- "little"
  ## Get the path to the file
  if(is.null(file.name))
    file.name <- myOpen(defaultextension = 'dcf', filetypes = list('dcf' = 'DIANA Compressed File'), initialfile = "", title= 'Load', multiple = FALSE)
  if(file.name == '')
    return(invisible())
  
  ## Open connection to binary file and check the file format
  con <- myFile(file.name, "rb")
  
  ## Get file information
  fileInfo <- file.info(file.name)
  fsize <- fileInfo$size
  mtime <- fileInfo$mtime	
  
  header <- readHeader(con, endFormat)
  
  closeAllConnections()
  
  ## Make a new header with extracted data
  head <- buildFileFolderPar()
  head$file.name <- file.name
  head$file.size <- fsize
  head$date.modified <- mtime
  head$endian <- endFormat
  head$user_title <- header$fileName
  head$version <- header$version
  head$species <- header$species
  
  head$speciesID <- header$speciesID
  head$genotype <- header$genotype
  head$matrix_location <- header$start_sparse_matrix
  head$species_location <- header$start_species
  head$match_list_location <- header$start_match_list
  
  head$noise_est <- header$g_noiseEst
  head$noise_sd <- header$g_noiseSD
  head$noise_max <- header$g_noiseMax
  head$noise_mean <- header$g_noiseMean
  head$noise_low_hinge <- header$g_noiseLowHinge
  head$noise_upper_hinge <- header$g_noiseUpperHinge
  head$noise_min <- header$g_noiseMin
  
  head$min_intensity <- header$g_min_int
  head$max_intensity <- header$g_max_int
  
  head$zero_offset <- header$g_zero	
  head$downfield_ppm <- header$aa_downfield
  head$upfield_ppm <- header$aa_upfield
  head$center_ppm <- header$aa_center
  
  head$aa_noise_est <- header$aa_noiseEst
  head$aa_noise_sd <- header$aa_noiseSD
  head$aa_noise_max <- header$aa_noiseMax
  head$aa_noise_mean <- header$aa_noiseMean
  head$aa_noise_low_hinge <- header$aa_noiseLowHinge
  head$aa_noise_upper_hinge <- header$aa_noiseUpperHinge
  head$aa_noise_min <- header$aa_noiseMin
  
  head$aa_min_intensity <- header$aa_min_int
  head$aa_max_intensity <- header$aa_max_int
  
  head$aa_zero_offset <- header$aa_zero
  head$files_merged <- header$files_merged
  head$genotype_override <- header$genotype_override
  head$title_override <- header$title_override
  head$appended_vectors <- header$appended_vectors
  head$vecInfo <- header$vecInfo
  
  
  if(header$noise_multiplier != -1)
  {
    head$noise_override <- header$noise_multiplier		
  }
  
  return(  list(file.par = head ))
}

## Internal function diana
## Reads diana (.dcf) format files and returns the species and match matrix information
## file.name - name of file to be read, NULL opens a file slection window
## w1Range = Min and Max chemical shifts to be read in the indirect dimension,
##           NULL will return all values.
## w2Range = Min and Max chemical shifts to be read in the direct dimension
##           NULL will return all values.
## file.par - file header (from dianaHead()), NULL will read the header first
## notes: providing a file.par from memory is used to speed up graphics
## returns the designated region of a spectrum and associated parameters
diana <- function (file.par = NULL)
{
  ## Get the path to the file and ucsf header
  
  if(is.null(file.par) || is.null(file.par$file.name))
  {
    file.name <- myOpen(defaultextension = 'dcf', filetypes = list('dcf' = 'DIANA Compressed File'), initialfile = "", title= 'Load', multiple = FALSE)
  }else
    file.name <- file.par$file.name
  
  if(is.null(file.name) || !nzchar(file.name))
    return(invisible())
  
  if(file.par$force_reload == TRUE )
  {
    file.par <- dianaHead(file.name = file.name, print.info = FALSE)[[1]]
    file.par$force_reload = FALSE
  }
  
  aaMap <- retrieveAAVector(fileName = file.name, endian = file.par$endian, offset = file.par$matrix_location)
  
  mData <- mapToGene(species, aaMap)
  
  outFolder <- list(file.par = file.par, mapData = mData,  species = species)
  
  return(outFolder)
}

prepGenomePlot <- function(w2Range, lSpecies, mapData, noise_est = -1)
{
  seqLen <- lSpecies$genes$seqLength
  cumLen <- cumsum(seqLen)
  
  if (!is.null(w2Range))
  {	
    firstGenomeIndex <- as.integer(min(w2Range))
    lastGenomeIndex <- as.integer(max(w2Range))
  }else
  {
    firstGenomeIndex <- 1
    lastGenomeIndex <- cumLen[length(cumLen)]
  }
  
  firstGeneIndex <- genomeToGene(firstGenomeIndex, lSpecies)
  lastGeneIndex <- genomeToGene(lastGenomeIndex, lSpecies)
  
  if((firstGeneIndex == -1)||(lastGeneIndex == -1))
  {
    ## invalid indices...
    
  }
  #print(firstGeneIndex)
  #print(lastGeneIndex)	
  geneNames = lSpecies$genes$name[firstGeneIndex:lastGeneIndex]
  #print(geneNames)
  chromData <- getChromPosition(lSpecies, firstGeneIndex, lastGeneIndex)
  
  chromDetails <- data.frame( chromName = chromData$name,
                              chromEnd = chromData$upper,
                              labelPos = chromData$mid,
                              stringsAsFactors = FALSE)
  
  w2Range <- c(firstGenomeIndex, lastGenomeIndex)
  
  geneSpan <- data.frame(x1 = cumLen[firstGeneIndex:lastGeneIndex] - seqLen[firstGeneIndex:lastGeneIndex] + 1,
                         x2 = cumLen[firstGeneIndex:lastGeneIndex])
  
  geneAvgData <- mapData$geneMap[firstGeneIndex:lastGeneIndex]
  
  if (noise_est == -1)
  {
    myMsg('Insufficient data to generate noise estimate.', type='ok', icon='error', title='DIANA')
  }
  
  gOI <- retrieveGenesOfInterest(lSpecies, geneAvgData, noise_est, firstGeneIndex, lastGeneIndex)
  
  gOI$genomeAxis <- cumLen[gOI$genomeAxis] - seqLen[gOI$genomeAxis] + 1 		# absolute value  
  gOI$geneAxis <- gOI$genomeAxis - firstGenomeIndex + 1 # to convert absolute index into relative
  xAxis <- firstGenomeIndex:lastGenomeIndex	# visual range (relative value)
  data <- mapData$aminoAcidMap[firstGenomeIndex:lastGenomeIndex] #y values for visual range
  geneData <- list(names = geneNames, span = geneSpan, genesOfInterest = gOI)
  
  results <- list(genes = geneData, 
                  w2Range = w2Range, 
                  w2 = xAxis, 
                  data = data, 
                  chromData = chromDetails)
  
  return(results)
}


prepGenePlot <- function(geneNames, lSpecies, mapData, noise_est = -1)
{
  
  gOI <- data.frame(genomeAxis = -1, geneAxis = -1, name = "", stringsAsFactors = FALSE)
  
  if(regexpr(',', geneNames)[[1]][1] > 0)
  {
    sGeneNameString <- unlist(strsplit(geneNames, ','))
    numGenes <- length(sGeneNameString)
  }else
  {
    sGeneNameString <- geneNames
    numGenes <- 1
  }
  
  vecStart <- vector(mode="integer", length=numGenes)
  vecLen <- vector(mode="integer", length=numGenes)
  vecCumSum <- vector(mode="integer", length=numGenes)
  
  for(i in 1:numGenes)
  {
    idx <- which(lSpecies$genes$name == sGeneNameString[i])
    
    if(length(idx)>0)
    {
      vecStart[i] <- lSpecies$genes$seqStartIdx[idx]
      vecLen[i] <- lSpecies$genes$seqLength[idx]
    }else
    {
      numGenes <- 0
      sMsg <- paste('Gene', sGeneNameString[i], 'not found.', sep=' ') 
      myMsg(sMsg, type='ok', icon='error', title='DIANA')
      results <- list(genesOfInterest = gOI, w2Range = c(0,0), geneNames = '', w2 = 0, data = 0, chromData = 0, xAxisLabels = 0)
      return(results)
    }
  }
  
  vIndices <- vector(mode="integer", length=sum(vecLen[1:numGenes]))
  xAxisLabels <- vector(mode="integer", length=numGenes)
  j <- 1
  k <- 0
  
  for(i in 1:numGenes)
  {
    k <- k + vecLen[i]
    vIndices[j:k] <- seq(from = vecStart[i], to = (vecStart[i] + vecLen[i] - 1), by = 1)
    xAxisLabels[i] <- substr(lSpecies$seq, vIndices[j], vIndices[k])
    j <- k + 1
  }
  
  vecCumSum <- cumsum(vecLen)
  vecCumStart <- c(0, vecCumSum)
  
  chromDetails <- data.frame('chromName' = sGeneNameString, 'chromEnd' = vecCumSum, 'labelPos' =  floor((vecCumSum + vecCumStart[1:numGenes]) / 2), stringsAsFactors = FALSE)
  
  w2Range <- c(1, length(vIndices))
  
  xAxis <- 1:length(vIndices)
  
  geneAvgData <- mapData$aminoAcidMap[vIndices] #y values
  
  if (noise_est == -1)
  {
    myMsg('Insufficient data to generate noise estimate.', type='ok', icon='error', title='DIANA')
  }
  
  return(list(genesOfInterest = gOI, w2Range = w2Range, geneNames = geneNames, w2 = xAxis, data = geneAvgData, 
              chromData = chromDetails, xAxisLabels = xAxisLabels, aaIndices = vIndices))
}

prepChromPlot <- function(w2Range, lSpecies, mapData, noise_est = -1)
{
  seqLen <- lSpecies$genes$seqLength
  cumLen <- cumsum(seqLen)
  
  if (!is.null(w2Range))
  {	
    firstGenomeIndex <- as.integer(min(w2Range))
    lastGenomeIndex <- as.integer(max(w2Range))
  }else
  {
    firstGenomeIndex <- 1
    lastGenomeIndex <- cumLen[length(cumLen)]
  }
  
  firstGeneIndex <- genomeToGene(firstGenomeIndex, lSpecies)
  lastGeneIndex <- genomeToGene(lastGenomeIndex, lSpecies)
  
  geneNames = lSpecies$genes$name[firstGeneIndex:lastGeneIndex]
  
  chromData <- getChromPosition(lSpecies, firstGeneIndex, lastGeneIndex)
  
  chromDetails <- data.frame( chromName = chromData$name,
                              chromEnd = chromData$upper,
                              labelPos = chromData$mid,
                              stringsAsFactors = FALSE)
  
  w2Range <- c(firstGenomeIndex, lastGenomeIndex)
  
  xAxis <- cumLen[firstGeneIndex:lastGeneIndex] - seqLen[firstGeneIndex:lastGeneIndex] + 1
  
  geneAvgData <- mapData$geneMap[firstGeneIndex:lastGeneIndex] #y values
  
  if (noise_est == -1)
  {
    myMsg('Insufficient data to generate noise estimate.', type='ok', icon='error', title='DIANA')
  }
  
  gOI <- retrieveGenesOfInterest(lSpecies, geneAvgData, noise_est, firstGeneIndex, lastGeneIndex)
  
  if((length(gOI$genomeAxis) == 1 ) && (gOI$genomeAxis == -1)) # no genes of interest found
  {		
    ## do nothing
  }else
  {
    gOI$genomeAxis <- cumLen[gOI$genomeAxis] - seqLen[gOI$genomeAxis] + 1 
  }
  results <- list(genesOfInterest = gOI, w2Range = w2Range, geneNames = geneNames, w2 = xAxis, data = geneAvgData, chromData = chromDetails)
  return(results)	
}

##Internal graphics wrapper function drawNMR
## note: this function implements all of the lower level draw functions
##Draws an NMR spectrum from a binary connection.
##in.folder - Header of the file to be plotted, default is current spectrum
##w1Range - Chemical shift range in the indirect dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)          
##w2Range - Chemical shift range in the direct dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)
##pos.zlim - Min and max positive intensities to be displayed, default is the 
##          most recent setting used with the file, the format is c(lower,upper)
##neg.zlim - Max and min negative intensities to be displayed, default is the 
##          most recent setting used with the file, the format is c(lower,upper)
##type  - specifies the type of plot to be generated: 2D data can be draw as
##        'auto', 'image', 'contour', 'filled', and 'persp'.
##        Image plots are the fastest, contour plots are more detailed, filled 
##        produces a filled contour plot, persp generates a 3D perspective plot,
##        and auto (the default) switches between image and contour depending on 
##        the amount of data being displayed.
##        1D values can also be passed as a type. Arguments include 'l', 'p',
##        and 'b' for line, point, and both line and points, respectively. 1D 
##        spectra default to 'l'. 2D data, when passes any of the 1D arguments 
##        will invoke proj1D and 'l', 'p' or 'b' will be passed to proj1D.
##        If type = NULL, type is taken from the spectrum's last setting. 
##pos.color - color of positive contours, default is the most recent setting
##         for the file, see colors() for the many color options
##neg.color - color of negative contours, default is the most recent setting
##         for the file, see colors() for the many color options
##col    - color for 1D data and 1D slices/projections of 2D data
##note:  All 2D plots use pos.color for positive intensities, and neg.color for 
##       negative intensities. 3D perspective plots use pos.color for all data,
##       and all 1D plots use col.
##nlevels - the number of contour intervals to be drawn, the default is the most
##         recent setting used for the file
##xlab    - x axis label, default is the direct detected nucleus in PPM
##ylab    - y axis label, default is the indirect detected nucleus in PPM
##main    - main title for plot, default is the spectrum name
##conDisp    - logical vector, c(TRUE, TRUE) plots positive and negative 
##           contours, c(TRUE, FALSE) plots only positive, c(FALSE, TRUE) plots
##           only the negative contours, c(FALSE, FALSE) plots no contours
## bg, fg, col.axis, col.lab, col.main, col.sub, col - see par()
##add       - logical argument, TRUE adds new data to an existing plot, FALSE 
##            generates a new plot
##p.window  -  The window to be used, can be 'main', 'sub', 'multi', or 'stats'
##axes      - logical argument, TRUE makes pretty labels
##offset    - Numeric argument expressing the % of total z range with which to 
##            displace a spectrum. This is used to create stacked 1D spectra
##...       - Additional graphics parameters can be passed to par()
##returns   - a plot of a 2D NMR spectrum
drawPeptides <- function ( in.folder = fileFolder[[wc()]], w1Range, w2Range, 
                           pos.zlim, neg.zlim, type, pos.color, neg.color, nlevels, conDisp,
                           bg, fg, col.axis, col.lab, col.main, col.sub, col, xlab = NULL, ylab = NULL,
                           main = in.folder$file.par$user_title, add = FALSE, p.window = 'main', 
                           axes = TRUE, offset = 0, ...){
  
  ## Supply graphics parameters to data missing graphics.par
  if ( is.null(in.folder$graphics.par) )
    in.folder$graphics.par <- defaultSettings
  
  ## Update graphics.par with any specified parameters
  if( !missing(w1Range) )
    in.folder$graphics.par$usr[3:4] <- w1Range
  if( !missing(w2Range) )
    in.folder$graphics.par$usr[1:2] <- w2Range
  if( !missing(type) )
    in.folder$graphics.par$type <- type
  if( !missing(pos.color) )
    in.folder$graphics.par$pos.color <- pos.color 
  if( !missing(neg.color) )
    in.folder$graphics.par$neg.color <- neg.color
  if( !missing(nlevels) )
    in.folder$graphics.par$nlevels <- nlevels
  if( !missing(conDisp) )
    in.folder$graphics.par$conDisp <- conDisp
  if( !missing(bg) )
    in.folder$graphics.par$bg <- bg
  if( !missing(fg) )
    in.folder$graphics.par$fg <- fg
  if( !missing(col.axis) )
    in.folder$graphics.par$col.axis <- col.axis 
  if( !missing(col.lab) )
    in.folder$graphics.par$col.lab <- col.lab	
  if( !missing(col.main) )
    in.folder$graphics.par$col.mainv <- col.main
  if( !missing(col.sub) )
    in.folder$graphics.par$col.sub <- col.sub
  if( !missing(col) )
    in.folder$graphics.par$col <- in.folder$graphics.par$proj.color <- col
  
  ## Calculate pos/neg zlims if missing	
  if( missing(pos.zlim) )
    pos.zlim <- c( in.folder$file.par$noise_est * in.folder$graphics.par$clevel,
                   in.folder$file.par$noise_est * in.folder$graphics.par$clevel * 
                     in.folder$graphics.par$nlevels )
  if( missing(neg.zlim) )
    neg.zlim <- -(rev(pos.zlim))
  
  ##Set plotting window
  setWindow( 
    p.window = p.window, 
    bg = in.folder$graphics.par$bg, 
    fg = in.folder$graphics.par$fg, 
    col.axis = in.folder$graphics.par$col.axis, 
    col.lab = in.folder$graphics.par$col.lab, 
    col.main = in.folder$graphics.par$col.main, 
    col.sub = in.folder$graphics.par$col.sub, 
    col = in.folder$graphics.par$col
  )
  
  ## Remove labels when appropriate
  if(p.window == 'multi')
    xlab <- ylab <- main <- ''
  if(p.window == 'sub')
    xlab <- ylab <- ''
  
  if( in.folder$file.par$number_dimensions == 1 )
  {
    if (in.folder$graphics.par$plotAA)
    {
      plotAA( in.folder = in.folder, xlab = xlab, ylab = ylab, 
              main = main, add = add, axes = axes, 
              offset = offset, ...)			
    }else
    {
      plotGenes( in.folder = in.folder, xlab = xlab, ylab = ylab, 
                 main = main, add = add, axes = axes, 
                 offset = offset, ...)
      #		print("manually call plot1D")
    }
    
  }
  bringFocus(-1) #return focus to console
}

plotGenes <- function(	
    in.folder = fileFolder[[wc()]],
    w1Range=in.folder$graphics.par$usr[3:4], ## w1Range is the y-axis
    w2Range=in.folder$graphics.par$usr[1:2], ## w2Range is the x-axis
    col = in.folder$graphics.par$proj.color, 
    type = in.folder$graphics.par$type,	
    xlab = NULL, 
    ylab = NULL, 
    main = in.folder$file.par$user_title, 
    roiMax = globalSettings$roiMax, 
    add = FALSE, 
    axes = TRUE, 
    offset = 0,
    geneNames = globalSettings$geneDisp,
    ... )
{
  Rprof("profile.txt")
  ## Remove any non line types
  if(!any(type == c('l', 'b', 'p')))
    type <- 'l'
  
  if(nchar(in.folder$file.par$title_override) > 0)
    main <- in.folder$file.par$title_override
  
  ## Read in dataset if we don't have data
  if( is.null(c(in.folder$data, in.folder$w2)) )
  {
    force_reload <- in.folder$file.par$force_reload
    
    new.folder <- diana(file.par = in.folder$file.par)
    
    in.folder$file.par 	<- new.folder$file.par
    species 			<- new.folder$species
    mapData 			<- new.folder$mapData
    
    if(force_reload)
    {## refresh the ranges for plot
      usr <- c( new.folder$file.par$downfield_ppm[1],
                new.folder$file.par$upfield_ppm[1], 
                new.folder$file.par$min_intensity,
                new.folder$file.par$max_intensity )			
      
      fileFolder[[wc()]]$file.par <- new.folder$file.par
      fileFolder[[wc()]]$graphics.par$usr <- usr
      
      myAssign("fileFolder", fileFolder, save.backup=TRUE)
      w1Range <- usr[3:4]
      w2Range <- usr[1:2]
    }
    ## data has been loaded
    
    if (in.folder$file.par$noise_est == -1) ## a noise estimate was not generated during data processing
    {
      print('Insufficient data to generate noise estimate using FDR.')
      flush.console()		
      noise_est <- sd(mapData$geneMap) * globalSettings$sd_noise_multiplier
      
    }else
    {
      noise_est <- in.folder$file.par$noise_est 
    }
    
    
    if (in.folder$file.par$noise_override != -1)
    {
      noise_est <- noise_est * in.folder$file.par$noise_override
    }else
    {
      noise_est <- noise_est * in.folder$file.par$noise_multiplier
    }		
    
    ## check if geneNames have been specified for plotting
    if(nchar(geneNames)>0)
    {
      bPlotGeneLevel <- TRUE
      results <- prepGenePlot(geneNames, species, mapData, noise_est)
      
      if(length(results$geneNames) == 0)
        return
      
      w2Range <- results$w2Range
      xLabels <- results$xAxisLabels
      
    }else
    {
      bPlotGeneLevel <- FALSE
      ## check range of w2Range... w2Range is the x-axis
      ## it should never exceed the min / max of the data indices
      if(min(w2Range)<=0)
      {
        w2Range[1] <- 1
      }
      
      if(max(w2Range)>nchar(species$seq))
      {
        w2Range[2] <- nchar(species$seq)
      }
      
      results <- prepChromPlot(w2Range, species, mapData, noise_est)
    }
    
    in.folder$w2 <- results$w2
    in.folder$data <- results$data
    in.folder$geneNames <- results$geneNames
    in.folder$chromDetails <- results$chromData
    in.folder$file.par$gOI <- results$genesOfInterest				
  }
  
  ## Erase old plot if add is false
  if(!add)
  {
    if(is.null(xlab))
    {
      if(bPlotGeneLevel)
      {
        xlab <- paste0(paste0(toupper(substr(species$name, 1, 1)), paste0(substr(species$name, 2, nchar(species$name))))," Gene(s)")
        
      }else
        xlab <- paste0(paste0(toupper(substr(species$name, 1, 1)), paste0(substr(species$name, 2, nchar(species$name))))," Chromosome")
    }
    
    if(is.null(ylab))
      ylab <- 'Intensity'
    
    plot(0,0, axes=FALSE, type= 'n', xlab = "", ylab="", main=main, xaxs="r", yaxs="r")
    
    title(main = main, sub=NULL, xlab = xlab, ylab=ylab, line=2)
    box(col=in.folder$graphics.par$fg, bty='o')
    
    w2Range <- sort(w2Range)
    
    winMin <- min(c(w1Range[1], min(in.folder$data)))
    winMax <- max(c(w1Range[2], max(in.folder$data)))		
    
    if(globalSettings$plotStyle$Variance && bPlotGeneLevel && (globalSettings$vectorType$Mean %in% in.folder$file.par$vecInfo$label))
    {
      dfPoly <- getVariancePolygon(in.folder$file.par, results$aaIndices, in.folder$w2, in.folder$data, offset)
      minPoly <- min(dfPoly$y)
      maxPoly <- max(dfPoly$y)
      
      if(winMax < maxPoly)
        winMax <- maxPoly
      
      if((minPoly < winMin) && (minPoly < 0))
      {
        winMin <- minPoly - ((max(in.folder$data) - min(in.folder$data)) * 0.05) + offset
      }
    }
    
    
    
    if (w1Range[2] == max(in.folder$data))
    {
      winMax <-((max(in.folder$data) - min(in.folder$data)) * 0.05) + max(in.folder$data) + offset
    }else
      winMax <- w1Range[2] + offset
    
    par(usr = c(w2Range[1], w2Range[2], winMin, winMax))
    
    if(axes)
    {
      if (!is.na(in.folder$graphics.par$xtck) &&	in.folder$graphics.par$xtck == 1)
        xlty <- 2
      else
        xlty <- 1
      
      if (!is.na(in.folder$graphics.par$ytck) && 	in.folder$graphics.par$ytck == 1)
        ylty <- 2
      else
        ylty <- 1
      
      if(bPlotGeneLevel)
        drawAndLabelAxisGenes(in.folder, in.folder$w2, in.folder$chromDetails$chromEnd, 
                              in.folder$chromDetails$labelPos, in.folder$chromDetails$chromName, xlty, xLabels)
      else
        drawAndLabelAxis(in.folder, in.folder$chromDetails$chromEnd, 
                         in.folder$chromDetails$labelPos, in.folder$chromDetails$chromName, xlty)		
      
      axis(side=2, lty=ylty, tck=in.folder$graphics.par$ytck, cex.axis=in.folder$graphics.par$cex.axis)				
    }			
  }
  
  # plot extra vector info
  if(bPlotGeneLevel)
  {
    if(globalSettings$plotStyle$EndPoints && (globalSettings$vectorType$EndPoints %in% in.folder$file.par$vecInfo$label))
      plot_endPoints(in.folder$file.par, results$aaIndices, offset, in.folder$w2, col)
    
    if(globalSettings$plotStyle$Peptides && ((globalSettings$vectorType$PepPoints %in% in.folder$file.par$vecInfo$label)))
      plot_peptides(in.folder$file.par, results$aaIndices, offset, col)
    
    if(globalSettings$plotStyle$Variance && ((globalSettings$vectorType$Mean %in% in.folder$file.par$vecInfo$label)))
    {
      if(!exists('dfPoly'))
        dfPoly <- getVariancePolygon(in.folder$file.par, results$aaIndices, in.folder$w2, in.folder$data, offset)			
      
      myCol <- col2rgb(col)
      #			myCol[1] <- myCol[1] * (1 - 0.25)  # darken
      #			myCol[2] <- myCol[2] * (1 - 0.25)
      #			myCol[3] <- myCol[3] * (1 - 0.25)
      myCol[1] <- myCol[1] + (255 - myCol[1]) * 0.4  ## lighten
      myCol[2] <- myCol[2] + (255 - myCol[2]) * 0.4
      myCol[3] <- myCol[3] + (255 - myCol[3])	* 0.4		
      
      polygon(dfPoly$x, dfPoly$y, col = rgb(myCol[1], myCol[2], myCol[3], maxColorValue = 255))
    }
    
    if(globalSettings$plotStyle$Accumulation)
      lines(x=in.folder$w2, y=in.folder$data + offset, col=col, type=type)
    
  }else
  {
    
    lines(x=in.folder$w2, y=in.folder$data + offset, col=col, type=type)
    
    if ((length(in.folder$file.par$gOI$genomeAxis) > 1))
    {
      text(in.folder$file.par$gOI$genomeAxis, jitter(in.folder$data[in.folder$file.par$gOI$geneAxis] + offset), 
           pos = 4, col = col, labels = in.folder$file.par$gOI$name, cex = in.folder$graphics.par$cex.axis * 0.8)
      
      points(in.folder$file.par$gOI$genomeAxis, in.folder$data[in.folder$file.par$gOI$geneAxis] + offset, col = col, pch = 16, type="p")
    }
  }
  Rprof(NULL)
}

## TSB Plot Gene level Detail
##Internal graphics wrapper function plot1D
##Draws a 1D NMR spectrum from a binary connection  
##in.folder - Header of the file to be plotted, default is current spectrum
##w1Range - Z limits for the plot          
##w2Range - Chemical shift range in the direct dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)
##col    - color for 1D data and 1D slices/projections of 2D data
##type  -  Can be 'l' (line), 'p' (points), or 'b' (both line and points)
##xlab    - x axis label, default is the appropriate nucleus in PPM
##ylab    - y axis label,default is 'intensity'
##main    - main title for plot, default is the file name
##roiMax    - logical argument, TRUE plots a point on the maximum
##            visible signal in the window
##add       - logical argument, TRUE adds new data to an existing plot, FALSE 
##            generates a new plot
##axes      - logical argument, TRUE makes pretty labels
##offset    - Numeric argument expressing the % of total z range with which to 
##            displace a spectrum. This is used to create stacked 1D spectra
##note: offset and vertical position (set by vp()) are not equivalent. 
##      vp() resets the zero point of the plot without affecting the max of  the 
##      zlimit. Offset shifts a given plot up/down from the vp() specified zero.
##...       - Additinal graphics paramaters can be passed to par()
##returns   - a plot of a 2D NMR spectrum  
plotAA <- function(	in.folder = fileFolder[[wc()]],
                    w1Range=in.folder$graphics.par$usr[3:4], 
                    w2Range=in.folder$graphics.par$usr[1:2], 
                    col = in.folder$graphics.par$proj.color, 
                    type = in.folder$graphics.par$type,	
                    xlab = NULL, ylab = NULL, 
                    main = in.folder$file.par$user_title, 
                    roiMax = globalSettings$roiMax, 
                    add = FALSE, axes = TRUE, offset = 0, ... )
{
  ## Remove any non line types
  if(!any(type == c('l', 'b', 'p')))
    type <- 'l'
  
  ## Redirect 2D NMR data to slice/projection function
  if(in.folder$file.par$number_dimensions == 1 )
  {	
    ## Read in dataset if we don't have data
    if( is.null(c(in.folder$data, in.folder$w2)) )
    {
      new.folder <- diana(file.par = in.folder$file.par)
      #new.folder <- ucsf2D(file.name = in.folder$file.par$file.name,	w2Range = w2Range, file.par = in.folder$file.par, getNames = TRUE)
      in.folder$file.par <- new.folder$file.par
      species <- new.folder$species
      mapData <- new.folder$mapData
      #####
      
      if (in.folder$file.par$noise_est == -1) ## a noise estimate was not generated during data processing
      {
        myMsg('Insufficient data to generate noise estimate using FDR.', type='ok', icon='error', title='DIANA')
        noise_est <- sd(mapData$geneMap) * globalSettings$sd_noise_multiplier
        
      }else
      {
        noise_est <- in.folder$file.par$noise_est 
      }
      
      
      if (in.folder$file.par$noise_override != -1)
      {
        noise_est <- noise_est * in.folder$file.par$noise_override
      }else
      {
        noise_est <- noise_est * in.folder$file.par$noise_multiplier
      }
      
      
      results <- prepGenomePlot(w2Range, species, mapData, noise_est)
      #####
      
      
      
      #			results <- prepGenomePlot(w2Range, species, mapData)
      in.folder$w2 <- results$w2
      in.folder$data <- results$data
      in.folder$geneNames <- results$genes$names
      in.folder$chromDetails <- results$chromData
      in.folder$file.par$gOI <- results$genes$genesOfInterest
      geneSpan <- results$genes$span
    }
    ## Erase old plot if add is false ## re-order this code so we only check for !add once ... do after loading new data if necessary
    if(!add)
    {
      w2Range <- sort(w2Range) 
      if(is.null(xlab))
      {
        if(length(geneSpan$x1) <= 10)
          xlab <- paste0(paste0(toupper(substr(species$name, 1, 1)), paste0(substr(species$name, 2, nchar(species$name))))," Gene")
        else
          xlab <- paste0(paste0(toupper(substr(species$name, 1, 1)), paste0(substr(species$name, 2, nchar(species$name))))," Chromosome")
      }
      
      if(is.null(ylab))
        ylab <- 'Intensity'
      
      plot(0,0, axes=FALSE, type= 'n', xlab = "", ylab="", main=main, xaxs="r", yaxs="r", ...)
      
      #			w1Range <- c(in.folder$file.par$zero_offset - (w1Range[2] - in.folder$file.par$zero_offset) * globalSettings$position.1D, w1Range[2])
      
      title(main = main, sub=NULL, xlab = xlab, ylab=ylab, line=2)
      box(col=in.folder$graphics.par$fg, bty='o')
      
      #			par(usr = c(w2Range[1], w2Range[2], 0, (max(in.folder$file.par$aa_max_intensity) * 1.10)))
      par(usr = c(w2Range[1], w2Range[2], 0, (max(in.folder$data) * 1.10)))
      
      #in.folder$data
      #			par(usr = c(w2Range[1], w2Range[2], 0, (max(in.folder$data) * 1.05) + offset))
      #			par(...)
      
      if(axes)
      {
        if (!is.na(in.folder$graphics.par$xtck) &&	in.folder$graphics.par$xtck == 1)
          xlty <- 2
        else
          xlty <- 1
        
        if (!is.na(in.folder$graphics.par$ytck) && 	in.folder$graphics.par$ytck == 1)
          ylty <- 2
        else
          ylty <- 1
        
        if (length(geneSpan$x1) <= 10)
        {			
          print(w2Range)
          drawAndLabelGenes(in.folder, w2Range[1:2], geneSpan, in.folder$geneNames, xlty)
          
        }else
        {
          drawAndLabelAxis(in.folder, in.folder$chromDetails$chromEnd, in.folder$chromDetails$labelPos, in.folder$chromDetails$chromName, xlty)
        }
        
        axis(side=2, lty=ylty, tck=in.folder$graphics.par$ytck, cex.axis=in.folder$graphics.par$cex.axis)				
      }			
    }
    
    lines(x=in.folder$w2, y=in.folder$data + offset, col=col, type=type, lwd=2)
    
    #		drawAndLabelGenes(in.folder, geneSpan, in.folder$geneNames, xlty)
    
    #		text(in.folder$file.par$gOI$genomeAxis, 
    #			jitter(in.folder$data[in.folder$file.par$gOI$geneAxis] + offset), 
    #			pos = 4, col = col, labels = in.folder$file.par$gOI$name, 
    #			cex = in.folder$graphics.par$cex.axis * 0.8)
    #	
    #
    #		points(in.folder$file.par$gOI$genomeAxis, 
    #				in.folder$data[in.folder$file.par$gOI$geneAxis] + offset, col = col, 
    #				pch = 16, type="p")
    
    if( roiMax &&  dev.cur() != 2 )
    {
      mShift <- maxShift( in.folder )
      points( mShift$w2, mShift$Height)
    }
  }
}

## converts gene index to genome index (amino acid space)
geneToGenome <- function(geneIndex, lSpecies)
{
  return(lSpecies$genes$seqStartIdx[geneIndex])	
}

genomeToGene <- function(genomeIndex, lSpecies)
{
  lBound <- lSpecies$genes$seqStartIdx
  uBound <- lBound + lSpecies$genes$seqLength - 1
  uBound[length(uBound)] <- uBound[length(uBound)] + 1  ## increment the upper bound of last entry
  ## so search includes entire range
  
  for (i in 1:length(lSpecies$genes$name))
  {
    if((genomeIndex >= lBound[i]) & (genomeIndex <= uBound[i]))
      return(i)
  }
  return(-1)
}

###############################

gl <- function()
{
  ##creates main window
  tclCheck()
  dlg <- myToplevel('gl')
  if (is.null(dlg))
    return(invisible())
  
  tkwm.title(dlg, 'Gene Name Threshold')
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  
  ##create file list box
  fileFrame <- ttklabelframe(dlg, text='Files')
  fileList <- tclVar()
  fileNames <- names(fileFolder)
  tclObj(fileList) <- getTitles(fileNames)
  fileBox <- tklistbox(fileFrame, height=10, width=34, listvariable=fileList,
                       selectmode='extended', active='dotbox',	exportselection=FALSE, bg='white', 
                       xscrollcommand=function(...) tkset(xscr, ...), yscrollcommand=function(...) tkset(yscr, ...))
  
  xscr <- ttkscrollbar(fileFrame, orient='horizontal', command=function(...) tkxview(fileBox, ...))
  yscr <- ttkscrollbar(fileFrame, orient='vertical',	command=function(...) tkyview(fileBox, ...))
  
  if (length(fileNames) > 2)
  {
    for (i in seq(0, length(fileNames) - 1, 2))
      tkitemconfigure(fileBox, i, background='#ececff')
  }
  
  tkselection.set(fileBox, wc() - 1)
  tcl(fileBox, 'see', wc() - 1)
  
  ##switches spectra on left-mouse double-click
  onDouble <- function()
  {
    usrSel <- 1 + as.integer(tkcurselection(fileBox))
    if (length(usrSel))
      usrFile <- fileNames[usrSel]
    else
      usrFile <- NULL
    if (!is.null(usrFile) && currentSpectrum != usrFile){
      currentSpectrum <- usrFile
      myAssign('currentSpectrum', currentSpectrum)
      refresh(multi.plot=FALSE)
      tkwm.deiconify(dlg)
      tkfocus(fileBox)
      
      if(fileFolder[[currentSpectrum]]$file.par$noise_override != -1)
      {
        tclObj(threshVal) <- fileFolder[[currentSpectrum]]$file.par$noise_override
      }else
      {
        tclObj(threshVal) <- fileFolder[[currentSpectrum]]$file.par$noise_multiplier	
      }			
      
    }
  }
  
  tkbind(fileBox, '<Double-Button-1>', onDouble)
  
  threshFrame <- ttklabelframe(fileFrame, text='Threshold', padding=2)
  entryFrame <- ttkframe(threshFrame)
  
  if(fileFolder[[currentSpectrum]]$file.par$noise_override != -1)
  {
    threshVal <- tclVar(fileFolder[[currentSpectrum]]$file.par$noise_override)
  }else
  {
    threshVal <- tclVar(fileFolder[[currentSpectrum]]$file.par$noise_multiplier)	
  }
  
  threshEntry <- ttkentry(entryFrame, width=6, textvariable=threshVal)
  entryLab <- ttklabel(entryFrame, text=' Threshold Multiplier')
  
  ##create apply button
  onApply <- function()
  {
    threshVal <- suppressWarnings(as.numeric(tclObj(threshVal)))
    
    if ((threshVal > 0)&&(threshVal != fileFolder[[currentSpectrum]]$file.par$noise_multiplier))
    {
      fileFolder[[currentSpectrum]]$file.par$noise_override <<- threshVal
    }else
    {
      tclObj(threshVal) <- fileFolder[[currentSpectrum]]$file.par$noise_multiplier
      fileFolder[[currentSpectrum]]$file.par$noise_override <<- -1
    }
    
    myAssign('fileFolder', fileFolder, save.backup = TRUE)
    refresh()
  }
  apply <- ttkbutton(threshFrame, text='Apply', width=10, command=onApply)
  
  ##create default button
  onDefault <- function()
  {
    tclObj(threshVal) <- fileFolder[[currentSpectrum]]$file.par$noise_multiplier
    fileFolder[[currentSpectrum]]$file.par$noise_override <<- -1
    refresh(multi.plot=FALSE)
  }
  default <- ttkbutton(threshFrame, text='Default', width=10, command=onDefault)
  
  ##add widgets to fileFrame
  tkgrid(fileFrame, column=1, row=1, sticky='nswe', pady=c(6, 0),	padx=8)
  tkgrid(fileBox, column=1, row=1, sticky='nswe')
  tkgrid(yscr, column=2, row=1, sticky='ns')
  tkgrid(xscr, column=1, row=2, sticky='we')
  
  ##make fileFrame stretch when window is resized
  tkgrid.columnconfigure(dlg, 1, weight=1)
  tkgrid.rowconfigure(dlg, 1, weight=10)
  tkgrid.columnconfigure(fileFrame, 1, weight=1)
  tkgrid.rowconfigure(fileFrame, 1, weight=1)
  
  
  tkgrid(threshFrame, column=1, columnspan=2, row=3, sticky='nwe', pady=2)
  tkgrid(entryFrame, column=1, columnspan=2, row=1, pady=c(0, 2))
  tkgrid(threshEntry, column=1, row=1)
  tkgrid(entryLab, column=2, row=1, sticky='w')
  tkgrid(apply, column=1, row=2, sticky='w')
  tkgrid(default, column=2, row=2, padx=4)
  
  tkgrid(ttksizegrip(dlg), column=2, row=2, sticky='se')
  
  ##Allows users to press the 'Enter' key to make selections
  onEnter <- function()
  {
    focus <- as.character(tkfocus())
    if (focus == '.pp.1.1')
      onDouble()
    else
      tryCatch(tkinvoke(focus), error=function(er){})
  }
  tkbind(dlg, '<Return>', onEnter)
  
  invisible()
}

pm <- function(species = globalSettings$speciesList)
{
  ##creates main window
  tclCheck()
  dlg <- myToplevel('pm')
  if (is.null(dlg))
    return(invisible())
  
  tkwm.title(dlg, 'Process Mascot Files')
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  
  speciesLabelFrame <- ttklabelframe(dlg, text='Species')
  
  speciesList <- tclVar()
  tclObj(speciesList) <- species
  
  tmp <- tclvalue(tkfont.actual(dlg))
  
  tmp2 <- strsplit(tmp, ' ', fixed=TRUE)
  fam <- tmp2[[1]][2]
  siz <- tmp2[[1]][4]
  wgt <- tmp2[[1]][6]
  slt <- tmp2[[1]][8]
  slt <- 'italic' ## change to italics
  
  myFont <- tkfont.create(family = fam, size = siz, weight = wgt, slant = slt)
  
  speciesBox <- tklistbox(speciesLabelFrame, height=10, width=34, listvariable=speciesList,
                          selectmode='extended', active='dotbox',	exportselection=FALSE, bg='white', font = myFont,
                          xscrollcommand=function(...) tkset(xscr, ...), yscrollcommand=function(...) tkset(yscr, ...))
  
  xscr <- ttkscrollbar(speciesLabelFrame, orient='horizontal', command=function(...) tkxview(speciesBox, ...))
  yscr <- ttkscrollbar(speciesLabelFrame, orient='vertical',	command=function(...) tkyview(speciesBox, ...))
  
  if (length(species) > 2)
  {
    for (i in seq(0, length(species) - 1, 2))
      tkitemconfigure(speciesBox, i, background='#ececff')
  }
  
  tkselection.set(speciesBox, 0)	
  #################
  ##create single vs. multiple file conversion radiobuttons
  fileConvFrame <- ttklabelframe(dlg, text='Convert File(s):')
  singleFileConvVal <- tclVar(TRUE)
  
  singleFileOnButton <- ttkradiobutton(fileConvFrame, variable=singleFileConvVal, 
                                       value=TRUE, text='Single File')
  
  singleFileOffButton <- ttkradiobutton(fileConvFrame, variable=singleFileConvVal, 
                                        value=FALSE, text='All Files in folder and sub folders')
  #################
  buttonFrame <- ttkframe(dlg)
  
  onApply <- function()
  {
    idx <- 1 + as.integer(tkcurselection(speciesBox))
    globalSettings$processSpeciesID <<- idx
    globalSettings$processSingleFile <<- (as.logical(tclObj(singleFileConvVal)) == TRUE)
    tkdestroy(dlg)
    batchConvert()
  }
  apply <- ttkbutton(buttonFrame, text='Apply', width=10, command=onApply)
  
  ##add widgets to speciesFrame
  tkgrid(speciesLabelFrame, column=1, row=1, sticky='nswe', pady=6, padx=6)
  tkgrid(speciesBox, column=1, row=1, sticky='nswe')
  tkgrid(yscr, column=2, row=1, sticky='ns')
  tkgrid(xscr, column=1, row=2, sticky='we')
  
  
  ##make fileFrame stretch when window is resized
  tkgrid.columnconfigure(dlg, 1, weight=1)
  tkgrid.rowconfigure(dlg, 1, weight=10)
  tkgrid.columnconfigure(speciesLabelFrame, 1, weight=1)
  tkgrid.rowconfigure(speciesLabelFrame, 1, weight=1)
  
  tkgrid(fileConvFrame, column=1, row=2, sticky='we', pady=6, padx=6)
  tkgrid(singleFileOnButton, column=1, row=1, sticky='nswe', pady=6, padx=6)
  tkgrid(singleFileOffButton, column=2, row=1, sticky='nswe', pady=6, padx=6)
  
  tkgrid(buttonFrame, column=1, row=3, sticky='nswe', pady=6, padx=6)
  tkgrid(apply, column=1, row=1, sticky='w')
  tkgrid.columnconfigure(buttonFrame, 1, weight=1)
  tkgrid.rowconfigure(buttonFrame, 1, weight=1)
  
  tkgrid(ttksizegrip(dlg), column=2, row=3, sticky='se')
  
  invisible()
}

sa <- function(saveFileName = '')
{
  if (saveFileName == '')
  {
    saveFileName <- mySave(title="Save As", defaultextension="dcf", filetypes=list('dcf'='DIANA Compressed File'))
    
    if (length(saveFileName) == 0 || !nzchar(saveFileName))
    {
      return(invisible())
    }
  }
  
  saveAs(saveFileName, fileFolder[[wc()]])
}


es <- function()
{	
  ##creates main window
  tclCheck()
  dlg <- myToplevel('es')
  if (is.null(dlg))
    return(invisible())
  
  tkwm.title(dlg, 'Edit Variables')
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  
  editStringFrame <- ttkframe(dlg)
  
  genotypeString <- tclVar()
  titleString <- tclVar()
  
  if (nchar(fileFolder[[wc()]]$file.par$genotype_override) > 0)
  {
    tclObj(genotypeString) <- fileFolder[[wc()]]$file.par$genotype_override
  }else
  {
    tclObj(genotypeString) <- fileFolder[[wc()]]$file.par$genotype
  }
  
  if(nchar(fileFolder[[wc()]]$file.par$title_override) > 0)
  {
    tclObj(titleString) <- fileFolder[[wc()]]$file.par$title_override
  }else
  {
    tclObj(titleString) <- fileFolder[[wc()]]$file.par$user_title
  }
  
  titleLabel <- ttklabel(editStringFrame, text='Title')
  genotypeLabel <- ttklabel(editStringFrame, text='Genotype')
  
  titleEntry <- ttkentry(editStringFrame, width=30, justify='left', 
                         textvariable=titleString)
  
  genotypeEntry <- ttkentry(editStringFrame, width=30, justify='left', 
                            textvariable=genotypeString)
  
  onApply <- function()
  {
    
    if(as.character(tclObj(titleString)) != fileFolder[[wc()]]$file.par$user_title)
      fileFolder[[wc()]]$file.par$title_override <- as.character(tclObj(titleString))
    
    
    if(as.character(tclObj(genotypeString)) != fileFolder[[wc()]]$file.par$genotype)
      fileFolder[[wc()]]$file.par$genotype_override <- as.character(tclObj(genotypeString))
    
    myAssign('fileFolder', fileFolder, save.backup = TRUE)
    
    tkdestroy(dlg)
    
    #	refresh()
  }
  apply <- ttkbutton(editStringFrame, text='Apply', width=10, command=onApply)
  
  ##create cancel button
  onCancel <- function() 
  {
    tkdestroy(dlg)
  }
  cancelButton <- ttkbutton(editStringFrame, text='Cancel', width=10, 
                            command=onCancel)	
  
  
  tkgrid(editStringFrame, column=1, row=1, sticky='nswe', pady=8, padx=2)
  tkgrid(titleLabel, column=1, row=1, padx=10, pady=4, sticky='nw')
  tkgrid(titleEntry, column=2, row=1, padx=10, pady=4, sticky='ne')
  tkgrid(genotypeLabel, column=1, row=2, padx=10, pady=4, sticky='nw')
  tkgrid(genotypeEntry, column=2, row=2, padx=10, pady=4, sticky='ne')
  tkgrid(apply, column=1, row=3, padx=10, pady=4, sticky='se')
  tkgrid(apply, column=2, row=3, padx=10, pady=4, sticky='se')
}

analyze_genes <- function(geneNames = '')
{
  globalSettings$geneDisp <- geneNames
  myAssign('globalSettings', globalSettings, save.backup = FALSE)
  zf()
  #	refresh()
}

plot_breaks <- function(geneNames = '')
{
  endian <- "little"
  loadFileName <- fileFolder[[wc()]]$file.par$file.name
  readCon <- file(loadFileName, "rb")
  header <- readHeader(readCon, endian)
  matches <- readMatches(readCon, endian, header$start_match_list)
}

plot_endPoints <- function(file.par, aaIndices, offset, xVals, col)
{
  vEndPoints <- getEndPoints(file.par, aaIndices)
  idx <- which(vEndPoints != 0)
  vEndPointsY <- vEndPoints[idx]
  vEndPointsX <- xVals[idx] 
  
  points(vEndPointsX, vEndPointsY + offset, col = col, pch = 16, type="p")	
}

plot_peptides <- function(file.par, aaIndices, offset, col)
{
  points <- getPeptidePoints(file.par, aaIndices)
  
  yLine <- 0.25
  prevX2 <- 0
  
  if(is.null(points))
  {
    return()
  }else
    print(length(points$x1))
  
  
  while(length(points$x1)>0)
  {
    i <- 1
    segments(points$x1[i], yLine, points$x2[i], yLine, col=col)
    prevX2 <- points$x2[i] + 1
    indices <- i
    
    j <- 2
    while(j <= length(points$x1))
    {
      if(points$x1[j] > prevX2)
      {
        segments(points$x1[j], yLine, points$x2[j], yLine, col=col)
        prevX2 <- points$x2[j] + 1
        indices <- c(indices, j)
      }
      j <- j + 1
    }
    points <- points[-indices,]
    yLine <- yLine + 0.75
  }	
}

getVariancePolygon <- function(file.par, aaIndices, xVals, yVals, offset)
{
  polygonYOffsets <- getVariancePoints(file.par, aaIndices)
  temp <- yVals + offset
  polygonY <- c(temp, rev(temp), temp[1])
  polyY <- polygonYOffsets + polygonY
  polygonX <- c(xVals, rev(xVals), xVals[1])
  
  return(data.frame(x = polygonX, y = polyY, stringsAsFactors = FALSE))
}


modalDialog <- function(parent, title, question, entryInit, entryWidth = 20,
                        returnValOnCancel = "ID_CANCEL") {
  dlg <- tktoplevel()
  tkwm.deiconify(dlg)
  tkgrab.set(dlg)
  tkfocus(dlg)
  tkwm.title(dlg, title)
  textEntryVarTcl <- tclVar(paste(entryInit))
  textEntryWidget <- ttkentry(dlg, width = paste(entryWidth),
                              textvariable = textEntryVarTcl)
  tkgrid(tklabel(dlg, text = question), textEntryWidget, padx = 10, pady = 15)
  returnVal <- returnValOnCancel
  
  onOK <- function() {
    returnVal <<- tclvalue(textEntryVarTcl)
    tkgrab.release(dlg)
    tkdestroy(dlg)
    tkfocus(parent)
  }
  
  onCancel <- function() {
    returnVal <<- returnValOnCancel
    tkgrab.release(dlg)
    tkdestroy(dlg)
    tkfocus(parent)
  }
  
  butOK <- ttkbutton(dlg, text = "OK", width = -6, command = onOK)
  butCancel <- ttkbutton(dlg, text = "Cancel", width = -6, command = onCancel)
  tkgrid(butCancel, butOK, padx = 10, pady = c(0, 15))
  
  tkfocus(dlg)
  tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg); tkfocus(parent)})
  tkbind(textEntryWidget, "<Return>", onOK)
  tkwait.window(dlg)
  
  returnVal
}

## Interactive GUI for manipulating plot settings
ps <- function(dispPane='co'){
  
  ##create main window
  current <- wc()
  tclCheck()
  dlg <- myToplevel('ps')
  if (is.null(dlg))
  {
    if (dispPane == 'co')
      tkselect('.ps.1', 0)
    else
      tkselect('.ps.1', 1)
    
    return(invisible())
  }
  
  tkwm.title(dlg, 'Plot Settings')
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  
  ##create paned notebook
  plotBook <- ttknotebook(dlg, padding=3)
  
  ##create plot settings panes
  coFrame <- ttkframe(plotBook, padding=c(0, 0, 2, 12)) 
  onedFrame <- ttkframe(plotBook, padding=c(0, 0, 2, 12))
  
  tkadd(plotBook, coFrame, text='Plot Colors')
  tkadd(plotBook, onedFrame, text='Spectra')
  
  ##add widgets to toplevel
  tkgrid(plotBook, column=1, row=1, sticky='nsew', padx=c(6, 0), pady=c(6, 0))
  tkgrid(ttksizegrip(dlg), column=2, row=2, sticky='se')
  tkgrid.columnconfigure(dlg, 1, weight=1)
  tkgrid.rowconfigure(dlg, 1, weight=1)		
  
  ##switch to the appropriate notebook pane
  if (dispPane == 'co')
    tkselect(plotBook, 0)
  else 
    tkselect(plotBook, 1)
  
  ####create widgets for coFrame
  ##create file list box
  coFileFrame <- ttklabelframe(coFrame, text='Files')
  coFileList <- tclVar()
  coFileNames <- names(fileFolder)
  tclObj(coFileList) <- getTitles(coFileNames)
  coFileBox <- tklistbox(coFileFrame, width=30, listvariable=coFileList, selectmode='extended', active='dotbox',	
                         exportselection=FALSE, bg='white', 	xscrollcommand=function(...) tkset(coXscr, ...), 
                         yscrollcommand=function(...) tkset(coYscr, ...))
  
  coXscr <- ttkscrollbar(coFileFrame, orient='horizontal', command=function(...) tkxview(coFileBox, ...))
  
  coYscr <- ttkscrollbar(coFileFrame, orient='vertical', command=function(...) tkyview(coFileBox, ...))
  
  if (length(coFileNames) > 2)
  {
    for (i in seq(0, length(coFileNames) - 1, 2))
      tkitemconfigure(coFileBox, i, background='#ececff')
  }
  
  tkselection.set(coFileBox, wc() - 1)
  tcl(coFileBox, 'see', wc() - 1)
  
  ##export fileBox selections to other tabs
  coSelect <- function(){
    usrSel <- 1 + as.integer(tkcurselection(coFileBox))
    onedFiles <- names(fileFolder)[which(sapply(fileFolder, function(x)
    {x$file.par$number_dimensions}) == 1)]
    
    if (!is.null(usrSel))
    {
      onedIndices <- na.omit(match(names(fileFolder)[usrSel], onedFiles))
      if (length(onedIndices))
      {
        tkselection.clear(onedFileBox, 0, 'end')
        
        for (i in onedIndices)
          tkselection.set(onedFileBox, i - 1)
      }
    }
    coConfigGui()
  }
  tkbind(coFileBox, '<<ListboxSelect>>', coSelect)
  
  ##switches spectra on left-mouse double-click
  coDouble <- function(){
    usrSel <- 1 + as.integer(tkcurselection(coFileBox))
    if (length(usrSel))
      usrFile <- coFileNames[usrSel]
    else
      usrFile <- NULL
    
    if (!is.null(usrFile) && currentSpectrum != usrFile)
    {
      myAssign('currentSpectrum', usrFile)
      refresh(multi.plot=FALSE)
      bringFocus()
      tkwm.deiconify(dlg)
      tkfocus(coFileBox)
    }
  }
  tkbind(coFileBox, '<Double-Button-1>', coDouble)
  
  ##create set axes color button
  coOptionFrame <- ttklabelframe(coFrame, text='Color options', padding=5)
  onAxes <- function()
  {
    usrSel <- 1 + as.integer(tkcurselection(coFileBox))
    if (length(usrSel))
      usrFiles <- coFileNames[usrSel]	
    else
      usrFiles <- currentSpectrum
    changeColor(dlg, 'axes', usrFiles)
  }
  axesColButton <- ttkbutton(coOptionFrame, text='Axes', width=11, command=onAxes)
  
  ##create set background color button
  onBg <- function()
  {
    usrSel <- 1 + as.integer(tkcurselection(coFileBox))
    if (length(usrSel))
      usrFiles <- coFileNames[usrSel]	
    else
      usrFiles <- currentSpectrum
    changeColor(dlg, 'bg', usrFiles)
  }
  bgColButton <- ttkbutton(coOptionFrame, text='BG', width=11, command=onBg)
  
  ##create set peak color button
  onPeak <- function()
  {
    usrSel <- 1 + as.integer(tkcurselection(coFileBox))
    if (length(usrSel))
      usrFiles <- coFileNames[usrSel]	
    else
      usrFiles <- currentSpectrum
    changeColor(dlg, 'peak', usrFiles)
  }
  peakColButton <- ttkbutton(coOptionFrame, text='Peak labels', width=11, command=onPeak)
  
  ##create set 1D color button
  onProj <- function()
  {
    usrSel <- 1 + as.integer(tkcurselection(coFileBox))
    if (length(usrSel))
      usrFiles <- coFileNames[usrSel]	
    else
      usrFiles <- currentSpectrum
    changeColor(dlg, 'proj', usrFiles)
  }
  projColButton <- ttkbutton(coOptionFrame, width=11, text='Plot', command=onProj)
  #	
  ##create print graphics button
  onContrast <- function(){
    usrSel <- 1 + as.integer(tkcurselection(coFileBox))
    if (length(usrSel))
      usrFiles <- coFileNames[usrSel]	
    else
      usrFiles <- currentSpectrum
    setGraphics(usrFiles, bg='white', line.color='black', pos.color='black', 
                neg.color='black', proj.color='black', peak.color='black', 
                roi.bcolor=c('red', 'black'), roi.tcolor=c('red', 'black'),
                refresh.graphics=TRUE)
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  contrastButton <- ttkbutton(coOptionFrame, text='High contrast', command=onContrast)
  
  ##create default colors button
  onDefaultColors <- function(){
    usrSel <- 1 + as.integer(tkcurselection(coFileBox))
    if (length(usrSel))
      usrFiles <- coFileNames[usrSel]	
    else
      usrFiles <- currentSpectrum
    setGraphics(usrFiles, bg=defaultSettings$bg, 
                line.color=defaultSettings$col.axis, 
                pos.color=defaultSettings$pos.color, 
                neg.color= defaultSettings$neg.color, 
                proj.color=defaultSettings$proj.color,
                peak.color=defaultSettings$peak.color,
                roi.bcolor=defaultSettings$roi.bcolor, 
                roi.tcolor=defaultSettings$roi.tcolor,
                refresh.graphics=TRUE)
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  defaultButton <- ttkbutton(coOptionFrame, text='Defaults', 
                             command=onDefaultColors)
  
  ##add widgets to fileFrame
  tkgrid(coFileFrame, column=1, row=1, sticky='nswe', pady=c(6, 4),	padx=8)
  tkgrid(coFileBox, column=1, row=1, sticky='nswe')
  tkgrid(coYscr, column=2, row=1, sticky='ns')
  tkgrid(coXscr, column=1, row=2, sticky='we')
  
  ##make fileFrame stretch when window is resized
  tkgrid.columnconfigure(coFrame, 1, weight=1)
  tkgrid.rowconfigure(coFrame, 1, weight=10)
  tkgrid.columnconfigure(coFileFrame, 1, weight=1)
  tkgrid.rowconfigure(coFileFrame, 1, weight=1)
  
  ##add widgets to coOptionFrame
  tkgrid(coOptionFrame, column=2, row=1, padx=c(4, 10))
  tkgrid(axesColButton, column=1, row=1, pady=c(3, 1), padx=1)
  tkgrid(bgColButton, column=2, row=1, pady=c(3, 1), padx=1)
  tkgrid(peakColButton, column=1, row=2, pady=1, padx=1)
  tkgrid(projColButton, column=2, row=2, pady=1, padx=1)
  
  tkgrid(contrastButton, column=1, columnspan=2, row=8, pady=c(16, 1), padx=25, sticky='we')
  tkgrid(defaultButton, row=8, column=1, columnspan=2, row=9, pady=6, padx=25, sticky='we')
  
  ##make optionFrame stretch when window is resized
  tkgrid.rowconfigure(coOptionFrame, 0, weight=1)
  tkgrid.rowconfigure(coOptionFrame, 9, weight=1)
  
  ##reconfigures widgets in GUI according to which spectra are open
  coConfigGui <- function()
  {
    usrSel <- 1 +	as.integer(tkcurselection(coFileBox))
    if (length(usrSel))
      usrFiles <- coFileNames[usrSel]
    else
      usrFiles <- currentSpectrum
  }
  coConfigGui()	
  
  ##resets file list and options whenever the mouse enters the GUI
  coMouse <- function()
  {
    reset(coFileList, coFileBox, coFileNames)
    coFileNames <<- names(fileFolder)
  }
  tkbind(coFrame, '<Enter>', coMouse)
  tkbind(coFrame, '<FocusIn>', coMouse)
  
  ####create widgets for onedFrame
  ##create file list box
  current <- wc()
  onedFileFrame <- ttklabelframe(onedFrame, text='Files')
  onedFileList <- tclVar()
  onedFileNames <- names(fileFolder)[which(sapply(fileFolder, 
                                                  function(x){x$file.par$number_dimensions}) == 1)]
  tclObj(onedFileList) <- getTitles(onedFileNames)
  onedFileBox <- tklistbox(onedFileFrame, width=30, listvariable=onedFileList, 
                           selectmode='extended', active='dotbox',	exportselection=FALSE, bg='white',
                           xscrollcommand=function(...) tkset(onedXscr, ...), 
                           yscrollcommand=function(...) tkset(onedYscr, ...))
  onedXscr <- ttkscrollbar(onedFileFrame, orient='horizontal',
                           command=function(...) tkxview(onedFileBox, ...))
  onedYscr <- ttkscrollbar(onedFileFrame, orient='vertical', 
                           command=function(...) tkyview(onedFileBox, ...))
  if (length(onedFileNames) > 2){
    for (i in seq(0, length(onedFileNames) - 1, 2))
      tkitemconfigure(onedFileBox, i, background='#ececff')
  }
  if (fileFolder[[current]]$file.par$number_dimensions == 1){
    tkselection.set(onedFileBox, match(currentSpectrum, onedFileNames) - 1)
    tcl(onedFileBox, 'see', match(currentSpectrum, onedFileNames) - 1)
  }
  
  ##export fileBox selections to other tabs
  onedSelect <- function(){
    usrSel <- 1 + as.integer(tkcurselection(onedFileBox))
    usrFile <- onedFileNames[usrSel]
    onedFiles <- names(fileFolder)[which(sapply(fileFolder, function(x)
    {x$file.par$number_dimensions}) == 1)] 
    selIndices <- match(onedFiles[usrSel], names(fileFolder))
    tkselection.clear(coFileBox, 0, 'end')
    for (i in selIndices)
      tkselection.set(coFileBox, i - 1)
    onedConfigGui()
  }
  tkbind(onedFileBox, '<<ListboxSelect>>', onedSelect)
  
  ##switches spectra on left-mouse double-click
  onedDouble <- function(){
    usrSel <- 1 + as.integer(tkcurselection(onedFileBox))
    if (length(usrSel))
      usrFile <- onedFileNames[usrSel]
    else
      usrFile <- NULL
    if (!is.null(usrFile) && currentSpectrum != usrFile){
      myAssign('currentSpectrum', usrFile)
      refresh(multi.plot=FALSE)
      bringFocus()
      tkwm.deiconify(dlg)
      tkfocus(onedFileBox)
      
      vecLabels <- fileFolder[[wc()]]$file.par$vecInfo$label
      numVecs <- length(vecLabels)
      
      statusEndPoints <- ifelse(globalSettings$vectorType$EndPoints %in% vecLabels, 'normal', 'disabled')
      statusPeptides <- ifelse(globalSettings$vectorType$PepPoints %in% vecLabels, 'normal', 'disabled')
      statusVariance <- ifelse(globalSettings$vectorType$Mean %in% vecLabels, 'normal', 'disabled')
      tkconfigure(cbPE, state=statusEndPoints)
      tkconfigure(cbPV, state=statusPeptides)
      tkconfigure(cbV, state=statusVariance)
      
    }
  }
  tkbind(onedFileBox, '<Double-Button-1>', onedDouble)
  
  ##creates switch number of dimensions radio buttons
  onedOptionFrame <- ttkframe(onedFrame)
  
  ##create plot type radiobuttons
  onedTypeFrame <- ttklabelframe(onedOptionFrame, text='Plot type')
  ptype <- switch(fileFolder[[current]]$graphics.par$type, 'auto'='line',
                  'p'='points', 'b'='both')
  if (is.null(ptype))
    ptype <- ''
  onedPlotType <- tclVar(ptype)
  onedType <- function(){
    usrSel <- 1 + as.integer(tkcurselection(onedFileBox))
    if (length(usrSel))
      usrFiles <- onedFileNames[usrSel]
    else{
      if (currentSpectrum %in% onedFileNames)
        usrFiles <- currentSpectrum
      else{
        tclObj(onedPlotType) <- 'line'
        err(paste('You must select a spectrum from the list before changing',
                  'plot options'), parent=dlg)
      }
    }
    pType <- switch(tclvalue(onedPlotType), 'line'='auto', 'points'='p', 
                    'both'='b')
    setGraphics(usrFiles, type=pType, refresh.graphics=TRUE)	
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  rbLine <- ttkradiobutton(onedTypeFrame, variable=onedPlotType, value='line',
                           text='Line', command=onedType)
  rbPoints <- ttkradiobutton(onedTypeFrame, variable=onedPlotType, 
                             value='points', text='Points', command=onedType)
  rbBoth <- ttkradiobutton(onedTypeFrame, variable=onedPlotType, value='both',
                           text='Both', command=onedType)
  
  ##create vertical position label
  vpFrame <- ttklabelframe(onedOptionFrame, text='Baseline')
  postn <- globalSettings$position.1D
  if(postn <= 1)
    postn <- (postn) / 2 * 100
  else
    postn <- 100 - (1 / postn) * 50 
  positionVal <- tclVar(postn)
  positionLab <- ttklabel(vpFrame, text='Position:')
  valLab <-	ttklabel(vpFrame, textvariable=positionVal, width=2)
  
  ##creates vertical position slider
  posSlider <- tkscale(vpFrame, from=99, to=0, variable=positionVal, 
                       orient='vertical', showvalue=F,	tickinterval=99, length=110, width=13, 
                       bg=as.character(tkcget(dlg, '-background')))
  onPosSlider <- function(){
    invisible(vp(as.numeric(tclObj(positionVal))))
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  tkbind(posSlider, '<ButtonRelease>', onPosSlider)	
  tkbind(posSlider, '<Return>', onPosSlider)	
  
  ##create default button
  onedDefault <- function()
  {
    usrSel <- 1 + as.integer(tkcurselection(onedFileBox))
    if (length(usrSel))
      usrFiles <- onedFileNames[usrSel]	
    else{
      if (currentSpectrum %in% onedFileNames)
        usrFiles <- currentSpectrum
      else
        err(paste('You must select a spectrum from the list before changing',
                  'plot options'), parent=dlg)
    }
    if (defaultSettings$type %in% c('p', 'b'))
      pType <- defaultSettings$type
    else
      pType <- 'l'		
    tclObj(onedPlotType) <- switch(pType, 'l'='line', 'p'='points', 'b'='both')
    postn <- defaultSettings$position.1D
    if (postn <= 1)
      postn <- (postn) / 2 * 100
    else
      postn <- 100 - (1 / postn) * 50 
    tclObj(positionVal) <- postn
    setGraphics(usrFiles, type=pType, proj.color=defaultSettings$proj.color, 
                save.backup=FALSE)
    invisible(vp(postn))
  }
  defaultButton <- ttkbutton(onedOptionFrame, text='Defaults', width=11, command=onedDefault)
  
  # create frame for gene level plot type
  genePlotTypeFrame <- ttklabelframe(onedOptionFrame, text='Gene Detail - Plot type')
  
  varCBPA <- tclVar(ifelse(globalSettings$plotStyle$Accumulation, 1, 0))
  varCBPE <- tclVar(ifelse(globalSettings$plotStyle$EndPoints, 1, 0))
  varCBPV <- tclVar(ifelse(globalSettings$plotStyle$Peptides, 1, 0))
  varCBV <- tclVar(ifelse(globalSettings$plotStyle$Variance, 1, 0))
  
  vecLabels <- fileFolder[[wc()]]$file.par$vecInfo$label
  numVecs <- length(vecLabels)
  
  statusEndPoints <- ifelse(globalSettings$vectorType$EndPoints %in% vecLabels, 'normal', 'disabled')
  statusPeptides <- ifelse(globalSettings$vectorType$PepPoints %in% vecLabels, 'normal', 'disabled')
  statusVariance <- ifelse(globalSettings$vectorType$Mean %in% vecLabels, 'normal', 'disabled')
  
  startCBPA <- tclvalue(varCBPA)
  startCBPE <- tclvalue(varCBPE)
  startCBPV <- tclvalue(varCBPV)
  startCBV <- tclvalue(varCBV)
  
  onPA <- function()
  {
    globalSettings$plotStyle$Accumulation <- !globalSettings$plotStyle$Accumulation
    myAssign( 'globalSettings', globalSettings, save.backup = FALSE )
  }
  cbPA <- ttkcheckbutton(genePlotTypeFrame, variable=varCBPA, command=onPA, text='Peptide Accumulation')
  
  onPE <- function()
  {
    globalSettings$plotStyle$EndPoints <- !globalSettings$plotStyle$EndPoints
    myAssign( 'globalSettings', globalSettings, save.backup = FALSE )
  }
  cbPE <- ttkcheckbutton(genePlotTypeFrame, variable=varCBPE, state=statusEndPoints, command=onPE, text='Peptide End Points')
  
  onPV <- function()
  {
    globalSettings$plotStyle$Peptides <- !globalSettings$plotStyle$Peptides
    myAssign( 'globalSettings', globalSettings, save.backup = FALSE )
  }
  cbPV <- ttkcheckbutton(genePlotTypeFrame, variable=varCBPV, state=statusPeptides, command=onPV, text='Peptide Visualization')
  
  onV <- function()
  {
    globalSettings$plotStyle$Variance <- !globalSettings$plotStyle$Variance
    myAssign( 'globalSettings', globalSettings, save.backup = FALSE )		
  }
  cbV  <-	ttkcheckbutton(genePlotTypeFrame, variable=varCBV, state=statusVariance, command=onV, text='Variance')
  
  onApply <- function()
  {
    if((startCBPA == tclvalue(varCBPA))&&(startCBPE == tclvalue(varCBPE))&&
       (startCBPV == tclvalue(varCBPV))&&(startCBV == tclvalue(varCBV)))
    {
      ## do nothing, the status of plot type is equal to the original configuration
    }else
    {
      print('Refreshing plot...')
      startCBPA <<- tclvalue(varCBPA)
      startCBPE <<- tclvalue(varCBPE)
      startCBPV <<- tclvalue(varCBPV)
      startCBV <<- tclvalue(varCBV)			
      refresh()
    }
  }
  applyButton <- ttkbutton(genePlotTypeFrame, text='Apply', width=11, command=onApply)
  
  
  onDisplayGene <- function()
  {
    geneName <- modalDialog(dlg, 'Gene Name Entry', 'Enter the name of the gene you wish to view in detail:', '')
    if(geneName == 'ID_CANCEL')
    {
      return()
    }else
    {
      if(nchar(geneName) == 0)
      {
        analyze_genes('')
      }
      else if(geneName %in% species$genes$name)
      {
        analyze_genes(geneName)
      }else
      {
        print(paste0(geneName, ' is not a valid gene of ', species$name))
      }
    }
    
  }	
  displayGeneButton <- ttkbutton(genePlotTypeFrame, text='Display Single Gene', width=21, command=onDisplayGene)
  
  onDisplayProteome <- function()
  {
    analyze_genes('')
  }
  displayGenomeButton <- ttkbutton(genePlotTypeFrame, text='Display Full Proteome', width=21, command=onDisplayProteome)
  
  ##add widgets to fileFrame
  tkgrid(onedFileFrame, column=1, row=1, sticky='nswe', pady=c(6, 4),	padx=8)
  tkgrid(onedFileBox, column=1, row=1, sticky='nswe')
  tkgrid(onedYscr, column=2, row=1, sticky='ns')
  tkgrid(onedXscr, column=1, row=2, sticky='we')
  
  ##make fileFrame stretch when window is resized
  tkgrid.columnconfigure(onedFrame, 1, weight=1)
  tkgrid.rowconfigure(onedFrame, 1, weight=10)
  tkgrid.columnconfigure(onedFileFrame, 1, weight=1)
  tkgrid.rowconfigure(onedFileFrame, 1, weight=1)
  
  ##add widgets to optionFrame
  tkgrid(onedOptionFrame, column=2, row=1, sticky='nswe', pady=c(10, 2), padx=c(4, 0))
  tkgrid(onedTypeFrame, column=1, row=2, padx=3, sticky='nsew')
  tkgrid(rbLine, column=1, row=1, padx=3, pady=6, sticky='w')
  tkgrid(rbPoints, column=1, row=2, padx=3, pady=6, sticky='w')
  tkgrid(rbBoth, column=1, row=3, padx=3, pady=6, sticky='w')
  
  tkgrid(defaultButton, column=1, row=3, padx=2, sticky='s')
  
  tkgrid(vpFrame, column=2, row=2, rowspan=2, padx=c(10, 15), sticky='nsew')
  tkgrid(positionLab, column=1, row=1, sticky='e', padx=c(1, 0))
  tkgrid(valLab, column=2, row=1, sticky='w', padx=2)
  tkgrid(posSlider, column=1, row=2, columnspan=2, padx=c(0, 5), sticky='ns')
  
  ##make optionFrame stretch when window is resized
  tkgrid.rowconfigure(onedOptionFrame, 0, weight=1)
  tkgrid.rowconfigure(onedOptionFrame, 4, weight=1)
  
  tkgrid(genePlotTypeFrame, column=1, row=5, columnspan=2, sticky='nsew')
  tkgrid(cbPA, column=1, row=1, padx=3, pady=6, sticky='w')
  tkgrid(cbPE, column=1, row=2, padx=3, pady=6, sticky='w')
  tkgrid(cbPV, column=1, row=3, padx=3, pady=6, sticky='w')
  tkgrid(cbV, column=1, row=4, padx=3, pady=6, sticky='w')
  tkgrid(applyButton, column=1, row=5, padx=3, pady=6, sticky='w')
  tkgrid(displayGeneButton, column=2, row=4, padx=3, pady=6, sticky='e')
  tkgrid(displayGenomeButton, column=2, row=5, padx=3, pady=6, sticky='e')
  tkgrid.rowconfigure(genePlotTypeFrame, 0, weight=1)
  tkgrid.rowconfigure(genePlotTypeFrame, 4, weight=1)
  
  ##reconfigures widgets in GUI according to which spectra are open
  onedConfigGui <- function(){
    usrSel <- 1 + as.integer(tkcurselection(onedFileBox))
    if (length(usrSel))
      usrFile <- onedFileNames[usrSel]
    else{
      if (length(names(fileFolder)) && currentSpectrum %in% onedFileNames)
        usrFile <- currentSpectrum
      else
        return(invisible())
    }
    if (length(usrFile) == 1){
      ptype <- switch(fileFolder[[usrFile]]$graphics.par$type, 'auto'='line',
                      'l'='line', 'p'='points', 'b'='both')
      if (is.null(ptype))
        ptype <- ''
      tclObj(onedPlotType) <- ptype
    }else{
      allEqual <- TRUE
      for (i in 2:length(usrFile)){
        if (fileFolder[[usrFile[i]]]$graphics.par$type != 
            fileFolder[[usrFile[i - 1]]]$graphics.par$type){
          allEqual <- FALSE
          break
        }
      }
      if (allEqual)
        tclObj(onedPlotType) <- 
          switch(fileFolder[[usrFile[1]]]$graphics.par$type, 'auto'='line', 
                 'l'='line', 'p'='points', 'b'='both')
      else
        tclObj(onedPlotType) <- ''	
    }
    postn <- globalSettings$position.1D
    if (postn <= 1)
      postn <- (postn) / 2 * 100
    else
      postn <- 100 - (1 / postn) * 50 
    tclObj(positionVal) <- postn
  }
  onedConfigGui()
  
  ##resets file list and options whenever the mouse enters the GUI
  onedMouse <- function(){
    reset(onedFileList, onedFileBox, onedFileNames, dims='1D')
    onedFileNames <<- names(fileFolder)[which(sapply(fileFolder, 
                                                     function(x){x$file.par$number_dimensions}) == 1)]
    onedConfigGui()
  }
  tkbind(onedFrame, '<Enter>', onedMouse)
  tkbind(onedFrame, '<FocusIn>', onedMouse)
  
  
  ##Allows users to press the 'Enter' key to make selections
  onEnter <- function()
  {
    focus <- as.character(tkfocus())
    if (length(grep('.1.1.1.1$', focus)))
      coDouble()
    else if (length(grep('.1.2.1.1$', focus)))
      onedDouble()
    
    else
      tryCatch(tkinvoke(focus), error=function(er){})
  }
  tkbind(dlg, '<Return>', onEnter) 
  
  ##enable\disable panes depending on which files are open
  onMouse <- function()
  {
    if (any(which(sapply(fileFolder, 
                         function(x){x$file.par$number_dimensions}) == 1)))
      tcl(plotBook, 'tab', 1, state='normal')
    else
      tcl(plotBook, 'tab', 1, state='disabled')
    
  }
  tkbind(dlg, '<Enter>', onMouse)
  tkbind(dlg, '<FocusIn>', onMouse)
  
  invisible()
}

## User graphics function ol
## Overlay open spectra onto the current spectrum
## askUser - Logical argument, TRUE opens the overlay GUI
## offset    - Numeric argument expressing the % of total z range with which to 
##            displace each spectrum. This is used to create stacked 1D spectra
##            and is not passed to 2D plots
## note: offset and vertical position (set by vp()) are not equivalent. 
##      vp() resets the zero point of the plot without affecting the max of the 
##      zlimit. Offset shifts a given plot up/down from the vp() specified zero.
## ...  - Additional plotting options can be passed to drawPeptides and par()
ol <- function(askUsr = TRUE, offset = NULL, ...)
{
  ### Define default overlay palletes
  p4col <- c(rgb(215,48,39,maxColorValue=256), rgb(253,174,97,maxColorValue=256), rgb(171,217,233,maxColorValue=256), rgb(69,117,180,maxColorValue=256))
  p6col <- c(rgb(215,48,39,maxColorValue=256), rgb(244,109,67,maxColorValue=256), rgb(253,174,97,maxColorValue=256), rgb(171,217,233,maxColorValue=256), rgb(116,173,209,maxColorValue=256), rgb(69,117,180,maxColorValue=256))
  p8col <- c(rgb(165,0,38,maxColorValue=256), rgb(215,48,39,maxColorValue=256), rgb(244,109,67,maxColorValue=256), rgb(253,174,97,maxColorValue=256), rgb(171,217,233,maxColorValue=256), rgb(116,173,209,maxColorValue=256), rgb(69,117,180,maxColorValue=256), rgb(49,54,149,maxColorValue=256))
  ###	
  
  ## Define current spectrum
  current <- wc()
  current.par <- fileFolder[[ current ]]$graphics.par
  c.nDim <- fileFolder[[current]]$file.par$number_dimensions
  
  ## Open GUI for making overlay list
  if(!exists('overlayList') )	
    myAssign('overlayList', NULL )
  
  if(askUsr==TRUE || is.null(overlayList)){
    os('ol')
    return(invisible())
  } 
  
  ## Fetch the offset parameter
  if( is.null(offset) )
    offset <- globalSettings$offset
  
  ## Remove the current spectrum from the overlay list
  if(currentSpectrum %in% overlayList)
    overlayList <- overlayList[-(which( overlayList == currentSpectrum))]
  if(length(overlayList) == 0)
    return(invisible())		
  
  ## Plot the overlay list
  o.nDim <- NULL
  
  if(globalSettings$overlay.text && globalSettings$overlay.textSuppress && length(gregexpr('/', fileFolder[[current]]$file.par$user_title)[[1]]) > 0)
  {
    startIdx <- gregexpr('/', fileFolder[[wc()]]$file.par$user_title)[[1]][length(gregexpr('/', fileFolder[[wc()]]$file.par$user_title)[[1]])] + 1
    endIdx <- gregexpr('\\.', fileFolder[[wc()]]$file.par$user_title)[[1]][length(gregexpr('\\.', fileFolder[[wc()]]$file.par$user_title)[[1]])] - 1
    plot.list <- substr(fileFolder[[wc()]]$file.par$user_title, startIdx, endIdx)
    
  }else
    plot.list <- fileFolder[[current]]$file.par$user_title
  
  if(c.nDim > 1)
    col.list <- current.par$pos.color
  else
    col.list <- current.par$proj.color	
  newset <- offset
  for(i in overlayList)
  {
    o.nDim <- fileFolder[[i]]$file.par$number_dimensions
    
    ## Overlay spectra on main plot 
    if(o.nDim == c.nDim)
    {
      ## DIANA -> PLOT 1D
      drawPeptides(fileFolder[[i]], type=current.par$type, add = TRUE, w1Range = current.par$usr[3:4], w2Range=current.par$usr[1:2], offset = newset, ...)	
    }else if (o.nDim == 1 && c.nDim > 1)## shouldn't happen in DIANA
    {
      ## Read 1D file
      in.folder <- fileFolder[[i]]
      data.folder <- ucsf1D(fileFolder[[i]]$file.par$file.name)
      
      ## Setup plot range
      in.folder$data <- data.folder$data
      in.folder$w2 <- data.folder$w2	
      newRange <- c(current.par$usr[1:2], min(data.folder$data), 
                    max(data.folder$data))
      in.folder$graphics.par$usr <- newRange
      op <- par('usr')
      par(usr=newRange)
      
      ## Plot 1D overlay
      plot1D(in.folder, add=TRUE, offset=newset,
             col=fileFolder[[i]]$graphics.par$proj.color, 
             type=fileFolder[[i]]$graphics.par$type)
      par(usr=op)
    }
    
    ## Keep track of overlaid spectra
    newset <- newset + offset
    
    if(globalSettings$overlay.text && globalSettings$overlay.textSuppress && length(gregexpr('/', fileFolder[[i]]$file.par$user_title)[[1]]) > 0)
    {
      startIdx <- gregexpr('/', fileFolder[[i]]$file.par$user_title)[[1]][length(gregexpr('/', fileFolder[[i]]$file.par$user_title)[[1]])] + 1
      EndIdx <- gregexpr('\\.', fileFolder[[i]]$file.par$user_title)[[1]][length(gregexpr('\\.', fileFolder[[i]]$file.par$user_title)[[1]])] - 1
      plot.list <- c(plot.list, substr(fileFolder[[i]]$file.par$user_title, startIdx, EndIdx))			
      
    }else
    {
      plot.list <- c(plot.list, fileFolder[[i]]$file.par$user_title)
    }
    
    if(o.nDim == 1)
    {
      col.list <- c(col.list, fileFolder[[i]]$graphics.par$proj.color)
      
    }else
      col.list <- c(col.list, 
                    fileFolder[[i]]$graphics.par$pos.color)		
  }
  
  ## Add a legend if there are files other than the current spectrum
  if( length(plot.list) > 1 && globalSettings$overlay.text)
  {
    legend("topleft", rev(plot.list), pch=NULL, bty='n', text.col = rev(col.list))
  }
}


## Interactive GUI for manipulating overlays and shift referencing
os <- function(dispPane='ol'){
  
  ##create main window
  current <- wc()
  tclCheck()
  dlg <- myToplevel('os')
  if (is.null(dlg))
  {
    if (dispPane == 'ol')
    {
      tkwm.title('.os', 'Overlays')
      tkselect('.os.1', 0)
    }
    return(invisible())
  }
  
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  if (dispPane == 'ol')
    tkwm.title(dlg, 'Overlays')
  
  ##create paned notebook
  osBook <- ttknotebook(dlg, padding=3)
  
  ##create overlay and referencing panes
  olFrame <- ttkframe(osBook, padding=c(0, 0, 2, 12)) 
  
  tkadd(osBook, olFrame, text='   Overlays   ')
  
  ##add widgets to toplevel
  tkgrid(osBook, column=1, row=1, sticky='nsew', padx=c(6, 0), pady=c(6, 0))
  tkgrid(ttksizegrip(dlg), column=2, row=2, sticky='se')
  tkgrid.columnconfigure(dlg, 1, weight=1)
  tkgrid.rowconfigure(dlg, 1, weight=1)
  
  ##switch to the appropriate notebook pane
  if (dispPane == 'ol')
    tkselect(osBook, 0)
  
  ####create widgets for olFrame
  ##create file list box
  olFileFrame <- ttklabelframe(olFrame, text='Files')
  olFileList <- tclVar()
  olFileNames <- names(fileFolder)
  overlayMatches <- match(overlayList, olFileNames)
  if (length(overlayMatches))
    olFileNames <- olFileNames[-overlayMatches]
  tclObj(olFileList) <- getTitles(olFileNames)
  olFileBox <- tklistbox(olFileFrame, height=13, width=25, 
                         listvariable=olFileList, selectmode='extended', active='dotbox', 
                         exportselection=FALSE, bg='white', 
                         xscrollcommand=function(...) tkset(olXscr, ...), 
                         yscrollcommand=function(...) tkset(olYscr, ...))
  olXscr <- ttkscrollbar(olFileFrame, orient='horizontal',
                         command=function(...) tkxview(olFileBox, ...))
  olYscr <- ttkscrollbar(olFileFrame, orient='vertical', 
                         command=function(...) tkyview(olFileBox, ...))
  if (length(olFileNames) > 2)
  {
    for (i in seq(0, length(olFileNames) - 1, 2))
      tkitemconfigure(olFileBox, i, background='#ececff')
  }
  currMatch <- match(currentSpectrum, olFileNames)
  if (!is.na(currMatch)){
    tkselection.set(olFileBox, currMatch - 1)
    tcl(olFileBox, 'see', currMatch - 1)
  }
  
  ##create add button
  middleFrame <- ttkframe(olFrame)
  buttonFrame <- ttkframe(middleFrame)
  onAdd <- function(){
    
    ##get selection
    usrSel <- 1 + as.integer(tkcurselection(olFileBox))
    if (!length(usrSel))
      err('You must select a file from the files list to overlay')
    
    ##update global object overlayList
    overlayList <- c(overlayList, olFileNames[usrSel])
    myAssign('overlayList', overlayList)
    refresh(sub.plot=FALSE, multi.plot=FALSE)
    
    ##update contents of overlayBox
    if (is.null(overlayList))
      overlayNames <<- character(0)
    else
      overlayNames <<- overlayList
    tclObj(overlaysList) <- getTitles(overlayNames)
    for (i in seq_along(usrSel))
      tkselection.set(overlayBox, length(overlayNames) - i)
    tcl(overlayBox, 'see', length(overlayNames) - 1)
    
    ##update contents of olFileBox
    olFileNames <<- names(fileFolder)
    overlayMatches <- match(overlayNames, olFileNames)
    if (length(overlayMatches))
      olFileNames <<- olFileNames[-overlayMatches]
    tclObj(olFileList) <- getTitles(olFileNames)
    tkselection.clear(olFileBox, 0, 'end')
    
    ##reconfigure GUI
    if (length(overlayNames) > 2){
      for (i in seq(0, length(overlayNames) - 1, 2))
        tkitemconfigure(overlayBox, i, background='#ececff')
    }
    olConfigGui()
    tkfocus(olFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  addButton <- ttkbutton(buttonFrame, text='Add -->', width=11, command=onAdd)
  
  ##create remove button
  onRemove <- function(){
    
    ##get selection
    usrSel <- 1 +	as.integer(tkcurselection(overlayBox))
    if (!length(usrSel))
      err('You must select a file from the overlays list to remove')
    selNames <- overlayNames[usrSel]
    
    ##update global object overlayList
    overlayList <- overlayList[-usrSel]
    if (!length(overlayList))
      overlayList=NULL
    myAssign('overlayList', overlayList)
    refresh(sub.plot=FALSE, multi.plot=FALSE)
    
    ##update contents overlayBox
    if (is.null(overlayList))
      overlayNames <<- character(0)
    else
      overlayNames <<- overlayList
    if (length(overlayNames))
      tclObj(overlaysList) <- getTitles(overlayNames)
    else
      tclObj(overlaysList) <- overlayNames
    tkselection.clear(overlayBox, 0, 'end')
    
    ##update contents of olFileBox
    olFileNames <<- names(fileFolder)
    overlayMatches <- match(overlayNames, olFileNames)
    if (length(overlayMatches))
      olFileNames <<- olFileNames[-overlayMatches]
    tclObj(olFileList) <- getTitles(olFileNames)
    tkselection.clear(olFileBox, 0, 'end')
    for (i in selNames)
      tkselection.set(olFileBox, match(i, olFileNames) - 1)
    tcl(olFileBox, 'see', match(i, olFileNames) - 1)
    
    ##reconfigure GUI
    if (length(overlayNames) > 2){
      for (i in seq(0, length(overlayNames) - 1, 2))
        tkitemconfigure(overlayBox, i, background='#ececff')
    }
    olConfigGui()
    tkfocus(olFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  removeButton <- ttkbutton(buttonFrame, text='Remove', width=11, 
                            state='disabled', command=onRemove)
  
  ##create offset label
  offsetFrame <- ttklabelframe(middleFrame, text='Offset')
  offsetVal <- tclVar(globalSettings$offset)
  offsetLab <- ttklabel(offsetFrame, text='Offset:')
  valLab <-	ttklabel(offsetFrame, textvariable=offsetVal, width=4)
  
  ##creates offset slider
  offsetSlider <- tkscale(offsetFrame, from=5, to=-5,	variable=offsetVal,
                          orient='vertical', showvalue=F,	tickinterval=50, resolution = 0.1,
                          bg=as.character(tkcget(dlg, '-background')))
  
  onOffset <- function(){
    setGraphics(offset=as.numeric(tclvalue(offsetVal)), refresh.graphics=TRUE)
    tkfocus(olFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  tkbind(offsetSlider, '<ButtonRelease>', onOffset)	
  
  ##create overlays list box
  overlayFrame <- ttklabelframe(olFrame, text='Overlays')
  if (is.null(overlayList)){
    overlayNames <- character(0)
  }else{
    overlayNames <- overlayList
  }
  overlaysList <- tclVar()
  if (!is.null(overlayList))
    tclObj(overlaysList) <- getTitles(overlayNames)
  else
    tclObj(overlaysList) <- character(0)
  overlayBox <- tklistbox(overlayFrame,	height=7, width=25, 
                          exportselection=FALSE, listvariable=overlaysList, selectmode='extended', 
                          active='dotbox', bg='white', 
                          xscrollcommand=function(...) tkset(overlayXscr, ...), 
                          yscrollcommand=function(...) tkset(overlayYscr, ...))
  overlayXscr <- ttkscrollbar(overlayFrame, orient='horizontal', 
                              command=function(...) tkxview(overlayBox, ...))
  overlayYscr <- ttkscrollbar(overlayFrame, orient='vertical',
                              command=function(...) tkyview(overlayBox, ...))
  if (length(overlayNames) > 2){
    for (i in seq(0, length(overlayNames) - 1, 2))
      tkitemconfigure(overlayBox, i, background='#ececff')
  }
  currMatch <- match(currentSpectrum, overlayNames)
  if (!is.na(currMatch)){
    tkselection.set(overlayBox, currMatch - 1)
    tcl(overlayBox, 'see', currMatch - 1)
  }
  
  ##switches spectra on left-mouse double-click
  olDouble <- function(box){
    usrSel <- 1 + as.integer(tkcurselection(box))
    if (length(usrSel)){
      if (box$ID == '.os.1.1.1.1')
        usrFile <- olFileNames[usrSel]
      else
        usrFile <- overlayList[usrSel]
    }else
      usrFile <- NULL
    if (!is.null(usrFile) && currentSpectrum != usrFile){
      myAssign('currentSpectrum', usrFile)
      refresh(multi.plot = FALSE)
      olConfigGui()
      tkwm.deiconify(dlg)
      tkfocus(box)
    }
  }
  tkbind(olFileBox, '<Double-Button-1>', function(...) olDouble(olFileBox))
  tkbind(overlayBox, '<Double-Button-1>', function(...) olDouble(overlayBox))
  
  ##create set positive contour color button
  colFrame <- ttklabelframe(olFrame, text='Overlay colors')
  
  ##create set 1D color button
  onProj <- function(){
    usrSel <- 1 +	as.integer(tkcurselection(overlayBox))
    usrFiles <- overlayList[usrSel]
    if (!length(usrFiles))
      err(paste('You must select a spectrum from the Overlays list before',
                'changing the plot color'))
    changeColor(dlg, '1D', usrFiles)
  }
  projColButton <- ttkbutton(colFrame, text='Plot', width=13, command=onProj)
  
  ##create overlay text checkbutton
  textFrame <- ttkframe(olFrame)
  textVal <- tclVar(ifelse(globalSettings$overlay.text, 1, 0))
  onText <- function()
  {
    if (as.logical(as.integer(tclvalue(textVal))))
    {
      setGraphics(overlay.text=TRUE, save.backup=TRUE)
      refresh(sub.plot=FALSE, multi.plot=FALSE)
    }else
    {
      setGraphics(overlay.text=FALSE, save.backup=TRUE)
      refresh(sub.plot=FALSE, multi.plot=FALSE)
    }
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  textButton <- ttkcheckbutton(textFrame, variable=textVal, command=onText, 
                               text='Display names of overlaid spectrum on plot')
  
  ### TSB
  suppTextVal <- tclVar(ifelse(globalSettings$overlay.textSuppress, 1, 0))
  onSuppressText <- function()
  {
    if (as.logical(as.integer(tclvalue(suppTextVal))))
    {
      setGraphics(overlay.textSuppress=TRUE, save.backup=TRUE)
      refresh(sub.plot=FALSE, multi.plot=FALSE)
    }else
    {
      setGraphics(overlay.textSuppress=FALSE, save.backup=TRUE)
      refresh(sub.plot=FALSE, multi.plot=FALSE)
    }
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  textSuppressButton <- ttkcheckbutton(textFrame, variable=suppTextVal, command=onSuppressText, 
                                       text='Suppress path of overlaid spectrum on plot')	
  ### end TSB
  
  
  ##add widgets to fileFrame
  tkgrid(olFileFrame, column=1, row=1, rowspan=3, sticky='nswe', pady=c(5, 8), 
         padx=6)
  tkgrid(olFileBox, column=1, row=1, sticky='nswe')
  tkgrid(olYscr, column=2, row=1, sticky='ns')
  tkgrid(olXscr, column=1, row=2, sticky='we')
  
  ##make fileFrame stretch when window is resized
  tkgrid.columnconfigure(olFrame, 1, weight=1)
  tkgrid.rowconfigure(olFrame, 1, weight=1)
  tkgrid.columnconfigure(olFileFrame, 1, weight=1)
  tkgrid.rowconfigure(olFileFrame, 1, weight=1)
  
  ##add widgets to middleFrame
  tkgrid(middleFrame, column=2, row=1, rowspan=2, sticky='nswe', pady=5, padx=2)
  tkgrid(buttonFrame, row=1, pady=c(28, 0))
  tkgrid(addButton, row=1)
  tkgrid(removeButton, row=2, pady=4)
  
  tkgrid(offsetFrame, row=2, pady=c(10, 0))
  tkgrid(offsetLab, column=1, row=1, sticky='e', padx=c(1, 0))
  tkgrid(valLab, column=2, row=1, sticky='w')
  tkgrid(offsetSlider, column=1, row=2, columnspan=2, padx=c(0, 5))
  tkgrid.rowconfigure(middleFrame, 2, weight=1)
  
  ##add widgets to overlayFrame
  tkgrid(overlayFrame, column=3, row=1, sticky='nswe', pady=5, padx=6)
  tkgrid(overlayBox, column=1, row=1, sticky='nswe')
  tkgrid(overlayYscr, column=2, row=1, sticky='ns')
  tkgrid(overlayXscr, column=1, row=2, sticky='we')
  
  ##make overlayFrame stretch when window is resized
  tkgrid.columnconfigure(olFrame, 3, weight=1)
  tkgrid.columnconfigure(overlayFrame, 1, weight=1)
  tkgrid.rowconfigure(overlayFrame, 1, weight=1)
  
  ##add widgets to colFrame
  tkgrid(colFrame, column=3, row=2, sticky='we', padx=20)
  tkgrid(projColButton, column=2, row=3, pady=c(2, 4))
  tkgrid.columnconfigure(colFrame, 1, weight=1)
  tkgrid.columnconfigure(colFrame, 3, weight=1)
  
  ##add widgets to textFrame
  tkgrid(textFrame, column=2, columnspan=2, row=3, sticky='we')
  tkgrid(textButton, sticky='w', padx=6, pady=6)
  tkgrid(textSuppressButton, sticky='w', padx=6, pady=6)
  
  ##reconfigures widgets in GUI according to which spectra are open
  olConfigGui <- function()
  {
    if (length(as.integer(tkcurselection(olFileBox))))
      tkconfigure(addButton, state='normal')
    else
      tkconfigure(addButton, state='disabled')
    
    configList <- list(offsetSlider, offsetLab, valLab)
    
    if (length(overlayList))
    {
      for (i in configList)
        tkconfigure(i, state='normal')
      tkconfigure(offsetSlider, fg='black')
    }else
    {
      for (i in configList)
        tkconfigure(i, state='disabled')
      tkconfigure(offsetSlider, fg='grey')
    }
    usrSel <- 1 +	as.integer(tkcurselection(overlayBox))
    usrFiles <- overlayList[usrSel]
    configList <- list(projColButton, removeButton)
    
    if (!length(usrSel))
    {
      for (i in configList)
        tkconfigure(i, state='disabled')
    }else
    {
      for (i in configList)
        tkconfigure(i, state='normal')
      
    }
    
    tclObj(offsetVal) <- globalSettings$offset
  }
  tkbind(overlayBox, '<<ListboxSelect>>', olConfigGui)
  olConfigGui()	
  
  ##resets widgets whenever the mouse enters the GUI
  olMouse <- function(){
    reset(list(olFileList, overlaysList), list(olFileBox, overlayBox), 
          list(olFileNames, overlayNames), c('files', 'overlays'))
    olFileNames <<- names(fileFolder)
    overlayMatches <- match(overlayList, olFileNames)
    if (length(overlayMatches))
      olFileNames <<- olFileNames[-overlayMatches]
    if (is.null(overlayList))
      overlayNames <<- character(0)
    else
      overlayNames <<- overlayList
    olConfigGui()
  }
  tkbind(olFrame, '<Enter>', olMouse)
  tkbind(olFrame, '<FocusIn>', olMouse)
  
  ##Allows users to press the 'Enter' key to make selections
  onEnter <- function(){
    focus <- as.character(tkfocus())
    if (length(grep('.1.1.1.1$', focus)))
      olDouble(olFileBox)
    else if (length(grep('.1.1.3.1$', focus)))
      olDouble(overlayBox)
    else
      tryCatch(tkinvoke(focus), error=function(er){})
  }
  tkbind(dlg, '<Return>', onEnter) 
  
  invisible()
}

mf <- function()
{
  ##create main window
  current <- wc()
  tclCheck()
  dlg <- myToplevel('mf')
  
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  
  tkwm.title(dlg, 'Manipulate Files')
  
  fileFrame <- ttklabelframe(dlg, text='Files')
  opFrame <- ttklabelframe(dlg, text='Operations')
  
  tkgrid(fileFrame, column=1, row=1, sticky='nsew', padx=c(6, 0), pady=c(6, 0))
  tkgrid(opFrame, column=2, row=1, sticky='nse', padx=c(6, 0), pady=c(6, 0))
  
  tkgrid.columnconfigure(dlg, 1, weight=1)
  tkgrid.rowconfigure(dlg, 1, weight=1)
  tkgrid.columnconfigure(fileFrame, 1, weight=1)
  tkgrid.rowconfigure(fileFrame, 1, weight=1)
  
  fileList <- tclVar()
  fileNames <- names(fileFolder)
  
  tclObj(fileList) <- getTitles(fileNames)
  fileBox <- tklistbox(fileFrame, height=8, width=25, listvariable=fileList, selectmode='extended', active='dotbox', 
                       exportselection=FALSE, bg='white', xscrollcommand=function(...) tkset(fileXscr, ...),	yscrollcommand=function(...) tkset(fileYscr, ...))
  
  fileXscr <- ttkscrollbar(fileFrame, orient='horizontal', command=function(...) tkxview(fileBox, ...))
  fileYscr <- ttkscrollbar(fileFrame, orient='vertical', command=function(...) tkyview(fileBox, ...))
  
  if (length(fileNames) > 2)
  {
    for (i in seq(0, length(fileNames) - 1, 2))
      tkitemconfigure(fileBox, i, background='#ececff')
  }
  
  currMatch <- match(currentSpectrum, fileNames)
  if (!is.na(currMatch))
  {
    tkselection.set(fileBox, currMatch - 1)
    tcl(fileBox, 'see', currMatch - 1)
  }
  
  ##switches spectra on left-mouse double-click
  onDouble <- function(box)
  {
    usrSel <- 1 + as.integer(tkcurselection(box))
    if (length(usrSel))
    {
      usrFile <- fileNames[usrSel]
      
    }else
      usrFile <- NULL
    
    if (!is.null(usrFile) && currentSpectrum != usrFile)
    {
      myAssign('currentSpectrum', usrFile)
      refresh(multi.plot = FALSE)
      configGui()
      tkwm.deiconify(dlg)
      tkfocus(box)
    }
  }
  tkbind(fileBox, '<Double-Button-1>', function(...) onDouble(fileBox))
  
  onMergeRb <- function()
  {
    tkconfigure(factorLabel, state='disabled') # normal / disabled
    tkconfigure(factorEntry, state='disabled')
    tclvalue(factorVal) <- ''
  }
  
  onAddRb <- function()
  {
    tkconfigure(factorLabel, state='normal') # normal / disabled
    tkconfigure(factorEntry, state='normal')
    tclvalue(factorVal) <- 'Operand'
  }
  
  onSubRb <- function()
  {
    tkconfigure(factorLabel, state='normal') # normal / disabled
    tkconfigure(factorEntry, state='normal')
    tclvalue(factorVal) <- '(-) Operand'
  }
  
  onDivRb <- function()
  {
    tkconfigure(factorLabel, state='normal') # normal / disabled
    tkconfigure(factorEntry, state='normal')
    tclvalue(factorVal) <- 'Denominator'
  }
  
  onMultRb <- function()
  {
    tkconfigure(factorLabel, state='normal') # normal / disabled
    tkconfigure(factorEntry, state='normal')
    tclvalue(factorVal) <- 'Operand'
  }
  
  onMeanRb <- function()
  {
    tkconfigure(factorLabel, state='disabled') # normal / disabled
    tkconfigure(factorEntry, state='disabled')
    tclvalue(factorVal) <- ''
  }
  
  opType <- ''
  
  rbMerge <- ttkradiobutton(opFrame, variable=opType, value='merge', text='Merge', command=onMergeRb)
  
  rbAdd <- ttkradiobutton(opFrame, variable=opType, value='add', text='Add', command=onAddRb)
  
  rbSub <- ttkradiobutton(opFrame, variable=opType, value='sub', text='Subtract', command=onSubRb)
  
  rbDiv <- ttkradiobutton(opFrame, variable=opType, value='div', text='Divide', command=onDivRb)
  
  rbMult <- ttkradiobutton(opFrame, variable=opType, value='mult', text='Multiply', command=onMultRb)
  
  rbMean <- ttkradiobutton(opFrame, variable=opType, value='mean', text='Mean', command=onMeanRb)
  
  tclvalue(opType) <- 'add'
  
  factorVal <- tclVar('Operand')
  factorLabel <- ttklabel(opFrame, textvariable=factorVal)
  
  factorEntryVar <- tclVar(0)
  factorEntry <- ttkentry(opFrame, width=12, justify='right', textvariable=factorEntryVar)	
  
  ##create set 1D color button
  onExecute <- function()
  {
    usrSel <- 1 +	as.integer(tkcurselection(fileBox))
    selectedFiles <- length(usrSel)
    myFactorVal <- tclvalue(factorEntryVar)
    
    #		print(paste0('Selected Files: ', selectedFiles))
    #		print(paste0('myFactorVal: ', myFactorVal))
    
    if (!length(usrSel))
      err(paste('You must select one or more files before',
                'pressing the execute button'))
    else
    {
      for(i in usrSel)
        print(names(fileFolder)[i])
      
      switch(tclvalue(opType),
             'add' = {
               if(selectedFiles > 1)		
                 math_op(fnc = "add", fileFolder, usrSel, opp = NULL)
               else
               {
                 if(!is.null(myFactorVal) && myFactorVal != 0)
                   math_op(fnc = "add", fileFolder, usrSel, opp = myFactorVal)
               }
             },
             'sub' = {
               if(selectedFiles > 1)
               {
                 print('Please enter the number of the file you wish to act as the negative operand, in the operand entry field.')
                 
                 if((myFactorVal != 0) && (myFactorVal %in% usrSel))
                 {
                   sub <- myFactorVal
                   usrSel <- which(usrSel != sub)
                   math_op(fnc = 'sub', fileFolder, c(usrSel, sub), opp = NULL)
                 }
               }else if(selectedFiles == 1 && !is.null(myFactorVal) && myFactorVal != 0)
               {
                 math_op(fnc = 'sub', fileFolder, usrSel, opp = myFactorVal)
               }	
             },
             'mult' = {
               if(selectedFiles == 1 && !is.null(myFactorVal) && myFactorVal != 0)
               {
                 math_op(fnc = 'mult', fileFolder, usrSel, opp = myFactorVal)
               }
             },
             'div' = {
               if(selectedFiles == 1 && !is.null(myFactorVal) && myFactorVal != 0)
               {
                 math_op(fnc = 'div', fileFolder, usrSel, opp = myFactorVal)
               }						
             },
             'merge' = {
               if(selectedFiles > 1)
               {
                 mergeFiles(fileFolder, usrSel, species)
               }							
             },
             'mean' = {
               if(selectedFiles > 1)
                 math_op(fnc = 'mean', fileFolder, usrSel, opp = NULL)
             }				
      )
    }
    
    tkdestroy(dlg)
  }
  
  executeButton <- ttkbutton(dlg, text='Execute', width=13, command=onExecute)
  
  ##add widgets to fileFrame
  tkgrid(fileBox, column=1, row=1, sticky='nswe')
  tkgrid(fileYscr, column=2, row=1, sticky='ns')
  tkgrid(fileXscr, column=1, row=2, sticky='we')
  
  ##add widgets to colFrame
  tkgrid(rbAdd, column=1, row=3, sticky='nw')
  tkgrid(rbSub, column=2, row=3, sticky='nw', padx=16)
  tkgrid(rbMult, column=1, row=4, sticky='nw')
  tkgrid(rbDiv, column=2, row=4, sticky='nw', padx=16)
  tkgrid(rbMerge, column=1, row=5, sticky='nw')
  tkgrid(rbMean, column=2, row=5, sticky='nw', padx=16)
  tkgrid(factorLabel, column=1, row=6, sticky='nw')
  tkgrid(factorEntry, column=2, row=6, sticky='nw', padx=16)
  tkgrid(executeButton, column=1, columnspan=2, row=2, sticky='se', pady=c(2, 4))
  tkgrid(ttksizegrip(dlg), column=3, row=3, sticky='se')
  
  ##reconfigures widgets in GUI according to which spectra are open
  configGui <- function()
  {		
    if(length(as.integer(tkcurselection(fileBox))) == 1)
    {
      configList <- list(executeButton, rbAdd, rbSub, rbMult, rbDiv)
      for(i in configList)
        tkconfigure(i, state='normal')
      
      if(tclvalue(opType) == 'merge' || tclvalue(opType) == 'mean')
        tclvalue(opType) <- 'add'
      
      tkconfigure(rbMerge, state='disabled')
      tkconfigure(rbMean, state='disabled')
      
    }else if (length(as.integer(tkcurselection(fileBox))) == 2)
    {
      configList <- list(executeButton, rbAdd, rbSub, rbMerge, rbMean)
      for(i in configList)
        tkconfigure(i, state='normal')
      
      if(tclvalue(opType) == 'mult' || tclvalue(opType) == 'div')
        tclvalue(opType) <- 'add'
      
      tkconfigure(rbMult, state='disabled')
      tkconfigure(rbDiv, state='disabled')
      
    }else if(length(as.integer(tkcurselection(fileBox))) > 2)
    {
      configList <- list(executeButton, rbAdd, rbMerge, rbMean)
      for(i in configList)
        tkconfigure(i, state='normal')
      
      configList <- list(rbSub, rbMult, rbDiv)
      for(i in configList)
        tkconfigure(i, state='disabled')
      
      if(tclvalue(opType) == 'add' || tclvalue(opType) == 'mean')
        tclvalue(opType) <- tclvalue(opType) # do nothing
      else
        tclvalue(opType) <- 'add'
      
    }else
    {
      configList <- list(executeButton, rbAdd, rbSub, rbMult, rbDiv, rbMerge, rbMean)
      for(i in configList)
        tkconfigure(i, state='disabled')
    }	
  }
  #	tkbind(targetFrame, '<<ListboxSelect>>', olConfigGui)
  configGui()	
  
  #	##resets widgets whenever the mouse enters the GUI
  ##resets widgets whenever the mouse enters the GUI
  onMouse <- function()
  {
    #		reset(list(fileList), list(fileBox), 
    #			list(fileNames), c('files'))
    #		fileNames <<- names(fileFolder)
    
    configGui()
  }
  tkbind(fileFrame, '<Enter>', onMouse)
  tkbind(fileFrame, '<Leave>', onMouse)
  tkbind(fileFrame, '<FocusIn>', onMouse)
  
  ##Allows users to press the 'Enter' key to make selections
  onEnter <- function()
  {
    focus <- as.character(tkfocus())
    if (length(grep('.1.1.1.1$', focus)))
      onDouble(fileBox)
    
    else
      tryCatch(tkinvoke(focus), error=function(er){})
  }
  tkbind(dlg, '<Return>', onEnter) 
  
  invisible()
}


checkColumnList <- function(sStr)
{
  return(!grepl("[^[:space:]\\,\\:0-9+]", sStr))
}

im <- function()
{
  tclCheck()
  dlg <- myToplevel('im')
  if (is.null(dlg))
    return(invisible())
  
  tkwm.title(dlg, 'Import Maven File')
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  
  mavenLabelFrame <- ttklabelframe(dlg, text='Maven Details')
  
  fileLabelFrame <- ttklabelframe(mavenLabelFrame, text='Use the \'Browse\' button to select a Maven file to import.')
  
  fileNameLabel <- ttklabel(fileLabelFrame, text='File Name:')
  
  fnOutput <- tclVar(character(0))
  fileNameOutputLabel <- ttklabel(fileLabelFrame, width=32, justify='left', textvariable=fnOutput)
  
  
  columnFrame <- ttklabelframe(mavenLabelFrame, text='Use the following text entry boxes to provide corresponding 
														column numbers from the Maven output file.')
  
  descriptionLabel <- ttklabel(columnFrame, text='Separate column numbers with commas, i.e. 1,3,4  or 
													use a colon to include all column numbers between an upper
													and lower index, i.e. 7:12 would indicate columns 7, 8, 9, 10, 11 and 12')	
  
  
  compoundIdxLabel <- ttklabel(columnFrame, text='Compounds:')
  blankIdxLabel <- ttklabel(columnFrame, text='Blanks:')
  sampleIdxLabel <- ttklabel(columnFrame, text='Samples:')
  
  
  thresholdFrame <- ttklabelframe(mavenLabelFrame, text='Adjust this value to determine the threshold between blanks and data readings. 
															Higher values reduce chances of false positives.')
  thresholdMultLabel <- ttklabel(thresholdFrame, text='Threshold Multiplier:')
  
  
  
  ciEntry <- tclVar(character(0))
  compoundIdxEntry <- ttkentry(columnFrame, width=24, justify='left', textvariable=ciEntry)
  
  biEntry <- tclVar(character(0))
  blankIdxEntry <- ttkentry(columnFrame, width=24, justify='left', textvariable=biEntry)
  
  siEntry <- tclVar(character(0))
  sampleIdxEntry <- ttkentry(columnFrame, width=24, justify='left', textvariable=siEntry)
  
  tmEntry <- tclVar(6)
  thresholdMultEntry <- ttkentry(thresholdFrame, width=6, justify='left', textvariable=tmEntry)
  
  
  onImport <- function()
  {
    if(nchar(tclvalue(fnOutput)) == 0)
    {
      print('Please use the \'browse\' button to select a Maven file to import. ')	
    }else if(!file.exists(tclvalue(fnOutput)))
    {
      print('The selected file does not exist. Please use the \'browse\' button to select a new file.')
    }else
    {
      if((nchar(tclvalue(ciEntry)) == 0) || (nchar(tclvalue(biEntry)) == 0) || (nchar(tclvalue(siEntry)) == 0))
      {
        print('Please ensure you have input the column numbers for compounds, blanks and samples.')
      }else
      {
        compoundsIdx <- tclvalue(ciEntry)
        blanksIdx <- tclvalue(biEntry)
        samplesIdx <- tclvalue(siEntry)
        thresholdMult <- tclvalue(tmEntry)
        
        if(checkColumnList(compoundsIdx) && checkColumnList(blanksIdx) && checkColumnList(samplesIdx) && is.numeric(thresholdMult) )
        {
          print(paste0(compoundsIdx, '_', blanksIdx, '_', samplesIdx, '_', thresholdMult))
          parseMaven(sMavenFileName = tclvalue(fnOutput), compoundsIdx = compoundsIdx, vBlanksIdx = blanksIdx, vSamplesIdx = samplesIdx, kThresholdMult = thresholdMult)
          tkdestroy(dlg)			
        }else
        {
          if(!checkColumnList(compoundsIdx))
            print('Please ensure that the compound column list contains only integer values separated by commas or a colon.')
          
          if(!checkColumnList(blanksIdx))
            print('Please ensure that the blanks column list contains only integer values separated by commas or a colon.')
          
          if(!checkColumnList(samplesIdx))
            print('Please ensure that the samples column list contains only integer values separated by commas or a colon.')
          
          if(!is.numeric(thresholdMult))
          {
            print('Please ensure that the threshold multiplier is a numeric value.')
            print(paste0('threshold: ', thresholdMult))
          }
        }			
      }
    }
  }
  importButton <- ttkbutton(mavenLabelFrame, text='Import', width=12, command=onImport)
  
  onBrowse <- function()
  {
    loadFileName <- myOpen(defaultextension = 'csv', filetypes = list('csv' = 'Excel File'),	initialfile = "", title= 'Load')
    
    if(nchar(loadFileName))
      tclvalue(fnOutput) <- loadFileName
  }
  browseButton <- ttkbutton(mavenLabelFrame, text='Browse', width=12, command=onBrowse)
  
  ##add widgets to mavenLabelFrame
  tkgrid(mavenLabelFrame, column=1, row=1, sticky='nswe', pady=6, padx=6)
  tkgrid(fileLabelFrame, column=1, columnspan=2, row=1, sticky='nswe', pady=6, padx=6)
  tkgrid(fileNameLabel, column=1, row=1, pady=6, padx=6, sticky='nw')	
  tkgrid(fileNameOutputLabel, column=2, row=1, pady=6, padx=6, sticky='ne')
  tkgrid(columnFrame, column=1, columnspan=2, sticky='nswe', pady=6, padx=6)
  tkgrid(descriptionLabel, column=1, row=2, columnspan=2, pady=6, padx=6, sticky='we')	
  tkgrid(compoundIdxLabel, column=1, row=3, pady=6, padx=6, sticky='w')
  tkgrid(compoundIdxEntry, column=2, row=3, pady=6, padx=6, sticky='e')	
  tkgrid(blankIdxLabel, column=1, row=4, pady=6, padx=6, sticky='w')
  tkgrid(blankIdxEntry, column=2, row=4, pady=6, padx=6, sticky='e')
  tkgrid(sampleIdxLabel, column=1, row=5, pady=6, padx=6, sticky='w')
  tkgrid(sampleIdxEntry, column=2, row=5, pady=6, padx=6, sticky='e')
  tkgrid(thresholdFrame, column=1, columnspan=2, pady=6, padx=6, sticky='we')
  tkgrid(thresholdMultLabel, column=1, pady=6, padx=6, row=6, sticky='w')
  tkgrid(thresholdMultEntry, column=2, pady=6, padx=6, row=6, sticky='e')
  
  tkgrid(browseButton, column=1, row=7, pady=6, padx=6, sticky='sw')
  tkgrid(importButton, column=2, row=7, pady=6, padx=6, sticky='se')
  
  ##make fileFrame stretch when window is resized
  tkgrid.columnconfigure(dlg, 1, weight=1)
  tkgrid.rowconfigure(dlg, 1, weight=10)
  tkgrid.columnconfigure(mavenLabelFrame, 1, weight=1)
  tkgrid.rowconfigure(mavenLabelFrame, 1, weight=1)
  tkgrid.columnconfigure(fileLabelFrame, 1, weight=1)
  tkgrid.rowconfigure(fileLabelFrame, 1, weight=1)
  tkgrid.columnconfigure(columnFrame, 1, weight=1)
  tkgrid.rowconfigure(columnFrame, 1, weight=1)
  tkgrid.columnconfigure(thresholdFrame, 1, weight=1)
  tkgrid.rowconfigure(thresholdFrame, 1, weight=1)
  
  tkgrid(ttksizegrip(dlg), column=2, row=3, sticky='se')
  
  invisible()
}

exportFileListMascotOutput <- function()
{
  ## get user input on the file(s) for processing
  dir <- choose.dir()
  
  if (dir.exists(dir))
  {
    lFileList <- vector(mode="character")
    subDirs2 <- vector(mode="character")
    counter <- 1
    
    subDirs <- list.dirs(path = file.path(dir), full.names = TRUE, recursive = FALSE)
    #		subDirs2 <- vector(mode="character", length=length(subDirs))
    if (length(subDirs) == 0)
      subDirs <- dir
    else
    {
      counter2 <- 1
      for (i in 1:length(subDirs))
      {
        tmp <- paste0(subDirs[i], '\\Mascot Output')
        if(dir.exists(tmp))
        {
          subDirs2[counter2] <- tmp
          counter2 <- counter2 + 1
        }
      }
    }
    
    for (i in 1:length(subDirs2))
    {
      lFilesCSV <- list.files(path = subDirs2[i], pattern = '*.csv', all.files=TRUE, full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
      
      if (length(lFilesCSV) > 0) ## only process if there are files selected
      {
        for (j in 1:length(lFilesCSV))
        {
          lFileList[counter] <- lFilesCSV[j]
          counter <- counter + 1
          
        }
      }
      
    }
    return(lFileList)
  }else
  {
    return(NULL)			
  }
}

getFileListCSV <- function()
{
  ## get user input on the file(s) for processing
  dir <- choose.dir()
  
  if (dir.exists(dir))
  {
    lFileList <- vector(mode="character")
    counter <- 1
    
    subDirs <- list.dirs(path = file.path(dir), full.names = TRUE, recursive = FALSE)
    if (length(subDirs) == 0)
      subDirs <- dir
    
    for (i in 1:length(subDirs))
    {
      lFilesCSV <- list.files(path = subDirs[i], pattern = '*.csv', all.files=TRUE, full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
      
      if (length(lFilesCSV) > 0) ## only process if there are files selected
      {
        for (j in 1:length(lFilesCSV))
        {
          lFileList[counter] <- lFilesCSV[j]
          counter <- counter + 1
          
        }
      }
      
    }
    
    return(lFileList)
  }else
  {
    return(NULL)
  }
}



exportFileListRAW <- function()
{
  ## get user input on the file(s) for processing
  dir <- choose.dir()
  
  if (dir.exists(dir))
  {
    lFileList <- vector(mode="character")
    counter <- 1
    
    subDirs <- list.dirs(path = file.path(dir), full.names = TRUE, recursive = FALSE)
    if (length(subDirs) == 0)
      subDirs <- dir
    
    for (i in 1:length(subDirs))
    {
      lFilesCSV <- list.files(path = subDirs[i], pattern = '*.RAW', all.files=TRUE, full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
      
      if (length(lFilesCSV) > 0) ## only process if there are files selected
      {
        for (j in 1:length(lFilesCSV))
        {
          lFileList[counter] <- lFilesCSV[j]
          counter <- counter + 1
          
        }
      }
      
    }
  }
}

################################################################################
#
# Diana Processing
#
################################################################################

falcilysinA <- c(41, 74, 90)
falcilysinB <- c(40, 73, 78, 132)
plasmepsin1A <- c(33, 46, 98)
plasmepsin1B <- c(31, 41, 129)
plasmepsin2A <- c(33, 108, 136)
plasmepsin2B <- c(32)
falcipainA 	<- c(31, 33)
falcipainB	<- c(32, 69, 82)

# AH1  	A  	33/34  	RMF/LSF
#		  	46/47  	PHF/DLS
#		  	98/99	VNF/KLL
#
#	  	B	31/32  	GRL/LVV
#		  	41/42  	QRF/FES
#		  	129/130	QAA/YQK
#
# AH2	A	33/34	RMF/LSF
#			108/109	LVT/LAA
#			136/137 TVL/TSK
#		B	32/33	RLL/VVY
#
# Cysteine
#		A	31/32	LER/MFL
#			33/34	RMF/LSF
#		B	32/33	RLL/VVY
#			69/70	VLG/AFS
#			82/83	LSA/LSD

plotCutSites <- function()
{
  if(globalSettings$geneDisp == 'HBA1')
  {
    plotCutSite('RMFLSF', 'black')
    plotCutSite('PHFDLS', 'black')
    plotCutSite('VNFKLL', 'black')
    plotCutSite('RMFLSF', 'black')
    plotCutSite('LVTLAA', 'black')
    plotCutSite('TVLTSK', 'black')
    plotCutSite('LERMFL', 'black')
    plotCutSite('RLLVVY', 'black')
    
  }else if(globalSettings$geneDisp == 'HBB')
  {
    plotCutSite('GRLLVV', 'black')
    plotCutSite('QRFFES', 'black')
    plotCutSite('QAAYQK', 'black')
    plotCutSite('RLLVVY', 'black')
    plotCutSite('RLLVVY', 'black')
    plotCutSite('VLGAFS', 'black')
    plotCutSite('LSALSD', 'black')
  }else
  {
    
  }
  
}


plotCutSite <- function(seq, col)
{
  currGene <- globalSettings$geneDisp
  
  if(!is.null(currGene) && !is.null(seq) && !missing(seq) && !is.null(species))
  {
    idx <- which(species$genes$name == currGene)
    
    if(length(idx) > 0)
    {
      start <- species$genes$seqStartIdx[idx]
      end <- start + species$genes$seqLength[idx]
      currSeq <- substr(species$seq, start, end)
      loc <- gregexpr(seq, currSeq,fixed=TRUE)
      
      if(loc[[1]][1] != -1)
      {
        numMatches <- length(loc[[1]])
        
        for(i in 1:numMatches)
        {
          abline(v=loc[[1]][i]+3, col = col, lty=2)
        }
        msg <- 'Sequence(s) found.'
      }else
      {
        msg <- 'Sequence not found.'
      }
      
    }else
    {
      msg <- 'Current gene not found.'	
    }
    
  }else
  {
    if(is.null(currGene))
      msg <- 'Current gene is NULL. '
    
    if(is.null(seq) || empty(seq))
      msg <- paste0(msg, 'Provided message is NULL or empty. ')
    
    if(is.null(species))
      msg <- paste0(msg, 'Object \' species\' has not been initialized.')
  }
  print(msg)
  flush.console()
}


##############################################################################
# parseMascot()
# 
# Author: Travis
###############################################################################
parseMascot <- function(sMascotFileName = "")
{
  ## Read Mascot file
  if (sMascotFileName == "")
  {
    sMascotFileName <- myOpen("csv", list(csv = "Comma Separated Values File, xls = Excel File"), multiple = FALSE)
    dfMascot <- read.csv( sMascotFileName, head = TRUE, stringsAsFactors = FALSE, skip = 3)
  }else
    dfMascot <- read.csv( sMascotFileName, head = TRUE, stringsAsFactors = FALSE, skip = 3)
  
  
  if(length(dfMascot$pep_seq) == 0)
  {
    return(NULL)
  }
  
  ###
  bCols = NULL
  
  for(i in 1:length(dfMascot) )
  {
    bCols <- c(bCols, all(is.na(dfMascot[, i])))
  }
  pepSeqIdx <- which(names(dfMascot)=="pep_seq")
  pepScoreIdx <- which(names(dfMascot)=="pep_score")
  
  ## when the first row of dfMascot$prot_acc is filled with the indicated string, that indicates
  ## that there were no proteins matched to the peptides by Mascot --- due to a bug in the Mascot
  ## software, the peptide data is written to the data file in an offset manner, depending on the number
  ## of columns written before pep_query. Therefore, check which column "pep_query" falls in, and 
  ## subtract that many columns to get the new location of the desired data
  if(dfMascot$prot_acc[1] == "--------------------------------------------------------")
  {
    offset <- which(names(dfMascot)=="pep_query") - 1
    pepSeqIdx <- pepSeqIdx - offset
    pepScoreIdx <- pepScoreIdx - offset
  }
  ###
  
  sGenotype <- getGenotype(sMascotFileName)
  
  dfNew <- subset(dfMascot, select = c(pepSeqIdx, pepScoreIdx))
  
  names(dfNew)[1] = "sequence"
  names(dfNew)[2] = "score"	
  
  ## check for any records that have empty sequences
  indices <- which(dfNew$sequence == '')
  ## if there are records with empty peptide sequences -> remove them
  if (length(indices) > 0)
    dfNew <- dfNew[-indices,]
  
  indices <- which(!duplicated(dfNew$sequence))
  
  pep <- cbind(ID = 1:length(indices), dfNew[indices,], stringsAsFactors = FALSE)
  pep <- pep[order(pep$sequence),]	
  #dfMascot <- dfMascot[order(dfMascot$peptide),]
  #	return(preparePeptides(dfMascot, sMascotFileName))
  
  return(list("genotype"=sGenotype, "pep" = pep, "fileName" = sMascotFileName))
}


## DO NOT DELETE
## this function allows for the potential use of "frequency" of peptides
#preparePeptides <- function(dfMascot, fileName)
#{
#	sGenotype <- dfMascot$genotype[1]
#	vInit <- vector(mode = "integer", length = length(dfMascot$peptide))
#	
#	pepScore <- dfMascot[!duplicated(dfMascot$peptide),]$score
#	
#	dfPeptides<- as.data.frame(table(dfMascot$peptide), stringsAsFactors = FALSE)	
#	
#	names(dfPeptides)[1] = "sequence"
##	names(dfPeptides)[2] = "frequency"
#	
##	dfPeptides <- cbind(dfPeptides, ID = (1:length(unique(dfMascot$peptide))))
##	dfPeptides <- cbind(dfPeptides, score = score)
#	dfPeptides <- data.frame("ID" = (1:length(unique(dfMascot$peptide))), "sequence" = dfPeptides$sequence, "score" = pepScore, stringsAsFactors=FALSE)
#	
#	result <- list("genotype" = sGenotype, "pep" = dfPeptides, "fileName" = fileName)
#	
#	return(result)
#}

prepareSpecies <- function(name, ID, genes)
{
  genes <- genes[order(genes$start),]
  genes <- genes[order(genes$chrom),]
  genes <- genes[order(suppressWarnings(as.numeric(genes$chrom))),]
  genes$seqStartIdx <- as.integer(cumsum(nchar(genes$seq)) - nchar(genes$seq) + 1)
  genes$seqLength <- as.integer(nchar(genes$seq))
  
  lSpecies <- list()
  lSpecies$name <- name
  lSpecies$ID <- ID
  lSpecies$genes <- genes
  lSpecies$seq <- paste(genes$seq, collapse = '', sep = '')
  lSpecies$genes["seq"] <- NULL
  return(lSpecies)
}
# add gc()
findAlignments <- function (pep, targetSeq)
{
  totPeps <- length(pep$ID)
  indices <- vector(mode = "integer")
  lMatch <- list(indices = indices, pepLength = 0, score = 0)
  
  lMatches <- vector(mode = "list", length = totPeps)
  nonMatchIndices <- vector(mode = "integer", length=0)
  
  numMatches = 0
  
  for (i in 1:totPeps)
  {
    matches <- gregexpr(pep$sequence[i], targetSeq, fixed = TRUE) ## fixed = TRUE for performance
    if (matches[[1]][1] != -1)	# make sure valid index found by checking first returned index
    {						# make sure it is not -1 (not found) and a valid index
      
      numMatches <- numMatches + 1
      lMatch$indices <- matches[[1]][1:length(matches[[1]])]
      lMatch$pepLength <- attr(matches[[1]], "match.length")[1]
      lMatch$score <- pep$score[i]
      lMatches[[numMatches]] <- lMatch
    }
  }
  
  if(numMatches > 0)
  {
    lMatches <- lMatches[1:numMatches]	#remove any excess members of list that weren't found
    return(lMatches)
  }else
  {
    return(NULL)
  }
}



##	Searches each peptide sequence against the proteome of the species of interest
## the starting index of each location where the entire peptide matches the proteome is recorded
## along with the mascot score of that peptide and the length of the peptide (for later use)
#
# findAlignments <- function (pep, targetSeq)
# {
# 	totPeps <- length(pep$ID)
# 	indices <- vector(mode = "integer")
# 	lMatch <- list(indices = indices, pepLength = 0, score = 0)
# 	
# 	lMatches <- vector(mode = "list", length = totPeps)
# 	nonMatchIndices <- vector(mode = "integer", length=0)
# 	
# 	numMatches <- 0
# 
# 	for (i in 1:totPeps)
# 	{	
# 		matches <- cpp_gregexpr2(targetSeq, pep$sequence[i])
# 		if (length(matches) > 0)	
# 		{						
# 			numMatches <- numMatches + 1
# 			lMatch$indices <- matches
# 			lMatch$pepLength <- nchar(pep$sequence[i])
# 			lMatch$score <- pep$score[i]
# 			lMatches[[numMatches]] <- lMatch
# 		}
# 	}
# 	
# 	if(numMatches > 0)
# 	{
# 		lMatches <- lMatches[1:numMatches]	#remove any excess members of list that weren't found
# 		return(lMatches)
# 	}else
# 	{
# 		return(NULL)
# 	}
# }

##############################################################################
# getGenotype()
# Extracts genotype from filename
# Author: Travis
###############################################################################
getGenotype <- function(sPath)
{
  vIndices <- gregexpr('__', sPath, fixed = TRUE)
  firstPeriod <- gregexpr('.', sPath, fixed = TRUE)
  
  if((vIndices[[1]][1] != -1)&&(firstPeriod[[1]][1]!= -1))
  {
    sGeno <- substr(sPath, vIndices[[1]][1]+3, firstPeriod[[1]][1]-1)
    
    vIndices <- gregexpr('_', sGeno, fixed = TRUE)
    
    if(vIndices[[1]][1] != -1)
    {
      sGeno <- substr(sGeno, 1, vIndices[[1]][1]-1)
    }
    
  }else
  {
    sGeno <- "Not Available"	
  }
  
  return(sGeno)
}

mapToGene <- function(lSpecies, aaMap)
{
  geneMap <- vector(mode = "numeric", length = length(lSpecies$genes$name))
  seqLen <- lSpecies$genes$seqLength
  geneStart <- lSpecies$genes$seqStartIdx
  geneEnd <- lSpecies$genes$seqStartIdx + lSpecies$genes$seqLength - 1
  
  for (iGene in 1:length(lSpecies$genes$name))
  {
    geneMap[iGene] <- sum(aaMap[geneStart[iGene]:geneEnd[iGene]])/seqLen[iGene]		
  }	
  
  mapData <- 	list(geneMap, geneStart, geneEnd, aaMap) ### why are we passing back gene start and end -> is this necessary? is this for plotting?
  names(mapData) <- c('geneMap', 'geneStart', 'geneEnd', 'aminoAcidMap')
  
  return(mapData)
}

mapToAA_old <- function(lSpecies, lMatches)
{
  aaMap <- vector(mode = "integer", length = nchar(lSpecies$seq))
  
  for (i in 1:length(lMatches))
  {
    indices <- lMatches[[i]]$indices
    
    for (j in 1:lMatches[[i]]$pepLength)
    {
      aaMap[indices] <- aaMap[indices] + 1  ## if we want to use frequency in the future... add frequency instead of "1" here
      indices <- indices + 1
    }
  }
  
  idx <- which(aaMap > 0)
  val <- as.integer(aaMap[idx])
  
  mat <- matrix(c(idx, val), nrow=2, byrow=TRUE)
  
  return(mapToGene(lSpecies, mat))
}

mapToAAPeptides <- function(lMatches)
{
  indexCount <- sapply(lMatches, function(x){length(x$indices)})
  indicesTot <- sum(indexCount)	
  
  vecStart <- vector(mode = "integer", length = indicesTot)
  vecEnd <- vector(mode = "integer", length = indicesTot)
  index = 1
  
  for (i in 1:length(lMatches))
  {
    indexStarts <- lMatches[[i]]$indices
    len <- lMatches[[i]]$pepLength - 1
    
    for (j in 1:length(indexStarts))
    {
      vecStart[index] <- indexStarts[j]
      vecEnd[index] <- (indexStarts[j] + len)
      index <- index + 1
    }
  }
  
  return(c(vecStart, vecEnd))
}


mapToAAEndPoints <- function(lSpecies, lMatches)
{
  aaMap <- vector(mode = "integer", length = nchar(lSpecies$seq))
  
  for (i in 1:length(lMatches))
  {
    indexStarts <- lMatches[[i]]$indices
    len <- lMatches[[i]]$pepLength - 1
    
    for (j in 1:length(indexStarts))
    {
      indices <- c(indexStarts[j], (indexStarts[j] + len))
      aaMap[indices] <- aaMap[indices] + 1  ## if we want to use frequency in the future... add frequency instead of "1" here
    }
  }
  
  return(aaMap)
}


mapToAA <- function(lSpecies, lMatches)
{
  aaMap <- vector(mode = "integer", length = nchar(lSpecies$seq))
  
  for (i in 1:length(lMatches))
  {
    indexStarts <- lMatches[[i]]$indices		
    len <- lMatches[[i]]$pepLength - 1
    
    for (j in 1:length(indexStarts))
    {
      indices <- seq.int(from = indexStarts[j], to = indexStarts[j] + len, by = 1)
      aaMap[indices] <- aaMap[indices] + 1  ## if we want to use frequency in the future... add frequency instead of "1" here
    }
  }
  
  return(mapToGene(lSpecies, aaMap))
}

processFile <- function(fileName, species)
{
  msg <- paste0('Parsing Mascot file: ', fileName)
  print(msg)
  flush.console()
  resList <- parseMascot(fileName)
  
  if(is.null(resList))
  {
    print('No records found.')
    return(NULL)
  }else
  {
    res <- alignAndScorePeptides(resList, species)
    
    if(is.numeric(res) && res == -1)
    {
      return(-1)
    }else if(res$Gene == -1 && res$Matches == -1) ## no matches found
    {
      return(-2)
    }else
    {
      return(res)
    }
  }
}

alignAndScorePeptides <- function(resList, lSpecies)
{
  #Rprof("myProfile.txt")
  appended_vectors <- 0
  
  if(is.null(resList))
  {
    return(-1)
  }
  
  dfPeptides <- resList$pep
  sGenotype <- resList$genotype
  fileName <- resList$fileName
  
  if(is.null(lSpecies))
    lSpecies <- get("species", envir=.GlobalEnv)
  ### generate a decoy species proteome sequence to perform an alignment against for false discovery rate thresholding
  
  lDecoy <- lSpecies
  lDecoy$name <- "Decoy"
  lDecoy$ID <- lDecoy$ID + 1
  
  tmp <- lSpecies$seq
  tmp1 <- unlist(strsplit(tmp, "", fixed = TRUE)) ## convert string to vector
  tmp2 <- sample(tmp1) ## randomly scramble vector of characters
  lDecoy$seq <- paste(tmp2, collapse = '') ## convert vector back to string 
  
  ## perform real alignment	
  msg <- paste0('Aligning ', length(dfPeptides$ID), ' peptides along genome. This may take a few minutes.')
  print(msg)
  flush.console()	
  lMatches <- findAlignments(dfPeptides, lSpecies$seq)
  ## map matches
  if(!is.null(lMatches)) 
  {
    msg <- paste0('Generating coincidence map.')
    print(msg)
    flush.console()		
    mapData <- mapToAA(lSpecies, lMatches)
    vEndPoints <- mapToAAEndPoints(lSpecies, lMatches)
    vPepPoints <- mapToAAPeptides(lMatches)
    aminoAcidMap <- mapData$aminoAcidMap
    geneMap <- mapData$geneMap
  }else
  {
    vPepPoints <- NULL
    vEndPoints <- NULL
    aminoAcidMap <- NULL
    geneMap <- NULL
  }
  
  ## perform alignment against decoy proteome
  msg <- paste0('Aligning ', length(dfPeptides$ID), ' peptides along decoy genome for false discovery rate estimation. This may take a few minutes.')
  print(msg)
  flush.console()	
  lDecoyMatches <- findAlignments(dfPeptides, lDecoy$seq)
  ## map decoy matches	
  if(!is.null(lDecoyMatches)) 
  {
    mapDecoy <- mapToAA(lDecoy, lDecoyMatches)
    genomeSummary <- summary(mapDecoy$aminoAcidMap)
    geneSummary <- summary(mapDecoy$geneMap)
    genomeNoiseSD <- sd(mapDecoy$aminoAcidMap)
    geneNoiseSD <- sd(mapDecoy$geneMap)
    
  }else
  {
    genomeSummary <- -1
    genomeNoiseSD <- -1
    geneSummary <- -1
    geneNoiseSD <- -1
  }
  
  result <- vector(mode="list", length=12)
  
  names(result) <- c("Genotype", "File", "Matches", "Genome", "Gene", "GenomeSummary", 
                     "GenomeNoiseSD", "GeneSummary", "GeneNoiseSD", "appended_vectors", "lVectors", "sVecLabel")	
  
  result$Genotype <- sGenotype
  result$File <- substr(fileName, gregexpr('/', fileName)[[1]][length(gregexpr('/', fileName)[[1]])] + 1, nchar(fileName))
  result$GenomeSummary <- genomeSummary
  result$GenomeNoiseSD <- genomeNoiseSD
  result$GeneSummary <- geneSummary
  result$GeneNoiseSD <- geneNoiseSD
  
  if(is.null(lMatches))
    result$Matches <- -1
  else
    result$Matches <- lMatches
  
  if(is.null(aminoAcidMap))
    result$Genome <- -1
  else
    result$Genome <- aminoAcidMap
  
  if(is.null(geneMap))
    result$Gene <- -1
  else
    result$Gene <- geneMap
  
  
  if(!is.null(vEndPoints))
  {
    appended_vectors <- appended_vectors + 1
    result$appended_vectors <- appended_vectors
    result$lVectors[[result$appended_vectors]] <- vEndPoints
    result$sVecLabel[appended_vectors] <- globalSettings$vectorType$EndPoints		
  }
  
  if(!is.null(vPepPoints))
  {
    appended_vectors <- appended_vectors + 1
    result$appended_vectors <- appended_vectors
    result$lVectors[[result$appended_vectors]] <- vPepPoints
    result$sVecLabel[appended_vectors] <- globalSettings$vectorType$PepPoints
  }
  
  return(result)
}

math_op <- function(fnc = "add", fileFolder, ffIndices, opp = NULL, saveFileName = '')
{
  if((!is.null(fileFolder)) && (exists("species")))
  {
    numIndices <- length(ffIndices)
    vVals <- NULL
    sVecLabel <- NULL
    
    switch(fnc,
           'add' = ## add aa vectors, re-map gene map and combining match lists / any number of files 
             {
               if(numIndices == 1)
               {
                 aaResult <- vector(mode="integer", length=nchar(species$seq))	
                 par <- fileFolder[[ffIndices[1]]]$file.par
                 aa1 <- retrieveAAVector(par$file.name, par$endian, par$matrix_location)
                 aaResult <- aa1 + opp
                 noise_adj <- opp
                 
               }else if(numIndices > 1)
               {			
                 #							# load all aaMaps from indicated files
                 aaResult <- vector(mode="integer", length=nchar(species$seq))
                 for(i in 1:numIndices)
                 {
                   par <- fileFolder[[ffIndices[i]]]$file.par
                   aaResult <- aaResult + retrieveAAVector(par$file.name, par$endian, par$matrix_location)
                 }
                 appended_vectors <- 0
               }else
               {
                 print("Inappropriate number of files selected for addition operation")
                 return(-1)							
               }
             },
           'sub' = ## subtract: requires 2 files.. 1 aa vector from another, re-map gene map. Common matches can be eliminated... but non-common matches
             ## are problematic?? - negative coincidence?
             {
               if(numIndices == 2)
               {
                 aaResult <- vector(mode="integer", length=nchar(species$seq))
                 par <- fileFolder[[ffIndices[1]]]$file.par
                 aa1 <- retrieveAAVector(par$file.name, par$endian, par$matrix_location)
                 par <- fileFolder[[ffIndices[2]]]$file.par
                 aa2 <- retrieveAAVector(par$file.name, par$endian, par$matrix_location)
                 aaResult <- aa1 - aa2
               }else if(numIndices == 1 && !is.null(opp) && (is.numeric(opp) || is.integer(opp)))
               {
                 aaResult <- vector(mode="integer", length=nchar(species$seq))	
                 par <- fileFolder[[ffIndices[1]]]$file.par
                 aa1 <- retrieveAAVector(par$file.name, par$endian, par$matrix_location)
                 aaResult <- aa1 - opp
                 noise_adj <- opp
               }else
               {
                 print("Inappropriate number of files selected for subtraction operation")
                 return(-1)
               }
               
             },
           'mult'= ## single file? multiply aa vector by a scalar... re-map gene map... leave match list alone?
             {
               if(numIndices == 1 && !is.null(opp) && (is.numeric(opp) || is.integer(opp)))
               {
                 aaResult <- vector(mode="integer", length=nchar(species$seq))	
                 par <- fileFolder[[ffIndices[1]]]$file.par
                 aa1 <- retrieveAAVector(par$file.name, par$endian, par$matrix_location)
                 aaResult <- aa1 * opp							
                 noise_adj <- opp
               }else
               {
                 print("Inappropriate number of files selected for multiplication operation")
                 return(-1)
               }
             },
           'div' = ## single file? divide aa vector by a scalar... re-map gene map... leave match list alone?
             {
               if(numIndices == 1 && !is.null(opp) && (is.numeric(opp) || is.integer(opp)))
               {
                 aaResult <- vector(mode="integer", length=nchar(species$seq))	
                 par <- fileFolder[[ffIndices[1]]]$file.par
                 aa1 <- retrieveAAVector(par$file.name, par$endian, par$matrix_location)
                 aaResult <- aa1 / opp
                 noise_adj <- opp
               }else
               {
                 print("Inappropriate number of files selected for division operation")
                 return(-1)												
               }
               
             },
           'mean' =
             {
               if(numIndices > 1)
               {
                 # mResults + 1 to hold the mean values
                 # mResults + 2 to hold the standard deviation values							
                 mResults <- matrix(nrow=numIndices+2, ncol=nchar(species$seq), byrow=TRUE)
                 mResults[numIndices+1,] <- 0
                 mResults[numIndices+2,] <- 0
                 
                 for(i in 1:numIndices)
                 {
                   par <- fileFolder[[ffIndices[i]]]$file.par
                   mResults[i,] <- retrieveAAVector(par$file.name, par$endian, par$matrix_location)
                 }
                 # sum the previous rows
                 for(i in 1:numIndices)
                 {
                   mResults[numIndices+1,] <- mResults[numIndices+1,] + mResults[i,]
                 }
                 
                 # divide by numIndices to get the mean in row numIndices + 1
                 mResults[numIndices+1,] <- mResults[numIndices+1,] / numIndices
                 
                 # subtract mean from each of the the original values
                 for(i in 1:numIndices)
                 {
                   mResults[i,] <- mResults[i,] - mResults[numIndices+1,] 
                 }
                 
                 # square each of the resulting values
                 for(i in 1:numIndices)
                 {
                   mResults[i,] <- mResults[i,] * mResults[i,] 
                 }
                 
                 # sum the squares
                 for(i in 1:numIndices)
                 {
                   mResults[numIndices+2,] <- mResults[numIndices+2,] + mResults[i,] 
                 }							
                 # divide by the number of indices to get the mean
                 mResults[numIndices+2,] <- mResults[numIndices+2,] / numIndices
                 #take the sqrt
                 mResults[numIndices+2,] <- sqrt(mResults[numIndices+2,])
                 
                 # mResults + 1 holds the mean
                 # mResults + 2 holds the standard deviation values
                 aaResult <- mResults[numIndices+1,] ## mean at the aa level
                 vVals <- mResults[numIndices+2,]  ## sd at the aa level
                 sVecLabel <- globalSettings$vectorType$Mean
               }
             }
    )
    
    
    if(!is.null(aaResult))
    {
      mapData <- mapToGene(species, aaResult)
      new.file <- fileFolder[[ffIndices[1]]]
      
      ## after performing vector operations, any previously calculated 'bonus' vectors are invalid... cancel them out
      new.file$file.par$appended_vectors <- 0
      
      if(saveFileName == '')
        saveFileName <- mySave(defaultextension = '.dcf', filetypes = list('dcf' = 'DIANA Compressed File'),	initialfile = '', title= 'Save')
      
      if(is.null(saveFileName) || (nchar(saveFileName) <= 0))
      {
        return(-1)
      }
      
      new.file$file.par$file.name <- saveFileName
      new.file$file.par$user_title <- saveFileName		
      new.file$file.par$force_reload <- TRUE
      
      if(!is.null(vVals))
      {
        lVectors <- vVals
      }else
      {
        lVectors <- NULL
      }
      
      saveAs(saveFileName, in.folder = new.file, mapData, includeMatches = FALSE, lVectors, sVecLabel)
      
      fileFolder[[(length(fileFolder)+1)]] <- new.file
      names(fileFolder)[length(fileFolder)] <- new.file$file.par$file.name
      currentSpectrum <- new.file$file.par$file.name
      
      myAssign("fileFolder", fileFolder, save.backup = FALSE)
      myAssign("currentSpectrum", currentSpectrum, save.backup = FALSE)
      refresh()
      
    }
  }
}

mergeFiles <- function(fileFolder, ffIndices, species, name = '')
{
  
  resList <- mergeResultList(fileFolder, ffIndices, species)
  result <- alignAndScorePeptides(resList, species)
  
  sName <- writeDIANA(name, species, result) 
  print(paste0(sName, ' saved at ', format(Sys.time(), "%H:%M"), '.'))	
  
}


getStats <- function(myVector, noise_sd)
{
  statSummary <- myVector
  
  
  # if statSummary == -1 there were insufficient matches found to determine FDR
  if(as.vector(statSummary)[1] == -1)
  {
    ## use 6 * standard deviation of decoy data ( ~= noise) + mean of noise as our noise threshold
    noiseEst <- as.numeric(-1)
    
    noiseSD <- as.numeric(noise_sd) # 8 bytes
    noiseMax <- as.numeric(-1) # 8 bytes
    noiseMean <- as.numeric(-1) # 8 bytes
    noiseLowHinge <- as.numeric(-1) # 8 bytes
    noiseUpperHinge <- as.numeric(-1) # 8 bytes
    noiseMin <- as.numeric(-1) # 8 bytes
  }else
  {
    if (noise_sd == -1)
    {
      noiseEst <- as.numeric(-1) # 8 bytes
      noiseSD <- as.numeric(noise_sd) # 8 bytes
    }else
    {
      ## (use standard deviation of decoy data ( ~= noise) * multiplier constant) + mean of noise as our noise threshold
      noiseEst <- as.numeric(statSummary[4] + globalSettings$sd_noise_multiplier * (noise_sd)) # 8 bytes			
      noiseSD <- as.numeric(noise_sd) # 8 bytes
    }
    
    noiseMax <- as.numeric(statSummary[6]) # 8 bytes
    noiseMean <- as.numeric(statSummary[4]) # 8 bytes
    noiseLowHinge <- as.numeric(statSummary[2]) # 8 bytes
    noiseUpperHinge <- as.numeric(statSummary[5]) # 8 bytes
    noiseMin <- as.numeric(statSummary[1]) # 8 bytes
  }
  
  stats <- data.frame(est = noiseEst, sd = noiseSD, max = noiseMax, min = noiseMin, mean = noiseMean, 
                      lHinge = noiseLowHinge, uHinge = noiseUpperHinge, stringsAsFactors = FALSE)
  
  return(stats)
}

getPeptidePoints <- function(file.par, aaIndices)
{
  
  #	file.par <- fileFolder[[1]]$file.par
  #	idx <- which(species$genes$name == 'HBB')
  #	sIdx <- species$genes$seqStartIdx[idx]
  #	eIdx <- species$genes$seqLength[idx] + sIdx -1
  #	aaIndices <- seq.int(from = sIdx, to = eIdx, by = 1)	
  
  vStartResult <- NULL
  vEndResult <- NULL
  
  vecIdx <- which(file.par$vecInfo$label == globalSettings$vectorType$PepPoints)
  
  if(length(vecIdx) > 0)
  {
    aaVec <- retrieveVectorInt(fileName = file.par$file.name, endi = file.par$endian, offset = file.par$vecInfo$loc[vecIdx])
  }else
    return(aaResult)
  
  numIndices <- length(aaVec)
  mid <- numIndices / 2
  
  if(numIndices > 0)
  {
    vecStart <- aaVec[1:mid]
    vecEnd <- aaVec[(mid+1):numIndices]
  }
  
  tmp <- diff(aaIndices)
  idx <- which(tmp > 1)
  
  if(length(idx) > 0) ## if this is true more than 1 set of indices in vector
  {					## the indices holding values greater than 1 indicate non-consecutive indices
    ## since diff returns a vector of -1 in length to vector it was applied to
    ## must add 1 to the idx to get indices of interest
    vStart <- c(1, idx+1)
    vEnd <- c(idx+1, length(aaIndices))
  }else
  {
    vStart <- 1
    vEnd <- length(aaIndices)
  }
  
  # for each 'set' of indices or each gene
  for(i in 1:length(vStart))
  {
    ## need to figure out how to determine if line from start to end vector spans the area in the plot... if so include
    ## otherwise do not pass back for plotting. Rather than finding out which to include... figure out what to exclude:
    ## can exlude: those which have right sides < than left side of plot and those that have left sides > than right side
    
    # iMin contains the minimum index in the current gene
    # iMax contains the maximum index in the current gene
    iMin <- aaIndices[vStart[i]]
    iMax <-aaIndices[vEnd[i]]
    
    tooBig <- which(vecStart > iMax)
    tooSmall <- which(vecEnd < iMin)
    exclude <- c(tooBig, tooSmall)
    
    if (i == 0)
    {
      globalIDs <- vecStart[-exclude]
      vStartResult <- globalIDs - iMin + 1
      
      globalIDs <- vecEnd[-exclude]
      vEndResult <- globalIDs - iMin + 1
    }else
    {
      globalIDs <- vecStart[-exclude]
      vStartResult <- c(vStartResult, globalIDs - iMin + 1)
      globalIDs <- vecEnd[-exclude]
      vEndResult <- c(vEndResult, globalIDs - iMin + 1)
    }
  }
  df <- data.frame('x1' = vStartResult, 'x2' = vEndResult, stringsAsFactors = FALSE)
  df <- df[order(df$x1),]
  
  return(df)	
}
#
#indices <- list()
# for(i in 1:length(vStartResult))
# indices[[i]] <- seq.int(from = vStartResult[i], to = vEndResult[i], by = 1)
#
#peps <- list()
#for(i in 1:length(vStartResult))
#	peps[[i]] <- substr(species$seq, vStartResult[i], vEndResult[i])

getEndPoints <- function (file.par, aaIndices)
{
  aaResult <- 0
  
  vecIdx <- which(file.par$vecInfo$label == globalSettings$vectorType$EndPoints)
  
  if(length(vecIdx) > 0)
  {
    aaVec <- retrieveAAVectorInt(file.par$file.name, file.par$endian, file.par$vecInfo$loc[vecIdx])
    #		mat <- readMatrixInt()
    
    ## aaIndices are the indices of the genes of interest
    aaResult <- aaVec[aaIndices]
  }
  return(aaResult)	
}

getVariancePoints <- function (file.par, aaIndices)
{
  aaResult <- 0
  
  vecIdx <- which(file.par$vecInfo$label == globalSettings$vectorType$Mean)
  
  if(length(vecIdx) > 0)
  {
    aaVec <- retrieveAAVector(file.par$file.name, file.par$endian, file.par$vecInfo$loc[vecIdx])
    
    ## aaIndices are the indices of the genes of interest
    aaInterest <- aaVec[aaIndices]
    
    ## to plot variance, we need a negative and positive version of these values to add to the mean values
    ## to make a polygon, need to have the second set of values in reverse order, need to add the initial value again to complete the path
    aaResult <- c(aaInterest, rev(aaInterest)*-1, aaInterest[1])		
  }
  return(aaResult)	
}


parseMaven <- function(sMavenFileName = '', compoundsIdx = NULL, vBlanksIdx, vSamplesIdx, kThresholdMult = 6)
{
  ## Read Maven file
  if (sMavenFileName == "")
  {
    sMavenFileName <- myOpen("csv", list(csv = "Comma Separated Values File, xls = Excel File"), multiple = FALSE)
    dfMaven <- read.csv( sMavenFileName, head = TRUE, stringsAsFactors = FALSE, skip = 0)
  }else
    dfMaven <- read.csv( sMavenFileName, head = TRUE, stringsAsFactors = FALSE, skip = 0)
  
  ## get rid of crap first row
  dfMaven <- dfMaven[-1,]
  
  if(is.null(compoundsIdx))
  {
    compoundsIdx <- which(names(dfMaven) == 'compound')
  }
  
  if(length(compoundsIdx)>0)
  {
    vCompoundsIdx <- which(dfMaven[compoundsIdx] != '')
    vCompoundNames <- dfMaven[[compoundsIdx]][vCompoundsIdx]
  }
  
  vBlankAvgs <- vector(mode = "numeric", length = length(dfMaven[[compoundsIdx]]))
  # average blank values for each peptide
  for(i in vCompoundsIdx)
  {
    vBlankAvgs[i] <- mean(as.numeric(dfMaven[vBlanksIdx][i,]))
  }
  idx <- which(vBlankAvgs != 0)
  otherIdx <- setdiff(1:length(vBlankAvgs), idx)
  
  vBlankAvgs[otherIdx] <- mean(vBlankAvgs[idx])
  
  vThreshold <- vBlankAvgs * kThresholdMult
  vSamplesNames <- names(dfMaven[vSamplesIdx])
  
  mat <- matrix(data = NA, nrow= length(vCompoundsIdx), ncol=length(vSamplesIdx))
  
  for(i in 1:length(vCompoundsIdx))
  {
    mat[i,] <- 0 < (as.numeric(dfMaven[vSamplesIdx][i,]) - vThreshold[i])
  }
  
  # get target directory for saving files
  dir <- choose.dir()
  
  for(i in 1:length(vSamplesIdx))
  {
    fileName <- vSamplesNames[i]
    peps <- vCompoundNames[which(mat[,i])]
    
    if(length(peps) > 0)
    {
      df <- data.frame()
      
      for(i in 1:length(peps))
      {
        row <- list(prot_hit_num = 1, prot_acc = 'dummy', prot_mass = 1, pep_query = 1, pep_rank=1, pep_isbold=1, pep_isunique=0, pep_exp_mz=1,
                    pep_exp_mr=1, pep_delta=0.1, pep_score=1, pep_res_before='M', pep_seq = peps[i], pep_res_after='M')
        
        df <- rbind(df, row, stringsAsFactors = FALSE)
      }
      sPath <- paste0(dir, '\\', fileName, '.csv')
      
      ## to insert blank lines in csv
      a <-rep(NA,dim(df)[2])
      a <-rbind(a,a,a)
      
      write.table(a, file =sPath, row.names = FALSE, na="", col.names=FALSE, sep=',')
      
      write.table(df, file=sPath, row.names=FALSE, na='', col.names = TRUE, append=TRUE, sep=',')
    }
  }
}

mergeResultList <- function(fileFolder, indices, species)
{
  numIndices <- length(indices)
  
  for(i in 1:numIndices)
  {
    par <- fileFolder[[indices[i]]]$file.par
    print(paste0('Retrieving peptides from: ', par$file.name, '.'))
    flush.console()		
    if(i == 1)
    {
      pepListHolder <- retrievePeptidesList(par$file.name, par$endian, par$match_list_location, species)
      
    }else
    {
      pepListHolder <- append(pepListHolder, retrievePeptidesList(par$file.name, par$endian, par$match_list_location, species))
    }
  }
  
  print('Compiling master peptide list.')
  flush.console()
  sequences <- sapply(pepListHolder, function(x){x$sequence})
  idx <- which(duplicated(sequences))
  
  if(length(idx) > 0)
  {
    uniquePepList <- pepListHolder[-idx]	
  }else
  {
    uniquePepList <- pepListHolder
  }
  
  sequences <- sapply(uniquePepList, function(x){x$sequence})
  uniquePepList <- uniquePepList[order(sequences)]
  scores <- sapply(uniquePepList, function(x){x$score})
  
  pep <- data.frame("ID" = 1:length(uniquePepList), "sequence" = sequences, "score" = scores, stringsAsFactors = FALSE)
  resList <- list()
  resList$genotype <- fileFolder[[indices[1]]]$file.par$genotype
  resList$pep <- pep
  resList$fileName <- ''
  
  ## append fileName and genotype information before returning
  ## make sure format is identical to that returned from parseMascot
  return(resList)	
}

highlightZone <- function(seq, y, col, wid = 1)
{
  currGene <- globalSettings$geneDisp
  
  if(!is.null(currGene) && !is.null(seq) && !missing(seq) && !is.null(species))
  {
    idx <- which(species$genes$name == currGene)
    
    if(length(idx) > 0)
    {
      start <- species$genes$seqStartIdx[idx]
      end <- start + species$genes$seqLength[idx]
      currSeq <- substr(species$seq, start, end)
      loc <- gregexpr(seq, currSeq,fixed=TRUE)
      
      if(loc[[1]][1] != -1)
      {
        numMatches <- length(loc[[1]])
        
        for(i in 1:numMatches)
        {
          segments(loc[[1]][i], y, loc[[1]][i] + nchar(seq) - 1, y, col, lwd = wid)
        }
        msg <- 'Sequence(s) found.'
      }else
      {
        msg <- 'Sequence not found.'
      }
      
    }else
    {
      msg <- 'Current gene not found.'	
    }
    
  }else
  {
    if(is.null(currGene))
      msg <- 'Current gene is NULL. '
    
    if(is.null(seq) || empty(seq))
      msg <- paste0(msg, 'Provided message is NULL or empty. ')
    
    if(is.null(species))
      msg <- paste0(msg, 'Object \' species\' has not been initialized.')
  }
  print(msg)
  flush.console()
}

tabulateFilesForMerge <- function()
{
  dir <- choose.dir()
  subDirs <- dir
  i <- 1
  lFilesDCF <- list.files(path = subDirs[i], pattern = '*.dcf', all.files=TRUE, full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  
  vGenotype <- vector(mode = "character", length = length(lFilesDCF))
  vBatch <- vector(mode = "character", length = length(lFilesDCF))
  vLow <- vector(mode = "character", length = length(lFilesDCF))
  vHigh <- vector(mode = "character", length = length(lFilesDCF))
  vSuffix <- vector(mode = "character", length = length(lFilesDCF))
  
  for(i in 1:length(lFilesDCF))
  {
    tmp <- substr(lFilesDCF[i], gregexpr('___', lFilesDCF[i])[[1]][1] + 3, nchar(lFilesDCF[i]))
    tmp1 <- substr(tmp, 1, gregexpr('.mgf', tmp)[[1]][1]-1)
    tmp2 <- strsplit(tmp1, '_')
    vGenotype[i] <- tmp2[[1]][1]
    vBatch[i] <- tmp2[[1]][2]
    vLow[i] <- tmp2[[1]][3]
    vHigh[i] <- tmp2[[1]][4]
    vSuffix[i] <- tmp2[[1]][5]
  }
  
  df <- data.frame('file' = lFilesDCF, 'genotype' = vGenotype, 'batch' = vBatch, 'lowRange' = vLow, 'highRange' = vHigh, 'suffix' = vSuffix, stringsAsFactors = FALSE)
  return(df)
}

autoMergeFiles <- function()
{
  df <- tabulateFilesForMerge()
  gtCount <- length(unique(df$genotype))
  gTypes <- unique(df$genotype)
  bGoodToMerge <- FALSE
  sSuf <- ''
  sPath <- paste0(dirname(df$file[1]), '/')
  
  missedFiles <- 0
  lMissedFiles <- ''
  
  for(i in 1:gtCount)
  {
    indices <- which(df$genotype == gTypes[i]) 
    tmp <- df[indices,]
    while(length(tmp$genotype)>0)
    {
      cur <- tmp[1,]
      idx <- which(tmp$batch == cur$batch)
      
      if(length(idx) > 2)
      {
        if(is.na(cur$suffix))
        {
          idx <- which(is.na(tmp$suffix[idx]))
        }else
          idx <- which(tmp$suffix[idx] == cur$suffix)
      }
      
      if(length(idx) == 2)  # should only be the pair we want to merge, but let's make some checks to be sure
      {
        bGoodToMerge <- (tmp$lowRange[idx[1]] != tmp$lowRange[idx[2]]) && (tmp$highRange[idx[1]] != tmp$highRange[idx[2]]) 
        if(is.na(tmp$suffix[idx[1]]))
        {
          bGoodToMerge <- bGoodToMerge && is.na(tmp$suffix[idx[2]]) ## if both are na, that's fine
          sSuf <- ''
        }else### suffix at index 1 is NOT NA
        {
          if(!is.na(tmp$suffix[idx[2]]))## neither are NA, they should be equal then
          {
            bGoodToMerge <- bGoodToMerge && (tmp$suffix[idx[1]] == tmp$suffix[idx[2]])
            sSuf <- tmp$suffix[idx[2]]
          }else
          {
            bGoodToMerge <- FALSE
          }
        }
        
        if(bGoodToMerge)
        {
          ## close all current files
          if (!exists("fileFolder") || is.null(fileFolder))
          {
            #do nothing
          }else if(length(names(fileFolder))>0)
            fc(names(fileFolder))
          
          ## open target files
          fo(tmp$file[idx[1]])
          fo(tmp$file[idx[2]])
          ## merge files					
          mergeFiles(fileFolder, c(1,2), species, paste0(sPath, tmp$genotype[1], '_', tmp$batch[1], '_', sSuf, '_merged.dcf'))				
        }
        
      }else
      {
        for(j in 1:length(idx))
        {
          missedFiles <- missedFiles + 1
          lMissedFiles[missedFiles] <- tmp$file[idx[j]]
        }
      }
      tmp <- tmp[-idx,]
    }
  }
  if(missedFiles > 0)
  {
    write.csv(lMissedFiles, file=paste0(sPath, 'missedMerge.csv'))
  }
}

tabulateFilesForMean <- function()
{
  dir <- choose.dir()
  subDirs <- dir
  i <- 1
  lFilesDCF <- list.files(path = subDirs[i], pattern = '*.dcf', all.files=TRUE, full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  idx <- unlist(gregexpr('merge', lFilesDCF))
  idx <- which(idx != -1)
  lFilesDCF <- lFilesDCF[idx]
  
  vGenotype <- vector(mode = "character", length = length(lFilesDCF))
  vBatch <- vector(mode = "character", length = length(lFilesDCF))
  vSuffix <- vector(mode = "character", length = length(lFilesDCF))
  
  for(i in 1:length(lFilesDCF))
  {
    tmp <- substr(lFilesDCF[i], gregexpr('/', lFilesDCF[i])[[1]][1] + 1, nchar(lFilesDCF[i]))
    tmp1 <- substr(tmp, 1, gregexpr('.dcf', tmp)[[1]][1]-1)
    tmp2 <- strsplit(tmp1, '_')
    vGenotype[i] <- tmp2[[1]][1]
    vBatch[i] <- tmp2[[1]][2]
    vSuffix[i] <- tmp2[[1]][3]
  }
  
  df <- data.frame('file' = lFilesDCF, 'genotype' = vGenotype, 'batch' = vBatch, 'suffix' = vSuffix, stringsAsFactors = FALSE)
  df <- df[order(df$genotype, as.numeric(df$batch)),]
  return(df)
}

autoMeanFiles <- function()
{
  df <- tabulateFilesForMean()
  gtCount <- length(unique(df$genotype))
  gTypes <- unique(df$genotype)
  bGoodToMerge <- FALSE
  sSuf <- ''
  sPath <- paste0(dirname(df$file[1]), '/')
  
  missedFiles <- 0
  lMissedFiles <- ''
  
  for(i in 1:gtCount)
  {
    
    if (!exists("fileFolder") || is.null(fileFolder))
    {
      #do nothing
      print('fileFolder does not exist')
      
    }else if(length(names(fileFolder))>0)
      fc(names(fileFolder))
    
    indices <- which(df$genotype == gTypes[i]) 
    tmp <- df[indices,]
    
    tot <- length(tmp$genotype)
    
    if(tot > 1)
    {
      # open all files
      for(i in 1:tot)
        fo(tmp$file[i])
      # mean all the files
      math_op('mean', fileFolder, 1:tot, NULL, paste0(sPath, tmp$genotype[1], '_mean_all.dcf'))
      # check if there is a file with a different suffix
      iFirst <- which(tmp$suffix == unique(tmp$suffix)[1])
      iSecond <- which(tmp$suffix == unique(tmp$suffix)[2])
      
      if(tot > 2)
      {
        if(length(iFirst)<length(iSecond))
        {
          math_op('mean', fileFolder, iSecond, NULL, paste0(sPath, tmp$genotype[1], '_mean.dcf'))
        }else
        {
          math_op('mean', fileFolder, iFirst, NULL, paste0(sPath, tmp$genotype[1], '_mean.dcf'))
        }
      }
    }
  }
}


replace_FileListCSV <- function(sTarget = '', sReplacement = '', ext = '.csv')
{
  ## get user input on the file(s) for processing
  dir <- choose.dir()
  
  if (dir.exists(dir))
  {
    lFilesCSV <- list.files(path = dir, pattern = paste0('*', ext), all.files=TRUE, full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
    
    if (length(lFilesCSV) > 0) ## only process if there are files selected
    {
      for (j in 1:length(lFilesCSV))
      {
        idx <- gregexpr(sTarget, lFilesCSV[j], fixed=TRUE)
        if(attr(idx[[1]], "match.length") != -1)
        {
          tmp <- lFilesCSV[j]
          tmp <- gsub(sTarget, sReplacement, tmp)
          file.rename(from=lFilesCSV[j], to=tmp)
        }
      }
    }
  }
}

testGui <- function()
{
  ##creates main window
  tclCheck()
  
  for(i in 1:20)
  {
    sName <- paste0('tst', i)
    
    dlg <- myToplevel(sName)
    if (is.null(dlg))
      return(invisible())
    
    sName <- paste0('Test ', i)
    tkwm.title(dlg, sName)
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    
    testFrame <- ttkframe(dlg)
    titleString <- tclVar()
    
    sName <- paste0('Title ', i)
    
    testLabel <- ttklabel(testFrame, text=sName)
    
    
    titleEntry <- ttkentry(testFrame, width=30, justify='left', textvariable=titleString)
    
    onOk <- function()
    {	
      tkdestroy(dlg)
    }
    ok <- ttkbutton(testFrame, text='Ok', width=10, command=onOk)
    
    
    tkgrid(testFrame, column=1, row=1, sticky='nswe', pady=8, padx=2)
    tkgrid(testLabel, column=1, row=1, padx=10, pady=4, sticky='nw')
    tkgrid(titleEntry, column=2, row=1, padx=10, pady=4, sticky='ne')
    tkgrid(ok, column=1, row=2, padx=10, pady=4, sticky='se')				
  }
  
}

################################################################################
#
# Diane FileIO
#
################################################################################

###############################################################################
# loadMascot () read and order a mascot file
# 
# Author: Travis
###############################################################################
loadMascot <- function(sMascotFilename = "")
{
  ## Read Mascot file
  if (sMascotFilename == "")
    dfMascot <- read.csv( myOpen("csv", list(csv = "Comma Separated Values File, xls = Excel File"), multiple = FALSE), head = TRUE, stringsAsFactors = FALSE)
  else
    dfMascot <- read.csv( sMascotFilename, head = TRUE, stringsAsFactors = FALSE)	
  
  dfMascot <- subset(dfMascot, select = c("sample", "compound", "mScore"))
  dfMascot <- dfMascot[order(dfMascot$compound),]
  
  names(dfMascot)[1] = "genotype"
  names(dfMascot)[2] = "peptide"
  
  return(prepareMascot(dfMascot))
}
### clean this up

###############################################################################
# TODO: Add comment
# 
# Author: Travis
###############################################################################
writeDIANA <- function(saveFileName, lSpecies, result, fileStatus = NULL)
{
  endian <- "little"
  if(saveFileName == "")
  {
    saveFileName <- mySave(defaultextension = '.dcf', filetypes = list('dcf' = 'DIANA Compressed File'),	initialfile = '', title= 'Save')
    result[["File"]] <- saveFileName
  }else
  {
    result[["File"]] <- saveFileName
  }
  
  writeCon <- file(saveFileName, "w+b")
  
  writeHeader(writeCon, lSpecies, result, endian, fileStatus)
  writeMatrixNum(writeCon, result[["Genome"]], endian)
  writeSpecies_noLoop(writeCon, lSpecies, endian)
  writeMatches(writeCon, result$Matches, endian)
  
  if(any(names(result) == "appended_vectors"))
  {
    if(result$appended_vectors > 0)
    {
      for(i in 1:result$appended_vectors)
      {
        if(any(names(result) == "lVectors"))
        {
          if(result$sVecLabel[i] == globalSettings$vectorType$Mean)
          {
            writeMatrixNum(writeCon, result$lVectors[[i]], endian)
          }else if(result$sVecLabel[i] == globalSettings$vectorType$PepPoints)
          {
            writeVectorInt(writeCon, result$lVectors[[i]], endian)
            
          }else
            writeMatrixInt(writeCon, result$lVectors[[i]], endian)
        }
      }	
    }
    
  }
  
  close(writeCon)
  
  return(saveFileName)
}

readDIANA <- function(loadFileName)
{
  endian <- "little"
  if(loadFileName == "")
    loadFileName <- myOpen(defaultextension = 'dcf', filetypes = list('dcf' = 'DIANA Compressed File'),	initialfile = "", title= 'Load')
  
  readCon <- file(loadFileName, "rb")
  header <- readHeader(readCon, endian)
  mat <- readMatrixNum(readCon, endian)
  
  lSpecies <- readSpecies_noLoop(readCon, endian)
  myAssign(in.name = "species", in.object = lSpecies)	
  
  close(readCon)
  
  aaMap <- vector(mode = "numeric", length = nchar(lSpecies$seq))
  aaMap[mat[1,]] <- mat[2,]	
  
  result <- vector(mode="list", length=8)
  names(result) <- c("Genotype", "File", "Genome", "Gene", "Noise_Multiplier", "Genotype_Override", "Files_Merged", "Title_Override")
  
  result[["Genotype"]] <- header$genotype	
  
  result[["File"]] <- header$fileName	
  
  mapData <- mapToGene(lSpecies, aaMap)
  
  result[["Genome"]] <- aaMap
  result[["Gene"]] <- mapData$geneMap
  result[["Noise_Multiplier"]] <-header$noise_multiplier
  result[["Genotype_Override"]] <- header$genotype_override
  result[["Files_Merged"]] <- header$files_merged
  result[["Title_Override"]] <- header$title_override
  
  return(result)
}

writeMatrixInt <- function(writeCon, vec, endi)
{
  if (!is.null(vec))
  {
    idx <- as.integer(which(vec != 0))
    val <- as.integer(vec[idx])
    ## write number of matches
    writeLength <- as.integer(2 * length(idx))
    writeBin(writeLength, writeCon, endian = endi)
    writeBin(c(idx, val), writeCon, endian = endi)
  }else
  {
    writeBin(0, writeCon, endian = endi)
  }
}

writeMatrixNum <- function(writeCon, vec, endi)
{
  if (!is.null(vec))
  {
    idx <- as.numeric(which(vec != 0))
    val <- as.numeric(vec[idx])
    ## write number of matches
    writeLength <- as.integer(2 * length(idx))
    writeBin(writeLength, writeCon, endian = endi)
    writeBin(as.numeric(c(idx, val)), writeCon, endian = endi)
  }else
  {
    writeBin(0, writeCon, endian = endi)
  }
}

writeVectorInt <- function(writeCon, vec, endi)
{
  if (!is.null(vec))
  {
    val <- as.integer(vec)
    ## write number of matches
    writeLength <- as.integer(length(val))
    writeBin(writeLength, writeCon, endian = endi)
    writeBin(val, writeCon, endian = endi)
  }else
  {
    writeBin(0, writeCon, endian = endi)
  }
}

retrieveAAVector <- function(fileName, endian, offset = NULL)
{
  species_start <- NULL
  
  if (missing(fileName))
  {
    fileName <- sort(myOpen(defaultextension = 'dcf', filetypes = list('dcf' = 'DIANA Compressed File'),	initialfile = "", title= 'Load', multiple=FALSE))
    if(!nzchar(usrList) || !file.exists(fileName))
      return(invisible())
  }
  
  readCon <- file(fileName, "rb")
  
  if(is.null(offset))
  {
    header <- readHeader(readCon, endian)
    offset <- header$start_sparse_matrix
    species_start <- header$start_species
  }
  
  seek(readCon, where = offset, origin = "start")
  mat <- readMatrixNum(readCon, endian)
  
  ## only read species if species defined in header hasn't been loaded already
  if(!exists("species") || is.null(species) || is.null(species$ID))
  {
    if(is.null(species_start))
    {
      seek(readCon, where=0, origin="start")
      header <- readHeader(readCon, endian)
      species_start <- header$start_species
    }
    seek(readCon, where = species_start, origin = "start")
    lSpecies <- readSpecies_noLoop(readCon, endian)
    myAssign(in.name = "species", in.object = lSpecies)
  }else
  {
    lSpecies <- species
  }
  
  close(readCon)
  vecLen <- nchar(species$seq)
  aaMap <- vector(mode = "numeric", length = vecLen)
  aaMap[mat[1,]] <- mat[2,]
  return(aaMap)
}

retrieveAAVectorInt <- function(fileName, endian, offset = NULL)
{
  species_start <- NULL
  
  if (missing(fileName))
  {
    fileName <- sort(myOpen(defaultextension = 'dcf', filetypes = list('dcf' = 'DIANA Compressed File'),	initialfile = "", title= 'Load', multiple=FALSE))
    if(!nzchar(usrList) || !file.exists(fileName))
      return(invisible())
  }
  
  readCon <- file(fileName, "rb")
  
  if(is.null(offset))
  {
    header <- readHeader(readCon, endian)
    offset <- header$start_sparse_matrix
    species_start <- header$start_species
  }
  
  seek(readCon, where = offset, origin = "start")
  mat <- readMatrixInt(readCon, endian)
  
  ## only read species if species defined in header hasn't been loaded already
  if(!exists("species") || is.null(species) || is.null(species$ID))
  {
    if(is.null(species_start))
    {
      seek(readCon, where=0, origin="start")
      header <- readHeader(readCon, endian)
      species_start <- header$start_species
    }
    seek(readCon, where = species_start, origin = "start")
    lSpecies <- readSpecies_noLoop(readCon, endian)
    myAssign(in.name = "species", in.object = lSpecies)
  }else
  {
    lSpecies <- species
  }
  
  close(readCon)
  vecLen <- nchar(species$seq)
  aaMap <- vector(mode = "integer", length = vecLen)
  aaMap[mat[1,]] <- mat[2,]
  return(aaMap)
}

retrievePeptideVector <- function(fileName, endian, offset = NULL)
{
  species_start <- NULL
  
  if (missing(fileName))
  {
    fileName <- sort(myOpen(defaultextension = 'dcf', filetypes = list('dcf' = 'DIANA Compressed File'),	initialfile = "", title= 'Load', multiple=FALSE))
    if(!nzchar(usrList) || !file.exists(fileName))
      return(invisible())
  }
  
  readCon <- file(fileName, "rb")
  
  if(is.null(offset))
  {
    
  }
  
  seek(readCon, where = offset, origin = "start")
  mat <- readMatrixInt(readCon, endian)
  
  close(readCon)
  return(mat)
}


readMatrixInt <- function(readCon, endi)
{
  ## read length of input
  readLength <-	readBin(readCon, what="integer", endian=endi)
  vecLen <- readLength / 2
  
  # allocate vectors
  idx <- vector(mode = "integer", length = vecLen)
  val <- vector(mode = "integer", length = vecLen)
  
  idx <- readBin(readCon, what = "integer", n=vecLen, endian=endi)
  val <- readBin(readCon, what = "integer", n=vecLen, endian=endi)
  
  mat <- matrix(c(idx, val), nrow=2, byrow=TRUE)
  
  return(mat)
}

readMatrixNum <- function(readCon, endi)
{
  ## read length of input
  readLength <-	readBin(readCon, what="integer", endian=endi)
  vecLen <- readLength / 2
  
  # allocate vectors
  idx <- vector(mode = "numeric", length = vecLen)
  val <- vector(mode = "numeric", length = vecLen)
  
  idx <- readBin(readCon, what = "numeric", n=vecLen, endian=endi)
  val <- readBin(readCon, what = "numeric", n=vecLen, endian=endi)
  
  mat <- matrix(c(idx, val), nrow=2, byrow=TRUE)
  
  return(mat)
}

retrieveVectorInt <- function(fileName, endi, offset = NULL)
{
  
  if (missing(fileName))
  {
    fileName <- sort(myOpen(defaultextension = 'dcf', filetypes = list('dcf' = 'DIANA Compressed File'),	initialfile = "", title= 'Load', multiple=FALSE))
    if(!nzchar(usrList) || !file.exists(fileName))
      return(invisible())
  }
  
  readCon <- file(fileName, "rb")
  
  if(is.null(offset))
  {
    return(-1)
  }
  
  seek(readCon, where = offset, origin = "start")
  
  ## read length of input
  readLength <-	readBin(readCon, what="integer", endian=endi)
  
  val <- readBin(readCon, what = "integer", n=readLength, endian=endi)
  close(readCon)
  
  return(val)
}

writeMatches_old <- function(writeCon, lMatches, endi)
{
  ## write number of matches
  
  if (!is.null(lMatches))
  {
    numMatches <- length(lMatches)
    
    writeBin(numMatches, writeCon, endian = endi)
    
    for (i in 1:numMatches)
    {
      writeBin(length(lMatches[[i]]$indices), writeCon, endian = endi)
      writeBin(lMatches[[i]]$indices, writeCon, endian = endi)
      writeBin(lMatches[[i]]$pepLength, writeCon, endian = endi)
      writeBin(lMatches[[i]]$score, writeCon, endian = endi)
    }
  }else #no matches found
  {
    writeBin(0, writeCon, endian = endi)
  }
}

readMatches_old <- function(readCon, endi, startPos)
{
  seek(readCon, where=startPos, origin="start")
  ## read number of matches
  numMatches <-	readBin(readCon, what="integer", endian=endi)
  if (numMatches > 0)
  {
    lMatch <- list('indices' = 0L, 'pepLength'= 0L, 'score' = 0)	
    
    # allocate vectors
    lMatches <- vector(mode = "list", length = numMatches)
    
    for (i in 1:numMatches)
    {
      numIndices <- readBin(readCon, what = "integer", endian=endi)
      lMatch$indices <- readBin(readCon, what = "integer", n = numIndices, endian=endi)
      lMatch$pepLength <- readBin(readCon, what = "integer", endian=endi)
      lMatch$score <- readBin(readCon, what = "numeric", endian=endi)
      lMatches[[i]] <- lMatch
    }
    return(lMatches)	
  }else
  {
    return(NULL)
  }
}

matchesFileSize <- function(lMatches)
{
  ## write number of matches
  if (!is.null(lMatches))
  {
    numMatches <- length(lMatches)
    pepLength <- sapply(lMatches, function(x){x$pepLength})
    score <- sapply(lMatches, function(x){x$score})
    indexCount <- sapply(lMatches, function(x){length(x$indices)})
    indicesTot <- sum(indexCount)
    
    size <- numMatches * (8 + 4 + 4) + (indicesTot * 4)
    # for each match, numeric for score, integer for pepLength + an integer for each indexCount.
    # plus an integer to store each index position
    
  }else #no matches found
  {
    size <- 0
  }
  
  return(size + 4) ### add 4 for storing the integer size of matches whether 0 or other
}


writeMatches <- function(writeCon, lMatches, endi)
{
  ## write number of matches
  if (!is.null(lMatches))
  {
    numMatches <- length(lMatches)
    pepLength <- sapply(lMatches, function(x){x$pepLength})
    score <- sapply(lMatches, function(x){x$score})
    indexCount <- sapply(lMatches, function(x){length(x$indices)})
    indices <- unlist(sapply(lMatches, function(x){x$indices}))
    
    writeBin(as.integer(numMatches), writeCon, endian = endi) ## total number of matches
    writeBin(as.integer(pepLength), writeCon, endian = endi) ## the pepLength for each match
    writeBin(as.numeric(score), writeCon, endian = endi) ## the score for each match
    writeBin(as.integer(indexCount), writeCon, endian = endi) ## the number of indices per match
    writeBin(as.integer(indices), writeCon, endian = endi) ## all the indices
    
  }else #no matches found
  {
    writeBin(0L, writeCon, endian = endi)
  }
}

readMatches <- function(readCon, endi, startPos)
{
  seek(readCon, where=startPos, origin="start")
  ## read number of matches
  numMatches <-	readBin(readCon, what="integer", endian=endi)
  
  if (numMatches > 0)
  {
    pepLength <- readBin(readCon, what = "integer", n = numMatches, endian = endi)
    score <- readBin(readCon, what = "numeric", n = numMatches, endian = endi)
    indexCount <- readBin(readCon, what = "integer", n = numMatches, endian = endi)
    numIndices <- sum(indexCount)
    
    indices <- readBin(readCon, what = "integer", n = numIndices, endian = endi)
    cumIndices <- cumsum(indexCount)
    indexStart <- c(0, cumIndices) + 1
    
    lMatch <- list('indices' = 0L, 'pepLength'= 0L, 'score' = 0)	
    
    # allocate vectors
    lMatches <- vector(mode = "list", length = numMatches)
    
    for (i in 1:numMatches)
    {
      numIndex <- numIndices[i] 
      lMatch$indices <- indices[indexStart[i]:(indexStart[i]+indexCount[i]-1)]
      lMatch$pepLength <- pepLength[i]
      lMatch$score <- score[i]
      lMatches[[i]] <- lMatch
    }
    return(lMatches)	
  }else
  {
    return(NULL)
  }
}

retrievePeptidesList <- function(fileName, endi, startPos, lSpecies)
{
  
  if (missing(fileName))
  {
    fileName <- sort(myOpen(defaultextension = 'dcf', filetypes = list('dcf' = 'DIANA Compressed File'),	initialfile = "", title= 'Load', multiple=FALSE))
    if(!nzchar(usrList) || !file.exists(fileName))
      return(invisible())
  }
  
  readCon <- file(fileName, "rb")
  
  seek(readCon, where=startPos, origin="start")
  ## read number of matches
  numMatches <-	readBin(readCon, what="integer", endian=endi)
  
  if (numMatches > 0)
  {
    pepLength <- readBin(readCon, what = "integer", n = numMatches, endian = endi)
    score <- readBin(readCon, what = "numeric", n = numMatches, endian = endi)
    indexCount <- readBin(readCon, what = "integer", n = numMatches, endian = endi)
    numIndices <- sum(indexCount)
    
    indices <- readBin(readCon, what = "integer", n = numIndices, endian = endi)
    close(readCon)
    cumIndices <- cumsum(indexCount)
    indexStart <- c(0, cumIndices) + 1
    
    lMatch <- list('indices' = 0L, 'pepLength'= 0L)	
    pep <- list('ID' = 0L, 'sequence' = "", 'score' = 0.0)
    
    # allocate vectors
    lPeptides <- vector(mode = "list", length = numMatches)
    
    for (i in 1:numMatches)
    {
      numIndex <- numIndices[i] 
      lMatch$indices <- indices[indexStart[i]:(indexStart[i]+indexCount[i]-1)]
      lMatch$pepLength <- pepLength[i]
      pep$score <- score[i]
      pep$ID <- i
      pep$sequence <- substr(lSpecies$seq, lMatch$indices[1], lMatch$indices[1] + lMatch$pepLength -1)
      lPeptides[[i]] <- pep 
    }
    return(lPeptides)	
  }else
  {
    close(readCon)
    return(NULL)
  }
}
# readCon <- file(loadFileName, "rb")
#peps <- readPeptidesList(readCon, endian, fileFolder[[wc()]]$file.par$match_list_location, species)

writeSpecies_old <- function(writeCon, lSpecies, endi)
{
  writeBin(as.integer(nchar(lSpecies$name)), writeCon, endian=endi) ## write the length of the species name
  writeBin(lSpecies$name, writeCon, endian = endi)
  writeBin(as.integer(lSpecies$ID), writeCon, endian = endi)
  writeBin(as.integer(length(lSpecies$genes$name)), writeCon, endian = endi)
  
  for (i in 1: length(lSpecies$genes$name))
  {
    writeBin(as.integer(nchar(lSpecies$genes$name[i])), writeCon, endian=endi) ## write the length of the gene name
    writeBin(lSpecies$genes$name[i], writeCon, endian = endi) ## write gene name
    writeBin(as.integer(nchar(lSpecies$genes$chrom[i])), writeCon, endian=endi)		
    writeBin(lSpecies$genes$chrom[i], writeCon, endian = endi)
    writeBin(as.integer(lSpecies$genes$start[i]), writeCon, endian = endi)
    writeBin(as.integer(lSpecies$genes$seqStartIdx[i]), writeCon, endian = endi)
    writeBin(as.integer(lSpecies$genes$seqLength[i]), writeCon, endian = endi)
  }
  writeBin(as.integer(nchar(lSpecies$seq)), writeCon, endian = endi)
  writeBin(lSpecies$seq, writeCon, endian = endi)
}


calculateSpeciesSpace <- function(lSpecies)
{
  num_genes <- length(lSpecies$genes$name)
  space_species <- 0
  
  space_species <- 4 # to store size length of species name
  space_species <- space_species + nchar(lSpecies$name) + 1 # + 1 to include string terminator
  space_species <- space_species + 4 # + 4 to include room for integer to store species ID
  space_species <- space_species + 4 # store number of genes
  
  space_species <- space_species + (4 * num_genes) ## space for vector of gene name string lengths
  space_species <- space_species + sum(nchar(lSpecies$genes$name))+1 ## add space for concatonated string of gene names with a terminator at the end
  space_species <- space_species + (4 * num_genes) ## space for vector of chrom name string lengths
  space_species <- space_species + sum(nchar(lSpecies$genes$chrom))+1 ## add space for concatonated string of gene's chromosome name with a terminator at the end
  space_species <- space_species + (4 * num_genes) ## add space to store start vector
  space_species <- space_species + (4 * num_genes) ## add space to store sequence start index
  
  space_species <- space_species + 4 ## room to store an integer with length of sequence
  space_species <- space_species + nchar(lSpecies$seq) + 1
  return(space_species)
}

calculateSpeciesSpace_old <- function(lSpecies)
{
  space_species <- 0
  
  space_species <- 4 # to store size length of species name
  space_species <- space_species + nchar(lSpecies$name) + 1 # + 1 to include string terminator
  space_species <- space_species + 4 # + 4 to include room for integer to store species ID
  space_species <- space_species + 4 # store number of genes
  
  for (i in 1: length(lSpecies$genes$name))
  {
    space_species <- space_species + 8 #store length of gene name and chrom name in an integer
    space_species <- space_species + (nchar(lSpecies$genes$name[i])+1) # add a byte for each character and +1 for string terminator
    space_species <- space_species + (nchar(lSpecies$genes$chrom[i])+1) # add a byte for each character and +1 for string terminator		
    space_species <- space_species + (4 + 4 + 4) # add an integer to store: start, seq_start_idx, and seqLength
  }
  space_species <- space_species + 4 ## room to store an integer with length of sequence
  space_species <- space_species + nchar(lSpecies$seq) + 1
  return(space_species)
}

readSpecies_old <- function(readCon, endi)
{
  lSpecies <- vector(mode="list", length=4)
  names(lSpecies) <- c("name", "ID", "genes", "seq")	
  
  len <-	readBin(readCon, what= integer(), endian=endi)
  lSpecies$name <- readChar(readCon, len+1)
  lSpecies$ID	<- readBin(readCon, what= integer(), endian=endi)
  
  numGenes <- readBin(readCon, what= integer(), endian=endi)
  
  # allocate vectors
  name <- vector(mode = "character", length = numGenes)
  chrom <- vector(mode = "character", length = numGenes)	
  start <- vector(mode = "integer", length = numGenes)	
  seqStartIdx <- vector(mode = "integer", length = numGenes)		
  seqLength <- vector(mode = "integer", length = numGenes)	
  
  for (i in 1: numGenes)
  {
    len <-	readBin(readCon, what = integer(), endian=endi)
    name[i] <- readChar(readCon, len+1)
    len <-	readBin(readCon, what = integer(), endian=endi)
    chrom[i] <- readChar(readCon, len+1)
    start[i] <- readBin(readCon, what = integer(), endian=endi)
    seqStartIdx[i] <- readBin(readCon, what = integer(), endian=endi)
    seqLength[i] <- readBin(readCon, what = integer(), endian=endi)
  }
  
  lSpecies$genes <- data.frame("name" = name, "chrom" = chrom, 
                               "start" = start, "seqStartIdx" = seqStartIdx, 
                               "seqLength" = seqLength, stringsAsFactors = FALSE)	
  
  len <-	readBin(readCon, what = integer(), endian=endi)
  lSpecies$seq <- vector(mode = "character", length = len)
  lSpecies$seq <- readChar(readCon, len+1)  # > 4s for this operation ?....
  
  return(lSpecies)
}

writeSpecies_noLoop <- function(writeCon, lSpecies, endi)
{
  writeBin(as.integer(nchar(lSpecies$name)), writeCon, endian=endi) 				## write the length of the species name
  writeBin(lSpecies$name, writeCon, endian = endi) 								## write the species name
  writeBin(as.integer(lSpecies$ID), writeCon, endian = endi) 						## write the species ID
  writeBin(as.integer(length(lSpecies$genes$name)), writeCon, endian = endi) 		## write the number of genes
  writeBin(as.integer(nchar(lSpecies$genes$name)), writeCon, endian = endi) 		## write vector of gene name lengths (110152)
  writeBin(paste0(lSpecies$genes$name, collapse=""), writeCon, endian = endi) 	## write string of non-separated gene-names
  writeBin(as.integer(nchar(lSpecies$genes$chrom)), writeCon, endian = endi) 		## write vector of chrom name lengths
  writeBin(paste0(lSpecies$genes$chrom, collapse=""), writeCon, endian = endi) 	## write string of non-separated chrom-names
  writeBin(as.integer(lSpecies$genes$start), writeCon, endian = endi) 			## write vector of start indices
  writeBin(as.integer(lSpecies$genes$seqStartIdx), writeCon, endian = endi) 		## write vector of sequence start indices
  writeBin(as.integer(nchar(lSpecies$seq)), writeCon, endian = endi)				## write length of sequence
  writeBin(lSpecies$seq, writeCon, endian = endi)									## write character sequence
}


readSpecies_noLoop <- function(readCon, endi)
{
  lSpecies <- vector(mode="list", length=4) 
  names(lSpecies) <- c("name", "ID", "genes", "seq")
  
  len <-	readBin(readCon, what= integer(), endian=endi)			## read length of the species name
  lSpecies$name <- readChar(readCon, len+1)						## read the name of the species
  lSpecies$ID	<- readBin(readCon, what= integer(), endian=endi)	## read the species ID
  numGenes <- readBin(readCon, what= integer(), endian=endi)		## read the number of genes
  
  # allocate vectors
  names <- vector(mode = "character", length = numGenes)
  geneNameLengths <- vector(mode = "integer", length = numGenes)
  chromNameLengths <- vector(mode = "integer", length = numGenes)
  chrom <- vector(mode = "character", length = numGenes)	
  start <- vector(mode = "integer", length = numGenes)	
  seqStartIdx <- vector(mode = "integer", length = numGenes)		
  seqLength <- vector(mode = "integer", length = numGenes)	
  
  geneNameLengths <- readBin(readCon, what = integer(), n = numGenes, endian = endi) ## read the lengths of the gene names
  nameString <- readChar(readCon, sum(geneNameLengths)+1)			## read all the gene names
  geneCumLen <- cumsum(geneNameLengths)
  geneStart <- c(1, (geneCumLen + 1))
  geneStart <- geneStart[1:length(geneStart)-1]
  names <- substring(nameString, geneStart, geneCumLen)
  
  chromNameLengths <- readBin(readCon, what = integer(), n = numGenes, endian = endi)
  chromNameString <- readChar(readCon, sum(chromNameLengths)+1)
  chromCumLen <- cumsum(chromNameLengths)
  chromStart <- c(1, (chromCumLen + 1))
  chromStart <- chromStart[1:length(chromStart)-1]
  chrom <- substring(chromNameString, chromStart, chromCumLen)
  
  start <- readBin(readCon, what = integer(), n = numGenes, endian = endi)
  seqStartIdx <- readBin(readCon, what = integer(), n = numGenes, endian = endi)
  
  seqLen <-	readBin(readCon, what = integer(), endian=endi)
  
  ## use diff(seqStartIdx) instead of storing seqLength?
  ## calculate peptide sequence length
  tmp <- c(seqStartIdx, seqLen)
  seqLength <- diff(tmp)
  
  lSpecies$genes <- data.frame("name" = names, "chrom" = chrom, "start" = start, 
                               "seqStartIdx" = seqStartIdx, "seqLength" = seqLength, 
                               stringsAsFactors = FALSE)	
  
  lSpecies$seq <- vector(mode = "character", length = seqLen)
  lSpecies$seq <- readChar(readCon, seqLen+1) 
  
  return(lSpecies)
}

## files_merged = 0, genotype_override = ''
## always write version first for backwards compatability
writeHeader <- function(writeCon, lSpecies, result, endi, fileStatus = NULL)
{
  ## An integer takes 4 bytes of memory
  ## A double precision floating point number takes 8 bytes of memory
  ## a character takes one byte per character
  
  ### getStats returns a data.frame of 7 numeric values containing
  # the min, 1st quartile, median, mean, 3rd quartile, max and noise estimate information
  gNoise <- getStats(result$GeneSummary, result$GeneNoiseSD)
  matches_size <- as.integer(matchesFileSize(result$Matches))
  
  if(any(names(result) == "appended_vectors")&& !is.null(result$appended_vectors))
  {
    appended_vectors <- result$appended_vectors
    
    if(any(names(result) == "lVectors"))
    {
      vecLengths <- vector(mode = "integer", length = appended_vectors)
      vecLabelLengths <- vector(mode = "integer", length = appended_vectors)
      vecLabelLengths <- nchar(result$sVecLabel)
      
      for(i in 1:result$appended_vectors)
      {
        vecLengths[i] <- length(which(result$lVectors[[i]] != 0))	
      }
    }else
    {
      vecLengths <- -1
    }
    
  }else
  {
    appended_vectors <- 0 # by default there are no additional vectors	
  }
  
  if(as.vector(result$Gene)[1] == -1)			
  {	
    g_min_intensity <- as.numeric(-1)	# 8 bytes  
    g_max_intensity <- as.numeric(-1)	# 8 bytes
    g_zero_offset <- as.numeric(-1)		# 8 bytes
  }else
  {
    g_min_intensity <- as.numeric( min(result$Gene) )  # 8 bytes
    g_max_intensity <- as.numeric( max(result$Gene) )  # 8 bytes
    g_zero_offset <- as.numeric( median(result$Gene) ) # 8 bytes		
  }
  
  g_downfield_ppm <- as.integer(1) 						# 4 bytes
  g_upfield_ppm <- as.integer(length(lSpecies$genes$name)) # 4 bytes
  g_center_ppm <- as.integer( round( (g_upfield_ppm + g_downfield_ppm)/2 ) ) # 4 bytes
  
  ## aNoise contains another 56 bytes worth of data
  aNoise <- getStats(result$GenomeSummary, result$GenomeNoiseSD)
  
  if(as.vector(result$Genome)[1] == -1)
  {
    aa_min_intensity <- as.numeric(-1) # 8 bytes
    aa_max_intensity <- as.numeric(-1) # 8 bytes
    aa_zero_offset <- as.numeric(-1) # 8 bytes
    length_sparse_matrix <- as.integer(0)
  }else
  {
    aa_min_intensity <- as.numeric( min(result$Genome) ) # 8 bytes
    aa_max_intensity <- as.numeric( max(result$Genome) ) # 8 bytes
    aa_zero_offset <- as.numeric( round( median(result$Genome) ) ) # 8 bytes
    length_sparse_matrix <- as.integer( (length(which( result[["Genome"]] > 0 ))*2 * 8) + 4) # 4
  }
  
  if(!is.null(fileStatus))
  {
    files_merged <- fileStatus$files_merged
    genotype_override <- fileStatus$genotype_override
    noise_multiplier <- fileStatus$noise_multiplier
    title_override <- fileStatus$title_override
  }else
  {
    files_merged <- 0
    genotype_override <- ""
    noise_multiplier <- -1
    title_override <- ""
  }
  
  aa_downfield_ppm <- as.integer(1)
  aa_upfield_ppm <- as.integer(nchar(lSpecies$seq))
  aa_center_ppm <- as.integer( round( (aa_upfield_ppm + aa_downfield_ppm)/2 ) ) 
  
  version <- as.numeric(0.2)
  size_species <- as.integer(calculateSpeciesSpace(lSpecies))
  size_fileName <- as.integer(nchar(result$File))
  size_speciesName <- as.integer(nchar(lSpecies$name))
  size_genotype <- as.integer(nchar(result$Genotype))
  name_file <- result$File
  name_species <- lSpecies$name
  genotype <- result$Genotype
  id_species <- as.integer(lSpecies$ID)
  genotype_override_len <- nchar(genotype_override)
  title_override_len <- nchar(title_override)
  
  ## base size of header is 248 bytes, not including the strings
  size_header <- 248
  
  if (size_fileName > 0)
  { 
    size_header <- size_header + size_fileName + 1 ## add space for string terminator
  }	
  
  if (size_speciesName > 0)
  { 
    size_header <- size_header + size_speciesName + 1 ## add space for string terminator
  }
  
  if (size_genotype > 0)
  { 
    size_header <- size_header + size_genotype + 1 ## add space for string terminator
  }	
  
  if (genotype_override_len > 0)
  { 
    size_header <- size_header + genotype_override_len + 1 ## add space for string terminator
  }
  
  if(title_override_len > 0)
  {
    size_header <- size_header + title_override_len + 1
  }
  
  if(appended_vectors > 0)
  {
    
    size_header <- size_header + (appended_vectors * 4) # allow space for storing lengths of vectors
    size_header <- size_header + (appended_vectors * 4) + sum(vecLabelLengths) + appended_vectors # allow space for storing length of vector labels  and the labels
  }
  
  writeBin(as.numeric(version), writeCon, endian = endi) 						# 8
  writeBin(as.integer(size_header), writeCon, endian = endi)					# 4
  writeBin(length_sparse_matrix, writeCon, endian = endi)						# 4
  writeBin(size_species, writeCon, endian = endi)								# 4
  writeBin(matches_size, writeCon, endian = endi)								# 4
  writeBin(size_fileName, writeCon, endian = endi)							# 4
  writeBin(size_speciesName, writeCon, endian = endi)							# 4
  writeBin(size_genotype, writeCon, endian = endi)							# 4
  writeBin(c(name_file, name_species, genotype), writeCon, endian = endi)		# nchar(name_file, name_species, genotype) + 3
  ## append_vectors
  writeBin(as.integer(appended_vectors), writeCon, endian = endi)				# 4
  if(appended_vectors > 0)
  {
    writeBin(vecLengths, writeCon, endian = endi)
    writeBin(vecLabelLengths, writeCon, endian = endi)
    
    for(i in 1:appended_vectors)
      writeBin(result$sVecLabel[[i]] ,writeCon, endian = endi)
  }
  ##
  
  writeBin(as.integer(id_species), writeCon, endian = endi)					# 4
  
  writeBin(as.integer(files_merged), writeCon, endian = endi)					# 4
  
  writeBin(as.integer(nchar(genotype_override)), writeCon, endian = endi)		# 4
  if (genotype_override_len > 0)
    writeBin(genotype_override, writeCon, endian = endi)					#nchar(genotype_override) + 1
  
  writeBin(as.integer(nchar(title_override)), writeCon, endian = endi)		# 4
  if(title_override_len > 0)
    writeBin(title_override, writeCon, endian = endi)						#nchar(title_override) + 1
  
  writeBin(as.numeric(noise_multiplier), writeCon, endian = endi)				# 8
  writeBin(gNoise$est, writeCon, endian = endi)								# 8
  writeBin(gNoise$sd, writeCon, endian = endi)								# 8
  writeBin(gNoise$max, writeCon, endian = endi)								# 8
  writeBin(gNoise$mean, writeCon, endian = endi)								# 8
  writeBin(gNoise$lHinge, writeCon, endian = endi)							# 8
  writeBin(gNoise$uHinge, writeCon, endian = endi)							# 8
  writeBin(gNoise$min, writeCon, endian = endi)								# 8
  writeBin(g_min_intensity, writeCon, endian = endi)							# 8
  writeBin(g_max_intensity, writeCon, endian = endi)							# 8
  writeBin(g_zero_offset, writeCon, endian = endi)							# 8
  writeBin(g_downfield_ppm, writeCon, endian = endi)							# 4
  writeBin(g_upfield_ppm, writeCon, endian = endi)							# 4
  writeBin(g_center_ppm, writeCon, endian = endi)								# 4
  
  writeBin(aNoise$est, writeCon, endian = endi)								# 8
  writeBin(aNoise$sd, writeCon, endian = endi)								# 8
  writeBin(aNoise$max, writeCon, endian = endi)								# 8
  writeBin(aNoise$mean, writeCon, endian = endi)								# 8
  writeBin(aNoise$lHinge, writeCon, endian = endi)							# 8
  writeBin(aNoise$uHinge, writeCon, endian = endi)							# 8
  writeBin(aNoise$min, writeCon, endian = endi)								# 8
  writeBin(aa_min_intensity, writeCon, endian = endi)							# 8
  writeBin(aa_max_intensity, writeCon, endian = endi)							# 8
  writeBin(aa_zero_offset, writeCon, endian = endi)							# 8
  writeBin(aa_downfield_ppm, writeCon, endian = endi)							# 4
  writeBin(aa_upfield_ppm, writeCon, endian = endi)							# 4
  writeBin(aa_center_ppm, writeCon, endian = endi)							# 4
}
## always read version first to allow version specific processing
readHeader <- function(readCon, endi)
{
  ### need to figure out size of species!!
  ### File Layout as below  ###########################################
  ### |Header   | AA Match Matrix | Species Details | Match List |  ###
  ### A         B					C   		      D            E  ###
  #####################################################################
  ### Precalculate size of header and store as size_header allows us to skip to start of AA matrix (Position B)
  ### Precalculate size of  AA Match Matrix - therefore length_sparse_matrix + size_header allows
  ### us to skip from file start to start of Species Details (Position C)
  ### We can also precalculate Species Details to know where Match list starts
  ### providing a short-cut to the start of each section of data
  genotype_override <- ""
  title_override <- ""
  version_file	<- readBin(readCon, what = numeric(), endian = endi)	
  size_header 	<- readBin(readCon, what = integer(), endian = endi)
  length_sparse_matrix 	<- readBin(readCon, what = integer(), endian = endi)
  size_species 	<- readBin(readCon, what = integer(), endian = endi)
  
  speciesStartLoc <- size_header + length_sparse_matrix
  ## * 8 to convert length to a per unit length for storage
  
  if(version_file < 0.2) ## versions before 0.2 did not have size of matches 
  {	
    size_matches  <- 0
  }else
  {
    size_matches 	<- readBin(readCon, what = integer(), endian = endi)
  }
  
  size_fileName 	<- readBin(readCon, what = integer(), endian = endi)
  size_speciesName <- readBin(readCon, what = integer(), endian = endi)
  size_genotype <- readBin(readCon, what = integer(), endian = endi)
  
  name_file 		<- readChar(readCon, (size_fileName + 1))
  name_species 	<- readChar(readCon, (size_speciesName + 1))
  genotype		<- readChar(readCon, (size_genotype + 1))	
  
  if(version_file < 0.2) ## versions before 0.2 did not have appended_vectors field 
  {
    appended_vectors <- 0
    vecLengths <- 0
    vecLocs <- 0
    vecLabels <- ''
  }else
  {
    appended_vectors <- readBin(readCon, what = integer(), endian = endi)
  }
  
  if(appended_vectors > 0)
  {
    vecLengths <- readBin(readCon, what = integer(), n= appended_vectors, endian= endi)
    vecLabelLengths <- readBin(readCon, what = integer(), n= appended_vectors, endian= endi)
    
    vecLabels <- vector(mode="character", length=appended_vectors)
    
    for (i in 1:appended_vectors)
    {
      vecLabels[[i]] <- readChar(readCon, vecLabelLengths[i] + 1)
    }
    
    ## vecLocs needs to be calculated
    startLoc <- speciesStartLoc + size_species + size_matches
    vecLocs <- vector(mode="integer", length=appended_vectors)
    
    vecLocs[[1]] <- startLoc 
    ## calculate start locations of extra vectors
    if(appended_vectors > 1)
    {
      for (i in 2:appended_vectors)
      {
        if(vecLabels[[i-1]] == globalSettings$vectorType$PepPoints)
        {
          vLen <- vecLengths[[i-1]]  # PepPoints vector is not saved as a matrix.. no need to double size
        }else
          vLen <- vecLengths[[i-1]] * 2
        
        if(vecLabels[[i]] == globalSettings$vectorType$EndPoints || vecLabels[[i]] == globalSettings$vectorType$PepPoints)
          vSize <- vLen * 4
        else if(vecLabels[[i]] == globalSettings$vectorType$Mean)
          vSize <- vLen * 8
        
        vecLocs[[i]] <- vecLocs[[i-1]] + vSize + 4  	
      }
    }
  }	
  
  vecInfo <- data.frame(len = vecLengths, loc = vecLocs, label = vecLabels, stringsAsFactors = FALSE)
  
  id_species 		<- readBin(readCon, what = integer(), endian= endi)
  
  files_merged 	<- readBin(readCon, what = integer(), endian = endi)
  genotype_override_len <- readBin(readCon, what = integer(), endian = endi)
  
  if (genotype_override_len > 0)
    genotype_override <- readChar(readCon, (genotype_override_len + 1), endian = endi)	
  
  title_override_len <- readBin(readCon, what = integer(), endian = endi)
  
  if (title_override_len > 0)
    title_override <- readChar(readCon, (title_override_len + 1), endian = endi)	
  
  noise_multiplier <- readBin(readCon, what = numeric(), endian = endi)
  g_noiseEst 		<- readBin(readCon, what = numeric(), endian = endi)
  g_noiseSD 		<- readBin(readCon, what = numeric(), endian = endi)
  g_noiseMax		<- readBin(readCon, what = numeric(), endian = endi)
  g_noiseMean		<- readBin(readCon, what = numeric(), endian = endi)
  g_noiseLowHinge <- readBin(readCon, what = numeric(), endian = endi)
  g_noiseUpperHinge <- readBin(readCon, what = numeric(), endian = endi)
  g_noiseMin		<- readBin(readCon, what = numeric(), endian = endi)
  g_min_intensity <- readBin(readCon, what = numeric(), endian = endi)
  g_max_intensity <- readBin(readCon, what = numeric(), endian = endi)
  g_zero_offset 	<- readBin(readCon, what = numeric(), endian = endi)
  g_downfield_ppm <- readBin(readCon, what = integer(), endian = endi)
  g_upfield_ppm 	<- readBin(readCon, what = integer(), endian = endi)
  g_center_ppm 	<- readBin(readCon, what = integer(), endian = endi)
  
  aa_noiseEst 	<- readBin(readCon, what = numeric(), endian = endi)
  aa_noiseSD 		<- readBin(readCon, what = numeric(), endian = endi)
  aa_noiseMax 	<- readBin(readCon, what = numeric(), endian = endi)
  aa_noiseMean 	<- readBin(readCon, what = numeric(), endian = endi)
  aa_noiseLowHinge <- readBin(readCon, what = numeric(), endian = endi)
  aa_noiseUpperHinge <- readBin(readCon, what = numeric(), endian = endi)
  aa_noiseMin 	<- readBin(readCon, what = numeric(), endian = endi)		
  
  aa_min_intensity <- readBin(readCon, what = numeric(), endian = endi)
  aa_max_intensity <- readBin(readCon, what = numeric(), endian = endi)
  aa_zero_offset 	<- readBin(readCon, what = numeric(), endian = endi)
  aa_downfield_ppm <- readBin(readCon, what = integer(), endian = endi)
  aa_upfield_ppm 	<- readBin(readCon, what = integer(), endian = endi)
  aa_center_ppm 	<- readBin(readCon, what = integer(), endian = endi)
  
  header <- list(	
    start_sparse_matrix = size_header,	start_species = speciesStartLoc, 
    start_match_list = speciesStartLoc + size_species,			
    version = version_file,	speciesID = id_species,
    species = name_species,	fileName = name_file,
    g_noiseEst = g_noiseEst, g_noiseSD = g_noiseSD, g_noiseMax = g_noiseMax, g_noiseMean = g_noiseMean, 
    g_noiseLowHinge = g_noiseLowHinge, g_noiseUpperHinge = g_noiseUpperHinge, g_noiseMin = g_noiseMin,
    g_min_int = g_min_intensity, g_max_int = g_max_intensity, g_zero = g_zero_offset,
    g_downfield = g_downfield_ppm, g_upfield = g_upfield_ppm, g_center = g_center_ppm,
    aa_noiseEst = aa_noiseEst, aa_noiseSD = aa_noiseSD, aa_noiseMax = aa_noiseMax, aa_noiseMean = aa_noiseMean, 
    aa_noiseLowHinge = aa_noiseLowHinge, aa_noiseUpperHinge = aa_noiseUpperHinge, aa_noiseMin = aa_noiseMin,
    aa_min_int = aa_min_intensity, aa_max_int = aa_max_intensity, aa_zero = aa_zero_offset,
    aa_downfield = aa_downfield_ppm, aa_upfield = aa_upfield_ppm, aa_center = aa_center_ppm,
    genotype = genotype, genotype_override = genotype_override, files_merged = files_merged, 
    noise_multiplier = noise_multiplier, title_override = title_override, 
    appended_vectors = appended_vectors, vecInfo = vecInfo)
  
  return(header)	
}

saveAs <- function(saveFileName, in.folder, mapData = NULL, includeMatches = TRUE, lVectors = NULL, sVecLabel = NULL)
{
  
  files_merged <- in.folder$file.par$files_merged
  genotype_override <- in.folder$file.par$genotype_override
  noise_multiplier <- in.folder$file.par$noise_override
  title_override <- in.folder$file.par$title_override
  appended_vectors <- in.folder$file.par$appended_vectors
  
  if(includeMatches)
  {
    readCon <- file(in.folder$file.par$file.name, "rb")
    lMatches <- readMatches(readCon, "little", in.folder$file.par$match_list_location)
    close(readCon)
  }else
    lMatches <- NULL
  
  if(is.null(mapData))
  {
    mapData <- diana(file.name = in.folder$file.par$file.name, file.par = in.folder$file.par)$mapData 
  }
  
  result <- vector(mode="list", length=12)
  
  names(result) <- c("Genotype", "File", "Matches", "Genome", "Gene", "GenomeSummary", 
                     "GenomeNoiseSD", "GeneSummary", "GeneNoiseSD", "appended_vectors", "lVectors", "sVecLabel")
  statsum <- vector(mode="numeric", length=6)
  
  result$Genotype <- in.folder$file.par$genotype
  result$File <- saveFileName
  result$Matches <-  lMatches
  result$Genome <- mapData$aminoAcidMap
  result$Gene <- mapData$geneMap
  result$GenomeSummary <- c(in.folder$file.par$aa_noise_min, in.folder$file.par$aa_noise_low_hinge, 0, 
                            in.folder$file.par$aa_noise_mean, in.folder$file.par$aa_noise_upper_hinge, 
                            in.folder$file.par$aa_noise_max)
  
  result$GenomeNoiseSD <- in.folder$file.par$aa_noise_sd
  
  result$GeneSummary <- c(in.folder$file.par$noise_min, in.folder$file.par$noise_low_hinge, 0, 
                          in.folder$file.par$noise_mean, in.folder$file.par$noise_upper_hinge, 
                          in.folder$file.par$noise_max)
  
  result$GeneNoiseSD <- in.folder$file.par$noise_sd
  result$lVectors <- vector(mode = "list")
  
  if(!is.null(lVectors))
  {
    result$appended_vectors <- appended_vectors + 1
    result$lVectors[[result$appended_vectors]] <- lVectors
    if(!is.null(sVecLabel))
    {
      result$sVecLabel <- sVecLabel
    }
  }
  
  fileStatus <- vector(mode="list", length=4)
  names(fileStatus) <- c("files_merged", "genotype_override", "noise_multiplier", "title_override")
  
  fileStatus$files_merged <- files_merged
  fileStatus$genotype_override <- genotype_override
  fileStatus$noise_multiplier <- noise_multiplier
  fileStatus$title_override <- title_override
  
  return(writeDIANA(saveFileName, species, result, fileStatus))		
}

################################################################################
#
# Prepare species files
#
################################################################################
#
# Human proteome
# 
###############################################################################

initHuman <- function(sFileName = '')
{
  if (sFileName == '')
    sFileName <- paste(getwd(), '/chrHs.csv', collapse='', sep='')
  
  human <- loadHumanGenes(sFileName)
  
  return(human)
}

loadHumanGenes <- function(sHumanFilename = "")
{
  ## Read Human Proteome to Genome Map file
  if (sHumanFilename == "")
  {
    dfHuman <- read.csv( myOpen("csv", list(csv = "Comma Separated Values File, xls = Excel File"), multiple = FALSE), head = TRUE, stringsAsFactors = FALSE)
  }else{
    dfHuman <- read.csv( sHumanFilename, head = TRUE, stringsAsFactors = FALSE)
  }
  
  fileInfo <- file.info(sHumanFilename)
  humanID <- as.integer(fileInfo$mtime)
  
  ## rename "chromosome" to "chrom" for consistency and brevity
  names(dfHuman)[names(dfHuman) == 'chromosome'] <- "chrom"
  
  ## Fix chromosome name problems
  ## make sure all x and y chromosomes are in same case
  dfHuman$chrom[which(dfHuman$chrom == "Y")] <- "y"
  dfHuman$chrom[which(dfHuman$chrom == "X")] <- "x"
  ## make sure unknown designations share a tag
  dfHuman$chrom[grep("CHR", dfHuman$chrom )] <- "UKN"
  dfHuman$chrom[grep("GL", dfHuman$chrom )] <- "UKN"
  dfHuman$chrom[grep("KI", dfHuman$chrom )] <- "UKN"
  dfHuman$chrom[dfHuman$chrom == ""] <- "UKN"
  dfHuman <- subset(dfHuman, select = c("GeneName", "seq", "chrom", "start"))
  names(dfHuman)[1] <- "name"
  
  return(prepareSpecies("human", humanID, dfHuman))
}

################################################################################
# 
# Human Sickle
#
################################################################################
initHuman <- function(sFileName = '')
{
  if (sFileName == '')
    sFileName <- paste(getwd(), './sickle.csv', collapse='', sep='')
  
  human <- loadHumanGenes(sFileName)
  
  return(human)
}

loadHumanGenes <- function(sHumanFilename = "")
{
  ## Read Human Proteome to Genome Map file
  if (sHumanFilename == "")
  {
    dfHuman <- read.csv( myOpen("csv", list(csv = "Comma Separated Values File, xls = Excel File"), multiple = FALSE), head = TRUE, stringsAsFactors = FALSE)
  }else{
    dfHuman <- read.csv( sHumanFilename, head = TRUE, stringsAsFactors = FALSE)
  }
  
  fileInfo <- file.info(sHumanFilename)
  humanID <- as.integer(fileInfo$mtime)
  
  ## rename "chromosome" to "chrom" for consistency and brevity
  names(dfHuman)[names(dfHuman) == 'chromosome'] <- "chrom"
  
  ## Fix chromosome name problems
  ## make sure all x and y chromosomes are in same case
  dfHuman$chrom[which(dfHuman$chrom == "Y")] <- "y"
  dfHuman$chrom[which(dfHuman$chrom == "X")] <- "x"
  ## make sure unknown designations share a tag
  dfHuman$chrom[grep("CHR", dfHuman$chrom )] <- "UKN"
  dfHuman$chrom[grep("GL", dfHuman$chrom )] <- "UKN"
  dfHuman$chrom[grep("KI", dfHuman$chrom )] <- "UKN"
  dfHuman$chrom[dfHuman$chrom == ""] <- "UKN"
  dfHuman <- subset(dfHuman, select = c("GeneName", "seq", "chrom", "start"))
  names(dfHuman)[1] <- "name"
  
  return(prepareSpecies("human", humanID, dfHuman))
}

################################################################################
#
# BosTaurus Proteome
#
################################################################################

initTaurus <- function(sFileName = '')
{
  if (sFileName == '')
    sFileName <- paste(getwd(), '/chrBtaurus1.csv', collapse='', sep='')
  
  taurus <- loadTaurusGenes(sFileName)
  
  return(taurus)
}

loadTaurusGenes <- function(sTaurusFilename = "")
{
  ## Read Human Proteome to Genome Map file
  if (sTaurusFilename == "")
  {
    dfTaurus <- read.csv( myOpen("csv", list(csv = "Comma Separated Values File, xls = Excel File"), multiple = FALSE), head = TRUE, stringsAsFactors = FALSE)
  }else{
    dfTaurus <- read.csv( sTaurusFilename, head = TRUE, stringsAsFactors = FALSE)
  }
  
  fileInfo <- file.info(sTaurusFilename)
  TaurusID <- as.integer(fileInfo$mtime)
  
  ## rename "chromosome" to "chrom" for consistency and brevity
  names(dfTaurus)[names(dfTaurus) == 'chromosome'] <- "chrom"
  
  ## make sure all x chromosomes are in same case
  dfTaurus$chrom[which(dfTaurus$chrom == "X")] <- "x"
  
  ## make sure unknown designations share a tag
  dfTaurus$GeneName[dfTaurus$GeneName == ""] <- "UKN"
  dfTaurus$SwissID[dfTaurus$SwissID == ""] <- "UKN"
  dfTaurus <- subset(dfTaurus, select = c("GeneName", "seq", "chrom", "start"))
  names(dfTaurus)[1] <- "name"
  
  return(prepareSpecies("BosTaurus", TaurusID, dfTaurus))
}

################################################################################
#
# Plasmodium falciparum
#
################################################################################

initPf <- function(sFileName = '')
{
  if (sFileName == '')
    sFileName <- paste(getwd(), '/3D7Pf.csv', collapse='', sep='')	
  
  plasmodium <- loadPfGenes(sFileName)
  
  return(lPf)
}

loadPfGenes <- function(sPfFileName = "")
{
  if (sPfFileName == "")
  {
    dfPf <- read.csv( myOpen("csv", list(csv = "Comma Separated Values File, xls = Excel File"), multiple = FALSE), head = TRUE, stringsAsFactors = FALSE)
  }else{
    dfPf <- read.csv( sPfFileName, head = TRUE, stringsAsFactors = FALSE)
  }
  
  fileInfo <- file.info(sPfFileName)
  pfID <- as.integer(fileInfo$mtime)
  
  dfPf <- subset(dfPf, select = c("geneID", "seq", "chromosome", "start", "end"))
  names(dfPf)[1] <- "name"
  names(dfPf)[3] <- "chrom"
  
  return(prepareSpecies("pf", pfID, dfPf))
}

################################################################################
#
# Pseudomonas aerginosa 
#
################################################################################

###############################################################################
#	initPAO1Genes()
# 2017-11-08
# Author: Dani Kilani
###############################################################################
initPAO1 <- function(sFileName = '')
{
  if (sFileName == '')
    sFileName <- paste(getwd(), '/chrPAO1.csv', collapse='', sep='')
  
  dfPAO1 <- loadPAGenes(sFileName)
  dfPAO1$genes <- dfPAO1$genes[order(dfPAO1$genes$start),]
  dfPAO1$genes <- dfPAO1$genes[order(dfPAO1$genes$chrom),]
  
  lSpecies <- list()
  lSpecies$name <- "PAO1"
  lSpecies$ID <- dfPAO1$ID
  lSpecies$genes <- dfPAO1$genes
  lSpecies$seq <- paste(dfPAO1$genes$seq, collapse = '', sep = '')
  
  assign("lSpecies", lSpecies, envir=.GlobalEnv)
  
  return(lSpecies)
}

###############################################################################
#	loadPA14Genes()
# 2017-08-11
# Author: Dani Kilani
###############################################################################
initPA14 <- function(sFileName = '')
{
  if (sFileName == '')
    sFileName <- paste(getwd(), '/chrPA14.csv', collapse='', sep='')
  
  dfPA14 <- loadPAGenes(sFileName)
  
  dfPA14$genes <- dfPA14$genes[order(dfPA14$genes$start),]
  dfPA14$genes <- dfPA14$genes[order(dfPA14$genes$chrom),]	
  
  lSpecies <- list()
  lSpecies$name <- "PA14"
  lSpecies$ID <- dfPA14$ID
  lSpecies$genes <- dfPA14$genes
  lSpecies$seq <- paste(dfPA14$genes$seq, collapse = '', sep = '')
  
  assign("lSpecies", lSpecies, envir=.GlobalEnv)
  
  return(lSpecies)
}

###############################################################################
#	loadPAO1Genes()
# 2017-11-08
# Author: Dani Kilani
###############################################################################
loadPAO1Genes <- function(sPAO1Filename = "")
{
  ## Read PAO1 Proteome to Genome Map file
  if (sPAO1Filename == "")
  {
    dfPAO1 <- read.csv( myOpen("csv", list(csv = "Comma Separated Values File, xls = Excel File"), multiple = FALSE), head = TRUE, stringsAsFactors = FALSE)
  }else{
    dfPAO1 <- read.csv( sPAO1Filename, head = TRUE, stringsAsFactors = FALSE)
  }
  
  fileInfo <- file.info(sPAO1Filename)
  paO1ID <- as.integer(fileInfo$mtime)
  
  ## rename column labels  for consistency and brevity
  names(dfPAO1)[names(dfPAO1) == 'Chromosome'] <- "chrom"
  names(dfPAO1)[names(dfPAO1) == 'Stop'] <- "end"
  names(dfPAO1)[names(dfPAO1) == 'Start'] <- "start"
  names(dfPAO1)[names(dfPAO1) == 'Amino_Acid_Sequence'] <- "seq"
  names(dfPAO1)[names(dfPAO1) == 'Frame_Number'] <- "frame" 
  dfPAO1 <- subset(dfPAO1, select = c("Gene_Name", "seq", "chrom", "start", "end", "frame"))
  names(dfPAO1)[1] <- "name"
  dfPAO1 <- dfPAO1[order(dfPAO1$frame),]
  
  lPAO1 <- list(genes = dfPAO1, ID = paO1ID)
  return(lPAO1)
}

###############################################################################
#	loadPAGenes()
# 2017-11-08
# Author: Dani Kilani
###############################################################################
loadPAGenes <- function(sPAFilename = "")
{
  ## Read PA Proteome to Genome Map file
  if (sPAFilename == "")
  {
    dfPA <- read.csv( myOpen("csv", list(csv = "Comma Separated Values File, xls = Excel File"), multiple = FALSE), head = TRUE, stringsAsFactors = FALSE)
  }else{
    dfPA <- read.csv( sPAFilename, head = TRUE, stringsAsFactors = FALSE)
  }
  
  fileInfo <- file.info(sPAFilename)
  paID <- as.integer(fileInfo$mtime)	
  
  ## rename column labels  for consistency and brevity
  names(dfPA)[names(dfPA) == 'Chromosome'] <- "chrom"
  names(dfPA)[names(dfPA) == 'Stop'] <- "end"
  names(dfPA)[names(dfPA) == 'Start'] <- "start"
  names(dfPA)[names(dfPA) == 'Amino_Acid_Sequence'] <- "seq"
  names(dfPA)[names(dfPA) == 'Frame_Number'] <- "frame"
  
  #	dfPA <- subset(dfPA, select = c("Gene_Name", "seq", "chrom", "start", "end", "frame"))
  ## Substitute frame for chromosome value TSB 2017-15-08 	
  dfPA <- subset(dfPA, select = c("Gene_Name", "seq", "frame", "start", "end"))#, "frame"))
  
  names(dfPA)[1] <- "name"
  ## name frame as chrom TSB 2017-15-08	
  names(dfPA)[3] <- "chrom"
  #	dfPA <- dfPA[order(dfPA$frame),]
  dfPA <- dfPA[order(dfPA$chrom),]
  lPA <- list(genes = dfPA, ID = paID)
  return(lPA)
}

###############################################################################
#	loadPA14Genes()
# 2017-11-08
# Author: Dani Kilani
###############################################################################
loadPA14Genes <- function(sPA14Filename = "")
{
  ## Read PA14 Proteome to Genome Map file
  if (sPA14Filename == "")
  {
    dfPA14 <- read.csv( myOpen("csv", list(csv = "Comma Separated Values File, xls = Excel File"), multiple = FALSE), head = TRUE, stringsAsFactors = FALSE)
  }else{
    dfPA14 <- read.csv( sPA14Filename, head = TRUE, stringsAsFactors = FALSE)
  }
  
  fileInfo <- file.info(sPA14Filename)
  pa14ID <- as.integer(fileInfo$mtime)	
  
  ## rename column labels  for consistency and brevity
  names(dfPA14)[names(dfPA14) == 'Chromosome'] <- "chrom"
  names(dfPA14)[names(dfPA14) == 'Stop'] <- "end"
  names(dfPA14)[names(dfPA14) == 'Start'] <- "start"
  names(dfPA14)[names(dfPA14) == 'Amino_Acid_Sequence'] <- "seq"
  names(dfPA14)[names(dfPA14) == 'Frame_Number'] <- "frame" 
  dfPA14 <- subset(dfPA14, select = c("Gene_Name", "seq", "chrom", "start", "end", "frame"))
  names(dfPA14)[1] <- "name"
  dfPA14 <- dfPA14[order(dfPA14$frame),]
  
  lPA14 <- list(genes = dfPA14, ID = pa14ID)
  return(lPA14)
}

################################################################################
#                                                                              #
# New functions DD RA                                                          #
#                                                                              #
################################################################################

################################################################################
#
# Libraries to load
library(tcltk)
library(ggplot2)
library(ggridges)
library(dplyr)
#
################################################################################
#Remove duplicates GUI.
#Modify the code to allow for the selection of one or multiple files => check pm function by Travis
#Add close when clicked on the button. 

rd <- function() {
  ## creates main window
  tclCheck()
  dlg <- tktoplevel()
  if (is.null(dlg))
    return(invisible())
  
  tkwm.title(dlg, 'Prepare Mascot file for DigestR')
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  
  ## create Remove Duplicates button
  onRemoveDuplicates <- function() {
    source_file <- tclvalue(tcl("tk_getOpenFile", "-initialdir", ".", "-filetypes", "{{CSV Files} {.csv}}"))
    # Call the RemoveDuplicate() function passing the selected file
    result_file <- removeDuplicates(source_file)
    # Perform actions with the result_file, such as displaying or opening it
  }
  removeDuplicatesButton <- ttkbutton(dlg, text='Remove Duplicates', width=15, command=onRemoveDuplicates)
  
  ## add button to the main window
  tkgrid(removeDuplicatesButton, column=1, row=1, padx=5, pady=5)
  
  ## make main window stretch when resized
  tkgrid.columnconfigure(dlg, 1, weight=1)
  tkgrid.rowconfigure(dlg, 1, weight=1)
  
  tkwait.window(dlg)
}

removeDuplicates <- function(input_file) {
  # Read the input CSV file
  data <- read.csv(input_file, skip = 3, header = TRUE)
  
  # Check if 'pep_seq' column exists
  if (!("pep_seq" %in% colnames(data))) {
    stop("Error: 'pep_seq' column not found in the input file.")
  }
  
  # Remove duplicate values in the pep_seq column
  rd_data <- data[!duplicated(data$pep_seq), ]
  
  # Prompt the user to select the destination CSV file
  dest_file <- tclvalue(tcl("tk_getSaveFile", "-initialdir", ".", "-defaultextension", ".csv", "-filetypes", "{{CSV Files} {.csv}}"))
  
  # Write the unique data to the destination CSV file
  write.csv(rd_data, file = dest_file, row.names = FALSE)
  
  # Return the path of the destination CSV file
  return(dest_file)
}

# Test the RemoveDuplicate() function
#rd()

################################################################################
# Unique peptides function  GUI 
## Change code so that it uses the current directory to export the file to.

up <- function() {
  ## creates main window
  tclCheck()
  dlg <- tktoplevel()
  if (is.null(dlg))
    return(invisible())
  tkwm.title(dlg, 'Unique peptides')
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  
  ## create select Query button
  L1Var <- tclVar()  # Variable to store the path of the selected L1 file
  onSelectL1 <- function() {
    selectedL1 <- tclvalue(tcl("tk_getOpenFile", "-initialdir", ".", "-filetypes", "{{CSV Files} {.csv}}"))
    tclvalue(L1Var) <- selectedL1
  }
  selectL1Button <- ttkbutton(dlg, text='Query', width=15, command=onSelectL1)
  
  ## create select L2 button
  L2Var <- tclVar()  # Variable to store the path of the selected L2 file
  onSelectL2 <- function() {
    selectedL2 <- tclvalue(tcl("tk_getOpenFile", "-initialdir", ".", "-filetypes", "{{CSV Files} {.csv}}"))
    tclvalue(L2Var) <- selectedL2
  }
  selectL2Button <- ttkbutton(dlg, text='Experimental', width=15, command=onSelectL2)
  
  ## create Remove Duplicates button
  onuniquePeptides <- function() {
    L1 <- read.csv(tclvalue(L1Var), header=TRUE)
    L2_file <- tclvalue(L2Var)
    L2 <- read.csv(L2_file, header=TRUE)
    
    # Check if 'pep_seq' column exists in L1 and L2
    if (!("pep_seq" %in% colnames(L1)) || !("pep_seq" %in% colnames(L2))) {
      stop("Error: 'pep_seq' column not found in the input files.")
    }
    
    # Find unique peptides in L2 not present in L1
    L2_uniques <- L2[!(L2$pep_seq %in% L1$pep_seq), ]
    
    # Get the name of the L2 file without extension
    l2_filename <- tools::file_path_sans_ext(basename(L2_file))
    
    # Create the new filename with "_unique" attached
    dest_file <- paste0(l2_filename, "_Unique.csv")
    
    # Write the unique data to the destination CSV file
    write.csv(L2_uniques, file = dest_file, row.names = FALSE)
    
    tkmessageBox(message=paste("Unique peptides saved to:", dest_file))
  }
  uniquePeptidesButton <- ttkbutton(dlg, text='Experimental unique', width=20, command=onuniquePeptides)
  
  ## add widgets to the main window
  tkgrid(selectL1Button, column=1, row=1, padx=5, pady=5)
  tkgrid(selectL2Button, column=1, row=2, padx=5, pady=5)
  tkgrid(uniquePeptidesButton, column=1, row=3, padx=5, pady=5)
  
  ## make the main window stretch when resized
  tkgrid.columnconfigure(dlg, 1, weight=1)
  tkgrid.rowconfigure(dlg, 1, weight=1)
  
  tkwait.window(dlg)
}

# Test the up() function
#up()

################################################################################
# Density plot function GUI 

# pd <- function() {
#   tclCheck()
#   dlg <- myToplevel('pd')
#   if (is.null(dlg))
#     return(invisible())

#' A Function from digest.R
#'
#' This function does something even more interesting.
#'
#' @export
pd <- function() {
  # Define variables
  csv_dir <- ""
  plot_type <- ""
  
  # Function to select the CSV directory
  selectCSVDir1 <- function() {
    csv_dir <<- tclvalue(tcl("tk_chooseDirectory"))
    csv_dir <- tclVar(csv_dir)  # Update the csv_dir variable
    tkconfigure(csvDirEntry, text = csv_dir)  # Update the csvDirEntry widget
  }
  
  # Function to set the selected plot type
  setPlotType <- function(value) {
    plot_type <<- value
  }
  
  group_csv_files <- function(csv_dir) {
    # Get the list of CSV files in the directory
    csv_files <- list.files(path = csv_dir, pattern = "*.csv", full.names = TRUE)
    
    # Create an empty list to store the grouped dataframes
    grouped_df <- list()
    
    # Loop through the CSV files and group them based on the second string
    for (i in 1:length(csv_files)) {
      # Get the filename without the directory path
      filename <- basename(csv_files[i])
      # Get the second string in the filename by splitting on "_"
      second_str <- strsplit(filename, "_")[[1]][2]
      # Read the CSV file and extract the pep_seq column
      csv_data <- read.csv(csv_files[i], skip = 0, header = TRUE)
      pep_seq_col <- csv_data$pep_seq
      # Create a new dataframe with the pep_seq column and the group name
      pep_seq_df <- data.frame(group_name = second_str, pep_seq = pep_seq_col)
      # If the group already exists in the list, append the dataframe
      if (second_str %in% names(grouped_df)) {
        grouped_df[[second_str]] <- rbind(grouped_df[[second_str]], pep_seq_df)
      }
      # If the group doesn't exist in the list, create a new list element
      else {
        grouped_df[[second_str]] <- pep_seq_df
      }
    }
    
    # Return the list of grouped dataframes
    return(grouped_df)
  }
  
  # Function to create the density plot
  createDensityPlot <- function() {
    grouped_df <- group_csv_files(csv_dir)
    
    # Create an empty data frame to store the nchar values for all groups
    nchar_df <- data.frame(group_name = character(), nchar = integer())
    
    for (i in 1:length(grouped_df)) {
      group_name <- names(grouped_df)[i]
      pep_seq_col <- grouped_df[[i]]$pep_seq
      group_nchar_df <- data.frame(group_name = group_name, nchar = nchar(pep_seq_col))
      nchar_df <- rbind(nchar_df, group_nchar_df)
    }
    
    nchar_mean_df <- nchar_df %>%
      group_by(group_name) %>%
      summarise(mean_nchar = mean(nchar))
    
    if (plot_type == "Overlay") {
      plot <- ggplot(nchar_df, aes(x = nchar, group = group_name, fill = group_name)) +
        geom_density(alpha = 0.5) +
        labs(title = "Peptide Length Distribution", x = "Peptide length in AA", y = "", fill = "Group") +
        theme_bw() +
        geom_vline(data = nchar_mean_df, aes(xintercept = mean_nchar, color = group_name), size = 1)
    } else if (plot_type == "Ridges") {
      plot <- ggplot(nchar_df, aes(x = nchar, y = group_name, fill = group_name)) +
        geom_density_ridges2(alpha = 0.5, rel_min_height = 0.01, scale = 7, quantile_lines = TRUE, quantile_fun = function(x, ...) mean(x)) +
        labs(title = "Peptide Length Distribution", x = "Peptide length in AA", y = "", fill = "Group") +
        theme_bw()
    } else if (plot_type == "Colored Ridges") {
      plot <- ggplot(nchar_df, aes(x = nchar, y = group_name, fill = stat(x))) +
        geom_density_ridges_gradient(alpha = 0.5, rel_min_height = 0.01, scale = 3, quantile_lines = TRUE, quantile_fun = function(x, ...) mean(x)) +
        labs(title = "Peptide Length Distribution", x = "Peptide length in AA", y = "", fill = "Group") +
        scale_fill_viridis_c(name = "Length", option = "H")
    } else {
      stop("Invalid plot_type argument. Please choose either 'Overlay', 'Ridges', or 'Colored Ridges'.")
    }
    
    print(plot)
  }
  
  # Create the main window
  win <- tktoplevel()
  tkwm.title(win, "Density Plot Creator")
  
  # Create and configure the CSV directory selection frame
  csvDirFrame <- ttkframe(win)
  tkgrid(csvDirFrame, padx = 20, pady = 20)
  
  # Create and configure the CSV directory label
  csvDirLabel <- ttklabel(csvDirFrame, text = "CSV Directory:")
  tkgrid(csvDirLabel, column = 1, row = 1, sticky = "w")
  
  # Create and configure the CSV directory entry
  csvDirEntry <- ttkentry(csvDirFrame, textvariable = csv_dir)
  tkgrid(csvDirEntry, column = 2, row = 1)
  
  # Create and configure the CSV directory browse button
  csvDirButton <- ttkbutton(csvDirFrame, text = "Browse", command = selectCSVDir1)
  tkgrid(csvDirButton, column = 3, row = 1, padx = 5)
  
  # Create and configure the plot type selection frame
  plotTypeFrame <- ttkframe(win)
  tkgrid(plotTypeFrame, padx = 10, pady = 10)
  
  # Create and configure the plot type label
  plotTypeLabel <- ttklabel(plotTypeFrame, text = "Plot Type:")
  tkgrid(plotTypeLabel, column = 1, row = 1, sticky = "w")
  
  # Create and configure the plot type radiobuttons
  plotTypeRadioFrame <- ttkframe(plotTypeFrame)
  tkgrid(plotTypeRadioFrame, column = 2, row = 1)
  
  # Create and configure the Overlay radiobutton
  overlayRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Overlay", variable = plot_type, value = "Overlay", command = function() setPlotType("Overlay"))
  tkgrid(overlayRadio, column = 1, row = 1, sticky = "w")
  
  # Create and configure the Ridges radiobutton
  ridgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Ridges", variable = plot_type, value = "Ridges", command = function() setPlotType("Ridges"))
  tkgrid(ridgesRadio, column = 1, row = 2, sticky = "w")
  
  # Create and configure theColored Ridges radiobutton
  coloredRidgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Colored Ridges", variable = plot_type, value = "Colored Ridges", command = function() setPlotType("Colored Ridges"))
  tkgrid(coloredRidgesRadio, column = 1, row = 3, sticky = "w")
  
  # Create and configure the create plot button
  createPlotButton <- ttkbutton(win, text = "Create Plot", command = createDensityPlot)
  tkgrid(createPlotButton, pady = 10)
  
  # Start the event loop
  tkwait.visibility(win)
}
# call the function
# pd()

################################################################################ 
# Cter, Nter function GUI 

csp <- function () {
  # tclCheck()
  # dlg <- myToplevel('csp')
  # if (is.null(dlg))
  #   return(invisible())
  
  # Define variables
  csv_directory <- ""
  PlotType <- ""
  
  # Function to select the CSV directory
  selectCSVDir <- function() {
    csv_directory <<- tclvalue(tcl("tk_chooseDirectory"))
    csv_directory <- tclVar(csv_directory)  # Update the csv_directory variable
    tkconfigure(csvDirEntry, text = csv_directory)  # Update the csvDirEntry widget
  }
  
  # Function to set the selected plot type
  setPlotType <- function(value) {
    PlotType <<- value
  }
  
  # Function to create the bar plot
  createBarPlot <- function() 
    #plotCutSites <- function(csv_directory, PlotType) 
  {
    library(dplyr)
    library(ggplot2)
    
    extract_first_last_letters_from_csv <- function(file_path) {
      # Read the CSV file skipping the first 3 lines (assuming headers are in line 4)
      data <- read.csv(file_path, skip = 0, header = TRUE)
      
      # Check if the 'pep_seq' column exists
      if (!("pep_seq" %in% colnames(data))) {
        stop("Error: 'pep_seq' column not found in the CSV file.")
      }
      
      # Extract the first and last letters from each sequence
      data$first_letter <- substr(data$pep_seq, 1, 1)
      data$last_letter <- substr(data$pep_seq, nchar(data$pep_seq), nchar(data$pep_seq))
      
      # Count the occurrences of each first letter
      first_letter_counts <- table(factor(data$first_letter, levels = c("A", "C", "D", "E", "F", "G","H" ,"K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")))
      
      # Count the occurrences of each last letter
      last_letter_counts <- table(factor(data$last_letter, levels = c("A", "C", "D", "E", "F", "G","H", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")))
      
      # Combine the counts for first and last letters
      combined_counts <- merge(
        data.frame(letter = names(first_letter_counts), first_letter_count = as.numeric(first_letter_counts)),
        data.frame(letter = names(last_letter_counts), last_letter_count = as.numeric(last_letter_counts)),
        by = "letter", all = TRUE
      )
      combined_counts$combined_count <- rowSums(combined_counts[, c("first_letter_count", "last_letter_count")], na.rm = TRUE)
      
      # Replace NA values with 0
      combined_counts[is.na(combined_counts)] <- 0
      
      # Return the combined counts
      return(combined_counts)
    }
    
    group_csv_files <- function(csv_directory) {
      # Get the list of CSV files in the directory
      csv_files <- list.files(path = csv_directory, pattern = "*.csv", full.names = TRUE)
      
      # Create an empty list to store the grouped CSV files
      grouped_csv <- list()
      
      # Loop through the CSV files and group them based on the second string
      for (i in 1:length(csv_files)) {
        # Get the filename without the directory path
        filename <- basename(csv_files[i])
        # Get the second string in the filename by splitting on "_"
        second_str <- strsplit(filename, "_")[[1]][2]
        # If the group already exists in the list, append the CSV file
        if (second_str %in% names(grouped_csv)) {
          grouped_csv[[second_str]] <- c(grouped_csv[[second_str]], csv_files[i])
        }
        # If the group doesn't exist in the list, create a new list element
        else {
          grouped_csv[[second_str]] <- list(csv_files[i])
        }
      }
      
      # Return the list of grouped CSV files
      return(grouped_csv)
    }
    
    # Group the CSV files
    csv <- group_csv_files(csv_directory)
    
    # Initialize an empty dataframe to store the results
    dat <- data.frame()
    
    # Iterate over each group of CSV files
    for (groups in csv) {
      
      # Create an empty dataframe to store the results for the group
      group_results_df <- data.frame(
        Group = character(),
        Mean_First_Letter_Count = double(),
        SD_First_Letter_Count = double(),
        Mean_Last_Letter_Count = double(),
        SD_Last_Letter_Count = double(),
        Mean_Combined_Count = double(),
        Percentage = double(),
        SD_Combined_Count = double(),
        stringsAsFactors = FALSE
      )
      
      # Iterate over each CSV file in the group
      for (i in groups) {
        print(i)
        
        # Extract the first and last letters from the CSV file
        file_results <- extract_first_last_letters_from_csv(i)
        
        # Get the filename without the directory path
        filename <- basename(i)
        
        # Get the second string in the filename by splitting on "_"
        grp_str <- rep(strsplit(filename, "_")[[1]][2], nrow(file_results))
        
        # Store the label for the group
        file_results$group <- data.frame(grp_str)
        
        # Append the results to the dataframe
        dat <- rbind(dat, file_results)
      }
    }
    
    # Aggregate the results by letter and group
    agg_df <- dat %>%
      group_by(letter, group) %>%
      dplyr::summarize(
        mean_combined_count = mean(combined_count),
        sd_combined_count = sd(combined_count),
        se_combined_count = sd_combined_count / sqrt(n()),
        sum_combined_count = sum(combined_count),
        mean_first_letter_count = mean(first_letter_count),
        sd_first_letter_count = sd(first_letter_count),
        se_first_letter_count = sd_first_letter_count / sqrt(n()),
        sum_first_letter_count = sum(first_letter_count),
        mean_last_letter_count = mean(last_letter_count),
        sd_last_letter_count = sd(last_letter_count),
        se_last_letter_count = sd_last_letter_count / sqrt(n()),
        sum_last_letter_count = sum(last_letter_count)
      )
    
    # Calculate percentage of mean_Last_letter_count and mean_first_letter_count with SE
    agg_df_sum <- agg_df %>%
      group_by(group) %>%
      dplyr::summarize(
        total_sum_combined_count = sum(sum_combined_count),
        total_sum_first_letter_count = sum(sum_first_letter_count),
        total_sum_last_letter_count = sum(sum_last_letter_count)
      )
    
    agg_df <- agg_df %>%
      left_join(agg_df_sum, by = "group") %>%
      dplyr::mutate(
        perc_mean_first_letter_count = mean_first_letter_count / total_sum_first_letter_count * 100,
        perc_mean_last_letter_count = mean_last_letter_count / total_sum_last_letter_count * 100,
        SE_perc_mean_first_letter_count = sqrt(perc_mean_first_letter_count / 100 * (100 - perc_mean_first_letter_count) / total_sum_first_letter_count),
        SE_perc_mean_last_letter_count = sqrt(perc_mean_last_letter_count / 100 * (100 - perc_mean_last_letter_count) / total_sum_last_letter_count)
      )
    
    # Convert the grouped dataframe to a regular dataframe
    agg_df <- data.frame(agg_df)
    
    # Convert the group variable to a factor
    grp <- agg_df$group[[1]]
    grp <- as.factor(grp)
    agg_df$group <- grp
    
    # Create the appropriate density plot based on the PlotType argument
    
    # Nter Percentage Plot
    if (PlotType == "Nter") {
      cutsiteplots <- ggplot(agg_df, aes(x = letter, y = perc_mean_first_letter_count, fill = group)) +
        geom_bar(position = "dodge", stat = "identity") +
        geom_errorbar(
          aes(ymin = perc_mean_first_letter_count - SE_perc_mean_first_letter_count,
              ymax = perc_mean_first_letter_count + SE_perc_mean_first_letter_count),
          position = position_dodge(width = 0.9), width = 0.25
        ) +
        labs(title = "N-terminal cleavage sites",
             x = "", y = "Percentage") +
        theme_grey(base_size = 12) +
        theme(
          axis.text.x = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          title = element_text(size = 14, face = "bold")
        )
    }
    
    # Cter Percentage Plot
    else if (PlotType == "Cter") {
      cutsiteplots <- ggplot(agg_df, aes(x = letter, y = perc_mean_last_letter_count, fill = group)) +
        geom_bar(position = "dodge", stat = "identity") +
        geom_errorbar(
          aes(ymin = perc_mean_last_letter_count - SE_perc_mean_last_letter_count,
              ymax = perc_mean_last_letter_count + SE_perc_mean_last_letter_count),
          position = position_dodge(width = 0.9), width = 0.25
        ) +
        labs(title = "C-terminal cleavage sites",
             x = "", y = "Percentage") +
        theme_grey(base_size = 12) +
        theme(
          axis.text.x = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          title = element_text(size = 14, face = "bold")
        )
    }
    
    # Plot mean_first_letter_count and mean_last_letter_count as percentage with SE
    else if (PlotType == "Combined") {
      cutsiteplots <- ggplot(agg_df, aes(x = letter, y = perc_mean_first_letter_count, fill = group)) +
        geom_bar(position = "dodge", stat = "identity") +
        geom_errorbar(
          aes(ymin = perc_mean_first_letter_count - SE_perc_mean_first_letter_count,
              ymax = perc_mean_first_letter_count + SE_perc_mean_first_letter_count),
          width = 0.2, position = position_dodge(0.9)
        ) +
        ggtitle("Percentage of Mean First Letter Count by Group") +
        xlab("Letter") +
        ylab("Percentage of Mean First Letter Count") +
        theme_bw()
      
      ggplot(agg_df, aes(x = letter, y = perc_mean_last_letter_count, fill = group)) +
        geom_bar(position = "dodge", stat = "identity") +
        geom_errorbar(
          aes(ymin = perc_mean_last_letter_count - SE_perc_mean_last_letter_count,
              ymax = perc_mean_last_letter_count + SE_perc_mean_last_letter_count),
          width = 0.2, position = position_dodge(0.9)
        ) +
        ggtitle("Percentage of Mean Last Letter Count by Group") +
        xlab("Letter") +
        ylab("Percentage of Mean Last Letter Count") +
        theme_bw()
    }
    
    # If the PlotType argument is not valid
    else {
      stop("Error: Invalid PlotType argument. Valid options are 'Nter', 'Cter', and 'Combined'.")
    }
    
    # Return the plot
    print(cutsiteplots)
  }
  
  # Create the main window
  win <- tktoplevel()
  tkwm.title(win, "Cut sites distribution")
  
  # Create and configure the CSV directory selection frame
  csvDirFrame <- ttkframe(win)
  tkgrid(csvDirFrame, padx = 20, pady = 20)
  
  # Create and configure the CSV directory label
  csvDirLabel <- ttklabel(csvDirFrame, text = "CSV Directory:")
  tkgrid(csvDirLabel, column = 1, row = 1, sticky = "w")
  
  # Create and configure the CSV directory entry
  csvDirEntry <- ttkentry(csvDirFrame, textvariable = csv_directory)
  tkgrid(csvDirEntry, column = 2, row = 1)
  
  # Create and configure the CSV directory browse button
  csvDirButton <- ttkbutton(csvDirFrame, text = "Browse", command = selectCSVDir)
  tkgrid(csvDirButton, column = 3, row = 1, padx = 5)
  
  # Create and configure the plot type selection frame
  plotTypeFrame <- ttkframe(win)
  tkgrid(plotTypeFrame, padx = 10, pady = 10)
  
  # Create and configure the plot type label
  plotTypeLabel <- ttklabel(plotTypeFrame, text = "Plot Type:")
  tkgrid(plotTypeLabel, column = 1, row = 1, sticky = "w")
  
  # Create and configure the plot type radiobuttons
  plotTypeRadioFrame <- ttkframe(plotTypeFrame)
  tkgrid(plotTypeRadioFrame, column = 2, row = 1)
  
  # Create and configure the Overlay radiobutton
  overlayRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Cter", variable = PlotType, value = "Cter", command = function() setPlotType("Cter"))
  tkgrid(overlayRadio, column = 1, row = 1, sticky = "w")
  
  # Create and configure the Ridges radiobutton
  ridgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Nter", variable = PlotType, value = "Nter", command = function() setPlotType("Nter"))
  tkgrid(ridgesRadio, column = 1, row = 2, sticky = "w")
  
  # Create and configure the Colored Ridges radiobutton
  coloredRidgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Combined", variable = PlotType, value = "Combined", command = function() setPlotType("Combined"))
  tkgrid(coloredRidgesRadio, column = 1, row = 3, sticky = "w")
  
  # Create and configure the create plot button
  createPlotButton <- ttkbutton(win, text = "Create Plot", command = createBarPlot)
  tkgrid(createPlotButton, pady = 10)
  
  # Start the event loop
  tkwait.visibility(win)
  
}

# call the function
#csp ()

################################################################################
# Plotting cut sites GUI.

#' A Function from gui.R
#'
#' This function does something even more interesting.
#'
#' @export
cs <- function() {
  tclCheck()
  dlg <- myToplevel('cs')
  if (is.null(dlg))
    return(invisible())
  
  # Read protease data from CSV
  protCutSites <- read.csv('C:/Users/dimit/OneDrive/Bureau/DigestR/Proteasecutsiteslist.csv')
  
  # Create the main window
  mainWindow <- tktoplevel()
  tkwm.title(mainWindow, "Plot Cut Site")
  
  # Create protease label and dropdown menu
  proteaseLabel <- tklabel(mainWindow, text = "Select Protease:")
  tkgrid(proteaseLabel, padx = 10, pady = 10)
  
  proteaseDropdown <- tklistbox(mainWindow)
  tkgrid(proteaseDropdown, padx = 10, pady = 10)
  
  # Populate protease dropdown menu with protease names
  proteaseList <- unique(protCutSites$protease)
  for (prot in proteaseList) {
    tkinsert(proteaseDropdown, "end", prot)
  }
  
  # Create color label and entry field
  colorLabel <- tklabel(mainWindow, text = "Enter Color for Lines:")
  tkgrid(colorLabel, padx = 10, pady = 10)
  
  colorEntry <- tkentry(mainWindow)
  tkgrid(colorEntry, padx = 10, pady = 10)
  
  # Create plot button
  plotButton <- tkbutton(mainWindow, text = "Plot", command = function() {
    # Get the selected protease and color
    selectedProteaseIndex <- tclvalue(tkcurselection(proteaseDropdown))
    selectedProtease <- tkget(proteaseDropdown, selectedProteaseIndex)
    selectedColor <- tkget(colorEntry)
    
    # Call the plotCutSite function with the selected protease and color
    plotCutSite(selectedProtease, selectedColor, protCutSites, colorEntry)
    
    # Close the GUI window
    tkdestroy(mainWindow)
  })
  tkgrid(plotButton, padx = 10, pady = 10)
}

# Function to plot cut sites
plotCutSite <- function(prot, colour, protCutSites, colorEntry) {
  # Retrieve current gene information
  currGene <- globalSettings$geneDisp
  
  # Find the index of the current gene in the species data
  idx <- which(species$genes$name == currGene)
  
  if (length(idx) > 0) {
    # Retrieve the corresponding sequence based on the gene index
    start <- species$genes$seqStartIdx[idx]
    end <- start + species$genes$seqLength[idx]
    currSeq <- substr(species$seq, start, end)
  } else {
    stop("No current gene information found. Function doesn't work with proteome-wide view. Call a specific protein.")
  }
  
  # Set gene sequence and colour variables
  gene <- globalSettings$geneDisp
  seq <- currSeq
  colour <- tclvalue(tkget(colorEntry))
  
  # Find the index of the protease in the protCutSites data
  protIndex <- which(as.character(protCutSites$protease) == as.character(prot))
  col_num <- ncol(protCutSites)
  
  # Get the cut sites for the protease
  protase_cutsite <- protCutSites[protIndex, 3:col_num]
  loc <- 0
  
  # Iterate through the cut sites
  for (j in protase_cutsite) {
    if (j == "") {
      break
    }
    print('################################################################')
    print(j)
    print('################################################################')
    
    prev_index <- 0
    protcuts <- strsplit(seq, j)
    
    # Iterate through the sequence segments
    for (k in protcuts[[1]]) {
      k <- paste(k, j, sep = '')
      print(k)
      
      # Find the location of the segment in the sequence
      loc <- gregexpr(k, seq, fixed = TRUE)
      start <- loc[[1]][1]
      len <- nchar(k)
      
      # Check if the position is greater than 0 and the segment length is greater than 0
      if (start > 0 && len > 0) {
        # Calculate the plot site
        plot_site <- start + len - 1
        
        if (k == "K") {
          # Check if it is the last character in the sequence
          if (plot_site == nchar(seq)) {
            plot_site <- prev_index
          } else {
            plot_site <- prev_index + 1
          }
          print('IS IT HERE?')
        }
        
        print(plot_site)
        
        # Assuming you have initialized a plot or set up the plotting device
        abline(v = plot_site, col = colour, lty = 2)
        
        
        prev_index <- plot_site
      }
    }
  }
  return(invisible(prot))
}

# Call the GUI handler function to start the application
#cs()

################################################################################
# rd <- function() {
#   ## creates main window
#   tclCheck()
#   dlg <- tktoplevel()
#   tkwm.title(dlg, 'Prepare Mascot file for DigestR')
#   tkfocus(dlg)
#   tkwm.deiconify(dlg)
#   
#   ## create Remove Duplicates button
#   onRemoveDuplicates <- function() {
#     source_file <- tclvalue(tcl("tk_getOpenFile", "-initialdir", ".", "-filetypes", "{{CSV Files} {.csv}}"))
#     # Call the RemoveDuplicate() function passing the selected file
#     result_file <- removeDuplicates(source_file)
#     # Perform actions with the result_file, such as displaying or opening it
#   }
#   removeDuplicatesButton <- ttkbutton(dlg, text='Remove Duplicates', width=15, command=onRemoveDuplicates)
#   
#   ## add button to the main window
#   tkgrid(removeDuplicatesButton, column=1, row=1, padx=5, pady=5)
#   
#   ## make main window stretch when resized
#   tkgrid.columnconfigure(dlg, 1, weight=1)
#   tkgrid.rowconfigure(dlg, 1, weight=1)
#   
#   tkwait.window(dlg)
# }
# 
# removeDuplicates <- function(input_file) {
#   # Read the input CSV file
#   data <- read.csv(input_file, skip = 3, header = TRUE)
#   
#   # Check if 'pep_seq' column exists
#   if (!("pep_seq" %in% colnames(data))) {
#     stop("Error: 'pep_seq' column not found in the input file.")
#   }
#   
#   # Remove duplicate values in the pep_seq column
#   rd_data <- data[!duplicated(data$pep_seq), ]
#   
#   # Prompt the user to select the destination CSV file
#   dest_file <- tclvalue(tcl("tk_getSaveFile", "-initialdir", ".", "-defaultextension", ".csv", "-filetypes", "{{CSV Files} {.csv}}"))
#   
#   # Write the unique data to the destination CSV file
#   write.csv(rd_data, file = dest_file, row.names = FALSE)
#   
#   # Return the path of the destination CSV file
#   return(dest_file)
# }
# 
# # Test the RemoveDuplicate() function
# #rd()
# 
# ################################################################################
# # Unique peptides function  GUI 
# ## Change code so that it uses the current directory to export the file to.
# 
# up <- function() {
#   ## creates main window
#   tclCheck()
#   dlg <- tktoplevel()
#   tkwm.title(dlg, 'Unique peptides')
#   tkfocus(dlg)
#   tkwm.deiconify(dlg)
#   
#   ## create select Query button
#   L1Var <- tclVar()  # Variable to store the path of the selected L1 file
#   onSelectL1 <- function() {
#     selectedL1 <- tclvalue(tcl("tk_getOpenFile", "-initialdir", ".", "-filetypes", "{{CSV Files} {.csv}}"))
#     tclvalue(L1Var) <- selectedL1
#   }
#   selectL1Button <- ttkbutton(dlg, text='Query', width=15, command=onSelectL1)
#   
#   ## create select L2 button
#   L2Var <- tclVar()  # Variable to store the path of the selected L2 file
#   onSelectL2 <- function() {
#     selectedL2 <- tclvalue(tcl("tk_getOpenFile", "-initialdir", ".", "-filetypes", "{{CSV Files} {.csv}}"))
#     tclvalue(L2Var) <- selectedL2
#   }
#   selectL2Button <- ttkbutton(dlg, text='Experimental', width=15, command=onSelectL2)
#   
#   ## create Remove Duplicates button
#   onuniquePeptides <- function() {
#     L1 <- read.csv(tclvalue(L1Var), header=TRUE)
#     L2_file <- tclvalue(L2Var)
#     L2 <- read.csv(L2_file, header=TRUE)
#     
#     # Check if 'pep_seq' column exists in L1 and L2
#     if (!("pep_seq" %in% colnames(L1)) || !("pep_seq" %in% colnames(L2))) {
#       stop("Error: 'pep_seq' column not found in the input files.")
#     }
#     
#     # Find unique peptides in L2 not present in L1
#     L2_uniques <- L2[!(L2$pep_seq %in% L1$pep_seq), ]
#     
#     # Get the name of the L2 file without extension
#     l2_filename <- tools::file_path_sans_ext(basename(L2_file))
#     
#     # Create the new filename with "_unique" attached
#     dest_file <- paste0(l2_filename, "_Unique.csv")
#     
#     # Write the unique data to the destination CSV file
#     write.csv(L2_uniques, file = dest_file, row.names = FALSE)
#     
#     tkmessageBox(message=paste("Unique peptides saved to:", dest_file))
#   }
#   uniquePeptidesButton <- ttkbutton(dlg, text='Experimental unique', width=20, command=onuniquePeptides)
#   
#   ## add widgets to the main window
#   tkgrid(selectL1Button, column=1, row=1, padx=5, pady=5)
#   tkgrid(selectL2Button, column=1, row=2, padx=5, pady=5)
#   tkgrid(uniquePeptidesButton, column=1, row=3, padx=5, pady=5)
#   
#   ## make the main window stretch when resized
#   tkgrid.columnconfigure(dlg, 1, weight=1)
#   tkgrid.rowconfigure(dlg, 1, weight=1)
#   
#   tkwait.window(dlg)
# }
# 
# # Test the up() function
# #up()
# 
# ################################################################################
# # Density plot function GUI 
# pd <- function() {
#   # Define variables
#   csv_dir <- ""
#   plot_type <- ""
#   
#   # Function to select the CSV directory
#   selectCSVDir <- function() {
#     csv_dir <<- tclvalue(tcl("tk_chooseDirectory"))
#     csv_dir <- tclVar(csv_dir)  # Update the csv_dir variable
#     tkconfigure(csvDirEntry, text = csv_dir)  # Update the csvDirEntry widget
#   }
#   
#   # Function to set the selected plot type
#   setPlotType <- function(value) {
#     plot_type <<- value
#   }
#   
#   # Function to create the density plot
#   createDensityPlot <- function() {
#     grouped_df <- groupCSVFiles(csv_dir)
#     
#     # Create an empty data frame to store the nchar values for all groups
#     nchar_df <- data.frame(group_name = character(), nchar = integer())
#     
#     for (i in 1:length(grouped_df)) {
#       group_name <- names(grouped_df)[i]
#       pep_seq_col <- grouped_df[[i]]$pep_seq
#       group_nchar_df <- data.frame(group_name = group_name, nchar = nchar(pep_seq_col))
#       nchar_df <- rbind(nchar_df, group_nchar_df)
#     }
#     
#     nchar_mean_df <- nchar_df %>%
#       group_by(group_name) %>%
#       summarise(mean_nchar = mean(nchar))
#     
#     if (plot_type == "Overlay") {
#       plot <- ggplot(nchar_df, aes(x = nchar, group = group_name, fill = group_name)) +
#         geom_density(alpha = 0.5) +
#         labs(title = "Peptide Length Distribution", x = "Peptide length in AA", y = "", fill = "Group") +
#         theme_bw() +
#         geom_vline(data = nchar_mean_df, aes(xintercept = mean_nchar, color = group_name), size = 1)
#     } else if (plot_type == "Ridges") {
#       plot <- ggplot(nchar_df, aes(x = nchar, y = group_name, fill = group_name)) +
#         geom_density_ridges2(alpha = 0.5, rel_min_height = 0.01, scale = 7, quantile_lines = TRUE, quantile_fun = function(x, ...) mean(x)) +
#         labs(title = "Peptide Length Distribution", x = "Peptide length in AA", y = "", fill = "Group") +
#         theme_bw()
#     } else if (plot_type == "Colored Ridges") {
#       plot <- ggplot(nchar_df, aes(x = nchar, y = group_name, fill = stat(x))) +
#         geom_density_ridges_gradient(alpha = 0.5, rel_min_height = 0.01, scale = 3, quantile_lines = TRUE, quantile_fun = function(x, ...) mean(x)) +
#         labs(title = "Peptide Length Distribution", x = "Peptide length in AA", y = "", fill = "Group") +
#         scale_fill_viridis_c(name = "Length", option = "H")
#     } else {
#       stop("Invalid plot_type argument. Please choose either 'Overlay', 'Ridges', or 'Colored Ridges'.")
#     }
#     
#     print(plot)
#   }
#   
#   # Create the main window
#   win <- tktoplevel()
#   tkwm.title(win, "Density Plot Creator")
#   
#   # Create and configure the CSV directory selection frame
#   csvDirFrame <- ttkframe(win)
#   tkgrid(csvDirFrame, padx = 20, pady = 20)
#   
#   # Create and configure the CSV directory label
#   csvDirLabel <- ttklabel(csvDirFrame, text = "CSV Directory:")
#   tkgrid(csvDirLabel, column = 1, row = 1, sticky = "w")
#   
#   # Create and configure the CSV directory entry
#   csvDirEntry <- ttkentry(csvDirFrame, textvariable = csv_dir)
#   tkgrid(csvDirEntry, column = 2, row = 1)
#   
#   # Create and configure the CSV directory browse button
#   csvDirButton <- ttkbutton(csvDirFrame, text = "Browse", command = selectCSVDir)
#   tkgrid(csvDirButton, column = 3, row = 1, padx = 5)
#   
#   # Create and configure the plot type selection frame
#   plotTypeFrame <- ttkframe(win)
#   tkgrid(plotTypeFrame, padx = 10, pady = 10)
#   
#   # Create and configure the plot type label
#   plotTypeLabel <- ttklabel(plotTypeFrame, text = "Plot Type:")
#   tkgrid(plotTypeLabel, column = 1, row = 1, sticky = "w")
#   
#   # Create and configure the plot type radiobuttons
#   plotTypeRadioFrame <- ttkframe(plotTypeFrame)
#   tkgrid(plotTypeRadioFrame, column = 2, row = 1)
#   
#   # Create and configure the Overlay radiobutton
#   overlayRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Overlay", variable = plot_type, value = "Overlay", command = function() setPlotType("Overlay"))
#   tkgrid(overlayRadio, column = 1, row = 1, sticky = "w")
#   
#   # Create and configure the Ridges radiobutton
#   ridgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Ridges", variable = plot_type, value = "Ridges", command = function() setPlotType("Ridges"))
#   tkgrid(ridgesRadio, column = 1, row = 2, sticky = "w")
#   
#   # Create and configure the Colored Ridges radiobutton
#   coloredRidgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Colored Ridges", variable = plot_type, value = "Colored Ridges", command = function() setPlotType("Colored Ridges"))
#   tkgrid(coloredRidgesRadio, column = 1, row = 3, sticky = "w")
#   
#   # Create and configure the create plot button
#   createPlotButton <- ttkbutton(win, text = "Create Plot", command = createDensityPlot)
#   tkgrid(createPlotButton, pady = 10)
#   
#   # Start the event loop
#   tkwait.visibility(win)
#   
# }
# 
# #call the function
# #pd()
# 
# ################################################################################ 
# # Cter, Nter function GUI 
# csp <- function () {
#   
#   # Define variables
#   csv_directory <- ""
#   PlotType <- ""
#   
#   # Function to select the CSV directory
#   selectCSVDir <- function() {
#     csv_directory <<- tclvalue(tcl("tk_chooseDirectory"))
#     csv_directory <- tclVar(csv_directory)  # Update the csv_directory variable
#     tkconfigure(csvDirEntry, text = csv_directory)  # Update the csvDirEntry widget
#   }
#   
#   # Function to set the selected plot type
#   setPlotType <- function(value) {
#     PlotType <<- value
#   }
#   
#   # Function to create the density plot
#   createBarPlot <- function() 
#     #plotCutSites <- function(csv_directory, PlotType) 
#   {
#     library(dplyr)
#     library(ggplot2)
#     
#     extract_first_last_letters_from_csv <- function(file_path) {
#       # Read the CSV file skipping the first 3 lines (assuming headers are in line 4)
#       data <- read.csv(file_path, skip = 0, header = TRUE)
#       
#       # Check if the 'pep_seq' column exists
#       if (!("pep_seq" %in% colnames(data))) {
#         stop("Error: 'pep_seq' column not found in the CSV file.")
#       }
#       
#       # Extract the first and last letters from each sequence
#       data$first_letter <- substr(data$pep_seq, 1, 1)
#       data$last_letter <- substr(data$pep_seq, nchar(data$pep_seq), nchar(data$pep_seq))
#       
#       # Count the occurrences of each first letter
#       first_letter_counts <- table(factor(data$first_letter, levels = c("A", "C", "D", "E", "F", "G","H" ,"K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")))
#       
#       # Count the occurrences of each last letter
#       last_letter_counts <- table(factor(data$last_letter, levels = c("A", "C", "D", "E", "F", "G","H", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")))
#       
#       # Combine the counts for first and last letters
#       combined_counts <- merge(
#         data.frame(letter = names(first_letter_counts), first_letter_count = as.numeric(first_letter_counts)),
#         data.frame(letter = names(last_letter_counts), last_letter_count = as.numeric(last_letter_counts)),
#         by = "letter", all = TRUE
#       )
#       combined_counts$combined_count <- rowSums(combined_counts[, c("first_letter_count", "last_letter_count")], na.rm = TRUE)
#       
#       # Replace NA values with 0
#       combined_counts[is.na(combined_counts)] <- 0
#       
#       # Return the combined counts
#       return(combined_counts)
#     }
#     
#     group_csv_files <- function(csv_directory) {
#       # Get the list of CSV files in the directory
#       csv_files <- list.files(path = csv_directory, pattern = "*.csv", full.names = TRUE)
#       
#       # Create an empty list to store the grouped CSV files
#       grouped_csv <- list()
#       
#       # Loop through the CSV files and group them based on the second string
#       for (i in 1:length(csv_files)) {
#         # Get the filename without the directory path
#         filename <- basename(csv_files[i])
#         # Get the second string in the filename by splitting on "_"
#         second_str <- strsplit(filename, "_")[[1]][2]
#         # If the group already exists in the list, append the CSV file
#         if (second_str %in% names(grouped_csv)) {
#           grouped_csv[[second_str]] <- c(grouped_csv[[second_str]], csv_files[i])
#         }
#         # If the group doesn't exist in the list, create a new list element
#         else {
#           grouped_csv[[second_str]] <- list(csv_files[i])
#         }
#       }
#       
#       # Return the list of grouped CSV files
#       return(grouped_csv)
#     }
#     
#     # Group the CSV files
#     csv <- group_csv_files(csv_directory)
#     
#     # Initialize an empty dataframe to store the results
#     dat <- data.frame()
#     
#     # Iterate over each group of CSV files
#     for (groups in csv) {
#       
#       # Create an empty dataframe to store the results for the group
#       group_results_df <- data.frame(
#         Group = character(),
#         Mean_First_Letter_Count = double(),
#         SD_First_Letter_Count = double(),
#         Mean_Last_Letter_Count = double(),
#         SD_Last_Letter_Count = double(),
#         Mean_Combined_Count = double(),
#         Percentage = double(),
#         SD_Combined_Count = double(),
#         stringsAsFactors = FALSE
#       )
#       
#       # Iterate over each CSV file in the group
#       for (i in groups) {
#         print(i)
#         
#         # Extract the first and last letters from the CSV file
#         file_results <- extract_first_last_letters_from_csv(i)
#         
#         # Get the filename without the directory path
#         filename <- basename(i)
#         
#         # Get the second string in the filename by splitting on "_"
#         grp_str <- rep(strsplit(filename, "_")[[1]][2], nrow(file_results))
#         
#         # Store the label for the group
#         file_results$group <- data.frame(grp_str)
#         
#         # Append the results to the dataframe
#         dat <- rbind(dat, file_results)
#       }
#     }
#     
#     # Aggregate the results by letter and group
#     agg_df <- dat %>%
#       group_by(letter, group) %>%
#       dplyr::summarize(
#         mean_combined_count = mean(combined_count),
#         sd_combined_count = sd(combined_count),
#         se_combined_count = sd_combined_count / sqrt(n()),
#         sum_combined_count = sum(combined_count),
#         mean_first_letter_count = mean(first_letter_count),
#         sd_first_letter_count = sd(first_letter_count),
#         se_first_letter_count = sd_first_letter_count / sqrt(n()),
#         sum_first_letter_count = sum(first_letter_count),
#         mean_last_letter_count = mean(last_letter_count),
#         sd_last_letter_count = sd(last_letter_count),
#         se_last_letter_count = sd_last_letter_count / sqrt(n()),
#         sum_last_letter_count = sum(last_letter_count)
#       )
#     
#     # Calculate percentage of mean_Last_letter_count and mean_first_letter_count with SE
#     agg_df_sum <- agg_df %>%
#       group_by(group) %>%
#       dplyr::summarize(
#         total_sum_combined_count = sum(sum_combined_count),
#         total_sum_first_letter_count = sum(sum_first_letter_count),
#         total_sum_last_letter_count = sum(sum_last_letter_count)
#       )
#     
#     agg_df <- agg_df %>%
#       left_join(agg_df_sum, by = "group") %>%
#       dplyr::mutate(
#         perc_mean_first_letter_count = mean_first_letter_count / total_sum_first_letter_count * 100,
#         perc_mean_last_letter_count = mean_last_letter_count / total_sum_last_letter_count * 100,
#         SE_perc_mean_first_letter_count = sqrt(perc_mean_first_letter_count / 100 * (100 - perc_mean_first_letter_count) / total_sum_first_letter_count),
#         SE_perc_mean_last_letter_count = sqrt(perc_mean_last_letter_count / 100 * (100 - perc_mean_last_letter_count) / total_sum_last_letter_count)
#       )
#     
#     # Convert the grouped dataframe to a regular dataframe
#     agg_df <- data.frame(agg_df)
#     
#     # Convert the group variable to a factor
#     grp <- agg_df$group[[1]]
#     grp <- as.factor(grp)
#     agg_df$group <- grp
#     
#     # Create the appropriate density plot based on the PlotType argument
#     
#     # Nter Percentage Plot
#     if (PlotType == "Nter") {
#       cutsiteplots <- ggplot(agg_df, aes(x = letter, y = perc_mean_first_letter_count, fill = group)) +
#         geom_bar(position = "dodge", stat = "identity") +
#         geom_errorbar(
#           aes(ymin = perc_mean_first_letter_count - SE_perc_mean_first_letter_count,
#               ymax = perc_mean_first_letter_count + SE_perc_mean_first_letter_count),
#           position = position_dodge(width = 0.9), width = 0.25
#         ) +
#         labs(title = "N-terminal cleavage sites",
#              x = "", y = "Percentage") +
#         theme_grey(base_size = 12) +
#         theme(
#           axis.text.x = element_text(size = 14, face = "bold"),
#           axis.text.y = element_text(size = 12, face = "bold"),
#           title = element_text(size = 14, face = "bold")
#         )
#     }
#     
#     # Cter Percentage Plot
#     else if (PlotType == "Cter") {
#       cutsiteplots <- ggplot(agg_df, aes(x = letter, y = perc_mean_last_letter_count, fill = group)) +
#         geom_bar(position = "dodge", stat = "identity") +
#         geom_errorbar(
#           aes(ymin = perc_mean_last_letter_count - SE_perc_mean_last_letter_count,
#               ymax = perc_mean_last_letter_count + SE_perc_mean_last_letter_count),
#           position = position_dodge(width = 0.9), width = 0.25
#         ) +
#         labs(title = "C-terminal cleavage sites",
#              x = "", y = "Percentage") +
#         theme_grey(base_size = 12) +
#         theme(
#           axis.text.x = element_text(size = 14, face = "bold"),
#           axis.text.y = element_text(size = 12, face = "bold"),
#           title = element_text(size = 14, face = "bold")
#         )
#     }
#     
#     # Plot mean_first_letter_count and mean_last_letter_count as percentage with SE
#     else if (PlotType == "Combined") {
#       cutsiteplots <- ggplot(agg_df, aes(x = letter, y = perc_mean_first_letter_count, fill = group)) +
#         geom_bar(position = "dodge", stat = "identity") +
#         geom_errorbar(
#           aes(ymin = perc_mean_first_letter_count - SE_perc_mean_first_letter_count,
#               ymax = perc_mean_first_letter_count + SE_perc_mean_first_letter_count),
#           width = 0.2, position = position_dodge(0.9)
#         ) +
#         ggtitle("Percentage of Mean First Letter Count by Group") +
#         xlab("Letter") +
#         ylab("Percentage of Mean First Letter Count") +
#         theme_bw()
#       
#       ggplot(agg_df, aes(x = letter, y = perc_mean_last_letter_count, fill = group)) +
#         geom_bar(position = "dodge", stat = "identity") +
#         geom_errorbar(
#           aes(ymin = perc_mean_last_letter_count - SE_perc_mean_last_letter_count,
#               ymax = perc_mean_last_letter_count + SE_perc_mean_last_letter_count),
#           width = 0.2, position = position_dodge(0.9)
#         ) +
#         ggtitle("Percentage of Mean Last Letter Count by Group") +
#         xlab("Letter") +
#         ylab("Percentage of Mean Last Letter Count") +
#         theme_bw()
#     }
#     
#     # If the PlotType argument is not valid
#     else {
#       stop("Error: Invalid PlotType argument. Valid options are 'Nter', 'Cter', and 'Combined'.")
#     }
#     
#     # Return the plot
#     return(cutsiteplots)
#   }
#   
#   # Create the main window
#   win <- tktoplevel()
#   tkwm.title(win, "Cut sites distribution")
#   
#   # Create and configure the CSV directory selection frame
#   csvDirFrame <- ttkframe(win)
#   tkgrid(csvDirFrame, padx = 20, pady = 20)
#   
#   # Create and configure the CSV directory label
#   csvDirLabel <- ttklabel(csvDirFrame, text = "CSV Directory:")
#   tkgrid(csvDirLabel, column = 1, row = 1, sticky = "w")
#   
#   # Create and configure the CSV directory entry
#   csvDirEntry <- ttkentry(csvDirFrame, textvariable = csv_directory)
#   tkgrid(csvDirEntry, column = 2, row = 1)
#   
#   # Create and configure the CSV directory browse button
#   csvDirButton <- ttkbutton(csvDirFrame, text = "Browse", command = selectCSVDir)
#   tkgrid(csvDirButton, column = 3, row = 1, padx = 5)
#   
#   # Create and configure the plot type selection frame
#   plotTypeFrame <- ttkframe(win)
#   tkgrid(plotTypeFrame, padx = 10, pady = 10)
#   
#   # Create and configure the plot type label
#   plotTypeLabel <- ttklabel(plotTypeFrame, text = "Plot Type:")
#   tkgrid(plotTypeLabel, column = 1, row = 1, sticky = "w")
#   
#   # Create and configure the plot type radiobuttons
#   plotTypeRadioFrame <- ttkframe(plotTypeFrame)
#   tkgrid(plotTypeRadioFrame, column = 2, row = 1)
#   
#   # Create and configure the Overlay radiobutton
#   overlayRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Cter", variable = PlotType, value = "Cter", command = function() setPlotType("Cter"))
#   tkgrid(overlayRadio, column = 1, row = 1, sticky = "w")
#   
#   # Create and configure the Ridges radiobutton
#   ridgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Nter", variable = PlotType, value = "Nter", command = function() setPlotType("Nter"))
#   tkgrid(ridgesRadio, column = 1, row = 2, sticky = "w")
#   
#   # Create and configure the Colored Ridges radiobutton
#   coloredRidgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Combined", variable = PlotType, value = "Combined", command = function() setPlotType("Combined"))
#   tkgrid(coloredRidgesRadio, column = 1, row = 3, sticky = "w")
#   
#   # Create and configure the create plot button
#   createPlotButton <- ttkbutton(win, text = "Create Plot", command = createBarPlot)
#   tkgrid(createPlotButton, pady = 10)
#   
#   # Start the event loop
#   tkwait.visibility(win)
#   
# }
# 
# # call the function
# #csp()
# 
# ################################################################################
# # Plotting cut sites GUI.
# 
# cs <- function() {
#   # Read protease data from CSV
#   protCutSites <- read.csv('C:/Users/dimit/OneDrive/Bureau/DigestR/Proteasecutsiteslist.csv')
#   
#   # Create the main window
#   mainWindow <- tktoplevel()
#   tkwm.title(mainWindow, "Plot Cut Site")
#   
#   # Create protease label and dropdown menu
#   proteaseLabel <- tklabel(mainWindow, text = "Select Protease:")
#   tkgrid(proteaseLabel, padx = 10, pady = 10)
#   
#   proteaseDropdown <- tklistbox(mainWindow)
#   tkgrid(proteaseDropdown, padx = 10, pady = 10)
#   
#   # Populate protease dropdown menu with protease names
#   proteaseList <- unique(protCutSites$protease)
#   for (prot in proteaseList) {
#     tkinsert(proteaseDropdown, "end", prot)
#   }
#   
#   # Create color label and entry field
#   colorLabel <- tklabel(mainWindow, text = "Enter Color for Lines:")
#   tkgrid(colorLabel, padx = 10, pady = 10)
#   
#   colorEntry <- tkentry(mainWindow)
#   tkgrid(colorEntry, padx = 10, pady = 10)
#   
#   # Create plot button
#   plotButton <- tkbutton(mainWindow, text = "Plot", command = function() {
#     # Get the selected protease and color
#     selectedProteaseIndex <- tclvalue(tkcurselection(proteaseDropdown))
#     selectedProtease <- tkget(proteaseDropdown, selectedProteaseIndex)
#     selectedColor <- tkget(colorEntry)
#     
#     # Call the plotCutSite function with the selected protease and color
#     plotCutSite(selectedProtease, selectedColor, protCutSites, colorEntry)
#     
#     # Close the GUI window
#     tkdestroy(mainWindow)
#   })
#   tkgrid(plotButton, padx = 10, pady = 10)
# }
# 
# # Function to plot cut sites
# plotCutSite <- function(prot, colour, protCutSites, colorEntry) {
#   # Retrieve current gene information
#   currGene <- globalSettings$geneDisp
#   
#   # Find the index of the current gene in the species data
#   idx <- which(species$genes$name == currGene)
#   
#   if (length(idx) > 0) {
#     # Retrieve the corresponding sequence based on the gene index
#     start <- species$genes$seqStartIdx[idx]
#     end <- start + species$genes$seqLength[idx]
#     currSeq <- substr(species$seq, start, end)
#   } else {
#     stop("No current gene information found. Function doesn't work with proteome-wide view. Call a specific protein.")
#   }
#   
#   # Set gene sequence and colour variables
#   gene <- globalSettings$geneDisp
#   seq <- currSeq
#   colour <- tclvalue(tkget(colorEntry))
#   
#   # Find the index of the protease in the protCutSites data
#   protIndex <- which(as.character(protCutSites$protease) == as.character(prot))
#   col_num <- ncol(protCutSites)
#   
#   # Get the cut sites for the protease
#   protase_cutsite <- protCutSites[protIndex, 3:col_num]
#   loc <- 0
#   
#   # Iterate through the cut sites
#   for (j in protase_cutsite) {
#     if (j == "") {
#       break
#     }
#     print('################################################################')
#     print(j)
#     print('################################################################')
#     
#     prev_index <- 0
#     protcuts <- strsplit(seq, j)
#     
#     # Iterate through the sequence segments
#     for (k in protcuts[[1]]) {
#       k <- paste(k, j, sep = '')
#       print(k)
#       
#       # Find the location of the segment in the sequence
#       loc <- gregexpr(k, seq, fixed = TRUE)
#       start <- loc[[1]][1]
#       len <- nchar(k)
#       
#       # Check if the position is greater than 0 and the segment length is greater than 0
#       if (start > 0 && len > 0) {
#         # Calculate the plot site
#         plot_site <- start + len - 1
#         
#         if (k == "K") {
#           # Check if it is the last character in the sequence
#           if (plot_site == nchar(seq)) {
#             plot_site <- prev_index
#           } else {
#             plot_site <- prev_index + 1
#           }
#           print('IS IT HERE?')
#         }
#         
#         print(plot_site)
#         
#         # Assuming you have initialized a plot or set up the plotting device
#         abline(v = plot_site, col = colour, lty = 2)
#         
#         
#         prev_index <- plot_site
#       }
#     }
#   }
#   return(invisible(prot))
# }
# 
# # Call the GUI handler function to start the application
# #cs()
