
################################################################################
################################################################################
##                                                                            ##
##                                                                            ##
## DigestR version 1.0.0, Tools for viewing and analyzing protein catabolism. ##
##      Copyright (C) 2024, Dimitri Desmonts de Lamache, Raied Aburashed,     ##
##               Travis A. Bingemann, Sören Wacker and Ian A. Lewis           ##
##				under GPL-3                                   ##
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


################################################################################
##                                                                            ##
##  Internal functions for creating, saving and updating digestR objects      ##
##                                                                            ##
################################################################################

library(tcltk)
library(tcltk2)
library(Rcpp)
library(ggplot2)
library(ggridges)
library(dplyr)
library(magrittr)
library(biomaRt)
library(ggrepel)
library(ggseqlogo)


## Assigns objects to the global environment and creates an undo point
myAssign <- function(in.name = NULL, in.object, save.backup = TRUE){
  
  ## Make sure the file name is the correct format
  if(!is.character(in.name)) 
    stop( 'myAssign requires a character input', call. = FALSE)
  if( in.name == 'currentSpectrum' )
    in.object <- in.object[1]
  
  ## Assign object to global environment
  if(in.name != 'zoom')
    assign(in.name, in.object, inherits=FALSE, envir=.GlobalEnv)
  else  
    assign("fileFolder", in.object, inherits=FALSE, envir=.GlobalEnv)
  
  ## Backup copies of the environment for undo/redo 
  if( save.backup && (in.name == "fileFolder" || in.name == "roiTable" ||
                      in.name == "currentSpectrum" || in.name == "roiSummary" || 
                      in.name == "overlayList" || in.name == "zoom" || 
                      in.name == "globalSettings" )){
    
    ## Assign NA to non existing files as a structure placeholder
    if(!exists("fileFolder") || is.null(fileFolder))
      fileFolder <- list(NA, NA)
    if(!exists("roiTable") || is.null(roiTable))
      roiTable <- list(NA, NA)
    if(!exists("currentSpectrum") || is.null(currentSpectrum))
      currentSpectrum <- list(NA)
    if(!exists("roiSummary") || is.null(roiSummary))  
      roiSummary <- list(NA, NA)
    if(!exists("overlayList") || is.null(overlayList))
      overlayList <- list(NA, NA)
    if(!exists("globalSettings") || is.null(globalSettings))
      globalSettings <- list(NA, NA)
    
    ## Assign the current global state of the environment to the undo list
    if(!exists('oldFolder'))
      oldFolder <- list( undo.index = 0, assign.index = 0, fileFolder = NULL, 
                         roiTable = NULL, currentSpectrum = NULL, 
                         roiSummary = NULL, overlayList = NULL,
                         zoom.history = NULL, zoom.list = NULL)
    if(is.null(oldFolder$undo.index))
      oldFolder$undo.index <- 0
    oldFolder$undo.index <- oldFolder$undo.index + 1  
    if(is.null(oldFolder$assign.index))
      oldFolder$assign.index <- 0
    oldFolder$assign.index <- oldFolder$assign.index + 1  
    oldFolder$fileFolder[[oldFolder$undo.index]] <- fileFolder  
    oldFolder$roiTable[[oldFolder$undo.index]] <- roiTable
    oldFolder$currentSpectrum[[oldFolder$undo.index]] <- currentSpectrum  
    oldFolder$roiSummary[[oldFolder$undo.index]] <-
      roiSummary 
    oldFolder$overlayList[[oldFolder$undo.index]] <- overlayList
    oldFolder$globalSettings[[oldFolder$undo.index]] <- globalSettings
    
    ## Keep a seperate zoom list for zoom previous command
    if( in.name == 'currentSpectrum' )
    {
      oldFolder$zoom.list <- list()
      oldFolder$zoom.list[[1]] <- fileFolder[[currentSpectrum]]$graphics.par$usr
    }
    if((in.name == 'zoom' || is.null(oldFolder$zoom.list)) && 
       (!is.na(fileFolder[1]) && !is.na(currentSpectrum[1])) )
    {
      oldFolder$zoom.list[[(length(oldFolder$zoom.list) + 1)]] <- 
        fileFolder[[currentSpectrum]]$graphics.par$usr
      oldFolder$zoom.history[[oldFolder$undo.index]] <- TRUE  
    }else
      oldFolder$zoom.history[[oldFolder$undo.index]] <- FALSE      
    
    ## Trim oldFolder after undo
    if(oldFolder$undo.index < length(oldFolder$fileFolder))
    {
      for(i in 1:length(oldFolder)){
        if(names(oldFolder)[i] != 'undo.index' && 
           names(oldFolder)[i] != 'assign.index' &&
           names(oldFolder)[i] != 'zoom.list') 
          oldFolder[[i]] <- oldFolder[[i]][1:oldFolder$undo.index]
      }      
    }
    
    ## Limit oldFolder to 10 entries
    if(oldFolder$undo.index > 10){
      for(i in 1:length(oldFolder)){
        if(names(oldFolder)[i] != 'undo.index'&& 
           names(oldFolder)[i] != 'assign.index')
          oldFolder[[i]] <- rev(rev(oldFolder[[i]])[1:10])
      }
      oldFolder$undo.index <- 10  
    }
    
    ## Save changes to oldFolder
    assign("oldFolder", oldFolder, inherits=FALSE, envir=.GlobalEnv) 
    
    ## Save a backup copy of the workspace
    if (defaultSettings$autoBackup && (oldFolder$assign.index == 1 || 
                                       !oldFolder$assign.index %% 10))
    {
      cat('\nPerforming automatic backup . . . ')
      tryCatch(invisible(save(list=ls(envir=.GlobalEnv, all.names=TRUE), 
                              file=file.path('~', '.digestRbackup'), version=NULL, ascii=FALSE, 
                              compress=FALSE, envir=.GlobalEnv, eval.promises=FALSE, 
                              precheck=FALSE)), error=function(){})
      cat('complete\n')
    }
  }         
}

## Internal function for checking values in defaultSettings
## newDef - object to check
## returns newDef with the correct formatting, all invalid entries will be set
##	to their default values
checkDef <- function(newDef){
  if (missing(newDef))
    newDef <- defaultSettings
  
  ##ensure newDef is not missing any values
  defSet <- createObj('defaultSettings', returnObj=TRUE)
  defLength <- length(defSet)
  defNames <- names(defSet)
  newNames <- names(newDef)
  if (length(newDef) < defLength){
    missingNames <- which(is.na(match(defNames, newNames)))
    newDef <- c(newDef, defSet[missingNames])
  }
  
  ##reformat values
  for (i in defNames){
    if (length(grep('pch$', i))){
      if (!is.na(suppressWarnings(as.numeric(newDef[[i]]))))
        defMode <- 'numeric'
      else
        defMode <- 'character'
    }
    else if (length(grep('tck$', i)))
      defMode <- 'numeric'
    else
      defMode <- storage.mode(defSet[[i]]) 
    if (defMode != 'function')
      tryCatch(suppressWarnings(storage.mode(newDef[[i]]) <- defMode), 
               error=function(er) newDef[[i]] <- NULL)
    
    ##ensure values are the correct length
    if (!i %in% c('libLocs', 'searchLibs') && 
        length(newDef[[i]]) != length(defSet[[i]]))
      newDef[[i]] <- defSet[[i]]
    
    ##check for NA values in digestR specific settings
    digestRnames <- defNames[66:length(defSet)]
    digestRnames <- digestRnames[-match(c('xtck', 'ytck'), digestRnames)]
    if (i %in% digestRnames && suppressWarnings(is.na(newDef[[i]])))
      newDef[[i]] <- defSet[[i]]
    
    ##check for valid colors
    colorNames <- defNames[grep('color', defNames)]
    if (i %in% colorNames){
      colTest <- try(col2rgb(newDef[[i]]), silent=TRUE)
      if (class(colTest) == 'try-error')
        newDef[[i]] <- defSet[[i]]
    }
  }	
  
  ##check for invalid "par" values and reset to default
  dev.new()
  par(defSet[1:65])
  parNames <- defNames[1:65]
  for (i in parNames)
    tryCatch(par(newDef[i]), error=function(er) 
      newDef[i] <<- defSet[i])
  dev.off()
  
  ##check for invalid digestR-specific parameters
  if (!newDef$type %in% c('auto', 'image', 'contour', 'filled', 'l', 'p', 'b'))
    newDef$type <- defSet$type
  if (newDef$position.1D < 0 || newDef$position.1D > 99)
    newDef$position.1D <- defSet$position.1D
  if (newDef$offset < -100 || newDef$offset > 100)
    newDef$offset <- defSet$offset
  if (!newDef$proj.type %in% c('l', 'p', 'b'))
    newDef$proj.type <- defSet$proj.type
  if (newDef$proj.direct != 1 && newDef$proj.direct != 2)
    newDef$proj.direct <- defSet$proj.direct	
  if (!newDef$peak.noiseFilt %in% c(0, 1, 2))
    newDef$peak.noiseFilt <- defSet$peak.noiseFilt
  if (!newDef$peak.labelPos %in% c('top', 'bottom', 'left', 'right', 
                                   'center'))
    newDef$peak.labelPos <- defSet$roi.labelPos
  if (newDef$clevel <= 0)
    newDef$clevel <- defSet$clevel
  if (newDef$nlevels < 1 || newDef$nlevels > 1000)
    newDef$nlevels <- defSet$nlevels
  if (any(newDef$roi.lwd <= 0))
    newDef$roi.lwd <- defSet$roi.lwd
  if (!newDef$roi.labelPos %in% c('top', 'bottom', 'left', 'right', 'center'))
    newDef$roi.labelPos <- defSet$roi.labelPos
  if (!newDef$roi.noiseFilt %in% c(0, 1, 2))
    newDef$roi.noiseFilt <- defSet$peak.noiseFilt
  
  ##check for valid library locations
  newDef$libLocs <- newDef$libLocs[file.exists(newDef$libLocs)]
  if (!length(newDef$libLocs))
    newDef$libLocs <- defSet$libLocs
  newDef$searchLibs <- newDef$searchLibs[file.exists(newDef$searchLibs)]
  if (!length(newDef$searchLibs))
    newDef$searchLibs <- defSet$searchLibs
  
  return(newDef)
}

## Internal function for writing defaultSettings out to file
## configFile - character string; file path to save the default settings to
## defSet - list; the object containing a list of default settings for digestR
writeDef <- function(configFile, defSet){
  if (missing(configFile))
    configFile <- file.path(path.expand('~'), '.digestR')
  if (missing(defSet))
    defSet <- defaultSettings
  dput(defSet, file=configFile)
  
  return(configFile)
}

## Internal function for reading in defaultSettings from file
## This function is depracted and is used only for backward-compatibility with
##	the previous defaultSettings format
## configFile - character string; full path to the file where defaultSettings
## 	is saved
## returns values from configFile as a list, in the standard defaultSettings
##	format
oldReadDef <- function(configFile, check=TRUE){
  
  ##read from file
  if (missing(configFile))
    configFile <- file.path(path.expand('~'), '.digestR')
  if (!file.exists(configFile))
    return(NULL)
  defText <- readLines(configFile, warn=FALSE)
  
  ##create list object from file text
  defNames <- defText[grep('#$', defText, fixed=TRUE)]
  defNames <- sapply(strsplit(defNames, '#$', fixed=TRUE), function(x) x[2])
  newDef <- as.list(rep("", length(defNames)))
  names(newDef) <- defNames
  valName <- NULL
  defVal <- NULL
  for (i in seq_along(defText)){
    
    ##increment loop if text is a blank line
    if (!nzchar(defText[i]))
      next
    
    ##increment loop and store value if text is a parameter name
    if (length(grep('#$', defText[i], fixed=TRUE))){
      valName <- unlist(strsplit(defText[i], '#$', fixed=TRUE))[2]
      next
    }
    
    ##get parameter value(s)
    defVal <- c(defVal, defText[i])
    nextVal <- defText[i + 1]
    if (i != length(defText) && nzchar(nextVal))
      next
    
    ##insert values into newDef
    if (!is.null(valName)){
      if (valName == 'filter' && defVal[1] != '0'){
        newDef[[valName]] <- function(x){}
        body(newDef[[valName]]) <- parse(text=paste(defVal, collapse='\n'))
      }else
        newDef[[valName]] <- defVal
    }
    defVal <- NULL
  }
  if (check)
    newDef <- checkDef(newDef)
  
  return(newDef)
}

## Internal function for reading in defaultSettings from file
## configFile - character string; full path to the file where defaultSettings
## 	is saved
## returns values from configFile as a list, in the standard defaultSettings
##	format
readDef <- function(configFile, check=TRUE){
  
  ##check configFile path and format
  if (missing(configFile))
    configFile <- file.path(path.expand('~'), '.digestR')
  if (!file.exists(configFile)){
    newDef <- NULL
  }else if (!length(readLines(configFile, warn=FALSE))){
    newDef <- NULL
  }else{	
    newDef <- tryCatch(dget(configFile), error=function(er)
      oldReadDef(configFile, check))
  }
  
  ##check values
  if (check)
    newDef <- checkDef(newDef)
  
  return(newDef)
}

## Internal utility function for creating digestR objects
## objList - character vector, a list of objects to create
## overwrite - logical, overwrites old data if TRUE
## returnObj - logical, returns a created object if TRUE
createObj <- function(objList, overwrite=FALSE, returnObj=FALSE){
  
  ## Turn off locator bell
  options(locatorBell=FALSE)
  
  ## Make a list of all digestR objects if objList is not provided
  if (missing(objList))
    objList <- c('currentSpectrum', 'defaultSettings', 'fileFolder',
                 'globalSettings', 'oldFolder', 'overlayList', 'pkgVar', 
                 'roiSummary', 'roiTable', 'assignments')
  forceCurrent <- forceDefault <- forceFile <- forceGlobal <- forceOld <-
    forceOverlay <- forcePkg <- forceSum <- forceRoi <- forceAssign <- FALSE
  
  ## Create defaultSettings
  defSet <- list(adj=0.5, ann=TRUE, ask=FALSE, bg="white", bty="7", cex=1, 
                 cex.axis=.95, cex.lab=1, cex.main=1, cex.sub=1, col="black",
                 col.axis="black", col.lab="black", col.main="black", col.sub="black", 
                 crt=0, err=0, family="", fg="black", fig=c(0, 1, 0, 1),
                 fin=c(10, 6.98958333333333), font=1, font.axis=1, font.lab=1, font.main=2, 
                 font.sub=1, lab=c(5, 5, 7), las=0, lend="round", lheight=1, ljoin="round", 
                 lmitre=10, lty="solid", lwd=1, mai=c(1.02, 0.82, 0.82, 0.42), 
                 ##TSB##			mar=c(2.25, 2.25, 1.5, 1), mex=1, mfcol=c(1, 1), mfg=c(1, 1, 1, 1), 
                 mar=c(3.20, 3.20, 1.5, 1), mex=1, mfcol=c(1, 1), mfg=c(1, 1, 1, 1),
                 mfrow=c(1, 1),	mgp=c(3, 1, 0), mkh=0.001, new=FALSE, oma=c(0, 0, 0, 0), 
                 omd=c(0, 1, 0, 1), omi=c(0, 0, 0, 0), pch=1, 
                 pin=c(8.76, 5.14958333333333), 
                 plt=c(0.082, 0.958, 0.145931445603577, 0.882682563338301), ps=12, pty="m", 
                 smo=1, srt=0, tck=NA, tcl=-0.5, usr=c(0, 10, 0,	10), xaxp=c(2, 8, 3), 
                 xaxs="r", xaxt="s", xlog=FALSE, xpd=FALSE, yaxp=c(20, 80, 6), yaxs="r", 
                 yaxt="s", ylog=FALSE,
                 pos.color="black", neg.color="green", conDisp=c(TRUE, TRUE), nlevels=20,
                 clevel=6, type="auto", theta=10,	phi=10,	asp=4, position.1D=.1, offset=0, 
                 proj.color='blue', proj.type="l",	proj.mode=FALSE, proj.direct=1, 
                 filter=function(x){range(x)[which.max(abs(range(x)))]}, peak.disp=FALSE, 
                 peak.color='black', peak.cex=.7, peak.pch='x', peak.labelPos='top', 
                 peak.noiseFilt=0, thresh.1D=6, roiMain=TRUE, roi.multi=TRUE, roiMax=TRUE, 
                 roi.bcolor=c('red', 'grey'), roi.tcolor=c('red', 'grey'),	
                 roi.lwd=c(2, 1), roi.lty=c('solid', 'dashed'), roi.labelPos='center', 
                 roi.cex=c(.95, .95), roi.noiseFilt=2, roi.w1=0, roi.w2=0, roi.pad=15, 
                 xtck=NA, ytck=NA, size.sub=c(6, 1.5), size.main=c(12.0, 6), #size.main=c(6, 4.25),	
                 mar.sub=c(.1, .1, 1.6, .1),size.multi=c(8, 6.25), #size.multi=c(6, 4.25), 
                 mar.multi=c(0, 0, 0, 0), cex.roi.multi=1, cex.files.multi=1, 
                 cex.roi.sub=1, overlay.text=TRUE, overlay.textSuppress = FALSE, autoBackup=TRUE, sdi=TRUE, update=TRUE, 
                 wd=path.expand('~'),
                 libLocs=gsub('\\', '/', system.file('Libraries/1H_13C_HSQC_pH7.4', 
                                                     package='digestR'), fixed=TRUE),
                 searchLibs=gsub('\\', '/', system.file('Libraries/1H_13C_HSQC_pH7.4', 
                                                        package='digestR'), fixed=TRUE), libUpdate=TRUE,
                 plotAA = FALSE, scaleSNR = FALSE, processSpeciesID = '', processSingleFile = TRUE,
                 speciesList = c('Homo sapiens',
                                 'Homo sapiens sickle', 
                                 'Plasmodium falciparum (3D7)', 
                                 'Pseudomonas aeruginosa 01', 
                                 'Pseudomonas aeruginosa 14', 
                                 'BosTaurus'
                                 ), 
                 speciesFiles = c('chrHs.csv', 
                                  'sickle.csv', 
                                  'chr3D7.csv', 
                                  'pa01.csv', 
                                  'pa14.csv',
                                  'chrBtaurus1.csv'
                                  ),
                 sd_noise_multiplier = 6, geneDisp = '', 
                 vectorType = data.frame(EndPoints = 'EndPoints', PepPoints = 'PepPoints', Mean = 'Mean', stringsAsFactors = FALSE),
                 plotStyle = data.frame(EndPoints = FALSE, Peptides = FALSE, Variance = FALSE, Accumulation = TRUE, stringsAsFactors = FALSE))  ##TSB##
  ## speciesFiles holds the file name that corresponds to the sequence information for the corresponding species in speciesList
  ## **for now** -- all point to the chrHs.csv which holds homo sapiens data.'
  
  ## Check defaultSettings
  if ('defaultSettings' %in% objList)
  {
    
    ## Replace defaultSettings if it is not the correct format
    if (!exists('defaultSettings', envir=.GlobalEnv) || 
        !is.list(defaultSettings) || length(defaultSettings) < length(defSet))
      forceDefault <- TRUE
    
    if (!returnObj)
    {
      
      ## Assign defaultSettings 
      #			if (overwrite || forceDefault)
      #			{
      assign('defaultSettings', defSet, envir=.GlobalEnv)
      #			}
      
      ## Read defaultSettings from file	
      #			configFile <- file.path(path.expand('~'), '.digestR')
      #			if (length(configFile) && file.exists(configFile))
      #			{
      #				defSet <- readDef(configFile)
      #				assign('defaultSettings', defSet, envir=.GlobalEnv)
      #			}
    }
  }
  
  ## Create globalSettings
  globalPars <- c('offset', 'position.1D', 'filter', 'proj.direct', 'proj.mode', 
                  'proj.type', 'peak.disp', 'peak.noiseFilt', 'thresh.1D', 'peak.pch', 
                  'peak.cex', 'peak.labelPos', 'roiMain', 'roiMax', 'roi.bcolor', 
                  'roi.tcolor', 'roi.lwd', 'roi.lty', 'roi.cex', 'roi.labelPos', 
                  'roi.noiseFilt', 'roi.w1', 'roi.w2', 'roi.pad', 'cex.roi.multi', 
                  'cex.files.multi', 'cex.roi.sub', 'size.main', 'size.sub', 'size.multi', 
                  'mar', 'mar.sub', 'mar.multi', 'overlay.text', 'overlay.textSuppress', 'processSpeciesID', 
                  'processSingleFile', 'speciesList', 'speciesFiles', 'sd_noise_multiplier', 'geneDisp', 'vectorType', 'plotStyle')
  
  defGlobal <- defSet[globalPars]
  
  ## Create fileFolder
  defFile <- NULL
  
  ## Create oldFolder
  defOld <- list( undo.index = 0, assign.index = 0, fileFolder = NULL, 
                  roiTable = NULL, currentSpectrum = NULL,	roiSummary = NULL, 
                  overlayList = NULL, zoom.history = NULL, zoom.list = NULL)
  
  ## Create currentSpectrum
  if (exists('fileFolder'))
    defCurrent <- names(fileFolder)[length(fileFolder)]
  else	
    defCurrent <- NULL
  
  ## Create overlaylist
  defOverlay <- NULL
  
  ## Create roiTable
  defRoi <- NULL
  
  ## Create roiSummary
  defSummary <- NULL
  
  ## Create pkgVar
  defPkg <- list()
  defPkg$prevDir <- defSet$wd
  
  # Remove this to see if the error disappears
  #defPkg$version <- suppressWarnings(
  #  paste(packageDescription('digestR', fields='Version'), ' (',	packageDescription('digestR', fields='Date'), ')', sep='')
  #)
  
  ## Create assignments
  defAssign <- NULL
  
  ## Returns the newly created object
  if (returnObj)
    return(switch(objList[1], 'currentSpectrum'=defCurrent, 
                  'defaultSettings'=defSet, 
                  'fileFolder'=defFile,
                  'globalSettings'=defGlobal, 
                  'oldFolder'=defOld, 
                  'overlayList'=defOverlay, 
                  'pkgVar'=defPkg, 
                  'roiSummary'=defSummary, 
                  'roiTable'=defRoi, 
                  'assignments'=defAssign))	
  
  ## Check globalSettings
  if ('globalSettings' %in% objList){
    
    ## Replace globalSettings if it is not the correct format
    if (!exists('globalSettings', envir=.GlobalEnv) || !is.list(globalSettings) 
        || length(globalSettings) != length(defGlobal))
      forceGlobal <- TRUE
    
    ## Assign globalSettings
    if (overwrite || forceGlobal)
      assign('globalSettings', defGlobal, envir=.GlobalEnv)
  }
  
  ## Check fileFolder
  if ('fileFolder' %in% objList){
    
    ## Replace fileFolder if it is not the correct format
    if (!exists('fileFolder', envir=.GlobalEnv) || !is.null(fileFolder) && 
        !is.list(fileFolder))
      forceFile <- TRUE
    
    ## Assign fileFolder
    if (overwrite || forceFile)
      assign('fileFolder', defFile, envir=.GlobalEnv)
  }
  
  ## Check oldFolder
  if ('oldFolder' %in% objList){
    
    ## Replace oldFolder if it is not the correct format
    if (!exists('oldFolder', envir=.GlobalEnv) || !is.list(oldFolder))
      forceOld <- TRUE
    
    ## Assign oldFolder
    if (overwrite || forceOld)
      assign('oldFolder', defOld, envir=.GlobalEnv)
  }
  
  ## Check currentSpectrum
  if ('currentSpectrum' %in% objList){
    
    ## Replace currentSpectrum if it is not the correct format
    if (!exists('currentSpectrum', envir=.GlobalEnv) || 
        !is.null(currentSpectrum) && (!nzchar(currentSpectrum) || 
                                      length(currentSpectrum) != 1 || !is.character(currentSpectrum)))
      forceCurrent <- TRUE
    
    ## Assign currentSpectrum
    if (overwrite || forceCurrent)
      assign('currentSpectrum', defCurrent, envir=.GlobalEnv)
  }
  
  ## Check overlaylist
  if ('overlayList' %in% objList){
    
    ## Replace overlayList if it is not the correct format
    if (!exists('overlayList', envir=.GlobalEnv) || !is.null(overlayList) && 
        !is.character(overlayList))
      forceOverlay <- TRUE
    
    ## Assign overlayList
    if (overwrite || forceOverlay)
      assign('overlayList', defOverlay, envir=.GlobalEnv)
  }
  
  ## Check roiTable
  if ('roiTable' %in% objList){
    
    ## Replace roiTable if it is not the correct format
    if (!exists('roiTable', envir=.GlobalEnv) || !is.null(roiTable) && 
        !is.data.frame(roiTable))
      forceRoi <- TRUE
    
    ## Assign roiTable
    if (overwrite || forceRoi)
      assign('roiTable', defRoi, envir=.GlobalEnv)
  }
  
  ## Check roiSummary
  if ('roiSummary' %in% objList){
    
    ## Replace roiSummary if it is not the correct format
    if (!exists('roiSummary', envir=.GlobalEnv) || !is.null(roiSummary) && 
        !is.list(roiSummary))
      forceSum <- TRUE
    
    ## Assign roiSummary
    if (overwrite || forceSum)
      assign('roiSummary', defSummary, envir=.GlobalEnv)
  }
  
  ## Check pkgVar
  if ('pkgVar' %in% objList){
    
    ## Replace pkgVar if it is not the correct format
    if (exists('pkgVar', envir=.GlobalEnv) && !is.null(pkgVar$prevDir) && 
        !overwrite && pkgVar$prevDir != path.expand('~'))
      defPkg$prevDir <- pkgVar$prevDir
    
    ## Assign pkgVar
    assign('pkgVar', defPkg, envir=.GlobalEnv)
  }
  
  ## Check assignments
  if ('assignments' %in% objList){
    
    ## Replace assignments if it is not the correct format
    if (!exists('assignments', envir=.GlobalEnv) || 
        !is.null(assignments) && !is.data.frame(assignments))
      forceAssign <- TRUE
    
    ## Assign roiTable
    if (overwrite || forceAssign)
      assign('assignments', defAssign, envir=.GlobalEnv)
  }
  
  invisible()
}

## Patch for code compatibilty with older digestR workspaces
## delete - logical, deletes old objects if TRUE
patch <- function(delete=TRUE){
  
  ## Do not apply the patch if digestR version is up to date
  if (exists('pkgVar') && identical(pkgVar, createObj('pkgVar', 
                                                      returnObj=TRUE)$version))
    return(invisible())
  
  ## Check for missing values in defaultSettings
  defSet <- createObj('defaultSettings', returnObj=TRUE)
  defaultSettings <- readDef(check=FALSE)
  if (!is.null(defaultSettings)){
    if (length(defaultSettings) == length(defSet)){
      defSet <- checkDef(defaultSettings)
    }else{
      
      ## Append missing values
      missingNames <- names(defSet)[!names(defSet) %in% names(defaultSettings)]
      defaultSettings[missingNames] <- defSet[missingNames]
      
      ## Write out defaultSettings
      defaultSettings <- checkDef(defaultSettings)
      writeDef(defSet=defaultSettings)
      defSet <- defaultSettings
    }
    myAssign("defaultSettings", defaultSettings, save.backup=FALSE)
  }
  
  ## Check fileFolder for correct structure
  if (exists('fileFolder') && !is.null(fileFolder) && length(fileFolder)){
    for (i in seq_along(fileFolder)){
      if (length(fileFolder[[i]]$graphics.par) < length(defSet)){
        missingItems <- defSet[!names(defSet) %in% 
                                 names(fileFolder[[i]]$graphics.par)]
        fileFolder[[i]]$graphics.par <- c(fileFolder[[i]]$graphics.par, 
                                          missingItems)
        fileFolder[[i]]$graphics.par$mar <- defSet$mar
        fileFolder[[i]]$graphics.par$cex.main <- defSet$cex.main
        fileFolder[[i]]$graphics.par$cex.axis <- defSet$cex.axis
      }
      if (is.null(fileFolder[[i]]$file.par$file.size) || 
          is.null(fileFolder[[i]]$file.par$date.modified)){
        fileInfo <- tryCatch(file.info(fileFolder[[i]]$file.par$file.name), 
                             error=function(er) NULL)
        fileFolder[[i]]$file.par$file.size <- fileInfo$size
        fileFolder[[i]]$file.par$date.modified <- fileInfo$mtime
      }
      if (is.null(fileFolder[[i]]$file.par$user_title)){
        userTitle <- basename(names(fileFolder)[i])
        userTitles <- sapply(fileFolder, function(x) x$file.par$user_title)
        if (userTitle %in% userTitles)
          userTitle <- names(fileFolder)[i]
        fileFolder[[i]]$file.par$user_title <- userTitle
      }
    }
    myAssign("fileFolder", fileFolder, save.backup=FALSE)
  }
  if (exists('file.folder') && !exists('fileFolder')){
    createObj(c('defaultSettings', 'globalSettings'), overwrite=TRUE)		
    
    ## Update fileFolder
    fileFolder <- file.folder
    if( is.list(fileFolder) ){
      for( i in 1:length(fileFolder) ){
        usr <- fileFolder[[i]]$graphics.par$usr
        fileFolder[[i]]$graphics.par <- defaultSettings
        fileFolder[[i]]$graphics.par$usr <- 
          c(rev(sort(usr[1:2])), rev(sort(usr[3:4])))
      }
    } 
    myAssign("fileFolder", fileFolder, save.backup = FALSE )
  }
  
  ## Rename digestR objects
  if(exists('roi.table') && (!exists('roiTable') || (exists('roiTable') && 
                                                     is.null(roiTable))))
    myAssign("roiTable", roi.table, save.backup=FALSE)
  if(exists('overlay.list') && (!exists('overlayList') || (exists('overlayList') 
                                                           && is.null(overlayList))))
    myAssign("overlayList", overlay.list, save.backup=FALSE)
  if(exists('old.folder') && (!exists('oldFolder') || (exists('oldFolder') &&
                                                       is.null(oldFolder$fileFolder))))
    myAssign("oldFolder", createObj('oldFolder', returnObj=TRUE), 
             save.backup=FALSE)
  if(exists('current.roi.summary') && (!exists('roiSummary') || 
                                       (exists('roiSummary') && is.null(roiSummary))))
    myAssign("roiSummary", current.roi.summary, save.backup=FALSE)
  if(exists('prevDir') && (!exists('pkgVar') || (exists('pkgVar') && 
                                                 is.null(pkgVar$prevDir)))){
    pkgVar <- createObj('pkgVar', returnObj=TRUE)
    pkgVar$prevDir <- prevDir
    myAssign('pkgVar', pkgVar)
  }
  
  ## Create any missing digestR objects
  createObj()
  
  ## Update roiTable format
  if( !is.null(roiTable) ){
    if(length(which(names(roiTable) == 'nDim')) == 0){
      roiTable$nDim <- rep(2, nrow(roiTable))
      myAssign("roiTable", roiTable, save.backup = FALSE) 
    }
  }
  
  ## Update roiSummary format
  if( !is.null(roiSummary) ){
    if(!is.list(roiSummary))
      myAssign("roiSummary", NULL, save.backup = FALSE)
    else if (!is.null(roiSummary$data) && 'GROUP' %in% names(roiSummary$data)){
      newSum <- NULL
      newSum$data <- roiSummary$data[2:ncol(roiSummary$data)]
      newSum$summary.par$summary.type <- 'maximum'
      newSum$summary.par$norm.data.source <- NA
      if (!is.null(roiSummary$summary.par$normalization.ROIs)){
        if (is.na(roiSummary$summary.par$normalization.ROIs[1]))
          newSum$summary.par$normalization <- 'none'
        else if (roiSummary$summary.par$normalization.ROIs[1] == 
                 'Signal to noise')
          newSum$summary.par$normalization <- 'signal/noise'
        else{
          newSum$summary.par$normalization <- 'internal'
          newSum$summary.par$norm.data.source <- 
            roiSummary$summary.par$normalization.ROIs
        }
      }else
        newSum$summary.par$normalization <- 'none'
      myAssign("roiSummary", newSum, save.backup=FALSE)
    }
  } 
  
  ## Update version
  if (exists('pkgVar')){
    pkgVar$version <- createObj('pkgVar', returnObj=TRUE)$version
    myAssign('pkgVar', pkgVar, save.backup=FALSE)	
  }
  
  ## Remove all of the old objects
  oldObj <- c('about', 'addGui', 'ANOVA', 'appendPeak', 'load', 'pReg',
              'assignGroups', 'autoRef', 'bringFocus', 'changeColor', 'changeRoi', 'co', 
              'buttonDlg', 'createObj', 'ct', 'ct1D', 'ct2D', 'ctd', 'ctu', 'cw', 'da', 
              'dd', 'devGui', 'di', 'dp', 'dr', 'draw2D', 'drawPeptides', 'drf', 'ed', 'err', 
              'export', 'fancyPeak2D', 'fc', 'ff', 'fillCon', 'findTiles', 'fo', 'gui', 
              'hideGui', 'import', 'isNoise', 'loc', 'localMax', 'matchShift', 
              'maxShift', 'mmcd', 'more', 'myAssign', 'myDialog', 'myMsg', 'myOpen', 
              'mySave', 'mySelect', 'newRange', 'nf', 'obs2List', 'ol', 'orderROI', 
              'overlays', 'pa', 'paAll', 'pan', 'patch', 'pd', 'pDel', 'pDelAll', 
              'pdisp', 'pe', 'peakDel', 'peakPick', 'peakPick1D', 'peakPick2D', 
              'peakVolume', 'per', 'persp2D', 'ph', 'pj', 'pjv', 'pl', 'plot1D', 
              'plot2D', 'popupGui', 'pp', 'pr', 'proj1D', 'pseudo1D', 'pu', 'pv', 'pw', 
              'pwAll', 'pz', 'ra', 'randomColors', 'rc', 'rcd', 'rci', 'rd', 'rdAll', 
              'rDel', 're', 'recall', 'red', 'refresh', 'regionMax', 'rei', 'reset', 
              'rmd', 'rml', 'rmr', 'rmu', 'rn', 'roi', 'roi.anova', 'roi.pca', 
              'roi.var', 'roi.xy', 'roi.ztdist', 'rotc', 'rotcc', 'rotd', 'rotu', 'rp', 
              'rpAll', 'rs', 'rsAll', 'rsf', 'rSum', 'rv', 'rvm', 'rvs', 'selList', 
              'selMain', 'selMulti', 'selSub', 'setGraphics', 'setWindow', 'shiftToROI',
              'showGui', 'showRoi', 'soon', 'spin', 'sr', 'sr1D', 'sr2D', 'ss', 
              'tclCheck', 'trans2Peak', 'ucsf1D', 'ucsf2D', 'ucsfHead', 'ucsfTile', 
              'ud', 'vp', 'vpd', 'vpu', 'vs', 'wc', 'wl', 'ws', 'xy.plot', 'zc', 'zf', 
              'zi', 'zm', 'zo', 'zp', 'zz', 'anova.plot', 'assign.groups', 'auto.ref', 
              'auto.roi', 'change.roi', 'clickZoom', 'crd', 'cri', 'ct.all', 
              'current.draw.all', 'current.graphics.all', 'default.graphics', 
              'default.settings', 'delete.roi', 'deselect.all', 'draw.roi', 'draw2D', 
              'draw2d', 'drf', 'drff', 'e', 'edit.roi', 'erd', 'eri', 'export.roi', 
              'get.file.name', 'hp', 'import.roi', 'import.summary', 'load.groups', 
              'match.shifts', 'max.shift', 'modify.plot', 'mrd', 'mrl', 'mrr', 'mru', 
              'my.assign', 'my.biplot', 'my.filled', 'neg.only', 'p.clear', 
              'p.clear.all', 'p.edit', 'p.off', 'p.on', 'p.print', 'p.roi', 'p.roi.all',
              'p.save', 'pa.all', 'pc', 'PCA', 'pca.plot', 'peak.pick', 'peak.pick.1D', 
              'peak.pick.2D', 'peak.volume', 'plot.colors', 'plot.gui', 'plot.label', 
              'plot.popup.gui', 'pos.neg', 'pos.only', 'print.data', 'print.graphics', 
              'pw.all', 'pz', 'random.colors', 'renumber.rois', 'replot', 'roi.files', 
              'roi.plot', 'roi.summary', 'save.data', 'select.all', 'select.roi', 
              'set.graphics', 'set.stat.graphics', 'show.roi', 'stat.default.graphics', 
              'stat.print.graphics', 'v.off', 'v.on', 'v1d', 'vabs', 
              'variance.stabilize', 'vd', 'vi', 'view.plot', 'vmax', 'vmin', 'vsi', 
              'z.plot', 'z.test', 'roi.summary', 'roi.table', 'file.folder', 'prevDir',
              'old.folder', 'overlay.list', 'current.roi.summary', 'pdisp.all', 
              'pdisp.off', 'pdisp.off.all', 'plot.roi.summary', 'pm', 'pm.abs', 'pmv', 
              'pmv.abs', 'save.roi.summary', 'ucsfData', 'vsv', 'sdi', 'mdi', 'cf', 
              'tkGuis', 'convList', 'rd', 'up', 'cs', 'csd', 'pd', 'gp' )
  
  if (delete)
    suppressWarnings(rm(list=oldObj, envir=.GlobalEnv))
  updateFiles(halt=FALSE)
}


################################################################################
##                                                                            ##
##   Internal functions for showing, hiding and working with digestR GUIs     ##
##                                                                            ##
################################################################################

## Internal function for checking that the tcltk package is loaded
tclCheck <- function(){
  
  ##load tcltk package
  tryCatch(suppressMessages(library(tcltk)), error=function(er)
    stop("digestR requires Tcl/Tk version 8.5 or greater", call.=FALSE))
  
  ##make sure tcl 8.5 or later is installed
  tclVers <- as.numeric(tcl("info", "tclversion"))
  if (tclVers < 8.5)
    stop("digestR requires Tcl/Tk version 8.5 or greater", call.=FALSE)
}

## Internal function for creating tcl images from files included with digestR
## imageName - character string; the name for the tcl image to create
## path - character sting; path to the image file (GIF only) to use, must be 
##	relative to the digestR package directory
createTclImage <- function(imageName, path=NULL){
  tclImages <- as.character(tcl('image', 'names'))
  if (imageName %in% tclImages)
    return(invisible())
  if (is.null(path))
    path <- paste(imageName, 'gif', sep='.')
  imagePath <- system.file(path, package='digestR')
  if (!file.exists(imagePath))
    err(paste('Image file:', imagePath, 'does not exist.'))
  tcl('image', 'create', 'photo', imageName, '-file', imagePath)
}

## Internal replacement function for tktoplevel
## Creates a toplevel widget with a given name
## id - character; the pathName for the toplevel, must begin with '.'
## parent - character; parent for the toplevel
## Note:  If a toplevel with the given id already exists the toplevel will be 
##   deiconified (redisplayed)
# myToplevel <- function (id, parent, ...){
  
#   ##check arguments
#   tclCheck()
#   if (missing(parent))
#     parent <- .TkRoot
#   if (missing(id)){
#     tryCatch(id <- paste(parent$ID, evalq(num.subwin <- num.subwin + 1, 
#                                           parent$env), sep = "."),
#              error=function(er){
#                assign('num.subwin', parent$env$num.subwin + 1, parent$env)
#                id <<- paste(parent$ID, parent$env$num.subwin, sep = ".")
#              })
#   }else if (length(unlist(strsplit(id, '.', fixed=TRUE))) == 1)
#     id <- paste(parent$ID, '.', id, sep='')
  
#   ##if a toplevel with the same id already exists, display it
#   if (as.logical(tcl('winfo', 'exists', id))){
#     hideGui(id)
#     showGui(id)
#     tkfocus(id)
#     return(NULL)
#   }
  
#   ##create a new window environment
#   win <- .Tk.newwin(id)
#   assign(id, win, envir=parent$env)
#   assign("parent", parent, envir=win$env)
  
#   ##create the new toplevel
#   win$ID <- id
#   tcl("toplevel", id, ...)
  
#   ##configure the window to be displayed on top of its parent
#   if (parent$ID != ""){
#     parentTop <- as.logical(tcl('wm', 'attributes', parent, '-topmost'))
#     if (parentTop)
#       tcl('wm', 'attributes', parent, topmost=FALSE)
#     if (as.logical(tkwinfo('viewable', parent)))
#       tkwm.transient(win, parent)
#     tkwm.withdraw(win)
#     tkwm.deiconify(win)
#     tkbind(id, "<Destroy>", function(){
#       if (parentTop)
#         tryCatch(tcl('wm', 'attributes', parent, topmost=TRUE), 
#                  error=function(er){})
#       if (exists(id, envir=parent$env, inherits=FALSE)) 
#         rm(list=id, envir=parent$env)
#       tkbind(id, "<Destroy>", "")})
#   }else{
#     tkbind(id, "<Destroy>", function(){
#       if (exists(id, envir=parent$env, inherits=FALSE)) 
#         rm(list=id, envir=parent$env)
#       tkbind(id, "<Destroy>", "")})
#   }
  
#   return(win)
# }
####################################################################
myToplevel <- function (id, parent, ...) {
  
  ## Check arguments
  tclCheck()
  if (missing(parent)) {
    parent <- .TkRoot
  }
  
  if (missing(id)) {
    tryCatch(id <- paste(parent$ID, evalq(num.subwin <- num.subwin + 1, 
                                          parent$env), sep = "."),
             error = function(er) {
               assign('num.subwin', parent$env$num.subwin + 1, parent$env)
               id <<- paste(parent$ID, parent$env$num.subwin, sep = ".")
             })
  } else if (length(unlist(strsplit(id, '.', fixed = TRUE))) == 1) {
    id <- paste(parent$ID, '.', id, sep = '')
  }
  
  ## If a toplevel with the same id already exists, reuse it
  if (as.logical(tcl('winfo', 'exists', id))) {
    hideGui(id)
    showGui(id)
    
    ## Ensure the window is raised and focused
    tkwm.deiconify(id)   # Ensure the window is visible
    tkraise(id)          # Bring the window to the front
    tkfocus(id)          # Focus on the window
    
    ## Force an update to ensure the window is drawn immediately
    tcl("update")
    
    return(get(id, envir = parent$env))
  }
  
  ## Create a new window environment
  win <- .Tk.newwin(id)
  assign(id, win, envir = parent$env)
  assign("parent", parent, envir = win$env)
  
  ## Create the new toplevel
  win$ID <- id
  tcl("toplevel", id, ...)
  
  ## Force the window to be visible and bring it to the front
  tkwm.deiconify(id)
  tkraise(id)
  tkfocus(id)
  
  ## Force the update to ensure the window is drawn
  tcl("update")
  
  ## Configure the window to be displayed on top of its parent
  if (parent$ID != "") {
    parentTop <- as.logical(tcl('wm', 'attributes', parent, '-topmost'))
    if (parentTop)
      tcl('wm', 'attributes', parent, topmost = FALSE)
    if (as.logical(tkwinfo('viewable', parent)))
      tkwm.transient(win, parent)
    
    ## Withdraw and deiconify the window to force redrawing
    tkwm.withdraw(win)
    tkwm.deiconify(win)
    
    ## Set up destroy callback
    tkbind(id, "<Destroy>", function() {
      if (parentTop)
        tryCatch(tcl('wm', 'attributes', parent, topmost = TRUE), 
                 error = function(er) {})
      if (exists(id, envir = parent$env, inherits = FALSE)) 
        rm(list = id, envir = parent$env)
      tkbind(id, "<Destroy>", "")
    })
  } else {
    tkbind(id, "<Destroy>", function() {
      if (exists(id, envir = parent$env, inherits = FALSE)) 
        rm(list = id, envir = parent$env)
      tkbind(id, "<Destroy>", "")
    })
  }
  
  ## Return the window object
  return(win)
}

## Internal function for destroying a tk GUI that is currently open 
## guiName - character string: the name for the GUI
closeGui <- function(guiName){
  if (missing(guiName)){
    for (i in ls(envir=.TkRoot$env, all.names=TRUE))
      tryCatch(tkdestroy(i), error=function(er){})
    guiList <- c('per', 'ps', 'co', 'ct', 'os', 'sr', 'pj', 'roi', 'zm', 'pp', 
                 'fs', 'ep', 'ca', 'cf', 'aa', 'pm', 'gl', 'vd', 'up', 'csd', 'pd', 'cs','gp') ### TSB -- added pm and gl
    for (i in guiList)
      tryCatch(tkdestroy(paste('.', i, sep='')), error=function(er){})
  }else{
    for (i in guiName){
      tryCatch(tkdestroy(i), error=function(er){})
      tryCatch(tkdestroy(paste('.', i, sep='')), error=function(er){})
    }
  }
  
  bringFocus()
}

## Internal function for hiding open tk GUIs
## guiName - character string: the name for the GUI to hide
hideGui <- function(guiName){
  if (missing(guiName)){
    for (i in ls(envir=.TkRoot$env, all.names=TRUE))
      tryCatch(tkwm.iconify(i), error=function(er){})
    guiList <- c('per', 'ps', 'co', 'ct', 'os', 'sr', 'pj', 'roi', 'zm', 'pp', 
                 'fs', 'ep', 'ca', 'cf', 'aa', 'pm', 'gl', 'vd', 'up', 'csd', 'pd', 'cs', 'gp') ### TSB -- added pm and gl
    for (i in guiList)
      tryCatch(tkwm.iconify(paste('.', i, sep='')), error=function(er){})
  }else{
    for (i in guiName){
      tryCatch(tkwm.iconify(i), error=function(er){})
      tryCatch(tkwm.iconify(paste('.', i, sep='')), error=function(er){})
    }
  }
  bringFocus()
}

## Internal function for hiding open tk GUIs
## guiName - character string: the name for the GUI to show
showGui <- function(guiName){
  if (missing(guiName)){
    for (i in ls(envir=.TkRoot$env, all.names=TRUE))
      tryCatch(tkwm.deiconify(i), error=function(er){})
    guiList <- c('per', 'ps', 'co', 'ct', 'os', 'sr', 'pj', 'roi', 'zm', 'pp', 
                 'fs', 'ep', 'ca', 'cf', 'aa', 'pm', 'gl', 'vd', 'up', 'csd', 'pd', 'cs', 'gp') ### TSB -- added pm and gl
    for (i in guiList)
      tryCatch(tkwm.deiconify(paste('.', i, sep='')), error=function(er){})
  }else{
    for (i in guiName){
      tryCatch(tkwm.deiconify(i), error=function(er){})
      tryCatch(tkwm.deiconify(paste('.', i, sep='')), error=function(er){})
    }
  }
  #	bringFocus()
}

## Internal graphics function gui
## Sets window popup and dropdown windows 
## top - specifies a tk toplevel to send the tk menus to
#' Launches the digestR GUI Menu
#'
#' This function creates a graphical user interface (GUI) menu for interacting
#' with the digestR package. The GUI menu provides various options for manipulating
#' CSV and DCF files, performing editing tasks, managing graphics, and accessing help
#' topics related to the digestR package.
#'
#' @param top An optional parameter indicating the top-level window to use for the GUI.
#'            If not provided, a new top-level window will be created.
#'
#' @details On Windows systems with Rgui, the function adds menu items to the Rgui
#'          interface. On other systems or when a top-level window is provided, the
#'          function creates a GUI using the Tcl/Tk package to display menu options.
#'
#' @seealso \code{\link{myToplevel}}, \code{\link{tclCheck}}, \code{\link{tkmenu}},
#'          \code{\link{tkadd}}, \code{\link{tkconfigure}}, \code{\link{tkfocus}},
#'          \code{\link{tkwm.deiconify}}
#'
#'
#' @examples
#' # Launch the digestR GUI
#' gui()
#'
#' # Use a pre-existing top-level window for the GUI
#' top <- myToplevel('menu', width = 400, height = 300)
#' gui(top)
#'
#' @export
gui <- function(top=NULL){
  
  #splashScreen()

  if (.Platform$OS.type == 'windows' && .Platform$GUI == 'Rgui' && 
      is.null(top)){
    if ("  digestR -->  " %in% winMenuNames())
      return(invisible())
    winMenuAdd("  digestR -->  ")
    
    winMenuAdd("Manipulate .csv files")
    
    #winMenuAddItem("Manipulate csv", 'Unique peptides				unique_peptides()/up()',"up()")
    winMenuAddItem("Manipulate .csv files", 'Process Mascot Files			pm()', "pm()")
    winMenuAddItem("Manipulate .csv files", 'Generate New Proteome			gp()', "gp()")
    winMenuAddItem("Manipulate .csv files", 'Plot cut site distribution		        csd()', "csd()")
    winMenuAddItem("Manipulate .csv files", 'Plot peptide distribution			pd()', "pd()")
    winMenuAddItem("Manipulate .csv files", 'Venn Diagram			        vd()', "vd()")
    #winMenuAddItem("Manipulate csv", 'Import Maven Files			im()', "im()")
    
    winMenuAdd("Manipulate .dcf files")
    #winMenuAddItem("Manipulate dcf", 'Open/Close files		fs()', "fs()")
    winMenuAddItem('Manipulate .dcf files', 'Manipulate Files		mf()', "mf()")
    winMenuAddItem('Manipulate .dcf files', 'Plot settings	        ct()', "ct()")
    winMenuAddItem('Manipulate .dcf files', 'Overlays			ol()', "ol()")
    winMenuAddItem("Manipulate .dcf files", 'Save as			sa()', "sa()")   
    
    winMenuAdd("Edit")
    winMenuAddItem("Edit", 'Undo				ud()', "ud()")
    winMenuAddItem("Edit", 'Redraw spectrum 				dd()', "dd()")
    #winMenuAddItem("Edit", 'Redo                  redo/rd()', "redo/rd()") 
    #winMenuAddItem("Edit", 'Genotype            genotype/gt()', "genotype/gt()")
    #winMenuAddItem("Edit", 'Title                 title/ti()', "title/ti()")		
    #winMenuAddItem("Edit", '--', "none")
    #winMenuAddItem("Edit", 'Preferences                       ep()', "ep()")
    
    winMenuAdd("Graphics")
    winMenuAddItem("Graphics", 'Plot colors			co()', "co()")
    #winMenuAddItem("Graphics", 'Plot settings			ct()', "ct()")
    
    winMenuAdd("View")
    winMenuAddItem("View", 'Zoom					zm()', "zm()")
    #winMenuAddItem("View", 'Overlays					overlay/ol()', "ol()")
    winMenuAddItem("View", 'Display protease cut site			cs()', "cs()")
    winMenuAddItem("View", 'Gene labeling				glab()', "glab()")		
    #winMenuAddItem("View", 'Redraw spectrum				dd()', "dd()")
    
    winMenuAdd("Help")
    winMenuAddItem("Help", 'Help topics', "?digestR")
    winMenuAddItem("Help", 'List functions', "?more")
    winMenuAddItem("Help", 'User manual', "digestR('user_manual', TRUE)")
    winMenuAddItem("Help", 'Developer\'s guide', 
                   "digestR:::myHelp('developers_guide/developers_guide', TRUE)")
    winMenuAddItem("Help", '--', "none")
    winMenuAddItem("Help", 'Homepage', 
                   "browseURL(PROJECT_URL)")
    winMenuAddItem("Help", 'Update digestR', "digestR:::updater()")
    winMenuAddItem("Help", 'About digestR', "digestR:::about()")
    
  } else {
    tclCheck()
    if (is.null(top)){
      if (.Platform$OS.type == 'windows')
        top <- myToplevel('menu', width=255, height=30)
      else
        top <- myToplevel('menu', width=285, height=1)
      if (is.null(top))
        return(invisible())
      tkwm.title(top, 'digestR Menu')
      tcl('wm', 'attributes', top, topmost=TRUE)
    }
    topMenu <- tkmenu(top)
    tkconfigure(top, menu=topMenu)
    
    fileMenu <- tkmenu(topMenu, tearoff=FALSE)
    tkadd(fileMenu, 'command', label='Open/Clofo()se files', accelerator='fs()', 
          command=function() fs())
    tkadd(fileMenu, 'separator') 
    tkadd(topMenu, 'cascade', label='File', menu=fileMenu)
    
    editMenu <- tkmenu(topMenu, tearoff=FALSE)
    tkadd(editMenu, 'command', label='Undo', accelerator='ud()', 
          command=function() ud())
    tkadd(editMenu, 'command', label='Redo', accelerator='rd()', 
          command=function() rd())
    
    tkadd(editMenu, 'command', label='Edit Genotype', accelerator='gt()', 
          command=function() gt())
    
    tkadd(editMenu, 'command', label='Edit Title', accelerator='ti()', 
          command=function() ti())												
    tkadd(editMenu, 'separator')
    tkadd(editMenu, 'command', label='Preferences', accelerator='ep()', 
          command=function() ep())
    tkadd(topMenu, 'cascade', label='Edit', menu=editMenu)
    
    graphicsMenu <- tkmenu(topMenu, tearoff=FALSE)
    tkadd(graphicsMenu, 'command', label='Plot colors', 
          accelerator='co()', command=function() co())
    tkadd(graphicsMenu, 'command', label='Plot settings', 
          accelerator='ct()',	command=function() ct())	
    tkadd(graphicsMenu, 'command', label='Perspective', accelerator='per()', 
          command=function() per())
    tkadd(topMenu, 'cascade', label='Graphics', menu=graphicsMenu)
    
    viewMenu <- tkmenu(topMenu, tearoff=FALSE)
    tkadd(viewMenu, 'command', label='Zoom', accelerator='zm()', 
          command=function() zm())
    tkadd(viewMenu, 'command', label='Overlays', accelerator='ol()', 
          command=function() ol())
    
    tkadd(viewMenu, 'command', label='Gene labeling', accelerator='gl()', 
          command=function() gl())				
    
    tkadd(viewMenu, 'command', label='Redraw main spectrum', 
          accelerator='dd()', command=function() dd())

    helpMenu <- tkmenu(topMenu, tearoff=FALSE)
    tkadd(helpMenu, 'command', label='Help topics',	
          command=function(...) digestR:::myHelp('digestR-package'))
    tkadd(helpMenu, 'command', label='List functions', 
          command=function(...) digestR:::myHelp('more'))
    tkadd(helpMenu, 'command', label='Update digestR', 
          command=function() digestR:::updater())
    tkadd(helpMenu, 'command', label='User manual', 
          command=function(...) digestR:::myHelp('user_manual', TRUE))
    tkadd(helpMenu, 'command', label='Developer\'s guide', command=function(...) 
      digestR:::myHelp('developers_guide/developers_guide', TRUE))
    tkadd(helpMenu, 'command', label='Homepage', 
          command=function(...) browseURL(PROJECT_URL))
    tkadd(helpMenu, 'command', label='About digestR',
          command=function() digestR:::about())
    tkadd(topMenu, 'cascade', label='Help', menu=helpMenu)
    
    tkfocus(top)
    tkwm.deiconify(top)
    
    invisible()
  }
}

#########
## TSB function 
## Toggle state of plotAA logical flag
#########
pAA <- function(){	
  current <- wc()
  in.folder <- fileFolder
  
  in.folder[[current]]$graphics.par$plotAA <- !(in.folder[[current]]$graphics.par$plotAA) 
  
  myAssign('fileFolder', in.folder, save.backup = TRUE)
  refresh()
}

## Internal graphics function popupGui
## Sets window popup and dropdown windows 
popupGui <- function(dev){
  
  if (.Platform$GUI != 'Rgui')
    return(invisible())
  
  ##creates context menu for the main plot window
  if (dev == 'main' && !'$Graph2Popup/2D Plot type' %in% winMenuNames()){   	
    #		winMenuAdd("$Graph2Popup/2D Plot type")
    #		winMenuAddItem("$Graph2Popup/2D Plot type", 'auto     da()', "da()")
    #		winMenuAddItem("$Graph2Popup/2D Plot type", 'contour  dr()', "dr()")
    #		winMenuAddItem("$Graph2Popup/2D Plot type", 'filled   drf()', "drf()")
    #		winMenuAddItem("$Graph2Popup/2D Plot type", 'image    di()', "di()")
    #		winMenuAdd("$Graph2Popup/1D Plot type")
    #		winMenuAddItem("$Graph2Popup/1D Plot type", 'line', 
    #				"setGraphics(type='auto', refresh.graphics=TRUE)")
    #		winMenuAddItem("$Graph2Popup/1D Plot type", 'points', 
    #				"setGraphics(type='p', refresh.graphics=TRUE)")
    #		winMenuAddItem("$Graph2Popup/1D Plot type", 'both', 
    #				"setGraphics(type='b', refresh.graphics=TRUE)")
    #######TSB
    winMenuAdd("$Graph2Popup/Amino Acid Detail")
    winMenuAddItem("$Graph2Popup/Amino Acid Detail", 'Toggle Amino Acid Plotting    pAA()', "pAA()") 
    
    #######TSB
    #		winMenuAdd("$Graph2Popup/2D Contours")
    #		winMenuAddItem("$Graph2Popup/2D Contours", 'raise   ctu()', "ctu()")
    #		winMenuAddItem("$Graph2Popup/2D Contours", 'lower  ctd()', "ctd()")
    #		winMenuAdd("$Graph2Popup/1D Slices")
    #		winMenuAddItem("$Graph2Popup/1D Slices", 'direct slice    vs(1)', "vs(1)")
    #		winMenuAddItem("$Graph2Popup/1D Slices", 'indirect slice  vs(2)', "vs(2)")
    #		winMenuAddItem("$Graph2Popup/1D Slices", 'projection      pjv()', "pjv()")		
    #		winMenuAdd("$Graph2Popup/1D position")
    #		winMenuAddItem("$Graph2Popup/1D position", 'raise   vpu()', "vpu()")
    #		winMenuAddItem("$Graph2Popup/1D position", 'lower  vpd()', "vpd()")
    winMenuAdd("$Graph2Popup/Zoom")
    winMenuAddItem("$Graph2Popup/Zoom", 'in                zi()', "zi()")
    winMenuAddItem("$Graph2Popup/Zoom", 'out            zo()', "zo()")
    winMenuAddItem("$Graph2Popup/Zoom", '-', "none")
    winMenuAddItem("$Graph2Popup/Zoom", 'center       zc()', "zc()")
    winMenuAddItem("$Graph2Popup/Zoom", 'full             zf()', "zf()")
    winMenuAddItem("$Graph2Popup/Zoom", '--', "none")
    winMenuAddItem("$Graph2Popup/Zoom", 'hand         zz()', "zz()")
    winMenuAddItem("$Graph2Popup/Zoom", 'point         pz()', "pz()")
    winMenuAddItem("$Graph2Popup/Zoom", '---', "none")
    winMenuAddItem("$Graph2Popup/Zoom", 'previous    zp()', "zp()")
    #		winMenuAddItem("$Graph2Popup/Zoom", 'get shifts  loc()', "loc()")
  }
  
  ##creates context menu for the sub plot window
  if (dev == 'sub' && !'$Graph3Popup/Select' %in% winMenuNames()){   	
    winMenuAdd("$Graph3Popup/Select")
    winMenuAddItem("$Graph3Popup/Select", 'from list  rs(1)', "rs(1)")
    winMenuAddItem("$Graph3Popup/Select", 'from window  rs(3)', "rs(3)")
    winMenuAdd("$Graph3Popup/Edit")
    winMenuAddItem("$Graph3Popup/Edit", 'edit table  re()', "re()")
    winMenuAddItem("$Graph3Popup/Edit", 'delete ROI rDel()', "rDel()")
    winMenuAdd("$Graph3Popup/Plot")
    winMenuAddItem("$Graph3Popup/Plot", 'Replot rvs()', "rvs()")
  }
  
  ##creates context menu for the multiple file window
  if (dev == 'multi' && !'$Graph4Popup/Select' %in% winMenuNames()){   	
    winMenuAdd("$Graph4Popup/Select")
    winMenuAddItem("$Graph4Popup/Select", 'from list  rs(1)', "rs(1)")
    winMenuAddItem("$Graph4Popup/Select", 'from window  rs(4)', "rs(4)")
    winMenuAddItem("$Graph4Popup/Select", 'files  rsf()', "rsf()")
    winMenuAdd("$Graph4Popup/Edit")
    winMenuAddItem("$Graph4Popup/Edit", 'sort files  fs()', "fs()")
    winMenuAddItem("$Graph4Popup/Edit", 'edit table  re()', "re()")
    winMenuAddItem("$Graph4Popup/Edit", 'delete ROI  rDel()', "rDel()")
    winMenuAdd("$Graph4Popup/Summary")
    winMenuAddItem("$Graph4Popup/Summary", 'rSum()', "rSum()")
    winMenuAdd("$Graph4Popup/Plot")
    winMenuAddItem("$Graph4Popup/Plot", 'Replot rvm()', "rvm()")
  }
}

# Displays the digestR splash screen
splashScreen <- function(){

  par(mar=defaultSettings$mar, cex.axis=defaultSettings$cex.axis, 
      cex.main=defaultSettings$cex.main, bg='black')
  colMain <- '#b4d0f3'
  colBack <- '#0065ca'
  plot(0, 0, type='n', xlab='', ylab='', col.axis='black')
  # Define your color palette
  #colMain <- '#FF5733'  # A warm main color
  #colBack <- '#2E86C1'  # A contrasting background color
  
  # Letter positions and colors
  letters <- c('D', 'I', 'G', 'E', 'S', 'T', 'R')
  colors <- c(colMain, colMain, colMain, colMain, colMain, colMain, colBack)
  cex_values <- c(7, 7, 7, 7, 7, 7, 6.5)
  offset_values <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.4)
  
  # Loop to create and position letters
  for (i in 1:length(letters)) {
    text(-0.75 + (i-1)*0.25, 0.2, letters[i], col=colors[i], cex=cex_values[i], pos=3, offset=offset_values[i])
  }
  
  # Add a decorative 'R' in a different color
  text(0.70, 0.2, 'R', col='#E74C3C', cex=5.5, pos=3, offset=0.4)
  
  # Other text elements (modify as needed)
  text(0, 0.08, 'Digestomics Analyzer', col=colMain, cex=2.5, font=1)
  text(0, -0.15, paste('version 1.0.0', pkgVar$version), col=colMain, font=3)
  text(0, -0.25, 'gp() - Generate New Proteome', col=colMain)
  text(0, -0.35, 'pm() - Process Mascot files', col=colMain)
  text(0, -0.45, 'fo() - Open *.dcf files', col=colMain)

  # Force the graphics device to refresh
  dev.flush()
}

##################################################################################
# New splashScreen to display when loading DigestR


splashScreen1 <- function(){
library(png)
  # Load the image and get its dimensions
  img_path <- system.file("extdata", "DigestRpicture.png", package = "digestR")
  img <- readPNG(img_path)
  plot(0, 0, type = 'n', xlab = '', ylab = '', axes = FALSE, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))
  # Display the image properly (full window)
  rasterImage(img, -1.5, -1.5, 1.5, 1.5)

  # Letter positions and modern color scheme for "PANDAS"
  letters <- c('P', 'A', 'N', 'D', 'A', 'S')
  x_positions <- c(-0.9, -0.6, -0.3, 0.0, 0.3, 0.6)  # Custom X positions for each letter
  y_positions <- c(0.62, 0.62, 0.62, 0.62, 0.62, 0.62)  # Same Y position for all
  colors <- rep('#4da6ff', 6)  # Light blue for all letters
  cex_values <- rep(5, 6)  # Font size for letters
  
  # Add shadow effect to the letters for depth
  shadow_col <- '#00000050'  # Semi-transparent black shadow
  for (i in 1:length(letters)) {
    text(x_positions[i] + 0.02, y_positions[i] - 0.02, letters[i], col=shadow_col, cex=cex_values[i], pos=3)
    text(x_positions[i], y_positions[i], letters[i], col=colors[i], cex=cex_values[i], pos=3)
  }

  # Add black shadow for the title "Peptide Analyzer of Naturally Digested Amino acid Sequences"
  text(0, -0.4, 'Peptide Analyzer of', col='black', cex=1.80, font=2.3)  # Black shadow
  text(0, -0.4, 'Peptide Analyzer of', col='#4da6ff', cex=1.8, font=2)  # Main text
  
  # Shadow for "Naturally Digested Amino acid Sequences"
  text(0, -0.53, 'Naturally Digested Amino acid Sequences', col='black', cex=1.80, font=2.3)  # Black shadow
  text(0, -0.53, 'Naturally Digested Amino acid Sequences', col='#4da6ff', cex=1.8, font=2)  # Main text
  
  # Version and command text, shifted further downwards
  version <- '1.0.0'
  text(0, -0.75, paste('version', version), col = '#00b3b3', font = 2)  # Dynamic color # Dynamic color
  
  # Clean, minimal function list with better spacing and shifted downwards
  text(0, -0.86, 'gp() - Generate New Proteome', col='#00b3b3', cex=1.1, font = 2)  # Light blue
  text(0, -0.94, 'pm() - Process Mascot files', col='#00b3b3', cex=1.1, font = 2)  # Light blue
  text(0, -1.02, 'fo() - Open *.dcf files', col='#00b3b3', cex=1.1, font = 2)  # Light blue
  
  # Force the graphics device to refresh
  dev.flush()
}

##################################################################################

about <- function(){
  ##creates toplevel
  dlg <- tktoplevel()
  tkwm.title(dlg, 'About digestR')
  tkwm.resizable(dlg, FALSE, FALSE)
  tkwm.deiconify(dlg)
  tcl('wm', 'attributes', dlg, topmost=TRUE)
  
  ##display digestR package info
  msg <- paste('digestR version', pkgVar$version, '\n',
               'Copyright (C) 2024 Dimitri Desmonts de Lamache, Travis Bingeman, Raied Aburashed, Sören Wacker and Ian A. Lewis' )
  msgLabel <- ttklabel(dlg, text=msg)
  
  ##creates ok button
  okButton <- ttkbutton(dlg, text='OK', width=10, command=function(...)
    tkdestroy(dlg))
  
  ##creates release notes button
  onRelease <- function(){
    myHelp('release_notes', TRUE)
    tkdestroy(dlg)
  }
  relButton <- ttkbutton(dlg, text='Release Notes', width=13,
                         command=onRelease)
  
  ##add widgets to toplevel
  tkgrid(msgLabel, column=1, columnspan=2, row=1, sticky='w', pady=c(12, 0), 
         padx=10)
  tkgrid(okButton, column=1, row=2, padx=c(6, 4), pady=c(6, 10), sticky='e')
  tkgrid(relButton, column=2, row=2, padx=c(4, 6), pady=c(6, 10), sticky='w')
  
  ##allow users to press the enter key to make selections
  tkbind(dlg, '<Return>', function(...) tryCatch(tkinvoke(tkfocus()), 
                                                 error=function(er){}))
  tkfocus(okButton)
  
  invisible()
}

## Internal function for displaying HTML help pages
## page - character string; name of the help page to open
## docDir - logical; searches the digestR 'doc' directory for the HTML help file
## 					 if TRUE, otherwise the default HTML directory for digestR is used.
myHelp <- function(page, docDir=FALSE){
  
  if (page == 'user_manual')
    page <- paste(page, '.pdf', sep='', collapse='')
  else
    page <- paste(page, '.html', sep='', collapse='')
  if (docDir){
    fileDir <- system.file('doc', package='digestR')
    msg <- 'Could not find specified help page.'
  }else{
    fileDir <- system.file('html', package='digestR')
    msg <- paste('Could not find specified help page.\n', 
                 'Try entering "?function_name" in the R console.', sep='')
  }
  tryCatch(browseURL(paste('file:///', file.path(fileDir, page), sep='', 
                           collapse='')), error=function(er) myMsg(msg, icon='error'))
  
  return(invisible())
}


################################################################################
##                                                                            ##
##                  Spectral (binary) file reading utilities                  ##
##                                                                            ##
################################################################################

## Internal function ucsfHead
## Reads a sparky format file and returns a header and other file info 
## file.name - File to be read, NULL will open a file selection window
## print.info - prints a summary of the data if TRUE
## returns a header for a sparky format file
ucsfHead <- function (file.name = NULL, print.info = TRUE){
  
  ## Get the path to the file
  if(is.null(file.name))
    file.name <- myOpen(title = "Select a spectrum", multiple = FALSE)
  if(file.name == '')
    return(invisible())
  
  ## Open connection to binary file and check the file format
  con <- myFile(file.name, "rb")
  fileType <- readLines(con, n=1, warn=FALSE)
  if (fileType == "RSD_NMR"){
    close(con)
    return(rsdHead(file.name, print.info))
  }
  if(fileType != "UCSF NMR"){
    close(con)
    stop(paste('digestR only accepts UCSF (sparky) format data.',
               'Use cf() to convert spectra to the correct format.',
               'If you want to load an R workspace, use wl().', sep='\n'), 
         call.=FALSE)
  }
  
  ## Make sure data is modern sparky format and is 1D, 2D, or 3D
  seek(con, where=10, origin = "start")
  nDim <- readBin(con, size=1, what='integer') 
  if(nDim != 1 && nDim != 2 && nDim != 3){
    close(con)
    stop('digestR only supports 1D, 2D, and 3D sparky format data', call.=FALSE)
  }
  if(readBin(con, size = 1, what = 'integer') != 1){
    close(con)
    stop('digestR does not support complex data', call.=FALSE)
  }
  seek(con, where = 13, origin = "start")
  if(readBin(con, size=1, what='integer') != 2 ){ 
    close(con)
    stop('digestR only supports the current version of sparky formatting', 
         call.=FALSE)
  }
  
  ## Get file information
  fileInfo <- file.info(file.name)
  fsize <- fileInfo$size
  mtime <- fileInfo$mtime
  
  ## Set points for the start of each header
  headStart <- 180
  
  ## Set endian format -- sort of a bonehead way of doing it...
  endFormat <- 'big'
  seek(con, where=(headStart + 20), origin = "start")
  w2Frq <- readBin(con, size=4, what='double', endian=endFormat) 
  
  ## Switch endian format if the frequency is wierd
  if( is.na(w2Frq) || w2Frq < 1 || w2Frq > 1500 )
    endFormat <- 'little'
  
  ## Read noise estimate if present
  seek(con, where=14, origin="start")
  noiseText <- readBin(con, size=8, what='character')
  noiseEst <- NULL
  if (noiseText == 'noiseEst')
    noiseEst <- readBin(con, size=4, what='double', endian=endFormat)
  
  ## Read first header
  if( nDim == 1 ){
    axis <- 'w2'
  }else
    axis <- 'w1'
  
  ## Set nucleus name (1H, 13C, 15N, 31P, ...)
  seek(con, where=headStart, origin = "start")
  nuc <- readBin(con, size=6, what='character')
  
  ## number of data points in this axis
  seek(con, where=(headStart + 8), origin = "start")
  msize <- readBin(con, size=4, what='integer', endian=endFormat)
  
  ## tile size along this axis
  seek(con, where=(headStart + 16), origin = "start")
  bsize <- readBin(con, size=4, what='integer', endian=endFormat)
  
  ##spectrometer frequency (MHz)
  seek(con, where=(headStart + 20), origin = "start")
  sFrq <- readBin(con, size=4, what='double', endian=endFormat)
  
  ## spectral width (Hz)
  seek(con, where=(headStart + 24), origin = "start")
  swHz <- readBin(con, size=4, what='double', endian=endFormat)
  
  ## center of data (ppm)
  seek(con, where= (headStart + 28), origin = "start")
  cFrqPPM <- readBin(con, size=4, what='double', endian=endFormat)
  
  ## up and downfield shifts for this axis
  seek(con, where=(headStart + 32), origin = "start")
  upShift <- readBin(con, size=4, what='double', endian=endFormat)
  downShift <- readBin(con, size=4, what='double', endian=endFormat)
  
  ## Start of NMR data for 1D
  binLoc <- headStart + 128
  
  ##Read the second header
  if( nDim > 1){
    
    headStart <- headStart + 128
    
    ## Start of NMR data for 2D
    binLoc <- headStart + 128 
    
    axis <- c(axis, 'w2')
    
    ## Set nucleus name (1H, 13C, 15N, 31P, ...)
    seek(con, where=headStart, origin = "start")
    nuc <- c(nuc, readBin(con, size=6, what='character'))
    
    ## number of data points in this axis
    seek(con, where=(headStart + 8), origin = "start")
    msize <- c(msize, readBin(con, size=4, what='integer', endian=endFormat))
    
    ## tile size along this axis
    seek(con, where=(headStart + 16), origin = "start")
    bsize <- c(bsize, readBin(con, size=4, what='integer', endian=endFormat))
    
    #spectrometer frequency (MHz)
    seek(con, where=(headStart + 20), origin = "start")
    sFrq <- c(sFrq, readBin(con, size=4, what='double', endian=endFormat))
    
    ## spectral width (Hz)
    seek(con, where=(headStart + 24), origin = "start")
    swHz <- c(swHz, readBin(con, size=4, what='double', endian=endFormat))
    
    ## center of data (ppm)
    seek(con, where=(headStart + 28), origin = "start")
    cFrqPPM <- c(cFrqPPM, readBin(con, size=4, what='double', endian=endFormat))
    
    ## up and downfield shifts for this axis
    seek(con, where=(headStart + 32), origin = "start")
    upShift <- c(upShift, readBin(con, size=4, what='double', 
                                  endian=endFormat))
    downShift <- c(downShift, readBin(con, size=4, what='double', 
                                      endian=endFormat))
  }
  
  ##Read the third header
  if( nDim == 3 ){
    
    headStart <- headStart + 128
    
    ## Start of NMR data for 2D
    binLoc <- headStart + 128 
    
    axis <- c('w1', 'w2', 'w3')
    
    ## Set nucleus name (1H, 13C, 15N, 31P, ...)
    seek(con, where=headStart, origin = "start")
    nuc <- c(nuc[2], readBin(con, size=6, what='character'), nuc[1])
    
    ## number of data points in this axis
    seek(con, where=(headStart + 8), origin = "start")
    msize <- c(msize[2], readBin(con, size=4, what='integer', endian=endFormat), 
               msize[1])
    
    ## tile size along this axis
    seek(con, where=(headStart + 16), origin = "start")
    bsize <- c(bsize[2], readBin(con, size=4, what='integer', endian=endFormat),
               bsize[1])
    
    #spectrometer frequency (MHz)
    seek(con, where=(headStart + 20), origin = "start")
    sFrq <- c(sFrq[2], readBin(con, size=4, what='double', endian=endFormat),
              sFrq[1])
    
    ## spectral width (Hz)
    seek(con, where=(headStart + 24), origin = "start")
    swHz <- c(swHz[2], readBin(con, size=4, what='double', endian=endFormat),
              swHz[1])
    
    ## center of data (ppm)
    seek(con, where=(headStart + 28), origin = "start")
    cFrqPPM <- c(cFrqPPM[2], readBin(con, size=4, what='double', 
                                     endian=endFormat), cFrqPPM[1])
    
    ## up and downfield shifts for this axis
    seek(con, where=(headStart + 32), origin = "start")
    upShift <- c(upShift[2], readBin(con, size=4, what='double', 
                                     endian=endFormat), upShift[1])
    downShift <- c(downShift[2], readBin(con, size=4, what='double', 
                                         endian=endFormat), downShift[1])
  }
  
  ## Find the range, interquartile range, and median of spectrum for noise est
  ## For 2D/3D data, this is based on the first tile
  seek(con, where = binLoc, origin = "start")
  if(nDim > 1){          
    file.range <- fivenum(readBin(con, size=4, what='double',
                                  endian = endFormat, n=(bsize[1] * bsize[2])))                                                
  }else
    file.range <- fivenum(readBin(con, size=4, what='double',
                                  endian = endFormat, n= msize[1])) 
  
  ## Close binary connection
  closeAllConnections()
  
  ## Translate sweep width to ppm
  swPPM <- swHz / sFrq
  
  ## Calculate the up and downfield shifts
  if (all(downShift == 0) && all(upShift == 0)){
    upPPM <- cFrqPPM - (swPPM / 2)
    downPPM <- cFrqPPM + (swPPM / 2)
  }else{
    
    ## Use the up and downfield shifts from the header if present
    ## This is used for files that have been converted from ASCII
    upPPM <- round(upShift, 5)
    downPPM <- round(downShift, 5)
  }
  
  ## Calculate noise estimate, if not present in main header
  if (is.null(noiseEst))
    noiseEst <- (diff(file.range[c(2, 4)]) / 2) / .674 ## 1 sd
  
  ## Make a new header with extracted data
  head <- list(
    file.name = file.name,
    file.size = fsize,
    date.modified = mtime,
    axis = axis,
    nucleus = nuc,
    matrix_size = msize,
    block_size = bsize,
    upfield_ppm = upPPM,
    downfield_ppm = downPPM,
    spectrum_width_Hz = swHz,
    transmitter_MHz = sFrq,
    center_ppm = cFrqPPM,
    binary_location = binLoc,
    endian = endFormat,
    number_dimensions = nDim,
    noise_est = noiseEst,
    min_intensity = file.range[1],
    max_intensity = file.range[5],
    zero_offset = file.range[3],
    user_title = basename(file.name)
  )
  if (nDim == 3)
    head$z_value <- upPPM[3]
  
  ## Print useful file info if print.info = TRUE
  if(print.info)
    log_message(data.frame(head[5:12]))
  
  ## Modify upfield PPM for digestR format 
  ## NOTE: Fourier transform shifts the observed center
  ##       frequency to the right in the frequency domain.
  ##       This small defect can cause problems when the number of points 
  ##       collected is very small. The code below corrects the problem
  ##       for most datasets. However, if a different Fourier algorithm is used
  ##       then the correction may be applied in the wrong direction.
  ##       A more robust method for determining this correction would be 
  ##       an obvious method for improving the accuracy of peak picking.
  if (all(downShift == 0) && all(upShift == 0)){
    cor <- ((head$downfield_ppm - head$upfield_ppm) / (msize - 1))	* 
      -(msize %% 2 - 1)
    head$upfield_ppm <- head$upfield_ppm + cor
  }
  return( list( file.par = head ) )
}


## Internal function rsdHead
## Reads an RSD file and returns header information 
## fileName - File to be read, NULL will open a file selection window
## print.info - prints a summary of the data if TRUE
## returns a header for a sparky format file
rsdHead <- function (file.name, print.info=TRUE){
  
  ## Get the path to the file
  if (missing(file.name))
    file.name <- myOpen(title='Select a spectrum', multiple=FALSE, 
                        filetypes=list(rsd="RSD File"))
  if (!nzchar(file.name))
    return(invisible())
  
  ## Open connection to binary file and check the file format
  con <- myFile(file.name, 'rb')
  fileType <- readBin(con, size=10, what='character')
  if (fileType != 'RSD_NMR'){
    close(con)
    stop('Specified file is not in RSD format.', call.=FALSE)
  }
  
  ## Read top-level header
  seek(con, where=20)
  nDim <- readBin(con, size=4, what='integer', endian='big')
  nblocks <- readBin(con, size=4, what='integer', endian='big')
  noiseEst <- readBin(con, size=4, what='double', endian='big')
  
  ## Read header information for original spectrum
  origNp <- origDown <- origUp <- sw <- sf <- nuc <- NULL
  for (i in 1:nDim){
    seek(con, where=(100 + (i - 1) * 50 + 4))
    origNp <- c(origNp, readBin(con, size=4, what='integer', endian='big'))
    origDown <- c(origDown, readBin(con, size=4, what='double', endian='big'))
    origUp <- c(origUp, readBin(con, size=4, what='double', endian='big'))
    sw <- c(sw, readBin(con, size=4, what='double', endian='big'))
    sf <- c(sf, readBin(con, size=4, what='double', endian='big'))
    nuc <- c(nuc, readBin(con, size=10, what='character'))
  }
  
  ## Read block headers
  np <- downShifts <- upShifts <- as.list(1:nDim)
  seek(con, where=(100 + nDim * 50))
  for (i in 1:nDim){
    for (j in 1:nblocks){
      seek(con, where=4, origin='current')
      np[[i]][j] <- readBin(con, size=4, what='integer', endian='big')
      downShifts[[i]][j] <- readBin(con, size=4, what='double',	endian='big')
      upShifts[[i]][j] <- readBin(con, size=4, what='double', endian='big')
      if (!(j == nblocks && i == nDim))
        seek(con, where=14, origin='current')
    }
  }
  
  ## Get file information
  binLoc <- 100 + nDim * 50 + nDim * nblocks * 30
  fileInfo <- file.info(file.name)
  fsize <- fileInfo$size
  mtime <- fileInfo$mtime
  
  ## Format values
  if (nDim == 1){
    axis <- 'w2'
    totalPoints <- sum(np[[1]])
  }else{
    axis <- c('w1', 'w2')
    totalPoints <- sum(np[[1]] * np[[2]])
  }
  names(np) <- names(downShifts) <- names(upShifts) <- axis
  
  ## Find the range, interquartile range, and median of spectrum for noise est
  seek(con, where=binLoc, origin='start')
  file.range <- fivenum(readBin(con, size=4, what='double',	endian='big', 
                                n=totalPoints)) 
  close(con)
  
  ## Make a new header with extracted data
  head <- list(file.name=file.name,
               file.size=fsize,
               date.modified=mtime,
               user_title=basename(file.name),
               axis=axis,
               nucleus=nuc,
               matrix_size=origNp,
               upfield_ppm=origUp,
               downfield_ppm=origDown,
               block_size=np,
               block_upfield_ppms=upShifts,
               block_downfield_ppms=downShifts,
               spectrum_width_Hz=sw,
               transmitter_MHz=sf,
               binary_location=binLoc,
               endian='big',
               number_dimensions=nDim,
               noise_est=noiseEst,
               min_intensity=file.range[1],
               max_intensity=file.range[5],
               zero_offset=file.range[3]
  )
  
  ## Print useful file info if print.info = TRUE
  if (print.info)
    log_message(data.frame(head[5:9]))
  
  return(list(file.par=head))
}


## Internal function ucsf1D
## Reads sparky format 1D spectra and returns data within specified range
## file.name - File to be read, NULL will open a file selection window
## w2Range   - Numeric list of length 2, chemical shift range to be read
## file.par  - The header returned from ucsfHead can be passed to ucsf1D,
##             this speeds up graphics by eliminating noise estimates
## returns   - List object with file header, w2 shifts, and data
ucsf1D <- function(file.name = NULL, w2Range = NULL, file.par = NULL){
  
  ## Get the path to the file and ucsf header
  if(is.null(file.name))
    file.name <- myOpen(title = "Select a spectrum", multiple = FALSE)
  if(file.name == '')
    return(invisible())
  if( is.null(file.par) ){
    if (file.name %in% names(fileFolder))
      file.par <- fileFolder[[file.name]]$file.par
    file.par <- ucsfHead(file.name = file.name, print.info = FALSE)[[1]]
  }
  
  ## Redirect to rsd1D for RSD format spectra
  if (!is.null(file.par$block_upfield_ppms)){
    outFolder <- rsd1D(file.name, w2Range, file.par)
    return(outFolder)
  }
  
  ## Setup outgoing list 
  outFolder <- list( file.par = file.par, w2 = seq(file.par$upfield_ppm[1], 
                                                   file.par$downfield_ppm[1], length.out = file.par$matrix_size[1]), 
                     data = NULL )
  
  ## Make sure valid w2 ranges were submitted
  out <- list( w2 = NULL )
  if( is.null(w2Range) || length(w2Range) != 2 || !is.numeric(w2Range) ){
    w2Range <- c(file.par$upfield_ppm[1], file.par$downfield_ppm[1])
    out$w2 <- c(1,length(outFolder$w2))
  }else{
    w2Range <- sort(w2Range)	
    t1 <- findInterval(w2Range, outFolder$w2, all.inside = TRUE)	
    t2 <- t1 + 1
    for(i in 1:2)
      out$w2[i] <- switch( which.min(c(
        abs(w2Range[i] - outFolder$w2[t1[i]]),
        abs(w2Range[i] - outFolder$w2[t2[i]]))), t1[i], t2[i])
  }
  
  ## Read data from memory if entry exists in fileFolder
  if (file.name %in% names(fileFolder) && 
      !is.null(fileFolder[[file.name]]$data))
    outFolder$data <- fileFolder[[file.name]]$data
  else{
    
    ## Read binary data
    endFormat <- outFolder$file.par$endian 
    con <- myFile(file.name, "rb")
    seek(con, where=outFolder$file.par$binary_location, origin = "start")
    outFolder$data <- readBin(con, size=4, what='double',
                              endian = endFormat, n=outFolder$file.par$matrix_size[1])
    closeAllConnections()
  }
  
  ## Trim data to fit w2 ranges
  outFolder$w2 <- outFolder$w2[out$w2[1]:out$w2[2]]
  
  ## Invert selections to match binary
  out$w2 <- sort(outFolder$file.par$matrix_size[1] - out$w2) + 1	
  outFolder$data <- outFolder$data[out$w2[1]:out$w2[2]]
  
  return(outFolder)
}


## Internal function ucsf2D
## Reads sparky format spectra and returns all of the sparky tiles covered
##    by the chemical shift range provided.
## file.name - name of spectrum to be read, NULL opens a file slection window
## w1Range = Min and Max chemical shifts to be read in the indirect dimension,
##           NULL will return all values.
## w2Range = Min and Max chemical shifts to be read in the direct dimension
##           NULL will return all values.
## file.par - file header (from ucsfHead()), NULL will read the header first
## notes: providing a file.par from memory is used to speed up graphics
## returns the designated region of a spectrum and associated parameters
ucsf2D <- function ( file.name = NULL, w1Range=NULL, w2Range=NULL,
                     file.par=NULL){
  
  ## Get the path to the file and ucsf header
  if(is.null(file.name))
    file.name <- myOpen(title = "Select a spectrum", multiple = FALSE)
  if(file.name == '')
    return(invisible())
  if( is.null(file.par) ){
    if (file.name %in% names(fileFolder))
      file.par <- fileFolder[[file.name]]$file.par
    file.par <- ucsfHead(file.name = file.name, print.info = FALSE)[[1]]
  }
  
  ## Redirect to ucsf1D if 1D file is opened
  if( file.par$number_dimensions == 1 )
    return(ucsf1D( file.name=file.name, w2Range=w2Range, file.par=file.par ))
  
  ## Redirect to rsd2D if RSD file is opened
  if( !is.null(file.par$block_upfield_ppms) )
    return(rsd2D( file.name=file.name, w1Range=w1Range, w2Range=w2Range, 
                  file.par=file.par ))
  
  ## Redirect to ucsf3D if 3D file is opened
  if( file.par$number_dimensions == 3 )
    return(ucsf3D( file.name=file.name, w1Range=w1Range, w2Range=w2Range, 
                   file.par=file.par ))
  
  ## Setup outgoing list 
  outFolder <- list( file.par = file.par, 
                     w1 = seq(file.par$upfield_ppm[1], file.par$downfield_ppm[1], 
                              length.out = file.par$matrix_size[1]),
                     w2 = seq(file.par$upfield_ppm[2], file.par$downfield_ppm[2], 
                              length.out = file.par$matrix_size[2]), data = NULL )
  
  ## Make sure valid w2 ranges were submitted
  out <- list( w1 = NULL, w2 = NULL )
  if( is.null(w1Range) || length(w1Range) != 2 || !is.numeric(w1Range) )
    w1Range <- c(file.par$upfield_ppm[1], file.par$downfield_ppm[1])
  if( is.null(w2Range) || length(w2Range) != 2 || !is.numeric(w2Range) )
    w2Range <- c(file.par$upfield_ppm[2], file.par$downfield_ppm[2])
  
  ## Find best w1/w2 matches	
  w1Range <- sort(w1Range); w2Range <- sort(w2Range)
  t1 <- findInterval(w1Range, outFolder$w1, all.inside = TRUE)
  d1 <- findInterval(w2Range, outFolder$w2, all.inside = TRUE)
  t2 <- t1 + 1; d2 <- d1 + 1
  for(i in 1:2){
    out$w1[i] <- switch( which.min(c(
      abs(w1Range[i] - outFolder$w1[t1[i]]),
      abs(w1Range[i] - outFolder$w1[t2[i]]))), t1[i], t2[i])
    out$w2[i] <- switch( which.min(c(
      abs(w2Range[i] - outFolder$w2[d1[i]]),
      abs(w2Range[i] - outFolder$w2[d2[i]]))), d1[i], d2[i])
  }
  
  ## Trim w1/w2 Ranges
  outFolder$w1 <- outFolder$w1[ (out$w1[1] : out$w1[2] )]
  outFolder$w2 <- outFolder$w2[ (out$w2[1] : out$w2[2] )]
  
  ## Invert w1/w2 selection to match binary data format
  out$w1 <- sort(file.par$matrix_size[1] - out$w1) + 1
  out$w2 <- sort(file.par$matrix_size[2] - out$w2) + 1
  
  ## Read and trim data from memory if entry exists in fileFolder
  if (file.name %in% names(fileFolder) && 
      !is.null(fileFolder[[file.name]]$data)){
    outFolder$data <- fileFolder[[file.name]]$data[out$w2[1]:out$w2[2], 
                                                   out$w1[1]:out$w1[2]]
    return(outFolder)
  }
  
  ## Find sparky tiles that will be plotted
  w1Tiles <- (ceiling(out$w1[1] / file.par$block_size[1]):
                ceiling(out$w1[2] / file.par$block_size[1]))
  w2Tiles <- (ceiling(out$w2[1] / file.par$block_size[2]):
                ceiling(out$w2[2] / file.par$block_size[2]))
  tiles <- NULL
  for ( i in 1:length(w1Tiles))
    tiles <- c(tiles, ((w1Tiles[i]-1) * (ceiling(file.par$matrix_size[2] / 
                                                   file.par$block_size[2])) + (w2Tiles-1)))
  
  ## Reset the match index to the active tiles
  out$w1 <- out$w1 - (w1Tiles[1] -1 ) * file.par$block_size[1]
  out$w2 <- out$w2 - (w2Tiles[1] -1 ) * file.par$block_size[2]
  gc(FALSE)
  
  ## Open connection to binary file
  w2TileNum <- ceiling(file.par$matrix_size[2] / file.par$block_size[2] )
  endFormat <- file.par$endian
  outFolder$data <- matrix( rep( NA, file.par$block_size[1] * length(w1Tiles) * 
                                   file.par$block_size[2] * length(w2Tiles)), 
                            ncol = file.par$block_size[1] * length(w1Tiles))
  con <- myFile(file.name, "rb")
  seek(con, where = file.par$binary_location, origin = "start")
  
  ##Read binary data for each tile
  j <- -1
  w1Count <- w2Count <- 1
  for(i in tiles){
    w1R <- (1:file.par$block_size[1]) + file.par$block_size[1] * (w1Count -1)
    w2R <- (1:file.par$block_size[2]) + file.par$block_size[2] * (w2Count -1)
    w2Count <- w2Count + 1
    if(w2Count > length(w2Tiles)){
      w2Count <- 1
      w1Count <- w1Count + 1
    }
    
    ## read binary
    seek(con, where = (file.par$block_size[1] *
                         file.par$block_size[2] * 4 * (i - j - 1)), origin = "current")
    outFolder$data[w2R, w1R] <- readBin(con, size=4, what='double',
                                        endian = endFormat, 
                                        n=(file.par$block_size[1] * file.par$block_size[2]))
    j <- i
  }
  
  ## Close binary conection
  closeAllConnections()
  
  ## Trim data to fit w1/w2 ranges
  outFolder$data <- outFolder$data[out$w2[1]:out$w2[2], out$w1[1]:out$w1[2]]
  gc(FALSE)
  
  return(outFolder)
}


## Internal function ucsf3D
## Reads sparky format spectra and returns all of the sparky tiles covered
##    by the chemical shift range provided.
## file.name - name of spectrum to be read, NULL opens a file slection window
## w1Range - Min and Max chemical shifts to be read in the indirect dimension,
##           NULL will return all values.
## w2Range - Min and Max chemical shifts to be read in the direct dimension
##           NULL will return all values.
## w3Range - chemical shift in the z dimension, must be a single value
## file.par - file header (from ucsfHead()), NULL will read the header first
## notes: providing a file.par from memory is used to speed up graphics
## returns the designated region of a spectrum and associated parameters
ucsf3D <- function(file.name=NULL, w1Range=NULL, w2Range=NULL, w3Range=NULL, 
                   file.par=NULL){
  
  ## Get the path to the file and ucsf header
  if (is.null(file.name))
    file.name <- myOpen(title="Select a spectrum", multiple=FALSE)
  if (file.name == '')
    return(invisible())
  if( is.null(file.par) ){
    if (file.name %in% names(fileFolder))
      file.par <- fileFolder[[file.name]]$file.par
    file.par <- ucsfHead(file.name = file.name, print.info = FALSE)[[1]]
  }
  
  ## Define some local variables
  bs <- file.par$block_size
  ms <- file.par$matrix_size
  uf <- file.par$upfield_ppm
  df <- file.par$downfield_ppm
  endFormat <- file.par$endian
  binLoc <- file.par$binary_location
  
  ## Setup outgoing list 
  outFolder <- list(file.par=file.par, w1=seq(uf[1], df[1], length.out=ms[1]),
                    w2=seq(uf[2], df[2], length.out=ms[2]), w3=seq(uf[3], df[3], 
                                                                   length.out=ms[3]), data=NULL)
  
  ## Make sure valid w2 ranges were submitted
  if (is.null(w1Range) || length(w1Range) != 2 || !is.numeric(w1Range))
    w1Range <- c(file.par$upfield_ppm[1], file.par$downfield_ppm[1])
  if (is.null(w2Range) || length(w2Range) != 2 || !is.numeric(w2Range))
    w2Range <- c(file.par$upfield_ppm[2], file.par$downfield_ppm[2])
  if (is.null(w3Range) || length(w3Range) != 1 || !is.numeric(w3Range))
    w3Range <- rep(file.par$z_value, 2)
  
  ## Find best w1/w2 matches	
  w1Range <- sort(w1Range); w2Range <- sort(w2Range)
  t1 <- findInterval(w1Range, outFolder$w1, all.inside=TRUE)
  d1 <- findInterval(w2Range, outFolder$w2, all.inside=TRUE)
  z1 <- findInterval(w3Range, outFolder$w3, all.inside=TRUE)
  t2 <- t1 + 1; d2 <- d1 + 1; z2 <- z1 + 1
  out <- NULL
  for(i in 1:2){
    out$w1[i] <- switch(which.min(c(abs(w1Range[i] - outFolder$w1[t1[i]]),
                                    abs(w1Range[i] - outFolder$w1[t2[i]]))), t1[i], t2[i])
    out$w2[i] <- switch(which.min(c(abs(w2Range[i] - outFolder$w2[d1[i]]),
                                    abs(w2Range[i] - outFolder$w2[d2[i]]))), d1[i], d2[i])
    out$w3[i] <- switch(which.min(c(abs(w3Range[i] - outFolder$w3[z1[i]]),
                                    abs(w3Range[i] - outFolder$w3[z2[i]]))), z1[i], z2[i])
  }
  
  ## Trim w1/w2 Ranges
  outFolder$w1 <- outFolder$w1[(out$w1[1]:out$w1[2])]
  outFolder$w2 <- outFolder$w2[(out$w2[1]:out$w2[2])]
  outFolder$w3 <- outFolder$w3[(out$w3[1]:out$w3[2])]
  
  ## Invert w1/w2 selection to match binary data format
  out$w1 <- sort(ms[1] - out$w1) + 1
  out$w2 <- sort(ms[2] - out$w2) + 1
  
  ## Read and trim data from memory if entry exists in fileFolder
  if (file.name %in% names(fileFolder) && 
      !is.null(fileFolder[[file.name]]$data)){
    outFolder$data <- fileFolder[[file.name]]$data[out$w2[1]:out$w2[2], 
                                                   out$w1[1]:out$w1[2]]
    return(outFolder)
  }
  
  ## Find sparky tiles that will be plotted
  w1Tiles <- (ceiling(out$w1[1] / bs[1]):ceiling(out$w1[2] / bs[1]))
  w2Tiles <- (ceiling(out$w2[1] / bs[2]):ceiling(out$w2[2] / bs[2]))
  w3Tiles <- (ceiling(out$w3[1] / bs[3]):ceiling(out$w3[2] / bs[3]))
  tiles <- NULL
  numTiles <- ceiling(ms / bs)
  for (i in w1Tiles)
    tiles <- c(tiles, (i - 1) * numTiles[2] + (w2Tiles - 1))
  tiles <- tiles + (w3Tiles - 1) * numTiles[1] * numTiles[2]
  
  
  ## Reset the match index to the active tiles
  out$w1 <- out$w1 - (w1Tiles[1] - 1) * bs[1]
  out$w2 <- out$w2 - (w2Tiles[1] - 1) * bs[2]
  out$w3 <- out$w3 - (w3Tiles[1] - 1) * bs[3]
  
  ## Open connection to binary file
  w2TileNum <- ceiling(ms[2] / bs[2] )
  outFolder$data <- matrix(NA, nrow=bs[2] * length(w2Tiles), 
                           ncol=bs[1] * length(w1Tiles))
  con <- myFile(file.par$file.name, "rb")
  
  ##Read binary data for each tile
  w1Count <- w2Count <- 1
  tileSize <- bs[1] * bs[2] * bs[3] * 4
  for (i in tiles){
    
    ## Define current w1/w2 range
    w1R <- (1:bs[1]) + bs[1] * (w1Count - 1)
    w2R <- (1:bs[2]) + bs[2] * (w2Count - 1)
    w2Count <- w2Count + 1
    if (w2Count > length(w2Tiles)){
      w2Count <- 1
      w1Count <- w1Count + 1
    }
    
    ## Define data location for current tile
    tileLoc <- tileSize * i + binLoc
    if (bs[3] > 1)
      zPos <- (out$w3[1] - 1) * bs[3]
    else
      zPos <- 0
    dataLoc <- tileLoc + bs[1] * bs[2] * 4 * zPos
    
    ## read binary
    seek(con, dataLoc, origin='start')
    outFolder$data[w2R, w1R] <- readBin(con, size=4, what='double', 
                                        endian=endFormat,	n=bs[1] * bs[2])
  }
  close(con)
  
  ## Trim data to fit shift ranges
  outFolder$data <- outFolder$data[out$w2[1]:out$w2[2], out$w1[1]:out$w1[2]]
  
  return(outFolder)
} 


## Internal function rsd1D
## Reads RSD format 1D spectra and returns data within specified range
## file.name - name of spectrum to be read, NULL opens a file slection window
## w2Range - Numeric list of length 2, chemical shift range to be read
## file.par - file header (from rsdHead()), NULL will read the header first
## notes: providing a file.par from memory is used to speed up graphics
## returns the designated region of a spectrum and associated parameters
rsd1D <- function(file.name=NULL, w2Range=NULL, file.par=NULL){
  
  ## Get the path to the file and RSD header
  if (is.null(file.name))
    file.name <- myOpen(title='Select a spectrum', multiple=FALSE)
  if (!length(file.name) || !nzchar(file.name))
    return(invisible())
  if( is.null(file.par) )
    file.par <- rsdHead(file.name=file.name, print.info=FALSE)[[1]]
  
  ## Setup outgoing list 
  outFolder <- list(file.par=file.par, w2=seq(file.par$upfield_ppm[1], 
                                              file.par$downfield_ppm[1], length.out=file.par$matrix_size[1]), 
                    data=NULL)
  
  ## Make sure a valid w2 range was submitted
  if (is.null(w2Range) || length(w2Range) != 2 || !is.numeric(w2Range))
    w2Range <- c(file.par$upfield_ppm[1], file.par$downfield_ppm[1])
  
  ## Find best w2 matches
  w2Range <- sort(w2Range)	
  t1 <- findInterval(w2Range, outFolder$w2, all.inside=TRUE)	
  t2 <- t1 + 1
  out <- list(w2=NULL)
  for (i in 1:2)
    out$w2[i] <- switch(which.min(c(abs(w2Range[i] - outFolder$w2[t1[i]]),
                                    abs(w2Range[i] - outFolder$w2[t2[i]]))), t1[i], t2[i])
  winW2 <- outFolder$w2[(out$w2[1]:out$w2[2])]
  
  ## Find tiles that will be plotted
  tiles <- NULL
  for (tNum in seq_along(file.par$block_size$w2)){
    
    ## Get the chemical shifts for the current tile
    blockW2 <- seq(file.par$block_upfield_ppms$w2[tNum], 
                   file.par$block_downfield_ppms$w2[tNum], 
                   length.out=file.par$block_size$w2[tNum])
    
    ## Check the window for the presence of any shift in the current block
    if (any(round(blockW2, 3) %in% round(winW2, 3)))
      tiles <- c(tiles, tNum)
  }
  
  ## Define data locations for each block
  blockLoc <- file.par$binary_location
  for (i in seq_along(file.par$block_size$w2))
    blockLoc <- c(blockLoc, blockLoc[i] + 4 * file.par$block_size$w2[i])
  
  ## Read binary data for each tile
  outFolder$data <- rep(0, length(outFolder$w2))
  con <- myFile(file.name, 'rb')
  for (i in tiles){
    
    ## Find best w2 matches for current tile
    w2Range <- c(file.par$block_upfield_ppms$w2[i], 
                 file.par$block_downfield_ppms$w2[i])
    t1 <- findInterval(w2Range, outFolder$w2, all.inside=TRUE)
    t2 <- t1 + 1
    tile <- NULL
    for (j in 1:2)
      tile$w2[j] <- switch(which.min(c(abs(w2Range[j] - outFolder$w2[t1[j]]),
                                       abs(w2Range[j] - outFolder$w2[t2[j]]))), t1[j], t2[j])
    
    ## Read data for current tile
    seek(con, blockLoc[i], origin='start')
    outFolder$data[tile$w2[1]:tile$w2[2]] <- rev(readBin(con, size=4, 
                                                         what='double', endian='big',	n=file.par$block_size$w2[i]))
  }
  close(con)
  
  ## Trim data to fit w2 ranges
  outFolder$w2 <- winW2
  outFolder$data <- rev(outFolder$data[out$w2[1]:out$w2[2]])
  
  return(outFolder)  
}


## Internal function rsd2D
## Reads 2D RSD format spectra and returns all of the tiles covered by the 
##	chemical shift range provided.
## file.name - name of spectrum to be read, NULL opens a file slection window
## w1Range - Min and max chemical shifts to be read in the indirect dimension,
##           NULL will return all values.
## w2Range - Min and max chemical shifts to be read in the direct dimension
##           NULL will return all values.
## file.par - file header (from rsdHead()), NULL will read the header first
## notes: providing a file.par from memory is used to speed up graphics
## returns the designated region of a spectrum and associated parameters
rsd2D <- function(file.name=NULL, w1Range=NULL, w2Range=NULL,
                  file.par=NULL){
  
  ## Get the path to the file and RSD header
  if (is.null(file.name))
    file.name <- myOpen(title='Select a spectrum', multiple=FALSE, 
                        filetypes=list(rsd="RSD File"))
  if (!length(file.name) || !nzchar(file.name))
    return(invisible())
  if (is.null(file.par))
    file.par <- rsdHead(file.name=file.name, print.info=FALSE)[[1]]
  
  ## Setup outgoing list 
  outFolder <- list(file.par=file.par, w1=seq(file.par$upfield_ppm[1], 
                                              file.par$downfield_ppm[1], length.out=file.par$matrix_size[1]),
                    w2=seq(file.par$upfield_ppm[2], file.par$downfield_ppm[2], 
                           length.out=file.par$matrix_size[2]), data=NULL) 
  
  ## Make sure valid w2 ranges were submitted
  if (is.null(w1Range) || length(w1Range) != 2 || !is.numeric(w1Range))
    w1Range <- c(file.par$upfield_ppm[1], file.par$downfield_ppm[1])
  if (is.null(w2Range) || length(w2Range) != 2 || !is.numeric(w2Range))
    w2Range <- c(file.par$upfield_ppm[2], file.par$downfield_ppm[2])
  
  ## Find best w1/w2 matches	
  w1Range <- sort(w1Range)
  w2Range <- sort(w2Range)
  t1 <- findInterval(w1Range, outFolder$w1, all.inside=TRUE)
  d1 <- findInterval(w2Range, outFolder$w2, all.inside=TRUE)
  t2 <- t1 + 1
  d2 <- d1 + 1
  out <- list(w1=NULL, w2=NULL)
  for (i in 1:2){
    out$w1[i] <- switch(which.min(c(abs(w1Range[i] - outFolder$w1[t1[i]]),
                                    abs(w1Range[i] - outFolder$w1[t2[i]]))), t1[i], t2[i])
    out$w2[i] <- switch(which.min(c(abs(w2Range[i] - outFolder$w2[d1[i]]),
                                    abs(w2Range[i] - outFolder$w2[d2[i]]))), d1[i], d2[i])
  }
  winW1 <- outFolder$w1[(out$w1[1]:out$w1[2])]
  winW2 <- outFolder$w2[(out$w2[1]:out$w2[2])]
  
  ## Find tiles that will be plotted
  tiles <- NULL
  upShifts <- file.par$block_upfield_ppms
  downShifts <- file.par$block_downfield_ppms
  blockSizes <- file.par$block_size
  for (tNum in seq_along(file.par$block_size$w1)){
    
    ## Get the chemical shifts for the current tile
    blockW1 <- seq(upShifts$w1[tNum], downShifts$w1[tNum], 
                   length.out=blockSizes$w1[tNum])
    blockW2 <- seq(upShifts$w2[tNum], downShifts$w2[tNum],
                   length.out=blockSizes$w2[tNum])
    
    ## Check the window for the presence of any shift in the current block
    if (any(round(blockW1, 3) %in% round(winW1, 3)) && 
        any(round(blockW2, 3) %in% round(winW2, 3)))
      tiles <- c(tiles, tNum)
  }
  
  ## Define data locations for each block
  blockLoc <- file.par$binary_location
  for (i in seq_along(file.par$block_size$w1))
    blockLoc <- c(blockLoc, blockLoc[i] + 4 * file.par$block_size$w1[i] * 
                    file.par$block_size$w2[i])
  
  ## Read binary data for each tile
  outData <- matrix(0, nrow=length(outFolder$w2), ncol=length(outFolder$w1))
  con <- myFile(file.name, "rb")
  for (i in tiles){
    
    ## Find best w1/w2 matches for current tile
    w1Range <- c(upShifts$w1[i], downShifts$w1[i])
    w2Range <- c(upShifts$w2[i], downShifts$w2[i])
    t1 <- findInterval(w1Range, outFolder$w1, all.inside=TRUE)
    d1 <- findInterval(w2Range, outFolder$w2, all.inside=TRUE)
    t2 <- t1 + 1
    d2 <- d1 + 1
    tile <- NULL
    for (j in 1:2){
      tile$w1[j] <- switch(which.min(c(abs(w1Range[j] - outFolder$w1[t1[j]]),
                                       abs(w1Range[j] - outFolder$w1[t2[j]]))), t1[j], t2[j])
      tile$w2[j] <- switch(which.min(c(abs(w2Range[j] - outFolder$w2[d1[j]]),
                                       abs(w2Range[j] - outFolder$w2[d2[j]]))), d1[j], d2[j])
    }
    
    ## Read data for current tile and place in matrix
    seek(con, blockLoc[i], origin='start')
    outData[tile$w2[1]:tile$w2[2], tile$w1[1]:tile$w1[2]] <- 
      matrix(rev(readBin(con, size=4, what='double', endian='big', 
                         n=file.par$block_size$w1[i] * file.par$block_size$w2[i])), 
             ncol=file.par$block_size$w1[i])
  }
  close(con)
  
  ## Trim data to fit w1/w2 ranges
  outFolder$w1 <- winW1
  outFolder$w2 <- winW2
  outData <- outData[out$w2[1]:out$w2[2], out$w1[1]:out$w1[2]]
  outFolder$data <- matrix(rev(outData), ncol=length(outFolder$w1))
  
  return(outFolder)
}


## Internal function ucsfTile
## Reads a single sparky tile from a binary connection and returns a data folder
## file.par - digestR file header for file to be read (output from ucsfHead)
## con      - A seekable connection, if missing, a connection will be opened
##            to the path provided by file.par$file.name 
## tile     - The sparky tile to be read as a zero indexed integer
## w1       - If tile is missing, then w1 and w2 chemical shifts can be given
##            and the entire tile that corresponds to the shifts is returned
## w2       - Numeric value of w2 chemical shift in ppm
## returns  - A data folder with the file header, data, and chemical shifts
##            for the entire tile.
## Note:    - This function is used for plotting 2D data and 2D peak picking
##            and is used to avoid loading large chunks of data into memory.
ucsfTile <- function(file.par, con, tile, w1, w2){
  
  ## Check incoming arguments
  if( missing(tile) && missing(w1) && missing(w2) )
    stop('A tile or w1Range and w2Range must be provided')
  if( missing(file.par) )
    stop('A file.par must be provided')
  if( file.par$number_dimensions != 2 )
    stop('ucsfTile only supports 2D sparky data')
  if( missing(con) || !isSeekable(con) ){
    con <- myFile(file.par$file.name, "rb")
    closeCon <- TRUE
  }else
    closeCon <- FALSE
  
  ## Setup outgoing list 
  outFolder <- list( file.par = file.par, w1 = NULL, w2 = NULL, data = NULL )
  tTiles <- ceiling(file.par$matrix_size / file.par$block_size )
  
  if( missing(tile)){
    
    ## Make sure valid w2 ranges were submitted
    if( any( !is.numeric(c(w1, w2)) ) || length(c(w1, w2)) != 2 )
      stop('w1 and w2 must be numeric values of length 1, use ucsf2D')
    
    ## Find w1/w2 tiles and sparky tiles
    mShift <- matchShift(outFolder, w1 = w1, w2 = w2, 
                         return.inc = TRUE, return.seq = FALSE, invert = TRUE)		
    w1Tiles <- (mShift$w1 - 1)  %/% file.par$block_size[1]
    w2Tiles <- (mShift$w2 - 1)  %/% file.par$block_size[2]
    tile <- w1Tiles * tTiles[2] + w2Tiles[1]
    
  }else{
    ## Stop if more than one tile is provided
    if(length(tile) > 1 )
      stop('ucsfTile can only read 1 tile at a time, use ucsf2D')
    w1Tiles <- floor( tile %/% tTiles[2] )
    w2Tiles <- floor( tile %% tTiles[2] ) 	
  }
  
  ## Set the outgoing chemical shift range
  w1Range <- (1:file.par$block_size[1]) + (w1Tiles * file.par$block_size[1])
  w1Range <- w1Range[w1Range <= file.par$matrix_size[1]]
  w2Range <- (1:file.par$block_size[2]) + (w2Tiles * file.par$block_size[2])
  w2Range <- w2Range[w2Range <= file.par$matrix_size[2]]
  outFolder$w1 <- rev(seq(file.par$downfield_ppm[1], file.par$upfield_ppm[1],
                          length.out = file.par$matrix_size[1])[w1Range])	
  outFolder$w2 <- rev(seq(file.par$downfield_ppm[2],file.par$upfield_ppm[2],
                          length.out = file.par$matrix_size[2])[w2Range])
  
  ## Find the binary location an read the data
  seek(con, file.par$block_size[1] * file.par$block_size[2] * 4 * tile + 
         file.par$binary_location, origin = 'start')
  outFolder$data <- matrix(readBin(con, size=4, what='double', 
                                   endian = file.par$endian, 
                                   n=(file.par$block_size[1] * file.par$block_size[2])), 
                           ncol = file.par$block_size[1])
  outFolder$data <- outFolder$data[(1:length(w2Range)), (1:length(w1Range))]
  
  if(closeCon)
    closeAllConnections()	
  
  return(outFolder)
}


################################################################################
##                                                                            ##
##                    Chemical shift utilities                                ##
##                                                                            ##
################################################################################

## Internal function matchShift
## Finds the closest match for the provided chemical shifts in a spectrum.
##       If shifts are outside of the spectrum's range, the closest shift 
##       from the spectrum will be returned
## inFolder  - Any entry from the fileFolder, default is the current spectrum
## w1        - Numeric argument or vector of w1 chemical shifts to be matched
## w2        - Numeric argument or vector of w2 chemical shifts to be matched
## Note      - If a vector of shifts is provided, then this function returns 
##             the shifts matched from the range of the vector
## w1.pad    - Integer, number of padding points over the w1 range to return
## w2.pad    - Integer, number of padding points over the w2 range to return
## Note      - Padding is used in peak picking, it allows peaks occurring on the 
##             edges of the spectral window to be peak picked
## invert    - Logical argument, will return shifts from the far edge of the 
##             spectrum rather than the closest match. i.e. the shifts as 
##             reflected across the center of the spectrum. 
## return.inc- Logical argument, TRUE returns the file index (point), FALSE
##             (the default) returns the chemical shift
## return.seq- Logical argument, TRUE returns the entire sequence of shifts
##             covered by the range of w1 and w2 provided, FALSE returns 
##             the range of the matched shifts (FALSE is the default)
## overRange - Logical argument, TRUE returns NA if shifts are outside spectra
##             window; FALSE returns closest match to the shifts provided 
## returns the range of chemical shift, or indices , of the closest w1 and w2 
##             shifts found in a spectrum. 
matchShift <- function (inFolder = fileFolder[[wc()]], w1 = NULL, w2 = NULL, 
                        w1.pad = 0, w2.pad = 0, invert = FALSE, return.inc = FALSE, 
                        return.seq = FALSE, overRange = FALSE ){
  
  if(is.null(w1) && is.null(w2))
    stop('No shifts were entered')
  
  
  ## Mantain code format for 1D data
  if(inFolder$file.par$number_dimensions == 1){
    w1 <- NULL
    for( i in c('upfield_ppm', 'downfield_ppm', 'matrix_size'))
      inFolder$file.par[[i]] = c(inFolder$file.par[[i]],
                                 inFolder$file.par[[i]])
  }
  
  ## Find best chemical shift match for w1
  out <- list(w1 = NULL, w2 = NULL)
  if( !is.null(w1) ){
    inFolder$w1 <- seq(inFolder$file.par$upfield_ppm[1],
                       inFolder$file.par$downfield_ppm[1],
                       length.out = inFolder$file.par$matrix_size[1])
    w1Range <- range(w1)	
    
    ## Check for values outside spectral window
    if( overRange && 
        ( all(w1Range > max(inFolder$w1)) || all(w1Range < min(inFolder$w1))))
      naOut <- TRUE
    else
      naOut <- FALSE
    
    ## Find shifts in data that are the closest match to input data
    t1 <- findInterval(w1Range, inFolder$w1, all.inside = TRUE)	
    t2 <- t1 + 1
    for(i in 1:2)
      out$w1[i] <- switch( which.min(c(
        abs(w1Range[i] - inFolder$w1[t1[i]]),
        abs(w1Range[i] - inFolder$w1[t2[i]]))), t1[i], t2[i])
    
    ## Pad outgoing data
    out$w1 <- out$w1 + c(-w1.pad, w1.pad)
    out$w1[out$w1 < 1] <- 1
    out$w1[out$w1 > inFolder$file.par$matrix_size[1]] <- 
      inFolder$file.par$matrix_size[1]
    
    
    ## Invert, convert to sequence, and translate to shifts if requested
    if( invert )
      out$w1 <- sort(inFolder$file.par$matrix_size[1] - out$w1) + 1
    if(return.seq)
      out$w1 <- out$w1[1]:out$w1[2]
    if(!return.inc)
      out$w1 <- inFolder$w1[out$w1]
    if(naOut)
      out$w1 <- NA
  }
  
  ## Find best chemical shift match for w2
  if( !is.null(w2) ){
    inFolder$w2 <- seq(inFolder$file.par$upfield_ppm[2],
                       inFolder$file.par$downfield_ppm[2],
                       length.out = inFolder$file.par$matrix_size[2])
    w2Range <- range(w2)
    
    ## Check for values outside spectral window
    if( overRange && 
        ( all(w2Range > max(inFolder$w2)) || all(w2Range < min(inFolder$w2))))
      naOut <- TRUE
    else
      naOut <- FALSE
    
    ## Find shifts in data that are the closest match to input data
    d1 <- findInterval(w2Range, inFolder$w2, all.inside = TRUE)
    d2 <- d1 + 1
    for(i in 1:2)	
      out$w2[i] <- switch( which.min(c(
        abs(w2Range[i] - inFolder$w2[d1[i]]),
        abs(w2Range[i] - inFolder$w2[d2[i]]))), d1[i], d2[i])
    
    ## Pad outgoing data
    out$w2 <- out$w2 + c(-w2.pad, w2.pad)
    out$w2[ out$w2 < 1] <- 1
    out$w2[ out$w2 > inFolder$file.par$matrix_size[2]] <- 
      inFolder$file.par$matrix_size[2]
    
    
    ## Invert, convert to sequence, and translate to shifts if requested
    if( invert )
      out$w2 <- sort(inFolder$file.par$matrix_size[2] - out$w2) + 1
    if(return.seq)
      out$w2 <- out$w2[1]:out$w2[2]
    if(!return.inc)
      out$w2 <- inFolder$w2[out$w2]
    if(naOut)
      out$w2 <- NA
  }
  
  return(out)       
}

## Internal chemical shift table utility function shiftToROI
## Converts a peak list to an roi table
## shiftList - chemical shift file containing an w1 and a w2 column
##           a column named Code or Assignment will be used for the 
##           ROI names if present.
## w1Delta - Width of the desired ROI in the w1 dimension
## w2Delta - Width of the ROI in the w2 dimension
## returns a table in the same format as the roiTable format
shiftToROI <- function(shiftList = NULL, w1Delta = 1, w2Delta = .05 ){
  
  ## Check the input format
  if( is.null(shiftList) )
    stop( 'No chemical shift table was entered' )
  if( length(which(names(shiftList) == 'w2')) == 0 ){
    if( length(which(names(shiftList) == 'Height')) == 0 )
      stop( '1D shift tables must have columns labeled w2 and Height' )
    else
      stop( '1D shift tables must have an w2 column' )
  }
  if( !is.numeric(w1Delta) || !is.numeric(w2Delta) )
    stop( 'w1Delta and w2Delta must be numeric')
  w1Delta <- w1Delta / 2
  w2Delta <- w2Delta / 2
  
  ## Generate a 1D roi table
  if( length(shiftList$w1) == 0 ){
    
    outTable <- list( Name = NULL, 
                      w2_downfield = shiftList$w2 + w2Delta, 
                      w2_upfield = shiftList$w2 - w2Delta, 
                      w1_downfield = rep(0, length(shiftList$w2)), 
                      w1_upfield = shiftList$Height + 0.1*shiftList$Height, 
                      ACTIVE = rep(TRUE, length(shiftList$w2)), 
                      nDim = rep(1, length(shiftList$w2)) )		
  }else{
    
    ## Generate a 2D roi table
    outTable <- list( Name = NULL, 
                      w2_downfield = shiftList$w2 + w2Delta, 
                      w2_upfield = shiftList$w2 - w2Delta, 
                      w1_downfield = shiftList$w1 + w1Delta, 
                      w1_upfield = shiftList$w1 - w1Delta, 
                      ACTIVE = rep(TRUE, length(shiftList$w2)), 
                      nDim = rep(2, length(shiftList$w2)) )		
    
  }
  
  ## Add assignments if possible
  if( !is.null(shiftList$Code) || !is.null(shiftList$Assignment) ){
    if(!is.null(shiftList$Assignment) )
      outTable$Name <- shiftList$Assignment
    else
      outTable$Name <- shiftList$Code	
    
    ## Replace NA or "NA" with ROI names ("NA" is returned from digestR peak lists)
    outTable$Name[ which(outTable$Name == 'NA') ] <- NA
    if( any(is.na(outTable$Name)) )
      outTable$Name[ which( is.na(outTable$Name)) ] <- 
      paste('ROI', 1:length(which(is.na(outTable$Name))), sep='.')	
    
    ## Make sure ROI names are unique
    for( i in 1: length(outTable$Name) ){
      tName <- which(outTable$Name == outTable$Name[i])
      if( length (tName) != 1 )
        outTable$Name[tName] <- paste(outTable$Name[i], 
                                      1:length(tName), sep ='.' )
    }
    
  }else
    outTable$Name <- paste('ROI', 1:length(shiftList$w2), sep='.')
  
  return(data.frame( outTable, stringsAsFactors = FALSE ))      
}

## Internal function autoRef
## fileName  - File to be referenced, defaults to currentSpectrum
## w2Range  - Chemical shift range, in ppm, to be searched,
##            NULL searches the entire spectrum
## w2Shift  - The reference frequency in ppm
## Note: This function is is used to reference well phased spectra containing
##       DSS, TSP, etc. The function finds furthest up field peak above the 
##       1D threshold and defines it as w2Shift (0 by default). 
## Note: Currently autoRef only handles 1D data
## Note: The w2Range option allows some flexibility over which portion of the 
##       spectrum is searched. However, if the referencing is really off, then
##       setting the w2Range argument may lead to some undesirable behavior.
## Note: The chemical shift referencing in digestR does not alter the original data.
## Returns a referenced spectrum/spectra and refreshes the active plots
autoRef <- function(fileNames = currentSpectrum, w2Range=NULL, 
                    w2Shift = 0, all1D = FALSE ){   
  
  oneDs <- which(sapply(fileFolder[fileNames], function(x) 
    x$file.par$number_dimensions) == 1)
  if (length(oneDs) == 0)
    err('No 1D files are open, digestR can only auto reference 1D files')
  specList <- fileNames[oneDs]
  if (length(oneDs) != length(specList))
    cat('\n', 'The following spectra could not be referenced because they ', 
        'are not one-dimensional:\n  ', 
        paste(getTitles(specList[-oneDs], FALSE), '\n  '), '\n')
  
  ## Update chemical shifts of files
  for( i in  specList ){
    
    ## Load spectral data from binary 
    in.folder <- ucsf1D( file.name = fileFolder[[i]]$file.par$file.name,
                         file.par=fileFolder[[i]]$file.par, w2Range = w2Range)
    in.folder$data <- rev( in.folder$data ) ## This is needed to match binary
    
    ## Find Furthest downfield signal above the peak picking threshold
    thresh <- globalSettings$thresh.1D * 
      in.folder$file.par$noise_est + in.folder$file.par$zero_offset
    w2DSS <- min ( in.folder$w2[ in.folder$data >= thresh ] )
    
    ## Define the possible chemical shift range for DSS
    w2DSS <- c(w2DSS + 0.1, w2DSS)
    w2DSS <- which(in.folder$w2 >= w2DSS[2] & in.folder$w2 <= w2DSS[1])
    
    ## Define DSS as the biggest peak in the upfield sub region
    w2DSS <- in.folder$w2[ w2DSS ][which.max(in.folder$data[ w2DSS ])]
    
    if( length(w2DSS) == 0 ){
      log_message( paste( 'Auto referencing failed in:',
                    basename(in.folder$file.par$file.name)), quote = FALSE )
      next()
    }
    
    ## Move reference to specified location
    fileFolder[[i]]$file.par$downfield_ppm[1] <- 
      fileFolder[[i]]$file.par$downfield_ppm[1] - w2DSS + w2Shift
    fileFolder[[i]]$file.par$upfield_ppm[1] <- 
      fileFolder[[i]]$file.par$upfield_ppm[1] - w2DSS + w2Shift
    
    ## Update peaklist and graphics
    fileFolder[[i]]$graphics.par$usr[1:2] <- 
      fileFolder[[i]]$graphics.par$usr[1:2] - w2DSS + w2Shift
    if(length(fileFolder[[i]]$peak.list) > 0)
      fileFolder[[i]]$peak.list$w2 <- (
        fileFolder[[i]]$peak.list$w2 - w2DSS + w2Shift)
    
  }
  
  ## Assign changes to the file folder and refresh graphics 
  myAssign("fileFolder", fileFolder)
  refresh() 
  lineCol <- fileFolder[[wc()]]$graphics.par$fg
  abline ( v = w2Shift, lty = 2, col=lineCol )
  
}

## Internal peak picking function isNoise
## A simple filter for detecting noise signals in 1D data 
## x    - Numeric argument, a candidate signal
## data - Numeric vector, the field of data being tested
## thresh  - Numeric argument expressing the threshold. Values range from
##        0 (no filtering) to -1 no data returned. The default, .15 seems 
##        like a reasonable compromise between false discovery and sensitivity
## Note: This function is an obvious area for future improvement, the filter
##       currently excludes broad signals.
isNoise <- function( x, data, thresh = -.15 ){
  return( (min((data - x)/ x, na.rm = TRUE )) > thresh  )
} 

## Internal function for combining two peak lists
appendPeak <- function(newList, oldList){
  
  if( is.null(oldList) || length(oldList) == 0)
    return(newList)
  if( is.null(newList) || length(newList) == 0)
    return(oldList)
  
  ## Allow lists to differ in structure
  newName <- unique(c(names(newList), names(oldList)))
  
  ## Add/update the index column to the new list
  if (is.null(newList$Index))
    newList$Index <- 1:nrow(newList)
  if (is.null(oldList$Index))
    oldList$Index <- 1:nrow(oldList)
  newList$Index <- 1:nrow(newList) + max(oldList$Index)
  
  cList <- data.frame(matrix(rep(NA, length(newName) * 2), nrow=2))
  names(cList) <- newName
  
  noMatch <- which(is.na(match(names(cList), names(oldList))))
  if (length(noMatch)){
    oNames <- names(oldList)
    oldList <- cbind(oldList, matrix(rep(NA, 
                                         nrow(oldList) * length(noMatch)), ncol=length(noMatch)))
    names(oldList) <- c(oNames, names(cList)[noMatch])
  }
  
  noMatch <- which(is.na(match(names(cList), names(newList))))
  if (length(noMatch)){
    oNames <- names(newList)
    newList <- cbind(newList, matrix(rep(NA, 
                                         nrow(newList) * length(noMatch)), ncol=length(noMatch)))
    names(newList) <- c(oNames, names(cList)[noMatch])
  }
  
  ## Remove duplicate entries
  cList <- rbind(cList[-(1:2),], oldList, newList)
  checkVar <- match(c('w1', 'w2', 'Height'), names(cList))
  checkVar <- which(duplicated(cList[,checkVar]))
  if (length(checkVar))
    cList <- cList[-checkVar, ]
  
  row.names(cList) <- NULL
  cList$Index <- 1:nrow(cList)
  return(cList)
}

## Internal wrapper function for implementing 1D and 2D peak picking
## fileName    - name of the file to be peak picked, NULL will pick the current
## inFile			 - file data as returned by ucsf1D or ucsf2D, used instead of the
##							 fileName argument to peak pick a file without opening it
## w1Range     - w1 chemical shift range c(downfield,upfield) to be used
## w2Range     - w2 chemical shift range c(downfield,upfield) to be used
## append - logical argument, TRUE appends peaks to old list
## internal - logical argument, TRUE returns list without assigning it to
##            fileFolder, FALSE assigns the file to fileFolder
## ... - Additional peak picking parameters passed to trans2peak
## Returns the new peak list 
peakPick <- function( fileName = currentSpectrum, inFile = NULL, 
                      w1Range = NULL, w2Range = NULL, append = FALSE, internal = FALSE, ...){
  
  if( length(fileName) > 1 && internal )
    stop('Only one peak list can be returned internally')
  
  if (!is.null(inFile)){
    if(inFile$file.par$number_dimensions == 1)
      nList <- peakPick1D( inFolder = inFile, w2Range = w2Range, ...)
    else if(!is.null(inFile$file.par$block_upfield_ppms))
      nList <- peakPickRsd( inFolder = inFile, w1Range = w1Range, 
                            w2Range = w2Range, ...)
    else if(inFile$file.par$number_dimensions == 3)
      nList <- peakPick3D( inFolder = inFile, w1Range = w1Range, 
                           w2Range = w2Range, ...)
    else
      nList <- peakPick2D( inFolder = inFile, w1Range = w1Range, 
                           w2Range = w2Range, ...)
    row.names(nList) <- NULL
    return(nList)
  }
  
  for(i in fileName){			
    if(fileFolder[[i]]$file.par$number_dimensions == 1)
      nList <- peakPick1D( fileName = i, w2Range = w2Range, ...)
    else if(!is.null(fileFolder[[i]]$file.par$block_upfield_ppms))
      nList <- peakPickRsd( fileName = i, w1Range = w1Range, 
                            w2Range = w2Range, ...)
    else if(fileFolder[[i]]$file.par$number_dimensions == 3)	
      nList <- peakPick3D( fileName = i, w1Range = w1Range, 
                           w2Range = w2Range, ...)
    else
      nList <- peakPick2D( fileName = i, w1Range = w1Range, 
                           w2Range = w2Range, ...)
    row.names(nList) <- NULL
    if( internal )
      break()
    
    ## Create peak labels if not provided
    if (!is.null(nList) && (all(is.na(nList$Assignment)) || 
                            length(which(nList$Assignment == 'NA')) == length(nList$Assignment)))
      nList$Assignment <- paste('P', nList$Index, sep='')
    
    ## Append unique new peaks to old list
    oList <- fileFolder[[i]]$peak.list
    if( append )
      nList <- appendPeak(nList, oList)
    fileFolder[[i]]$peak.list <- nList
    
    if ( !internal ){
      cat(paste(fileFolder[[i]]$file.par$user_title, ':', '\n Total peaks: ', 
                nrow(nList), '\n', sep = ''))
      if( append ){
        if(is.null(oList))
          cat(paste(' New peaks:', nrow(nList), '\n'))
        else
          cat(paste(' New peaks:', nrow(nList) - nrow(oList), '\n'))
      }
      
      flush.console()
    }	
  }
  
  if( !internal )
    myAssign('fileFolder', fileFolder)
  
  return(nList)
}

## Internal peak picking function peakPick1D
## General local maximum (hill climbing method) for 1D peak picking
## fileName    - name of the file to be peak picked, NULL will pick the current
## inFile			 - file data as returned by ucsf1D or ucsf2D, used instead of the
##							 fileName argument to peak pick a file without opening it
## w2Range     - w2 chemical shift range c(downfield,upfield) to be used
## w2Gran      - integer controlling granularity of search space, 
##         			 smaller values are more exhaustive bigger values supress noise
## noiseFilt - Integer argument that can be set to 0, 1 or 2; 
##              0 does not apply a noise filter, 1 applies a mild filter
##              (adjacent points in the direct dimension must be above the 
##              noise threshold), 2 applies a strong filter (all adjacent points
##              must be above the noise threshold
## maxOnly - logical, if TRUE only the maximum peak is returned
## Note: currently the filter excludes broad peaks. 
## Note: the threshold used for peak picking taken from the global parameters
## ... - Additional peak picking parameters
## Returns a new peak list 
peakPick1D <- function( fileName = currentSpectrum, inFile = NULL, 
                        w2Range = NULL,	w2Gran = 2, noiseFilt = globalSettings$peak.noiseFilt, 
                        maxOnly = FALSE,  ... ){
  
  ## Define the current spectrum and w2 range 
  if (is.null(inFile)){
    inFile <- fileFolder[[fileName]]
    inFile <- ucsf1D(file.name = fileName, file.par=inFile$file.par)
  }
  if( inFile$file.par$number_dimensions != 1 )
    stop('You must use peakPick1D to peak pick 1D files')
  if(is.null(w2Range) || length(w2Range) != 2 )
    w2Range <- c(inFile$file.par$upfield_ppm[1],
                 inFile$file.par$downfield_ppm[1])
  
  ## Find all local maxes
  absData <- abs(inFile$data)
  allPos <- intersect(which(c(NA, absData) <  c(absData, NA)), 
                      which(c(NA, absData) >  c(absData, NA))-1)
  
  ## Remove local maxes below the threshold
  thresh <- globalSettings$thresh.1D * 
    inFile$file.par$noise_est + inFile$file.par$zero_offset 
  allPos <- allPos[which(absData[allPos] >= thresh )]
  
  ## Apply granularity and noise filters
  out <- NULL
  for(i in 1:length(allPos) ){
    subRan <- try(absData[(allPos[i] - w2Gran):(allPos[i] + w2Gran)], 
                  silent = TRUE)
    if(class(subRan) == "try-error")
      next()
    x <- subRan[ w2Gran + 1 ]
    if( max( subRan, na.rm = TRUE) == x  ){
      if( noiseFilt == 0  )
        out <- c(out, i)
      else{
        noise <- isNoise( x = x, data = subRan, ...)
        if( !noise )
          out <- c(out, i)
      }
    }
  }
  allPos <- allPos[out]
  
  ## Build a peak list
  mShift <- matchShift(inFolder = inFile, w2 = w2Range, return.inc = TRUE, 
                       invert = TRUE)$w2
  allPos <- allPos[ allPos <= mShift[2] & allPos >= mShift[1] ]
  if(length(allPos) == 0)
    return(NULL)
  n1 <- length(allPos)
  peak <- data.frame( list( Index = 1:n1, w1 = rep(NA, n1), 
                            w2 = rev(inFile$w2)[allPos], Height = inFile$data[allPos],
                            Assignment = rep('NA', n1)), stringsAsFactors = FALSE )
  if( maxOnly ){
    peak <- peak[which.max(abs(peak$Height)),]
    peak$Index <- 1
  }
  
  return(peak)
}

## Internal peak picking function peakPick2D
## General local maximum (hill climbing method) for 2D peak picking
## fileName    - name of the file to be peak picked, NULL will pick the current
## inFile			 - file data as returned by ucsf1D or ucsf2D, used instead of the
##							 fileName argument to peak pick a file without opening it
## w1Range     - w1 chemical shift range c(downfield,upfield) to be used
## w2Range     - w2 chemical shift range c(downfield,upfield) to be used
## fancy       - logical argument, FALSE implements a basic peak picker
##               that returns local maxima only, this is fastest; TRUE
##               determines chemical shifts of peaks, groups multiplets,
##               and measures line width and volume
## noiseFilt   - Integer argument that can be set to 0, 1 or 2; 
##               0 does not apply a noise filter, 1 applies a mild filter
##               (adjacent points in the direct dimension must be above the 
##               noise threshold), 2 applies a strong filter (all adjacent points
##               must be above the noise threshold
## maxOnly - logical, if TRUE only the absolute maximum peak is returned
## Note: the threshold used for peak picking taken from the graphics parameters
## ... - Additional peak picking parameters passed to trans2peak
## Returns a new peak list for the w1/w2 ranges provided
peakPick2D <- function(fileName = currentSpectrum, inFile = NULL, 
                       w1Range = NULL, w2Range = NULL, fancy = FALSE, 
                       noiseFilt = globalSettings$peak.noiseFilt, maxOnly = FALSE, ...){
  
  ## Error checking
  if (is.null(inFile))
    inFile <- fileFolder[[fileName]]	
  if(is.null(w1Range))
    w1Range <- c(inFile$file.par$upfield_ppm[1],
                 inFile$file.par$downfield_ppm[1])
  if(is.null(w2Range))
    w2Range <- c(inFile$file.par$upfield_ppm[2],
                 inFile$file.par$downfield_ppm[2])
  if(is.null(inFile$graphics.par$tiles))
    inFile <- findTiles(in.folder = inFile, internal=TRUE)
  
  ## Find best chemical shift match
  bs <- inFile$file.par$block_size
  w1Range <- sort(w1Range); w2Range <- sort(w2Range)
  mShift <- matchShift(inFile, w1 = w1Range, w2 = w2Range, 
                       return.inc = TRUE, return.seq = TRUE, invert = TRUE)
  
  ## Find the total possible sparky tiles for the range
  w1Tiles <- unique((mShift$w1 - 1)  %/% bs[1])
  w2Tiles <- unique((mShift$w2 - 1)  %/% bs[2])
  tTiles <- ceiling(inFile$file.par$matrix_size / bs )
  tiles <- NULL
  for( i in 1:length(w2Tiles) )
    tiles <- c(tiles, w1Tiles * tTiles[2] + w2Tiles[i])
  tiles <- sort(tiles)
  
  
  ## Set up for tile overlap
  w1Slice <- rep(NA, tTiles[1] * bs[1])
  w1Slice <- rbind(w1Slice, w1Slice, w1Slice, w1Slice, w1Slice)
  w2Slice <- rep(NA, tTiles[2] * bs[2] + 5)
  w2Slice <- cbind(w2Slice, w2Slice, w2Slice, w2Slice, w2Slice)
  
  
  ## Open connection to binary file
  con <- myFile(inFile$file.par$file.name, "rb")
  binLoc <- inFile$file.par$binary_location
  endForm <- inFile$file.par$endian
  conDisp <- inFile$graphics.par$conDisp
  thresh <- inFile$file.par$noise_est * inFile$graphics.par$clevel
  locMax <- allPos <- NULL
  for(i in 1:length(tiles)){
    
    ## Define the current w1/w2 increments
    w1B <- 1:bs[1] + bs[1] * floor( tiles[i] %/% tTiles[2] ) 
    w2B <- 1:bs[2] + bs[2] * floor( tiles[i] %% tTiles[2] )	
    
    ##skip tiles with no data above the threshold and set w1/w2 slices to NA
    if( length(which(inFile$graphics.par$tiles == (tiles[i]))) == 0 ){	
      w1Slice[ 1:5, w1B ] <- NA
      w2Slice[ w2B, 1:5 ] <- NA
      next() 			
    }
    
    ## Read binary
    seek(con, bs[1] * bs[2] * 4 * tiles[i] + binLoc, origin = 'start')
    tData <- matrix(rev(readBin(con, size=4, what='double', 
                                endian = endForm, n=(bs[1] * bs[2]))), ncol = bs[1])
    
    ## Make a list of all observable signals for fancy peak picking
    if(fancy){
      tObs <- obs2List(tData, conDisp = conDisp, thresh = thresh, bs = bs)
      tObs$tRow <- tObs$tRow + bs[2] * 
        floor( (tTiles[1] * tTiles[2] - tiles[i] -1) %% tTiles[2])
      tObs$tCol <- tObs$tCol + bs[1] * 
        floor( (tTiles[1] * tTiles[2] - tiles[i] -1) %/% tTiles[2])	
      allPos <- rbind(allPos, tObs)
    }
    
    ## Overlap tile, and save upfield slices for next tile
    w1Lap <- w1Slice[1:5, w1B]
    pad <- rev(w2B)[1] + c(1,2,3,4,5)
    w2Lap <- w2Slice[c(w2B, pad), 1:5]
    w1Slice[1:5, w1B] <- tData[1:5,] 
    w2Slice[w2B, 1:5] <- tData[,1:5]
    tData <- cbind(rbind(tData, w1Lap), w2Lap)
    
    ## Find local maxes
    if( conDisp[1] && conDisp[2] )
      oThresh <- localMax(abs(tData), thresh=thresh, noiseFilt=noiseFilt ) - 1
    else{
      if( conDisp[1])
        oThresh <- localMax(tData, thresh=thresh, noiseFilt=noiseFilt ) - 1
      else
        oThresh <- localMax(-tData, thresh=thresh, noiseFilt=noiseFilt ) - 1			
    }
    
    peak <-data.frame(list( 
      tRow = (oThresh %% (bs[2] +5)) + 1, 
      tCol = (oThresh %/% (bs[2] +5)) + 1, 
      Height = tData[oThresh + 1]))
    peak <- peak[ peak[,1] > 1 & peak[,1] < bs[2]+2, ]		
    peak <- peak[peak[,2] > 1 & peak[,2] < bs[1] + 2, ]
    if (!length(peak) || nrow(peak) == 0)
      next()
    
    ## Translate tile row/column to spectrum row/column
    peak$tRow <- peak$tRow + bs[2] * 
      floor( (tTiles[1] * tTiles[2] - tiles[i] -1) %% tTiles[2])
    peak$tCol <- peak$tCol + bs[1] * 
      floor( (tTiles[1] * tTiles[2] - tiles[i] -1) %/% tTiles[2])
    locMax <- rbind(locMax, peak)
    
  }
  
  ## Close binary conection
  closeAllConnections()
  
  ## Clean up peak list
  if(is.null(locMax))
    return(NULL)
  locMax <- unique(locMax)
  
  ## Match point locations with chemical shifts
  w1 <- seq(inFile$file.par$upfield_ppm[1], 
            inFile$file.par$downfield_ppm[1],
            length.out = inFile$file.par$matrix_size[1])
  if( inFile$file.par$matrix_size[1] %% bs[1] != 0 )
    w1 <- c(rep(NA, bs[1] - inFile$file.par$matrix_size[1] %% bs[1]), w1)
  w2 <- seq(inFile$file.par$upfield_ppm[2], 
            inFile$file.par$downfield_ppm[2],
            length.out = inFile$file.par$matrix_size[2])
  if( inFile$file.par$matrix_size[2] %% bs[2] != 0 )
    w2 <- c(rep(NA, bs[2] - inFile$file.par$matrix_size[2] %% bs[2] ), w2)
  w1Range <- sort(rev(w1)[range(mShift$w1)])
  w2Range <- sort(rev(w2)[range(mShift$w2)])
  
  if(fancy){
    peak <- trans2Peak( allPos = allPos, locMax = locMax, w1 = w1, w2 = w2, 
                        w1Range = w1Range, w2Range = w2Range, ...)
  }else{
    
    peak <- data.frame( list(	
      Index = 1:nrow(locMax),
      w1 = w1[locMax$tCol],
      w2 = w2[locMax$tRow],
      Height = locMax$Height,
      Assignment = rep("NA", length(locMax$tCol)) ), 
      stringsAsFactors = FALSE)
    
    ## Filter the outgoing list to match the data range
    peak <- peak[peak$w1 <= w1Range[2] & peak$w1 >= w1Range[1] & 
                   peak$w2 <= w2Range[2] & peak$w2 >= w2Range[1],]
    if (!length(peak) || nrow(peak) == 0)
      return(NULL)
    if( maxOnly )
      peak <- peak[which.max(abs(peak$Height)),]
    peak <- peak[order(peak$w1, peak$w2),]
    peak$Index <- 1:nrow(peak)
    
  }
  
  ## Return the new list
  if (!length(peak) || nrow(peak) == 0)
    return(NULL)
  
  return(peak)
}


## Internal peak picking function peakPick3D
## General local maximum (hill climbing method) for peak picking 2D slices of
## 3D files
## fileName    - name of the file to be peak picked, NULL will pick the current
## inFile			 - file data as returned by ucsf1D or ucsf2D, used instead of the
##							 fileName argument to peak pick a file without opening it
## w1Range     - w1 chemical shift range c(downfield,upfield) to be used
## w2Range     - w2 chemical shift range c(downfield,upfield) to be used
## fancy       - logical argument, FALSE implements a basic peak picker
##               that returns local maxima only, this is fastest; TRUE
##               determines chemical shifts of peaks, groups multiplets,
##               and measures line width and volume
## noiseFilt   - Integer argument that can be set to 0, 1 or 2; 
##               0 does not apply a noise filter, 1 applies a mild filter
##               (adjacent points in the direct dimension must be above the 
##               noise threshold), 2 applies a strong filter (all adjacent points
##               must be above the noise threshold
## maxOnly - logical, if TRUE only the absolute maximum peak is returned
## Note: the threshold used for peak picking taken from the graphics parameters
## ... - Additional peak picking parameters passed to trans2peak
## Returns a new peak list for the w1/w2 ranges provided
peakPick3D <- function(fileName = currentSpectrum, inFile = NULL, 
                       w1Range = NULL, w2Range = NULL, fancy = FALSE, 
                       noiseFilt = globalSettings$peak.noiseFilt, maxOnly = FALSE, ...){
  
  ## Error checking
  if (is.null(inFile))
    inFile <- fileFolder[[fileName]]	
  if (is.null(w1Range))
    w1Range <- c(inFile$file.par$upfield_ppm[1],
                 inFile$file.par$downfield_ppm[1])
  if(is.null(w2Range))
    w2Range <- c(inFile$file.par$upfield_ppm[2],
                 inFile$file.par$downfield_ppm[2])
  if (is.null(inFile$graphics.par$tiles))
    inFile <- findTiles(in.folder=inFile, internal=TRUE)
  
  ## Define some local variables
  filePar <- inFile$file.par
  graphicsPar <- inFile$graphics.par
  bs <- filePar$block_size
  ms <- filePar$matrix_size
  uf <- filePar$upfield_ppm
  df <- filePar$downfield_ppm
  endForm <- filePar$endian
  binLoc <- filePar$binary_location
  conDisp <- graphicsPar$conDisp
  thresh <- inFile$file.par$noise_est * graphicsPar$clevel
  
  ## Find best chemical shift match
  shifts <- NULL
  shifts$w1 <- seq(uf[1], df[1], length.out=ms[1])
  shifts$w2 <- seq(uf[2], df[2],	length.out=ms[2])	
  shifts$w3 <- seq(uf[3], df[3],	length.out=ms[3])	
  w1Range <- sort(w1Range); w2Range <- sort(w2Range)
  w3Range <- rep(filePar$z_value, 2)
  t1 <- findInterval(w1Range, shifts$w1, all.inside=TRUE)
  d1 <- findInterval(w2Range, shifts$w2, all.inside=TRUE)
  z1 <- findInterval(w3Range, shifts$w3, all.inside=TRUE)
  t2 <- t1 + 1; d2 <- d1 + 1; z2 <- z1 + 1
  out <- NULL
  for (i in 1:2){
    out$w1[i] <- switch(which.min(c(abs(w1Range[i] - shifts$w1[t1[i]]),
                                    abs(w1Range[i] - shifts$w1[t2[i]]))), t1[i], t2[i])
    out$w2[i] <- switch(which.min(c(abs(w2Range[i] - shifts$w2[d1[i]]),
                                    abs(w2Range[i] - shifts$w2[d2[i]]))), d1[i], d2[i])
    out$w3[i] <- switch(which.min(c(abs(w3Range[i] - shifts$w3[z1[i]]),
                                    abs(w3Range[i] - shifts$w3[z2[i]]))), z1[i], z2[i])
  }
  
  ## Invert w1/w2 selection to match binary data format
  out$w1 <- sort(ms[1] - out$w1) + 1
  out$w2 <- sort(ms[2] - out$w2) + 1
  
  ## Find sparky tiles that will be plotted
  w1Tiles <- (ceiling(out$w1[1] / bs[1]):ceiling(out$w1[2] / bs[1]))
  w2Tiles <- (ceiling(out$w2[1] / bs[2]):ceiling(out$w2[2] / bs[2]))
  w3Tiles <- (ceiling(out$w3[1] / bs[3]):ceiling(out$w3[2] / bs[3]))
  numTiles <- ceiling(ms / bs)
  tiles <- NULL
  for (i in w1Tiles)
    tiles <- c(tiles, (i - 1) * numTiles[2] + (w2Tiles - 1))
  tiles <- tiles + (w3Tiles - 1) * numTiles[1] * numTiles[2]
  
  ## Set up for tile overlap
  w1Slice <- rep(NA, numTiles[1] * bs[1])
  w1Slice <- rbind(w1Slice, w1Slice, w1Slice, w1Slice, w1Slice)
  w2Slice <- rep(NA, numTiles[2] * bs[2] + 5)
  w2Slice <- cbind(w2Slice, w2Slice, w2Slice, w2Slice, w2Slice)
  
  ## Open connection to binary file
  con <- myFile(filePar$file.name, "rb")
  w2TileNum <- ceiling(ms[2] / bs[2])
  tileSize <- bs[1] * bs[2] * bs[3] * 4
  locMax <- allPos <- NULL
  for(i in tiles){
    
    ## Define current w1/w2 range
    tile2D <- i - (w3Tiles - 1) * numTiles[1] * numTiles[2]
    w1 <- ceiling((tile2D + 1) / w2TileNum)
    w1B <- (1:bs[1]) + bs[1] * (w1 - 1) 
    w2 <- (tile2D + 1) - (w1 * w2TileNum ) + w2TileNum
    w2B <- (1:bs[2]) + bs[2] * (w2 - 1) 
    
    ##skip tiles with no data above the threshold and set w1/w2 slices to NA
    if (!i %in% graphicsPar$tiles){	
      w1Slice[1:5, w1B] <- NA
      w2Slice[w2B, 1:5] <- NA
      next() 			
    }
    
    ## Define data location for current tile
    tileLoc <- tileSize * i + binLoc
    if (bs[3] > 1)
      zPos <- (out$w3[1] - 1) * bs[3]
    else
      zPos <- 0
    dataLoc <- tileLoc + bs[1] * bs[2] * 4 * zPos
    
    ## Read binary
    seek(con, dataLoc, origin='start')
    tData <- matrix(rev(readBin(con, size=4, what='double', endian=endForm, 
                                n=bs[1] * bs[2])), ncol=bs[1])
    
    ## Make a list of all observable signals for fancy peak picking
    if (fancy){
      tObs <- obs2List(tData, conDisp=conDisp, thresh=thresh, bs=bs)
      tObs$tRow <- tObs$tRow + bs[2] * 
        floor((numTiles[1] * numTiles[2] - tile2D - 1) %% numTiles[2])
      tObs$tCol <- tObs$tCol + bs[1] * 
        floor((numTiles[1] * numTiles[2] - tile2D - 1) %/% numTiles[2])	
      allPos <- rbind(allPos, tObs)
    }
    
    ## Overlap tile, and save upfield slices for next tile
    w1Lap <- w1Slice[1:5, w1B]
    pad <- rev(w2B)[1] + c(1, 2, 3, 4, 5)
    w2Lap <- w2Slice[c(w2B, pad), 1:5]
    w1Slice[1:5, w1B] <- tData[1:5, ] 
    w2Slice[w2B, 1:5] <- tData[, 1:5]
    tData <- cbind(rbind(tData, w1Lap), w2Lap)
    
    ## Find local maxima
    if (conDisp[1] && conDisp[2])
      oThresh <- localMax(abs(tData), thresh=thresh, noiseFilt=noiseFilt) - 1
    else{
      if (conDisp[1])
        oThresh <- localMax(tData, thresh=thresh, noiseFilt=noiseFilt) - 1
      else
        oThresh <- localMax(-tData, thresh=thresh, noiseFilt=noiseFilt) - 1			
    }
    
    ## Format local maxima
    peak <- data.frame(list(tRow=(oThresh %% (bs[2] + 5)) + 1, 
                            tCol=(oThresh %/% (bs[2] + 5)) + 1, 
                            Height=tData[oThresh + 1]))
    peak <- peak[peak[, 1] > 1 & peak[, 1] < bs[2] + 2, ]		
    peak <- peak[peak[, 2] > 1 & peak[, 2] < bs[1] + 2, ]
    if (!nrow(peak))
      next()
    
    ## Translate tile row/column to spectrum row/column
    peak$tRow <- peak$tRow + bs[2] * 
      floor((numTiles[1] * numTiles[2] - tile2D - 1) %% numTiles[2])
    peak$tCol <- peak$tCol + bs[1] * 
      floor((numTiles[1] * numTiles[2] - tile2D - 1) %/% numTiles[2])
    locMax <- rbind(locMax, peak)
    
  }
  close(con)
  
  ## Clean up peak list
  if (is.null(locMax))
    return(NULL)
  locMax <- unique(locMax)
  
  ## Match point locations with chemical shifts
  w1 <- seq(uf[1], df[1],	length.out=ms[1])
  if (ms[1] %% bs[1])
    w1 <- c(rep(NA, bs[1] - ms[1] %% bs[1]), w1)
  w2 <- seq(uf[2], df[2], length.out=ms[2])
  if (ms[2] %% bs[2])
    w2 <- c(rep(NA, bs[2] - ms[2] %% bs[2] ), w2)
  w1Range <- sort(rev(w1)[range(out$w1)])
  w2Range <- sort(rev(w2)[range(out$w2)])
  
  if (fancy){
    peak <- trans2Peak(allPos=allPos, locMax=locMax, w1=w1, w2=w2, 
                       w1Range=w1Range, w2Range=w2Range, ...)
  }else{
    
    peak <- data.frame(list(Index=1:nrow(locMax), w1=w1[locMax$tCol],
                            w2=w2[locMax$tRow], Height=locMax$Height, 
                            Assignment=rep("NA", length(locMax$tCol))), stringsAsFactors=FALSE)
    
    ## Filter the outgoing list to match the data range
    peak <- peak[peak$w1 <= w1Range[2] & peak$w1 >= w1Range[1] & 
                   peak$w2 <= w2Range[2] & peak$w2 >= w2Range[1], ]
    if (!nrow(peak))
      return(NULL)
    if (maxOnly)
      peak <- peak[which.max(abs(peak$Height)), ]
    peak <- peak[order(peak$w1, peak$w2), ]
    peak$Index <- 1:nrow(peak)
  }
  
  ## Return the new list
  if (!nrow(peak))
    return(NULL)
  
  return(peak)
}


## Internal peak picking function peakPickRsd
## General local maximum (hill climbing method) for peak picking RSD files
## fileName    - name of the file to be peak picked, NULL will pick the current
## inFile			 - file data as returned by ucsf1D or ucsf2D, used instead of the
##							 fileName argument to peak pick a file without opening it
## w1Range     - w1 chemical shift range c(downfield,upfield) to be used
## w2Range     - w2 chemical shift range c(downfield,upfield) to be used
## fancy       - logical argument, FALSE implements a basic peak picker
##               that returns local maxima only, this is fastest; TRUE
##               determines chemical shifts of peaks, groups multiplets,
##               and measures line width and volume
## noiseFilt   - Integer argument that can be set to 0, 1 or 2; 
##               0 does not apply a noise filter, 1 applies a mild filter
##               (adjacent points in the direct dimension must be above the 
##               noise threshold), 2 applies a strong filter (all adjacent points
##               must be above the noise threshold
## maxOnly - logical, if TRUE only the absolute maximum peak is returned
## Note: the threshold used for peak picking taken from the graphics parameters
## ... - Additional peak picking parameters passed to trans2peak
## Returns a new peak list for the w1/w2 ranges provided
peakPickRsd <- function(fileName=currentSpectrum, inFile=NULL, w1Range=NULL, 
                        w2Range=NULL, fancy=FALSE, noiseFilt=globalSettings$peak.noiseFilt, 
                        maxOnly=FALSE, ...){
  
  ## Error checking
  if (is.null(inFile))
    inFile <- fileFolder[[fileName]]	
  if (is.null(w1Range))
    w1Range <- c(inFile$file.par$upfield_ppm[1],
                 inFile$file.par$downfield_ppm[1])
  if (is.null(w2Range))
    w2Range <- c(inFile$file.par$upfield_ppm[2],
                 inFile$file.par$downfield_ppm[2])
  if (is.null(inFile$graphics.par$tiles))
    inFile <- findTiles(in.folder=inFile, internal=TRUE)
  
  ## Define some local variables
  filePar <- inFile$file.par
  graphicsPar <- inFile$graphics.par
  conDisp <- graphicsPar$conDisp
  thresh <- filePar$noise_est * graphicsPar$clevel
  w1Range <- sort(w1Range)
  w2Range <- sort(w2Range)
  shifts <- matchShift(inFolder=inFile, w1=w1Range, w2=w2Range, return.seq=TRUE)
  winW1 <- shifts$w1
  winW2 <- shifts$w2
  
  ## Find tiles that will be plotted
  tiles <- NULL
  upShifts <- filePar$block_upfield_ppms
  downShifts <- filePar$block_downfield_ppms
  blockSizes <- filePar$block_size
  for (tNum in seq_along(filePar$block_size$w1)){
    
    ## Get the chemical shifts for the current tile
    blockW1 <- seq(upShifts$w1[tNum], downShifts$w1[tNum], 
                   length.out=blockSizes$w1[tNum])
    blockW2 <- seq(upShifts$w2[tNum], downShifts$w2[tNum],
                   length.out=blockSizes$w2[tNum])
    
    ## Check the window for the presence of any shift in the current block
    if (any(round(blockW1, 3) %in% round(winW1, 3)) && 
        any(round(blockW2, 3) %in% round(winW2, 3)))
      tiles <- c(tiles, tNum)
  }
  
  ## Define data locations for each block
  blockLoc <- filePar$binary_location
  for (i in seq_along(filePar$block_size$w1))
    blockLoc <- c(blockLoc, blockLoc[i] + 4 * filePar$block_size$w1[i] * 
                    filePar$block_size$w2[i])
  
  ## Peak pick the spectrum one tile at a time
  con <- myFile(filePar$file.name, "rb")
  peaks <- NULL
  for (tNum in tiles){
    
    ## Skip tiles with no data above the noise threshold
    if (!(tNum - 1) %in% graphicsPar$tiles)
      next
    
    ## Read data
    bs <- c(filePar$block_size$w1[tNum], filePar$block_size$w2[tNum])
    seek(con, blockLoc[tNum], origin='start')
    tData <- matrix(rev(readBin(con, size=4, what='double', 
                                endian=filePar$endian, n=bs[1] * bs[2])), 
                    ncol=filePar$block_size$w1[tNum])
    
    ## Find local maxima
    if (conDisp[1] && conDisp[2])
      oThresh <- localMax(abs(tData), thresh=thresh, noiseFilt=noiseFilt) - 1
    else{
      if (conDisp[1])
        oThresh <- localMax(tData, thresh=thresh, noiseFilt=noiseFilt) - 1
      else
        oThresh <- localMax(-tData, thresh=thresh, noiseFilt=noiseFilt) - 1			
    }
    
    ## Convert maxima indices to row number, column number, and height
    locMax <- data.frame(list(tRow=oThresh %% bs[2] + 1, 
                              tCol=oThresh %/% bs[2] + 1, 
                              Height=tData[oThresh + 1]))
    if (!nrow(locMax))
      next
    
    ## Define chemical shifts for current tile
    w1 <- seq(upShifts$w1[tNum], downShifts$w1[tNum], length.out=bs[1])
    w2 <- seq(upShifts$w2[tNum], downShifts$w2[tNum], length.out=bs[2])
    
    if (fancy){
      
      ## Pass a list of all observable signals to the fancy peak picker
      tObs <- obs2List(tData, conDisp=conDisp, thresh=thresh, bs=bs)		
      peak <- trans2Peak(allPos=tObs, locMax=locMax, w1=w1, w2=w2, 
                         w1Range=range(w1), w2Range=range(w2), ...)
    }else{
      
      ## Format local maxima in to a standard peak list
      peak <- data.frame(list(Index=1:nrow(locMax), 
                              w1=w1[locMax$tCol],
                              w2=w2[locMax$tRow],
                              Height=locMax$Height,
                              Assignment=rep("NA", length(locMax$tCol))), 
                         stringsAsFactors=FALSE)
    }
    peaks <- rbind(peaks, peak)
  }
  close(con)
  
  ## Filter the outgoing list to match the data range
  peaks <- peaks[peaks$w1 <= w1Range[2] & peaks$w1 >= w1Range[1] & 
                   peaks$w2 <= w2Range[2] & peaks$w2 >= w2Range[1], ]
  if (!nrow(peaks))
    return(NULL)
  
  ## Filter list to include only the peaks with the absolute maximum intensity
  if (!fancy && maxOnly)
    peaks <- peaks[which.max(abs(peaks$Height)), ]
  
  ## Format and return the new list
  peaks <- peaks[order(peaks$w1, peaks$w2), ]
  peaks$Index <- 1:nrow(peaks)
  return(peaks)
}


## Internal 2D fancy peak picking helper function
trans2Peak <- function( allPos, locMax, w1, w2, w1Range, w2Range, ...){
  
  ## Filter the outgoing list to match the current data range
  locMax <- data.frame( list(	
    tCol = locMax$tCol,
    tRow = locMax$tRow,
    w1 = w1[locMax$tCol],
    w2 = w2[locMax$tRow], 
    Height = locMax$Height))	
  locMax <- locMax[locMax$w1 <= w1Range[2] & locMax$w1 >= w1Range[1] & 
                     locMax$w2 <= w2Range[2] & locMax$w2 >= w2Range[1],]
  if(nrow(locMax) == 0)
    return(NULL)
  
  ## Group transitions together
  locMax <- fancyPeak2D( allPos = allPos, locMax = locMax, ...)
  locMax$rStart <- w2[locMax$rStart]
  locMax$rEnd <- w2[locMax$rEnd]
  locMax$cStart <- w1[locMax$cStart]
  locMax$cEnd <- w1[locMax$cEnd]
  
  peak <- NULL
  j <- 1
  for( i in unique(locMax$Index) ){
    tSub <- locMax[locMax$Index == i, ]
    tSub$Index <- j
    tSub$Multiplet <- 1:nrow(tSub)
    
    tSub <- rbind(c( mean(tSub$w1),  mean(tSub$w2), 
                     tSub$Height[which.max(abs(tSub$Height))], 
                     min(tSub$rStart), max(tSub$rEnd),
                     min(tSub$cStart), max(tSub$cEnd), j, NA), tSub[,-(1:2)])
    
    peak <- rbind(peak, tSub)
    j <- j + 1
  }		
  
  peak$w1D <- peak$cEnd - peak$cStart
  peak$w2D <- peak$rEnd - peak$rStart
  peak$Assignment <- rep("NA", nrow(peak))		
  peak <- peak[,match(c('Index', 'w1', 'w2', 'Height', 'Assignment', 
                        'Multiplet', 'w1D', 'w2D'), names(peak))]
  peak <- peak[order(peak$w1, peak$w2),]
  
  return(peak)
}

## Internal 2D peak finding function 
## allPos  - a data frame with columns labled tRow, tCol, and Height that define
##           the row, column and intensities of observed peaks in a 2D matrix
## w1Gran  - Integer 0 or greater, defines the number of sub threshold points a 
##           peak can cross in either direction of the indirect dimension
## w2Gran  - Integer 0 or greater, defines the number of sub threshold points a 
##           peak can cross in either direction of the direct dimension
## Note:   - Granularities of 0 mean that multiplets must be in the same 
##           row/column and none of the points are separated by points below
##           the noise threshold.
## returns - a data frame with the rows and columns defining each peak 
##           with each transition and chemical shift indicated
fancyPeak2D <- function( allPos, locMax, w1Gran = 1, w2Gran = 3 ){
  
  if( w1Gran < 0 )
    stop('w1Gran must be an integer greater than 0')
  if( w2Gran < 0 )
    stop('w2Gran must be an integer greater than 0')
  
  ## Find the w2 row breaks defining each peak 
  allPos <- allPos[order(allPos$tRow, decreasing = FALSE),]
  rowFilt <- NULL	
  for( i in unique(allPos$tCol) ){
    tmp <- allPos[allPos$tCol ==  i, ]
    breaks <-  which((c(tmp$tRow, NA) - c(NA, tmp$tRow + 1) != 0))
    if( length(breaks) == 0 ){
      sBr <- 1
      eBr <- nrow(tmp)
      mRow <- which.max(abs(tmp$Height))
      rowFilt <- rbind(rowFilt, cbind(tmp$tRow[sBr], tmp$tRow[eBr], 
                                      tmp$tRow[mRow], tmp$tCol[mRow], tmp$tCol[mRow],tmp$tCol[mRow], 
                                      tmp$Height[mRow]))
      next()
    }
    
    sBr <- sort(c(1, breaks - 1, breaks, nrow(tmp) ))
    eBr <- sBr[ gl( 2, 1, length = length(sBr) ) == 2 ]
    sBr <- sBr[ gl( 2, 1, length = length(sBr) ) == 1 ]
    
    for( j in 1:length(sBr) ){
      subTmp <- tmp[sBr[j]:eBr[j],]
      mRow <- which.max(abs(subTmp$Height))
      rowFilt <- rbind(rowFilt, cbind(tmp$tRow[sBr[j]], tmp$tRow[eBr[j]], 
                                      subTmp$tRow[mRow], subTmp$tCol[mRow], subTmp$tCol[mRow], 
                                      subTmp$tCol[mRow], subTmp$Height[mRow]))
    }		
  }
  
  ## Find the w1 column breaks defining each peak 
  allPos <- allPos[order(allPos$tCol, decreasing = FALSE),]
  colFilt <- NULL	
  for( i in unique(allPos$tRow) ){
    tmp <- allPos[allPos$tRow ==  i, ]
    breaks <-  which((c(tmp$tCol, NA) - c(NA, tmp$tCol + 1) != 0))
    if( length(breaks) == 0 ){
      sBr <- 1
      eBr <- nrow(tmp)
      mRow <- which.max(abs(tmp$Height))
      colFilt <- rbind(colFilt, cbind(tmp$tRow[mRow], tmp$tRow[mRow], 
                                      tmp$tRow[mRow], tmp$tCol[sBr], tmp$tCol[eBr], tmp$tCol[mRow], 
                                      tmp$Height[mRow]))
      next()
    }
    
    sBr <- sort(c(1, breaks - 1, breaks, nrow(tmp) ))
    eBr <- sBr[ gl( 2, 1, length = length(sBr) ) == 2 ]
    sBr <- sBr[ gl( 2, 1, length = length(sBr) ) == 1 ]
    
    for( j in 1:length(sBr) ){
      subTmp <- tmp[sBr[j]:eBr[j],]
      mRow <- which.max(abs(subTmp$Height))
      colFilt <- rbind(colFilt, cbind(subTmp$tRow[mRow], subTmp$tRow[mRow], 
                                      subTmp$tRow[mRow], tmp$tCol[sBr[j]], tmp$tCol[eBr[j]], 
                                      subTmp$tCol[mRow], subTmp$Height[mRow]))
    }		
  }
  
  rowFilt <- data.frame(rowFilt)
  colFilt <- data.frame(colFilt)
  names(rowFilt) <- names(colFilt) <- c('rStart', 'rEnd', 'rMax', 'cStart', 
                                        'cEnd', 'cMax', 'Height')
  
  ## Find the peak boundries for each transition
  trans <- data.frame(list( tRow = NA, tCol = NA, Height = NA, rStart = NA, 
                            rEnd = NA, cStart = NA, cEnd = NA))[-1,]
  for( i in 1:nrow(locMax) ){
    mSub <- locMax[i, ]
    rowSub <- rowFilt[ rowFilt$cStart == mSub$tCol, ]
    rowSub <- rowSub[ rowSub$rStart <= mSub$tRow &
                        rowSub$rEnd >= mSub$tRow, ]					
    colSub <- colFilt[ colFilt$rStart == mSub$tRow, ]
    colSub <- colSub[ colSub$cStart <= mSub$tCol &
                        colSub$cEnd >= mSub$tCol, ]
    tSub <- rbind(rowSub, colSub)
    trans <- rbind(trans, cbind(mSub, min(tSub$rStart), max(tSub$rEnd), 
                                min(tSub$cStart), max(tSub$cEnd)))
  }
  names(trans) <- c('tRow', 'tCol', 'w1', 'w2', 'Height', 'rStart', 
                    'rEnd', 'cStart','cEnd')
  trans$Index <- 1:nrow(trans)
  
  ## Group transitions together and generate a master peak summary 
  j <- 1
  for( i in trans$Index ){
    
    tSub <- trans[ trans$Index == i , ]
    if( nrow(tSub) == 0 )
      next()
    
    tmp <- trans[trans$rStart >= (tSub$rStart - w2Gran),]
    tmp <- tmp[tmp$rEnd <= (tSub$rEnd + w2Gran), ]
    tmp <- tmp[tmp$cStart >= (tSub$cStart - w1Gran), ]
    tSub <- tmp[tmp$cEnd <= (tSub$cEnd + w1Gran), ]
    trans[ which(!is.na(match(trans$Index, tSub$Index))), ]$Index <- j	
    
    j <- j + 1
  }
  
  return(trans)
}

## Internal helper function for fancy peak picking
## Converts a data matrix into a list with the rows and columns of 
##          signals that are above the threshold 
## x        - Numeric matrix of data to be searched
## Condisp  - Logical vector of the form c(TRUE, TRUE), allows independent 
##            thresholding of negative and positive values 
## thresh   - Numeric value specifying the absolute limit for observable signals
## bs       - Integer vector of the form c(columns, rows) defining the 
##            dimensions of the matrix
## returns a list of observable signals with the columns tRow, tCol, and Height
obs2List <- function(x, conDisp, thresh, bs = c(ncol(x), nrow(x))){
  
  ## Threshold data
  if( conDisp[1] && conDisp[2] )
    x[ which(abs(x) < thresh) ] <- NA
  else{
    if( conDisp[1])
      x[ x < thresh ] <- NA
    else
      x[ (-x) < thresh ] <- NA	
  }
  oThresh <- which(!is.na(x)) - 1
  if( length(oThresh) < 1 )
    return(NULL)
  
  ## Convert tile row/columns to overall matrix row/column numbers 
  return(data.frame(list(
    tRow = (oThresh %% bs[2]) + 1, 
    tCol = (oThresh %/% bs[2]) + 1, 
    Height = x[oThresh + 1])))	
}

## Internal 2D peak picking function
## Finds points in a matrix that are larger than all surrounding points
## x  -  A numeric matrix containing the range of data to be peak picked
## thresh    - Numeric value specifying the minimum level to be included
## noiseFilt - Integer argument that can be set to 0, 1 or 2; 
##              0 does not apply a noise filter, 1 applies a mild filter
##              (adjacent points in the direct dimension must be above the 
##              noise threshold), 2 applies a strong filter (all adjacent points
##              must be above the noise threshold
## Returns a vector of points defining the local maxima
localMax <- function(x, thresh, noiseFilt ){
  
  nC <- ncol(x)
  nR <- nrow(x)
  if( noiseFilt == 2 )
    x[x < thresh] <- NA  
  
  ## Find row/column local maxes
  if( noiseFilt == 1 ){
    y <- x
    y[ y < thresh ] <- NA
    vMax <- intersect(which(c(NA, y) < c(y, NA)), which(c(NA, y) > c(y, NA))-1)
  }else
    vMax <- intersect(which(c(NA, x) < c(x, NA)), which(c(NA, x) > c(x, NA))-1)
  x <- t(x)
  hMax <- intersect(which(c(NA, x) < c(x, NA)), which(c(NA, x) > c(x, NA))-1)-1
  hMax <- (hMax %% nC * nR) + hMax %/% nC + 1
  
  
  ## Find diagonal maxima
  x <- t(x)
  hvMax <- intersect(vMax, hMax)
  if( noiseFilt == 0 )
    hvMax <- hvMax[ x[hvMax] > thresh ]
  dMax <- cbind(hvMax, hvMax - nR + 1, hvMax - nR - 1, hvMax + nR + 1, 
                hvMax + nR - 1)
  dMax[dMax < 1 | dMax > nC*nR ] <- NA
  dMax <- which(max.col(cbind( x[dMax[,1]], x[dMax[,2]], x[dMax[,3]], 
                               x[dMax[,4]], x[dMax[,5]])) == 1)
  
  return(hvMax[dMax])
}

## Estimates volumes of 2D peaks using stacked elipsoids
## inFile  - file parameters and data for desired spectrum as returned by ed()
## gran       - Integer indicating granularity of contour fitting 
## c.vol      - Logical argument, TRUE returns stacked elipsoid volumes, 
##               FALSE returns sum of visible data
## note: the default, contour.vol = FALSE, seems to work best
## baselineCorr - local baseline correction for 1D 
## Returns volume for 2D ROIs and area of 1D ROIs
peakVolume <- function(inFile, gran = 200, c.vol = FALSE, 
                       baselineCorr = FALSE){
  
  ## Volume estimates for 2D data
  if(inFile$file.par$number_dimensions > 1){
    if(!c.vol ){
      volume <- sum(inFile$data[inFile$data >= 
                                  inFile$file.par$noise_est * inFile$graphics.par$clevel]) 
      if(length(volume) == 0)
        volume <- NA
    }else{
      ## Generate contour lines for data
      zlim <- c(inFile$file.par$noise_est *
                  inFile$graphics.par$clevel, max(inFile$data))
      c.levels <- seq(zlim[1], zlim[2], length.out=gran )
      c.int <- diff(c.levels[1:2])
      contour.file <- contourLines(z = inFile$data, levels = c.levels )
      
      ## Estimate volume as stacked elipsoids
      if(length(contour.file) > 0 && zlim[1] < zlim[2]){
        volume <- NULL 
        for(i in 1:length(contour.file))
          volume <- sum(volume, (4 / 3 * pi * diff(range(contour.file[[i]]$x)) *
                                   diff(range(contour.file[[i]]$y)) * c.int ))
      }else
        volume <- NA
    }  
    
    ## Area estimate for 1D data
  }else{
    
    ## Local baseline correction
    if( baselineCorr )
      inFile$data <- inFile$data - fivenum(inFile$data)[2]
    volume <- sum(inFile$data)
    
  }  
  return(volume)
}

################################################################################
##                                                                            ##
##                    Internal graphics functions                             ##
##                                                                            ##
################################################################################

## Internal graphics function bringFocus
## This is a platform independant version of bringToTop
bringFocus <- function(dev = -1){
  if( .Platform$OS.type == 'windows' && sdiCheck(FALSE) )
    # bringToTop( dev )
    return(invisible())
}

## Internal graphics function setWindow
## Makes a new plotting window with the correct title, width/height or sets
## an existing window as the active device
## p.window - The window type, can be 'main', 'sub', 'multi', or 'stats'
## ...      - Additional R graphics parameters, see par()
## returns a new graphics device 
setWindow <- function( p.window = 'main', ...){
  
  ## Check to see if window is open
  devNum <- which(c('main', 'sub', 'multi', 'stats') == p.window) + 1
  if(length(which(dev.list() == devNum)) == 0 ){
    odev <- dev.list()
    devTitle <- switch(devNum - 1, "Main Plot Window", "ROI Subplot Window",
                       "Multiple File Window")
    devWidth <- switch(devNum - 1, globalSettings$size.main[1], 
                       globalSettings$size.sub[1], globalSettings$size.multi[1])
    devHeight <- switch(devNum - 1, globalSettings$size.main[2], 
                        globalSettings$size.sub[2], globalSettings$size.multi[2])
    
    ## Open new graphics window 
    while(dev.cur()  != devNum){
      if (.Platform$OS == 'windows')
        dev.new(title = devTitle, width = devWidth, height = devHeight)
      else
        X11(title = devTitle, width = devWidth,	height = devHeight)
    }
    for(i in dev.list()){
      if(length(which(c(odev, devNum) == i)) == 0)
        dev.off(i)    
    }
    if (p.window == 'main'){
      par(defaultSettings[1:65])
      par(mar=globalSettings$mar)
    }
    par(...)
    #gui(p.window)
    popupGui(p.window)
  }else{
    dev.set(devNum)
    bringFocus(devNum)
    par(...)
  }
  invisible()
}

## Internal function to cycle through open graphics windows
## note: A device can specified by setting dev (used internally)
cw <- function(dev=NULL){
  if(is.null(dev))
    dev.set(which = dev.next())
  else
    dev.set(which = dev)
  bringFocus(dev.cur())
  
  ## Leave the console active
  bringFocus(-1)  
}

## Internal graphics function set.graphics
## General function for changing graphics settings used by lower level functions
## file.name - a list of file names to be modified, default is the current file,
##             names must match names from fileFolder
## all.files - logical argument, if TRUE all files will be updated
## save.backup - logical argument, TRUE will save a copy of the environment
##               for undo/redo, FALSE updates the file folder without saving
##               a backup copy. 
## refresh.graphics - logical argument, refreshes the active plots if TRUE
## par() parameters: bg, fg, col.axis, col.lab, col.main, col.sub, col, and usr 
##          are R arguments, see par for documentation
## line.color  - sets fg, col, col.axis, col.lab, col.main, and col.sub to 
##               a single input color
## drawPeptides parameters: pos.color, neg.color, proj.color, conDisp, nlevels,
##                     clevel, w1Range, w2Range, and type. See drawPeptides for 
##                     documentation
## perspective parameters: theta, phi, asp, see persp for documentation
## peak display parameters: peak.color, peak.disp, noiseFilt see pdisp
## 1D file settings: thresh.1D, position.1D, offset, and overlay.text, 
##									 see peakPick2D, vp, and overlay for documentation
## 1D projections: proj.mode, proj.type, proj.direct, filter, see proj1D
## roi parameters: roi.multi, roiMain, roiMax, roi.bcolor, roi.tcolor
## Note: All changes are applied to the files in fileFolder, except:
##       offset, peak.disp, thresh.1D, position.1D, roiMain, roiMax, 
##				proj.direct, and filter.
##       These exceptions are modified in the globalSettings list
setGraphics <- function (file.name = currentSpectrum, all.files = FALSE, 
                         save.backup = TRUE, refresh.graphics = FALSE, bg = NULL, fg = NULL, 
                         col.axis = NULL, col.lab = NULL, col.main = NULL,	col.sub = NULL, 
                         col = NULL, usr = NULL, line.color = NULL, pos.color = NULL, 
                         neg.color = NULL, proj.color = NULL, conDisp = NULL, nlevels = NULL,	
                         clevel = NULL, type = NULL, theta = NULL, phi = NULL,	asp = NULL,	
                         peak.color = NULL, peak.disp = NULL, peak.noiseFilt = NULL, peak.pch = NULL,
                         peak.cex = NULL, peak.labelPos = NULL, thresh.1D = NULL, position.1D = NULL, 
                         offset = NULL, proj.mode = NULL, proj.type = NULL,	proj.direct = NULL, 
                         filter = NULL, roi.multi = NULL, roiMain = NULL, roiMax = NULL, 
                         roi.bcolor = NULL, roi.tcolor = NULL, roi.lwd = NULL, roi.lty = NULL, 
                         roi.cex = NULL, roi.labelPos = NULL, roi.noiseFilt = NULL, roi.w1 = NULL, 
                         roi.w2 = NULL, roi.pad = NULL, w1Range = NULL, w2Range = NULL, 
                         overlay.text = NULL, overlay.textSuppress = NULL){
  
  ## Set global graphics changes
  if( !is.null(offset) ){
    if(is.numeric(offset))
      globalSettings$offset <- offset
    else
      log_message('offset must be a numeric value', quote = FALSE)		
  }
  if(!is.null(proj.mode)){
    if( is.logical(proj.mode))
      globalSettings$proj.mode <- proj.mode
    else
      log_message('proj.mode must be either TRUE or FALSE', quote = FALSE)		
  }
  if(!is.null(proj.type)){
    if(proj.type %in% c('l', 'p', 'b'))
      globalSettings$proj.type <- proj.type
    else
      log_message('proj.type must be either "l", "p", or "b"', quote = FALSE)			
  }
  if(!is.null(proj.direct)){
    if(proj.direct == 1 || proj.direct == 2 )
      globalSettings$proj.direct <- proj.direct
    else
      log_message('proj.direct must be either 1 or 2', quote = FALSE)				
  }
  if(!is.null(filter)){
    if(is.function(filter))
      globalSettings$filter <- filter
    else
      log_message('filter must be a function', quote = FALSE)	
  }
  if( !is.null(peak.disp) ){
    if(is.logical(peak.disp))
      globalSettings$peak.disp <- peak.disp
    else
      log_message('peak.disp must be either TRUE or FALSE', quote = FALSE)		
  }
  if( !is.null(peak.noiseFilt) ){
    if( any(peak.noiseFilt == c(0, 1, 2)) )
      globalSettings$peak.noiseFilt <- peak.noiseFilt
    else
      log_message('peak.noiseFilt must be 0, 1, or 2', quote = FALSE)			
  }
  if (!is.null(peak.pch)){
    if (is.numeric(peak.pch) || nchar(peak.pch) == 1)
      globalSettings$peak.pch <- peak.pch
    else
      log_message('peak.pch must be numeric or a single ASCII character', quote=FALSE)			
  }
  if (!is.null(peak.cex)){
    if (is.numeric(peak.cex))
      globalSettings$peak.cex <- peak.cex
    else
      log_message('peak.cex must be a numeric value', quote=FALSE)			
  }
  if (!is.null(peak.labelPos)){
    if (peak.labelPos %in% c('top', 'bottom', 'left', 'right', 'center'))
      globalSettings$peak.labelPos <- peak.labelPos
    else
      log_message(paste('peak.labelPos must be either "top", "bottom", "left",', 
                  '"right", or "center"'), quote=FALSE)			
  }
  if( !is.null(thresh.1D) ){
    if(is.numeric(thresh.1D))
      globalSettings$thresh.1D <- thresh.1D
    else
      log_message('thresh.1D must be a numeric value', quote = FALSE)			
  }
  if( !is.null(position.1D) ){
    if(is.numeric(position.1D))
      globalSettings$position.1D <- position.1D
    else
      log_message('position.1D must be a numeric value', quote = FALSE)		
  }
  if( !is.null(roiMain)){
    if(is.logical(roiMain))
      globalSettings$roiMain <- roiMain
    else
      log_message('roiMain must be either TRUE or FALSE', quote = FALSE)		
  }
  if( !is.null(roiMax)){
    if(is.logical(roiMax))
      globalSettings$roiMax <- roiMax
    else
      log_message('roiMax must be either TRUE or FALSE', quote = FALSE)			
  }
  if(!is.null(roi.bcolor)){
    colErr <- TRUE
    if(length(roi.bcolor) == 2){
      colErr <- FALSE
      for (i in roi.bcolor){
        colTest <- try(col2rgb(i), silent=TRUE)
        if (class(colTest) == 'try-error')
          colErr <- TRUE
      }
    }
    if (!colErr)
      globalSettings$roi.bcolor <- roi.bcolor
    else
      log_message('ROI box color must be a vector of colors of length 2', 
            quote = FALSE)	
  }
  if(!is.null(roi.tcolor)){
    colErr <- TRUE
    if(length(roi.tcolor) == 2){
      colErr <- FALSE
      for (i in roi.tcolor){
        colTest <- try(col2rgb(i), silent=TRUE)
        if (class(colTest) == 'try-error')
          colErr <- TRUE
      }
    }
    if (!colErr)
      globalSettings$roi.tcolor <- roi.tcolor
    else
      log_message('ROI text color must be a vector of colors of length 2', 
            quote = FALSE)	
  }
  if (!is.null(roi.lwd)){
    if (is.numeric(roi.lwd) && length(roi.lwd) == 2)
      globalSettings$roi.lwd <- roi.lwd
    else
      log_message('roi.lwd must be a numeric vector of length 2', quote=FALSE)			
  }
  if (!is.null(roi.lty)){
    if (roi.lty %in% c('solid', 'dashed', 'dotted', 'dotdash', 'longdash', 
                       'twodash', 'blank') && length(roi.lty) == 2)
      globalSettings$roi.lty <- roi.lty
    else
      log_message(paste('roi.lty must be a vector of length 2. Valid options include', 
                  '"solid", "dashed", "dotted", "dotdash", "longdash", "twodash",',
                  'or "blank"'), quote=FALSE)			
  }
  if (!is.null(roi.cex)){
    if (is.numeric(roi.cex) && length(roi.cex) == 2)
      globalSettings$roi.cex <- roi.cex
    else
      log_message('roi.cex must be a numeric vector of length 2', quote=FALSE)			
  }
  if (!is.null(roi.labelPos)){
    if (roi.labelPos %in% c('top', 'bottom', 'left', 'right', 'center'))
      globalSettings$roi.labelPos <- roi.labelPos
    else
      log_message(paste('roi.labelPos must be either "top", "bottom", "left",', 
                  '"right", or "center"'), quote=FALSE)					
  }
  if (!is.null(roi.noiseFilt)){
    if (roi.noiseFilt %in% c(0, 1, 2))
      globalSettings$roi.noiseFilt <- roi.noiseFilt
    else
      log_message('roi.noiseFilt must be 0, 1, or 2', quote = FALSE)			
  }
  if (!is.null(roi.w1)){
    if (is.numeric(roi.w1))
      globalSettings$roi.w1 <- roi.w1
    else
      log_message('roi.w1 must be a numeric value', quote=FALSE)			
  }
  if (!is.null(roi.w2)){
    if (is.numeric(roi.w2))
      globalSettings$roi.w2 <- roi.w2
    else
      log_message('roi.w2 must be a numeric value', quote=FALSE)			
  }
  if (!is.null(roi.pad)){
    if (is.numeric(roi.pad))
      globalSettings$roi.pad <- roi.pad
    else
      log_message('roi.pad must be a numeric value', quote=FALSE)			
  }
  if(!is.null(overlay.text))
  {
    if( is.logical(overlay.text))
      globalSettings$overlay.text <- overlay.text
    else
      log_message('overlay.text must be either TRUE or FALSE', quote = FALSE)		
  }
  if(!is.null(overlay.textSuppress))
  {
    if( is.logical(overlay.textSuppress))
      globalSettings$overlay.textSuppress <- overlay.textSuppress
    else
      log_message('overlay.textSuppress must be either TRUE or FALSE', quote = FALSE)		
  }
  
  ## Assign global parameters
  myAssign( 'globalSettings', globalSettings, save.backup = FALSE )
  
  
  ## Define the files to be modified
  if(all.files)
    current <- 1:length(fileFolder)
  else     
    current <- match(file.name, names(fileFolder))
  
  ## Update file settings
  for( i in current){  
    current.gpar <- fileFolder[[i]]$graphics.par
    current.fpar <- fileFolder[[i]]$file.par
    nDim <- fileFolder[[i]]$file.par$number_dimensions
    
    
    ## Set all line colors
    if(!is.null(line.color))
      fg <- col.axis <- col.lab <- col.main <- col.sub <- col <- line.color
    
    ## Change general graphics parameters related to par
    if(!is.null(bg))
      current.gpar$bg <- bg
    if(!is.null(fg))
      current.gpar$fg <- fg
    if(!is.null(col.axis))
      current.gpar$col.axis <- col.axis
    if(!is.null(col.lab))  
      current.gpar$col.lab <- col.lab
    if(!is.null(col.main))
      current.gpar$col.main <- col.main
    if(!is.null(col.sub))
      current.gpar$col.sub <- col.sub
    if(!is.null(col))
      current.gpar$col <- col
    if(!is.null(usr)){
      usr <- suppressWarnings(as.numeric(usr))
      if( !any(is.na(usr)) && length(usr) == 4){
        if( diff(usr[1:2]) == 0 || diff(usr[3:4]) == 0  )
          log_message('The plot ranges must be greater than zero', quote = FALSE)
        else
          current.gpar$usr <- usr
      }else
        log_message(paste('The plot vector "usr" must be a numeric vector with four', 
                    'elements'),	quote = FALSE)
    } 
    
    ## Change digestR graphics parameters
    if(!is.null(pos.color))
      current.gpar$pos.color <- pos.color
    if(!is.null(neg.color))
      current.gpar$neg.color <- neg.color
    if(!is.null(proj.color))
      current.gpar$proj.color <- proj.color
    if(!is.null(conDisp) && is.logical(conDisp) && length(conDisp) == 2 )
      current.gpar$conDisp <- conDisp	
    if(!is.null(nlevels)){
      nlevels <- suppressWarnings(as.integer(nlevels))
      if( is.na(nlevels) || nlevels < 1 || nlevels > 1000)
        err(paste('The number of contour levels must be an integer greater', 
                  'than 0 and less than or equal to 1000'))
      else
        current.gpar$nlevels <- nlevels
    }
    if(!is.null(clevel)){
      clevel <- suppressWarnings(as.numeric(clevel))
      if( is.na(clevel) || clevel < 0 )
        err('The minimum contour level must be greater than 0')
      else{
        if( clevel != current.gpar$clevel)
          current.gpar$tiles <- NULL
        current.gpar$clevel <- clevel
      }
    }
    if(!is.null(type)){
      if (nDim == 1){
        if( any (type == c('auto', 'l', 'p', 'b' )))
          current.gpar$type <- type
        else{
          cat(paste('Plot type must be a string (in quotes) named:', '\n', 
                    'auto, l, p, or b', '\n'))
        }
      }else{
        if( any (type == c('auto', 'image', 'contour', 'filled', 'persp' )))
          current.gpar$type <- type
        else{
          cat(paste('Plot type must be a string (in quotes) named:', '\n', 
                    'auto, image, filled, or persp', '\n'))
        }
      }
    }		
    if(!is.null(theta)){
      if(is.numeric(theta))
        current.gpar$theta <- theta
      else
        log_message('theta must be a numeric value', quote = FALSE)	
      
    }
    if(!is.null(phi)){
      if( is.numeric(phi) )
        current.gpar$phi <- phi 	
      else
        log_message('phi must be a numeric value', quote = FALSE)	
    }                         
    if(!is.null(asp)){
      if( is.numeric(asp) )
        current.gpar$asp <- asp 	
      else
        log_message('asp must be a numeric value', quote = FALSE)			
    }
    if(!is.null(peak.color))
      current.gpar$peak.color <- peak.color
    if(!is.null(roi.multi)){
      if(is.logical(roi.multi) && length(roi.multi) == 1)
        current.gpar$roi.multi <- roi.multi
      else
        log_message('roi.multi must be either TRUE or FALSE', quote = FALSE)				
    }
    
    ## Change file parameters
    if (!is.null(w1Range)){
      w1Range <- suppressWarnings(as.numeric(w1Range))
      if (any(is.na(w1Range)))
        log_message('Chemical shift values must be numeric', quote = FALSE)
      else{
        if (nDim == 1)
          current.gpar$usr[3:4] <- sort(w1Range)
        else
          current.gpar$usr[3:4] <- sort(w1Range)
      }
    }
    if (!is.null(w2Range)){
      w2Range <- suppressWarnings(as.numeric(w2Range))
      if (any(is.na(w2Range)))
        log_message('Chemical shift values must be numeric', quote = FALSE)
      else
        current.gpar$usr[1:2] <- sort(w2Range)
    }
    
    ## Update the file folder
    fileFolder[[i]]$graphics.par <- current.gpar
    fileFolder[[i]]$file.par <- current.fpar
  }       
  
  ## Update graphics parameters
  if( !is.null(usr) && save.backup )
    myAssign("zoom", fileFolder )
  else
    myAssign("fileFolder", fileFolder, save.backup = save.backup)
  
  ## Refresh plots
  if(refresh.graphics)
    refresh()
}

##Refresh active windows 
refresh <- function(main.plot = TRUE, overlay = TRUE, sub.plot = TRUE, 
                    multi.plot = TRUE, ...){
  
  current <- wc()
  
  ## Refresh the main plot
  if(main.plot){
    
    drawPeptides(...)
    if( globalSettings$roiMain  && !is.null(roiTable) && nrow(roiTable) > 0)
      showRoi()  		
    
    ## Display peaks if appropriate
    pList <- fileFolder[[current]]$peak.list
    if( globalSettings$peak.disp )
      pdisp()       
    
    ## Display 1D projection if appropriate
    if(globalSettings$proj.mode && 
       (fileFolder[[current]]$file.par$number_dimensions > 1))
      proj1D() 
    
    ## Add overlays
    if(exists('overlayList') && overlay && !is.null(overlayList) &&
       fileFolder[[current]]$graphics.par$type != 'persp')
      ol(askUsr=FALSE)		
    
    ## Add menus if not present
    if (.Platform$OS == 'windows' && .Platform$GUI == 'Rgui' &&
        !length(grep('$Graph2', winMenuNames(), fixed=TRUE))){
      #gui('main')
      popupGui('main')
    }
  }
  
  ## Refresh the roi sub plot
  if(length(which(dev.list() == 3)) == 1 && sub.plot)
    rvs() 
  
  ## Refresh the roi sub plot
  if(length(which(dev.list() == 4)) == 1 && multi.plot)
    rvm() 
}


## Internal graphics function findTiles
## Finds 2D NMR sparky tiles above user defined noise threshold
## in.folder  - digestR spectral header file to be searched
## internal - if TRUE, changes are made to in.folder and returned, rather than
##						fileFolder
## returns (invisible) an updated file folder list of tiles with viewable 
findTiles <- function(in.folder=fileFolder[[wc()]], internal=FALSE){
  
  ## Define some local variables
  filePar <- in.folder$file.par
  nDim <- filePar$number_dimensions
  bs <- filePar$block_size
  upShifts <- filePar$block_upfield_ppms
  ms <- filePar$matrix_size
  fileName <- filePar$file.name
  
  ## Find total number of tiles
  if (nDim == 1)
    return(in.folder)
  if (is.null(upShifts)){
    tiles <- ceiling(ms / bs)
    if (nDim == 2){
      tiles[3] <- 1
      bs[3] <- 1
    }
    tiles <- (tiles[1] * tiles[2] * tiles[3]) - 1
    lRow <- bs[2]
    lCol <- bs[1]
    w1Above <- ceiling(ms[2] / bs[2]) 
    n <- bs[1] * bs[2] * bs[3]
  }else{
    tiles <- length(bs$w1) - 1
  }
  zlim <- filePar$noise_est * in.folder$graphics.par$clevel
  endFormat <- filePar$endian
  
  ## Open binary connection
  outTiles <- NULL
  con <- myFile(fileName, "rb")
  seek(con, where=filePar$binary_location, origin="start")
  
  ## Read binary	and find the viewable tiles 
  for (i in 0:tiles){
    if (!is.null(upShifts)){
      lCol <- bs$w1[i + 1]
      n <- bs$w1[i + 1] * bs$w2[i + 1]
    }
    cTile <- matrix(abs(readBin(con, size=4, what='double', endian=endFormat, 
                                n=n)), ncol=lCol) 
    if (any(cTile > zlim)){
      outTiles <- c(outTiles, i)
      if (is.null(upShifts)){
        if (any(cTile[, lCol] > zlim))
          outTiles <- c(outTiles, i + w1Above)
        if (any(cTile[, 1] > zlim))
          outTiles <- c(outTiles, i - w1Above) 			
        if (any(cTile[lRow, ] > zlim))
          outTiles <- c(outTiles, i + 1)
        if (any(cTile[1, ] > zlim))
          outTiles <- c(outTiles, i - 1)
      }
    }
  }
  close(con)
  
  ## Assign tiles to global environment
  tiles <- unique(outTiles[outTiles <= tiles])
  in.folder$graphics.par$tiles <- tiles
  if (internal)
    return(in.folder)
  fileNames <- sapply(fileFolder, function(x) x$file.par$file.name)
  fileFolder[[match(fileName, fileNames)]]$graphics.par$tiles <- tiles
  myAssign("fileFolder", fileFolder, save.backup=FALSE)
  
  invisible(in.folder)
}


################################################################################
##                                                                            ##
##                    Internal plotting functions                             ##
##                                                                            ##
################################################################################

## Internal plotting function pdisp
## Displays current peak list and prints any non NA labels
## col  - color for labels (see points)
## cex  - Character expansion (see points)
## pch  - Integer from 1:26 used to define label type (see points)
## ... - Additional arguments can be passed points
pdisp <- function(col, cex, pch, pos, offset, ...){
  
  ## Define current spectrum
  current <- wc()
  lineCol <- fileFolder[[current]]$graphics.par$fg
  if (length(which(dev.list() == 2)) == 0)
    refresh(sub.plot=FALSE, multi.plot=FALSE)
  inFolder <- fileFolder[[current]]
  nDim <- inFolder$file.par$number_dimensions
  peaklist <- inFolder$peak.list
  labels <- which(peaklist$Assignment != 'NA')
  
  ## Make sure something has been peak picked
  #if(length(peaklist$w2) == 0)
  #	return(invisible())
  
  ## Set peak color
  if( missing(col) )
    col <- inFolder$graphics.par$peak.color
  
  ## Set peak label size
  if (nDim == 1)
    w2Range <- inFolder$file.par$downfield_ppm[1] - 
    inFolder$file.par$upfield_ppm[1]
  else
    w2Range <- inFolder$file.par$downfield_ppm[2] - 
    inFolder$file.par$upfield_ppm[2]
  currRange <- inFolder$graphics.par$usr[1] - inFolder$graphics.par$usr[2]
  zmFactor <- w2Range / currRange
  zmFactor <- zmFactor^(1 / (.9 * zmFactor))
  if (missing(cex)){
    cex <- globalSettings$peak.cex
    cex <- cex * zmFactor
  }
  
  ## Set peak label alignment
  if (missing(pos)){
    if (globalSettings$peak.labelPos == 'top')
      pos <- 3
    else if (globalSettings$peak.labelPos == 'bottom')
      pos <- 1
    else if (globalSettings$peak.labelPos == 'left')
      pos <- 2
    else if (globalSettings$peak.labelPos == 'right')
      pos <- 4
    else
      pos <- NULL
  }
  if (missing(offset))
    offset <- .3
  
  ## Set graphics for peak display    
  if(nDim > 1){
    if( missing(pch) )
      pch <- globalSettings$peak.pch 
    if(is.null(peaklist$Multiplet))
      mainShifts <- NULL
    else
      mainShifts <- which( is.na(peaklist$Multiplet)  )
    if(length(mainShifts) != 0 )
      points(peaklist$w2[mainShifts], peaklist$w1[mainShifts], col = col, 
             cex = cex, pch = 7, ... )
    
    points( peaklist$w2, peaklist$w1, col = col, cex = cex, pch = pch, ... )
    y <- peaklist$w1[labels]
  }else{
    abline( h = inFolder$file.par$noise_est * globalSettings$thresh.1D +
              inFolder$file.par$zero_offset, col=lineCol )
    if( missing(pch) )
      pch = 25
    points(peaklist$w2, peaklist$Height, col = col, pch = pch, cex = cex, ... )
    y <- peaklist$Height[labels]
  }
  
  ## Add assignment labels
  if(length(labels) > 0)
    text(peaklist$w2[labels], y, labels = peaklist$Assignment[labels], 
         cex = cex * .8, col = col, pos = pos, offset = offset)
  
  ## Leave the console active
  bringFocus(-1) 
  
  invisible()	
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
drawNMR <- function ( in.folder = fileFolder[[wc()]], w1Range, w2Range, 
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
  
  ## Plot 1D data, slices, projections
  if( in.folder$file.par$number_dimensions == 1 ){
    #		if (!is.null(in.folder$file.par$block_upfield_ppms))
    #			plotRsd1D( in.folder = in.folder, xlab = xlab, ylab = ylab, main = main, 
    #					add = add, axes = axes, offset = offset, ...)
    #		else
    plot1D( in.folder = in.folder, xlab = xlab, ylab = ylab, main = main, 
            add = add, axes = axes, offset = offset, ...)
  }else{
    
    ## Set axis labels
    if( is.null(xlab) )
      xlab <- paste(in.folder$file.par$nucleus[2])#, 'PPM', sep=' ')
    if( is.null(ylab) )
      ylab <- paste(in.folder$file.par$nucleus[1])# , 'PPM', sep=' ')
    
    ## Plot data in the multiple file/subplot windows, or data in memory
    if (in.folder$graphics.par$type != 'persp' && 
        (p.window != 'main' || (!is.null(in.folder$data) && 
                                !is.null(in.folder$w1) && !is.null(in.folder$w2)))){
      draw2D( in.folder = in.folder, pos.zlim = pos.zlim, neg.zlim = neg.zlim,
              xlab = xlab, ylab = ylab, main = main, add = add, axes = axes, ...)
      
      
    }else if (!is.null(in.folder$file.par$block_upfield_ppms)){
      
      ## Plot RSD file
      if (in.folder$graphics.par$type != 'persp')
        plotRsd2D( in.folder = in.folder, pos.zlim = pos.zlim, 
                   neg.zlim = neg.zlim, xlab = xlab, ylab = ylab, main = main, 
                   add = add, axes = axes, ...)
      else
        perspRsd(in.folder = in.folder, xlab = xlab, ylab = ylab, main = main, 
                 ...) 
      
    }else if (in.folder$file.par$number_dimensions == 2){
      
      ## Plot 2D UCSF file
      if (in.folder$graphics.par$type != 'persp')
        plot2D( in.folder = in.folder, pos.zlim = pos.zlim, 
                neg.zlim = neg.zlim, xlab = xlab, ylab = ylab, main = main, 
                add = add, axes = axes, ...)
      else
        persp2D(in.folder = in.folder, xlab = xlab, ylab = ylab, main = main, 
                ...) 
      
    }else if (in.folder$file.par$number_dimensions == 3){
      
      ## Plot 3D UCSF file
      if (in.folder$graphics.par$type != 'persp')
        plot3D( in.folder = in.folder, pos.zlim = pos.zlim, 
                neg.zlim = neg.zlim, xlab = xlab, ylab = ylab, main = main, 
                add = add, axes = axes, ...)
      else
        persp2D(in.folder = in.folder, xlab = xlab, ylab = ylab, main = main, 
                ...) 
    }
  }
  
  bringFocus(-1) #return focus to console
}


##Internal graphics wrapper function plot2D
##Draws a 2D NMR spectrum from a binary connection
##in.folder - Header of the file to be plotted, default is current spectrum
##w1Range - Chemical shift range in the indirect dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)          
##w2Range - Chemical shift range in the direct dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)
##pos.zlim - Min and max poisitive intensities to be displayed, default is the 
##          most recent setting used with the file, the format is c(lower,upper)
##neg.zlim - Max and min negitive intensities to be displayed, default is the 
##          most recent setting used with the file, the format is c(lower,upper)
##type  - specifies the type of plot to be generated: 'image' is the fastest, 
##        'contour' produces a contour plot, 'filled' is a filled contour plot 
##         and 'persp' produces a 3D perspective plot. Default is the most 
##         recent type used with the file
##pos.color - color of positive contours, default is the most recent setting
##         for the file, see colors() for the many color options
##neg.color - color of negative contours, default is the most recent setting
##         for the file, see colors() for the many color options
##nlevels - the number of contour intervals to be drawn, the default is the most
##         recent setting used for the file
##xlab    - x axis label, default is the direct detected nucleus in PPM
##ylab    - y axis label, default is the indirect detected nucleus in PPM
##main    - main title for plot, default is the file name
##conDisp    - logical vector, c(TRUE, TRUE) plots positive and negative 
##           contours, c(TRUE, FALSE) plots only positive, c(FALSE, TRUE) plots
##           only the negative contours, c(FALSE, FALSE) plots no contours
##add       - logical argument, TRUE adds new data to an existing plot, FALSE 
##            generates a new plot
##p.window  -  The window to be used, can be 'main', 'sub', 'multi', or 'stats'
##axes      - logical argument, TRUE draws axes
##...       - Additinal graphics paramaters can be passed to par()
##returns   - a plot of a 2D NMR spectrum  
plot2D <- function(     
    in.folder = fileFolder[[wc()]],
    w1Range = in.folder$graphics.par$usr[3:4],
    w2Range = in.folder$graphics.par$usr[1:2],
    pos.zlim = c(in.folder$file.par$noise_est * in.folder$graphics.par$clevel,
                 in.folder$file.par$noise_est * in.folder$graphics.par$clevel *
                   in.folder$graphics.par$nlevels),
    neg.zlim = -(rev(pos.zlim)),
    type = in.folder$graphics.par$type,
    pos.color = in.folder$graphics.par$pos.color,
    neg.color = in.folder$graphics.par$neg.color,
    nlevels = in.folder$graphics.par$nlevels,
    conDisp = in.folder$graphics.par$conDisp,
    xlab = paste(in.folder$file.par$nucleus[2]), 
    ylab = paste(in.folder$file.par$nucleus[1]), 
    main = in.folder$file.par$user_title, add = FALSE, axes = TRUE, ...){
  
  ## Find total usable tiles if none exists
  if(is.null(in.folder$graphics.par$tiles))
    in.folder <- findTiles(in.folder = in.folder)
  
  ## Define some local variables
  bs <- in.folder$file.par$block_size
  ms <- in.folder$file.par$matrix_size
  uf <- in.folder$file.par$upfield_ppm
  df <- in.folder$file.par$downfield_ppm
  endFormat <- in.folder$file.par$endian
  binLoc <- in.folder$file.par$binary_location
  
  ## Find best chemical shift match
  in.folder$w1 <- seq(uf[1], df[1], length.out = ms[1])
  in.folder$w2 <- seq(uf[2], df[2],	length.out = ms[2])	
  w1Range <- sort(w1Range); w2Range <- sort(w2Range)
  out <- NULL
  t1 <- findInterval(w1Range, in.folder$w1, all.inside = TRUE)
  d1 <- findInterval(w2Range, in.folder$w2, all.inside = TRUE)
  t2 <- t1 + 1; d2 <- d1 + 1
  for(i in 1:2){
    out$w1[i] <- switch( which.min(c(
      abs(w1Range[i] - in.folder$w1[t1[i]]),
      abs(w1Range[i] - in.folder$w1[t2[i]]))), t1[i], t2[i])
    out$w2[i] <- switch( which.min(c(
      abs(w2Range[i] - in.folder$w2[d1[i]]),
      abs(w2Range[i] - in.folder$w2[d2[i]]))), d1[i], d2[i])
  }
  
  ## Invert w1/w2 selection to match binary data format
  out$w1 <- sort(ms[1] - out$w1) + 1
  out$w2 <- sort(ms[2] - out$w2) + 1
  
  ## Find sparky tiles that will be plotted
  w1Tiles <- (ceiling(out$w1[1] / bs[1]): ceiling(out$w1[2] / bs[1]))
  w2Tiles <- (ceiling(out$w2[1] / bs[2]): ceiling(out$w2[2] / bs[2]))
  tiles <- NULL
  for ( i in 1:length(w1Tiles))
    tiles <- c(tiles, ((w1Tiles[i]-1) * (ceiling(ms[2] / bs[2])) + (w2Tiles-1)))
  
  if(!add){
    
    ## Set new axes and 'clear' old plot with a box
    plot(0, 0, axes=FALSE, type = 'n', xlab=xlab, ylab=ylab, main=main, 
         xaxs='i', yaxs='i', cex.main=in.folder$graphics.par$cex.main)
    par(usr = c(w2Range[2:1], w1Range[2:1]))
    rect(c(w2Range[2], w2Range[2], w2Range[2], uf[2] ),  
         c(w1Range[2], uf[1], w1Range[1], w1Range[1] ),
         c(w2Range[1], w2Range[1], df[2], w2Range[1] ),
         c(df[1], w1Range[1], w1Range[2], w1Range[2]), 
         col = "grey", border = NA)
    
    ## Draw plot labels
    box(col=in.folder$graphics.par$fg, bty='o')
    if(axes){
      par(usr = c(rev(par('usr')[1:2]), rev(par('usr')[3:4]))) 
      xlabvals <- pretty(w2Range, 10)    
      xlabvals <- c(xlabvals, par('usr')[2])
      n1 <- length(xlabvals)
      if (!is.na(in.folder$graphics.par$xtck) && 
          in.folder$graphics.par$xtck == 1)
        xlty <- 2
      else
        xlty <- 1
      if (!is.na(in.folder$graphics.par$ytck) && 
          in.folder$graphics.par$ytck == 1)
        ylty <- 2
      else
        ylty <- 1
      axis(side=1, at=min(w2Range) + (max(w2Range)-xlabvals),  
           labels = c( xlabvals[-n1], paste('    ', xlab)), lty=xlty,
           tck=in.folder$graphics.par$xtck, 
           cex.axis=in.folder$graphics.par$cex.axis)   
      ylabvals <- pretty(w1Range, 10)
      ylabvals <- c(ylabvals, par('usr')[4] )
      n1 <- length(ylabvals)
      axis(side=2, at=min(w1Range) + (max(w1Range) - ylabvals),  
           labels = c( ylabvals[-n1], paste('    ', ylab)), lty=ylty, 
           tck=in.folder$graphics.par$ytck, 
           cex.axis=in.folder$graphics.par$cex.axis) 
      par(usr = c(rev(par('usr')[1:2]), rev(par('usr')[3:4]))) 
    }
  }
  
  ## Set plot function for auto mode
  if (!type %in% c('auto', 'image', 'contour', 'filled'))
    type <- 'auto'
  if(type == 'auto'){
    if((length(w1Tiles) + length(w2Tiles)) < 17)
      type <- 'contour'
    else
      type <- 'image'
  }
  
  ## Set plotting function
  plotCur <- switch( which(c('image', 'contour', 'filled') == type),
                     function(...){ image(add = TRUE, ...)},
                     function(zlim, ...){ contour(zlim, add = TRUE, nlevels = nlevels, 
                                                  levels = seq(zlim[1], zlim[2], length.out = nlevels),  
                                                  drawlabels = FALSE, ...)},
                     function(zlim, ...){ fillCon(
                       levels = seq(zlim[1], zlim[2], length.out = nlevels), ...)}
  )
  
  
  in.folder$w1 <- 1:((ms[1] %/% bs[1] + 1) * bs[1]) * 
    ((df[1] - uf[1]) / (ms[1] - 1))
  in.folder$w1 <- in.folder$w1 + (df[1] - rev(in.folder$w1)[1] ) 
  in.folder$w2 <- 1:((ms[2] %/% bs[2] + 1) * bs[2]) * 
    ((df[2] - uf[2]) / (ms[2] - 1))
  in.folder$w2 <- in.folder$w2 + (df[2] - rev(in.folder$w2)[1] ) 
  
  ## Open connection to binary file
  con <- myFile(in.folder$file.par$file.name, "rb")
  outData <- NULL
  
  ##Read binary data for each tile
  w2TileNum <- ceiling(ms[2] / bs[2] )
  if(type != 'image'){
    w2Slice <- rep(NA, length(in.folder$w2))
    w1Slice <- rep(NA, length(in.folder$w1))
  } 
  
  for(i in tiles){
    ## define current w1/w2 range
    w1 <- ceiling( (i + 1) / w2TileNum )
    w1B <- (1:bs[1]) + bs[1] * (w1 - 1) 
    w2 <- (i + 1) - (w1 * w2TileNum ) + w2TileNum
    w2B <- (1:bs[2]) + bs[2] * (w2 - 1) 
    
    ##skip tiles with no data above the threshold
    if( length( which(in.folder$graphics.par$tiles == (i))) == 0 ){
      if(type != 'image'){
        w2Slice[w2B] <- rep(NA, length(w2B))
        w1Slice[w1B] <- rep(NA, length(w1B))
      }
      next() 
    }
    
    ## read binary
    w1 <- rev(in.folder$w1)[w1B]
    w2 <- rev(in.folder$w2)[w2B]
    w1 <- sort(w1); w2=sort(w2)
    seek(con, bs[1] * bs[2] * 4 * i + binLoc, origin = 'start')
    outData <- matrix(rev(readBin(con, size=4, what='double',
                                  endian = endFormat, n=(bs[1] * bs[2]))), ncol=bs[1])
    
    ## Overlap tile with previous tile to correct edge effects
    if(type != 'image'){   
      edge <- in.folder$w1[(which(in.folder$w1 == rev(w1)[1])) + 1]
      if(length(edge) != 0 && !is.na(edge)){    
        w1 <- c(w1, edge)
        outData <- cbind(outData, w2Slice[w2B])
      }  
      w2Slice[w2B] <- outData[, 1]
      
      edge <- in.folder$w2[(which(in.folder$w2 == rev(w2)[1])) + 1]
      if(length(edge) != 0 && !is.na(edge)){    
        w2 <- c(w2, edge)
        outData <- suppressWarnings(rbind(outData, w1Slice[w1B]))
      }  
      w1Slice[w1B] <- outData[1, (1:length(w1B))]             
    }
    
    ## Plot the current tile
    if(conDisp[1])
      plotCur(x = w2, y = w1, z = outData, zlim = pos.zlim, 
              col = pos.color, ...)
    if(conDisp[2])
      plotCur(x = w2, y = w1, z= outData, zlim = neg.zlim, col = neg.color, ...)
  }
  
  ## Close binary conection
  closeAllConnections()
}


##Internal graphics wrapper function plot2D
##Draws a single 2D slice from a 3D NMR spectrum from binary connection
##in.folder - Header of the file to be plotted, default is current spectrum
##w1Range - Chemical shift range in the indirect dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)          
##w2Range - Chemical shift range in the direct dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)
##pos.zlim - Min and max poisitive intensities to be displayed, default is the 
##          most recent setting used with the file, the format is c(lower,upper)
##neg.zlim - Max and min negitive intensities to be displayed, default is the 
##          most recent setting used with the file, the format is c(lower,upper)
##type  - specifies the type of plot to be generated: 'image' is the fastest, 
##        'contour' produces a contour plot, 'filled' is a filled contour plot 
##         and 'persp' produces a 3D perspective plot. Default is the most 
##         recent type used with the file
##pos.color - color of positive contours, default is the most recent setting
##         for the file, see colors() for the many color options
##neg.color - color of negative contours, default is the most recent setting
##         for the file, see colors() for the many color options
##nlevels - the number of contour intervals to be drawn, the default is the most
##         recent setting used for the file
##xlab    - x axis label, default is the direct detected nucleus in PPM
##ylab    - y axis label, default is the indirect detected nucleus in PPM
##main    - main title for plot, default is the file name
##conDisp    - logical vector, c(TRUE, TRUE) plots positive and negative 
##           contours, c(TRUE, FALSE) plots only positive, c(FALSE, TRUE) plots
##           only the negative contours, c(FALSE, FALSE) plots no contours
##add       - logical argument, TRUE adds new data to an existing plot, FALSE 
##            generates a new plot
##p.window  -  The window to be used, can be 'main', 'sub', 'multi', or 'stats'
##axes      - logical argument, TRUE draws axes
##...       - Additinal graphics paramaters can be passed to par()
##returns   - a plot of a 2D NMR spectrum  
plot3D <- function(in.folder=fileFolder[[wc()]],
                   w1Range=in.folder$graphics.par$usr[3:4],
                   w2Range=in.folder$graphics.par$usr[1:2],
                   pos.zlim=c(in.folder$file.par$noise_est * in.folder$graphics.par$clevel,
                              in.folder$file.par$noise_est * in.folder$graphics.par$clevel *
                                in.folder$graphics.par$nlevels),
                   neg.zlim=-(rev(pos.zlim)), type=in.folder$graphics.par$type, 
                   pos.color=in.folder$graphics.par$pos.color,
                   neg.color= in.folder$graphics.par$neg.color,
                   nlevels=in.folder$graphics.par$nlevels, 
                   conDisp=in.folder$graphics.par$conDisp,
                   xlab=paste(in.folder$file.par$nucleus[2]), 
                   ylab=paste(in.folder$file.par$nucleus[1]), 
                   main=in.folder$file.par$user_title, add=FALSE, axes=TRUE, ...){
  
  ## Find total usable tiles if none exists
  if (is.null(in.folder$graphics.par$tiles))
    in.folder <- findTiles(in.folder=in.folder)
  
  ## Define some local variables
  filePar <- in.folder$file.par
  graphicsPar <- in.folder$graphics.par
  bs <- filePar$block_size
  ms <- filePar$matrix_size
  uf <- filePar$upfield_ppm
  df <- filePar$downfield_ppm
  endFormat <- filePar$endian
  binLoc <- filePar$binary_location
  
  ## Find best chemical shift match
  in.folder$w1 <- seq(uf[1], df[1], length.out=ms[1])
  in.folder$w2 <- seq(uf[2], df[2],	length.out=ms[2])	
  in.folder$w3 <- seq(uf[3], df[3],	length.out=ms[3])	
  w1Range <- sort(w1Range); w2Range <- sort(w2Range)
  w3Range <- rep(filePar$z_value, 2)
  t1 <- findInterval(w1Range, in.folder$w1, all.inside=TRUE)
  d1 <- findInterval(w2Range, in.folder$w2, all.inside=TRUE)
  z1 <- findInterval(w3Range, in.folder$w3, all.inside=TRUE)
  t2 <- t1 + 1; d2 <- d1 + 1; z2 <- z1 + 1
  out <- NULL
  for (i in 1:2){
    out$w1[i] <- switch(which.min(c(abs(w1Range[i] - in.folder$w1[t1[i]]),
                                    abs(w1Range[i] - in.folder$w1[t2[i]]))), t1[i], t2[i])
    out$w2[i] <- switch(which.min(c(abs(w2Range[i] - in.folder$w2[d1[i]]),
                                    abs(w2Range[i] - in.folder$w2[d2[i]]))), d1[i], d2[i])
    out$w3[i] <- switch(which.min(c(abs(w3Range[i] - in.folder$w3[z1[i]]),
                                    abs(w3Range[i] - in.folder$w3[z2[i]]))), z1[i], z2[i])
  }
  
  ## Invert w1/w2 selection to match binary data format
  out$w1 <- sort(ms[1] - out$w1) + 1
  out$w2 <- sort(ms[2] - out$w2) + 1
  
  ## Find sparky tiles that will be plotted
  w1Tiles <- (ceiling(out$w1[1] / bs[1]):ceiling(out$w1[2] / bs[1]))
  w2Tiles <- (ceiling(out$w2[1] / bs[2]):ceiling(out$w2[2] / bs[2]))
  w3Tiles <- (ceiling(out$w3[1] / bs[3]):ceiling(out$w3[2] / bs[3]))
  numTiles <- ceiling(ms / bs)
  tiles <- NULL
  for (i in w1Tiles)
    tiles <- c(tiles, (i - 1) * numTiles[2] + (w2Tiles - 1))
  tiles <- tiles + (w3Tiles - 1) * numTiles[1] * numTiles[2]
  
  if (!add){
    
    ## Set new axes and 'clear' old plot with a box
    plot(0, 0, axes=FALSE, type='n', xlab=xlab, ylab=ylab, main=main,	xaxs='i', 
         yaxs='i', cex.main=graphicsPar$cex.main)
    par(usr=c(w2Range[2:1], w1Range[2:1]))
    rect(c(w2Range[2], w2Range[2], w2Range[2], uf[2]), 
         c(w1Range[2], uf[1], w1Range[1], w1Range[1] ),
         c(w2Range[1], w2Range[1], df[2], w2Range[1]), 
         c(df[1], w1Range[1], w1Range[2], w1Range[2]), col="grey", border=NA)
    
    ## Draw plot labels
    box(col=graphicsPar$fg, bty='o')
    if (axes){
      par(usr=c(rev(par('usr')[1:2]), rev(par('usr')[3:4]))) 
      xlabvals <- pretty(w2Range, 10)    
      xlabvals <- c(xlabvals, par('usr')[2])
      n1 <- length(xlabvals)
      if (!is.na(graphicsPar$xtck) && 
          graphicsPar$xtck == 1)
        xlty <- 2
      else
        xlty <- 1
      if (!is.na(graphicsPar$ytck) && 
          graphicsPar$ytck == 1)
        ylty <- 2
      else
        ylty <- 1
      axis(side=1, at=min(w2Range) + (max(w2Range) - xlabvals),  
           labels=c( xlabvals[-n1], paste('    ', xlab)), lty=xlty,
           tck=graphicsPar$xtck, 
           cex.axis=graphicsPar$cex.axis)   
      ylabvals <- pretty(w1Range, 10)
      ylabvals <- c(ylabvals, par('usr')[4] )
      n1 <- length(ylabvals)
      axis(side=2, at=min(w1Range) + (max(w1Range) - ylabvals),  
           labels = c( ylabvals[-n1], paste('    ', ylab)), lty=ylty, 
           tck=graphicsPar$ytck, 
           cex.axis=graphicsPar$cex.axis) 
      par(usr=c(rev(par('usr')[1:2]), rev(par('usr')[3:4]))) 
    }
  }
  
  ## Set plot function for auto mode
  if (!type %in% c('auto', 'image', 'contour', 'filled'))
    type <- 'auto'
  if (type == 'auto'){
    if ((length(w1Tiles) + length(w2Tiles)) < 17)
      type <- 'contour'
    else
      type <- 'image'
  }
  
  ## Set plotting function
  plotCur <- switch(which(c('image', 'contour', 'filled') == type),
                    function(...) image(add=TRUE, ...),
                    function(zlim, ...) contour(zlim, add=TRUE, nlevels=nlevels, 
                                                levels=seq(zlim[1], zlim[2], length.out=nlevels), drawlabels=FALSE, 
                                                ...),
                    function(zlim, ...) fillCon(levels=seq(zlim[1], zlim[2], 
                                                           length.out=nlevels), ...))
  
  
  in.folder$w1 <- 1:((ms[1] %/% bs[1] + 1) * bs[1]) * 
    ((df[1] - uf[1]) / (ms[1] - 1))
  in.folder$w1 <- in.folder$w1 + (df[1] - rev(in.folder$w1)[1]) 
  in.folder$w2 <- 1:((ms[2] %/% bs[2] + 1) * bs[2]) * 
    ((df[2] - uf[2]) / (ms[2] - 1))
  in.folder$w2 <- in.folder$w2 + (df[2] - rev(in.folder$w2)[1]) 
  
  ## Read binary data for each tile
  outData <- NULL
  tileSize <- bs[1] * bs[2] * bs[3] * 4
  w2TileNum <- ceiling(ms[2] / bs[2])
  if (type != 'image'){
    w2Slice <- rep(NA, length(in.folder$w2))
    w1Slice <- rep(NA, length(in.folder$w1))
  } 
  con <- myFile(filePar$file.name, "rb")
  for (i in tiles){
    
    ## Define current w1/w2 range
    tile2D <- i - (w3Tiles - 1) * numTiles[1] * numTiles[2]
    w1 <- ceiling((tile2D + 1) / w2TileNum)
    w1B <- (1:bs[1]) + bs[1] * (w1 - 1) 
    w2 <- (tile2D + 1) - (w1 * w2TileNum ) + w2TileNum
    w2B <- (1:bs[2]) + bs[2] * (w2 - 1) 
    
    ## Skip tiles with no data above the threshold
    if (!i %in% graphicsPar$tiles){
      if (type != 'image'){
        w2Slice[w2B] <- rep(NA, length(w2B))
        w1Slice[w1B] <- rep(NA, length(w1B))
      }
      next() 
    }
    
    ## Define data location for current tile
    tileLoc <- tileSize * i + binLoc
    if (bs[3] > 1)
      zPos <- (out$w3[1] - 1) * bs[3]
    else
      zPos <- 0
    dataLoc <- tileLoc + bs[1] * bs[2] * 4 * zPos
    
    ## Read binary
    seek(con, dataLoc, origin='start')
    outData <- matrix(rev(readBin(con, size=4, what='double', endian=endFormat, 
                                  n=bs[1] * bs[2])), ncol=bs[1])
    
    ## Overlap tile with previous tile to correct edge effects
    w1 <- rev(in.folder$w1)[w1B]
    w2 <- rev(in.folder$w2)[w2B]
    w1 <- sort(w1); w2=sort(w2)
    if (type != 'image'){   
      edge <- in.folder$w1[(which(in.folder$w1 == rev(w1)[1])) + 1]
      if (length(edge) != 0 && !is.na(edge)){    
        w1 <- c(w1, edge)
        outData <- cbind(outData, w2Slice[w2B])
      }  
      w2Slice[w2B] <- outData[, 1]
      edge <- in.folder$w2[(which(in.folder$w2 == rev(w2)[1])) + 1]
      if (length(edge) != 0 && !is.na(edge)){    
        w2 <- c(w2, edge)
        outData <- suppressWarnings(rbind(outData, w1Slice[w1B]))
      }  
      w1Slice[w1B] <- outData[1, (1:length(w1B))]             
    }
    
    ## Plot the current tile
    if (conDisp[1])
      plotCur(x=w2, y=w1, z=outData, zlim=pos.zlim, col=pos.color, ...)
    if (conDisp[2])
      plotCur(x=w2, y=w1, z=outData, zlim=neg.zlim, col=neg.color, ...)
  }
  
  ## Close binary conection
  close(con)
}


##Internal graphics wrapper function plotRsd2D
##Draws a 2D RSD NMR spectrum from a binary connection
##in.folder - Header of the file to be plotted, default is current spectrum
##w1Range - Chemical shift range in the indirect dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)          
##w2Range - Chemical shift range in the direct dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)
##pos.zlim - Min and max poisitive intensities to be displayed, default is the 
##          most recent setting used with the file, the format is c(lower,upper)
##neg.zlim - Max and min negitive intensities to be displayed, default is the 
##          most recent setting used with the file, the format is c(lower,upper)
##type  - specifies the type of plot to be generated: 'image' is the fastest, 
##        'contour' produces a contour plot, 'filled' is a filled contour plot 
##         and 'persp' produces a 3D perspective plot. Default is the most 
##         recent type used with the file
##pos.color - color of positive contours, default is the most recent setting
##         for the file, see colors() for the many color options
##neg.color - color of negative contours, default is the most recent setting
##         for the file, see colors() for the many color options
##nlevels - the number of contour intervals to be drawn, the default is the most
##         recent setting used for the file
##xlab    - x axis label, default is the direct detected nucleus in PPM
##ylab    - y axis label, default is the indirect detected nucleus in PPM
##main    - main title for plot, default is the file name
##conDisp    - logical vector, c(TRUE, TRUE) plots positive and negative 
##           contours, c(TRUE, FALSE) plots only positive, c(FALSE, TRUE) plots
##           only the negative contours, c(FALSE, FALSE) plots no contours
##add       - logical argument, TRUE adds new data to an existing plot, FALSE 
##            generates a new plot
##p.window  -  The window to be used, can be 'main', 'sub', 'multi', or 'stats'
##axes      - logical argument, TRUE draws axes
##...       - Additinal graphics paramaters can be passed to par()
##returns   - a plot of a 2D NMR spectrum  
plotRsd2D <- function(in.folder=fileFolder[[wc()]], 
                      w1Range=in.folder$graphics.par$usr[3:4],
                      w2Range=in.folder$graphics.par$usr[1:2],
                      pos.zlim=c(in.folder$file.par$noise_est * in.folder$graphics.par$clevel,
                                 in.folder$file.par$noise_est * in.folder$graphics.par$clevel *
                                   in.folder$graphics.par$nlevels),
                      neg.zlim=-(rev(pos.zlim)), type=in.folder$graphics.par$type,
                      pos.color=in.folder$graphics.par$pos.color,
                      neg.color=in.folder$graphics.par$neg.color,
                      nlevels=in.folder$graphics.par$nlevels,
                      conDisp=in.folder$graphics.par$conDisp,
                      xlab=paste(in.folder$file.par$nucleus[2]), 
                      ylab=paste(in.folder$file.par$nucleus[1]), 
                      main=in.folder$file.par$user_title, add=FALSE, axes=TRUE, ...){
  
  ## Find total usable tiles if none exists
  if (is.null(in.folder$graphics.par$tiles))
    in.folder <- findTiles(in.folder=in.folder)
  
  ## Define some local variables
  filePar <- in.folder$file.par
  graphicsPar <- in.folder$graphics.par
  uf <- filePar$upfield_ppm
  df <- filePar$downfield_ppm
  binLoc <- filePar$binary_location
  
  ## Find best w1/w2 matches	
  w1Range <- sort(w1Range)
  w2Range <- sort(w2Range)
  in.folder$w1 <- seq(filePar$upfield_ppm[1], filePar$downfield_ppm[1], 
                      length.out=filePar$matrix_size[1])
  in.folder$w2 <- seq(filePar$upfield_ppm[2], filePar$downfield_ppm[2], 
                      length.out=filePar$matrix_size[2])
  t1 <- findInterval(w1Range, in.folder$w1, all.inside=TRUE)
  d1 <- findInterval(w2Range, in.folder$w2, all.inside=TRUE)
  t2 <- t1 + 1
  d2 <- d1 + 1
  out <- list(w1=NULL, w2=NULL)
  for (i in 1:2){
    out$w1[i] <- switch(which.min(c(abs(w1Range[i] - in.folder$w1[t1[i]]),
                                    abs(w1Range[i] - in.folder$w1[t2[i]]))), t1[i], t2[i])
    out$w2[i] <- switch(which.min(c(abs(w2Range[i] - in.folder$w2[d1[i]]),
                                    abs(w2Range[i] - in.folder$w2[d2[i]]))), d1[i], d2[i])
  }
  w1Range <- c(in.folder$w1[(out$w1[1])], in.folder$w1[(out$w1[2])])
  w2Range <- c(in.folder$w2[(out$w2[1])], in.folder$w2[(out$w2[2])])
  winW1 <- in.folder$w1[(out$w1[1]:out$w1[2])]
  winW2 <- in.folder$w2[(out$w2[1]:out$w2[2])]
  
  ## Find tiles that will be plotted
  tiles <- NULL
  upShifts <- filePar$block_upfield_ppms
  downShifts <- filePar$block_downfield_ppms
  blockSizes <- filePar$block_size
  for (tNum in seq_along(filePar$block_size$w1)){
    
    ## Get the chemical shifts for the current tile
    blockW1 <- seq(upShifts$w1[tNum], downShifts$w1[tNum], 
                   length.out=blockSizes$w1[tNum])
    blockW2 <- seq(upShifts$w2[tNum], downShifts$w2[tNum],
                   length.out=blockSizes$w2[tNum])
    
    ## Check the window for the presence of any shift in the current block
    if (any(round(blockW1, 3) %in% round(winW1, 3)) && 
        any(round(blockW2, 3) %in% round(winW2, 3)))
      tiles <- c(tiles, tNum)
  }
  
  if (!add){
    
    ## Set new axes and 'clear' old plot with a box
    plot(0, 0, axes=FALSE, type='n', xlab=xlab, ylab=ylab, main=main, 
         xaxs='i', yaxs='i', cex.main=graphicsPar$cex.main)
    par(usr=c(w2Range[2:1], w1Range[2:1]))
    rect(c(w2Range[2], w2Range[2], w2Range[2], uf[2]),  
         c(w1Range[2], uf[1], w1Range[1], w1Range[1]),
         c(w2Range[1], w2Range[1], df[2], w2Range[1]),
         c(df[1], w1Range[1], w1Range[2], w1Range[2]), 
         col="grey", border=NA)
    
    ## Draw plot labels
    box(col=graphicsPar$fg, bty='o')
    if (axes){
      par(usr=c(rev(par('usr')[1:2]), rev(par('usr')[3:4]))) 
      xlabvals <- pretty(w2Range, 10)    
      xlabvals <- c(xlabvals, par('usr')[2])
      n1 <- length(xlabvals)
      if (!is.na(graphicsPar$xtck) && graphicsPar$xtck == 1)
        xlty <- 2
      else
        xlty <- 1
      if (!is.na(graphicsPar$ytck) && graphicsPar$ytck == 1)
        ylty <- 2
      else
        ylty <- 1
      axis(side=1, at=min(w2Range) + (max(w2Range) - xlabvals),  
           labels=c(xlabvals[-n1], paste('    ', xlab)), lty=xlty,
           tck=graphicsPar$xtck, cex.axis=graphicsPar$cex.axis)   
      ylabvals <- pretty(w1Range, 10)
      ylabvals <- c(ylabvals, par('usr')[4])
      n1 <- length(ylabvals)
      axis(side=2, at=min(w1Range) + (max(w1Range) - ylabvals),  
           labels = c( ylabvals[-n1], paste('    ', ylab)), lty=ylty, 
           tck=graphicsPar$ytck, cex.axis=graphicsPar$cex.axis) 
      par(usr=c(rev(par('usr')[1:2]), rev(par('usr')[3:4]))) 
    }
  }
  
  ## Set plot function for auto mode
  if (!type %in% c('auto', 'image', 'contour', 'filled'))
    type <- 'auto'
  if (type == 'auto'){
    if (length(tiles) < 64)
      type <- 'contour'
    else
      type <- 'image'
  }
  
  ## Set plotting function
  plotCur <- switch(which(c('image', 'contour', 'filled') == type),
                    function(...) image(add=TRUE, ...),
                    function(zlim, ...) contour(zlim, add=TRUE, nlevels=nlevels, 
                                                levels=seq(zlim[1], zlim[2], length.out=nlevels),  
                                                drawlabels=FALSE, ...),
                    function(zlim, ...) fillCon(levels=seq(zlim[1], zlim[2], 
                                                           length.out=nlevels), ...))
  
  ## Define data locations for each block
  blockLoc <- binLoc
  for (i in seq_along(filePar$block_size$w1))
    blockLoc <- c(blockLoc, blockLoc[i] + 4 * filePar$block_size$w1[i] * 
                    filePar$block_size$w2[i])
  
  ## Read binary data for each tile
  con <- myFile(filePar$file.name, "rb")
  outData <- NULL
  for (tNum in tiles){
    
    ## Skip tiles with no data above the threshold
    if (!(tNum - 1) %in% graphicsPar$tiles)
      next()
    
    ## Define current w1/w2 range
    w1 <- seq(upShifts$w1[tNum], downShifts$w1[tNum], 
              length.out=filePar$block_size$w1[tNum])
    w2 <- seq(upShifts$w2[tNum], downShifts$w2[tNum], 
              length.out=filePar$block_size$w2[tNum])
    
    ## Read data
    seek(con, blockLoc[tNum], origin='start')
    outData <- matrix(rev(readBin(con, size=4, what='double', 
                                  endian=filePar$endian, 
                                  n=filePar$block_size$w1[tNum] * filePar$block_size$w2[tNum])), 
                      ncol=filePar$block_size$w1[tNum])
    
    ## Plot the current tile
    if (conDisp[1])
      plotCur(x=w2, y=w1, z=outData, zlim=pos.zlim, col=pos.color, ...)
    if (conDisp[2])
      plotCur(x=w2, y=w1, z=outData, zlim=neg.zlim, col=neg.color, ...)
  }
  
  ## Close binary conection
  close(con)
}


##Internal graphics wrapper function draw2D
##Draws a 2D NMR spectrum from a binary connection
##note:    this replaces plot2D for the roi subplots and multiple file windows
##         for batch operations, this function is more efficient.       
##in.folder - Header of the file to be plotted, default is current spectrum
##w1Range - Chemical shift range in the indirect dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)          
##w2Range - Chemical shift range in the direct dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)
##pos.zlim - Min and max poisitive intensities to be displayed, default is the 
##          most recent setting used with the file, the format is c(lower,upper)
##neg.zlim - Max and min negitive intensities to be displayed, default is the 
##          most recent setting used with the file, the format is c(lower,upper)
##type  - specifies the type of plot to be generated: 'image' is the fastest, 
##        'contour' produces a contour plot, 'filled' is a filled contour plot 
##         and 'persp' produces a 3D perspective plot. Default is the most 
##         recent type used with the file
##pos.color - color of positive contours, default is the most recent setting
##         for the file, see colors() for the many color options
##neg.color - color of negative contours, default is the most recent setting
##         for the file, see colors() for the many color options
##nlevels - the number of contour intervals to be drawn, the default is the most
##         recent setting used for the file
##xlab    - x axis label, default is the direct detected nucleus in PPM
##ylab    - y axis label, default is the indirect detected nucleus in PPM
##main    - main title for plot, default is the spectrum name
##conDisp    - logical vector, c(TRUE, TRUE) plots positive and negative 
##           contours, c(TRUE, FALSE) plots only positive, c(FALSE, TRUE) plots
##           only the negative contours, c(FALSE, FALSE) plots no contours
##add       - logical argument, TRUE adds new data to an existing plot, FALSE 
##            generates a new plot
##axes      - logical argument, TRUE makes pretty labels
##roiMax    - logical argument, TRUE plots a point on the maximum
##            visible signal in the window
##...       - Additinal graphics paramaters can be passed to par()
##returns   - a plot of a 2D NMR spectrum  
draw2D <- function (
    in.folder = fileFolder[[wc()]],
    w1Range = in.folder$graphics.par$usr[3:4],
    w2Range = in.folder$graphics.par$usr[1:2],
    pos.zlim = c(in.folder$file.par$noise_est * in.folder$graphics.par$clevel,
                 in.folder$file.par$noise_est * in.folder$graphics.par$clevel *
                   in.folder$graphics.par$nlevels),
    neg.zlim = -(rev(pos.zlim)),
    type = in.folder$graphics.par$type,
    pos.color = in.folder$graphics.par$pos.color,
    neg.color = in.folder$graphics.par$neg.color,
    nlevels = in.folder$graphics.par$nlevels,
    conDisp = in.folder$graphics.par$conDisp,
    xlab = paste(in.folder$file.par$nucleus[2], 'PPM', sep=' '),
    ylab = paste(in.folder$file.par$nucleus[1], 'PPM', sep=' '),
    main = in.folder$file.par$user_title, 
    roiMax = globalSettings$roiMax,	add = FALSE, axes = TRUE, ...){
  
  ## Open dataset
  w1Range <- sort(w1Range); w2Range <- sort(w2Range)
  if(is.matrix(in.folder$data) && !is.null(in.folder$w1) && 
     !is.null(in.folder$w2)){
    plotFile <- in.folder
    standardPlot <- TRUE
    if (!(xlab == '' && ylab == ''))
      roiMax <- FALSE
  }else{
    plotFile <- ucsf2D(file.name = in.folder$file.par$file.name,
                       w1Range = w1Range, w2Range = w2Range, file.par = in.folder$file.par)
    standardPlot <- FALSE
  }
  if(is.matrix(plotFile$data))
    plotFile$data <- matrix(rev(plotFile$data), ncol = ncol(plotFile$data))
  
  ## Set window to match the data 
  if(is.matrix(plotFile$data)){
    if (standardPlot){
      shiftRange <- matchShift(in.folder, w1=w1Range, w2=w2Range)
      minW1 <- which.min(abs(plotFile$w1 - shiftRange$w1[1]))
      minW2 <- which.min(abs(plotFile$w2 - shiftRange$w2[1]))
      maxW1 <- which.min(abs(plotFile$w1 - shiftRange$w1[2]))
      maxW2 <- which.min(abs(plotFile$w2 - shiftRange$w2[2]))
      plotFile$w1 <- plotFile$w1[minW1:maxW1]
      plotFile$w2 <- plotFile$w2[minW2:maxW2]
      plotFile$data <- plotFile$data[minW2:maxW2, minW1:maxW1]
    }else{
      w1Range <- range(plotFile$w1)
      w2Range <- range(plotFile$w2)
    }
  }
  
  ## Set new axes and 'clear' old plot with a box
  if(!add){
    uf <- in.folder$file.par$upfield_ppm
    df <- in.folder$file.par$downfield_ppm
    plot(0, 0, axes=FALSE, type = 'n', xlab=xlab, ylab=ylab, main=main, 
         xaxs='i', yaxs='i')
    par(usr = c(w2Range[2:1], w1Range[2:1]))
    rect(c(w2Range[2], w2Range[2], w2Range[2], uf[2] ),  
         c(w1Range[2], uf[1], w1Range[1], w1Range[1] ),
         c(w2Range[1], w2Range[1], df[2], w2Range[1] ),
         c(df[1], w1Range[1], w1Range[2], w1Range[2]), 
         col = "grey", border = NA)
  }
  
  ## Plot 'NA' for out of range ROIs
  if(!is.matrix(plotFile$data)){
    par(usr = c(-1,1,-1,1))
    text( 0, 0, "NA")	
    box(col=in.folder$graphics.par$fg, bty='o')
    return(invisible())
  }
  
  ## Set plot function for auto mode
  if(!type %in% c('auto', 'image', 'contour', 'filled'))
    type <- 'auto'
  if(type == 'auto'){
    if (standardPlot){
      if(length(plotFile$data) < 1100000)
        type <- 'contour'
      else
        type <- 'image'
    }else{
      if(length(plotFile$data) > 150000)
        type <- 'image'
      else
        type <- 'contour'  
    }
  }
  
  ## Set plotting function
  plotCur <- switch( which(c('image', 'contour', 'filled') == type),
                     function(...){ image(add = TRUE, ...)},
                     function(zlim, ...){ contour(zlim, add = TRUE, nlevels = nlevels, 
                                                  levels = seq(zlim[1], zlim[2], length.out = nlevels),  
                                                  drawlabels = FALSE, ...)},
                     function(zlim, ...){ fillCon(
                       levels = seq(zlim[1], zlim[2], length.out = nlevels), ...)}
  )
  
  ## Plot the data
  if( conDisp[1] )
    plotCur(x=plotFile$w2, y=plotFile$w1, z=plotFile$data, zlim = pos.zlim, 
            col = pos.color, ...)
  if( conDisp[2] )
    plotCur(x=plotFile$w2, y=plotFile$w1, z=plotFile$data, zlim = neg.zlim, 
            col = neg.color, ...)
  if( roiMax ){
    mShift <- maxShift(plotFile, invert = TRUE, conDisp = conDisp)
    points( mShift$w2, mShift$w1)
  }
  
  ## Draw plot labels
  box(col=in.folder$graphics.par$fg, bty='o')
  if(axes){
    par(usr = c(rev(par('usr')[1:2]), rev(par('usr')[3:4]))) 
    xlabvals <- pretty(w2Range, 10)    
    xlabvals <- c(xlabvals, par('usr')[2])
    n1 <- length(xlabvals)
    if (!is.na(in.folder$graphics.par$xtck) && 
        in.folder$graphics.par$xtck == 1)
      xlty <- 2
    else
      xlty <- 1
    if (!is.na(in.folder$graphics.par$ytck) && 
        in.folder$graphics.par$ytck == 1)
      ylty <- 2
    else
      ylty <- 1
    axis(side=1, at=min(w2Range) + (max(w2Range)-xlabvals),  
         labels = c( xlabvals[-n1], paste('    ', xlab)), lty=xlty,
         tck=in.folder$graphics.par$xtck, 
         cex.axis=in.folder$graphics.par$cex.axis)   
    ylabvals <- pretty(w1Range, 10)
    ylabvals <- c(ylabvals, par('usr')[4] )
    n1 <- length(ylabvals)
    axis(side=2, at=min(w1Range) + (max(w1Range) - ylabvals),  
         labels = c( ylabvals[-n1], paste('    ', ylab)), lty=ylty, 
         tck=in.folder$graphics.par$ytck, 
         cex.axis=in.folder$graphics.par$cex.axis) 
    par(usr = c(rev(par('usr')[1:2]), rev(par('usr')[3:4]))) 
  }
}

## Internal graphics function persp2D
## Generates a 3D perspective plot of a 2D NMR spectrum
##in.folder - Header of the file to be plotted, default is current spectrum
##w1Range - Chemical shift range in the indirect dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)          
##w2Range - Chemical shift range in the direct dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)
##color   - main color for the plot 
##xlab    - x axis label, default is the direct detected nucleus in PPM
##ylab    - y axis label, default is the indirect detected nucleus in PPM
##main    - main title for plot, default is the file name
##axes      - logical argument, TRUE plots axis labels
##...       - Additinal graphics paramaters can be passed to par()
##returns   - a 3D perspective plot of a 2D NMR spectrum   
persp2D <- function (
    in.folder = fileFolder[[wc()]], 
    w1Range = in.folder$graphics.par$usr[3:4],
    w2Range = in.folder$graphics.par$usr[1:2], 
    color = in.folder$graphics.par$pos.color,
    xlab = paste(in.folder$file.par$nucleus[2], 'PPM', sep=' '),
    ylab = paste(in.folder$file.par$nucleus[1], 'PPM', sep=' '),
    main = in.folder$file.par$user_title,                     
    axes = TRUE, expand = .25, r = 10, ...){
  
  ## Find the total number tiles that will have to be read
  in.folder$w1 <- seq(in.folder$file.par$upfield_ppm[1],
                      in.folder$file.par$downfield_ppm[1],
                      length.out = in.folder$file.par$matrix_size[1])
  in.folder$w2 <- seq(in.folder$file.par$upfield_ppm[2],
                      in.folder$file.par$downfield_ppm[2],
                      length.out = in.folder$file.par$matrix_size[2])
  w1Range <- sort(w1Range); w2Range <- sort(w2Range)
  out <- NULL
  t1 <- findInterval(w1Range, in.folder$w1, all.inside = TRUE)
  d1 <- findInterval(w2Range, in.folder$w2, all.inside = TRUE)
  t2 <- t1 + 1; d2 <- d1 + 1
  for(i in 1:2){
    out$w1[i] <- switch( which.min(c(
      abs(w1Range[i] - in.folder$w1[t1[i]]),
      abs(w1Range[i] - in.folder$w1[t2[i]]))), t1[i], t2[i])
    out$w2[i] <- switch( which.min(c(
      abs(w2Range[i] - in.folder$w2[d1[i]]),
      abs(w2Range[i] - in.folder$w2[d2[i]]))), d1[i], d2[i])
  }
  w1Tiles <- (ceiling(out$w1[1] / in.folder$file.par$block_size[1]):
                ceiling(out$w1[2] / in.folder$file.par$block_size[1]))
  w2Tiles <- (ceiling(out$w2[1] / in.folder$file.par$block_size[2]):
                ceiling(out$w2[2] / in.folder$file.par$block_size[2]))
  
  ## Redirect to contour/image plot if too much data is involved
  if((length(w1Tiles) + length(w2Tiles)) > 6 ){
    log_message('Perspective plots need a small spectral window, use zi()', 
          quote= FALSE)
    plot2D(in.folder = in.folder, w1Range = w1Range, w2Range = w2Range, 
           type = 'auto')
  }else{  
    ## Turn axes off 
    if( !axes ){
      xlab <- ylab <- zlab <- ''
      box <- FALSE
    }else{
      box <- TRUE 
      zlab <- 'Intensity'
    }
    
    ## Read data from binary
    plotFile <- ucsf2D(file.name = in.folder$file.par$file.name,
                       w1Range =w1Range, w2Range = w2Range, file.par = in.folder$file.par)
    if(!is.matrix (plotFile$data) ){
      ud()
      err('Perspective plots can not be scrolled outside the plot region')			
    }
    
    ## Generate the plot       
    persp(x = plotFile$w2, y = plotFile$w1, z = plotFile$data, 
          col = color, theta = in.folder$graphics.par$theta, 
          phi = in.folder$graphics.par$phi, asp = in.folder$graphics.par$asp, 
          box = box, main = main, xlab = xlab, ylab = ylab, zlab = zlab, 
          expand = expand, r = 10, ...)
  }     
}      


## Internal graphics function perspRsd
## Generates a 3D perspective plot of an RSD spectrum
##in.folder - Header of the file to be plotted, default is current spectrum
##w1Range - Chemical shift range in the indirect dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)          
##w2Range - Chemical shift range in the direct dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)
##color   - main color for the plot 
##xlab    - x axis label, default is the direct detected nucleus in PPM
##ylab    - y axis label, default is the indirect detected nucleus in PPM
##main    - main title for plot, default is the file name
##axes      - logical argument, TRUE plots axis labels
##...       - Additinal graphics paramaters can be passed to par()
##returns   - a 3D perspective plot of a 2D NMR spectrum   
perspRsd <- function(in.folder=fileFolder[[wc()]], 
                     w1Range=in.folder$graphics.par$usr[3:4],
                     w2Range=in.folder$graphics.par$usr[1:2],
                     color=in.folder$graphics.par$pos.color,
                     xlab=paste(in.folder$file.par$nucleus[2], 'PPM', sep=' '),
                     ylab=paste(in.folder$file.par$nucleus[1], 'PPM', sep=' '),
                     main=in.folder$file.par$user_title, axes=TRUE, expand=.25, r=10, ...){
  
  ## Find best w1/w2 matches	
  w1Range <- sort(w1Range)
  w2Range <- sort(w2Range)
  filePar <- in.folder$file.par
  in.folder$w1 <- seq(filePar$upfield_ppm[1], filePar$downfield_ppm[1], 
                      length.out=filePar$matrix_size[1])
  in.folder$w2 <- seq(filePar$upfield_ppm[2], filePar$downfield_ppm[2], 
                      length.out=filePar$matrix_size[2])
  t1 <- findInterval(w1Range, in.folder$w1, all.inside=TRUE)
  d1 <- findInterval(w2Range, in.folder$w2, all.inside=TRUE)
  t2 <- t1 + 1
  d2 <- d1 + 1
  out <- NULL
  for (i in 1:2){
    out$w1[i] <- switch(which.min(c(abs(w1Range[i] - in.folder$w1[t1[i]]),
                                    abs(w1Range[i] - in.folder$w1[t2[i]]))), t1[i], t2[i])
    out$w2[i] <- switch(which.min(c(abs(w2Range[i] - in.folder$w2[d1[i]]),
                                    abs(w2Range[i] - in.folder$w2[d2[i]]))), d1[i], d2[i])
  }
  w1Range <- c(in.folder$w1[(out$w1[1])], in.folder$w1[(out$w1[2])])
  w2Range <- c(in.folder$w2[(out$w2[1])], in.folder$w2[(out$w2[2])])
  winW1 <- in.folder$w1[(out$w1[1]:out$w1[2])]
  winW2 <- in.folder$w2[(out$w2[1]:out$w2[2])]
  
  ## Find tiles that will be plotted
  tiles <- NULL
  upShifts <- filePar$block_upfield_ppms
  downShifts <- filePar$block_downfield_ppms
  blockSizes <- filePar$block_size
  for (tNum in seq_along(filePar$block_size$w1)){
    
    ## Get the chemical shifts for the current tile
    blockW1 <- seq(upShifts$w1[tNum], downShifts$w1[tNum], 
                   length.out=blockSizes$w1[tNum])
    blockW2 <- seq(upShifts$w2[tNum], downShifts$w2[tNum],
                   length.out=blockSizes$w2[tNum])
    
    ## Check the window for the presence of any shift in the current block
    if (any(round(blockW1, 3) %in% round(winW1, 3)) && 
        any(round(blockW2, 3) %in% round(winW2, 3)))
      tiles <- c(tiles, tNum)
  }
  
  ## Redirect to contour/image plot if too much data is involved
  if (length(tiles) > 6){
    log_message('Perspective plots need a small spectral window, use zi()', 
          quote=FALSE)
    plotRsd2D(in.folder=in.folder, w1Range=w1Range, w2Range=w2Range, 
              type='auto')
  }else{  
    
    ## Turn axes off 
    if (!axes){
      xlab <- ylab <- zlab <- ''
      box <- FALSE
    }else{
      box <- TRUE 
      zlab <- 'Intensity'
    }
    
    ## Read data from binary
    plotFile <- rsd2D(file.name=in.folder$file.par$file.name,	w1Range=w1Range, 
                      w2Range=w2Range, file.par=in.folder$file.par)
    if (!is.matrix(plotFile$data)){
      ud()
      err('Perspective plots can not be scrolled outside the plot region')			
    }
    
    ## Generate the plot       
    persp(x=plotFile$w2, y=plotFile$w1, z=plotFile$data, col=color, 
          theta=in.folder$graphics.par$theta, phi=in.folder$graphics.par$phi, 
          asp=in.folder$graphics.par$asp,	box=box, main=main, xlab=xlab, 
          ylab=ylab, zlab=zlab, expand=expand, r=10, ...)
  }     
} 


##Internal graphics function fillCon
##x  - Acending numeric vector
##y  - Acending numeric vector
##z  - Data matrix with x rows and y columns
##levels - Numeric vector with the breaks for the contour plot
##col - Used as a place holder for compatibility
##returns a filled contour plot
fillCon <- function(x, y, z, levels, col, ...){
  col <- NULL
  par(...)
  if(min(levels) >= 0)
    col <- topo.colors(n = length(levels)) 
  else 
    col <- cm.colors( n = length(levels))	
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), col = col)
}


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
##      vp() resets the zero point of the plot without affecting the max of the 
##      zlimit. Offset shifts a given plot up/down from the vp() specified zero.
##...       - Additinal graphics paramaters can be passed to par()
##returns   - a plot of a 2D NMR spectrum  
plot1D <- function(
    in.folder = fileFolder[[wc()]],
    w1Range=in.folder$graphics.par$usr[3:4],
    w2Range=in.folder$graphics.par$usr[1:2],
    col = in.folder$graphics.par$proj.color, 
    type = in.folder$graphics.par$type,
    xlab = NULL, ylab = NULL,
    main = in.folder$file.par$user_title,
    roiMax = globalSettings$roiMax,
    add = FALSE, axes = TRUE, offset = 0, ... ){
  
  ## Remove any non line types
  if(!any(type == c('l', 'b', 'p')))
    type <- 'l'
  
  ## Redirect 2D NMR data to slice/projection function
  if(in.folder$file.par$number_dimensions > 1 )
    proj1D(in.folder = in.folder, w1Range = w1Range, w2Range = w2Range,
           col = col, main = main, add = add, axes = axes, type = type, ...)
  else{
    
    ## Force vertical position for all 1D plots
    w1Range <- c(in.folder$file.par$zero_offset - 
                   (w1Range[2] - in.folder$file.par$zero_offset) * 
                   globalSettings$position.1D, w1Range[2])
    
    ## Expand w1Range by offset
    col <- col
    offset <- diff(w1Range) * (offset / 100)
    
    ## Erase old plot if add is false
    if(!add){
      w2Range <- sort(w2Range) 
      if(is.null(xlab))
        xlab <- paste(in.folder$file.par$nucleus[1]) #, 'PPM', sep=' ')
      if(is.null(ylab))
        ylab <- 'Intensity'
      plot(0,0, axes=FALSE, type = 'n', xlab = xlab, ylab = ylab,
           main=main, xaxs='i',yaxs='i')
      par(usr = c(w2Range[2:1],w1Range[1:2]))
      rect(w2Range[2], w1Range[1], w2Range[1], w1Range[2], 
           col = 'grey', border = NA)
      rect(in.folder$file.par$downfield_ppm[1], w1Range[1], 
           in.folder$file.par$upfield_ppm[1], w1Range[2], 
           col= in.folder$graphics.par$bg, border = NA)	
    }
    
    ## Draw plot labels
    if(!(add)){
      box(col=in.folder$graphics.par$fg, bty='o')
      if(axes){
        par(usr = c(rev(par('usr')[1:2]),par('usr')[3:4])) 
        xlabvals <- pretty(w2Range, 10)
        xlabvals <- c(xlabvals, par('usr')[2])
        n1 <- length(xlabvals)
        if (!is.na(in.folder$graphics.par$xtck) && 
            in.folder$graphics.par$xtck == 1)
          xlty <- 2
        else
          xlty <- 1
        if (!is.na(in.folder$graphics.par$ytck) && 
            in.folder$graphics.par$ytck == 1)
          ylty <- 2
        else
          ylty <- 1
        axis(side=1, at=min(w2Range) + (max(w2Range)-xlabvals),  
             labels = c( xlabvals[-n1], paste('    ', xlab)), lty=xlty, 
             tck=in.folder$graphics.par$xtck, 
             cex.axis=in.folder$graphics.par$cex.axis)     
        axis(side=2, lty=ylty, tck=in.folder$graphics.par$ytck, 
             cex.axis=in.folder$graphics.par$cex.axis)
        par(usr = c(rev(par('usr')[1:2]),par('usr')[3:4]))
      }
    }
    
    ## Read new dataset and plot
    if( is.null(c(in.folder$data, in.folder$w2)) ){
      new.folder <- ucsf1D(file.name = in.folder$file.par$file.name,
                           w2Range = w2Range, file.par = in.folder$file.par)
      in.folder$file.par <- new.folder$file.par
      in.folder$w2 <- new.folder$w2
      in.folder$data <- new.folder$data
    }
    lines(x=in.folder$w2, y=rev(in.folder$data) + offset, col=col, type=type)
    if( roiMax &&  dev.cur() != 2 ){
      mShift <- maxShift( in.folder )
      points( mShift$w2, mShift$Height)
    }
  }
}

## Internal graphics function proj1D
## Projects 2D spectrum into a single dimension
## in.folder - Header from which projection data should be derived, by default 
##            data is derived from the current spectrum
## w1Range - Chemical shift range in the indirect dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)          
## w2Range - Chemical shift range in the direct dimension, default is the most
##          recent setting used with the file, the format is c(lower,upper)
## col    - color for 1D data and 1D slices/projections of 2D data
## proj.direct - 1 projects data across the direct axis, 2 across the indirect
## filter    - A vector capable function (e.g. min, max, sd) to filter data, 
##            non function arguments prompt users to select a slice 
## type - The type of plot, 'l' = lines, 'p' = points, 'b' = lines and points
## xy   - x and y coordinates can be passed from locator(1)
## ... - additional arguments can be passed to higher level functions (i.e - 
##			drawNMR), but are only present here to prevent argument mismatches
## returns the 1D trace/projection of the current spectrum
proj1D <- function (
    in.folder = fileFolder[[wc()]], 
    w1Range = in.folder$graphics.par$usr[3:4],
    w2Range = in.folder$graphics.par$usr[1:2],
    col = in.folder$graphics.par$proj.color,
    filter = globalSettings$filter,
    proj.direct = globalSettings$proj.direct,
    type = globalSettings$proj.type,
    xy = NULL, ...){
  
  ## Set chemical shift range from current spectrum
  if( in.folder$file.par$number_dimensions < 2 )
    err("1D projections can only be generated from 2D data")
  current.par <- c(in.folder$file.par, in.folder$graphics.par)
  lineCol <- fileFolder[[wc()]]$graphics.par$fg
  
  ## Generate data for 1D slices
  if( !is.function(filter) ){
    
    ## Have user select the slice and match click location to data
    if( is.null(xy) || length(unlist(xy)) != 2)
      xy <- locator(1)
    
    if(proj.direct  == 1 ){
      abline( h = xy$y, col=lineCol )
      w1Range <- rep(xy$y, 2)
      w2Range <- current.par$usr[1:2]
    }
    
    if(proj.direct == 2 ){
      abline(v = xy$x, col=lineCol )
      w1Range <- current.par$usr[3:4]
      w2Range <- rep(xy$x, 2)
    }
    
    ## Find 1D slice
    data.folder <- ucsf2D(in.folder$file.par$file.name, w1Range = w1Range,
                          w2Range=w2Range, file.par = in.folder$file.par)
    
    ## Generate data for projections
  }else{
    
    ## Read 2D file
    data.folder <- ucsf2D(in.folder$file.par$file.name, w1Range = w1Range,
                          w2Range=w2Range, file.par = in.folder$file.par)
    
    ## Generate projection
    data.folder$data <- apply(data.folder$data, proj.direct, filter)
  }
  
  ## Set z limits for new plot
  zlim <- c( 0, current.par$clevel * current.par$nlevels *
               current.par$noise_est )
  zlim[1] <- current.par$zero_offset - (zlim[2] - current.par$zero_offset) * 
    globalSettings$position.1D
  
  ## Setup outgoing file folder
  if( proj.direct == 1){
    in.folder$data <- data.folder$data
    in.folder$w2 <- data.folder$w2	
    newRange <- c( sort(current.par$usr[1:2], decreasing = TRUE), zlim )
    
  }else{
    in.folder$w2 <- data.folder$data
    in.folder$data <- data.folder$w1	
    newRange <- c( zlim, sort(current.par$usr[3:4], decreasing = TRUE) )
  }
  in.folder$file.par$number_dimensions <- 1
  in.folder$graphics.par$usr <- newRange
  op <- par('usr')
  par( usr = newRange )
  
  ## Plot the slice/projection
  plot1D(in.folder, add = TRUE, col = col, type = type)
  par(usr = op)
  bringFocus(-1) 
  return(invisible(data.folder))
  
}


## Finds the absolute max of a vector but returns original sign
## This function is the default for proj1D
pseudo1D <- function(x){range(x)[which.max(abs(range(x)))]}


################################################################################
##                                                                            ##
##                   Plotting functions for users                             ##
##                                                                            ##
################################################################################

## User file function fo
## Opens a sparky format spectrum or multiple spectra and generates a plot
## fileName - character string or vector; full pathname for file(s) to be opened
## ...  - Additional plotting options can be passed to drawNMR and par()

#' Open and Load Spectrum Files for `digestR`
#'
#' The `fo` function is responsible for loading one or more spectrum files, processing them, 
#' and integrating them into the `digestR` environment. If the user does not specify a 
#' filename, a GUI-based file selector will be triggered. The function will then read the 
#' specified files, verify their format and content, and make necessary updates to the global 
#' state. This function also provides error handling, ensuring that files with unsupported or 
#' corrupted formats do not crash the application.
#' 
#' @param fileName A character string specifying the path to a spectrum file to be loaded. 
#' If `fileName` is missing, the user will be prompted to select files via a GUI-based selector.
#' @param ... Additional arguments passed to other methods or functions.
#' 
#' @return Returns an invisibly the user-selected or provided file paths.
#' 
#' @examples 
#' \dontrun{
#' # Load a single spectrum file (assumes valid path and format)
#' fo("path_to_spectrum_file.txt")
#' 
#' # Trigger GUI-based file selector for user to select files
#' fo()
#' }
#' 
#' @seealso `createObj`, `dianaHead`, `log_message`, `myAssign`, `refresh`
#'
#' @export
# file_open <- function(fileName, ...){
  
#   ## Create any/all of the digestR objects that are missing
#   createObj()
  
#   ## Have user select all files they wish to open
#   if (missing(fileName)){
#     usrList <- sort(myOpen())
#     if(!length(usrList) || !nzchar(usrList))
#       return(invisible())
#   }else
#     usrList <- fileName
  
#   ## Read selected files
#   errors <- FALSE
#   fileNames <- names(fileFolder)
#   userTitles <- NULL
#   if (!is.null(fileFolder))
#     userTitles <- sapply(fileFolder, function(x) x$file.par$user_title)
#   for( i in 1:length(usrList) ){
    
#     ##Read Sparky Header and file info from binary
#     if( length(usrList) == 1 )
#     {
#       new.file <- tryCatch(dianaHead(file.name=usrList[i], print.info=TRUE), 
#                            error=function(er){
#                              errors <<- TRUE
#                              return(er$message)})
#       if (errors)
#         err(new.file)
#     }else
#     {
#       new.file <- tryCatch(dianaHead(file.name=usrList[i], print.info=FALSE), 
#                            error=function(er){
#                              errors <<- TRUE
#                              paste('\nOpening file "', basename(usrList[i]), '" produced an', 
#                                    ' error:\n"', er$message, '"', sep='')
#                            })
#       if (!is.list(new.file)){
#         cat(new.file, '\n\n')
#         flush.console()
#         next()
#       }
#     }
    
#     ## Make sure input files are of the correct format
#     if( length(new.file$file.par) == 0 )
#     {
#       log_message( paste('ERROR:', basename(usrList)[i], "is unreadable" ), 
#              quote=FALSE)
#       flush.console()
#       next()			
#     }
    
#     ## Fetch the default graphics settings 
#     new.file$graphics.par <- defaultSettings
    
#     ## Set initial plotting range
#     if( new.file$file.par$number_dimensions == 1 ) 
#     {
#       #		3rd parameter math basically vertically "centers" the origin
#       #		not desired functionality for DIANA	
#       #			new.file$graphics.par$usr <- c( new.file$file.par$downfield_ppm[1],
#       #					new.file$file.par$upfield_ppm[1], 
#       #					new.file$file.par$zero_offset - 
#       #							(new.file$file.par$max_intensity - new.file$file.par$zero_offset) 
#       #							* globalSettings$position.1D,
#       #					new.file$file.par$max_intensity )
      
#       new.file$graphics.par$usr <- c( new.file$file.par$downfield_ppm[1],
#                                       new.file$file.par$upfield_ppm[1], new.file$file.par$min_intensity,
#                                       new.file$file.par$max_intensity )			
#     }else
#     {         
#       new.file$graphics.par$usr <- c( new.file$file.par$downfield_ppm[2],
#                                       new.file$file.par$upfield_ppm[2], 
#                                       new.file$file.par$downfield_ppm[1],
#                                       new.file$file.par$upfield_ppm[1] )
#     }    
    
#     ## Make a new entry in the file folder if file is not already present 
#     filePar <- new.file$file.par
#     if(!new.file$file.par$file.name %in% fileNames)
#     {
      
#       ## Add 1D/2D spectra to the file folder
#       if (new.file$file.par$number_dimensions < 3)
#       {
#         if (new.file$file.par$user_title %in% userTitles)
#           new.file$file.par$user_title <- new.file$file.par$file.name
#         fileFolder[[(length(fileFolder) + 1)]] <- new.file
#         names(fileFolder)[length(fileFolder)] <- new.file$file.par$file.name
#       }else
#       {
        
#         ## Make duplicate entries in fileFolder for each z-slice in 3D spectra
#         w3 <- seq(filePar$upfield_ppm[3], filePar$downfield_ppm[3], 
#                   length.out=filePar$matrix_size[3])
#         for (j in seq_along(w3))
#         {
#           userTitle <- paste(basename(filePar$file.name), ' (z=', w3[j], ')', 
#                              sep='')
#           new.file$file.par$user_title <- userTitle
#           new.file$file.par$z_value <- w3[j]
#           fileFolder[[length(fileFolder) + 1]] <- new.file					
#           names(fileFolder)[length(fileFolder)] <- userTitle
#         }
#       }
#     }else
#     {
      
#       ## Update fileFolder entry if file is already present in fileFolder
#       fLoc <- match(new.file$file.par$file.name, fileNames)
#       if (new.file$file.par$number_dimensions < 3)
#       {
#         fileFolder[[fLoc]] <- new.file
#         if (new.file$file.par$user_title %in% userTitles)
#           new.file$file.par$user_title <- new.file$file.par$file.name
#       }else
#       {
#         for (j in fLoc){
#           zVal <- fileFolder[[j]]$file.par$z_value
#           new.file$file.par$user_title <- paste(basename(filePar$file.name), 
#                                                 ' (z=', zVal, ')', sep='')
#           new.file$file.par$z_value <- zVal
#           fileFolder[[j]] <- new.file
#         }
#       }
#     }
    
#     ## Reassign currentSpectrum
#     if (new.file$file.par$number_dimensions < 3)
#       currentSpectrum <- new.file$file.par$file.name
#     else
#       currentSpectrum <- userTitle
    
#     ## Tell user which files have been loaded
#     log_message( basename(usrList)[i], quote = FALSE )
#     flush.console()
#   }
  
#   ## Assign the new objects to the global environment
#   myAssign("fileFolder", fileFolder, save.backup = FALSE)
#   myAssign("currentSpectrum", currentSpectrum, save.backup = FALSE)
  
#   ## Save an undo point and refresh the active graphics
#   if( !is.null(fileFolder) ){
#     myAssign("currentSpectrum", currentSpectrum, save.backup = TRUE)
    
#     refresh(...)   
#   }
  
#   ##display error dialog
#   if (errors)
#     myMsg(paste('Errors occurred while opening files ',
#                 'Check the R console for details.', sep='\n'), icon='error')
  
#   return(invisible(usrList))
# }
#################################################################################
# Updated fo function: DDLM
#################################################################################

# # # Working Fo but opening errors still happening.

# file_open <- function(fileName, ...) {
  
#   ## Create any/all of the digestR objects that are missing
#   createObj()
  
#   ## Have user select all files they wish to open
#   if (missing(fileName)) {
#     usrList <- sort(myOpen())
#     if(!length(usrList) || !nzchar(usrList))
#       return(invisible())
#   } else {
#     usrList <- fileName
#   }
  
#   ## Read selected files
#   errors <- FALSE
#   fileNames <- names(fileFolder)
#   userTitles <- NULL
#   if (!is.null(fileFolder))
#     userTitles <- sapply(fileFolder, function(x) x$file.par$user_title)
  
#   ## Temporarily suppress all warnings
#   old_warn <- options(warn = -1)
  
#   for (i in 1:length(usrList)) {
    
#     ## Try to read the file while suppressing warnings
#     new.file <- tryCatch(
#       dianaHead(file.name = usrList[i], print.info = TRUE), 
#       error = function(cond) handleFoErrors(cond, usrList[i])
#     )
    
#     if (!is.list(new.file)) {
#       # If not a list, the operation failed, log the error and skip to the next iteration
#       log_message(paste("file opened", basename(usrList[i]), ":", new.file))
#       next
#     }
    
#     ## Make sure input files are of the correct format
#     if (length(new.file$file.par) == 0) {
#       log_message(paste('ERROR:', basename(usrList)[i], "is unreadable"), quote = FALSE)
#       flush.console()
#       next
#     }
    
#     ## Fetch the default graphics settings 
#     new.file$graphics.par <- defaultSettings
    
#     ## Set initial plotting range
#     if (new.file$file.par$number_dimensions == 1) {
#       new.file$graphics.par$usr <- c(new.file$file.par$downfield_ppm[1],
#                                      new.file$file.par$upfield_ppm[1], 
#                                      new.file$file.par$min_intensity,
#                                      new.file$file.par$max_intensity)
#     } else {
#       new.file$graphics.par$usr <- c(new.file$file.par$downfield_ppm[2],
#                                      new.file$file.par$upfield_ppm[2], 
#                                      new.file$file.par$downfield_ppm[1],
#                                      new.file$file.par$upfield_ppm[1])
#     }
    
#     ## Make a new entry in the file folder if file is not already present 
#     filePar <- new.file$file.par
#     if (!new.file$file.par$file.name %in% fileNames) {
      
#       ## Add 1D/2D spectra to the file folder
#       if (new.file$file.par$number_dimensions < 3) {
#         if (new.file$file.par$user_title %in% userTitles)
#           new.file$file.par$user_title <- new.file$file.par$file.name
#         fileFolder[[(length(fileFolder) + 1)]] <- new.file
#         names(fileFolder)[length(fileFolder)] <- new.file$file.par$file.name
#       } else {
        
#         ## Make duplicate entries in fileFolder for each z-slice in 3D spectra
#         w3 <- seq(filePar$upfield_ppm[3], filePar$downfield_ppm[3], 
#                   length.out = filePar$matrix_size[3])
#         for (j in seq_along(w3)) {
#           userTitle <- paste(basename(filePar$file.name), ' (z=', w3[j], ')', sep = '')
#           new.file$file.par$user_title <- userTitle
#           new.file$file.par$z_value <- w3[j]
#           fileFolder[[length(fileFolder) + 1]] <- new.file
#           names(fileFolder)[length(fileFolder)] <- userTitle
#         }
#       }
#     } else {
      
#       ## Update fileFolder entry if file is already present in fileFolder
#       fLoc <- match(new.file$file.par$file.name, fileNames)
#       if (new.file$file.par$number_dimensions < 3) {
#         fileFolder[[fLoc]] <- new.file
#         if (new.file$file.par$user_title %in% userTitles)
#           new.file$file.par$user_title <- new.file$file.par$file.name
#       } else {
#         for (j in fLoc) {
#           zVal <- fileFolder[[j]]$file.par$z_value
#           new.file$file.par$user_title <- paste(basename(filePar$file.name), 
#                                                 ' (z=', zVal, ')', sep = '')
#           new.file$file.par$z_value <- zVal
#           fileFolder[[j]] <- new.file
#         }
#       }
#     }
    
#     ## Reassign currentSpectrum
#     if (new.file$file.par$number_dimensions < 3) {
#       currentSpectrum <- new.file$file.par$file.name
#     } else {
#       currentSpectrum <- userTitle
#     }
    
#     ## Tell user which files have been loaded
#     log_message(paste("File", basename(usrList[i]), "opened successfully."))
#     flush.console()
#   }
  
#   ## Restore warning options
#   options(old_warn)
  
#   ## Assign the new objects to the global environment
#   myAssign("fileFolder", fileFolder, save.backup = FALSE)
#   myAssign("currentSpectrum", currentSpectrum, save.backup = FALSE)
  
#   ## Save an undo point and refresh the active graphics
#   if (!is.null(fileFolder)) {
#     myAssign("currentSpectrum", currentSpectrum, save.backup = TRUE)
#     refresh(...)   
#   }
  
#   ## Display error dialog
#   if (errors) {
#     myMsg(paste('Errors occurred while opening files', 'Check the R console for details.', sep = '\n'), icon = 'error')
#   }
  
#   return(invisible(usrList))

#   # Call the splash screen after loading the file
#   splashScreen()
  
#   # Ensure the graphical device is updated immediately
#   dev.flush()  # Forces a redraw of the graphics device
  
#   # If necessary, call refresh to update the main plot window
#   refresh(...)

# }
#######################################################################
file_open <- function(fileName, ...) {
  
  ## Create any/all of the digestR objects that are missing
  createObj()
  browser()  # Breakpoint 1 - after creating objects
  
  ## Have user select all files they wish to open
  if (missing(fileName)) {
    usrList <- sort(myOpen())
    if(!length(usrList) || !nzchar(usrList))
      return(invisible())
  } else {
    usrList <- fileName
  }
  browser()  # Breakpoint 2 - after getting the user list
  
  ## Read selected files
  errors <- FALSE
  fileNames <- names(fileFolder)
  userTitles <- NULL
  if (!is.null(fileFolder))
    userTitles <- sapply(fileFolder, function(x) x$file.par$user_title)
  browser()  # Breakpoint 3 - after reading fileFolder
  
  ## Temporarily suppress all warnings
  old_warn <- options(warn = -1)
  
  for (i in 1:length(usrList)) {
    browser()  # Breakpoint 4 - inside the loop before reading the file
    
    ## Try to read the file while suppressing warnings
    new.file <- tryCatch(
      dianaHead(file.name = usrList[i], print.info = TRUE), 
      error = function(cond) handleFoErrors(cond, usrList[i])
    )
    
    if (!is.list(new.file)) {
      log_message(paste("file opened", basename(usrList[i]), ":", new.file))
      next
    }
    
    ## Make sure input files are of the correct format
    if (length(new.file$file.par) == 0) {
      log_message(paste('ERROR:', basename(usrList)[i], "is unreadable"), quote = FALSE)
      flush.console()
      next
    }
    
    ## Fetch the default graphics settings 
    new.file$graphics.par <- defaultSettings
    
    ## Set initial plotting range
    if (new.file$file.par$number_dimensions == 1) {
      new.file$graphics.par$usr <- c(new.file$file.par$downfield_ppm[1],
                                     new.file$file.par$upfield_ppm[1], 
                                     new.file$file.par$min_intensity,
                                     new.file$file.par$max_intensity)
    } else {
      new.file$graphics.par$usr <- c(new.file$file.par$downfield_ppm[2],
                                     new.file$file.par$upfield_ppm[2], 
                                     new.file$file.par$downfield_ppm[1],
                                     new.file$file.par$upfield_ppm[1])
    }
    browser()  # Breakpoint 5 - after setting plotting range
    
    ## Make a new entry in the file folder if file is not already present 
    filePar <- new.file$file.par
    if (!new.file$file.par$file.name %in% fileNames) {
      
      ## Add 1D/2D spectra to the file folder
      if (new.file$file.par$number_dimensions < 3) {
        if (new.file$file.par$user_title %in% userTitles)
          new.file$file.par$user_title <- new.file$file.par$file.name
        fileFolder[[(length(fileFolder) + 1)]] <- new.file
        names(fileFolder)[length(fileFolder)] <- new.file$file.par$file.name
      } else {
        
        ## Make duplicate entries in fileFolder for each z-slice in 3D spectra
        w3 <- seq(filePar$upfield_ppm[3], filePar$downfield_ppm[3], 
                  length.out = filePar$matrix_size[3])
        for (j in seq_along(w3)) {
          userTitle <- paste(basename(filePar$file.name), ' (z=', w3[j], ')', sep = '')
          new.file$file.par$user_title <- userTitle
          new.file$file.par$z_value <- w3[j]
          fileFolder[[length(fileFolder) + 1]] <- new.file
          names(fileFolder)[length(fileFolder)] <- userTitle
        }
      }
    } else {
      
      ## Update fileFolder entry if file is already present in fileFolder
      fLoc <- match(new.file$file.par$file.name, fileNames)
      if (new.file$file.par$number_dimensions < 3) {
        fileFolder[[fLoc]] <- new.file
        if (new.file$file.par$user_title %in% userTitles)
          new.file$file.par$user_title <- new.file$file.par$file.name
      } else {
        for (j in fLoc) {
          zVal <- fileFolder[[j]]$file.par$z_value
          new.file$file.par$user_title <- paste(basename(filePar$file.name), 
                                                ' (z=', zVal, ')', sep = '')
          new.file$file.par$z_value <- zVal
          fileFolder[[j]] <- new.file
        }
      }
    }
    
    ## Reassign currentSpectrum
    if (new.file$file.par$number_dimensions < 3) {
      currentSpectrum <- new.file$file.par$file.name
    } else {
      currentSpectrum <- userTitle
    }
    
    ## Tell user which files have been loaded
    log_message(paste("File", basename(usrList[i]), "opened successfully."))
    flush.console()
  }
  
  ## Restore warning options
  options(old_warn)
  browser()  # Breakpoint 6 - after file processing
  
  ## Assign the new objects to the global environment
  myAssign("fileFolder", fileFolder, save.backup = FALSE)
  myAssign("currentSpectrum", currentSpectrum, save.backup = FALSE)
  
  ## Save an undo point and refresh the active graphics
  if (!is.null(fileFolder)) {
    myAssign("currentSpectrum", currentSpectrum, save.backup = TRUE)
    refresh(...)   
  }
  
  ## Display error dialog
  if (errors) {
    myMsg(paste('Errors occurred while opening files', 'Check the R console for details.', sep = '\n'), icon = 'error')
  }
  
  return(invisible(usrList))

  # Call the splash screen after loading the file
  splashScreen()
  
  # Ensure the graphical device is updated immediately
  dev.flush()  # Forces a redraw of the graphics device
  
  # If necessary, call refresh to update the main plot window
  refresh(...)

}


handleFoErrors <- function(cond, fileName = NULL) {
  # Set errors flag to TRUE to indicate that an error occurred
  errors <<- TRUE
  
  # Customize the error message based on the error type or condition
  if (grepl("unused argument \\(cond\\)", cond$message)) {
    log_message <- "An unused argument error occurred in the function. Please check the function arguments."
  } else if (grepl("truncating string with embedded nuls", cond$message)) {
    # Specific handling for "truncating string with embedded nuls" warnings/errors
    log_message <- paste("Warning: A truncation error occurred due to embedded null characters in the file:", 
                         if (!is.null(fileName)) basename(fileName) else "", sep = " ")
  } else {
    # Generic error handling for other types of errors
    log_message <- paste("An error occurred while processing", if (!is.null(fileName)) basename(fileName) else "", ":", cond$message)
  }
  
  # Log the error message to a log file or suppress the output
  write(log_message, file = "fo_error_log.txt", append = TRUE)
  
  # Optionally print the message for debugging purposes (comment out to suppress completely)
  cat(log_message, "\n")
  
  # Return NULL or an appropriate value to allow the function to continue
  return(NULL)
}

#######################################################################################
## User file function fc
## Closes a user defined file from a list
## usrList - character string/vector; files to close
fc <- function(usrList=NULL){
  
  ## Make list of all files
  current <- wc()
  allFiles <- names(fileFolder)    
  fileNames <- getTitles(allFiles, FALSE)
  
  ## Have user select files to close   
  if (is.null(usrList)){
    usrList <- mySelect(fileNames, multiple=TRUE, index=TRUE,
                        title='Select files to close:', preselect=fileNames[current])
    if ( length(usrList) == 0 || !nzchar(usrList) )
      return(invisible())
  }else{
    usrList <- as.vector(na.omit(match(usrList, allFiles)))
    if (!length(usrList))
      return(invisible())
  }
  
  ## Reassign current spectrum if current is being deleted
  redraw <- FALSE
  if( !is.na(match(current, usrList)) ){
    redraw <- TRUE
    
    allFiles <- 1:length(allFiles)
    allFiles <- allFiles[-usrList]
    
    new.current <- rev(allFiles[allFiles < current])[1]
    if(is.na(new.current))
      new.current <- allFiles[allFiles > current][1]
    if(is.na(new.current))
      myAssign("currentSpectrum", NULL, save.backup = FALSE)
    else
      myAssign("currentSpectrum", names(fileFolder)[new.current],  
               save.backup = FALSE)
  }
  
  ## Remove selected files from the overlay list
  if( !is.null(overlayList) ){
    oldOl <- match( names(fileFolder)[usrList], overlayList )
    if( any(is.na(oldOl)) ) 
      oldOl <- oldOl[ -which(is.na(oldOl)) ]
    if( length(oldOl) > 0 ){
      overlayList <- overlayList[ -oldOl ]
      if( length( overlayList) == 0 )
        overlayList <- NULL
      redraw <- TRUE    
      myAssign('overlayList', overlayList , save.backup = FALSE)			
    }
  }
  
  
  ## Delete user selected files from file folder      
  fileFolder <- fileFolder[-usrList]
  if(length(fileFolder) == 0){
    fileFolder <- NULL
    redraw <- FALSE
    log_message('The file folder is now empty', quote = FALSE )
    
    if(length(which(dev.list() == 2)) == 1){
      dev.set(which = 2)
      plot(1, 1, col='transparent', axes=FALSE, xlab='', ylab='')
      text(1, 1, 'No files are open.\nUse fo()', cex=1)            
      
    }
    
    if(length(which(dev.list() == 3)) == 1){
      dev.set(which = 3)
      plot(1, 1, col='transparent', axes=FALSE, xlab='', ylab='')
      text(1, 1, 'ROIs will apear here \nuse roi() to designate an ROI', 
           cex=1)       
    }
    
    if(length(which(dev.list() == 4)) == 1){
      dev.set(which = 4)          
      plot(1, 1, col='transparent', axes=FALSE, xlab='', ylab='')
      text(1, 1, 'Active ROIs will appear here \nuse rs()', cex=1) 
    }
    
    ## Leave main plot window and console active
    if(length(which(dev.list() == 2)) == 1)          
      dev.set(which = 2)
    bringFocus(-1)
  }
  
  ## Assign fileFolder and save backup
  if ( !redraw )
    myAssign('fileFolder', fileFolder)
  else{
    myAssign('fileFolder', fileFolder, save.backup = FALSE)
    myAssign('currentSpectrum', currentSpectrum, save.backup = TRUE)
  }
  
  ## Refresh open plots
  if(exists('fileFolder') && !is.null(fileFolder)){
    if( redraw )
      refresh()
    else
      refresh(main.plot = FALSE, sub.plot = FALSE)
  }        
}

## User file function ss
## Switch the active spectrum to another file in memory
## ...  - Additional plotting options can be passed to drawPeptides and par()
ss <- function(...){
  
  ##Generate a list of file names    
  current <- wc()
  fileNames <- getTitles(names(fileFolder), FALSE)
  
  ##Have the user select a spectrum
  usrList <- mySelect(fileNames, multiple=FALSE, preselect = fileNames[current], 
                      title='Select a spectrum:', index=TRUE)
  
  ## Do nothing if user selects cancel
  if ( length(usrList) == 0 || !nzchar(usrList) )
    return(invisible())
  
  ##Set new file as current spectrum and refresh the plots
  myAssign('currentSpectrum', names(fileFolder)[usrList])
  refresh(multi.plot = FALSE, ...)  
}

#' Undo Last Action
#'
#' This function allows users to undo the most recent action performed in the application. It restores the state of various global variables to their previous values, effectively reversing the most recent changes made.
#'
#' @return NULL (The function primarily works with global variables and their modifications.)
#'
#' @details
#' The \code{undo} function works by decrementing the undo index in the \code{oldFolder} variable, which stores the history of changes. It resets each of the global files and adds new zoom changes to the \code{oldFolder} index.
#' 
#' @importFrom tcltk2
#' @export
ud <- function(){
  
  current <- wc()
  
  if( !exists('oldFolder') || oldFolder$undo.index < 2 ){
    cat('Cannot undo \n')
    return(invisible())
  }
  
  ## Set undo index
  oldFolder$undo.index <- oldFolder$undo.index - 1 
  
  ## Reset each of the global files    
  save.list <- names(oldFolder)
  for(i in 1:length(save.list)){
    if(names(oldFolder)[i] %in% c('undo.index', 'assign.index', 'zoom.list', 
                                  'zoom.history'))
      next()
    out.file <- oldFolder[[i]][[oldFolder$undo.index]]
    suppressWarnings( if(!length(out.file[[1]][1]) || 
                         is.na( out.file[[1]][1] )) out.file <- NULL )   
    myAssign( save.list[i], out.file, save.backup = FALSE)
    
  }
  
  ## Add new zoom changes to oldFolder index
  if(oldFolder$zoom.history[[oldFolder$undo.index]])
    oldFolder$zoom.list[[(length(oldFolder$zoom.list) + 1)]] <- 
    fileFolder[[wc()]]$graphics.par$usr
  
  ## Save oldFolder to global environment 
  myAssign( "oldFolder", oldFolder, save.backup = FALSE) 
  refresh()  
}

#' Refreshes the main plot without changing settings
#' ...  - Additional plotting options can be passed to drawPeptides and par() 
#' @import tcltk2
#'
#' @export

redraw <- function (...) {
  
  ##Redraw the open spectra
  refresh(sub.plot = FALSE, multi.plot = FALSE, ...)
}

## Redraws the plot with the graphic type automatically set
## ...  - Additional plotting options can be passed to drawPeptides and par() 
## note: this requires fo() to be used first  
da <- function (...) {
  
  setGraphics(type = 'auto')
  
  ##Redraw the open spectra
  refresh(sub.plot = FALSE, multi.plot = FALSE, ...)
}

## Redraws the plot as an image map 
## ...  - Additional plotting options can be passed to drawPeptides and par() 
## note: this requires fo() to be used first  
di <- function (...) {
  
  setGraphics(type = 'image')
  
  ##Redraw the open spectra
  refresh(sub.plot = FALSE, multi.plot = FALSE, ...)
}

## Redraws the plot as a contour plot
## ...  - Additional plotting options can be passed to drawPeptides and par() 
## note: this requires fo() to be used first
dr <- function (nlevels=20, ...) {
  
  setGraphics(type = 'contour', nlevels=nlevels)
  
  ## redraw the current spectrum
  refresh(sub.plot = FALSE, multi.plot = FALSE, ...)
}

## Redraws the plot as a filled contour map
## ...  - Additional plotting options can be passed to drawPeptides and par() 
## note: this requires fo() to be used first
drf <- function (...) {
  
  setGraphics(type='filled')
  
  ## Draw the new spectrum                    
  refresh(sub.plot = FALSE, multi.plot = FALSE, ...)
}

## Draws the image as a perspective plot
## ...  - Additional plotting options can be passed to drawPeptides and par() 
## note: this requires fo() to be used first
dp <- function (...) {
  
  setGraphics(type='persp')
  
  ## Draw the new spectrum                    
  refresh(overlay = FALSE, sub.plot = FALSE, multi.plot = FALSE, ...)
  
}

## Rotates a perspective plot clockwise
## ...  - Additional plotting options can be passed to drawPeptides and par()
rotc <- function(degrees = 10, ...){
  
  current <- wc()
  newRot <- fileFolder[[current]]$graphics.par$theta + degrees
  setGraphics(type='persp', theta = newRot)
  
  ## Draw the new spectrum                    
  refresh(overlay = FALSE, sub.plot = FALSE, multi.plot = FALSE, ...)
}

## Rotates a perspective plot counter clockwise
## ...  - Additional plotting options can be passed to drawPeptides and par()
rotcc <- function(degrees = 10, ...){
  rotc (degrees = -degrees, ...)
}

## Rotates a perspective plot up
## ...  - Additional plotting options can be passed to drawPeptides and par()
rotu <- function(degrees = 10, ...){
  
  current <- wc()
  newRot <- fileFolder[[current]]$graphics.par$phi + degrees
  setGraphics(type='persp', phi = newRot)
  
  ## Draw the new spectrum                    
  refresh(overlay = FALSE, sub.plot = FALSE, multi.plot = FALSE, ...)
}

## Rotates a perspective plot down
## ...  - Additional plotting options can be passed to drawPeptides and par()  
rotd <- function(degrees = 10, ...){
  rotu(degrees = -degrees, ...)
}

## Rotates a perspective plot 360 degrees in steps of 10 degrees
## ...  - Additional plotting options can be passed to drawPeptides and par() 
spin <- function (...){
  degrees <- 10
  current <- wc()
  for( i in 1:36){
    newRot <- fileFolder[[current]]$graphics.par$theta + degrees
    setGraphics(type='persp', theta = newRot, save.backup = FALSE )
    refresh(overlay = FALSE, sub.plot = FALSE, multi.plot = FALSE, ...)
  }
}

## User function vp
## Sets the vertical position of all 1D spectra and slices
## position - numeric value between 0 and 100 that expresses the percentage of 
##          the visible range shown below the mean spectral noise level.
##					The default value 5, creates plots with 5% of the intensity range
##          below the average noise.
## note: if no offset is given, then the current vertical offset is returned
## note2: This function modifies all of the open spectra 
## returns refreshed plots with new vertical positions if offset is not NULL or
##         the current vertical offset if offset is NULL
vp <- function(position = NULL){
  
  
  ## Print the current vertical position if offset is NULL
  if(is.null(position)){
    
    position <- globalSettings$position.1D	
    if(position <= 1)
      position <- (position) / 2 * 100
    else
      position <- 100 - (1 / position) * 50 
    log_message('Current vertical position:', quote = FALSE)
    return(position)
  }
  
  ## Check for valid offset values
  if(is.na(suppressWarnings(as.integer(position))) || position < 0 || 
     position >= 100)
    err('Vertical position value must be between 0 and 100 (5 is the default)')
  
  ## Update vertical position with new offset values
  if(position <= 50)
    position <- position / 50
  else
    position <- .5 / ((100 - position ) / 100) 
  
  setGraphics( position.1D = position, refresh.graphics = TRUE)					
} 

## User function vpu
## Increases vertical position of all 1D spectra and slices
## p - numeric value between 0 and 100; the amount to increase the 
##	vertical position by
vpu <- function(p=5){
  
  ## Check for valid offset values
  if(is.na(suppressWarnings(as.integer(p))) || p < 0 || 
     p >= 100)
    err('Vertical position value must be between 0 and 100 (5 is the default)')
  
  ## Update vertical position with new offset values
  if (p <= 50)
    p <- p / 50
  else
    p <- .5 / ((100 - p ) / 100) 
  
  prevPos <- globalSettings$position.1D		
  setGraphics(position.1D=prevPos + p, refresh.graphics=TRUE)					
} 

## User function vpd
## Decreases vertical position of all 1D spectra and slices
## p - numeric value between 0 and 100; the amount to decrease the 
##	vertical position by
vpd <- function(p=5){
  
  ## Check for valid offset values
  if(is.na(suppressWarnings(as.integer(p))) || p < 0 || 
     p >= 100)
    err('Vertical position value must be between 0 and 100 (5 is the default)')
  
  ## Update vertical position with new offset values
  if (p <= 50)
    p <- p / 50
  else
    p <- .5 / ((100 - p ) / 100) 
  
  prevPos <- globalSettings$position.1D		
  setGraphics(position.1D=prevPos - p, refresh.graphics=TRUE)					
}


################################################################################
#                                                                              #
#                        Internal zoom functions                               #
#                                                                              #
################################################################################  

## Internal utility function newRange
## Returns an expanded or contracted range
## x  - a optional numeric list from which range is calculated
## r  - A numeric range (length 2) to be extended or contracted
## f  - Numeric expansion/contraction factor, positive values are treated
##      as an expansion factor, negative values are treated as a contraction
##      factor. 
## checkF - logical; if TRUE values less than -1 will be reset to .9999
## returns a new numeric range of length two
newRange <- function(x, r = range(x, na.rm = TRUE), f = 0.05, checkF = TRUE){
  
  ## Error checking
  f <- suppressWarnings(as.numeric(f))
  if( is.na(f) )
    stop('The new range factor must be numeric')
  r <- suppressWarnings(as.numeric(r))
  if( any(is.na(r)) || length(r) != 2 )
    stop('Only numeric ranges can be used')	
  f <- f/2
  
  ## Change range
  if( f > 0 )
    return(r + c(-f, f) * diff(r))
  if( f < 0 ){
    f <- -f
    if( checkF && f > 1 )
      f <- 0.4999
    return( r - c(-f, f) * diff(r) )		
  }
  
  return(r)
}

## Internal function pan
## direction - can be set to 'u', 'd', 'l' or 'r' for up, down, left or right,
##             these valuse can also be combined in a list input: c('l', 'u')
## p         - Numeric argument, percentage of the current window to move 
## save.backup - Logical argument, TRUE saves an undo point
## ... - Additional parameters can be passed to drawPeptides
pan <- function( direction, p = 5, save.backup = TRUE, ...){
  
  ## Error checking
  if(missing(direction))
    err("No pan direction was provided")
  if( all( is.na(match(direction, c('u', 'd', 'l', 'r'))) ))
    err("The pan direction must be set to u, d, l, or r")
  p <- suppressWarnings(as.numeric(p))
  if( is.na(p) )
    err('The panning percentage must be numeric')
  
  ## Define current
  current <- wc()
  nDim <- fileFolder[[current]]$file.par$number_dimensions
  usr <- fileFolder[[current]]$graphics.par$usr
  p <- p/100
  
  ## Set new window
  if( any(direction == 'r' ))
    usr[1:2] <- usr[1:2] - (-p) * diff( sort(usr[1:2]) )
  
  if( any(direction == 'l' ))
    usr[1:2] <- usr[1:2] - p * diff( sort(usr[1:2]) )
  
  if( any(direction == 'u' )){
    if(nDim > 1)
      usr[3:4] <- usr[3:4] - p * diff( sort(usr[3:4]) )
    else
      usr[4] <- max(usr[3:4]) + p * max(usr[3:4])
  }
  if( any(direction == 'd' )){
    if(nDim > 1)
      usr[3:4] <- usr[3:4] - (-p) * diff( sort(usr[3:4]) )
    else
      usr[4] <- max(usr[3:4]) + (-p) * max(usr[3:4])
  }
  
  ## Assign the new window and refresh
  setGraphics( usr = usr, save.backup = save.backup )
  refresh(  sub.plot = FALSE, multi.plot = FALSE, ...)
  
} 


################################################################################
#                                                                              #
#                     2D zoom functions for users                              #
#                                                                              #
################################################################################  

## User zoom function zp
## Zoom the plot to previous zoom set
zp <- function(){
  current <- wc()
  zoom.par <- oldFolder$zoom.list
  usr <- unlist(zoom.par[2])
  
  ## Check for valid previous zoom
  if(length( zoom.par ) > 1 && !is.null(usr) && !is.na(usr) ){
    fileFolder[[current]]$graphics.par$usr <- usr
    myAssign('zoom', fileFolder)
    
    oldFolder$zoom.list <- zoom.par[-(length(zoom.par))]   
    myAssign("oldFolder", oldFolder, save.backup = FALSE)
    refresh(sub.plot = FALSE, multi.plot = FALSE)
  }else
    log_message('Cannot zoom previous', quote=FALSE)
}

## User zoom function zi
## Zooms the plot in to a narrower chemical shift range
## p  - The percentage reduction in range to be applied
## ...  - Additional plotting options can be passed to drawPeptides and par()
zi <- function ( p = 25, ...) {
  
  p <- suppressWarnings(as.numeric(p))
  if( is.na(p) )
    err('The zoom factor must be numeric')
  if (p > 100)
    err('The zoom factor may not exceed 100')
  
  ## Find current parameters
  current <- wc()
  usr <- fileFolder[[current]]$graphics.par$usr 
  
  ## Set new Ranges
  usr[1:2] <- sort(newRange(usr[1:2], f = -p/100))
  usr[3:4] <- sort(newRange(usr[3:4], f = -p/100))	
  
  ## Draw the new spectrum                    
  setGraphics( usr = usr )
  refresh(sub.plot = FALSE, multi.plot = FALSE, ...)
}

## User zoom function zo
## Zooms the plot in to a wider chemical shift range
## p  - The fractional reduction in range to be applied
## ...  - Additional plotting options can be passed to drawPeptides and par()
zo <- function (p = 25, ...) {
  zi( p = -p, ...)
}

## User pan function pr
## Pans spectrum to the right
## p   - The percentage of the current window to move the window
## ... - Additional parameters can be passed to drawPeptides
pr <- function( p = 5, ... ){
  
  pan( direction = 'r', p = p, ...)
  
}

## User pan function pl
## Pans spectrum to the left
## p   - The percentage of the current window to move the window
## ... - Additional parameters can be passed to drawPeptides
pl <- function( p = 5, ...){
  
  pan( direction = 'l', p = p, ...)
  
}

## User pan function pu
## Pans spectrum up
## p   - The percentage of the current window to move the window
## ... - Additional parameters can be passed to drawPeptides
pu <- function( p = 5, ... ){
  
  pan( direction = 'u', p = p, ...)
  
}

## User pan function pd
## Pans spectrum down
## p   - The percentage of the current window to move the window
## ... - Additional parameters can be passed to drawPeptides
pd <- function( p = 5, ... ){
  
  pan( direction = 'd', p = p, ...)
  
}

## Redraws spectrum at full chemical shift range as an image
## ...  - Additional plotting options can be passed to drawPeptides and par()
ff <- zf <- function (...) {
  ## Find current spectrum
  current <- wc()
  nDim <- fileFolder[[current]]$file.par$number_dimensions
  
  ## Reset chemical shift range
  if(nDim > 1){
    fileFolder[[current]]$graphics.par$usr = c(
      fileFolder[[current]]$file.par$downfield_ppm[2],
      fileFolder[[current]]$file.par$upfield_ppm[2],
      fileFolder[[current]]$file.par$downfield_ppm[1],
      fileFolder[[current]]$file.par$upfield_ppm[1])
  }else{
    
    if(nchar(globalSettings$geneDisp)>0)
    {
      fileFolder[[current]]$graphics.par$usr <- c( 
        fileFolder[[current]]$file.par$downfield_ppm[1],
        fileFolder[[current]]$file.par$upfield_ppm[1], 
        fileFolder[[current]]$file.par$zero_offset - 
          (fileFolder[[current]]$file.par$aa_max_intensity - 
             fileFolder[[current]]$file.par$zero_offset) 
        * globalSettings$position.1D,
        fileFolder[[current]]$file.par$aa_max_intensity )			
    }else
    {
      fileFolder[[current]]$graphics.par$usr <- c( 
        fileFolder[[current]]$file.par$downfield_ppm[1],
        fileFolder[[current]]$file.par$upfield_ppm[1], 
        fileFolder[[current]]$file.par$zero_offset - 
          (fileFolder[[current]]$file.par$max_intensity - 
             fileFolder[[current]]$file.par$zero_offset) 
        * globalSettings$position.1D,
        fileFolder[[current]]$file.par$max_intensity )
    }
  }                   
  
  ## Draw the new spectrum                    
  myAssign('zoom', fileFolder)
  refresh(sub.plot = FALSE, multi.plot = FALSE, ...) 
}

## User zoom/pan function zz
## Interactive click zooming and panning
## ...  - Additional plotting options can be passed to drawPeptides and pan 
zz <- function (...) {
  
  current <- wc() 
  lineCol <- fileFolder[[current]]$graphics.par$fg
  nDim <- fileFolder[[current]]$file.par$number_dimensions
  
  ## Opens the main plot window if not currently opened
  if (is.na(match(2, dev.list())))
    refresh(multi.plot=FALSE, sub.plot=FALSE)
  cw(dev=2)
  
  ## Give the user some instructions
  cat(paste('In the main plot window:\n',  
            ' Left-click two points inside the plot to zoom\n',  
            ' Left-click outside the plot to pan\n', 
            ' Right-click to exit\n'))
  flush.console()
  hideGui()
  
  while(TRUE){
    
    ## Tell the user what to do
    op <- par('font')
    par( font = 2 )
    legend("topleft", c('LEFT CLICK TO ZOOM','RIGHT CLICK TO EXIT'), 
           pch=NULL, bty='n', text.col=lineCol)
    par(font = op)
    
    xy <- data.frame(locator(1))
    if(length(xy) == 0)
      break()
    
    ## Pan if click is outside the plotting range
    direction <- NULL
    if( xy$x < fileFolder[[current]]$graphics.par$usr[1] )
      direction <- c(direction, 'l')
    if( xy$x > fileFolder[[current]]$graphics.par$usr[2] )
      direction <- c(direction, 'r')		
    if( xy$y > max(par('usr')[3:4]) )
      direction <- c(direction, 'd')
    if( xy$y < min(par('usr')[3:4]) )
      direction <- c(direction, 'u')
    if( !is.null(direction) ){
      pan(direction = direction, save.backup = FALSE, ...)
      next()	
    }
    
    ## Zoom if first click inside the plotting area	
    abline(v=xy$x, col=lineCol)
    if( nDim != 1 )
      abline(h=xy$y, col=lineCol )    
    xy2 <- data.frame(locator(1))
    if(length(xy2) == 0)
      break()
    abline(v=xy2$x, col=lineCol)
    if (nDim != 1)
      abline(h=xy2$y, col=lineCol)   
    xy <- rbind(xy, xy2)
    
    ## Do not rescale lower bound on 1D data
    if (nDim == 1){
      xy$y <- sort(xy$y)
      abline( h = xy$y[2], col=lineCol )
      xy$y[1] <- fileFolder[[current]]$file.par$zero_offset - 
        (xy$y[2] - fileFolder[[current]]$file.par$zero_offset) * 
        globalSettings$position.1D
    }else
      xy$y <- rev(sort(xy$y))
    
    ## Assign the new plotting regions and refresh
    setGraphics( usr = c(sort(xy$x), xy$y), save.backup = FALSE )	
    refresh(sub.plot = FALSE, multi.plot = FALSE, ...)
  }
  
  ## Assign the new plotting regions to the correct file
  showGui()
  myAssign("zoom", fileFolder, save.backup = TRUE)	
  refresh( sub.plot = FALSE, multi.plot = FALSE )
  
}

## Automatically zoom in on a point using mouse interface
## w1Delta  - Chemical shift range for new w1 window, NULL sets window to 
##            2.5ppm for all nuclei other than 1H and 0.25 ppm for 1H
## w2Delta  - Chemical shift range for new w2 window, NULL sets window to 
##            2.5ppm for all nuclei other than 1H and 0.25 ppm for 1H
## ...  - Additional plotting options can be passed to drawPeptides and par()
pz <- function (w1Delta = NULL, w2Delta = NULL, ...) {
  
  ## Find current spectrum and establish chemical shift range
  current <- wc()
  lineCol <- fileFolder[[current]]$graphics.par$fg
  nDim <- fileFolder[[current]]$file.par$number_dimensions
  
  ## Opens the main plot window if not currently opened
  if (is.na(match(2, dev.list())))
    refresh(multi.plot=FALSE, sub.plot=FALSE)
  cw(dev=2)
  
  ## Give the user some instructions
  cat(paste('\nIn the main plot window:\n',  
            ' Left-click to view chemical shifts\n',  
            ' Right-click to exit\n\n')) 
  flush.console()
  hideGui()
  
  ## Tell the user what to do
  op <- par('font')
  par(font=2)
  legend("topleft", c('LEFT CLICK TO ZOOM','RIGHT CLICK TO EXIT'), 
         pch=NULL, bty='n', text.col=lineCol)
  par(font=op)
  
  ## Set the chemical shift ranges of the point zoom window
  if(is.null(w1Delta)){
    if(fileFolder[[current]]$file.par$nucleus[1] != '1H')
      w1Delta <- 2.5 / 2
    else
      w1Delta <- .1 / 2
  }
  if(nDim == 1 && is.null(w2Delta))
    w2Delta <- w1Delta
  if(is.null(w2Delta)){
    if(fileFolder[[current]]$file.par$nucleus[2] != '1H')
      w2Delta <- 2.5 / 2
    else
      w2Delta <- .1 / 2
  }
  
  ## Have user locate the place to zoom
  xy=data.frame(locator(1, type='n'))
  if(length(xy) == 0){
    showGui()
    return(invisible())
  }
  abline(v=xy$x, col=lineCol )
  abline(h=xy$y, col=lineCol )
  
  ## Set new chemical shift ranges
  w2Range <- c(xy$x + w2Delta, xy$x - w2Delta)
  if(nDim != 1)
    w1Range <- c(xy$y + w1Delta, xy$y - w1Delta) 
  else
    w1Range <- c(fileFolder[[current]]$file.par$zero_offset - 
                   (xy$y - fileFolder[[current]]$file.par$zero_offset) * 
                   globalSettings$position.1D, xy$y)
  fileFolder[[current]]$graphics.par$usr <- c(w2Range, w1Range) 
  myAssign('zoom', fileFolder)
  showGui()
  
  ## Draw the new spectrum                    
  refresh(sub.plot = FALSE, multi.plot = FALSE, ...)
}

## Centers the window on the largest visable peak
## massCenter  - logical argument; TRUE centers peaks by center of mass,
##               false centers peaks by maximum signal observed
## note: graphics settings are used for choosing between negative and positive 
##       signals.
zc <- function ( massCenter = TRUE ){
  
  if( !is.logical (massCenter) ){
    massCenter = FALSE
    err('massCenter must be set to either TRUE or FALSE')
  }
  
  ## Find the biggest peak in the window
  current <- wc()
  w1Range <- fileFolder[[current]]$graphics.par$usr[3:4]
  w2Range <- fileFolder[[current]]$graphics.par$usr[1:2]
  conDisp <- fileFolder[[current]]$graphics.par$conDisp
  maxRes <- maxShift( ucsf2D( file.name = currentSpectrum, 
                              w1Range = w1Range, w2Range = w2Range, 
                              file.par = fileFolder[[current]]$file.par ), conDisp = conDisp, 
                      massCenter = massCenter )
  if(is.null(maxRes))
    return(invisible())
  
  ## Set new range
  if( fileFolder[[current]]$file.par$number_dimensions != 1 )
    w1Range <- c(maxRes$w1 + abs(diff(w1Range))/2,  
                 maxRes$w1 - abs(diff(w1Range))/2  )
  w2Range <- c(maxRes$w2 + abs(diff(w2Range))/2,  
               maxRes$w2 - abs(diff(w2Range))/2  ) 
  
  ## Print warning or reset the window
  fileFolder[[current]]$graphics.par$usr <- c(w2Range, w1Range ) 
  myAssign('zoom', fileFolder)
  
  ## Draw the new spectrum                    
  refresh(sub.plot = FALSE, multi.plot = FALSE)
}    

## Changes contours to a higher level
## n  - Number of standard deviations to raise the plotting threshold
## ...  - Additional plotting options can be passed to drawPeptides and par()
ctu <- function (n = 1, ...) {
  
  ## Update contour levels and recalculate the viewable tiles
  current <- wc()
  uClevel <- fileFolder[[ current ]]$graphics.par$clevel + n
  
  ## Reset the graphics and refresh the open plots 
  setGraphics(clevel = uClevel)
  refresh(...)  
  
}

## Changes contours to a lower level
## n  - Number of standard deviations to decrease plotting threshold
## ...  - Additional plotting options can be passed to drawPeptides and par()
ctd <- function (n = 1, ...) {
  ctu(n = -n, ...)  
}

## Find the 2D location of the pointer
## Returns chemical shifts rounded to four places
loc <- function(){
  
  current <- wc()
  lineCol <- fileFolder[[current]]$graphics.par$fg
  
  ## Opens the main plot window if not currently opened
  if (is.na(match(2, dev.list())))
    refresh(multi.plot=FALSE, sub.plot=FALSE)
  cw(dev=2)
  
  ## Give the user some instructions
  cat(paste('\nIn the main plot window:\n',  
            ' Left-click to view chemical shifts\n',  
            ' Right-click to exit\n\n')) 
  flush.console()
  hideGui()
  
  ##Print chemical shifts in the console 
  if(fileFolder[[current]]$file.par$number_dimensions == 1)
    locName <- c(fileFolder[[current]]$file.par$nucleus[1], 'Intensity')
  else    
    locName <- c(fileFolder[[current]]$file.par$nucleus[2], 
                 fileFolder[[current]]$file.par$nucleus[1]) 
  
  i <- 1
  out <- NULL
  while(TRUE){
    
    ## Tell the user what to do
    op <- par('font')
    par( font = 2 )
    legend("topleft", c('LEFT CLICK FOR SHIFTS','RIGHT CLICK TO EXIT'), 
           pch=NULL, bty='n', text.col=lineCol)
    par(font = op)
    
    xy = data.frame(locator(1, type='n'))
    if(length(xy) == 0)		
      break()
    if( i > 1 )
      refresh( sub.plot = FALSE, multi.plot = FALSE)
    abline( v = xy$x, h = xy$y )
    
    if( i == 1)
      cat('\n', locName, '\n')
    
    cat(unlist(round(xy, 4)), '\n')
    i <- 2
    out <- xy
    
    flush.console()
  }
  
  #return focus to console and print	
  showGui()
  if(!is.null(out))
    names(out) <- locName
  refresh( sub.plot = FALSE, multi.plot = FALSE)
  invisible(out)
}

## Find the chemical shift range in Hz and PPM between user-defined points
## Returns the difference between points in Hz and PPM rounded to four places
delta <- function(){
  
  current <- wc()
  lineCol <- fileFolder[[current]]$graphics.par$fg
  
  ## Opens the main plot window if not currently opened
  if (is.na(match(2, dev.list())))
    refresh(multi.plot=FALSE, sub.plot=FALSE)
  cw(dev=2)
  
  ## Give the user some instructions
  cat(paste('\nIn the main plot window:\n',  
            ' Left-click on two points\n',  
            ' Right-click to exit\n\n')) 
  flush.console()
  hideGui()
  
  ##Print chemical shifts in the console
  nDim <- fileFolder[[current]]$file.par$number_dimensions
  if( nDim == 1 ){
    locName <- c(fileFolder[[current]]$file.par$nucleus[1], 'Intensity')
    
    
  }else{
    locName <- c(fileFolder[[current]]$file.par$nucleus[2], 
                 fileFolder[[current]]$file.par$nucleus[1]) 		
  }    
  
  out <- NULL
  i <- j <- 1
  while(TRUE){
    
    ## Tell the user what to do
    op <- par('font')
    par( font = 2 )
    legend("topleft", c('LEFT CLICK FOR DELTA','RIGHT CLICK TO EXIT'), 
           pch=NULL, bty='n', text.col=lineCol)
    par(font = op)
    
    ## Let user select locations on spectrum
    xy <- data.frame(locator(1, type='n'))
    if(length(xy) == 0)		
      break()
    if( i == 1 && j > 1)
      refresh( sub.plot = FALSE, multi.plot = FALSE)
    
    ## Show selected points with ablines
    abline( v = xy$x, h = xy$y)
    j <- 2
    
    ## Create the output table for printing chemical shift deltas
    if( i == 1 ){
      out <- xy
      i <- 2
      next()
    }
    
    out <- rbind(out, xy)
    names(out) <- locName	
    
    out <- round(rbind( abs(out[1,]-out[2,]), abs(out[1,]-out[2,]) * 
                          rev(fileFolder[[wc()]]$file.par$transmitter_MHz) ), 4)
    row.names(out) <- c('PPM', 'Hz')
    if( nDim == 1 )
      out[2,2] <- NA
    log_message(out)
    flush.console()
    i <- 1
  }
  
  #return focus to console 	
  showGui()
  if(!is.null(out))
    names(out) <- locName
  refresh( sub.plot = FALSE, multi.plot = FALSE)
  invisible(out)
}

################################################################################
##                                                                            ##
##                   1D projection functions for users                        ##
##                                                                            ##
################################################################################

## View 1D slice of 2D spectrum
## proj.direct  - Integer of value 1 or 2 indicating the direction of the slice;
##                1 returns direct slices, 2 returns indirect slices
## ...  - Additional plotting options can be passed to draw1D, par, and proj1D
vs <- function (proj.direct = NULL, ...) {
  
  ## Stop if current spectrum is not 2D
  current <- wc()
  lineCol <- fileFolder[[current]]$graphics.par$fg
  in.folder <- fileFolder[[current]]
  if( in.folder$file.par$number_dimensions < 2 )
    err("Only 2D data can be projected into 1D")
  
  ## Have user define slice dimension 
  if (!any(proj.direct == c(1,2))){
    usrSel <- mySelect(c('Direct dimension', 'Indirect dimension'), 
                       title = 'View slices in:')
    if( length(usrSel) == 0 || !nzchar(usrSel) )
      return(invisible())
    if(usrSel == 'Direct dimension')
      proj.direct <- 1
    else
      proj.direct <- 2
  }
  
  ## Opens the main plot window if not currently opened
  if (is.na(match(2, dev.list())))
    refresh(multi.plot=FALSE, sub.plot=FALSE)
  cw(dev=2)
  
  ## Give the user some instructions
  cat(paste('In the main plot window:\n',  
            ' Left-click inside the plot to view slice\n',  
            ' Right-click to exit\n'))
  flush.console()
  hideGui()
  
  ## Generate slices
  xy <- outFile <- NULL
  
  while( TRUE ){
    
    ## Tell the user what to do
    op <- par('font')
    par( font = 2 )
    legend("topleft", c('LEFT CLICK FOR SLICE','RIGHT CLICK TO EXIT'), 
           pch=NULL, bty='n', text.col=lineCol)
    par(font = op)
    
    xy <- locator(1)
    if(is.null(xy))
      break()
    refresh( sub.plot = FALSE, multi.plot = FALSE)
    outFile <- proj1D( in.folder, filter = 0, 
                       proj.direct = proj.direct, xy = xy, ...)
    
  }
  showGui()
  refresh( sub.plot = FALSE, multi.plot = FALSE)
  invisible(outFile)
}

## User function for toggling the 1D projection display
pjv <- function(){
  if (globalSettings$proj.mode){
    setGraphics(proj.mode=FALSE)
    refresh(multi.plot=FALSE, sub.plot=FALSE)
    log_message('1D projection display off', quote=FALSE)
  }else{
    setGraphics(proj.mode=TRUE)
    refresh(multi.plot=FALSE, sub.plot=FALSE)
    log_message('1D projection display on', quote=FALSE)
  }
}

################################################################################
##                                                                            ##
##               Peak picking functions for users                             ##
##                                                                            ##
################################################################################

## Peak picks full sweep width of current spectrum
pa <- function(...){
  
  ## Error checking
  wc()
  
  ## Peak pick active window    
  peakPick(append = FALSE, ...)
  
  ## Refresh graphics and save a backup
  setGraphics(peak.disp = TRUE,  save.backup = FALSE)
  refresh(sub.plot = FALSE, multi.plot = FALSE)
  
}

## Peak picks full sweep width in all of the spectra in the file folder 
paAll <- function(...){
  
  ## Error checking
  wc()
  
  ## Peak pick each spectrum
  peakPick(fileName = names(fileFolder), append = FALSE, ...)
  
  ## Refresh graphics and save backup copy  
  setGraphics(peak.disp = TRUE,  save.backup = FALSE)
  refresh(sub.plot = FALSE, multi.plot = FALSE)
  
}

## Peak pick a region in the current spectrum
## fileName - character string or vector, names for the files to peak pick
## append  - logical, TRUE apppends new peaks to old list
pReg <- function(fileName = currentSpectrum, append = TRUE, ...){
  
  ## Define the current spectrum
  current <- wc()
  nDim <- fileFolder[[current]]$file.par$number_dimensions
  lineCol <- fileFolder[[current]]$graphics.par$fg
  
  ## Opens the main plot window if not currently opened
  if (is.na(match(2, dev.list())))
    refresh(multi.plot=FALSE, sub.plot=FALSE)
  cw(dev=2)
  
  ## Give the user some instructions
  hideGui()
  cat(paste('In the main plot window:\n',  
            ' Left-click two points inside the plot to define region\n'))
  flush.console()
  op <- par('font')
  par(font=2)
  legend("topleft", c('LEFT CLICK TO DEFINE REGION', 'RIGHT CLICK TO CANCEL'), 
         pch=NULL, bty='n', text.col=lineCol)
  par(font=op)
  
  ##define the first boundary for the region
  xy <- data.frame(locator(1))
  if (length(xy) == 0){
    refresh(multi.plot=FALSE, sub.plot=FALSE)
    showGui()
    return(invisible())
  }
  abline(v=xy$x, col=lineCol )
  if (nDim != 1)
    abline(h=xy$y, col=lineCol )    
  
  ##define other boundary for the region
  xy2 <- data.frame(locator(1))
  if (length(xy2) == 0){
    refresh(multi.plot=FALSE, sub.plot=FALSE)
    showGui()
    return(invisible())
  }
  abline(v=xy2$x, col=lineCol )
  if (nDim != 1)
    abline(h=xy2$y, col=lineCol )   
  xy <- rbind(xy, xy2)
  
  ## Peak pick region
  if (nDim == 1)
    peakPick(fileName, w2Range=c(rev(sort(xy$x))), append=append, ...)
  else
    peakPick(fileName, w1Range=c(rev(sort(xy$y))), 
             w2Range=c(rev(sort(xy$x))), append=append, ...)	
  
  ## Refresh graphics
  setGraphics(peak.disp=TRUE, save.backup=FALSE)
  refresh(sub.plot = FALSE, multi.plot = FALSE)
  showGui()
}

## User function regionMax
## fileName - character string or vector; spectrum name(s) as returned by 
##						names(fileFolder)
## redraw - logical, TRUE refreshes the spectrum before exiting 
## noiseCheck - logical, TRUE excludes data below the noise threshold
## Returns chemical shifts at absolute max intensity in a user-defined region 
regionMax <- function( fileName=currentSpectrum, redraw=TRUE, noiseCheck=TRUE ){
  
  ## Define the current spectrum
  current <- wc()
  nDim <- fileFolder[[current]]$file.par$number_dimensions
  lineCol <- fileFolder[[current]]$graphics.par$fg
  
  ## Opens the main plot window if not currently opened
  if (is.na(match(2, dev.list())))
    refresh(multi.plot=FALSE, sub.plot=FALSE)
  cw(dev=2)
  
  ## Give the user some instructions
  hideGui()
  cat(paste('In the main plot window:\n',  
            ' Left-click two points inside the plot to define region\n'))
  flush.console()
  op <- par('font')
  par(font=2)
  legend("topleft", c('LEFT CLICK TO DEFINE REGION', 'RIGHT CLICK TO CANCEL'), 
         pch=NULL, bty='n', text.col=lineCol)
  par(font=op)
  
  ##define the first boundary for the region
  xy <- data.frame(locator(1))
  if (length(xy) == 0){
    refresh(multi.plot=FALSE, sub.plot=FALSE)
    showGui()
    return(invisible())
  }
  abline(v=xy$x, col=lineCol )
  if (nDim != 1)
    abline(h=xy$y, col=lineCol )    
  
  ##define other boundary for the region
  xy2 <- data.frame(locator(1))
  if (length(xy2) == 0){
    refresh(multi.plot=FALSE, sub.plot=FALSE)
    showGui()
    return(invisible())
  }
  abline(v=xy2$x, col=lineCol )
  if (nDim != 1)
    abline(h=xy2$y, col=lineCol )   
  xy <- rbind(xy, xy2)
  showGui()
  
  ## Find max observable signal for each spectrum 
  outFile <- NULL
  outNames <- NULL
  for( i in fileName ){
    
    ## Define variables for the current spectrum
    cNum <- match(i, names(fileFolder))
    cFile <- fileFolder[[cNum]]
    cDim <- cFile$file.par$number_dimensions
    outNames <- c(outNames, i)
    cRange <- matchShift(cFile, w2 = xy$x, w1 = xy$y, w1.pad = 1, w2.pad = 1, 
                         overRange = TRUE)
    
    ## Skip spectra outside the user-defined chemical shift ranges
    if( any(is.na(unlist(cRange))) ){
      cData <- rep(NA, 3)
      names(cData) <- c('w1', 'w2', 'Height')
      outFile <- rbind(outFile, cData ) 
      next
    }
    
    ## Find the absolute maximum point in the range
    cData <- maxShift(inFile=ucsf2D(file.name=i, w1Range=cRange$w1, 
                                    w2Range=cRange$w2, file.par=cFile$file.par),
                      conDisp=cFile$graphics.par$conDisp)
    
    ## Exclude data outside the spectral window
    if( is.null(cData) ){
      cData <- rep(NA, 3)
      names(cData) <- c('w1', 'w2', 'Height')
      outFile <- rbind(outFile, cData ) 
      next
    }
    
    ## Exclude data below the noise threshold
    if(noiseCheck){
      if( cDim > 1 && abs(cData$Height) < 
          cFile$file.par$noise_est * cFile$graphics.par$clevel )
        cData <- rep(NA, 3)
      if( cDim == 1 && abs(cData$Height) < 
          globalSettings$thresh.1D * cFile$file.par$noise_est + 
          cFile$file.par$zero_offset )
        cData <- rep(NA, 3)
      names(cData) <- c('w1', 'w2', 'Height')
    }
    outFile <- rbind(outFile, cData ) 
  }
  
  ## Format output data
  outFile <- suppressWarnings(data.frame(outFile))
  row.names(outFile) <- outNames
  
  ## Refresh the graphics
  if(redraw)
    refresh(multi.plot=FALSE, sub.plot=FALSE)
  cInc <- match(names(fileFolder[current]), row.names(outFile))
  abline(v = outFile$w2[cInc], h = outFile$w1[cInc], lty = 2 )
  
  return(outFile)
}

## Peak pick the maximum intensity within a region
## fileName - character string or vector; spectrum name(s) as returned by 
##						names(fileFolder)
## append - logical; TRUE apppends new peaks to old list
pm2 <- function(fileName=currentSpectrum, append = TRUE, ...){
  log_message('pm() digestR_backend')
  pReg(fileName = fileName, append = append, maxOnly = TRUE, ...)
}

## Peak pick current window in current spectrum
## append  - logical argument, TRUE apppends new peaks to old list
pw <- function(append = TRUE, ...){
  
  ## Define the current spectrum
  current <- wc()
  w1Range <- fileFolder[[current]]$graphics.par$usr[3:4]
  w2Range <- fileFolder[[current]]$graphics.par$usr[1:2]
  
  
  ## Peak pick active window
  peakPick(w1Range = w1Range, w2Range= w2Range, append = append, ...)
  setGraphics(peak.disp = TRUE,  save.backup = FALSE)
  if( !append )
    refresh(sub.plot = FALSE, multi.plot = FALSE)
  else
    pdisp()
}

## Peak pick current window in all spectra
## append  - logical argument, TRUE apppends new peaks to old list
pwAll <- function( append = TRUE, ... ){
  
  ##Checks that spectra have the same nuclei on the same axis
  wc()
  fileName <- names(fileFolder)
  if (length(fileName) > 1){
    for (i in fileName){	
      if (fileFolder[[i]]$file.par$number_dimensions != 
          fileFolder[[1]]$file.par$number_dimensions){
        if (fileFolder[[i]]$file.par$number_dimensions == 1){
          if (match(fileFolder[[i]]$file.par$nucleus,
                    fileFolder[[1]]$file.par$nucleus) != 2)
            stop('All spectra must have the same nuclei, on the same axis', 
                 quote=FALSE)
        }else if (match(fileFolder[[i]]$file.par$nucleus,
                        fileFolder[[1]]$file.par$nucleus)[2] != 1)
          stop('All spectra must have the same nuclei, on the same axis', 
               quote=FALSE)
      }else if (!all(fileFolder[[i]]$file.par$nucleus == 
                     fileFolder[[1]]$file.par$nucleus))
        stop('All spectra must have the same nuclei, on the same axis', 
             quote=FALSE)
    }
  }
  
  ## Define the current spectrum
  current <- wc()
  w1Range <- fileFolder[[current]]$graphics.par$usr[3:4]
  w2Range <- fileFolder[[current]]$graphics.par$usr[1:2]
  
  ## Peak pick each spectrum
  peakPick(fileName = names(fileFolder), w1Range = w1Range, 
           w2Range = w2Range, append = append, ...)
  
  ## Refresh graphics and save backup copy  
  setGraphics(peak.disp = TRUE,  save.backup = FALSE)
  if( !append )
    refresh(sub.plot = FALSE, multi.plot = FALSE)
  else
    pdisp()
}

## User peak picking function rp
## Peak picks inside or outside ROIs 
## fileName - string or character vector; spectrum name(s) as returned by 
##						names(fileFolder)
## append  - logical argument, TRUE apppends new peaks to old list
## ...     - arguments can be passed to internal peak picking functions
## Saves the new peak list to the file folder and displays new peaks
rp <- function( fileName = currentSpectrum, append = TRUE, parent=NULL, ...){
  
  ## Error checking
  current <- wc()
  if(!exists('roiTable') || is.null(roiTable) || nrow(roiTable) == 0)
    err('No ROIs have been designated, use roi()')
  
  ## Prompt user for ranges to be picked  
  usrSel <- mySelect(c('Inside ROIs', 'Maxima of ROIs', 'Outside ROIs'), 
                     multiple=FALSE, preselect='Inside', title='Peak pick:', parent=parent)
  if ( length(usrSel) == 0 || !nzchar(usrSel) )
    return(invisible())
  
  ## Have user select the type of peak picking	
  if( usrSel != 'Outside ROIs' ){
    
    ## Have user select rois to peakPick 
    usrROI <- mySelect(roiTable$Name, multiple = TRUE, index = TRUE,
                       preselect = roiTable$Name[which(roiTable$ACTIVE == TRUE)], 
                       title = 'Select ROIs to peak pick', parent=parent)
    if ( length(usrROI) == 0 || !nzchar(usrROI) )
      return(invisible())
    usrROI <- roiTable[usrROI,]
    
  }else
    usrROI <- roiTable
  
  ## Peak pick each spectrum
  for(j in fileName ){
    
    ## Set parameters
    cNum <- which( names(fileFolder) == j )
    oList <- fileFolder[[cNum]]$peak.list
    cDim <- fileFolder[[cNum]]$file.par$number_dimensions
    
    ## Peak pick entire spectrum
    tList <- peakPick( fileName = j, internal = TRUE, ...)
    
    ## Find subset of peaks inside/outside ROIs
    out <- NULL
    for( i in 1:nrow(usrROI) ){		
      
      ## Find w2 subset
      subList <- tList[tList$w2 < usrROI$w2_downfield[i], ]
      subList <- subList[subList$w2 > usrROI$w2_upfield[i], ]
      if( usrROI$nDim[i] == 1 || cDim == 1 ){
        if (usrSel == 'Maxima of ROIs')
          subList <- subList[which.max(abs(subList$Height)), ]
        out <- rbind(out, subList)
        next()
      }
      
      ## Find w1 subset		
      subList <- subList[subList$w1 < usrROI$w1_downfield[i], ]
      subList <- subList[subList$w1 > usrROI$w1_upfield[i], ]
      if (usrSel == 'Maxima of ROIs')
        subList <- subList[which.max(abs(subList$Height)), ]	
      out <- rbind(out, subList )
    }
    
    ## Invert peak selection for outside roi option was selected
    if( usrSel == 'Outside ROIs' )
      nList <- tList[which(is.na(match(tList$Index, out$Index))),]
    else
      nList <- out
    
    ## Tell user which file was peak picked
    cat(paste(basename(j), '\n'))
    
    ## update appended lists
    if( append ){
      
      if( length(oList) && nrow(nList) ){
        nList <- appendPeak( newList = nList, oldList = oList )
        cat(paste(' Total peaks:', nrow(nList), '\n', 
                  'New peaks:', nrow(nList) - nrow(oList), '\n'))				
      }
      
      if( !length(oList) && nrow(nList) ){
        row.names(nList) <- NULL
        nList$Index <- 1:nrow(nList)
        cat(paste(' Total peaks:', nrow(nList), '\n', 
                  'New peaks:', nrow(nList), '\n'))	
      }
      
      if( length(oList) && !nrow(nList) ){
        nList <- oList
        cat(paste(' Total peaks:', nrow(nList), '\n', 
                  'New peaks:', 0, '\n'))	
      }
      
      if( !length(oList) && !nrow(nList) ){
        cat(paste(' Total peaks:', 0, '\n'))
        nList <- NULL	
      }
      
      ## update non-appended lists
    }else{
      if( nrow(nList) ){
        cat(paste(' Total peaks:', nrow(nList), '\n'))
        row.names(nList) <- NULL
        nList$Index <- 1:nrow(nList)
      }else{
        cat(paste(' Total peaks:', 0, '\n'))
        nList <- NULL
      }
    }
    
    ## Save new peak list to file folder
    fileFolder[[cNum]]$peak.list <- nList
    flush.console()
  }
  
  ## Assign file folder and refresh the graphics
  myAssign("fileFolder", fileFolder, save.backup = TRUE)
  setGraphics(peak.disp = TRUE,  save.backup = FALSE)
  if( !append )
    refresh(sub.plot = FALSE, multi.plot = FALSE)
  else
    pdisp()
}



## User function rpAll
## Peak pick ROIs in all files
## append  - logical argument, TRUE apppends new peaks to old list
## ...     - arguments can be passed to internal peak picking functions
rpAll <- function( append = TRUE, ...){
  
  wc()
  rp( fileName = names(fileFolder), append = append, ... )
  
}

## User peak function
## Allow users to pick peaks by hand
## forcePoint - Logical argument, when TRUE, chemical shifts and intensity
##              will be taken from the closest data point in the spectrum;
##              when FALSE the intensity is taken from the closest data point
##              but chemical shifts will not be corrected to match the spectrum
## Note: digestR does not interpolate between data points. The default behavior
##       of this function is to allows users to specify any chemical shift, but 
##       the intensities of these locations are derived from the closest 
##       neighboring point. The most reliable method for finding the maximum
##       intensity of a peak is to use ROI summary.
## returns a list of the unique hand picked peaks and appends these data
##       to the peak list for the active spectrum.
ph <- function( forcePoint = FALSE){
  
  ## Check to make sure the requisite files are present
  current <- wc()
  lineCol <- fileFolder[[current]]$graphics.par$fg
  
  ## Bring main into focus and show any active peaks
  if(length(which(dev.list() == 2)) != 1)
    drawPeptides()
  cw(dev=2);	pdisp()
  
  ## Give the user some instructions
  cat(paste('In the main plot window:\n',  
            ' Left-click inside the plot to peak pick\n',  
            ' Right-click to exit\n'))
  flush.console()
  hideGui()
  
  ## Have user select peaks
  log_message('New peaks:', quote = FALSE)
  outList <- NULL
  while( TRUE ){
    
    ## Tell the user what to do
    op <- par('font')
    par( font = 2 )
    legend("topleft", c('LEFT CLICK TO PICK','RIGHT CLICK TO EXIT'), 
           pch=NULL, bty='n', text.col=lineCol)
    par(font = op)
    
    xy <- locator(1)
    if(is.null(xy))
      break()
    
    ## Force chemical shifts to the closest datapoint
    if( forcePoint ){
      xy <- matchShift( w1 = xy$y, w2 =xy$x )
      peakFile <- ucsf2D( currentSpectrum, w1Range = xy$w1, 
                          w2Range = xy$w2, file.par = fileFolder[[current]]$file.par )
      
      ## Do not force chemical shifts, Height is taken from closest datapoint	
    }else{
      peakFile <- list( w1 = xy$y, w2 = xy$x )
      xy <- matchShift( w1 = xy$y, w2 =xy$x )
      peakFile$data <- ucsf2D(currentSpectrum, w1Range = xy$w1, 
                              w2Range = xy$w2, file.par = fileFolder[[current]]$file.par )$data
    }
    
    
    ## Make a new peak list or append to the old list
    if( fileFolder[[current]]$file.par$number_dimensions  > 1 )
      nPeak <- data.frame( list(Index = 1,  w1 = peakFile$w1, w2 = peakFile$w2, 
                                Height = peakFile$data, Assignment = "NA" ), 
                           stringsAsFactors = FALSE)
    else
      nPeak <- data.frame(list( Index = 1, w1 = NA, w2 = peakFile$w2, 
                                Height = peakFile$data, 
                                Assignment = "NA" ), stringsAsFactors = FALSE)
    
    fileFolder[[current]]$peak.list <- appendPeak(nPeak, 
                                                  fileFolder[[current]]$peak.list)
    
    ## Assign new peak list
    log_message(nPeak[2:4], quote=FALSE)
    flush.console()
    outList <- unique (rbind( outList, nPeak ))
    myAssign("fileFolder", fileFolder, save.backup = FALSE)	
    pdisp()
  }
  
  ## Assign the final version of the file and print the new peaks
  showGui()
  setGraphics(peak.disp = TRUE,  save.backup = TRUE)
  refresh( sub.plot = FALSE, multi.plot = FALSE)
  row.names( outList ) <- NULL
  invisible( outList ) 
}

## Turn on/off peak display
pv <- function(){
  
  ## Check to make sure the requisite files are present
  wc()
  
  ## Turn peak view on/off
  if( globalSettings$peak.disp ){
    setGraphics( peak.disp = FALSE )		
    refresh(sub.plot = FALSE, multi.plot = FALSE)
    cat('Peak display off \n')
  }else{
    setGraphics( peak.disp = TRUE )   
    pdisp()	
    cat('Peak display on \n')
  }
}

## Clear the peak list for the current spectrum 
pDel <- function(){
  current <- wc()
  if ( is.null(fileFolder[[current]]$peak.list) )
    log_message('The current peak list is empty', quote=FALSE)
  else
    peakDel()	
}

## Internal function for peak delete
## fileName - A list of file names from fileFolder
peakDel <- function( fileName = currentSpectrum ){
  for (i in fileName)
    fileFolder[[i]]$peak.list <- NULL
  myAssign('fileFolder', fileFolder)
  refresh(sub.plot=FALSE, multi.plot=FALSE)
}

## Clear peak list for all open spectra 
pDelAll <- function(){
  
  ## Check to make sure the requisite files are present
  current <- wc()
  peak.list <- fileFolder[[current]]$peak.list
  
  ## Clear peak lists and turn off peak display
  for(i in 1:length(fileFolder))
    fileFolder[[i]]$peak.list <- NULL   
  
  myAssign("fileFolder", fileFolder, save.backup = FALSE)
  if( nrow(peak.list) > 0 )
    refresh(sub.plot = FALSE, multi.plot = FALSE) 
}

## Open Madison Metabolomics Consortium Database in a web browser
mmcd <- function(){browseURL("http://mmcd.nmrfam.wisc.edu")}

## Changes the noise filter for peak picking
## filt - the desired filter level; either 0, 1, or 2 for off, mild, and strong
nf <- function(noiseFilt){
  if (missing(noiseFilt)){
    filtLevel <- switch(globalSettings$peak.noiseFilt + 1, 'off', 'mild', 
                        'strong')
    cat(paste('Current noise filter setting:', filtLevel, '\n' ))
    return(invisible())	
  }
  
  if (!any(noiseFilt == c(0, 1, 2)))
    stop('Invalid "noiseFilt" argument, use: 0 (off), 1 (mild) or 2 (strong)')
  
  filtLevel <- switch(noiseFilt + 1, 'off', 'mild', 'strong')
  setGraphics(peak.noiseFilt = noiseFilt)
  cat(paste('Peak picking noise filter has been set to:', filtLevel, '\n'))
}

################################################################################
##                                                                            ##
##               Defining and viewing Regions of interst                      ##
##                                                                            ##
################################################################################

## Internal utility function maxShift
## Finds chemical shifts of the maximum absolute intensity in the file folder
## inFile  - file parameters and data for desired spectrum as returned by ed()
## invert       - Logical argument, inverts the data before finding the max
##                shift. This is used in the graphics functions
## conDisp      - Logical vector of length 2; c(TRUE, TRUE) returns the 
##                absolute max, c(TRUE, FALSE) returns the max, c(TRUE, FALSE)
##                returns the min
## returns chemical shifts and intensity of the maximum intensity signal
maxShift <- function( inFile, invert = FALSE, conDisp = c(TRUE, TRUE), 
                      massCenter = FALSE ){
  
  
  ## Set function
  if( inFile$file.par$number_dimensions == 1 || all(conDisp )){
    if( invert )
      FUN <- function(x){ which.max(rev(abs(x)))}
    else
      FUN <- function(x){ which.max(abs(x))}
    
  }else{
    if( conDisp[1] ){
      if( invert )
        FUN <- function(x){ which.max(rev(x))}
      else
        FUN <- function(x){ which.max((x))}
      
    }else{
      if( invert )
        FUN <- function(x){ which.min(rev(x))}	
      else
        FUN <- function(x){ which.min((x))}		
      
    }	
  }
  
  if(inFile$file.par$number_dimensions > 1){ 
    if( !is.matrix(inFile$data) )
      return(NULL)
    
    ## Find the location of the absolute max intensity in the roi 
    bs <- nrow(inFile$data)
    sMax <- FUN(inFile$data) - 1
    rNum <- sMax %% bs 
    cNum <- sMax %/% bs 
    w1 <- rev(inFile$w1)[cNum + 1]
    w2 <- rev(inFile$w2)[rNum + 1]
    Height <- inFile$data[sMax + 1]
    
  }else{
    if( length(inFile$data) < 2 )
      return(NULL)
    
    ## Center w2 data by smoothed spectrum
    if( massCenter ){
      rNum <- try( FUN(smooth.spline(inFile$data, df = 5)$y ), silent = TRUE )
      if(class(rNum) == "try-error")
        rNum <- FUN(inFile$data)
    }else
      rNum <- FUN(inFile$data)
    w1 <- NA
    w2 <- rev(inFile$w2)[rNum]
    Height <- inFile$data[rNum]
    
  }
  return(data.frame(w1, w2, Height))
} 

## Internal roi utility function orderROI
## Removes any duplicates entries, renumbers ROIs, and sorts by ROI name
## roiTable  - An roi table object
## Returns an updated roi table
orderROI <- function( roiTable = NULL ){
  
  if( is.null(roiTable) )
    stop('No ROI table was entered')
  if( length(nrow(roiTable)) == 0 )
    return(roiTable)
  
  ## Remove duplicate ROIs
  tName <- try(as.vector(sapply(roiTable[,1], 
                                function(x){unlist(strsplit(x, ".", fixed = TRUE))[1]})),
               silent = TRUE)
  if(class(tName) == "try-error")
    tName <- roiTable$Name
  roiTable$Name <- tName
  roiTable <- unique(roiTable) 	
  
  ## Number any duplacate roi names
  for( i in 1: length(roiTable$Name) ){
    tName <- which(roiTable$Name == roiTable$Name[i])
    if( roiTable$Name[i] == 'ROI' || length (tName) > 1 )
      roiTable$Name[tName] <- paste(roiTable$Name[i], 
                                    1:length(tName), sep ='.' )
  }	
  row.names(roiTable) <- NULL	
  
  return(roiTable)
}

## Internal graphics function showRoi
## Plots boxes (2D) or lines (1D) to indicate the location of ROIs
## rTable  - roiTable objects can be passed from other functions
## col     - line color vector [active, inactive]
## text.col	- text color vector [active, inactive]
## lw      - line width vector [active, inactive]
## lty     - line type vector [active, inactive] (sudee par)
## cex     - Text size (see par)
## Note: this function is used to show ROIs on the main plot
showRoi <- function ( rTable = roiTable, col = globalSettings$roi.bcolor, 
                      text.col = globalSettings$roi.tcolor, lw = globalSettings$roi.lwd, 
                      lty = globalSettings$roi.lty, cex = globalSettings$roi.cex) {
  
  ## Check to make sure the requisite files are present
  if( !exists('roiTable') && is.null(rTable) )
    return('The ROI table is empty, use roi() to make a new roi')
  if( is.null(rTable) )
    rTable <- roiTable
  if( is.null(rTable$Name) || nrow(rTable) == 0 )
    return('The ROI table is empty, use roi() to make a new roi')	
  
  ## Set ROI label alignment
  if (globalSettings$roi.labelPos == 'top')
    pos <- 3
  else if (globalSettings$roi.labelPos == 'bottom')
    pos <- 1
  else if (globalSettings$roi.labelPos == 'left')
    pos <- 2
  else if (globalSettings$roi.labelPos == 'right')
    pos <- 4
  else
    pos <- NULL
  
  ## Define the current spectrum
  current <- wc()   
  nDim <- fileFolder[[current]]$file.par$number_dimensions    
  current.par <- fileFolder[[current]]$file.par 
  usr <- fileFolder[[current]]$graphics.par$usr
  if (is.na(match(2, dev.list())))
    refresh(sub.plot=FALSE, multi.plot=FALSE)
  
  for(rowNum in 1:nrow(rTable)){
    ##Find roi PPM locations
    w2Range <- sort(unlist(rTable[rowNum, 2:3]))
    w1Range <- sort(unlist(rTable[rowNum, 4:5]))
    
    ## Fudge the ROI table for 1D/2D compatibility
    if( rTable$nDim[rowNum] == 2 && nDim == 1)
      w1Range <- c(current.par$min_intensity, current.par$max_intensity)      
    if( rTable$nDim[rowNum] == 1 && nDim > 1)
      w1Range <- c(current.par$upfield_ppm[1], current.par$downfield_ppm[1]) 
    
    ## Set color parameters for ROIs
    if( rTable[rowNum, 6] )
      i <- 1
    else
      i <- 2
    
    ## Draw 1D rois as lines, 2D rois as boxes
    if(nDim > 1){
      rect(w2Range[1], w1Range[2], w2Range[2], w1Range[1], lwd = lw[i], 
           border = col[i], lty = lty[i], cex = cex[i])
      if (is.null(pos))
        textCoord <- c(mean(w2Range), mean(w1Range))
      else if (pos == 3)
        textCoord <- c(mean(w2Range), min(w1Range))
      else if (pos == 1)
        textCoord <- c(mean(w2Range), max(w1Range))
      else if (pos == 2)
        textCoord <- c(max(w2Range), mean(w1Range))
      else if (pos == 4)
        textCoord <- c(min(w2Range), mean(w1Range))
      else
        textCoord <- c(mean(w2Range), mean(w1Range))	
      text(textCoord[1], textCoord[2], rTable$Name[rowNum], cex = cex[i], 
           col = text.col[i], pos = pos, offset = .3)				
    }else{
      lines(x = rep( mean(w2Range), 2), y = c(w1Range[2], usr[4]* .75), 
            lty = lty[i], lwd = lw[i], col = col[i])
      lines(x = w2Range, y = rep(w1Range[2], 2), lty = lty[i], lwd = lw[i], 
            col = col[i])
      text(mean(w2Range), usr[4]*.88, rTable$Name[rowNum], cex = cex[i], 
           col = text.col[i], srt = 90, pos = pos, offset = .3)					
    }
  }  
}

## Plots all of the defined regions of interest (ROIs) in a new window
## ...  - Some plotting options can be passed to drawPeptides and par()
rvs <- function ( ... ){
  
  ## Define the current spectrum
  current <- wc()
  nDim <- fileFolder[[current]]$file.par$number_dimensions 
  current.par <- c(fileFolder[[current]]$file.par,  
                   fileFolder[[current]]$graphics.par) 
  
  ## Make a new window with some user instructions if no ROIs are active
  if( is.null(roiTable) || nrow(roiTable) == 0 ){
    setWindow('sub', bg = current.par$bg, fg = current.par$fg, 
              col.axis = current.par$col.axis, col.lab = current.par$col.lab,
              col.main = current.par$col.main, col.sub = current.par$col.sub,
              col = current.par$col, mfrow = c(1,1), pty = 'm', 
              mar = globalSettings$mar.sub , xpd = FALSE )
    plot(1, 1, type = 'n', axes=FALSE, xlab='', ylab='' )
    text(1, 1, 'ROIs will apear here \nuse roi() to designate an ROI', cex=1)
    dev.set(which = 2)  
    bringFocus(-1)
    return(invisible())
  }
  
  ## Set up a grid of roi sub plots for the roi plotting window   
  plot.grid <- c(ceiling(nrow(roiTable) / 10), 10)   
  setWindow('sub', bg=current.par$bg, fg=current.par$fg, 
            col.axis=current.par$col.axis, col.lab = current.par$col.lab, 
            col.main=current.par$col.main, col.sub =current.par$col.sub, 
            col=current.par$col, mfrow = plot.grid, pty = 'm', 
            mar = globalSettings$mar.sub, xpd = FALSE)
  
  
  ## Make sure that plots are possible
  answer <- try(plot(1,1, type = 'n', axes=FALSE, xlab='', ylab=''), 
                silent = TRUE)
  if(class(answer) == "try-error"){
    par( mfrow = c(1,1) )
    plot(1, 1, type = 'n', axes=FALSE, xlab='', ylab='' )
    text(1, 1, paste('Too many ROIs for this window:\nresize the window and ', 
                     'type refresh() or use rDel()', sep=''), cex=1)
    dev.set(which = 2)  
    bringFocus(-1)
    return(invisible())			
  }else
    par( mfrow = plot.grid )	
  
  
  ## Generate plots for each roi in a sub window of the roi plot window           
  for (i in 1:nrow(roiTable)){
    
    ##Find roi PPM locations
    w1Range <- unlist(roiTable[i, 4:5])
    w2Range <- unlist(roiTable[i, 2:3]) 
    
    ## Fudge the ROI table for 1D/2D compatibility
    if( roiTable$nDim[i] == 2 && nDim == 1)
      w1Range <- c(current.par$min_intensity, current.par$max_intensity)
    if( roiTable$nDim[i] == 1 && nDim > 1)
      w1Range <- c(current.par$upfield_ppm[1], current.par$downfield_ppm[1])
    
    
    ## Set the file parameters
    fileFolder[[current]]$graphics.par$usr <- c(w2Range, w1Range)    
    roi.name <- roiTable$Name[i]
    
    ##Plot the current roi
    if( dev.cur() != 3 ){
      dd()
      break()	
    }else{
      drawPeptides(in.folder = fileFolder[[current]],  main=roi.name, xlab='', 
                   ylab='', p.window = 'sub', axes = FALSE, 
                   cex = globalSettings$cex.roi.sub, ...)
      
      ##Draw active roi red
      if(roiTable[i, 6]){
        rect(par('usr')[1], par('usr')[4], par('usr')[2], par('usr')[3],  
             lwd=2, border='red')
      }
    }                   
  }
  
  dev.set(which = 2) 
  bringFocus(-1)
}

## Create a plot of ROIs from multiple files
## ...  - Some plotting options can be passed to drawPeptides and par()
rvm <- function ( ... ){
  
  ## Define the current spectrum
  current <- wc()  
  cPar <- c(fileFolder[[current]]$file.par, 
            fileFolder[[current]]$graphics.par) 
  
  ## Plot a black box if nothing is active
  if(length(which(roiTable[, 6] == TRUE)) == 0){
    setWindow('multi', bg = cPar$bg, fg = cPar$fg, 
              col.axis = cPar$col.axis, col.lab = cPar$col.lab,
              col.main = cPar$col.main, col.sub = cPar$col.sub,
              col = cPar$col, mfrow = c(1,1), pty = 'm', 
              mar = globalSettings$mar.multi, oma=c(.5, 0, 0, .5), xpd = FALSE )
    plot(1, 1, type = 'n', axes=FALSE, xlab='', ylab='')
    text(1, 1, 'Active ROIs will appear here \nuse rs()', cex=1)
    dev.set(which = 2)
    bringFocus(-1) 
    return(invisible())
  }   
  
  ## Find the active roi files
  actFile <- which(sapply(fileFolder, function(x){x$graphics.par$roi.multi}))
  
  ## Suppress rvm if no files are active     
  if( length(actFile) == 0 ){
    setWindow('multi', bg = cPar$bg, fg = cPar$fg, 
              col.axis = cPar$col.axis, col.lab = cPar$col.lab,
              col.main = cPar$col.main, col.sub = cPar$col.sub,
              col = cPar$col, mfrow = c(1,1), pty = 'm', 
              mar = globalSettings$mar.multi, oma=c(.5, 0, 0, .5),  xpd = FALSE) 
    plot(1, 1, type='n', axes=FALSE, xlab='', ylab='')
    text(1, 1, 'Active spectra will appear here \nuse rsf()', cex=1)
    dev.set(which = 2)
    bringFocus(-1)
    return(invisible())
  }
  
  
  ## Create grid for roi subplots   
  plot.grid =c(length(actFile) + 1, length(which(roiTable[, 6] == TRUE)) + 1)
  setWindow('multi', bg=cPar$bg, fg=cPar$fg, 
            col.axis=cPar$col.axis, col.lab = cPar$col.lab, 
            col.main=cPar$col.main, col.sub =cPar$col.sub, 
            col=cPar$col, mfrow = plot.grid, pty = "m", 
            mar = globalSettings$mar.multi, xpd=FALSE, oma=c(.5, 0, 0, .5)) 
  
  ## Make sure that plots are possible
  answer <- try(plot(1,1, type = 'n', axes=FALSE, xlab='', ylab=''), 
                silent = TRUE)
  if(class(answer) == "try-error"){
    par( mfrow = c(1,1) )
    plot(1, 1, type = 'n', axes=FALSE, xlab='', ylab='' )
    text(1, 1, paste('Too many ROIs for this window:\nresize the window and ',
                     'type refresh() or deactivate rois with rs()'), cex=1)
    dev.set(which = 2)  
    bringFocus(-1)
    return(invisible())			
  }else
    par( mfrow = plot.grid )	
  
  
  ## Print ROI names in first row
  for(j in  c(0, which(roiTable[, 6] == TRUE))){
    plot(0, 0, type='n', axes=FALSE, xlab='', ylab='')
    if(j == 0)
      next()
    text(0, -.5, roiTable$Name[j], cex=globalSettings$cex.roi.multi, 
         col = cPar$fg)            
  }    
  
  ## Plot active ROIs from each file
  for(i in 1:(length(actFile)) ){
    
    ## Print File name in first column
    if( dev.cur() != 4 ){
      dd()
      break()
    }else{
      plot(0, 0, type = 'n', axes=FALSE, xlab='', ylab='')
      text(0, 0, basename(names(actFile)[i]), xpd=FALSE, 
           cex=globalSettings$cex.files.multi, col = cPar$fg) 			
    }
    
    ## Plot the data from each ROI
    for(j in  which(roiTable[, 6] == TRUE)){
      w1Range <- as.numeric(roiTable[j, 4:5])
      w2Range <- as.numeric(roiTable[j, 2:3])
      rPar <- c(fileFolder[[actFile[i]]]$file.par,
                fileFolder[[actFile[i]]]$graphics.par)
      nDim <- rPar$number_dimensions
      
      ## Fudge the ROI table for 1D/2D compatibility
      if( roiTable$nDim[j] == 2 && nDim == 1)
        w1Range <- c(rPar$min_intensity, rPar$max_intensity)
      if( roiTable$nDim[j] == 1 && nDim > 1)
        w1Range <- c(rPar$downfield_ppm[1], rPar$upfield_ppm[1])
      
      ## Plot the roi
      if( dev.cur() != 4 ){
        dd()
        break()				
      }else
        drawPeptides (fileFolder[[actFile[i]]], w1Range = w1Range, 
                      w2Range = w2Range, fg = cPar$fg, col.axis = cPar$col.axis, 
                      xlab = '', ylab = '', main = '', axes = FALSE,  
                      p.window = 'multi', ...) 
    }         
  }
  #return focus to console  
  dev.set(which = 2) 
  bringFocus(-1) #return focus to console       
  
}

## Internal function changeRoi
## Modifies selected rois (use to shift, expand, contract etc.)
## ppmInc - logical argument, TRUE will interpret w1Inc and w2Inc in ppm,
##          FALSE will interpret w1Inc and w2Inc in number of points
## w1Inc - Numeric value expressing the ppm (ppmInc = TRUE) or percent change 
##         (ppmInc = FALSE) in an roi in the format c(downfield, upfield).
##          Positive increment values shift ROIs downfield (or up in 1D rois).
## Note:  - In 1D rois, only the second element of w1Inc is used 
## w2Inc - Numeric value expressing the ppm (ppmInc = TRUE) or percent change 
##         (ppmInc = FALSE) in an roi in the format c(downfield, upfield)
##          Positive increment values shift ROIs downfield
changeRoi <- function (ppmInc = FALSE, w1Inc = c(0, 0), w2Inc = c(0, 0), ...){
  
  ## Check incoming data
  if( is.null(roiTable) || nrow(roiTable) == 0)
    err('No ROIs have been designated, use roi()')	
  if( nrow(roiTable[roiTable$ACTIVE,]) == 0 )
    err( 'No ROIs have been selected' )
  if(length(c(w1Inc, w2Inc)) != 4  || !is.numeric(c(w1Inc, w2Inc)))
    err( 'Increments must be numeric and in the form c(downfield, upfield)')
  
  ## Define the current spectrum and find chemical shifts for each point
  actTable <- roiTable[roiTable$ACTIVE,]
  current <- wc()
  nDim <- fileFolder[[current]]$file.par$number_dimensions
  if (nDim == 3)
    nDim <- 2
  df <- fileFolder[[current]]$file.par$downfield_ppm
  uf <- fileFolder[[current]]$file.par$upfield_ppm
  ms <- fileFolder[[current]]$file.par$matrix_size
  
  ## Set default for cases where downfield/upfield roiTables are identical
  w1Min <- 0 * (df[1] - uf[1]) / (ms[1] - 1)
  if( nDim == 1 ){
    w2Min <- w1Min
    w1Min <- fileFolder[[current]]$file.par$noise_est
  }else
    w2Min <- 0 * (df[2] - uf[2]) / (ms[2] - 1)
  
  ## Update roi nDim to current spectrum if nDim does not match roiTable$nDim
  if(any(actTable$nDim != nDim)){
    usrInput <- myMsg(type="yesno", message = 
                        paste( "Do you want to convert the active ROIs to ", 
                               nDim, "D ROIs?", sep=''))
    if( usrInput == 'yes' ){
      
      ## Generate new list of w1 ranges
      if( nDim == 2 )
        newW1 <- c( df[1], uf[1] )
      else
        newW1 <- c(fileFolder[[current]]$file.par$min_intensity, 
                   fileFolder[[current]]$file.par$max_intensity )
      
      ## Update with the new ranges
      for(i in which(actTable$nDim != nDim))
        actTable[i,4:5] <- newW1
      actTable$nDim <- nDim
    }
  }
  
  ## Update w1/w2 ranges with ppm offsets
  if(ppmInc){
    actTable[ , 2] <- actTable[ , 2] + w2Inc[1]
    actTable[ , 3] <- actTable[ , 3] + w2Inc[2]
    for(i in which(actTable$nDim == 2)){
      actTable[ i , 4 ] <- actTable[i  , 4] + w1Inc[1]
      actTable[ i , 5] <- actTable[ i  , 5] + w1Inc[2]
    }
  }else{
    
    ## Update w2 ranges with percent offsets
    w1Inc <- w1Inc/50; w2Inc <- w2Inc/50
    for( i in 1:nrow(actTable) ){
      actTable[i,2:3] <- c(newRange(actTable[i, 2:3], f = w2Inc[1], ...)[2], 
                           newRange(actTable[i, 2:3], f = -w2Inc[2], ...)[1])
      if( actTable$nDim[i] != nDim )
        next()
      
      ## Update w1 ranges with percent offsets for 2D rois
      if (actTable$nDim[i] == 2){	
        actTable[i,4:5] <- c(newRange(actTable[i, 4:5], f = w1Inc[1], ...)[2],
                             newRange(actTable[i, 4:5], f = -w1Inc[2], ...)[1])
        
        ## Update w1 values for 1D rois in 1D files
      }else{
        ## This allows the user functions mru() and mrd() to behave properly
        if(w1Inc[1] == w1Inc[2])
          w1tmpInc <- -w1Inc
        else
          w1tmpInc <- -rev(w1Inc)
        
        ## Do not allow negative max and remove upfield/downfield difference 
        actTable[i, 5] <- actTable[i, 5] + (actTable[i, 5] * w1tmpInc[1])	
        if(actTable[i, 5] <= 0)
          actTable[i, 5] <- w1Min
      }															
    }
  }
  
  ## Assign new roi table and refresh the plots
  roiTable[match(row.names(actTable), row.names(roiTable)),] <- actTable
  myAssign('roiTable', roiTable)
  refresh() 
  
}

## Define region of interest from the active plot
## ...  - Additional plotting options can be passed to drawPeptides and par()
rn <- function( ... ){
  
  current <- wc() 
  lineCol <- fileFolder[[current]]$graphics.par$fg
  nDim <- fileFolder[[current]]$file.par$number_dimensions
  
  ##Open main and ROI subplot window if closed
  if(length(which(dev.list() == 2)) != 1)
    drawPeptides()
  
  ## Bring main into focus and show the rois 
  cw(dev=2);	showRoi()
  
  if(!exists('roiTable') || is.null(roiTable))
    roiTable <- list( Name = NULL, w2_downfield = NULL, w2_upfield = NULL,  
                      w1_downfield = NULL,  w1_upfield = NULL, ACTIVE = NULL, nDim = NULL)
  
  ## Give the user some instructions
  cat(paste('In the main plot window:\n',  
            ' Left-click two points inside the plot to designate an ROI\n',   
            ' Right-click to exit\n'))
  flush.console()
  hideGui()
  
  ## Have user define the plotting range by selecting a region 
  repeat{
    
    ## Tell the user what to do
    op <- par('font')
    par( font = 2 )
    legend("topleft", c('LEFT CLICK FOR ROI','RIGHT CLICK TO EXIT'), 
           pch=NULL, bty='n', text.col=lineCol)
    par(font = op)
    
    ## Have user define the plotting range by selecting a region 
    xy <- data.frame(locator(1, type = 'p', col = 'red', pch = 3))
    if(length(xy) == 0)
      break()   
    xy2 <- data.frame(locator(1, type = 'p', col = 'red', pch = 3))
    if(length(xy2) == 0)
      break()
    xy <- rbind(xy, xy2)
    xy$x <- sort(xy$x); xy$y <- sort(xy$y)
    
    if( nDim == 1 ){
      
      ## Show the new ROI and append it to the roi table
      abline( v = xy$x, col = 'red', lwd = 2 )
      newRoi <- list( Name = 'ROI', w2_downfield = xy$x[2], 
                      w2_upfield = xy$x[1], 
                      w1_downfield = fileFolder[[current]]$file.par$zero_offset - 
                        (xy$y[2] - fileFolder[[current]]$file.par$zero_offset) * 
                        globalSettings$position.1D, w1_upfield = xy$y[2], ACTIVE = TRUE, 
                      nDim = nDim)
      if (!is.null(ncol(roiTable)))
        newRoi <- c(newRoi, rep(list(NA), ncol(roiTable) - 7))
      names(newRoi) <- names(roiTable)
      roiTable <- orderROI(rbind(roiTable, data.frame(newRoi, 
                                                      stringsAsFactors = FALSE) ))
      
      drawPeptides(); showRoi(rTable = roiTable)
      
    }else{
      
      ## Show the new ROI and append it to the roi table
      rect(xleft = xy$x[2], ybottom = xy$y[2], xright = xy$x[1], 
           ytop = xy$y[1], border = 'red', lwd = 2 )
      newRoi <- list( Name = 'ROI', w2_downfield = xy$x[2], 
                      w2_upfield = xy$x[1],	w1_downfield = xy$y[2],  w1_upfield = xy$y[1], 
                      ACTIVE = TRUE, nDim = nDim)
      if (!is.null(ncol(roiTable)))
        newRoi <-  c(newRoi, rep(list(NA), ncol(roiTable) - 7))
      names(newRoi) <- names(roiTable)
      roiTable <- orderROI(rbind(roiTable, data.frame(newRoi, 
                                                      stringsAsFactors = FALSE) ))
      showRoi(rTable = roiTable)			
    }	
  }
  showGui()
  
  ## Assign the new roi table 	
  if (is.data.frame(roiTable))
    myAssign("roiTable", roiTable)
  
  ## Redraw the spectrum   
  refresh(... )
  
}

## Automatically draw rois
## w1Delta - w1 chemical shift window for each roi
## w2Delta - w2 chemical shift window for each roi
## p - padding percentage to be added to peak widths in 2D spectra when creating 
##     ROIs, if w1Delta and w2Delta are not provided.
## noiseFilt - Integer argument that can be set to 0, 1 or 2; 
##              0 does not apply a noise filter, 1 applies a mild filter
##              (adjacent points in the direct dimension must be above the 
##              noise threshold), 2 applies a strong filter (all adjacent points
##              must be above the noise threshold
## ... - Additional arguments can be passed peakPick()
ra <- function(w1Delta=globalSettings$roi.w1, w2Delta=globalSettings$roi.w2, 
               p=globalSettings$roi.pad, noiseFilt=globalSettings$roi.noiseFilt, ...) {
  
  ## Find current spectrum and establish chemical shift range
  p <- p/100 + .5
  current <- wc()
  usr <- fileFolder[[current]]$graphics.par$usr
  filePar <- fileFolder[[current]]$file.par
  nuc <- filePar$nucleus
  res <- (filePar$downfield_ppm - filePar$upfield_ppm ) / filePar$matrix_size	
  nDim <- filePar$number_dimensions 
  
  if( nDim > 1 ){
    
    ## Tell the user that something is happening
    if( w1Delta == 0 && w2Delta == 0 )
      cat('Measuring peak widths...\n')
    else
      cat('Generating peak list...\n')
    flush.console()
    
    ## Calculate the granularity constants
    gran <- NULL
    for( i in 1:2 ){
      if(nuc[i] == '1H')
        gran <- c(gran, .05)
      else
        gran <- c(gran, .5)
    }
    w1Gran <- gran[1] %/% res[1] + 1
    w2Gran <- gran[2] %/% res[2] + 1
    
    ## Generate a fancy peak list 
    peakList <- peakPick( w1Range= usr[3:4], w2Range = usr[1:2], append = FALSE, 
                          internal = TRUE, fancy = TRUE, w1Gran = w1Gran, w2Gran = w2Gran, 
                          noiseFilt = noiseFilt, ... )	
    if( is.null(peakList) ){
      cat('No peaks were detected \n')
      return(invisible())
    }
    peakList <- peakList[which( is.na(peakList$Multiplet)),]
    
    ## Set the peak widths
    if(w1Delta == 0)
      peakList$w1D <- (peakList$w1D * p) + res[1]
    else
      peakList$w1D <- w1Delta/2
    if(w2Delta == 0)
      peakList$w2D <- (peakList$w2D * p) + res[2]
    else
      peakList$w2D <- w2Delta/2
    
    ## Translate fancy peak list to roi table
    peakList <- data.frame(cbind( peakList$w2 + peakList$w2D, 
                                  peakList$w2 - peakList$w2D, peakList$w1 + peakList$w1D, 
                                  peakList$w1 - peakList$w1D), stringsAsFactors = FALSE)
    names(peakList) <- c('w2_downfield', 'w2_upfield', 'w1_downfield', 
                         'w1_upfield')	
    peakList$Name <- rep('ROI', nrow(peakList))
    peakList$ACTIVE <- TRUE
    peakList$nDim <- 2
    peakList <- peakList[, match(c('Name', 'w2_downfield', 'w2_upfield', 
                                   'w1_downfield', 'w1_upfield', 'ACTIVE', 'nDim'), 
                                 names(peakList))]
  }else{
    
    ## Calculate the granularity constant		
    if( w2Delta == 0 && nuc[1] == '1H')
      w2Delta <- .05
    if( w2Delta == 0 && nuc[1] != '1H')
      w2Delta <- .5
    w2Gran <- w2Delta %/% res[1] + 1
    
    ## Generate a peak list
    peakList <- peakPick(w2Range = usr[1:2], append = FALSE, internal = TRUE, 
                         w2Gran = w2Gran, noiseFilt = noiseFilt, ...)
    if( is.null(peakList) ){
      cat('No peaks were detected \n')
      return(invisible())
    }
    
    ## Translate peak list to roi table
    peakList <- data.frame(cbind( peakList$w2 + w2Delta/2, 
                                  peakList$w2 - w2Delta/2, filePar$zero_offset, 
                                  peakList$Height + peakList$Height * .15), stringsAsFactors = FALSE)
    names(peakList) <- c('w2_downfield', 'w2_upfield', 'w1_downfield', 
                         'w1_upfield')	
    peakList$Name <- rep('ROI', nrow(peakList))
    peakList$ACTIVE <- TRUE
    peakList$nDim <- 1
    peakList <- peakList[, match(c('Name', 'w2_downfield', 'w2_upfield', 
                                   'w1_downfield', 'w1_upfield', 'ACTIVE', 'nDim'), 
                                 names(peakList))]
    
    ## Center the new rois
    while( TRUE ){
      peakList <- rc(inTable = peakList)
      dup <- which(duplicated(peakList[,2:3]))
      if( length(dup) == 0 )
        break()	
      peakList <- peakList[-(dup),]	
    }
  }
  
  ## Save a backup copy and refresh graphics
  if (!is.null(roiTable) && ncol(roiTable) > 7){
    for (i in 1:(ncol(roiTable) - 7))
      peakList <- cbind(peakList, rep(NA, nrow(peakList)))
    names(peakList) <- names(roiTable)
  }
  peakList$clevel <- fileFolder[[current]]$graphics.par$clevel
  roiTable <- orderROI(rbind(roiTable, peakList))
  myAssign("roiTable", roiTable )
  refresh( main.plot = FALSE )
  showRoi()
  
  invisible(roiTable)	
}

## User wrapper function for lower level plot selection functions
## Allows users to select or deselect ROIs in any active plotting window
## preSel - integer indicating the plot from which ROIs will be selected,
##             1 = list, 2 = main, 3 = sub plot, 4 = multiple file plot, 
##             5 = region in main
## ...  - Plotting options can be passed to drawPeptides 
rs <- function( preSel = NULL, parent=NULL, ... ){
  
  ##Checks appropriate objects
  wc()
  if(length(roiTable) < 1 )
    err('No ROIs have been designated (use: roi())' )
  
  ##Build a readable table of open devices
  graphicsList <- c('List of ROIs', 'Main plot window', 
                    'ROI subplot window', 'Multiple file window', 'Region in main plot')
  
  ##Have user select the plotting window for selection    
  if(!any(preSel == 1:5) )
    usrList <- mySelect(graphicsList, multiple=FALSE, title='Select ROIs from:', 
                        parent=parent)
  else
    usrList <- switch(which(1:5 == preSel), 'List of ROIs', 'Main plot window', 
                      'ROI subplot window', 'Multiple file window', 'Region in main plot' )
  
  ## Invoke the correct selection function
  if(usrList == 'List of ROIs')
    selList(...)
  if(usrList == 'Main plot window')
    selMain(...)
  if(usrList == 'ROI subplot window')
    selSub(...)
  if(usrList == 'Multiple file window')
    selMulti(...)
  if(usrList == 'Region in main plot')
    selRegion(...)
  
}

## Internal helper function selList
## Allows users to activate/deactivate ROIs from a list 
## ...  - Plotting options can be passed to drawPeptides 
selList <- function(...){
  
  ## Have user select the ROIs to activate
  usr.sel <- mySelect(c('NONE', roiTable$Name), multiple = TRUE, index = TRUE,  
                      preselect = roiTable$Name[which(roiTable[, 6] == TRUE)], 
                      title = 'Select ROIs to activate')
  if ( length(usr.sel) == 0 || !nzchar(usr.sel) )
    return(invisible())
  
  ## Update ROI table
  roiTable$ACTIVE <- FALSE
  usr.sel <- usr.sel[ usr.sel != 1 ]
  if( length(usr.sel) != 0 )
    roiTable[usr.sel - 1,]$ACTIVE <- TRUE
  
  ## Refresh the plots
  myAssign('roiTable', roiTable)
  refresh(...)
  showRoi()
  invisible(usr.sel)
}

## Internal helper function selMain
## Allows users to activate/deactivate ROIs from the main plot
## ...  - Plotting options can be passed to drawPeptides 
selMain <- function(...){
  
  ## Set the current spectrum
  current <- wc()
  lineCol <- fileFolder[[current]]$graphics.par$fg
  nDim <- fileFolder[[current]]$file.par$number_dimensions
  if( is.null(roiTable) || nrow(roiTable) == 0 )
    err('No ROIs have been designated (use: roi())' )
  
  ##Open main and ROI subplot window if closed
  if(length(which(dev.list() == 2)) != 1)
    drawPeptides()
  cw(dev=2);	showRoi()
  
  
  cat(paste('\nIn the main plot window:\n',  
            ' Left-click to activate/deactivate ROIs\n',  
            ' Right-click to exit\n\n')) 
  flush.console()
  hideGui()
  
  ## Have user select rois from the main plot
  while( TRUE ){
    
    ## Tell the user what to do
    op <- par('font')
    par( font = 2 )
    legend("topleft", c('LEFT CLICK TO SELECT','RIGHT CLICK TO EXIT'), 
           pch=NULL, bty='n', text.col=lineCol)
    par(font = op)
    
    xy <- locator(1)		
    if(is.null(xy))
      break()		
    
    ## Locate the selected ROI/s
    if(nDim == 1)
      index <- which(roiTable$w2_downfield >= xy$x & 
                       roiTable$w2_upfield <= xy$x )
    else
      index <- which(roiTable$w2_downfield >= xy$x & 
                       roiTable$w2_upfield <= xy$x & roiTable$w1_downfield >= xy$y &
                       roiTable$w1_upfield <= xy$y )		
    if(length(index) == 0 )
      next()
    
    ## Update selection
    for( i in 1:length(index) ){
      if( roiTable[index[i],]$ACTIVE )
        roiTable[index[i],]$ACTIVE <- FALSE
      else
        roiTable[index[i],]$ACTIVE <- TRUE
    }
    
    ## Refresh the plot
    showRoi(rTable = roiTable)
  }
  
  ## Save the final copy of the list and update the plot
  showGui()
  myAssign('roiTable', roiTable)
  refresh( ... )
  showRoi()
}	

## Internal helper function selRegion
## Allows users to activate ROIs from a region in the main plot window
## ...  - Plotting options can be passed to drawPeptides 
selRegion <- function(...){
  
  ## Define the current spectrum
  current <- wc()
  nDim <- fileFolder[[current]]$file.par$number_dimensions
  lineCol <- fileFolder[[current]]$graphics.par$fg
  if( is.null(roiTable) || nrow(roiTable) == 0 )
    err('No ROIs have been designated (use: roi())' )
  
  ## Open main and ROI subplot window if closed
  if (!2 %in% dev.list())
    drawPeptides()
  cw(dev=2)	
  showRoi()
  
  ## Give the user some instructions
  hideGui()
  cat(paste('In the main plot window:\n',  
            ' Left-click two points inside the plot to define region\n'))
  flush.console()
  op <- par('font')
  par(font=2)
  legend("topleft", c('LEFT CLICK TO DEFINE REGION', 'RIGHT CLICK TO CANCEL'), 
         pch=NULL, bty='n', text.col=lineCol)
  par(font=op)
  
  ##define the first boundary for the region
  usrCoord1 <- data.frame(locator(1))
  if (length(usrCoord1) == 0){
    refresh(multi.plot=FALSE, sub.plot=FALSE)
    showGui()
    return(invisible())
  }
  abline(v=usrCoord1$x, col=lineCol )
  if (nDim != 1)
    abline(h=usrCoord1$y, col=lineCol )    
  
  ##define other boundary for the region
  usrCoord2 <- data.frame(locator(1))
  if (length(usrCoord2) == 0){
    refresh(multi.plot=FALSE, sub.plot=FALSE)
    showGui()
    return(invisible())
  }
  abline(v=usrCoord2$x, col=lineCol )
  if (nDim != 1)
    abline(h=usrCoord2$y, col=lineCol )   
  usrCoord <- rbind(usrCoord1, usrCoord2)
  usrCoord$x <- sort(usrCoord$x)
  usrCoord$y <- sort(usrCoord$y)
  showGui()
  
  ## Identify the selected ROIs
  includeUp <- which(roiTable$w2_upfield <= usrCoord$x[2] & 
                       roiTable$w2_upfield >= usrCoord$x[1])
  includeDown <- which(roiTable$w2_downfield <= usrCoord$x[2] & 
                         roiTable$w2_downfield >= usrCoord$x[1])
  if (nDim == 1)
    selRois <- unique(c(includeUp, includeDown))
  else{
    w2Matches <- unique(c(includeUp, includeDown))
    includeUp <- which(roiTable$w1_upfield <= usrCoord$y[2] & 
                         roiTable$w1_upfield >= usrCoord$y[1])
    includeDown <- which(roiTable$w1_downfield <= usrCoord$y[2] & 
                           roiTable$w1_downfield >= usrCoord$y[1])
    w1Matches <- unique(c(includeUp, includeDown))
    selRois <- w1Matches[w1Matches %in% w2Matches]
  }
  
  ## Update selection and refresh
  if (length(selRois) == 0){
    refresh(multi.plot=FALSE, sub.plot=FALSE)
    showRoi()
    return(invisible())
  }
  newTable <- roiTable
  newTable[selRois, 'ACTIVE'] <- TRUE
  myAssign('roiTable', newTable)
  refresh(...)
  showRoi()
}	

## Internal helper function selSub
## Allows users to activate/deactivate ROIs from the sub plot
## ...  - Plotting options can be passed to drawPeptides 
selSub <- function(...){
  
  ##Open main and ROI subplot window if closed
  devClose <- FALSE
  if(length(which(dev.list() == 3)) != 1){
    rvs()
    devClose = TRUE
  }
  
  ## Bring main into focus and reset par
  cw(dev=3) 	
  op <- par(mfg = c(1,1,par('mfg')[3:4])) 
  current <- wc()
  nDim <- fileFolder[[current]]$file.par$number_dimensions
  cat(paste('\nIn the ROI subplot window:\n',  
            ' Left-click to activate/deactivate ROIs\n',  
            ' Right-click to exit\n\n')) 
  flush.console()
  hideGui()
  
  ## Have user select rois from the sub plot
  while( TRUE ){
    
    ## Have user select the roi or exit
    xy <- locator(1)
    if(is.null(xy))
      break()		
    xy$x <- trunc( grconvertX(xy$x, from = 'user', to = 'nic') * op$mfg[4] ) + 1
    xy$y <- trunc( op$mfg[3] - 
                     grconvertY( xy$y, from = 'user', to = 'nic') * op$mfg[3] + 1 ) 
    userSel <- xy$x + (op$mfg[4] * (xy$y-1))
    if(userSel > nrow(roiTable) )
      next()
    
    ## Update selection
    if(roiTable[userSel, 6]){
      col = 'white'
      roiTable[userSel, 6] <- FALSE
    }else{
      col = 'red'
      roiTable[userSel, 6] <- TRUE
    }
    
    ## Refresh the plot
    par(mfg = c(unlist(rev(xy)), op$mfg[3:4]))
    box( which = 'plot', col = col)
  }
  
  ## Save the final copy of the list and update the plot
  showGui()
  par( op )
  myAssign('roiTable', roiTable)
  if(devClose)
    dev.off(3)
  refresh( ... )
  showRoi()
}	

## Internal helper function selMulti
## Allows users to activate/deactivate ROIs from the ROI multi plot
## ...  - Plotting options can be passed to drawPeptides 
selMulti <- function(...){
  
  ## Open select list if no ROIs are active
  if(length(which(roiTable$ACTIVE)) < 1 ){
    actRoi <- selList() 
    if(length( actRoi ) < 1 )
      return(invisible())
  }
  
  ## Open select list if no files are active
  if( length(which(sapply(fileFolder, 
                          function(x){x$graphics.par$roi.multi}))) < 1 )
    rsf()
  
  
  ##Open main and ROI subplot window if closed
  devClose <- FALSE
  if(length(which(dev.list() == 4)) != 1){
    rvm()
    devClose <- TRUE
  }
  
  ## Bring main into focus and reset par
  cw( dev = 4 ) 	
  op <- par(mfg = c(1,1,par('mfg')[3:4]))
  if (.Platform$OS.type == 'windows')
    cat(paste('\nIn the multiple file window:\n',  
              ' Left-click on file names to change files\n',  
              ' Left-click on ROI names to change ROIs\n',
              ' Left-click on spectra to view ROI in main window\n',
              ' Right-click to exit\n\n'))
  else
    cat(paste('\nIn the multiple file window:\n',  
              ' Left-click on spectra to view ROI in main window\n',
              ' Right-click to exit\n\n'))
  flush.console()	
  hideGui()
  
  ## Have user select rois from the sub plot
  while( TRUE ){
    
    ## Have user select the roi or exit
    xy <- locator(1)
    if(is.null(xy) )
      break()		
    xy$x <- trunc( grconvertX(xy$x, from = 'user', to = 'nic') * op$mfg[4] ) 
    xy$y <- trunc( op$mfg[3] - 
                     grconvertY( xy$y, from = 'user', to = 'nic') * op$mfg[3]  ) 
    
    ## Skip cases where nothing is selected
    if (xy$x == 0 && xy$y == 0 )
      next()
    
    ## Open ROI list if user selects an roi name
    if (xy$y == 0 && xy$x != 0){
      if (.Platform$OS.type == 'windows'){
        actRoi <- selList()
        if(length( actRoi ) < 1 ){
          cw(2)
          showGui()
          return(invisible())
        }
      }else
        next()
      cw(4)
      op <- par(mfg = c(1,1,par('mfg')[3:4]))
      break()
    }
    
    ## Open active file list if user selects a file name
    if( xy$x ==0 && xy$y !=0 ){
      if (.Platform$OS.type == 'windows'){
        if (is.null(rsf())){
          showGui()
          return(invisible())
        }
      }else
        next()
      cw(4)
      op <- par(mfg = c(1,1,par('mfg')[3:4]))
      break()
    }
    
    
    ## Change current spectrum 
    currentSpectrum <- names(which(sapply(fileFolder, 
                                          function(x){x$graphics.par$roi.multi})))[xy$y]
    myAssign('currentSpectrum', currentSpectrum, save.backup = FALSE)
    current <- wc()
    nDim <- fileFolder[[current]]$file.par$number_dimensions
    if (nDim == 3)
      nDim <- 2
    
    ## Update zoom
    tusr <- as.numeric(roiTable[roiTable$ACTIVE,][xy$x,2:5])
    if ( roiTable[roiTable$ACTIVE,][1,]$nDim != nDim && nDim == 2)
      tusr[3:4] <- c(fileFolder[[current]]$file.par$downfield_ppm[1], 
                     fileFolder[[current]]$file.par$upfield_ppm[1])
    if ( roiTable[roiTable$ACTIVE,][1,]$nDim != nDim && nDim == 1 )
      tusr[3:4] <- c(fileFolder[[current]]$file.par$min_intensity, 
                     fileFolder[[current]]$file.par$max_intensity)
    setGraphics(usr = tusr)
    
    ## Refresh the graphics
    refresh( multi.plot = FALSE, ...)
    cw(4)
  }
  showGui()
  par( op )
  if(devClose)
    dev.off(4)
  dev.set(2)
  showRoi()
}

## Activate all ROIs 
rsAll <- function(){
  if(is.null(roiTable))
    err('There are no ROIs to select')
  roiTable$ACTIVE <- TRUE
  myAssign('roiTable', roiTable)
  refresh()
}

## Deactivate all ROIs
rdAll <- function(){
  if(is.null(roiTable))
    err('There are no ROIs to select' )
  roiTable$ACTIVE <- FALSE
  myAssign('roiTable', roiTable)
  refresh()
}

## Allows user to delete ROIs 
##roiNames - character vector; names of ROIs to delete
rDel <- function(roiNames){
  if(is.null(roiTable))
    stop('There are no ROIs to delete', call.=FALSE)
  
  ## Present users with possible ROIs to delete
  if (missing(roiNames)){
    del.list <- mySelect(roiTable$Name, 
                         roiTable[roiTable$ACTIVE == TRUE, ]$Name, multiple = TRUE, index=TRUE, 
                         title = 'ROIs to be deleted')
    if( length(del.list) == 0 || !nzchar(del.list) )
      return(invisible())  
  }else{
    del.list <- unlist(na.omit(match(roiNames, roiTable$Name)))
    if (!length(del.list))
      return(invisible())
  }
  
  ## Double check that the user wants to delete the ROIs
  userChoice <- myMsg(type = "okcancel",
                      message=paste('Confirm deletion of ', length(del.list), ' ROI(s)?'))
  if(userChoice != 'ok')
    return(invisible())		
  
  ## Remove selected rois and duplicates and renumber rois
  roiTable <- roiTable[- del.list, ]
  roiTable$Name <- as.vector(sapply(roiTable[,1], 
                                    function(x){unlist(strsplit(x, ".", fixed = TRUE))[1]}))
  roiTable <- unique(roiTable)
  for( i in 1: length(roiTable$Name) ){
    tName <- which(roiTable$Name == roiTable$Name[i])
    if( length (tName) != 1 )
      roiTable$Name[tName] <- paste(roiTable$Name[i], 
                                    1:length(tName), sep ='.' )
  }
  
  ## Update file folder and refresh the graphics	
  row.names(roiTable) <- NULL
  if (nrow(roiTable) == 0)
    roiTable <- NULL
  myAssign('roiTable', roiTable)
  refresh()       
} 

## User ROI function rr
## Rename active ROIs
rr <- function(){
  
  ##display in error if there are no ROIs
  if (is.null(roiTable) || !nrow(roiTable))
    err('No ROIs have been designated, use roi()')	
  
  ##get active ROI names
  activeNames <- roiTable[roiTable$ACTIVE,]$Name
  if (all(activeNames == activeNames[1]))
    defName <- activeNames[1]
  else
    defName <- ''
  
  ##ask user for new name
  newName <- myDialog('ROI name:', defName, 'Rename ROIs')
  if (is.null(newName))
    return(invisible())
  
  ##assign name to active ROIs
  newTable <- roiTable
  newTable[roiTable$ACTIVE, 'Name'] <- newName
  newTable <- orderROI(newTable)
  myAssign('roiTable', newTable)
  refresh()
}

## Allow users to select the roi files in memory they wish to view in the
## multiple file window
rsf <- function(){
  
  ##Check for open files
  wc()
  
  ## Find the active roi files
  actFiles <- names(which(sapply(fileFolder, 
                                 function(x){x$graphics.par$roi.multi})))
  allFiles <- getTitles(names(fileFolder), FALSE)
  
  ## Have user select files to display    
  usrList <- mySelect(allFiles, multiple = TRUE, preselect = actFiles,
                      title = 'Select files for the ROI plot', index = TRUE )
  if( length(usrList) == 0 || !nzchar(usrList) )
    return(invisible())
  
  ## Reset the file selection an refresh
  setGraphics(all.files = TRUE, roi.multi = FALSE, save.backup = FALSE)
  setGraphics(file.name = names(fileFolder)[usrList], roi.multi = TRUE, 
              save.backup = TRUE)
  refresh( main.plot = FALSE, overlay = FALSE, sub.plot = FALSE )
  invisible(names(fileFolder)[usrList])
}

## User function rSum
## Summerizes NMR spectral ROIs by intensity, area, or chemical shift
## ask - logical; if TRUE, a series of dialogs are displayed to obtain input 
##	parameters from the user and all other arguments will be ignored
## sumFiles - character string/vector; spectrum name(s) as returned by 
##						names(fileFolder)
## sumRois - character string/vector; ROIs to be included in the summary
## sumType - character string; indicates how the ROI data should be summarized,
##	must be one of "maximum", "minimum", "absMax", "area", "absArea", "w1", or 
##	"w2"
## normType - character string; indicates how the ROI data should be normalized,
##	must be one of "none", "internal", "crossSpec", "signal/noise", or "sum"
## normList - characters string or vector; one or more files or ROIs to be used
##	when normalizing the ROI data (only applicable if normType is set to 
##	"internal" or "crossSpec")
## Returns ROI summary data table
rSum <- function(ask=TRUE, sumFiles, sumRois, sumType, normType='none', 
                 normList=NA){
  
  ## Check to make sure the requisite files/ROIs are present
  wc()
  if (!exists('roiTable') || is.null(roiTable))
    err( 'No ROIs have been designated, use roi()')
  
  ## Get summary parameters from user if ask is set to TRUE
  actFiles <- getTitles(names(which(sapply(fileFolder, function(x) 
    x$graphics.par$roi.multi))))
  
  if (ask){
    
    ## Have user select the files to use for the summary
    sumFiles <- mySelect(getTitles(names(fileFolder), FALSE), 
                         title='Select summary files:', preselect=actFiles, multiple=TRUE, 
                         index=TRUE)
    if (!length(sumFiles) || !nzchar(sumFiles))
      return(invisible())	
    sumFiles <- names(fileFolder)[sumFiles]
    
    ## Have user select the ROIs to use for the summary
    sumRois <- mySelect(roiTable$Name, title='Select summary ROIs:', 
                        preselect=roiTable$Name[roiTable$ACTIVE], multiple=TRUE)
    if (!length(sumRois) || !nzchar(sumRois))
      return(invisible())	
    
    ## Have user select the summary type 
    sumType <- mySelect(c('maximum', 'absolute max', 'minimum', 'area', 
                          'absolute area', 'w1 shift', 'w2 shift', 'custom'), 
                        preselect='Maximum', multiple=FALSE, title='Summarize ROI by:', 
                        index=TRUE)
    if (!length(sumType) || !nzchar(sumType))
      return(invisible())	
    sumType <- c('maximum', 'absMax', 'minimum', 'area', 'absArea', 'w1', 
                 'w2', 'custom')[sumType]
    if (sumType == 'custom'){
      funList <- NULL
      for (i in ls('.GlobalEnv')){
        if (is.function(get(i)))
          funList <- c(funList, i)
      }
      if (is.null(funList))
        err(paste('You must create a function in the global\nenvironment to',
                  'use for the summary.'))
      usrFun <- NULL
      while (is.null(usrFun)){
        usrFun <- mySelect(funList, multiple=FALSE, 
                           title='Select function name:')
        if (!nzchar(usrFun))
          return(invisible())
        if (is.null(formals(get(usrFun)))){
          usrSel <- myMsg(paste('You must select the name of a function from',
                                'the global\nenvironment that takes at least one argument.'), 
                          'okcancel', 'error')
          if (usrSel == 'cancel')
            return(invisible())
          else
            usrFun <- NULL
        }
      }
    }
    
    ## Have user select normilization mode
    if (sumType == 'w1' || sumType == 'w2'){
      normType <- 'none'
    }else{
      normType <- mySelect(c('NONE', 'Internal standard', 'Across spectra', 
                             'Signal to noise', 'Constant sum'), preselect='NONE',	
                           title='Select normalization:', index=TRUE)
      if (!length(normType) || !nzchar(normType))
        return(invisible())
      normType <- c('none', 'internal', 'crossSpec', 'signal/noise', 
                    'sum')[normType]
      if (normType == 'internal'){
        normList <- mySelect(sumRois, multiple=TRUE, 
                             title='Identify Internal Standards:')
        if (!length(normList) || !nzchar(normList))
          return(invisible())
      }else if (normType == 'crossSpec'){
        normList <- mySelect(getTitles(sumFiles, FALSE), title='Select spectra:', 
                             multiple=TRUE, index=TRUE)
        if (!length(normList) || !nzchar(normList))
          return(invisible())
        normList <- sumFiles[normList]
      }else{
        normList <- NA
      }
    }
  }else{
    if (missing(sumType))
      err('The "sumType" argument must be provided')
    if (sumType == 'custom'){
      funList <- NULL
      for (i in ls('.GlobalEnv')){
        if (is.function(get(i)))
          funList <- c(funList, i)
      }
      if (is.null(funList))
        err(paste('You must create a function in the global\nenvironment to',
                  'use for the summary.'))
      usrFun <- NULL
      while (is.null(usrFun)){
        usrFun <- mySelect(funList, multiple=FALSE, 
                           title='Select function name:')
        if (!nzchar(usrFun))
          return(invisible())
        if (is.null(formals(get(usrFun)))){
          usrSel <- myMsg(paste('You must select the name of a function from',
                                'the global\nenvironment that takes at least one argument.'), 
                          'okcancel', 'error')
          if (usrSel == 'cancel')
            return(invisible())
          else
            usrFun <- NULL
        }
      }
    }
    if (sumType == 'w1' || sumType == 'w2'){
      normType <- 'none'
    }else{
      if (normType == 'internal' || normType == 'crossSpec'){
        if (is.na(normList[1]))
          err('The "normList" argument is required')
      }else{
        normList <- NA
      }
    }
    if (missing(sumFiles))
      sumFiles <- actFiles
  }
  
  ## Convert user selections into functions
  sumFun <- switch(match(sumType, c('maximum', 'absMax', 'minimum', 'area', 
                                    'absArea', 'w1', 'w2', 'custom')), 
                   function(x, parm){max(x$data)}, 
                   function(x, parm){x$data[which.max(abs(x$data))]},
                   function(x, parm){min(x$data)},    
                   function(x, parm){x$graphics.par=parm; return(peakVolume(x))}, 
                   function(x, parm){x$data=abs(x$data); x$graphics.par=parm; 
                   return(peakVolume(x))}, 
                   function(x, parm){maxShift(x, conDisp=parm$conDisp)$w1},  
                   function(x, parm){maxShift(x, conDisp=parm$conDisp)$w2}, 
                   get(usrFun))
  if (is.null(sumFun))
    err(paste('Summary type must be either "maximum", "absMax", "minimum",',
              '"area", "absArea", "w1", "w2", or "custom"', sep=''))
  if (length(formals(sumFun)) < 2)
    formals(sumFun) <- c(formals(sumFun), alist(...=))
  
  ## Generate summary for each active ROI
  fileData <- outData <- NULL
  cat('Processing files: \n')
  for (i in sumFiles){
    cat(paste(basename(i), '. . . ' ))
    flush.console()
    currPar <- fileFolder[[i]]$file.par
    nDim <- currPar$number_dimensions
    if (nDim == 3)
      nDim <- 2
    fileData <- NULL
    
    ## Read data from active ROIs
    for (j in  match(sumRois, roiTable$Name)){
      
      w1Range <- as.numeric(roiTable[j, 4:5])
      w2Range <- as.numeric(roiTable[j, 2:3])
      
      ## Fudge the ROI table for 1D/2D compatibility
      if (roiTable$nDim[j] == 1 && nDim == 2)
        w1Range <- c(currPar$downfield_ppm[1], currPar$upfield_ppm[1])
      
      ## Read roi from binary and reject ROIs outside of the spectral window
      roiData <- ucsf2D(i, w1Range=w1Range, w2Range=w2Range, file.par=currPar)
      if ((nDim == 2 && (length(roiData$w1) < 2 || length(roiData$w2) < 2)) ||
          (nDim == 1 && length(roiData$data) < 2 ))
        fileData <- c(fileData, as.numeric(NA))
      else
        fileData <- c(fileData, sumFun(roiData, fileFolder[[i]]$graphics.par))      
    }
    
    ## Normalize data
    names(fileData) <- sumRois
    if (normType == 'internal'){
      
      ## Normalize by ROIs
      fileData <- fileData / mean(fileData[normList])
    }else if (normType == 'signal/noise'){
      
      ## Normalize by noise level
      fileData <- fileData / fileFolder[[i]]$file.par$noise_est
    }else if (normType == 'sum'){
      
      ## Normalize by sum of spectral data
      fileData <- fileData / sum(ucsf2D(i, 
                                        file.par=fileFolder[[i]]$file.par)$data)
    }
    outData <- rbind(outData, fileData)
    cat('done \n' )
    flush.console()
  }
  if (normType == 'crossSpec'){
    
    ## Normalize by spectra
    rownames(outData) <- sumFiles
    if (length(normList) > 1)
      normData <- mean(data.frame(outData[normList, ]))
    else
      normData <- unlist(outData[normList, ])
    if (ncol(outData) > 1)
      outData <- t(apply(outData, 1, function(x) x / normData))
    else{
      outData <- apply(outData, 1, function(x) x / normData)
    }
    normList <- getTitles(names(fileFolder), FALSE)[match(normList, 
                                                          names(fileFolder))]
  }
  
  ## Format data
  outData <- data.frame(outData, row.names=NULL)
  colnames(outData) <- sumRois
  fileNames <- getTitles(names(fileFolder), FALSE)[match(sumFiles, 
                                                         names(fileFolder))]
  outData <- data.frame(fileNames, outData, stringsAsFactors=FALSE)
  names(outData)[1] <- 'File'
  
  ## Set up outgoing file structure
  sumPar <- list(sumType, normType, normList)
  names(sumPar) <- c('summary.type', 'normalization', 'norm.data.source')
  newSum <- list(outData, sumPar)
  names(newSum) <- c('data', 'summary.par')
  
  ## Update summary table and print ROI data
  myAssign('roiSummary', newSum)
  se()
  
  return(invisible(newSum))
}

## User ROI function rc 
## Centers slected rois about the max peak in each active ROI
## massCenter  - logical argument; TRUE centers peaks by center of mass,
##               false centers peaks by maximum signal observed
## inTable		-  Used by internal functions only
## note: graphics settings are used for choosing between negative and positive 
##       signals.
rc <- function ( massCenter = TRUE, inTable ){
  
  if( !is.logical (massCenter) ){
    err('massCenter must argument must be either TRUE or FALSE')
  }
  
  if( !missing(inTable) )
    roiTable <- inTable
  if( is.null(roiTable) || nrow(roiTable) == 0)
    err('No ROIs have been designated, use roi()')	
  if( length(roiTable$ACTIVE) < 1 )
    err('No ROIs have been selected')
  
  ## Define the current spectrum
  current <- wc()
  nDim <- fileFolder[[current]]$file.par$number_dimensions
  conDisp <- fileFolder[[current]]$graphics.par$conDisp
  
  ##Open main and ROI subplot window if closed
  if(length(which(dev.list() == 2)) != 1)
    drawPeptides()
  
  ## Bring main into focus and show the rois 
  cw(dev=2);	showRoi()
  
  
  ## Apply function to all active ROIs  
  for(i in which(roiTable$ACTIVE)){
    
    ## Find current chemical shift range
    w1Range <- sort(as.numeric(roiTable[i, 4:5]))
    w2Range <- sort(as.numeric(roiTable[i, 2:3])) 
    
    ## Find the max shift for the current ROI
    current.roi <- maxShift(ucsf2D(currentSpectrum, w1Range = w1Range, 
                                   w2Range = w2Range, file.par = fileFolder[[current]]$file.par), 
                            conDisp = conDisp, massCenter = massCenter )
    if(is.null(current.roi))
      next()
    
    ## Update ROI table
    roiTable[i, 2:3 ] <- c(current.roi$w2 + diff(w2Range)/2, 
                           current.roi$w2 - diff(w2Range)/2)
    if( nDim > 1 && roiTable$nDim[i] > 1 ){
      roiTable[i, 4:5 ] <- c(current.roi$w1 + diff(w1Range)/2, 
                             current.roi$w1 - diff(w1Range)/2)			
    }    
  }
  
  ## Assign roi table to global environment and refresh the plot
  if(missing(inTable)){
    myAssign( 'roiTable', roiTable )
    refresh() 		
  }
  invisible(roiTable)   
}

## move ROI upfield in direct dimension (right)
##p - percentage to move ROI by
rmr <- function (p=1){
  changeRoi(w2Inc = c(-p, -p))
}

## move ROI downfield in direct dimension (left)
##p - percentage to move ROI by
rml <- function (p=1){
  changeRoi(w2Inc = c(p, p))
}

## move ROI upfield in indirect dimension (right)
##p - percentage to move ROI by
rmu <- function (p=1){
  changeRoi(w1Inc = c(-p, -p))
}

## move ROI downfield in indirect dimension (right)
##p - percentage to move ROI by
rmd <- function (p=1){
  changeRoi(w1Inc = c(p, p))
}

## Expand roi in direct dimension
##p - percentage to expand ROI by
red <- function (p=1){
  changeRoi(w2Inc = c(p, -p))
}

## Contract roi in direct dimension
##p - percentage to contract ROI by
rcd <- function (p=1){
  changeRoi(w2Inc = c(-p, p))
}

## Expand roi in the indirect dimension
##p - percentage to expand ROI by
rei <- function (p=1){
  changeRoi(w1Inc = c(p, -p))
}

## Contract roi in the indirect dimension
##p - percentage to contract ROI by
rci <- function (p=1){
  changeRoi(w1Inc = c(-p, p))
}

## User function for toggling the ROI display within the main plot window
rv <- function(){
  if (globalSettings$roiMain){
    setGraphics(roiMain=FALSE)
    refresh(multi.plot=FALSE, sub.plot=FALSE)
    cat('ROI display off \n')
  }else{
    setGraphics(roiMain=TRUE)
    showRoi()
    cat('ROI display on \n')
  }
}

############################################################
#                                                          #
#    User functions to Load/save objects and workspaces    #
#                                                          #
############################################################

## Import data from file
## object - character string; name of digestR object to import
import <- function(object, parent=NULL){
  
  ## Have user select type of file to import
  tclCheck()
  current <- wc()
  if (missing(object)){
    usrSel <- mySelect(c('ROI table', 'ROI summary', 'Peak list'), 
                       title='Import file type:', parent=parent)
    if ( length(usrSel) == 0 || !nzchar(usrSel) )
      return(invisible())
  }else
    usrSel <- object
  title <- switch(usrSel, 'ROI table'='Import ROI table', 
                  'ROI summary'='Import ROI summary', 'Peak list'='Import peak list')
  
  ## Have user select a file
  fileName <- myOpen(initialfile='', title=title, defaultextension='txt',
                     filetypes=list('xls'='Excel Files', 'txt'='Text Files'))
  if (length(fileName) == 0)
    return(invisible())
  overAppend <- 'Append'
  
  ## Import an ROI table
  if (usrSel == 'ROI table'){
    
    ## Checks for the correct columns
    newRoi <- read.table(fileName, header=FALSE, nrows=1, sep='\t', 
                         stringsAsFactors=FALSE)
    if (sum(newRoi ==  c('Name', 'w2_downfield', 'w2_upfield', 'w1_downfield', 
                         'w1_upfield', 'ACTIVE','nDim')) < 7)
      err(paste('ROI tables must be tab delimited text files with the',
                'following columns:\n\n', 
                '                           1. Name', '\n', 
                '                           2. w2_downfield', '\n', 
                '                           3. w2_upfield', '\n', 
                '                           4. w1_downfield', '\n', 
                '                           5. w1_upfield', '\n', 
                '                           6. ACTIVE', '\n', 
                '                           7. nDim'))
    
    ## Read ROI table and check with user for overwrite
    newRoi <- read.table(fileName, header=TRUE, sep='\t',	
                         stringsAsFactors=FALSE)
    if (exists('roiTable') && length(roiTable) > 0)
      overAppend <- buttonDlg('An ROI table already exists.', c('Overwrite', 
                                                                'Append', 'Cancel'), default='Cancel', parent=parent)
    if (overAppend == 'Cancel')
      return(invisible())
    
    ## Checks that the ROI table has the correct format
    newRoi <- re(newRoi)
    if (is.null(newRoi))
      return(invisible())
    
    ## Append new ROIs to existing table
    if (overAppend == 'Append'){
      
      ##check for matching columns in the two tables
      newNames <- colnames(newRoi)
      oldNames <- colnames(roiTable)
      if (!is.null(newNames) && !identical(newNames, oldNames)){
        
        ##add any columns from the new table missing from the existing table
        for (i in newNames[-match(oldNames, newNames, nomatch=0)]){
          roiTable <- cbind(roiTable, rep(NA, nrow(roiTable)))
          colnames(roiTable)[ncol(roiTable)] <- i
        }
        
        ##add any columns from the existing table missing from the new table
        for (i in oldNames[-match(newNames, oldNames, nomatch=0)]){
          newRoi <- cbind(newRoi, rep(NA, nrow(newRoi)))
          colnames(newRoi)[ncol(newRoi)] <- i
        }
      }
    }
    
    ## Assigns the ROI table to the global environment
    myAssign('roiTable', newRoi)		
    refresh()
    
    ## Open roi subplot window if it is closed
    if(length(which(dev.list() == 2)) == 1)
      showRoi()       
  }
  
  ## Import an ROI summary
  if (usrSel == 'ROI summary'){
    
    ## Checks for the correct columns
    newSum <- read.table(fileName, header=FALSE, nrows=1, sep='\t', 
                         stringsAsFactors = FALSE)
    fileIndex <- which(newSum == 'FILE')
    groupIndex <- which(newSum == 'GROUP')
    if (length(fileIndex) != 1 || length(groupIndex) > 1 ||
        any(!is.na(suppressWarnings(as.numeric(newSum[2:length(newSum)])))))
      err(paste('ROI summaries must be tab delimited text files with:\n',
                '1. A single "GROUP" column', '\n', 
                '2. A single "FILE" column and', '\n',
                '3. Columns for each ROI', '\n',
                ' (with headings beginning with a character)'))
    
    ## Checks that the columns are in the correct order
    newSum <- read.table(fileName, header=TRUE, sep='\t',	
                         stringsAsFactors=FALSE)
    if (length(groupIndex) == 0){
      tmp <- rep('Group 1', length(newSum$FILE))
      newSum <- data.frame(tmp, newSum, stringsAsFactors=FALSE)
      names(newSum)[1] <- 'GROUP'
      groupIndex <- which(names(newSum) == 'GROUP')
      fileIndex <- which(names(newSum) == 'FILE')		
    }
    if (groupIndex != 1){
      tmp <- newSum[, groupIndex]
      newSum <- newSum[, -groupIndex]
      newSum <- data.frame(tmp, newSum, stringsAsFactors=FALSE)
      names(newSum)[1] <- 'GROUP'
      groupIndex <- which(names(newSum) == 'GROUP')
      fileIndex <- which(names(newSum) == 'FILE')
    }
    if (fileIndex != 2){
      tmp <- newSum[, fileIndex]
      newSum <- newSum[, -fileIndex]
      newSum <- data.frame(newSum$GROUP, tmp, newSum[, -1], 
                           stringsAsFactors=FALSE)
      names(newSum)[1:2] <- c('GROUP', 'FILE')
    }
    
    ## Checks that the ROI summary has the correct format
    newSum <- re(newSum)
    if (is.null(newSum))
      return(invisible())
    
    ## Assigns the ROI summary to the global environment
    roiSummary <- NULL
    roiSummary$data <- newSum 
    roiSummary$summary.par$summary.type <- NA
    roiSummary$summary.par$normalization <- NA
    roiSummary$summary.par$norm.data.source <- NA
    myAssign('roiSummary', roiSummary)		     
  }
  
  ## Import a peak list
  if (usrSel == 'Peak list'){
    
    ## Checks for the correct columns
    newPeak <- read.table(fileName, header=FALSE, nrows=1, sep='\t', 
                          stringsAsFactors=FALSE)
    w2Index <- which(newPeak == 'w2')
    w1Index <- which(newPeak == 'w1')
    heightIndex <- which(newPeak == 'Height')
    assignIndex <- which(newPeak == 'Assignment')
    inIndex <- which(newPeak == 'Index')
    if (length(inIndex) > 1 || length(w1Index) > 1 || length(w2Index) > 1 || 
        length(heightIndex) > 1 || length(assignIndex) > 1)
      err(paste('Peak lists must contain only one of each of the following',
                'columns:\n\n', 
                '                          1. w2', '\n', 
                '                          2. w1', '\n', 
                '                          3. Height', '\n', 
                '                          4. Assignment', '\n', 
                '                          5. Index'))
    
    if (fileFolder[[current]]$file.par$number_dimensions > 1){
      if (length(w1Index) == 0 || length(w2Index) == 0)
        err('Peak lists for 2D spectra must contain a w1 and w2 column.')
    }else{
      if (length(w2Index) == 0)
        err('Peak lists for 1D spectra must contain a w2 column.')
    }
    
    ## Read peak list and check with user for overwrite
    newPeak <- read.table(fileName, header=TRUE, sep='\t', 
                          stringsAsFactors=FALSE)
    if (length(fileFolder[[current]]$peak.list) > 0)
      overAppend <- buttonDlg('A peak list already exists.', c('Overwrite', 
                                                               'Append', 'Cancel'), default='Cancel', parent=parent)
    if (overAppend == 'Cancel')
      return(invisible())
    
    ## Adds missing columns
    if (length(w1Index) == 0){
      w1 <- rep(NA, length(newPeak$w2))
      newPeak <- data.frame(newPeak, w1, stringsAsFactors=FALSE)
    }
    if (length(heightIndex) == 0){
      Height <- rep(NA, length(newPeak$w2))
      newPeak <- data.frame(newPeak, Height, stringsAsFactors=FALSE)
    }
    if (length(assignIndex) == 0){
      Assignment <- rep(NA, length(newPeak$w2))
      newPeak <- data.frame(newPeak, Assignment, stringsAsFactors=FALSE)
    }
    if (length(inIndex) == 0){
      Index <- 1:length(newPeak$w2)
      newPeak <- data.frame(newPeak, Index, stringsAsFactors=FALSE)
    }
    
    ## open peak list editor
    newPeak <- pe(newPeak)
    if (is.null(newPeak))
      return(invisible())
    
    ## Assigns the peak list to the global environment		
    newPeak$w1 <- as.numeric(newPeak$w1)
    newPeak$w2 <- as.numeric(newPeak$w2)
    if (overAppend == 'Append')
      newPeak <- appendPeak(newPeak, fileFolder[[current]]$peak.list)
    fileFolder[[current]]$peak.list <- newPeak
    setGraphics(peak.disp=TRUE, save.backup=FALSE)
    pdisp()
    myAssign('fileFolder', fileFolder)	
    refresh(sub.plot=FALSE, multi.plot=FALSE)
  }
}

## Export a data table to tab delimited file
## object - character string; name of digestR object to export
export <- function(object, parent=NULL){
  
  ## Have user select type of file to export
  current <- wc()
  usrList <- NULL
  if (!is.null(roiTable))
    usrList <- c(usrList, 'ROI table')
  if (exists('roiSummary') && !is.null(roiSummary$data))
    usrList <- c(usrList, 'ROI summary')
  if (!is.null(fileFolder[[current]]$peak.list))
    usrList <- c(usrList, 'Peak list')
  if (is.null(usrList))
    err('There are currently no data tables to export.')
  if (missing(object)){
    usrSel <- mySelect(usrList, title='Export file type:', parent=parent)
    if (length(usrSel) == 0 || !nzchar(usrSel))
      return(invisible())
  }else
    usrSel <- object
  
  ## Have user select file name
  initFile <- switch(usrSel, 'ROI table'='roiTable', 
                     'ROI summary'='roiSummary', 'Peak list'='peakList')
  fileName <- mySave(defaultextension='txt', initialfile=initFile, 
                     title='Export', filetypes=list('xls'='Excel Files', 'txt'='Text Files'))
  if (length(fileName) == 0 || !nzchar(fileName))
    return(invisible())
  
  ##writes the data table to the given file name
  dataTable <- switch(usrSel, 'ROI table'=roiTable, 
                      'ROI summary'=roiSummary$data, 
                      'Peak list'=fileFolder[[current]]$peak.list)
  write.table(dataTable, file=fileName, quote=FALSE, sep='\t', row.names=FALSE, 
              col.names=TRUE)
  log_message(paste('The data were saved to: ', fileName), quote=FALSE)
}

## Extract data from digestR objects
## Returns selected object and prints a summary to the console
ed <- function(){
  
  ## Have user select object to extract data from
  current <- wc()
  usrList <- c('Main plot window', 'Slice')
  if (!is.null(roiTable))
    usrList <- c(usrList, 'ROI table')
  if (exists('roiSummary') && !is.null(roiSummary$data))
    usrList <- c(usrList, 'ROI summary')
  if (!is.null(fileFolder[[current]]$peak.list))
    usrList <- c(usrList, 'Peak list')
  usrSel <- mySelect(usrList, title='Extract from:')
  if ( length(usrSel) == 0 || !nzchar(usrSel) )
    return(invisible())
  
  ## Extract data from the main plot window
  if (usrSel == 'Main plot window'){
    nDim <- fileFolder[[current]]$file.par$number_dimensions
    w1Range <- fileFolder[[current]]$graphics.par$usr[3:4]
    w2Range <- fileFolder[[current]]$graphics.par$usr[1:2]	
    return(ucsf2D(file.name = currentSpectrum, w1Range = w1Range, 
                  w2Range = w2Range, file.par = fileFolder[[current]]$file.par))
  }
  
  ## Extract data from a 1D slice
  if (usrSel == 'Slice'){
    tmp <- vs()
    return(tmp)
  }
  
  ## Extract data from an ROI table
  if (usrSel == 'ROI table')
    return(roiTable)
  
  ## Extract data from an ROI summary
  if (usrSel == 'ROI summary'){
    return(roiSummary[1:2])
  }
  
  ## Extract data from a peak list
  if (usrSel == 'Peak list'){
    return(fileFolder[[current]]$peak.list)
  }
}

## Load an R workspace
## fileName - character string, the file path for the workspace to load
## plot - logical, replots the currentSpectrum if TRUE
## clearAll - logical, clears all previous objects if TRUE, otherwise only digestR
##	          objects are cleared
load <- wl <- function(fileName, plot=TRUE, clearAll=TRUE){
  
  if (missing(fileName))
    fileName <- myOpen(filetypes=list('RData'='R image'), 
                       defaultextension='RData', multiple=FALSE, title='Load Workspace')
  if (length(fileName) == 0 || !nzchar(fileName))
    stop('Load cancelled')
  backup <- file.path(Sys.getenv('HOME'), 'zal3waozq')
  save.image(backup)
  tryCatch({if (clearAll)
    suppressWarnings(rm(list=ls(envir=.GlobalEnv), envir=.GlobalEnv))
    else{
      digestRob <- c('fileFolder', 'currentSpectrum', 'oldFolder', 'roiTable', 
                  'roiSummary', 'pkgVar', 'globalSettings', 'overlayList')
      suppressWarnings(rm(digestRob, envir=.GlobalEnv))
    }
    suppressWarnings(base::load(file=fileName, envir=.GlobalEnv))
    digestR:::patch()
    # gui()
    if (exists('fileFolder') && !is.null(fileFolder) && plot)
      dd()
    cat('"', fileName, '"', ' successfully loaded','\n', sep='')}, 
    error=function(er){
      if (length(grep('magic', er$message))){
        base::load(file=backup, envir=.GlobalEnv)
        suppressWarnings(file.remove(backup))
        err('Invalid workspace, no data loaded.')
      }
      suppressWarnings(file.remove(backup))
      return(er)
    })
  invisible(suppressWarnings(file.remove(backup)))
}

## Save an R workspace
## fileName - character string, the file path for the workspace to save to
ws <- function(fileName){
  
  if (missing(fileName))
    fileName <- mySave(defaultextension='RData', title='Save Workspace',
                       filetypes=list('RData'='R image'))
  if (length(fileName) == 0 || !nzchar(fileName))
    return(invisible())
  base::save.image(file=fileName)
  
  cat(paste('Workspace saved to ','"', fileName, '"', '\n', sep=''))
}

## Restore an R workspace
rb <- function(){
  backupFile <- file.path(path.expand('~'), '.digestRbackup')
  if (!file.exists(backupFile))
    return('No backup to restore, load cancelled.')
  tryCatch(wl(backupFile), error=function(er) 
    log_message('Unable to restore workspace.'))
}

############################################################
#                                                          #
#      General digestR fileFolder object utility functions    #
#                                                          #
############################################################

## Internal function for checking existence of files within fileFolder
## halt - logical, stops execution after updating files if TRUE
updateFiles <- function(halt=TRUE){
  
  if (is.null(fileFolder))
    return(invisible())
  
  ## Create list of files that need to be updated
  fileNames <- names(fileFolder)
  updateList <- file.access(unique(fileNames))
  updateList <- names(updateList[which(updateList == -1)])
  if (!length(updateList))
    return(invisible())
  
  ## Have the user select files to update
  usr <- myMsg(paste('Could not find previously opened file(s), press OK to', 
                     'update file locations'), type='okcancel', icon='error')
  
  ## Do nothing if user selects cancel
  if (usr == 'cancel')
    stop('File location update cancelled.', call.=FALSE)
  
  ## Allow user to select files to update
  prevPaths <- mySelect(updateList, multiple=TRUE, 
                        title='Select files to update')
  
  ## Do nothing if user selects cancel
  if (!length(prevPaths) || !nzchar(prevPaths))
    stop('File location update cancelled.', call.=FALSE)
  
  ## Close spectra the user chose not to update
  newFolder <- fileFolder
  newCS <- NULL
  newOL <- overlayList
  if (length(prevPaths) != length(updateList)){
    closeList <- updateList[-as.vector(na.omit(match(prevPaths, updateList)))]
    message <- paste('The following files will be closed:', paste(closeList, 
                                                                  collapse=' \n  '), sep='\n  ')
    usr <- myMsg(message, icon='info', type='okcancel')
    
    ## Do nothing if user selects cancel
    if (usr == 'cancel')
      stop('File location update cancelled.', call.=FALSE)
    
    newFolder <- newFolder[-match(closeList, fileNames)]
  }
  
  ## Get new file location
  newPath <- myOpen(title=paste('Update location for "', basename(prevPaths[1]),
                                '"', sep=''), multiple=FALSE)
  
  ## Do nothing if user selects cancel
  if (!length(newPath) || !nzchar(newPath))
    stop('File location update cancelled.', call.=FALSE)
  
  ## Updates selected file
  i <- match(prevPaths[1], fileNames)
  names(newFolder)[i] <- newPath
  prevUpfield <- newFolder[[i]]$file.par$upfield_ppm
  prevDownfield <- newFolder[[i]]$file.par$downfield_ppm
  newFolder[[i]]$file.par <- ucsfHead(newPath, FALSE)$file.par
  newFolder[[i]]$file.par$upfield_ppm <- prevUpfield
  newFolder[[i]]$file.par$downfield_ppm <- prevDownfield
  newDir <- dirname(newPath)
  if (currentSpectrum == prevPaths[1])
    newCS <- newPath
  if (!is.null(newOL)){
    overMatch <- match(prevPaths[1], newOL)
    if (!is.na(overMatch))
      newOL[overMatch] <- newPath
  }
  prevPaths <- prevPaths[-1]
  
  ## Updates remaining files
  if (length(prevPaths)){
    remaining <- TRUE
    while(remaining){
      
      ## Looks for file matches in provided directory
      newDirFiles <- list.files(newDir, full.names=TRUE)
      fileMatches <- na.omit(match(basename(prevPaths), basename(newDirFiles)))
      if (length(fileMatches)){
        newPaths <- newDirFiles[fileMatches]
        prevPathMatches <- match(basename(newPaths), basename(prevPaths))
        dirPrevPaths <- prevPaths[prevPathMatches]
        prevPaths <- prevPaths[-prevPathMatches]
        
        ## Updates files that match
        for (newPath in newPaths){
          i <- match(dirPrevPaths[1], fileNames)
          names(newFolder)[i] <- newPath
          prevUpfield <- newFolder[[i]]$file.par$upfield_ppm
          prevDownfield <- newFolder[[i]]$file.par$downfield_ppm
          newFolder[[i]]$file.par <- ucsfHead(newPath, FALSE)$file.par
          newFolder[[i]]$file.par$upfield_ppm <- prevUpfield
          newFolder[[i]]$file.par$downfield_ppm <- prevDownfield
          if (currentSpectrum == newPath)
            newCS <- newPath
          if (!is.null(newOL)){
            overMatch <- match(newPath, newOL)
            if (!is.na(overMatch))
              newOL[overMatch] <- newPath
          }
          dirPrevPaths <- dirPrevPaths[-1]
        }
        
        ## Tells the user which files were updated using previosly provided dir.
        cat(paste('The following files were found in "', newDir, '",\n', 
                  'and were automatically updated:\n  ', paste(basename(newPaths), 
                                                               collapse=' \n  '), '\n',	sep=''))
        flush.console()
        
        ## Checks for files that still need to be updated
        if (!length(prevPaths))
          remaining <- FALSE
      }else{
        
        ## Allows the user to select a new directory
        newPath <- myOpen(title=paste('Update location for "', prevPaths[1],
                                      '"', sep=''), multiple=FALSE)
        
        ## Do nothing if user selects cancel
        if (!length(newPath) || !nzchar(newPath))
          stop('File location update cancelled.', call.=FALSE)
        
        ## Updates selected file
        i <- match(prevPaths[1], fileNames)
        names(newFolder)[i] <- newPath
        prevUpfield <- newFolder[[i]]$file.par$upfield_ppm
        prevDownfield <- newFolder[[i]]$file.par$downfield_ppm
        newFolder[[i]]$file.par <- ucsfHead(newPath, FALSE)$file.par
        newFolder[[i]]$file.par$upfield_ppm <- prevUpfield
        newFolder[[i]]$file.par$downfield_ppm <- prevDownfield
        newDir <- dirname(newPath)
        if (currentSpectrum == newPath)
          newCS <- newPath
        if (!is.null(newOL)){
          overMatch <- match(newPath, newOL)
          if (!is.na(overMatch))
            newOL[overMatch] <- newPath
        }
        
        ## Checks for files that still need to be updated
        prevPaths <- prevPaths[-1]
        if (!length(prevPaths))
          remaining <- FALSE
      }
    }
  }	
  
  ## Assigns new objects to global environment
  if (length(newFolder) && is.null(newCS))
    newCS <- names(newFolder)[length(newFolder)]
  myAssign('currentSpectrum', newCS, save.backup=FALSE)
  myAssign('overlayList', newOL, save.backup=FALSE)
  myAssign('fileFolder', newFolder)
  refresh()
  if (halt)
    stop('Files updated, previous function may need to be recalled.', 
         call.=FALSE)
}

## Internal function for getting the current spectrum
## fileName - logical; returns the file name for the current spectrum if TRUE
## Returns the current spectrum's index in the file folder
wc <- function(fileName=FALSE){
  
  ##Checks for open files
  if (!exists("fileFolder") || is.null(fileFolder))
    err('The file folder is empty, use fo()')
  
  ##Checks for currentSpectrum
  if (!exists("currentSpectrum") || is.null(currentSpectrum))
    err('The file folder is empty, use fo()')
  
  ##Find the current spectrum in fileFolder
  current <- match(currentSpectrum, names(fileFolder))
  
  ##Return the full file path for the current spectrum
  if (fileName)
    return(fileFolder[[current]]$file.par$file.name)
  
  ##Return the fileFolder index for the current spectrum
  return(current)
}


############################################################
#                                                          #
#       Internal replacement functions and dialogs         #
#                                                          #
############################################################

## Internal function getTitles
## Given a vector of file names (as in names(fileFolder), returns user titles
## inNames - numeric or character vector, the file names to retrieve titles for
## index - logical argument, if TRUE, appends indices to the items in the list
## returns a vector of length inNames containing user titles
getTitles <- function( inNames, index = TRUE ){
  
  ## Get the user_title field for each file name
  inNames <- as.vector(sapply(fileFolder[inNames], function(x) 
    x$file.par$user_title))
  
  ## Append index number to list
  if (index)
    inNames <- paste(seq_along(inNames), inNames, sep=') ')
  
  return(inNames)
}

##Internal function 'mySelect'
##platform independant version of select.list()
##uses a modified version of tk_select.list
##index - logical; if TRUE returns the index for the selected item rather than
##        the list item itself
##parent - specifies a tktoplevel to be the parent window for the dialog 
mySelect <- function(list, preselect=NULL, multiple=FALSE, title=NULL, 
                     index=FALSE, parent=NULL){
  
  ##Append indices to list
  inList <- paste(seq_along(list), list, sep=') ')
  
  ##Default to R's list selection dialog on Windows platforms
  if (.Platform$OS.type == 'windows'){
    if( !is.null(preselect) ){
      preselect <- match(preselect, list)
      preselect <- inList[ preselect ]
    }
    usrList <- select.list( inList, preselect = preselect, multiple=multiple, 
                            title=title)
    usrList <- match(usrList, inList)
    if( length(usrList) == 0 || is.na(usrList) )
      return("")
    if( index )
      return( usrList )
    return( list[usrList] )
  }
  
  ##Checks for valid arguments
  if (is.null(preselect))
    preselect <- 1
  else{
    if (!is.character(preselect))
      stop('Invalid preselect argument')
    preselect <- as.numeric(na.exclude(match(preselect, list)))
    if (length(preselect) == 0)
      preselect <- 1	
  }
  if (!is.logical(multiple))
    stop('Invalid multiple argument')
  if (!is.null(title) && !is.character(title))
    stop('Invalid title argument')
  
  ##creates main window
  tclCheck()
  if (is.null(parent))
    dlg <- tktoplevel()
  else
    dlg <- myToplevel('dlg', parent=parent)
  tkwm.title(dlg, 'Selection')
  tcl('wm', 'attributes', dlg, topmost=TRUE)
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  
  ##determine size of list box
  vscr <- hscr <- FALSE
  ht <- length(inList)
  if (ht > 35){
    ht <- 35
    vscr <- TRUE
  }
  wd <- max(nchar(inList))
  if (wd > 40){
    wd <- 40
    hscr <- TRUE
  }
  wd <- wd + 5
  if (wd < 15)
    wd <- 15
  
  ##create font for listFrame
  fonts <- tcl('font', 'name')
  if (!'listFont' %in% as.character(fonts)){
    listFont <- as.character(tcl('font', 'configure', 'TkDefaultFont'))
    tcl('font', 'create', 'listFont', listFont[1], listFont[2], listFont[3], 
        listFont[4], listFont[5], 'bold', listFont[7], listFont[8], listFont[9], 
        listFont[10], listFont[11], listFont[12])
  }
  
  ##create list box
  if (is.null(title))
    title <- 'Select:'
  label <- ttklabel(dlg, text=title)
  listFrame <- ttklabelframe(dlg, padding=3, labelwidget=label)
  tcl(label, 'configure', '-font', 'listFont')
  lvar <- tclVar()
  tclObj(lvar) <- inList
  listBox <- tklistbox(listFrame, height=ht, width=wd, listvariable=lvar, 
                       selectmode=ifelse(multiple,'extended', 'browse'), active='dotbox', 
                       exportselection=FALSE,  bg='white',
                       xscrollcommand=function(...) tkset(xscr, ...), 
                       yscrollcommand=function(...) tkset(yscr, ...))
  xscr <- ttkscrollbar(listFrame, orient='horizontal', 
                       command=function(...) tkxview(listBox, ...))
  yscr <- ttkscrollbar(listFrame, command=function(...) tkyview(listBox, ...))
  if (length(inList) > 2){
    for (i in seq(0, length(inList) - 1, 2))
      tkitemconfigure(listBox, i, background='#ececff')
  }
  for (i in preselect)
    tkselection.set(listBox, i - 1)
  tcl(listBox, 'see', i - 1)
  
  ##create ok button
  bottomFrame <- ttkframe(dlg)
  returnVal <- ''
  onOK <- function() {
    usrSel <- as.integer(tkcurselection(listBox))
    if (length(usrSel) != 0){
      if (index)
        returnVal <<- 1 + usrSel
      else
        returnVal <<- list[1 + usrSel]
    }
    tkgrab.release(dlg)
    tkdestroy(dlg)
  }
  okButton <- ttkbutton(bottomFrame, text='OK', width=10, command=onOK)
  
  ##create cancel button
  onCancel <- function() {
    tkgrab.release(dlg)
    tkdestroy(dlg)
  }
  cancelButton <- ttkbutton(bottomFrame, text='Cancel', width=10, 
                            command=onCancel)
  tkbind(dlg, '<Destroy>', onCancel)
  
  ##add widgets to listFrame
  tkgrid(listFrame, column=1, row=1, sticky='nswe', pady=10, padx=c(14, 0))
  tkgrid(listBox, column=1, row=1, sticky='nswe')
  if (vscr)
    tkgrid(yscr, column=2, row=1, sticky='ns')
  if (hscr)
    tkgrid(xscr, column=1, row=2, sticky='we')
  
  ##make listFrame stretch when window is resized
  tkgrid.columnconfigure(dlg, 1, weight=1)
  tkgrid.rowconfigure(dlg, 1, weight=10)
  tkgrid.columnconfigure(listFrame, 1, weight=1)
  tkgrid.rowconfigure(listFrame, 1, weight=1)
  
  ##add buttons to bottom of toplevel
  tkgrid(bottomFrame, column=1, row=2, padx=c(22, 0))
  tkgrid(okButton, column=1, row=1, padx=4)
  tkgrid(cancelButton, column=2, row=1, padx=4)
  tkgrid(ttksizegrip(dlg), column=3, row=3, sticky='se')
  
  ##Allows users to press the 'Enter' key to make selections
  onEnter <- function(){
    focus <- as.character(tkfocus())
    if (length(grep('.2.1$', focus)))
      onOK()
    else
      tryCatch(tkinvoke(focus), error=function(er){})
  }
  tkbind(dlg, '<Return>', onEnter) 
  tkbind(listBox, '<Double-Button-1>', onOK)
  tkactivate(listBox, max(preselect) - 1)
  
  ##configure dialog window
  tkwm.deiconify(dlg)
  if (as.logical(tkwinfo('viewable', dlg)))
    tkgrab.set(dlg)
  tkfocus(listBox)
  tkwait.window(dlg)
  
  return(returnVal)
}

## Internal function 'myDialog'
## Tk version of winDialogString, creates a dialog with ok\cancel buttons and a 
##   text entry widget
## message - character string; message to display in the dialog
## default - character string; default text in the entry widget
## title - character string; title for the window
## entryWidth - positive integer; horizontal length of the entry widget
## parent - specifies a tktoplevel to be the parent window for the dialog 
myDialog <- function(message='', default='', title='digestR', entryWidth=20, 
                     parent=NULL){
  
  ##creates main window
  tclCheck()
  if (is.null(parent))
    dlg <- tktoplevel()
  else
    dlg <- myToplevel('dlg', parent=parent)
  tkwm.title(dlg, title)
  tkwm.resizable(dlg, FALSE, FALSE)
  tcl('wm', 'attributes', dlg, topmost=TRUE)
  if (.Platform$OS == 'windows')
    tcl('wm', 'attributes', dlg, toolwindow=TRUE)
  tkfocus(dlg)
  
  ##creates message label
  msgLabel <- ttklabel(dlg, text=message)
  
  ##creates text entry widget
  usrEntry <- tclVar(default)
  textEntry <- ttkentry(dlg, width=entryWidth, justify='center', 
                        textvariable=usrEntry)
  tkselection.range(textEntry, 0, nchar(default))
  
  ##creates ok button
  returnVal <- NULL
  onOK <- function(){
    returnVal <<- tclvalue(usrEntry)
    tkgrab.release(dlg)
    tkdestroy(dlg)
  }
  okButton <- ttkbutton(dlg, text="OK", width=8, command=onOK)
  
  ##creates cancel button
  onCancel <- function(){
    tkgrab.release(dlg)
    tkdestroy(dlg)
  }
  cancelButton <- ttkbutton(dlg, text="Cancel", width=8, command=onCancel)
  
  ##add widgets to toplevel
  tkgrid(msgLabel, column=1, columnspan=2, row=1, pady=c(8, 5), padx=6, 
         sticky='w')
  tkgrid(textEntry, column=1, columnspan=2, row=2, pady=5, padx=20)
  tkgrid(okButton, column=1, row=3, pady=8, padx=c(6, 1))
  tkgrid(cancelButton, column=2, row=3, pady=8, padx=c(1, 6))
  
  ##selects the text in the entry widget when Ctrl+A is pressed
  onCtrlA <- function(){
    tkfocus(textEntry)
    tkselection.range(textEntry, 0, nchar(tclvalue(usrEntry)))
  }
  tkbind(dlg, '<Control-a>', onCtrlA)
  
  ##allows users to press the 'Enter' key to make selections
  onEnter <- function(){
    focus <- as.character(tkfocus())
    if (length(grep('.2$', focus)))
      onOK()
    else
      tryCatch(tkinvoke(focus), error=function(er){})
  }
  tkbind(dlg, '<Return>', onEnter)
  
  ##configure dialog window
  tkwm.deiconify(dlg)
  if (as.logical(tkwinfo('viewable', dlg)))
    tkgrab.set(dlg)
  tkfocus(textEntry)
  tkwait.window(dlg)
  
  if (!is.null(returnVal))
    return(returnVal)
  invisible(returnVal)
}

## Internal function 'myDir'
## Modified version of tkchooseDirectory that saves and uses the last directory
## initialdir - specifies the initial directory displayed in the dialog
## parent - specifies a tktoplevel to be the parent window for the dialog 
## title - specifies a string to display as the title of the dialog
## mustexist - specifies whether or not the user must select an existing dir.
## see tcl/tk manual for additional documentation
myDir <- function(initialdir='', parent=NULL, title='', mustexist=TRUE){
  
  ## Get saved directory if no initial directory is entered
  if( initialdir == '' && !is.null(pkgVar$prevDir) )
    initialdir <- pkgVar$prevDir
  
  ## Let the user choose a directory using tk dialog
  tclCheck()
  if (!is.null(parent))
    returnVal <- tclvalue(tkchooseDirectory(initialdir=initialdir, title=title,
                                            parent=parent, mustexist=mustexist))
  else
    returnVal <- tclvalue(tkchooseDirectory(initialdir=initialdir, title=title, 
                                            mustexist=mustexist))
  
  ## Save the selected (not canceled) directory
  if (length(returnVal) && nzchar(returnVal)){
    returnVal <- gsub('\\', '/', returnVal, fixed=TRUE)
    pkgVar$prevDir <- returnVal
    myAssign("pkgVar", pkgVar, save.backup = FALSE)
  }
  
  return(returnVal)
}

## Internal function 'myOpen'
## Modified version of tkgetOpenFile
## defaultextension - a string specifying the file extension (without ".") for 
##										the default filetype to be displayed in	the dialog (must 
##										match one of the extensions provided in filetypes)
## filetypes - adds the provided file types to the file types listbox in the
##             open dialog if supported by the platform, must be in list format,
##             as such: list(txt = "Text File", xls = "Excel File")
## initialfile - specifies a filename to be displayed initially in the dialog
## initialdir - specifies the initial directory displayed in the dialog
## multiple - logical, should users be able to select more than one file?
## title - specifies a string to display as the title of the dialog
## parent - specifies a tktoplevel to be the parent window for the dialog 
## see tcl/tk manual for additional documentation
myOpen <- function( defaultextension='', filetypes='', initialfile='', 
                    initialdir='', multiple=TRUE, title='', parent=NULL){
  
  ## Get saved directory if no initial directory is entered
  if( initialdir == '' && !is.null(pkgVar$prevDir) )
    initialdir <- pkgVar$prevDir
  
  ## Reformat filetypes parameters to work with tkgetOpenFile
  if (is.list(filetypes)){
    tkParse <- NULL
    for (i in seq_along(filetypes)){
      if (names(filetypes)[i] == defaultextension)
        next
      tkParse <- paste(c(tkParse, paste('{{', filetypes[[i]],'} ', '{.',
                                        names(filetypes)[i], '}}', sep='')), collapse=' ')
    }
    defIn <- match(defaultextension, names(filetypes))
    #		if (!is.na(defIn)){
    sysInfo <- Sys.info()
    if (sysInfo['release'] == 'Vista' || sysInfo['release'] == '7'){
      tkParse <- paste(tkParse, paste('{{', filetypes[[defIn]],'} ', '{.',
                                      names(filetypes)[i], '}}', sep=''))
      tkParse <- paste('{{All Files} *}', tkParse)
    }else{
      tkParse <- paste(paste('{{', filetypes[[defIn]],'} ', '{.', 
                             names(filetypes)[i], '}}', sep=''), tkParse)
      tkParse <- paste(tkParse, '{{All Files} *}')
    }
    #		}
  }else
    tkParse <- ''
  
  ## Get the path to the file/files the user enters
  tclCheck()
  if (multiple){
    if (!is.null(parent))
      returnVal <- as.character(tkgetOpenFile(filetypes=tkParse, title=title,
                                              initialdir=initialdir, initialfile=initialfile, parent=parent, 
                                              multiple=multiple))
    else
      returnVal <- as.character(tkgetOpenFile(filetypes=tkParse, title=title,
                                              initialdir=initialdir, initialfile=initialfile, 
                                              multiple=multiple))
  }else{
    if (!is.null(parent))
      returnVal <- tclvalue(tkgetOpenFile(filetypes=tkParse, title=title,
                                          initialdir=initialdir, initialfile=initialfile, parent=parent, 
                                          multiple=multiple))
    else
      returnVal <- tclvalue(tkgetOpenFile(filetypes=tkParse, title=title,
                                          initialdir=initialdir, initialfile=initialfile, 
                                          multiple=multiple))
  }
  
  ## Save the selected (not canceled) directory
  if (length(returnVal) > 0 && nzchar(returnVal)){
    pkgVar$prevDir <- dirname(returnVal[1])
    myAssign("pkgVar", pkgVar, save.backup=FALSE)
  }
  return( returnVal )
}

## Internal function 'mySave'
## Modified version of tkgetSaveFile
## defaultextension - a string specifying the file extension to be appended to
##                    the file name if one is not provided
## filetypes - adds the provided file types to the file types listbox in the
##             open dialog if supported by the platform, must be in list format,
##             as such: list(txt = "Text File", xls = "Excel File")
## initialfile - specifies a filename to be displayed initially in the dialog
## initialdir - specifies the initial directory displayed in the dialog
## title - specifies a string to display as the title of the dialog
## parent - specifies a tktoplevel to be the parent window for the dialog 
## see tcl/tk manual for additional documentation
mySave <- function(defaultextension='', filetypes='', initialfile='', 
                   initialdir='', title='', parent=NULL){
  
  ## Get saved directory if no initial directory is entered
  if (initialdir == '' && !is.null(pkgVar$prevDir))
    initialdir <- pkgVar$prevDir
  
  ## Make sure a '.' is included in defaultextension argument (Linux issue)
  if (nzchar(defaultextension) && 
      !length(grep('.', defaultextension, fixed=TRUE)))
    defaultextension <- paste('.', defaultextension, sep='')
  
  ## Reformat filetypes parameters to work with tkgetSaveFile
  if (is.list(filetypes)){
    tkParse <- NULL
    for (i in seq_along(filetypes)){
      if (names(filetypes)[i] == defaultextension)
        next
      tkParse <- paste(c(tkParse, paste('{{', filetypes[[i]],'} ', '{.',
                                        names(filetypes)[i], '}}', sep='')), collapse=' ')
    }
    defIn <- match(defaultextension, names(filetypes))
    if (!is.na(defIn)){
      sysInfo <- Sys.info()
      if (sysInfo['release'] == 'Vista' || sysInfo['release'] == '7' || sysInfo['release'] == '10 x64')
        tkParse <- paste(tkParse, paste('{{', filetypes[[defIn]],'} ', '{.',
                                        names(filetypes)[i], '}}', sep=''))
      else
        tkParse <- paste(paste('{{', filetypes[[defIn]],'} ', '{.', 
                               names(filetypes)[i], '}}', sep=''), tkParse)
    }
    tkParse <- paste('{{All Files} *}', tkParse)
  }else
    tkParse <- ''
  
  ## Get the path to the file/files the user enters
  tclCheck()
  if (!is.null(parent))
    returnVal <- tclvalue(tkgetSaveFile(filetypes=tkParse, title=title,
                                        defaultextension=defaultextension, initialdir=initialdir,
                                        initialfile=initialfile, parent=parent))
  else
    returnVal <- tclvalue(tkgetSaveFile(filetypes=tkParse, title=title,
                                        defaultextension=defaultextension, initialdir=initialdir,
                                        initialfile=initialfile))
  
  ## Save the selected (not canceled) directory
  if (length(returnVal) > 0 && nzchar(returnVal)){
    pkgVar$prevDir <- dirname(returnVal[1])
    myAssign("pkgVar", pkgVar, save.backup=FALSE)
  }
  return(returnVal)
}

## Internal function 'myMessageBox'
## message - character string, the message to display in the dialog
## type - character string, the type of buttons to display in the dialog
## icon - character string, the type of icon to display in the dialog
## title - character string, the title for the dialog box
## parent - specifies a tktoplevel to be the parent window for the dialog 
## tk version of winDialog, tcl/tk manual for additional documentation
myMsg <- function(message='', type='ok', icon='question', title='digestR', 
                  parent=NULL){
  
  tclCheck()
  if (!is.null(parent))
    return(tclvalue(tkmessageBox(message=message, type=type, icon=icon, 
                                 title='digestR', parent=parent)))
  else
    return(tclvalue(tkmessageBox(message=message, type=type, icon=icon, 
                                 title='digestR')))
}

## Internal utility function for returning errors, invokes stop and 
## opens the error message in a tk window
## message - the desired error message
## parent - specifies a tktoplevel to be the parent window for the dialog 
## halt - logical, stops code execution if TRUE
err <- function(message, parent=NULL, halt=TRUE){
  
  ##display error message
  myMsg(message, icon='error', parent=parent)
  
  ##return focus to console
  bringFocus(-1)
  
  ##halt code execution
  if (halt)
    stop(message, call.=FALSE)
}

## Internal version of file
## Opens the updateFile function if a file path has changed
## fileName - character string; full path name for the spectrum to open a file
##						connection for
## open - character string; see R documentation for "file" function
myFile <- function(fileName, open){
  suppressWarnings(tryCatch(file(fileName, open), error=function(er){
    updateFiles()}))
}

## Internal function for creating a Tk dialog with up to three buttons
## message - character string; message to display in the dialog
## buttons - character vector; names for the buttons to display
## default - character string; specifies the default button and return value
## checkBox - logical; TRUE indicates that the last button should be a checkbox
## title - character string, the title for the dialog box
## parent - specifies a tktoplevel to be the parent window for the dialog 
buttonDlg <- function(message, buttons, checkBox=FALSE, default=buttons[1], 
                      title='digestR', parent=NULL){	
  
  ##creates dialog window
  if (is.null(parent))
    dlg <- tktoplevel()
  else
    dlg <- myToplevel('dlg', parent=parent)
  tkwm.title(dlg, title)
  tkwm.resizable(dlg, FALSE, FALSE)
  tcl('wm', 'attributes', dlg, topmost=TRUE)
  if (checkBox)
    checkVal <- tclVar(0)
  returnVal <- default
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  
  ##displays message
  msgLabel <- ttklabel(dlg, text=message, justify='left')
  tkgrid(msgLabel, column=1, row=1, pady=c(15, 0), padx=15)
  
  ##determines width of buttons
  if (checkBox)
    butWidth <- max(nchar(buttons[1:(length(buttons) - 1)]))
  else
    butWidth <- max(nchar(buttons))
  if (butWidth < 8)
    butWidth <- 8
  
  ##creates buttons
  buttonFrame <- ttkframe(dlg)
  tkgrid(buttonFrame, column=1, row=2, pady=10)
  buttonList <- as.list(buttons)
  for (i in seq_along(buttons)){
    if (i != length(buttons) || !checkBox){
      onButton <- function(){
        usrSel <- tclvalue(tkcget(tkfocus(), '-text'))
        if (checkBox){
          returnVal <<- data.frame(usrSel, 
                                   as.logical(as.integer(tclvalue(checkVal))), 
                                   stringsAsFactors=FALSE)
          names(returnVal) <<- c('button', 'checked')
        }else
          returnVal <<- usrSel
        tkgrab.release(dlg)
        tkdestroy(dlg)
      }
      buttonList[[i]] <- ttkbutton(buttonFrame, text=buttons[i], width=butWidth, 
                                   command=onButton)
    }else
      buttonList[[i]] <- ttkcheckbutton(buttonFrame, text=buttons[i], 
                                        variable=checkVal)
    if (i == 1)
      tkgrid(buttonList[[i]], column=i, row=1, padx=c(12, 3))
    else if (i == length(buttons)){
      if (checkBox)
        tkgrid(buttonList[[i]], column=i, row=1, padx=c(10, 12))
      else
        tkgrid(buttonList[[i]], column=i, row=1, padx=c(3, 12))
    }else
      tkgrid(buttonList[[i]], column=i, row=1, padx=3)
  }
  
  ##configure dialog window
  defButton <- buttonList[[match(default, buttons)]]
  tkconfigure(defButton, state='active')
  tkbind(dlg, '<Destroy>', function(...) return(returnVal))
  tkbind(dlg, '<Return>', function(...) tryCatch(tkinvoke(tkfocus()), 
                                                 error=function(er){}))
  tkwm.deiconify(dlg)
  if (as.logical(tkwinfo('viewable', dlg)))
    tkgrab.set(dlg)
  tkfocus(defButton)
  tkwait.window(dlg)
  
  return(returnVal)
}


############################################################
#                                                          #
#       User Functions for changing R GUI settings         #
#                                                          #
############################################################

## Changes R console settings to SDI mode
sdi <- function(){
  
  ## Doesn't run if Rgui isn't running
  if (.Platform$OS.type != 'windows' || .Platform$GUI != 'Rgui')
    return(invisible())
  
  ## Reads in the Rconsole file from the R user directory
  conPath <- file.path(Sys.getenv('R_USER'), 'Rconsole')
  if (file.exists(file.path(Sys.getenv('R_USER'), 'Rconsole'))){
    readCon <- file(conPath)
    conText <- readLines(readCon, warn=FALSE)
    close(readCon)
    file.remove(conPath)
    
    ## Reads in the Rconsole file from the R home directory
  }else if (file.access(file.path(R.home('etc'), 'Rconsole'), 2) == 0){
    conPath <- file.path(R.home('etc'), 'Rconsole')
    readCon <- file(conPath)
    conText <- readLines(readCon, warn=FALSE)
    close(readCon)
    file.remove(conPath)
    
    ## Copies Rconsole file from R home directory to R user directory
  }else{
    conPath <- file.path(R.home('etc'), 'Rconsole')
    readCon <- file(conPath)
    conText <- readLines(readCon, warn=FALSE)
    close(readCon)
    file.copy(conPath, file.path(Sys.getenv('R_USER'), 'Rconsole'))
  }
  
  ## Writes out a new Rconsole file
  outFile <- conText
  file.create(conPath)
  writeCon <- file(conPath, 'w')
  matches <- NULL
  for (i in c('MDI = yes', 'MDI= yes', 'MDI =yes', 'MDI=yes'))
    matches <- c(matches, length(grep(i, outFile)) != 0)
  if (any(matches)){
    for (i in c('MDI = yes', 'MDI= yes', 'MDI =yes', 'MDI=yes'))
      outFile <- gsub(i, 'MDI = no', outFile)
    writeLines(outFile, writeCon)
  }else
    writeLines(outFile, writeCon)
  close(writeCon)	
  invisible(myMsg(paste('             R console settings updated.', '\n', 
                        'R must be restarted for changes to take effect.'), icon='info'))
}

## Checks that the R console is in SDI mode
## dispMsg - if FALSE, returns TRUE if Rgui is in SDI mode, returns FALSE if 
##	in MDI mode, and does not display a message
sdiCheck <- function(dispMsg=TRUE){
  
  ## Doesn't check Rconsole if Rgui isn't running
  if (.Platform$OS.type != 'windows' || .Platform$GUI != 'Rgui')
    return(TRUE)
  
  ## Doesn't check if sdi is set to FALSE in defaultSettings
  if (!defaultSettings$sdi)
    return(FALSE)
  
  ## Check R user directory for Rconsole file
  if (file.exists(file.path(Sys.getenv('R_USER'), 'Rconsole')))
    conPath <- file.path(Sys.getenv('R_USER'), 'Rconsole')
  
  ## Uses Rconsole file in R home directory
  else
    conPath <- file.path(R.home('etc'), 'Rconsole')
  
  ## Reads in the Rconsole file
  readCon <- file(conPath)
  conText <- readLines(readCon, warn=FALSE)
  close(readCon)
  
  ## Checks for lines in Rconsole file that set MDI to yes
  mdiYes <- 0
  for (i in c('MDI = yes', 'MDI= yes', 'MDI =yes', 'MDI=yes')){
    matchedLines <- grep(i, conText)
    for (j in matchedLines){
      matchedText <- unlist(strsplit(conText[j], '#'))
      if (length(grep(i, matchedText[1])) != 0)
        mdiYes <- j
    }
  }
  
  ## Checks for lines in Rconsole file that set MDI to no
  mdiNo <- 0
  for (i in c('MDI = no', 'MDI= no', 'MDI =no', 'MDI=no')){
    matchedLines <- grep(i, conText)
    for (j in matchedLines){
      matchedText <- unlist(strsplit(conText[j], '#'))
      if (length(grep(i, matchedText[1])) != 0)
        mdiNo <- j
    }
  }
  if (mdiNo <= mdiYes)
    mdiMode <- TRUE
  else
    mdiMode <- FALSE
  if (!dispMsg)
    return(!mdiMode)
  
  ## Warns user about running R in MDI mode
  if (mdiMode){
    usr <- buttonDlg(paste('R is currently running in MDI mode.  For digestR, we',
                           ' suggest\nconfiguring R to display windows separately (SDI mode).',
                           '\n\nWould you like to switch to SDI mode?', sep=''), 
                     buttons=c('Yes', 'No', 'Don\'t display this message again'), TRUE, 
                     default='No')
    if (usr[1] == 'Yes')
      sdi()
    
    ## Edit defaultSettings if user doesn't want message to be displayed again
    if (as.logical(usr[2])){
      defaultSettings$sdi <- FALSE
      writeDef(defSet=defaultSettings)
      myAssign('defaultSettings', defaultSettings)
    }
  }
}

############################################################
#                                                          #
#                       Tk GUIs                            #
#                                                          #
############################################################

## Internal function used on file lists in Tk GUIs
## resets file and overlay lists to match any changes made to fileFolder
reset <- function(lists, boxes, prevPaths, update='files', dims='both'){
  
  ## Restructure inputs if only one list (the files list) is being reset
  if (length(lists) == 1){
    lists <- list(lists, NULL)
    boxes <- list(boxes, NULL)
    prevPaths <- list(prevPaths, NULL)
    update <- c('files', NULL)
  }
  
  ## Assign names to inputs
  names(lists) <- names(boxes) <- names(prevPaths) <- update
  
  ## Reset file lists
  if ('files' %in% update){
    
    ## Update file list using the names in fileFolder 
    if (dims == '1D')
      newPaths <- names(fileFolder)[which(sapply(fileFolder, 
                                                 function(x){x$file.par$number_dimensions}) == 1)]
    else if (dims == '2D')
      newPaths <- names(fileFolder)[which(sapply(fileFolder, 
                                                 function(x){x$file.par$number_dimensions}) > 1)]
    else
      newPaths <- names(fileFolder)
    if (!is.null(newPaths)){
      if ('overlays' %in% update){
        overlayMatches <- match(overlayList, newPaths)
        if (length(overlayMatches))
          newPaths <- newPaths[-overlayMatches]
      }
      tclObj(lists$files) <- getTitles(newPaths)
      
      ## Get previous selection and reset
      prevSel <- prevPaths$files[as.integer(tkcurselection(boxes$files)) + 1]
      if (length(prevSel)){
        tkselection.clear(boxes$files, 0, 'end')
        curSel <- match(prevSel, newPaths, nomatch=0) - 1
        for (i in curSel)
          tkselection.set(boxes$files, i)
      }
      
      ## Alternate colors in listbox
      if (length(newPaths) > 2){
        for (i in seq(0, length(newPaths) - 1, 2))
          tkitemconfigure(boxes$files, i, background='#ececff')
      }
    }else
      tclObj(lists$files) <- ''
  }
  
  ## Reset overlay lists
  if ('overlays' %in% update){
    
    ## Update overlay list using overlayList 
    if (!is.null(overlayList)){
      tclObj(lists$overlays) <- getTitles(overlayList)
      
      ## Get previous selection and reset
      prevSel <- prevPaths$overlays[as.integer(tkcurselection(boxes$overlays)) + 
                                      1]
      if (length(prevSel)){
        tkselection.clear(boxes$overlays, 0, 'end')
        curSel <- match(prevSel, overlayList, nomatch=0) - 1
        for (i in curSel)
          tkselection.set(boxes$overlays, i)
      }
      
      ## Alternate colors in listbox
      if (length(overlayList) > 2){
        for (i in seq(0, length(overlayList) - 1, 2))
          tkitemconfigure(boxes$overlays, i, background='#ececff')
      }
    }else
      tclObj(lists$overlays) <- character(0)
  }
}	

## Internal function for changing colors within Tk GUIs
## parent - the tktoplevel to be used as the parent window for the color widget
## type - the type of color change (eg. 'peak')
## usrFiles - list of files to apply color changes to
changeColor <- function(parent, type, usrFiles=NULL){
  
  ## Display color selection dialogue
  initCol <- switch(type, 'peak'=defaultSettings$peak.color,
                    'bg'=defaultSettings$bg,
                    'axes'=defaultSettings$col.axis,
                    'pos'=defaultSettings$pos.color,
                    'neg'=defaultSettings$neg.color,
                    'proj'=defaultSettings$proj.color,
                    '1D'=defaultSettings$proj.color,
                    'abox'=defaultSettings$roi.bcolor[1],
                    'ibox'=defaultSettings$roi.bcolor[2],
                    'atext'=defaultSettings$roi.tcolor[1],
                    'itext'=defaultSettings$roi.tcolor[2])
  usrColor <- tclvalue(tcl("tk_chooseColor", parent=parent, 
                           initialcolor=initCol))
  
  ## Set color
  if (nzchar(usrColor)){
    switch(type, 
           'peak'=setGraphics(usrFiles, peak.color=usrColor),
           'bg'=setGraphics(usrFiles, bg=usrColor),
           'axes'=setGraphics(usrFiles, line.color=usrColor),
           'pos'=setGraphics(usrFiles, pos.color=usrColor),
           'neg'=setGraphics(usrFiles, neg.color=usrColor),
           'proj'=setGraphics(usrFiles, proj.color=usrColor),
           '1D'=setGraphics(usrFiles, proj.color=usrColor),
           'abox'=setGraphics(roi.bcolor=c(usrColor, 
                                           globalSettings$roi.bcolor[2])),
           'ibox'=setGraphics(roi.bcolor=c(globalSettings$roi.bcolor[1], 
                                           usrColor)),
           'atext'=setGraphics(roi.tcolor=c(usrColor, 
                                            globalSettings$roi.tcolor[2])),
           'itext'=setGraphics(roi.tcolor=c(globalSettings$roi.tcolor[1], 
                                            usrColor)))
    refresh()
    tkfocus(parent)
    tkwm.deiconify(parent)
    bringFocus()
  }
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
  tkadd(plotBook, onedFrame, text='1D Spectra')
  
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
  projColButton <- ttkbutton(coOptionFrame, width=11, text='1D', command=onProj)
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
  onedDefault <- function(){
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
  defaultButton <- ttkbutton(onedOptionFrame, text='Defaults', width=11, 
                             command=onedDefault)
  
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
  tkgrid(onedOptionFrame, column=2, row=1, sticky='nswe', pady=c(10, 2), 
         padx=c(4, 0))
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
    #		else if (length(grep('.1.3.1.1$', focus)))
    #			twodDouble()
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

#' Wrapper function, displays the plot colors pane in ps()
#' @param dispPane A character specifying the default displayed pane in the GUI. Options are "co" (default) for
#'                 the "Plot Colors" pane, and "sp" for the "Spectra" pane.
#'
#' @return None (invisible return).
#'
#' @importFrom tcltk tkwm.title tkfocus tkwm.deiconify tkadd tkgrid ttksizegrip tkselect tkcurselection tkitemconfigure tkselection.set tkbind tkdestroy
#' @export
plot_colors <- function(){
  ps('co')
}

#' Wrapper function, displays the plot settings panes in ps()
#' Plot Settings GUI
#'
#' This function creates a graphical user interface (GUI) window for adjusting plot settings for DIANA spectra.
#' The GUI allows users to change plot colors, configure color options for various components, and switch
#' between different settings panels.
#'
#' @param dispPane A character specifying the default displayed pane in the GUI. Options are "co" (default) for
#'                 the "Plot Colors" pane, and "sp" for the "Spectra" pane.
#'
#' @return None (invisible return).
#'
#' @importFrom tcltk tkwm.title tkfocus tkwm.deiconify tkadd tkgrid ttksizegrip tkselect tkcurselection tkitemconfigure tkselection.set tkbind tkdestroy
#' @export
plot_settings <- function(){
  current <- wc()

  if (fileFolder[[current]]$file.par$number_dimensions == 1)
    ps('ct1D')
  else
    ps('ct2D')
}

## Interactive GUI for manipulating overlays and shift referencing
os <- function(dispPane='ol'){
  
  ##create main window
  current <- wc()
  tclCheck()

  if (as.logical(tcl('winfo', 'exists', '.os'))) { 
  tkdestroy('.os')
  } # New

  # Create a new window
  dlg <- tktoplevel() # New
  tkwm.title(dlg, if (dispPane == 'ol') 'Overlays' else 'Shift Referencing')  # Set window title
  tkwm.geometry(dlg, "800x500+200+200") # New

  # Withdraw the window to prevent flickering while setting up elements
  tkwm.withdraw(dlg) # New
  tkfocus(dlg) # New

  #dlg <- myToplevel('os')
  # if (is.null(dlg))
  # {
  #   if (dispPane == 'ol')
  #   {
  #     tkwm.title('.os', 'Overlays')
  #     tkselect('.os.1', 0)
  #   }else
  #   {
  #     tkwm.title('.os', 'Shift Referencing')
  #     if (dispPane == 'sr1D')
  #       tkselect('.os.1', 1)
  #     else
  #       tkselect('.os.1', 2)
  #   }
  #   return(invisible())
  # }
      if (dispPane == 'ol') {
      tkwm.title('.os', 'Overlays')
      tkselect('.os.1', 0)
    } else {
      tkwm.title('.os', 'Shift Referencing')
      if (dispPane == 'sr1D') {
        tkselect('.os.1', 1)
      } else {
        tkselect('.os.1', 2)
    }
    return(invisible())
  }
    
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  if (dispPane == 'ol')
    tkwm.title(dlg, 'Overlays')
  else
    tkwm.title(dlg, 'Shift Referencing')
  
  ##create paned notebook
  osBook <- ttknotebook(dlg, padding=3)
  
  ##create overlay and referencing panes
  olFrame <- ttkframe(osBook, padding=c(0, 0, 2, 12)) 
  onedFrame <- ttkframe(osBook, padding=c(0, 0, 12, 12))
  twodFrame <- ttkframe(osBook, padding=c(0, 0, 12, 12))
  tkadd(osBook, olFrame, text='   Overlays   ')
  tkadd(osBook, onedFrame, text='1D Referencing')
  tkadd(osBook, twodFrame, text='2D Referencing')
  
  ##add widgets to toplevel
  tkgrid(osBook, column=1, row=1, sticky='nsew', padx=c(6, 0), pady=c(6, 0))
  tkgrid(ttksizegrip(dlg), column=2, row=2, sticky='se')
  tkgrid.columnconfigure(dlg, 1, weight=1)
  tkgrid.rowconfigure(dlg, 1, weight=1)
  
  ##switch to the appropriate notebook pane
  if (dispPane == 'ol')
    tkselect(osBook, 0)
  else if (dispPane == 'sr1D')
    tkselect(osBook, 1)
  else
    tkselect(osBook, 2)
  
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
                         exportselection=FALSE, 
                         xscrollcommand=function(...) tkset(olXscr, ...), 
                         yscrollcommand=function(...) tkset(olYscr, ...),
			 font = "Helvetica 12 bold",
                         bg = "white", fg = "black")

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
  
  ##export fileBox selections to other tabs
  olSelect <- function()
  {
    usrSel <- 1 + as.integer(tkcurselection(olFileBox))
    onedFiles <- names(fileFolder)[which(sapply(fileFolder, function(x)
    {x$file.par$number_dimensions}) == 1)]
    twodFiles <- names(fileFolder)[which(sapply(fileFolder, function(x)
    {x$file.par$number_dimensions}) > 1)]
    if (!is.null(usrSel)){
      onedIndices <- na.omit(match(olFileNames[usrSel], onedFiles))
      if (length(onedIndices)){
        tkselection.clear(onedFileBox, 0, 'end')
        for (i in onedIndices)
          tkselection.set(onedFileBox, i - 1)
      }
      twodIndices <- na.omit(match(olFileNames[usrSel], twodFiles))
      if (length(twodIndices)){
        tkselection.clear(twodFileBox, 0, 'end')
        for (i in twodIndices)
          tkselection.set(twodFileBox, i - 1)
      }
    }
    olConfigGui()
  }
  tkbind(olFileBox, '<<ListboxSelect>>', olSelect)
  
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
  offsetFrame <- ttklabelframe(middleFrame, text='1D Offset')
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
  onPos <- function(){
    usrSel <- 1 +	as.integer(tkcurselection(overlayBox))
    usrFiles <- overlayList[usrSel]
    if (!length(usrFiles))
      err(paste('You must select a spectrum from the Overlays list before',
                'changing the positive contour color'))
    changeColor(dlg, 'pos', usrFiles)
  }
  posColButton <- ttkbutton(colFrame, text='+  Contour', width=13, command=onPos)
  
  ##create set negative contour color button
  onNeg <- function(){
    usrSel <- 1 +	as.integer(tkcurselection(overlayBox))
    usrFiles <- overlayList[usrSel]
    if (!length(usrFiles))
      err(paste('You must select a spectrum from the Overlays list before',
                'changing the negative contour color'))
    changeColor(dlg, 'neg', usrFiles)
  }
  negColButton <- ttkbutton(colFrame, text='-  Contour', width=13, 
                            command=onNeg)
  
  ##create set 1D color button
  onProj <- function(){
    usrSel <- 1 +	as.integer(tkcurselection(overlayBox))
    usrFiles <- overlayList[usrSel]
    if (!length(usrFiles))
      err(paste('You must select a spectrum from the Overlays list before',
                'changing the 1D color'))
    changeColor(dlg, '1D', usrFiles)
  }
  projColButton <- ttkbutton(colFrame, text='1D', width=13, command=onProj)
  
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
  tkgrid(posColButton, column=2, row=1, pady=2)
  tkgrid(negColButton, column=2, row=2, pady=2)
  tkgrid(projColButton, column=2, row=3, pady=c(2, 4))
  tkgrid.columnconfigure(colFrame, 1, weight=1)
  tkgrid.columnconfigure(colFrame, 3, weight=1)
  
  ##add widgets to textFrame
  tkgrid(textFrame, column=2, columnspan=2, row=3, sticky='we')
  tkgrid(textButton, sticky='w', padx=6, pady=6)
  tkgrid(textSuppressButton, sticky='w', padx=6, pady=6)
  
  ##reconfigures widgets in GUI according to which spectra are open
  olConfigGui <- function(){
    if (length(as.integer(tkcurselection(olFileBox))))
      tkconfigure(addButton, state='normal')
    else
      tkconfigure(addButton, state='disabled')
    configList <- list(offsetSlider, offsetLab, valLab)
    dims <- sapply(fileFolder, function(x) x$file.par$number_dimensions)
    if (length(overlayList) && any(dims[match(overlayList, 
                                              names(fileFolder))] == 1)){
      for (i in configList)
        tkconfigure(i, state='normal')
      tkconfigure(offsetSlider, fg='black')
    }else{
      for (i in configList)
        tkconfigure(i, state='disabled')
      tkconfigure(offsetSlider, fg='grey')
    }
    usrSel <- 1 +	as.integer(tkcurselection(overlayBox))
    usrFiles <- overlayList[usrSel]
    configList <- list(posColButton, negColButton, projColButton, 
                       removeButton)
    if (!length(usrSel)){
      for (i in configList)
        tkconfigure(i, state='disabled')
    }else{
      oneD <- FALSE
      for (i in usrFiles){
        if (fileFolder[[i]]$file.par$number_dimensions == 1){
          oneD <- TRUE
          break
        }
      }
      twoD <- FALSE
      for (i in usrFiles){
        if (fileFolder[[i]]$file.par$number_dimensions > 1){
          twoD <- TRUE
          break
        }
      }
      tkconfigure(removeButton, state='normal')
      if (oneD)
        tkconfigure(projColButton, state='normal')
      else
        tkconfigure(projColButton, state='disabled')
      if (twoD){
        tkconfigure(posColButton, state='normal')
        tkconfigure(negColButton, state='normal')
      }else{
        tkconfigure(posColButton, state='disabled')
        tkconfigure(negColButton, state='disabled')
      }
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
  
  ####create widgets for onedFrame
  ##create file list box
  onedFileFrame <- ttklabelframe(onedFrame, text='Files')
  onedFileList <- tclVar()
  onedFileNames <- names(fileFolder)[which(sapply(fileFolder, 
                                                  function(x){x$file.par$number_dimensions}) == 1)]
  tclObj(onedFileList) <- getTitles(onedFileNames)
  onedFileBox <- tklistbox(onedFileFrame, height=11, width=30, bg='white', 
                           listvariable=onedFileList, selectmode='extended', active='dotbox',
                           exportselection=FALSE, xscrollcommand=function(...) tkset(onedXscr, ...), 
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
    selIndices <- match(onedFileNames[usrSel], names(fileFolder))
    tkselection.clear(olFileBox, 0, 'end')
    for (i in selIndices)
      tkselection.set(olFileBox, i - 1)
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
      currentSpectrum <- usrFile
      myAssign('currentSpectrum', currentSpectrum)
      refresh(multi.plot = FALSE)
      tkwm.deiconify(dlg)
      tkfocus(onedFileBox)
    }
  }
  tkbind(onedFileBox, '<Double-Button-1>', onedDouble)
  
  ##creates new shift textbox
  onedOptionFrame <- ttkframe(onedFrame)
  onedRefValFrame <- ttklabelframe(onedOptionFrame, 
                                   text='Reference value (ppm)', padding=3)
  onedRefVal <- tclVar(0)
  onedRefEntry <- ttkentry(onedRefValFrame, width=6, justify='center', 
                           textvariable=onedRefVal)
  
  ##creates get shift button
  onedLoc <- function(){
    
    ##prompt user for type of shift selection
    usr <- mySelect(c('Designated point', 'Region maximum'), multiple=FALSE, 
                    title='Get shifts at:',	preselect='Designated point', parent=dlg)
    if (length(usr) == 0 || !nzchar(usr))
      return(invisible())
    else if (usr == 'Region maximum'){
      tryCatch(shift <- regionMax(currentSpectrum)$w2, 
               error=function(er){
                 showGui()
                 refresh(multi.plot=FALSE, sub.plot=FALSE)
                 stop('Shift not defined', call.=FALSE)})
      if (is.null(shift)){
        showGui()
        refresh(multi.plot=FALSE, sub.plot=FALSE)
        stop('Shift not defined', call.=FALSE)
      }
    }else{
      
      ## Opens the main plot window if not currently opened
      if (is.na(match(2, dev.list())))
        refresh(multi.plot=FALSE, sub.plot=FALSE)
      cw(dev=2)
      
      ##gives the user instructions
      hideGui()
      cat(paste('In the main plot window:\n',  
                ' Left-click a point inside the plot to designate position\n'))
      flush.console()
      op <- par('font')
      par( font = 2 )
      legend("topleft", c('LEFT CLICK TO DESIGNATE POSITION', 
                          'RIGHT CLICK TO EXIT'),	pch=NULL, bty='n', 
             text.col=fileFolder[[wc()]]$graphics.par$fg)
      par(font = op)
      
      ##get the chemical shift at designated postion
      tryCatch(shift <- locator(1), error=function(er){
        showGui()
        refresh(multi.plot=FALSE, sub.plot=FALSE)
        stop('Shift not defined', call.=FALSE)})
      if (is.null(shift)){
        showGui()
        refresh(multi.plot=FALSE, sub.plot=FALSE)
        stop('Shift not defined', call.=FALSE)
      }
      refresh(multi.plot=FALSE, sub.plot=FALSE)
      abline(v=shift$x, lty=2, col=fileFolder[[wc()]]$graphics.par$fg)
      shift <- shift$x
    }
    if (is.na(shift[1])){
      showGui()
      err('Specified region does not contain peaks above the noise level')
    }
    tclObj(onedRefVal) <<- round(unlist(shift), 4)
    showGui()
    tkfocus(onedFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  onedLocButton <- ttkbutton(onedRefValFrame, text='Get Shift', command=onedLoc)
  
  ##creates point reference button
  onedDefRefFrame <- ttklabelframe(onedOptionFrame, text='Define reference', 
                                   padding=3)
  onedPoint <- function(){
    
    ##checks for correct input
    usrSel <- 1 + as.integer(tkcurselection(onedFileBox))
    usrFiles <- onedFileNames[usrSel]
    if (length(usrSel) == 0)
      err(paste('You must select a spectrum from the list before designating',
                ' a position'), parent=dlg)
    newShift <- suppressWarnings(as.numeric(tclvalue(onedRefVal)))
    if (is.na(newShift))
      err('You must provide a numeric value for the chemical shift reference', 
          parent=dlg)
    
    ##make sure input files are similar to each other and the current spectrum
    lineCol <- fileFolder[[wc()]]$graphics.par$fg
    for (i in usrFiles){
      if (!identical(fileFolder[[wc()]]$file.par$number_dimensions, 
                     fileFolder[[i]]$file.par$number_dimensions))
        err(paste('All files must have the same number of dimensions as the', 
                  'current spectrum'), parent=dlg)
      if (!identical(fileFolder[[wc()]]$file.par$nucleus, 
                     fileFolder[[i]]$file.par$nucleus))
        err('All files must have the same nuclei as the current spectrum', 
            parent=dlg)
      if (!identical(fileFolder[[wc()]]$file.par$matrix_size, 
                     fileFolder[[i]]$file.par$matrix_size))
        err(paste('All files must have the same number of points in all', 
                  'dimensions as the current spectrum'), parent=dlg)
    }
    
    ##makes sure the current spectrum is included in the user's selection
    if (!currentSpectrum %in% usrFiles){
      usrSel <- myMsg(paste('The current spectrum was not included in your', 
                            ' selection and will be added automatically.\nDo you wish to', 
                            ' proceed?', sep=''), 'yesno', parent=dlg)
      if (usrSel == 'no'){
        return(invisible())
      }else{
        usrFiles <- c(currentSpectrum, usrFiles)
        tkselection.set(onedFileBox, match(currentSpectrum, 
                                           onedFileNames) - 1)
      }
    }
    
    ## Opens the main plot window if not currently opened
    if (is.na(match(2, dev.list())))
      refresh(multi.plot=FALSE, sub.plot=FALSE)
    cw(dev=2)
    
    ##gives the user instructions
    hideGui()
    cat(paste('In the main plot window:\n',  
              ' Left-click a point inside the plot to define the reference\n'))
    flush.console()
    op <- par('font')
    par( font = 2 )
    legend("topleft", c('LEFT CLICK TO DEFINE REFERENCE', 
                        'RIGHT CLICK TO EXIT'),	pch=NULL, bty='n', text.col=lineCol)
    par(font = op)
    tryCatch(pointVal <- locator(1)[[1]], error=function(er){
      showGui()
      stop('Point not defined', call.=FALSE)})
    if (length(pointVal) == 0 || is.null(pointVal)){
      showGui()
      refresh(multi.plot=FALSE, sub.plot=FALSE)
      stop('Point not defined', call.=FALSE)
    }
    showGui()
    pointVal <- pointVal - newShift
    
    ##sets up and downfield shifts to the current spectrum's
    currUp <- fileFolder[[wc()]]$file.par$upfield_ppm
    currDown <- fileFolder[[wc()]]$file.par$downfield_ppm
    keepFolder <- fileFolder
    for (i in usrFiles){
      fileFolder[[i]]$file.par$upfield_ppm <- currUp
      fileFolder[[i]]$file.par$downfield_ppm <- currDown
    }
    myAssign('fileFolder', fileFolder, save.backup=FALSE)
    
    ##references selected spectra
    for (i in usrFiles){
      fileFolder[[i]]$file.par$upfield_ppm <- currUp - pointVal
      fileFolder[[i]]$file.par$downfield_ppm <- currDown - pointVal
      totShiftChange <- keepFolder[[i]]$file.par$upfield_ppm - 
        fileFolder[[i]]$file.par$upfield_ppm
      newUsr <- c(fileFolder[[i]]$graphics.par$usr[1:2] - totShiftChange, 
                  fileFolder[[i]]$graphics.par$usr[3:4])
      fileFolder[[i]]$graphics.par$usr <- newUsr
      if (!is.null(fileFolder[[i]]$peak.list)){
        fileFolder[[i]]$peak.list$w2 <- fileFolder[[i]]$peak.list$w2 - 
          totShiftChange
      }
    }
    myAssign('fileFolder', fileFolder)
    refresh()
    abline(v=newShift, lty=2, col=lineCol)
    tkfocus(onedFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  onedPointButton <- ttkbutton(onedDefRefFrame, text='Point', width=8, 
                               command=onedPoint)
  
  ##creates region reference button
  onedRegion <- function(){
    
    ##checks for correct input
    usrSel <- 1 + as.integer(tkcurselection(onedFileBox))
    usrFiles <- onedFileNames[usrSel]
    if (length(usrSel) == 0)
      err('You must select a spectrum from the list before defining a region', 
          parent=dlg)
    newShift <- suppressWarnings(as.numeric(tclvalue(onedRefVal)))
    if (is.na(newShift))
      err('You must provide a numeric value for the chemical shift reference', 
          parent=dlg)
    
    ##make sure input files are similar to each other and the current spectrum
    for (i in usrFiles){
      if (!identical(fileFolder[[wc()]]$file.par$number_dimensions, 
                     fileFolder[[i]]$file.par$number_dimensions))
        err(paste('All files must have the same number of dimensions as the', 
                  'current spectrum'), parent=dlg)
      if (!identical(fileFolder[[wc()]]$file.par$nucleus, 
                     fileFolder[[i]]$file.par$nucleus))
        err('All files must have the same nuclei as the current spectrum', 
            parent=dlg)
      if (!identical(fileFolder[[wc()]]$file.par$matrix_size, 
                     fileFolder[[i]]$file.par$matrix_size))
        err(paste('All files must have the same number of points in all', 
                  'dimensions as the current spectrum'), parent=dlg)
    }
    
    ##makes sure the current spectrum is included in the user's selection
    if (!currentSpectrum %in% usrFiles){
      usrSel <- myMsg(paste('The current spectrum was not included in your', 
                            ' selection and will be added automatically.\nDo you wish to', 
                            ' proceed?', sep=''), 'yesno', parent=dlg)
      if (usrSel == 'no'){
        return(invisible())
      }else{
        usrFiles <- c(currentSpectrum, usrFiles)
        tkselection.set(onedFileBox, match(currentSpectrum, onedFileNames) - 1)
      }
    }
    
    ##sets up and downfield shifts to the current spectrum's
    currUp <- fileFolder[[wc()]]$file.par$upfield_ppm
    currDown <- fileFolder[[wc()]]$file.par$downfield_ppm
    keepFolder <- fileFolder
    for (i in usrFiles){
      fileFolder[[i]]$file.par$upfield_ppm <- currUp
      fileFolder[[i]]$file.par$downfield_ppm <- currDown
    }
    myAssign('fileFolder', fileFolder, save.backup=FALSE)
    
    ##references selected spectra
    tryCatch(ms <- regionMax(usrFiles, redraw=FALSE), error=function(er){
      myAssign('fileFolder', keepFolder, save.backup=FALSE)
      showGui()
      stop('Region not defined', call.=FALSE)})
    if (is.null(ms)){
      myAssign('fileFolder', keepFolder, save.backup=FALSE)
      showGui()
      stop('Region not defined', call.=FALSE)
    }
    for (i in usrFiles){
      if (is.na(ms[i, 'w2']))
        next
      regionVal <- ms[i, 'w2'] - newShift
      fileFolder[[i]]$file.par$upfield_ppm <- currUp - regionVal
      fileFolder[[i]]$file.par$downfield_ppm <- currDown - regionVal
      totShiftChange <- keepFolder[[i]]$file.par$upfield_ppm - 
        fileFolder[[i]]$file.par$upfield_ppm
      newUsr <- c(fileFolder[[i]]$graphics.par$usr[1:2] - totShiftChange, 
                  fileFolder[[i]]$graphics.par$usr[3:4])
      fileFolder[[i]]$graphics.par$usr <- newUsr
      if (!is.null(fileFolder[[i]]$peak.list)){
        fileFolder[[i]]$peak.list$w2 <-	fileFolder[[i]]$peak.list$w2 - 
          totShiftChange
      }
    }
    myAssign('fileFolder', fileFolder)
    refresh()
    
    ##draw lines to indicate region maximum
    lineCol <- fileFolder[[wc()]]$graphics.par$fg
    abline(v=newShift, lty=2, col=lineCol)		
    tkfocus(onedFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  onedRegionButton <- ttkbutton(onedDefRefFrame, text='Region', width=8, 
                                command=onedRegion)
  
  ##manual shift adjustment functions
  onedAdjFrame <- ttklabelframe(onedOptionFrame, text='Man. adjustment (ppm)', 
                                padding=3)
  onedAmountVal <- tclVar(1)
  onedArrow <- function(direct, n){
    
    ##checks for correct inputs
    tryCatch({n <- as.numeric(tclvalue((onedAmountVal)))
    if (n < 0)
      warning()
    }, warning = function(w){
      err('Increment value must be a positive number', parent=dlg)
    })
    usrSel <- 1 + as.integer(tkcurselection(onedFileBox))
    usrFiles <- onedFileNames[usrSel]
    fileIndices <- match(usrFiles, names(fileFolder))
    if (length(usrSel) == 0)
      err('You must select a spectrum from the list before adjusting shifts', 
          parent=dlg)
    
    ##adjust shifts to the left
    if (direct == 'left'){
      for (i in fileIndices){
        fileFolder[[i]]$file.par$upfield_ppm <- 
          fileFolder[[i]]$file.par$upfield_ppm + n
        fileFolder[[i]]$file.par$downfield_ppm <- 
          fileFolder[[i]]$file.par$downfield_ppm + n
        fileFolder[[i]]$graphics.par$usr[1:2] <- 
          fileFolder[[i]]$graphics.par$usr[1:2] + n
        if (!is.null(fileFolder[[i]]$peak.list))
          fileFolder[[i]]$peak.list$w2 <- 
          fileFolder[[i]]$peak.list$w2 + n
      }
      
      ##adjust shifts to the right
    }else{
      for (i in fileIndices){
        fileFolder[[i]]$file.par$upfield_ppm <- 
          fileFolder[[i]]$file.par$upfield_ppm - n
        fileFolder[[i]]$file.par$downfield_ppm <- 
          fileFolder[[i]]$file.par$downfield_ppm - n
        fileFolder[[i]]$graphics.par$usr[1:2] <- 
          fileFolder[[i]]$graphics.par$usr[1:2] - n
        if (!is.null(fileFolder[[i]]$peak.list))
          fileFolder[[i]]$peak.list$w2 <- 
          fileFolder[[i]]$peak.list$w2 - n
      }
    }
    myAssign('fileFolder', fileFolder)
    refresh()	
    tkfocus(onedFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  
  ##creates manual shift adjustment arrows
  onedLeftButton <- ttkbutton(onedAdjFrame, text='<', width=4, 
                              command=function(...) onedArrow('left', onedAmountVal))
  onedAmountEntry <- ttkentry(onedAdjFrame, textvariable=onedAmountVal, width=5, 
                              justify='center')
  onedRightButton <- ttkbutton(onedAdjFrame, text='>', width=4, 
                               command=function(...) onedArrow('right', onedAmountVal))
  
  ##creates auto referencing button
  onedAuto <- function(){
    usrSel <- 1 + as.integer(tkcurselection(onedFileBox))
    usrFiles <- onedFileNames[usrSel]
    if (length(usrSel) == 0)
      err('You must select a spectrum from the list before referencing shifts.', 
          parent=dlg)
    autoRef(usrFiles)
    tclObj(onedRefVal) <<- 0
    tkfocus(onedFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  onedAutoButton <- ttkbutton(onedOptionFrame, text='Auto Ref', width=10, 
                              command=onedAuto)
  
  ##creates default button
  onedDefault <- function(){
    
    ##checks for correct input
    usrSel <- 1 + as.integer(tkcurselection(onedFileBox))
    usrFiles <- onedFileNames[usrSel]
    if (length(usrSel) == 0)
      err(paste('You must select a spectrum from the list before restoring',
                'defaults'), parent=dlg)
    
    ##restores defaults
    for (i in usrFiles){
      filePar <- ucsfHead(i, print.info=FALSE)
      filePar <- filePar[[1]]
      prevFullUsr <- c(fileFolder[[i]]$file.par$downfield_ppm, 
                       fileFolder[[i]]$file.par$upfield_ppm, 
                       fileFolder[[i]]$file.par$zero_offset - 
                         (fileFolder[[i]]$file.par$max_intensity - 
                            fileFolder[[i]]$file.par$zero_offset) * 
                         globalSettings$position.1D, 
                       fileFolder[[i]]$file.par$max_intensity)
      defUsr <- c(filePar$downfield_ppm, filePar$upfield_ppm, 
                  filePar$zero_offset - (filePar$max_intensity - filePar$zero_offset) * 
                    globalSettings$position.1D, filePar$max_intensity)
      usrDiff <- prevFullUsr - defUsr
      fileFolder[[i]]$graphics.par$usr <- 
        fileFolder[[i]]$graphics.par$usr - usrDiff
      fileFolder[[i]]$file.par$upfield_ppm <- filePar$upfield_ppm
      fileFolder[[i]]$file.par$downfield_ppm <- filePar$downfield_ppm
      if (!is.null(fileFolder[[i]]$peak.list))
        fileFolder[[i]]$peak.list$w2 <- 
        fileFolder[[i]]$peak.list$w2 - usrDiff[1]
    }
    myAssign('fileFolder', fileFolder)
    refresh()
    tkfocus(onedFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  onedDefaultButton <- ttkbutton(onedOptionFrame, text='Default', width=10, 
                                 command=onedDefault)
  
  ##creates undo button
  onedUndo<- function(){
    ud()
    tkfocus(onedFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  onedUndoButton <- ttkbutton(onedOptionFrame, text='Undo', width=10, 
                              command=onedUndo)
  
  ##creates redo button
  onedRedo<- function(){
    rd()
    tkfocus(onedFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  onedRedoButton <- ttkbutton(onedOptionFrame, text='Redo', width=10, 
                              command=onedRedo)
  
  ##add widgets to fileFrame
  tkgrid(onedFileFrame, column=1, row=1, sticky='nswe', pady=c(6, 0),	padx=8)
  tkgrid(onedFileBox, column=1, row=1, sticky='nswe')
  tkgrid(onedYscr, column=2, row=1, sticky='ns')
  tkgrid(onedXscr, column=1, row=2, sticky='we')
  
  ##make fileFrame stretch when window is resized
  tkgrid.columnconfigure(onedFrame, 1, weight=1)
  tkgrid.rowconfigure(onedFrame, 1, weight=10)
  tkgrid.columnconfigure(onedFileFrame, 1, weight=1)
  tkgrid.rowconfigure(onedFileFrame, 1, weight=1)
  
  ##add widgets to optionFrame
  tkgrid(onedOptionFrame, column=2, row=1, sticky='nswe', pady=c(10, 2), 
         padx=c(4, 0))
  tkgrid(onedRefValFrame, column=1, columnspan=2, row=2, sticky='we', pady=4)
  tkgrid(onedRefEntry, column=1, row=1, padx=4)
  tkgrid(onedLocButton, column=2, row=1)
  
  tkgrid(onedDefRefFrame, column=1, columnspan=2, row=3, sticky='we', pady=4)
  tkgrid(onedPointButton, column=1, row=1, padx=4)
  tkgrid(onedRegionButton, column=2, row=1)
  
  tkgrid(onedAdjFrame, column=1, columnspan=2, row=4, sticky='we', pady=c(0, 6))
  tkgrid(onedLeftButton, column=1, row=1, padx=c(12, 0))
  tkgrid(onedAmountEntry, column=2, row=1)
  tkgrid(onedRightButton, column=3, row=1)
  
  tkgrid(onedAutoButton, column=1, row=5, pady=c(6, 2))
  tkgrid(onedDefaultButton, column=2, row=5, pady=c(6, 2))
  tkgrid(onedUndoButton, column=1, row=6)
  tkgrid(onedRedoButton, column=2, row=6)
  
  ##make optionFrame stretch when window is resized
  tkgrid.rowconfigure(onedOptionFrame, 0, weight=1)
  tkgrid.rowconfigure(onedOptionFrame, 7, weight=1)
  
  ##resets file list whenever the mouse enters the GUI
  onedMouse <- function(){
    reset(onedFileList, onedFileBox, onedFileNames, dims='1D')
    onedFileNames <<- names(fileFolder)[which(sapply(fileFolder, 
                                                     function(x){x$file.par$number_dimensions}) == 1)]
  }
  tkbind(onedFrame, '<Enter>', onedMouse)
  tkbind(onedFrame, '<FocusIn>', onedMouse)
  
  ####create widgets for twodFrame
  ##create file list box
  twodFileFrame <- ttklabelframe(twodFrame, text='Files')
  twodFileList <- tclVar()
  twodFileNames <- names(fileFolder)[which(sapply(fileFolder, 
                                                  function(x){x$file.par$number_dimensions}) > 1)]
  tclObj(twodFileList) <- getTitles(twodFileNames)
  twodFileBox <- tklistbox(twodFileFrame, height=12, width=30, bg='white', 
                           listvariable=twodFileList, selectmode='extended', active='dotbox',	
                           exportselection=FALSE, xscrollcommand=function(...) tkset(twodXscr, ...), 
                           yscrollcommand=function(...) tkset(twodYscr, ...))
  twodXscr <- ttkscrollbar(twodFileFrame, orient='horizontal', 
                           command=function(...) tkxview(twodFileBox, ...))
  twodYscr <- ttkscrollbar(twodFileFrame, orient='vertical', 
                           command=function(...) tkyview(twodFileBox, ...))
  if (length(twodFileNames) > 2){
    for (i in seq(0, length(twodFileNames) - 1, 2))
      tkitemconfigure(twodFileBox, i, background='#ececff')
  }
  if (fileFolder[[current]]$file.par$number_dimensions > 1){
    tkselection.set(twodFileBox, match(currentSpectrum, twodFileNames) - 1)
    tcl(twodFileBox, 'see', match(currentSpectrum, twodFileNames) - 1)
  }
  
  ##export fileBox selections to other tabs
  twodSelect <- function(){
    usrSel <- 1 + as.integer(tkcurselection(twodFileBox))
    usrFile <- twodFileNames[usrSel]
    selIndices <- match(twodFileNames[usrSel], names(fileFolder))
    tkselection.clear(olFileBox, 0, 'end')
    for (i in selIndices)
      tkselection.set(olFileBox, i - 1)
  }
  tkbind(twodFileBox, '<<ListboxSelect>>', twodSelect)
  
  ##switches spectra on left-mouse double-click
  twodDouble <- function(){
    usrSel <- 1 + as.integer(tkcurselection(twodFileBox))
    if (length(usrSel))
      usrFile <- twodFileNames[usrSel]
    else
      usrFile <- NULL
    if (!is.null(usrFile) && currentSpectrum != usrFile){
      currentSpectrum <- usrFile
      myAssign('currentSpectrum', currentSpectrum)
      refresh(multi.plot = FALSE)
      tkwm.deiconify(dlg)
      tkfocus(twodFileBox)
    }
  }
  tkbind(twodFileBox, '<Double-Button-1>', twodDouble)
  
  ##creates reference value textboxes
  twodOptionFrame <- ttkframe(twodFrame)
  twodRefValFrame <- ttklabelframe(twodOptionFrame, text='Reference value (ppm)', 
                                   padding=3)
  w1RefVal <- tclVar(0)
  w1Entry <- ttkentry(twodRefValFrame, width=6, justify='center', 
                      textvariable=w1RefVal)
  w2RefVal <- tclVar(0)
  w2Entry <- ttkentry(twodRefValFrame, width=6, justify='center', 
                      textvariable=w2RefVal)
  
  ##creates get shift button
  twodLoc <- function(){
    reset(twodFileList, twodFileBox, twodFileNames, dims='2D')
    twodFileNames <<- names(fileFolder)[which(sapply(fileFolder, 
                                                     function(x){x$file.par$number_dimensions}) > 1)]
    
    ##prompt user for type of shift selection
    usr <- mySelect(c('Designated point', 'Region maximum'), multiple=FALSE, 
                    title='Get shifts at:',	preselect='Designated point', parent=dlg)
    if (length(usr) == 0 || !nzchar(usr))
      return(invisible())
    else if (usr == 'Region maximum'){
      tryCatch(shift <- regionMax(currentSpectrum)[c('w2', 
                                                     'w1')], error=function(er){
                                                       showGui()
                                                       refresh(multi.plot=FALSE, sub.plot=FALSE)
                                                       stop('Shift not defined', call.=FALSE)})
      if (is.null(shift)){
        showGui()
        refresh(multi.plot=FALSE, sub.plot=FALSE)
        stop('Shift not defined', call.=FALSE)
      }
    }else{
      
      ## Opens the main plot window if not currently opened
      if (is.na(match(2, dev.list())))
        refresh(multi.plot=FALSE, sub.plot=FALSE)
      cw(dev=2)
      
      ##gives the user instructions
      hideGui()
      cat(paste('In the main plot window:\n',  
                ' Left-click a point inside the plot to designate position\n'))
      flush.console()
      op <- par('font')
      par(font=2)
      legend("topleft", c('LEFT CLICK TO DESIGNATE POSITION', 
                          'RIGHT CLICK TO EXIT'),	pch=NULL, bty='n', 
             text.col=fileFolder[[wc()]]$graphics.par$fg)
      par(font=op)
      
      ##get the chemical shift at designated postion
      tryCatch(shift <- locator(1), error=function(er){
        showGui()
        refresh(multi.plot=FALSE, sub.plot=FALSE)
        stop('Shift not defined', call.=FALSE)})
      if (is.null(shift)){
        showGui()
        refresh(multi.plot=FALSE, sub.plot=FALSE)
        stop('Shift not defined', call.=FALSE)
      }
    }
    refresh(multi.plot=FALSE, sub.plot=FALSE)
    if (is.na(shift[2])){
      showGui()
      err('Specified region does not contain peaks above the noise level')
    }
    points(shift[1], shift[2], pch='*', cex=1.5, 
           col=fileFolder[[wc()]]$graphics.par$fg)
    tclObj(w1RefVal) <<- round(unlist(shift[2]), 4)
    tclObj(w2RefVal) <<- round(unlist(shift[1]), 4)
    showGui()
    tkfocus(twodFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  twodLocButton <- ttkbutton(twodRefValFrame, text='Get Shifts', width=14, 
                             command=twodLoc)
  
  ##creates point reference button
  twodDefRefFrame <- ttklabelframe(twodOptionFrame, text='Define reference', 
                                   padding=3)
  twodPoint <- function(){
    
    ##checks for correct input
    usrSel <- 1 + as.integer(tkcurselection(twodFileBox))
    usrFiles <- twodFileNames[usrSel]
    if (length(usrSel) == 0)
      err(paste('You must select a spectrum from the list before designating a',
                'position'), parent=dlg)
    w1RefVal <- suppressWarnings(as.numeric(tclvalue(w1RefVal)))
    w2RefVal <- suppressWarnings(as.numeric(tclvalue(w2RefVal)))
    if (is.na(w1RefVal) || is.na(w2RefVal))
      err('You must provide numeric values for the chemical shift reference', 
          parent=dlg)
    
    ##make sure input files are similar to each other and the current spectrum
    for (i in usrFiles){
      if (!identical(fileFolder[[wc()]]$file.par$number_dimensions, 
                     fileFolder[[i]]$file.par$number_dimensions))
        err(paste('All files must have the same number of dimensions as the', 
                  'current spectrum'), parent=dlg)
      if (!identical(fileFolder[[wc()]]$file.par$nucleus, 
                     fileFolder[[i]]$file.par$nucleus))
        err('All files must have the same nuclei as the current spectrum', 
            parent=dlg)
      if (!identical(fileFolder[[wc()]]$file.par$matrix_size, 
                     fileFolder[[i]]$file.par$matrix_size))
        err(paste('All files must have the same number of points in all', 
                  'dimensions as the current spectrum'), parent=dlg)
    }
    
    ##makes sure the current spectrum is included in the user's selection
    if (!currentSpectrum %in% usrFiles){
      usrSel <- myMsg(paste('The current spectrum was not included in your', 
                            ' selection and will be added automatically.\nDo you wish to', 
                            ' proceed?', sep=''), 'yesno', parent=dlg)
      if (usrSel == 'no'){
        return(invisible())
      }else{
        usrFiles <- c(currentSpectrum, usrFiles)
        tkselection.set(twodFileBox, match(currentSpectrum, twodFileNames) - 1)
      }
    }
    
    ## Opens the main plot window if not currently opened
    if (is.na(match(2, dev.list())))
      refresh(multi.plot=FALSE, sub.plot=FALSE)
    cw(dev=2)
    
    ##gives the user instructions
    hideGui()
    cat(paste('In the main plot window:\n',  
              ' Left-click a point inside the plot to define the reference\n'))
    flush.console()
    op <- par('font')
    par(font=2)
    legend("topleft", c('LEFT CLICK TO DEFINE REFERENCE', 
                        'RIGHT CLICK TO EXIT'), pch=NULL, bty='n', 
           text.col=fileFolder[[wc()]]$graphics.par$fg)
    par(font = op)
    tryCatch(pointVal <- locator(1), error=function(er){
      showGui()
      stop('Point not defined', call.=FALSE)})
    if (length(pointVal) == 0 || is.null(pointVal)){
      showGui()
      refresh(multi.plot=FALSE, sub.plot=FALSE)
      stop('Point not defined', call.=FALSE)
    }
    showGui()
    pointVal <- c(pointVal[[2]] - w1RefVal, pointVal[[1]] - w2RefVal)
    
    ##sets up and downfield shifts to the current spectrum's
    currUp <- fileFolder[[wc()]]$file.par$upfield_ppm
    currDown <- fileFolder[[wc()]]$file.par$downfield_ppm
    keepFolder <- fileFolder
    for (i in usrFiles){
      fileFolder[[i]]$file.par$upfield_ppm <- currUp
      fileFolder[[i]]$file.par$downfield_ppm <- currDown
    }
    myAssign('fileFolder', fileFolder, save.backup=FALSE)
    
    ##references selected spectra
    for (i in usrFiles){
      fileFolder[[i]]$file.par$upfield_ppm <- currUp - pointVal
      fileFolder[[i]]$file.par$downfield_ppm <- currDown - pointVal
      totShiftChange <- keepFolder[[i]]$file.par$upfield_ppm - 
        fileFolder[[i]]$file.par$upfield_ppm
      newUsr <- c(fileFolder[[i]]$graphics.par$usr[1:2] - totShiftChange[2], 
                  fileFolder[[i]]$graphics.par$usr[3:4] - totShiftChange[1])
      fileFolder[[i]]$graphics.par$usr <- newUsr
      if (!is.null(fileFolder[[i]]$peak.list)){
        fileFolder[[i]]$peak.list$w1 <- fileFolder[[i]]$peak.list$w1 - 
          totShiftChange[1]
        fileFolder[[i]]$peak.list$w2 <- fileFolder[[i]]$peak.list$w2 - 
          totShiftChange[2]
      }
    }
    myAssign('fileFolder', fileFolder)
    refresh()
    points(w2RefVal, w1RefVal, pch='*', cex=1.5, 
           col=fileFolder[[wc()]]$graphics.par$fg)
    showGui()
    tkfocus(twodFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  twodPointButton <- ttkbutton(twodDefRefFrame, text='Point', width=8, 
                               command=twodPoint)
  
  ##creates region reference button
  twodRegion <- function(){
    
    ##checks for correct input
    reset(twodFileList, twodFileBox, twodFileNames, dims='2D')
    twodFileNames <<- names(fileFolder)[which(sapply(fileFolder, 
                                                     function(x){x$file.par$number_dimensions}) > 1)]
    usrSel <- 1 + as.integer(tkcurselection(twodFileBox))
    usrFiles <- twodFileNames[usrSel]
    if (length(usrSel) == 0)
      err('You must select a spectrum from the list before defining a region', 
          parent=twodFrame)
    w1RefVal <- suppressWarnings(as.numeric(tclvalue(w1RefVal)))
    w2RefVal <- suppressWarnings(as.numeric(tclvalue(w2RefVal)))
    if (is.na(w1RefVal) || is.na(w2RefVal))
      err('You must provide numeric values for the chemical shift reference', 
          parent=dlg)
    
    ##make sure input files are similar to each other and the current spectrum
    for (i in usrFiles){
      if (!identical(fileFolder[[wc()]]$file.par$number_dimensions, 
                     fileFolder[[i]]$file.par$number_dimensions))
        err(paste('All files must have the same number of dimensions as the', 
                  'current spectrum'), parent=dlg)
      if (!identical(fileFolder[[wc()]]$file.par$nucleus, 
                     fileFolder[[i]]$file.par$nucleus))
        err('All files must have the same nuclei as the current spectrum',
            parent=dlg)
      if (!identical(fileFolder[[wc()]]$file.par$matrix_size, 
                     fileFolder[[i]]$file.par$matrix_size))
        err(paste('All files must have the same number of points in all', 
                  'dimensions as the current spectrum'), parent=dlg)
    }
    
    ##makes sure the current spectrum is included in the user's selection
    if (!currentSpectrum %in% usrFiles){
      usrSel <- myMsg(paste('The current spectrum was not included in your', 
                            ' selection and will be added automatically.\nDo you wish to', 
                            ' proceed?', sep=''), 'yesno', parent=dlg)
      if (usrSel == 'no'){
        return(invisible())
      }else{
        usrFiles <- c(currentSpectrum, usrFiles)
        tkselection.set(twodFileBox, match(currentSpectrum, twodFileNames) - 1)
      }
    }
    
    ##sets up and downfield shifts to the current spectrum's
    currUp <- fileFolder[[wc()]]$file.par$upfield_ppm
    currDown <- fileFolder[[wc()]]$file.par$downfield_ppm
    keepFolder <- fileFolder
    for (i in usrFiles){
      fileFolder[[i]]$file.par$upfield_ppm <- currUp
      fileFolder[[i]]$file.par$downfield_ppm <- currDown
    }
    myAssign('fileFolder', fileFolder, save.backup=FALSE)
    
    ##references selected spectra
    tryCatch(ms <- regionMax(usrFiles, redraw=FALSE), error=function(er){
      myAssign('fileFolder', keepFolder, save.backup=FALSE)
      showGui()
      stop('Region not defined', call.=FALSE)})
    if (is.null(ms)){
      myAssign('fileFolder', keepFolder, save.backup=FALSE)
      showGui()
      stop('Region not defined', call.=FALSE)
    }
    for (i in usrFiles){
      if (is.na(ms[i, 'w1']) || is.na(ms[i, 'w2']))
        next
      regionVal <- c(ms[i, 'w1'] - w1RefVal,	ms[i, 'w2'] - w2RefVal)
      fileFolder[[i]]$file.par$upfield_ppm <- currUp - regionVal
      fileFolder[[i]]$file.par$downfield_ppm <- currDown - regionVal
      totShiftChange <- keepFolder[[i]]$file.par$upfield_ppm - 
        fileFolder[[i]]$file.par$upfield_ppm
      newUsr <- c(fileFolder[[i]]$graphics.par$usr[1:2] - totShiftChange[2], 
                  fileFolder[[i]]$graphics.par$usr[3:4] - totShiftChange[1])
      fileFolder[[i]]$graphics.par$usr <- newUsr
      if (!is.null(fileFolder[[i]]$peak.list)){
        fileFolder[[i]]$peak.list$w1 <- fileFolder[[i]]$peak.list$w1 - 
          totShiftChange[1]
        fileFolder[[i]]$peak.list$w2 <-	fileFolder[[i]]$peak.list$w2 - 
          totShiftChange[2]
      }
    }
    myAssign('fileFolder', fileFolder)
    refresh()
    
    ##draw lines to indicate region maximum
    lineCol <- fileFolder[[wc()]]$graphics.par$fg
    abline(h=w1RefVal,	v=w2RefVal, lty=2,	col=lineCol)
    tkfocus(twodFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }	
  twodRegionButton <- ttkbutton(twodDefRefFrame, text='Region', width=8, 
                                command=twodRegion)
  
  ##create manual shift adjustment arrows
  twodAdjFrame <- ttklabelframe(twodOptionFrame, text='Man. adjustment (ppm)', 
                                padding=3)
  twodAmountVal <- tclVar(1)
  twodArrow <- function(direct, n){
    reset(twodFileList, twodFileBox, twodFileNames, dims='2D')
    twodFileNames <<- names(fileFolder)[which(sapply(fileFolder, 
                                                     function(x){x$file.par$number_dimensions}) > 1)]
    
    ##checks for valid inputs
    tryCatch({n <- as.numeric(tclvalue((twodAmountVal)))
    if (n < 0)
      warning()
    }, warning = function(w){
      err('Increment value must be a positive number', parent=dlg)
    })
    usrSel <- 1 + as.integer(tkcurselection(twodFileBox))
    usrFiles <- twodFileNames[usrSel]
    if (length(usrSel) == 0)
      err('You must select a spectrum from the list before adjusting shifts', 
          parent=dlg)
    
    ##adjust shifts to the left
    if (direct == 'left'){
      ppmInc <- c(0, n)
      usrInc <- c(n, n, 0, 0)
      peakInc <- c('w2', n)
      
      ##adjust shifts to the right
    }else if (direct == 'right'){
      ppmInc <- c(0, -n, 0, -n)
      usrInc <- c(-n, -n, 0, 0)
      peakInc <- c('w2', -n)
      
      ##adjust shifts downward
    }else if (direct == 'down'){
      ppmInc <- c(n, 0)
      usrInc <- c(0, 0, n, n)
      peakInc <- c('w1', n)
      
      ##adjust shifts upward
    }else if (direct == 'up'){
      ppmInc <- c(-n, 0)
      usrInc <- c(0, 0, -n, -n)
      peakInc <- c('w1', -n)
    }
    
    ##set new shifts
    for (i in usrFiles){
      fileFolder[[i]]$file.par$upfield_ppm <- 
        fileFolder[[i]]$file.par$upfield_ppm + ppmInc
      fileFolder[[i]]$file.par$downfield_ppm <- 
        fileFolder[[i]]$file.par$downfield_ppm + ppmInc
      fileFolder[[i]]$graphics.par$usr <- 
        fileFolder[[i]]$graphics.par$usr + usrInc
      if (!is.null(fileFolder[[i]]$peak.list))
        fileFolder[[i]]$peak.list[, peakInc[1]] <- 
        fileFolder[[i]]$peak.list[, peakInc[1]] + as.numeric(peakInc[2])
    }
    
    myAssign('fileFolder', fileFolder)
    refresh()	
    tkfocus(twodFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }	
  twodUpButton <- ttkbutton(twodAdjFrame, text='^', width=5, 
                            command=function(...) twodArrow('up', twodAmountVal))
  twodDownButton <- ttkbutton(twodAdjFrame, text='v', width=5, 
                              command=function(...) twodArrow('down', twodAmountVal))
  twodLeftButton <- ttkbutton(twodAdjFrame, text='<', width=4, 
                              command=function(...) twodArrow('left', twodAmountVal))
  twodRightButton <- ttkbutton(twodAdjFrame, text='>', width=4, 
                               command=function(...) twodArrow('right', twodAmountVal))
  twodAmountEntry <- ttkentry(twodAdjFrame, textvariable=twodAmountVal, width=5, 
                              justify='center')
  
  ##creates default button
  twodDefault <- function(){
    
    ##checks for correct input
    reset(twodFileList, twodFileBox, twodFileNames, dims='2D')
    twodFileNames <<- names(fileFolder)[which(sapply(fileFolder, 
                                                     function(x){x$file.par$number_dimensions}) > 1)]
    usrSel <- 1 + as.integer(tkcurselection(twodFileBox))
    usrFiles <- twodFileNames[usrSel]
    if (length(usrSel) == 0)
      err(paste('You must select a spectrum from the list before restoring', 
                'defaults'), parent=dlg)
    
    ##restores defaults
    for (i in usrFiles){
      filePar <- ucsfHead(i, print.info=FALSE)[[1]]
      prevFullUsr <- c(fileFolder[[i]]$file.par$downfield_ppm[2],	
                       fileFolder[[i]]$file.par$upfield_ppm[2], 
                       fileFolder[[i]]$file.par$downfield_ppm[1], 
                       fileFolder[[i]]$file.par$upfield_ppm[1])
      defUsr <- c(filePar$downfield_ppm[2],	filePar$upfield_ppm[2], 
                  filePar$downfield_ppm[1], filePar$upfield_ppm[1])
      usrDiff <- prevFullUsr - defUsr
      fileFolder[[i]]$graphics.par$usr <- 
        fileFolder[[i]]$graphics.par$usr - usrDiff
      fileFolder[[i]]$file.par$upfield_ppm <- filePar$upfield_ppm
      fileFolder[[i]]$file.par$downfield_ppm <- filePar$downfield_ppm
      if (!is.null(fileFolder[[i]]$peak.list)){
        fileFolder[[i]]$peak.list$w1 <- 
          fileFolder[[i]]$peak.list$w1 - usrDiff[3]
        fileFolder[[i]]$peak.list$w2 <- 
          fileFolder[[i]]$peak.list$w2 - usrDiff[1]
      }
    }
    myAssign('fileFolder', fileFolder)
    refresh()
    tkfocus(twodFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  twodDefaultButton <- ttkbutton(twodOptionFrame, text='Default', width=18,
                                 command=twodDefault)
  
  ##creates undo button
  twodUndo<- function(){
    ud()
    tkfocus(twodFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  twodUndoButton <- ttkbutton(twodOptionFrame, text='Undo', width=10, 
                              command=twodUndo)
  
  ##creates redo button
  twodRedo<- function(){
    rd()
    tkfocus(twodFrame)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  twodRedoButton <- ttkbutton(twodOptionFrame, text='Redo', width=10, 
                              command=twodRedo)
  
  ##add widgets to fileFrame
  tkgrid(twodFileFrame, column=1, row=1, sticky='nswe', pady=c(6, 0),	padx=8)
  tkgrid(twodFileBox, column=1, row=1, sticky='nswe')
  tkgrid(twodYscr, column=2, row=1, sticky='ns')
  tkgrid(twodXscr, column=1, row=2, sticky='we')
  
  ##make fileFrame stretch when window is resized
  tkgrid.columnconfigure(twodFrame, 1, weight=1)
  tkgrid.rowconfigure(twodFrame, 1, weight=10)
  tkgrid.columnconfigure(twodFileFrame, 1, weight=1)
  tkgrid.rowconfigure(twodFileFrame, 1, weight=1)
  
  ##add widgets to optionFrame
  tkgrid(twodOptionFrame, column=2, row=1, sticky='nswe', pady=c(10, 2), 
         padx=c(4, 0))
  tkgrid(twodRefValFrame, column=1, columnspan=2, row=2, sticky='we', pady=4)
  tkgrid(ttklabel(twodRefValFrame, text='w1:'), column=1, row=1)
  tkgrid(w1Entry, column=2, row=1, padx=1)
  tkgrid(ttklabel(twodRefValFrame, text='w2:'), column=3, row=1)
  tkgrid(w2Entry, column=4, row=1, padx=c(1, 3))
  tkgrid(twodLocButton, column=1, columnspan=4, row=3, pady=c(3, 0))
  
  tkgrid(twodDefRefFrame, column=1, columnspan=2, row=3, sticky='we', pady=4)
  tkgrid(twodPointButton, column=1, row=1, padx=c(6, 4))
  tkgrid(twodRegionButton, column=2, row=1)
  
  tkgrid(twodAdjFrame, column=1, columnspan=2, row=4, sticky='we', pady=c(0, 6))
  tkgrid(twodUpButton, column=2, row=1, sticky='s', pady=c(2, 0))
  tkgrid(twodLeftButton, column=1, row=2, sticky='e', padx=c(11, 0))
  tkgrid(twodAmountEntry, column=2, row=2)
  tkgrid(twodRightButton, column=3, row=2, sticky='w')
  tkgrid(twodDownButton, column=2, row=3, sticky='n', pady=c(0, 4))
  
  tkgrid(twodDefaultButton, column=1, columnspan=2, row=5, pady=c(0, 2))
  tkgrid(twodUndoButton, column=1, row=6)
  tkgrid(twodRedoButton, column=2, row=6)
  
  ##make optionFrame stretch when window is resized
  tkgrid.rowconfigure(twodOptionFrame, 0, weight=1)
  tkgrid.rowconfigure(twodOptionFrame, 7, weight=1)
  
  ##resets file list whenever the mouse enters the GUI
  twodMouse <- function(){
    reset(twodFileList, twodFileBox, twodFileNames, dims='2D')
    twodFileNames <<- names(fileFolder)[which(sapply(fileFolder, 
                                                     function(x){x$file.par$number_dimensions}) > 1)]
  }
  tkbind(twodFrame, '<Enter>', twodMouse)
  tkbind(twodFrame, '<FocusIn>', twodMouse)
  
  ##Allows users to press the 'Enter' key to make selections
  onEnter <- function(){
    focus <- as.character(tkfocus())
    if (length(grep('.1.1.1.1$', focus)))
      olDouble(olFileBox)
    else if (length(grep('.1.1.3.1$', focus)))
      olDouble(overlayBox)
    else if (length(grep('.1.2.1.1$', focus)))
      onedDouble()
    else if (length(grep('.1.3.1.1$', focus)))
      twodDouble()
    else
      tryCatch(tkinvoke(focus), error=function(er){})
  }
  tkbind(dlg, '<Return>', onEnter) 
  
  ##Change the window title depending on which pane is displayed
  onSwitch <- function(){
    if (length(grep('.1.1', as.character(tkselect(osBook)))))
      tkwm.title(dlg, 'Overlays')
    else
      tkwm.title(dlg, 'Referencing')
  }
  
  ##enable\disable panes depending on which files are open
  onMouse <- function(){
    if (length(names(fileFolder)) && any(which(sapply(fileFolder, 
                                                      function(x){x$file.par$number_dimensions}) == 1)))
      tcl(osBook, 'tab', 1, state='normal')
    else
      tcl(osBook, 'tab', 1, state='disabled')
    if (length(names(fileFolder)) && any(which(sapply(fileFolder, 
                                                      function(x){x$file.par$number_dimensions}) > 1)))
      tcl(osBook, 'tab', 2, state='normal')
    else
      tcl(osBook, 'tab', 2, state='disabled')
  }
  tkbind(dlg, '<Enter>', onMouse)
  tkbind(dlg, '<FocusIn>', onMouse)
  tkbind(dlg, '<<NotebookTabChanged>>', onSwitch)

# Now that the window is set up, deiconify it to show the window properly
  tkwm.deiconify(dlg) # New
  tkraise(dlg) # New
  tcl("update") # New
  
  invisible()
}

#' Interactive GUI for zooming and scrolling
#' Zoom and Manipulate View using zoom
#'
#' This function opens a graphical user interface (GUI) that allows users to manipulate the zoom and view of a plot.
#' Users can choose between zooming or scrolling mode, adjust the zoom increment, and perform various zooming actions,
#' such as zooming in, zooming out, and focusing on specific regions of the plot. The GUI also provides options for
#' resetting the view to full, centering the view, and more.
#'
#' @return NULL (The function primarily works with GUI elements.)
#'
#' @import tcltk2
#'
#' @export
zoom <- function(){
  
  ##Checks for open files
  wc()
  
  ##creates main window
  tclCheck()

 # Destroy any existing window with the same ID # New
  if (as.logical(tcl('winfo', 'exists', '.zm'))) {
    tkdestroy('.zm')  # Destroy previous window to avoid conflicts
  }

  #dlg <- myToplevel('zm', pady=2, padx=4)
  #if (is.null(dlg))
  #  return(invisible())
  dlg <- tktoplevel() # New
  tkwm.title(dlg, 'Zoom')
  tkwm.geometry(dlg, "250x100+200+200")
  tkwm.resizable(dlg, FALSE, FALSE)

  # Withdraw the window to prevent flickering while setting up # New
  tkwm.withdraw(dlg)
	
  tcl('wm', 'attributes', dlg, topmost=TRUE)
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  
  ##create radiobuttons
  leftFrame <- ttkframe(dlg)
  radioFrame <- ttkframe(leftFrame)
  rbVal <- tclVar('zoom')
  onZoom <- function(){
    tkconfigure(upButton, text='In')
    tkconfigure(downButton, text='Out')
    tkconfigure(leftButton, state='disabled')
    tkconfigure(rightButton, state='disabled')
    tkfocus(arrowFrame)
  }
  zoomButton <- ttkradiobutton(radioFrame, variable=rbVal, value='zoom', 
                               text='zoom', command=onZoom)
  
  onScroll <- function(){
    tkconfigure(upButton, text='^')
    tkconfigure(downButton, text='v')
    tkconfigure(leftButton, text='<', state='normal')
    tkconfigure(rightButton, text='>', state='normal')
    tkfocus(arrowFrame)
  }
  scrollButton <- ttkradiobutton(radioFrame, variable=rbVal, value='scroll', 
                                 text='scroll', command=onScroll)
  
  ##create arrow buttons
  arrowFrame <- ttkframe(leftFrame)
  inc <- tclVar(20)
  onArrow <- function(type, direct, n){
    type=tclvalue(type)
    tryCatch({n <- as.numeric(tclvalue((inc)))
    if (n < 0)
      warning()
    }, warning = function(w){
      err('Increment value must be a positive number', parent=dlg)
    })
    if (type == 'scroll'){
      switch(direct, 'up'=pu(n), 
             'down'=pd(n), 
             'left'=pl(n), 
             'right'=pr(n))
    }else{
      switch(direct, 'up'=zi(n), 'down'=zo(n)) 
    }	
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  upButton <- ttkbutton(arrowFrame, text='In', width=5, command=function(...) 
    onArrow(rbVal, 'up', inc))
  downButton <- ttkbutton(arrowFrame, text='Out', width=5, command=function(...)
    onArrow(rbVal, 'down', inc))
  leftButton <- ttkbutton(arrowFrame, text='<', width=4, state='disabled', 
                          command=function(...) onArrow(rbVal, 'left', inc))
  rightButton <- ttkbutton(arrowFrame, text='>', width=4, state='disabled', 
                           command=function(...) onArrow(rbVal, 'right', inc))
  editEntry <- ttkentry(arrowFrame, textvariable=inc, width=5, justify='center')
  
  ##create other zoom buttons
  rightFrame <- ttkframe(dlg)
  onFf <- function(){
    ff()
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  ffButton <- ttkbutton(rightFrame, text='Full', width=8, command=onFf)
  
  onZc <- function(){
    zc()
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  zcButton <- ttkbutton(rightFrame, text='Center', width=8, command=onZc)
  
  onZz <- function(){
    zz()
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  zzButton <- ttkbutton(rightFrame, text='Hand', width=8, command=onZz)
  
  onPz <- function(){
    pz()
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  pzButton <- ttkbutton(rightFrame, text='Point', width=8, command=onPz)
  
  onZp <- function(){
    zp()
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  zpButton <- ttkbutton(rightFrame, text='Prev', width=8, command=onZp)
  
  onLoc <- function(){
    hideGui()
    usrLoc <- mySelect(c('Chemical Shift', 'Region Maximum', 'Delta in PPM/Hz'), 
                       title = 'Measure:', parent=dlg)
    if (!nzchar(usrLoc)){
      showGui()
      return(invisible())
    }
    if (usrLoc == 'Chemical Shift')
      loc()
    else if (usrLoc == 'Region Maximum'){
      tryCatch(shift <- regionMax(currentSpectrum), 
               error=function(er){
                 showGui()
                 refresh(multi.plot=FALSE, sub.plot=FALSE)
                 stop('Shift not defined', call.=FALSE)})
      showGui()
      if (is.null(shift)){
        refresh(multi.plot=FALSE, sub.plot=FALSE)
        stop('Shift not defined', call.=FALSE)
      }
      rownames(shift) <- NULL
      log_message(shift)
    }else
      delta()
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    bringFocus()
  }
  locButton <- ttkbutton(rightFrame, text='Get shifts', width=8, command=onLoc)
  
  ##add widgets to leftFrame
  tkgrid(leftFrame, column=1, row=1, padx=c(0, 4))
  tkgrid(radioFrame, column=1, columnspan=3, row=1)
  tkgrid(zoomButton, column=1, row=1, padx=3)
  tkgrid(scrollButton, column=2, row=1)
  
  tkgrid(arrowFrame, column=1, columnspan=3, row=2)
  tkgrid(upButton, column=2, row=1, sticky='s', pady=c(2, 0))
  tkgrid(leftButton, column=1, row=2, sticky='e')
  tkgrid(editEntry, column=2, row=2)
  tkgrid(rightButton, column=3, row=2, sticky='w')
  tkgrid(downButton, column=2, row=3, sticky='n', pady=c(0, 4))
  
  ##add widgets to rightFrame
  tkgrid(rightFrame, column=2, row=1)
  tkgrid(ffButton, column=1, row=1, pady=c(4, 2), padx=2)
  tkgrid(pzButton, column=2, row=1, pady=c(4, 2))
  tkgrid(zcButton, column=1, row=2, pady=2, padx=2)
  tkgrid(zpButton, column=2, row=2, pady=2)
  tkgrid(zzButton, column=1, row=3, pady=2, padx=2)
  tkgrid(locButton, column=2, row=3, pady=2)
  tkbind(dlg, '<Return>', function(...) tryCatch(tkinvoke(tkfocus()), 
                                                 error=function(er){}))
  
  ## Now that the window is set up, deiconify and show the window
  tkwm.deiconify(dlg) # New
  tcl("update")  # New, Ensure the window is fully drawn and updated
	 
  invisible()
}

## Internal function openStored
## Add an entry to fileFolder for a spectrum stored in the global environment
## inFolder - list; data and file parameters for input spectrum.  This should 
##	match the	output format of ucsf2D()
## fileName - character string; name for the new entry in fileFolder
openStored <- function(inFolder, fileName='storedSpec'){
  
  ## Check input folder format
  if (is.null(inFolder$file.par))
    inFolder <- inFolder[[1]]
  if (is.null(inFolder$file.par) || is.null(inFolder$w2) || 
      is.null(inFolder$data))
    err('Object must contain a spectrum matching the format returned by ed().')
  if (inFolder$file.par$number_dimensions == 2 && is.null(inFolder$w1))
    err('Object must contain a spectrum matching the format returned by ed().')
  
  ## Look for fileName in fileFolder
  while(fileName %in% names(fileFolder)){
    fileName <- myDialog(paste('Name already exists in the file folder.', 
                               'Please provide a unique name:', sep='\n'), fileName)
    if (!length(fileName) || !nzchar(fileName))
      return(invisible())
  }
  
  ## Modify file parameters
  inFolder$file.par$user_title <- fileName
  inFolder$file.par$file.name <- fileName
  if (is.null(inFolder$graphics.par))
    inFolder$graphics.par <- defaultSettings
  
  ## Add entry to fileFolder
  n <- length(fileFolder) + 1
  fileFolder[[n]] <- inFolder
  names(fileFolder)[n] <- fileName
  
  ## Save changes and plot spectrum
  currentSpectrum <- fileName
  myAssign('fileFolder', fileFolder, FALSE)
  myAssign('currentSpectrum', currentSpectrum, FALSE)
  zf()
  
  return(fileName)
}

###################################################################################################
# Old fs function
###################################################################################################

## Displays an interactive GUI for sorting files
#' Open and Manage Spectra Files using fs
#'
#' This function provides a graphical user interface (GUI) for opening and managing spectra files.
#' Users can view information about open files, rearrange their order, open new files, and close
#' existing ones.
#'
#' @return NULL (The function primarily works with GUI elements.)
#'
#' @import tcltk2
#'
# ' @export
fs <- function(){
  
  ##call fo() if there are no open files
  if (!exists("fileFolder") || is.null(fileFolder) || !exists("currentSpectrum") 
      || is.null(currentSpectrum)){
    usrSel <- fo()
    if (is.null(usrSel))
      return(invisible())
  }
  ##creates main window
  tclCheck()
  
  dlg <- myToplevel('fs')
  if (is.null(dlg))
  {
    return(invisible())
  }
  tkwm.title(dlg, 'Files')
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  
  #	log_message('d')
  flush.console()
  
  ##create tablelist widget
  tableFrame <- ttklabelframe(dlg, 
                              text='Double-click on a file path to switch spectra:')
  
  xscr <- ttkscrollbar(tableFrame, orient='horizontal', command=function(...) 
    tkxview(tableList, ...))
  
  yscr <- ttkscrollbar(tableFrame, orient='vertical', command=function(...) 
    tkyview(tableList, ...))
  
  tableList <- tkwidget(tableFrame, 'tablelist::tablelist', bg='white',  
                        columns=c('0', 'Spectrum', '0', 'File Path', '0', 'Size', 'right', '0', 'Date Modified', 'center'), 
                        height=11, width=110, labelcommand=function(...) onSort(...), selectmode='extended', spacing=3, stretch='all', 
                        activestyle='underline', exportselection=FALSE, editselectedonly=TRUE, xscrollcommand=function(...) tkset(xscr, ...),
                        yscrollcommand=function(...) tkset(yscr, ...))
  
  #	log_message('e')
  flush.console()
  
  for (i in 0:3)
    tcl(tableList, 'columnconfigure', i, sortmode='dictionary')
  tcl(tableList, 'columnconfigure', 0, editable=TRUE)
  
  #	log_message('f')
  flush.console()	
  
  ##format the size column
  formatSize <- function(size)
  {
    if (length(grep('KB', size)))
      return(tclVar(size))
    size <- as.numeric(size)
    if (size < 2^20)
      return(tclVar(paste(signif(size / 2^10, 3), 'KB')))
    else
      return(tclVar(paste(signif(size / 2^20, 3), 'MB')))
  }
  tcl(tableList, 'columnconfigure', 2, formatcommand=function(...) 
    formatSize(...))
  
  #	log_message('g')
  flush.console()	
  
  ##selects all rows Ctrl+A is pressed
  tkbind(dlg, '<Control-a>', function(...) 
    tkselection.set(tableList, 0, 'end'))
  
  ##get file information and add to tablelist
  getFileInfo <- function(fileNames){
    if (!length(fileNames) || !nzchar(fileNames))
      return(invisible())
    tmpFolder <- fileFolder[fileNames]
    userTitles <- sapply(tmpFolder, function(x) x$file.par$user_title)
    
    sizes <- as.character(sapply(tmpFolder, function(x) 
      x$file.par$file.size))
    
    mods <- sapply(tmpFolder, function(x) 
      as.character(x$file.par$date.modified))
    paths <- sapply(tmpFolder, function(x) x$file.par$file.name)
    #		fileData <- cbind(userTitles, dims, paths, sizes, mods, deparse.level=0)
    fileData <- cbind(userTitles, paths, sizes, mods, deparse.level=0)
    for (i in 1:nrow(fileData))
      tkinsert(tableList, 'end', fileData[i, ])
  }
  getFileInfo(names(fileFolder))
  if (!is.null(currentSpectrum)){
    tkselection.set(tableList, wc() - 1)
    tcl(tableList, 'see', wc() - 1)
  }
  
  #	log_message('h')
  flush.console()	
  
  ##switches spectra on left-mouse double-click
  onDouble <- function(W){
    if (W != '.fs.1.3.body')
      return(invisible())
    usrSel <- as.numeric(tcl(tableList, 'curselection')) + 1
    usrFile <- names(fileFolder)[usrSel]
    if (!is.null(usrFile) && !is.na(usrFile) && currentSpectrum != usrFile){
      currentSpectrum <- usrFile
      myAssign('currentSpectrum', currentSpectrum)
      refresh(multi.plot=FALSE)
      tkwm.deiconify(dlg)
      tkfocus(tableList)
    }
  }
  tkbind(dlg, '<Double-Button-1>', onDouble)
  
  ##updates fileFolder
  updateFileFolder <- function(newOrder){
    newFolder <- fileFolder[newOrder]
    myAssign('fileFolder', newFolder)
  }
  
  ##sort files when user clicks on column headers 
  onSort <- function(tbl, col){
    prevOrder <- as.character(tcl(tableList, 'getcolumns', 0))
    tcl('tablelist::sortByColumn', tbl, col)
    newOrder <- as.character(tcl(tableList, 'getcolumns', 0))
    updateFileFolder(match(newOrder, prevOrder))
  }
  
  ##allow spectrum names to be edited
  onEdit <- function(widget, rowNum, colNum, newVal){
    rowNum <- as.numeric(rowNum) + 1
    userTitles <- sapply(fileFolder, function(x) x$file.par$user_title)
    if (newVal %in% userTitles){
      myMsg(paste('A spectrum with that name is currently open.', 
                  'Please enter a unique name.', sep='\n'), icon='error', 
            parent=dlg)
      tcl(tableList, 'cancelediting')
    }else{
      fileFolder[[rowNum]]$file.par$user_title <- newVal
      myAssign('fileFolder', fileFolder)
      refresh()
    }
    
    return(tclVar(as.character(newVal)))
  }
  tkconfigure(tableList, editendcommand=function(...) onEdit(...))
  
  ##create top button
  optionFrame <- ttkframe(dlg)
  onTop <- function(){
    usrSel <- as.numeric(tcl(tableList, 'curselection'))
    if (!length(usrSel) || usrSel == 0)
      return(invisible())
    prevOrder <- as.character(tcl(tableList, 'getcolumns', 0))
    for (i in seq_along(usrSel))
      tkmove(tableList, usrSel[i], i - 1)
    tcl(tableList, 'see', 0)
    newOrder <- as.character(tcl(tableList, 'getcolumns', 0))
    updateFileFolder(match(newOrder, prevOrder))
  }
  topButton <- ttkbutton(optionFrame, text='Top', width=8, command=onTop)
  
  ##create up button
  onUp <- function(){
    usrSel <- as.numeric(tcl(tableList, 'curselection'))
    if (!length(usrSel) || usrSel == 0)
      return(invisible())
    prevOrder <- as.character(tcl(tableList, 'getcolumns', 0))
    for (selItem in usrSel)
      tkmove(tableList, selItem, selItem - 1)
    tcl(tableList, 'see', min(usrSel) - 1)
    newOrder <- as.character(tcl(tableList, 'getcolumns', 0))
    updateFileFolder(match(newOrder, prevOrder))
  }
  upButton <- ttkbutton(optionFrame, text='^', width=8, command=onUp)
  
  ##create down button
  onDown <- function(){
    usrSel <- rev(as.numeric(tcl(tableList, 'curselection')))
    if (!length(usrSel) || usrSel == length(fileFolder) - 1)
      return(invisible())
    prevOrder <- as.character(tcl(tableList, 'getcolumns', 0))
    for (selItem in usrSel)
      tkmove(tableList, selItem, selItem + 2)
    tcl(tableList, 'see', max(usrSel) + 1)
    newOrder <- as.character(tcl(tableList, 'getcolumns', 0))
    updateFileFolder(match(newOrder, prevOrder))
  }
  downButton <- ttkbutton(optionFrame, text='v', width=8, command=onDown)
  
  ##create bottom button
  onBottom <- function(){
    usrSel <- rev(as.numeric(tcl(tableList, 'curselection')))
    if (!length(usrSel) || usrSel == length(fileFolder) - 1)
      return(invisible())
    prevOrder <- as.character(tcl(tableList, 'getcolumns', 0))
    for (i in seq_along(usrSel))
      tkmove(tableList, usrSel[i], length(fileFolder) - i + 1)
    tcl(tableList, 'see', length(fileFolder) - 1)
    newOrder <- as.character(tcl(tableList, 'getcolumns', 0))
    updateFileFolder(match(newOrder, prevOrder))
  }
  bottomButton <- ttkbutton(optionFrame, text='Bottom', width=8, 
                            command=onBottom)
  
  ##create file open button
  onOpen <- function(){
    
    ##open files
    prevFiles <- NULL
    if (length(fileFolder))
      prevFiles <- names(fileFolder)
    newFiles <- fo()
    if (is.null(newFiles))
      return(invisible())
    fileMatches <- as.vector(na.omit(match(prevFiles, newFiles)))
    if (length(fileMatches))
      newFiles <- newFiles[-fileMatches]
    if (!length(newFiles))
      return(invisible())
    
    ##get file information for newly opened files
    getFileInfo(newFiles)
    
    ##select newly opened files
    tkselection.clear(tableList, 0, 'end')
    tkselection.set(tableList, length(fileFolder) - length(newFiles), 'end')
    tkyview.moveto(tableList, 1)
  }
  openButton <- ttkbutton(optionFrame, text='Open file', width=11, 
                          command=onOpen)
  
  ##create file close button
  onClose <- function(){
    
    ##get user selection
    usrSel <- as.numeric(tcl(tableList, 'curselection'))
    if (!length(usrSel))
      return(invisible())
    usrFiles <- names(fileFolder)[usrSel + 1]
    
    ##close selected files
    fc(usrFiles)
    tkdelete(tableList, usrSel)
    tkselection.set(tableList, usrSel[1])
  }
  closeButton <- ttkbutton(optionFrame, text='Close file', width=11, 
                           command=onClose)
  
  ##create openStored label
  #	openStoredLabel <- ttklabel(dlg, text=paste('*Please type "?fs" in', 
  #					'the R console for more details on this feature.'))
  
  ##add widgets to tableFrame
  tkgrid(tableFrame, column=1, row=1, sticky='nswe', pady=6, padx=6)
  tkgrid(tableList, column=1, row=1, sticky='nswe')
  tkgrid(xscr, column=1, row=2, sticky='we')
  tkgrid(yscr, column=2, row=1, sticky='ns')
  
  ##make tableFrame stretch when window is resized
  tkgrid.columnconfigure(dlg, 1, weight=1)
  tkgrid.rowconfigure(dlg, 1, weight=1)
  tkgrid.columnconfigure(tableFrame, 1, weight=1)
  tkgrid.rowconfigure(tableFrame, 1, weight=1)
  
  ##add widgets to optionFrame
  tkgrid(optionFrame, column=1, row=2, pady=c(6, 0))
  tkgrid(topButton, column=1, row=1, padx=c(0, 4), pady=2)
  tkgrid(upButton, column=2, row=1, pady=2)
  tkgrid(downButton, column=3, row=1, padx=c(0, 4), pady=2)
  tkgrid(bottomButton, column=4, row=1, pady=2)
  tkgrid(openButton, column=5, row=1, padx=c(30, 4), pady=2)
  #	tkgrid(openStoredButton, column=6, row=1, padx=4, pady=2)
  tkgrid(closeButton, column=7, row=1, pady=2)
  #	tkgrid(openStoredLabel, column=1, row=3, padx=c(10, 0), pady=c(10, 0), sticky='w')
  tkgrid(ttksizegrip(dlg), column=2, row=3, sticky='se')
  
  ##allows users to press the 'Enter' key to make selections
  onEnter <- function(){
    focus <- as.character(tkfocus())
    if (focus == '.fs.1.3.body')
      onDouble('.fs.1.3.body')
    else
      tryCatch(tkinvoke(focus), error=function(er){})
  }
  tkbind(dlg, '<Return>', onEnter)
  
  ##updates widgets whenever the mouse enters the GUI
  onMouse <- function()
  {
    if (!length(fileFolder))
    {
      tkdelete(tableList, '0', 'end')
      return(invisible())
    }
    filePaths <- as.character(tcl(tableList, 'getcolumns', 1))
    folderPaths <- as.vector(sapply(fileFolder, function(x) x$file.par$file.name))
    userTitles <- as.character(tcl(tableList, 'getcolumns', 0))
    folderTitles <- getTitles(names(fileFolder), FALSE)		
    if (identical(folderPaths, filePaths) && identical(folderTitles, userTitles))
      return(invisible())
    
    tkdelete(tableList, '0', 'end')
    getFileInfo(names(fileFolder))
  }
  tkbind(dlg, '<Enter>', onMouse)
  tkbind(dlg, '<FocusIn>', onMouse)
  
  invisible()
}
#########################################################################################
# New fs function:
# fs <- function() {
#   # Check if Tcl/Tk is available
#   if (!inherits(tclRequire("Tcl"), "tclObj")) {
#     stop("Tcl/Tk is not available.")
#   }
  
#   # Ensure the tablelist package is available
#   tclRequire("tablelist")
  
#   # Call fo() if there are no open files
#   if (!exists("fileFolder") || is.null(fileFolder) || !exists("currentSpectrum") || is.null(currentSpectrum)) {
#     usrSel <- fo()
#     if (is.null(usrSel)) return(invisible())
#   }
  
#   # Create the main window
#   dlg <- tktoplevel()        # Use tktoplevel() directly
#   tkwm.title(dlg, 'Files')   # Set window title
#   tkfocus(dlg)
#   tkwm.deiconify(dlg)
  
#   flush.console()
  
#   # Create the tablelist widget
#   tableFrame <- ttklabelframe(dlg, text = 'Double-click on a file path to switch spectra:')
  
#   xscr <- ttkscrollbar(tableFrame, orient = 'horizontal', command = function(...) tkxview(tableList, ...))
#   yscr <- ttkscrollbar(tableFrame, orient = 'vertical', command = function(...) tkyview(tableList, ...))
  
#   tableList <- tkwidget(tableFrame, 'tablelist::tablelist', 
#                         bg = 'white',  
#                         columns = c('0', 'Spectrum', '0', 'File Path', '0', 'Size', 'right', '0', 'Date Modified', 'center'), 
#                         height = 11, width = 110, labelcommand = function(...) onSort(...), selectmode = 'extended', 
#                         spacing = 3, stretch = 'all', activestyle = 'underline', exportselection = FALSE, 
#                         editselectedonly = TRUE, xscrollcommand = function(...) tkset(xscr, ...),
#                         yscrollcommand = function(...) tkset(yscr, ...))
  
#   flush.console()
  
#   # Make columns sortable
#   for (i in 0:3) {
#     tcl(tableList, 'columnconfigure', i, sortmode = 'dictionary')
#   }
#   tcl(tableList, 'columnconfigure', 0, editable = TRUE)
  
#   # Format the size column
#   formatSize <- function(size) {
#     if (grepl('KB', size)) return(tclVar(size))
#     size <- as.numeric(size)
#     if (size < 2^20) {
#       return(tclVar(paste(signif(size / 2^10, 3), 'KB')))
#     } else {
#       return(tclVar(paste(signif(size / 2^20, 3), 'MB')))
#     }
#   }
#   tcl(tableList, 'columnconfigure', 2, formatcommand = function(...) formatSize(...))
  
#   flush.console()
  
#   # Select all rows when Ctrl+A is pressed
#   tkbind(dlg, '<Control-a>', function(...) tkselection.set(tableList, 0, 'end'))
  
#   # Debugging: Add test entry to ensure table is populated
#   tkinsert(tableList, 'end', c("Test Spectrum", "Test Path", "123KB", "2024-09-11"))
  
#   # Function to populate table with file info
#   getFileInfo <- function(fileNames) {
#     if (length(fileNames) == 0 || all(nzchar(fileNames) == FALSE)) return(invisible())
    
#     tmpFolder <- fileFolder[fileNames]
#     userTitles <- sapply(tmpFolder, function(x) x$file.par$user_title)
#     sizes <- as.character(sapply(tmpFolder, function(x) x$file.par$file.size))
#     mods <- sapply(tmpFolder, function(x) as.character(x$file.par$date.modified))
#     paths <- sapply(tmpFolder, function(x) x$file.par$file.name)
    
#     fileData <- cbind(userTitles, paths, sizes, mods, deparse.level = 0)
#     for (i in 1:nrow(fileData)) {
#       tkinsert(tableList, 'end', fileData[i, ])
#     }
#   }
  
#   # Populate table with current file information
#   getFileInfo(names(fileFolder))
  
#   # Make the current spectrum selected
#   if (!is.null(currentSpectrum)) {
#     tkselection.set(tableList, wc() - 1)
#     tcl(tableList, 'see', wc() - 1)
#   }
  
#   flush.console()
  
#   # Switch spectra on left mouse double-click
#   # onDouble <- function(W) {
#   #   #if (W != '.fs.1.3.body') return(invisible())
    
#   #   usrSel <- as.numeric(tcl(tableList, 'curselection')) + 1
#   #   usrFile <- names(fileFolder)[usrSel]
    
#   #   if (!is.null(usrFile) && !is.na(usrFile) && currentSpectrum != usrFile) {
#   #     currentSpectrum <- usrFile
#   #     myAssign('currentSpectrum', currentSpectrum)
#   #     refresh(multi.plot = FALSE)
#   #     tkwm.deiconify(dlg)
#   #     tkfocus(tableList)
#   #   }
#   # }
#   # tkbind(dlg, '<Double-Button-1>', onDouble)
		    
#   onDouble <- function(box) {
#   # Get the selected index from the tableList widget
#   usrSel <- 1 + as.integer(tkcurselection(box))
  
#   # Determine which file is selected based on the box ID
#   if (length(usrSel)) {
#     if (box$ID == '.fs.1.1.1.1') {
#       usrFile <- olFileNames[usrSel]
#     } else {
#       usrFile <- names(fileFolder)[usrSel]  # Use fileFolder in this case
#     }
#   } else {
#     usrFile <- NULL
#   }
  
#   # If a valid file is selected and it's different from the current spectrum
#   if (!is.null(usrFile) && currentSpectrum != usrFile) {
#     # Update the current spectrum
#     myAssign('currentSpectrum', usrFile)
    
#     # Refresh the plot window to display the selected file's spectrum
#     refresh(multi.plot = FALSE)
    
#     # Additional GUI configuration (if needed)
#     olConfigGui()
    
#     # Ensure the plot window is brought to the foreground
#     tkwm.deiconify(dlg)
#     tkfocus(box)
#   }
# }

#   # Add widgets to tableFrame
#   tkgrid(tableFrame, column = 1, row = 1, sticky = 'nswe', pady = 6, padx = 6)
#   tkgrid(tableList, column = 1, row = 1, sticky = 'nswe')
#   tkgrid(xscr, column = 1, row = 2, sticky = 'we')
#   tkgrid(yscr, column = 2, row = 1, sticky = 'ns')
  
#   # Make tableFrame stretch when the window is resized
#   tkgrid.columnconfigure(dlg, 1, weight = 1)
#   tkgrid.rowconfigure(dlg, 1, weight = 1)
#   tkgrid.columnconfigure(tableFrame, 1, weight = 1)
#   tkgrid.rowconfigure(tableFrame, 1, weight = 1)
  
#   ##updates fileFolder
#   updateFileFolder <- function(newOrder){
#     newFolder <- fileFolder[newOrder]
#     myAssign('fileFolder', newFolder)
#   }
  
#   ##sort files when user clicks on column headers 
#   onSort <- function(tbl, col){
#     prevOrder <- as.character(tcl(tableList, 'getcolumns', 0))
#     tcl('tablelist::sortByColumn', tbl, col)
#     newOrder <- as.character(tcl(tableList, 'getcolumns', 0))
#     updateFileFolder(match(newOrder, prevOrder))
#   }
  
#   ##allow spectrum names to be edited
#   onEdit <- function(widget, rowNum, colNum, newVal){
#     rowNum <- as.numeric(rowNum) + 1
#     userTitles <- sapply(fileFolder, function(x) x$file.par$user_title)
#     if (newVal %in% userTitles){
#       myMsg(paste('A spectrum with that name is currently open.', 
#                   'Please enter a unique name.', sep='\n'), icon='error', 
#             parent=dlg)
#       tcl(tableList, 'cancelediting')
#     }else{
#       fileFolder[[rowNum]]$file.par$user_title <- newVal
#       myAssign('fileFolder', fileFolder)
#       refresh()
#     }
    
#     return(tclVar(as.character(newVal)))
#   }
#   tkconfigure(tableList, editendcommand=function(...) onEdit(...))
  
#   # Create the top button
#   optionFrame <- ttkframe(dlg)
  
#   onTop <- function() {
#     usrSel <- as.numeric(tcl(tableList, 'curselection'))
#     if (!length(usrSel) || usrSel == 0) return(invisible())
    
#     prevOrder <- as.character(tcl(tableList, 'getcolumns', 0))
#     for (i in seq_along(usrSel)) {
#       tkmove(tableList, usrSel[i], i - 1)
#     }
#     tcl(tableList, 'see', 0)
#     newOrder <- as.character(tcl(tableList, 'getcolumns', 0))
#     updateFileFolder(match(newOrder, prevOrder))
#   }
#   topButton <- ttkbutton(optionFrame, text = 'Top', width = 8, command = onTop)
  
#   # Create the up button
#   onUp <- function() {
#     usrSel <- as.numeric(tcl(tableList, 'curselection'))
#     if (!length(usrSel) || usrSel == 0) return(invisible())
    
#     prevOrder <- as.character(tcl(tableList, 'getcolumns', 0))
#     for (selItem in usrSel) {
#       tkmove(tableList, selItem, selItem - 1)
#     }
#     tcl(tableList, 'see', min(usrSel) - 1)
#     newOrder <- as.character(tcl(tableList, 'getcolumns', 0))
#     updateFileFolder(match(newOrder, prevOrder))
#   }
#   upButton <- ttkbutton(optionFrame, text = '^', width = 8, command = onUp)
  
#   # Create the down button
#   onDown <- function() {
#     usrSel <- rev(as.numeric(tcl(tableList, 'curselection')))
#     if (!length(usrSel) || usrSel == length(fileFolder) - 1) return(invisible())
    
#     prevOrder <- as.character(tcl(tableList, 'getcolumns', 0))
#     for (selItem in usrSel) {
#       tkmove(tableList, selItem, selItem + 2)
#     }
#     tcl(tableList, 'see', max(usrSel) + 1)
#     newOrder <- as.character(tcl(tableList, 'getcolumns', 0))
#     updateFileFolder(match(newOrder, prevOrder))
#   }
#   downButton <- ttkbutton(optionFrame, text = 'v', width = 8, command = onDown)
  
#   # Create the bottom button
#   onBottom <- function() {
#     usrSel <- rev(as.numeric(tcl(tableList, 'curselection')))
#     if (!length(usrSel) || usrSel == length(fileFolder) - 1) return(invisible())
    
#     prevOrder <- as.character(tcl(tableList, 'getcolumns', 0))
#     for (i in seq_along(usrSel)) {
#       tkmove(tableList, usrSel[i], length(fileFolder) - i + 1)
#     }
#     tcl(tableList, 'see', length(fileFolder) - 1)
#     newOrder <- as.character(tcl(tableList, 'getcolumns', 0))
#     updateFileFolder(match(newOrder, prevOrder))
#   }
#   bottomButton <- ttkbutton(optionFrame, text = 'Bottom', width = 8, command = onBottom)
  
#   ## Add option buttons to the optionFrame
#   tkgrid(optionFrame, column = 1, row = 2, pady = c(6, 0))
#   tkgrid(topButton, column = 1, row = 1, padx = c(0, 4), pady = 2)
#   tkgrid(upButton, column = 2, row = 1, pady = 2)
#   tkgrid(downButton, column = 3, row = 1, padx = c(0, 4), pady = 2)
#   tkgrid(bottomButton, column = 4, row = 1, pady = 2)
   
#  ##create file open button
#   onOpen <- function(){
    
#     ##open files
#     prevFiles <- NULL
#     if (length(fileFolder))
#       prevFiles <- names(fileFolder)
#     newFiles <- fo()
#     if (is.null(newFiles))
#       return(invisible())
#     fileMatches <- as.vector(na.omit(match(prevFiles, newFiles)))
#     if (length(fileMatches))
#       newFiles <- newFiles[-fileMatches]
#     if (!length(newFiles))
#       return(invisible())
    
#     ##get file information for newly opened files
#     getFileInfo(newFiles)
    
#     ##select newly opened files
#     tkselection.clear(tableList, 0, 'end')
#     tkselection.set(tableList, length(fileFolder) - length(newFiles), 'end')
#     tkyview.moveto(tableList, 1)
#   }
#   openButton <- ttkbutton(optionFrame, text='Open file', width=11, 
#                           command=onOpen)
  
#   ##create file close button
#   onClose <- function(){
    
#     ##get user selection
#     usrSel <- as.numeric(tcl(tableList, 'curselection'))
#     if (!length(usrSel))
#       return(invisible())
#     usrFiles <- names(fileFolder)[usrSel + 1]
    
#     ##close selected files
#     fc(usrFiles)
#     tkdelete(tableList, usrSel)
#     tkselection.set(tableList, usrSel[1])
#   }
#   closeButton <- ttkbutton(optionFrame, text='Close file', width=11, 
#                            command=onClose)
  
#   ## Additional configuration for resizing
#   tkgrid(ttksizegrip(dlg), column = 2, row = 3, sticky = 'se')
  
#   ## Bind Enter key to select or invoke commands
#   onEnter <- function() {
#     focus <- as.character(tkfocus())
#     if (focus == '.fs.1.3.body') {
#       onDouble('.fs.1.3.body')
#     } else {
#       tryCatch(tkinvoke(focus), error = function(er) {})
#     }
#   }
#   tkbind(dlg, '<Return>', onEnter)
  
#   ## Updates the window when mouse enters the GUI
#   onMouse <- function() {
#     if (!length(fileFolder)) {
#       tkdelete(tableList, '0', 'end')
#       return(invisible())
#     }
#     filePaths <- as.character(tcl(tableList, 'getcolumns', 1))
#     folderPaths <- as.vector(sapply(fileFolder, function(x) x$file.par$file.name))
#     userTitles <- as.character(tcl(tableList, 'getcolumns', 0))
#     folderTitles <- getTitles(names(fileFolder), FALSE)
#     if (identical(folderPaths, filePaths) && identical(folderTitles, userTitles)) return(invisible())
    
#     tkdelete(tableList, '0', 'end')
#     getFileInfo(names(fileFolder))
#   }
#   tkbind(dlg, '<Enter>', onMouse)
#   tkbind(dlg, '<FocusIn>', onMouse)
  
#   invisible()
# }
#########################################################################################################

## Interactive GUI for editing tables
## data - data.frame; the table to edit
## editable - logical vector; indicates whether entries in each column in the 
##	table should be editable, should be the same length as the number of columns
## title - character string; title for the GUI
## colVer - function list; functions, one for each column, used to verify the 
##	entries in each column.  Functions should return TRUE or FALSE.
## errMsgs - character vector; error messages to display if a function in colVer
##	returns FALSE, should be the same length as as the number of columns.  If 
##	NULL, no error messages are displayed
tableEdit <- function(data, editable=rep(TRUE, ncol(data)),	title='digestR', 
                      colVer=NULL, errMsgs=rep(paste('Data type for new entry must',  
                                                     'match the data type for a given column.'), ncol(data))){
  
  ##check colVer argument
  if (is.null(colVer)){
    verFun <- function(x) return(TRUE)
    colVer <- rep(list(verFun), ncol(data))
  }
  if (length(colVer) != ncol(data))
    stop('length of colVer argument must equal the number of columns in data')
  
  ##check errMsgs argument
  if (!is.null(errMsgs) && length(errMsgs) != ncol(data))
    stop('length of errMsgs argument must equal the number of columns in data')
  
  ##creates main window
  dlg <- tktoplevel()
  tcl('wm', 'attributes', dlg, topmost=TRUE)
  returnVal <- NULL
  tkwm.title(dlg, title)
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  colNames <- colnames(data)
  
  ##create tablelist widget
  tableFrame <- ttklabelframe(dlg, text='Data Table:')
  xscr <- ttkscrollbar(tableFrame, orient='horizontal', command=function(...) 
    tkxview(tableList, ...))
  yscr <- ttkscrollbar(tableFrame, orient='vertical', command=function(...) 
    tkyview(tableList, ...))
  colVals <- NULL
  for (i in colNames)
    colVals <- c(colVals, '0', i, 'center')
  tableList <- tkwidget(tableFrame, 'tablelist::tablelist', columns=colVals, 
                        activestyle='underline', height=11, width=110, exportselection=FALSE,
                        labelcommand='tablelist::sortByColumn', selectmode='extended', bg='white', 
                        spacing=3, stretch='all', editselectedonly=TRUE, selecttype='cell',
                        xscrollcommand=function(...) tkset(xscr, ...),
                        yscrollcommand=function(...) tkset(yscr, ...))
  
  ##add data to tablelist widget
  for (i in 1:nrow(data))
    tkinsert(tableList, 'end', unlist(data[i, ]))
  for (i in 1:ncol(data))
    tcl(tableList, 'columnconfigure', i - 1, sortmode='dictionary', 
        editable=editable[i], width=0, align='left')
  
  ##get the data types for each column
  colTypes <- NULL
  for (i in 1:ncol(data))
    colTypes <- c(colTypes, storage.mode(data[, i]))
  
  ##selects all rows Ctrl+A is pressed
  tkbind(tableList, '<Control-a>', function(...) 
    tkselect(tableList, 'set', 0, 'end'))
  
  ##rewrites data after GUI is updated
  writeData <- function(){
    
    ##get the data from the GUI
    newData <- NULL
    numRows <- as.numeric(tcl(tableList, 'index', 'end'))
    if (numRows == 0){
      data <<- newData
      return(invisible())
    }
    for (i in 0:numRows)
      newData <- rbind(newData, as.character(tcl(tableList, 'get', i)))
    
    ##format data
    colnames(newData) <- colNames
    newData <- as.data.frame(newData, stringsAsFactors=FALSE)
    for (i in 1:ncol(newData))
      suppressWarnings(storage.mode(newData[, i]) <- colTypes[i])
    data <<- newData
  }
  
  ##save tableList data after table is sorted
  tkbind(dlg, '<<TablelistColumnSorted>>', writeData)
  
  ##create top button
  optionFrame <- ttkframe(tableFrame)
  moveFrame <- ttklabelframe(optionFrame, text='Move selected rows', padding=6)
  onTop <- function(){
    usrSel <- as.numeric(tcl(tableList, 'curselection'))
    if (!length(usrSel) || usrSel == 0)
      return(invisible())
    tkselection.set(tableList, usrSel)
    for (i in seq_along(usrSel))
      tkmove(tableList, usrSel[i], 0 + i - 1)
    tcl(tableList, 'see', 0)
    writeData()
  }
  topButton <- ttkbutton(moveFrame, text='Top', width=11, command=onTop)
  
  ##create up button
  onUp <- function(){
    usrSel <- as.numeric(tcl(tableList, 'curselection'))
    if (!length(usrSel) || usrSel == 0)
      return(invisible())
    tkselection.set(tableList, usrSel)
    for (selItem in usrSel)
      tkmove(tableList, selItem, selItem - 1)
    tcl(tableList, 'see', min(usrSel) - 1)
    writeData()
  }
  upButton <- ttkbutton(moveFrame, text='^', width=9, command=onUp)
  
  ##create down button
  onDown <- function(){
    usrSel <- rev(as.numeric(tcl(tableList, 'curselection')))
    if (!length(usrSel) || usrSel == nrow(data) - 1)
      return(invisible())
    tkselection.set(tableList, usrSel)
    for (selItem in usrSel)
      tkmove(tableList, selItem, selItem + 2)
    tcl(tableList, 'see', max(usrSel) + 1)
    writeData()
  }
  downButton <- ttkbutton(moveFrame, text='v', width=9, command=onDown)
  
  ##create bottom button
  onBottom <- function(){
    usrSel <- rev(as.numeric(tcl(tableList, 'curselection')))
    if (!length(usrSel) || usrSel == nrow(data) - 1)
      return(invisible())
    tkselection.set(tableList, usrSel)
    for (i in seq_along(usrSel))
      tkmove(tableList, usrSel[i], nrow(data) - i + 1)
    tcl(tableList, 'see', nrow(data) - 1)
    writeData()
  }
  bottomButton <- ttkbutton(moveFrame, text='Bottom', width=11, 
                            command=onBottom)
  
  ##create sig. fig. spinbox
  sigFigFrame <- ttklabelframe(optionFrame, text='Display', padding=6)
  onSigFig <- function(){
    if (tclvalue(sigFigVal) == 'max'){
      for (i in 1:nrow(data))
        tcl(tableList, 'rowconfigure', i - 1, text=unlist(data[i, ]))
      return(invisible())
    }
    sigFig <- as.numeric(tclvalue(sigFigVal))
    for (i in seq_along(data[1, ])){
      if (any(is.logical(data[, i])))
        next
      newData <- tryCatch(signif(data[, i], sigFig), 
                          error=function(er) return(data[, i]))
      newData[is.na(newData)] <- 'NA'
      tcl(tableList, 'columnconfigure', i - 1, text=newData)
    }
  }
  sigFigVal <- tclVar('max')
  sigFigBox <- tkwidget(sigFigFrame, 'spinbox', width=6, wrap=TRUE,
                        textvariable=sigFigVal, values=c('max', 1:9), command=onSigFig)
  sigFigLabel <- ttklabel(sigFigFrame, text='significant figures')
  
  ##create table edit widgets
  if (any(editable)){
    
    ##check interactively edited cells using functions provided in colVer
    onEdit <- function(widget, rowNum, colNum, newVal, tclReturn=TRUE){
      rowNum <- as.numeric(rowNum) + 1
      colNum <- as.numeric(colNum) + 1
      if (newVal == 'NA'){
        if (tclReturn)
          return(tclVar(as.character(newVal)))
        else
          return(TRUE)
      }
      suppressWarnings(storage.mode(newVal) <- colTypes[colNum])
      if (!colVer[[colNum]](newVal)){
        if (!is.null(errMsgs))
          myMsg(errMsgs[colNum], icon='error', parent=dlg)
        if (tclReturn)
          tcl(tableList, 'cancelediting')
        else
          return(FALSE)
      }else
        data[rowNum, colNum] <<- newVal
      if (tclReturn)
        return(tclVar(as.character(newVal)))
      else
        return(TRUE)
    }
    tkconfigure(tableList, editendcommand=function(...) onEdit(...))
    
    ##create cell editing textbox
    ceditFrame <- ttklabelframe(optionFrame, text='Edit selected cells', 
                                padding=6)
    usrEntry <- tclVar(character(0))
    textEntry <- ttkentry(ceditFrame, width=13, justify='center', 
                          textvariable=usrEntry)
    
    ##update cell editing textbox with current cell selection value
    onCellSel <- function(){
      usrSel <- as.character(tcl(tableList, 'curcellselection'))
      if (!length(usrSel))
        tclObj(usrEntry) <- character(0)
      selVals <- as.character(tcl(tableList, 'getcells', usrSel))
      if (length(grep(selVals[1], selVals, fixed=TRUE)) == length(selVals))
        tclObj(usrEntry) <- selVals[1]
      else
        tclObj(usrEntry) <- character(0)
    }
    tkbind(tableList, '<<TablelistSelect>>', onCellSel)
    
    ##create apply button
    onApply <- function(){
      tcl(tableList, 'finishediting')
      newVal <- tclvalue(usrEntry)
      usrSel <- as.character(tcl(tableList, 'curcellselection'))
      for (i in usrSel){
        rowNum <- unlist(strsplit(i, ','))[1]
        colNum <- unlist(strsplit(i, ','))[2]
        isValid <- onEdit(rowNum=rowNum, colNum=colNum, newVal=newVal, 
                          tclReturn=FALSE)
        if (isValid)
          tcl(tableList, 'cellconfigure', i, text=newVal)
        else
          return(invisible())
      }
      writeData()
    }
    applyButton <- ttkbutton(ceditFrame, text='Apply', width=8, command=onApply)
    
    ##create copy button
    reditFrame <- ttklabelframe(optionFrame, text='Edit rows', padding=6)
    clipboard <- NULL
    onCopy <- function(){
      usrSel <- as.numeric(tcl(tableList, 'curselection'))
      if (!length(usrSel) || usrSel == 0)
        return(invisible())
      tkselection.set(tableList, usrSel)
      selVals <- NULL
      for (i in usrSel)
        selVals <- rbind(selVals, as.character(tcl(tableList, 'get', i)))
      clipboard <<- selVals
    }
    copyButton <- ttkbutton(reditFrame, text='Copy', width=10, command=onCopy)
    
    ##create paste button
    onPaste <- function(){
      if (is.null(clipboard))
        return(invisible())
      for (i in 1:nrow(clipboard))
        tkinsert(tableList, 'end', unlist(clipboard[i, ]))
      writeData()
      tcl(tableList, 'see', nrow(data) - 1)
    }
    pasteButton <- ttkbutton(reditFrame, text='Paste', width=10, 
                             command=onPaste)
    
    ##create insert button
    onInsert <- function(){
      tkinsert(tableList, 'end', as.character(rep(NA, ncol(data))))
      writeData()
      tcl(tableList, 'see', nrow(data) - 1)
    }
    insertButton <- ttkbutton(reditFrame, text='Insert', width=10, 
                              command=onInsert)
    
    ##create delete button
    onDelete <- function(){
      usrSel <- as.numeric(tcl(tableList, 'curselection'))
      if (!length(usrSel))
        return(invisible())
      tkselection.set(tableList, usrSel - 1)
      tkdelete(tableList, usrSel)
      writeData()
    }
    deleteButton <- ttkbutton(reditFrame, text='Delete', width=10, 
                              command=onDelete)
  }
  
  ##create ok button
  bottomFrame <- ttkframe(dlg)
  onOk <- function(){
    
    ##verify data
    tcl(tableList, 'finishediting')
    if (!is.null(data)){
      for (i in 1:ncol(data)){
        suppressWarnings(storage.mode(data[, i]) <- colTypes[i])
        if (!colVer[[i]](data[, i])){
          if (!is.null(errMsgs))
            myMsg(errMsgs[i], icon='error', parent=dlg)
          return(invisible())
        }
      }
    }
    
    ##return the data and close the GUI
    returnVal <<- data
    tkgrab.release(dlg)
    tkdestroy(dlg)
    return(returnVal)
  }
  okButton <- ttkbutton(bottomFrame, text='OK', width=10, command=onOk)
  
  ##create cancel button
  onCancel <- function(){
    tkgrab.release(dlg)
    tkdestroy(dlg)
    return(returnVal)
  }
  cancelButton <- ttkbutton(bottomFrame, text='Cancel', width=10, 
                            command=onCancel)
  
  ##create export button
  onExport <- function(){
    tkwm.iconify(dlg)
    if ('ACTIVE' %in% names(data))
      initFile <- 'roiTable'
    else if ('Index' %in% names(data))
      initFile <- 'peakList'
    else
      initFile <- 'roiSummary'
    fileName <- mySave(initialfile=initFile, defaultextension='txt', 
                       title='Export', filetypes=list('xls'='Excel Files', 'txt'='Text Files'))
    if (length(fileName) == 0 || !nzchar(fileName)){
      tkwm.deiconify(dlg)
      return(invisible())
    }
    write.table(data, file=fileName, quote=FALSE, sep='\t', row.names=FALSE, 
                col.names=TRUE)
    tkwm.deiconify(dlg)
  }
  exportButton <- ttkbutton(bottomFrame, text='Export', width=10, 
                            command=onExport)
  
  ##add widgets to treeFrame
  tkgrid(tableFrame, column=1, row=1, sticky='nswe', pady=6, padx=6)
  tkgrid(tableList, column=1, row=1, sticky='nswe')
  tkgrid(xscr, column=1, row=2, sticky='we')
  tkgrid(yscr, column=2, row=1, sticky='ns')
  
  ##make treeFrame stretch when window is resized
  tkgrid.columnconfigure(dlg, 1, weight=1)
  tkgrid.rowconfigure(dlg, 1, weight=1)
  tkgrid.columnconfigure(tableFrame, 1, weight=1)
  tkgrid.rowconfigure(tableFrame, 1, weight=1)
  
  ##add widgets to moveFrame
  tkgrid(optionFrame, column=1, columnspan=2, row=3, pady=8)
  tkgrid(moveFrame, column=1, row=1, padx=8)
  tkgrid(topButton, column=1, row=1, pady=2, padx=c(0, 4))
  tkgrid(upButton, column=2, row=1, pady=2, padx=1)
  tkgrid(downButton, column=3, row=1, padx=1, pady=2)
  tkgrid(bottomButton, column=4, row=1, pady=2, padx=c(4, 0))
  
  ##add widgets to sigFigFrame
  tkgrid(sigFigFrame, column=2, row=1, padx=8)
  tkgrid(sigFigBox, column=1, row=1, padx=c(4, 2), pady=c(2, 4))
  tkgrid(sigFigLabel, column=2, row=1, padx=c(0, 4), pady=c(2, 4))
  
  ##add editing widgets
  if (any(editable)){
    
    ##add widgets to rowFrame
    tkgrid(reditFrame, column=1, row=2, pady=4, padx=8)
    tkgrid(copyButton, column=1, row=1, padx=c(0, 2))
    tkgrid(pasteButton, column=2, row=1, padx=c(0, 8))
    tkgrid(insertButton, column=3, row=1, padx=c(0, 2))
    tkgrid(deleteButton, column=4, row=1, padx=c(0, 0))
    
    ##add widgets to ceditFrame
    tkgrid(ceditFrame, column=2, row=2, pady=4, padx=8)
    tkgrid(textEntry, column=1, row=1, padx=2)
    tkgrid(applyButton, column=3, row=1, padx=2)
  }
  
  ##add widgets to bottomFrame
  tkgrid(bottomFrame, column=1, row=2, pady=c(6, 0))
  tkgrid(okButton, column=1, row=1, padx=4)
  tkgrid(cancelButton, column=2, row=1, padx=4)
  tkgrid(exportButton, column=3, row=1, padx=c(20, 4))
  tkgrid(ttksizegrip(dlg), column=1, row=3, sticky='se')
  
  ##Allows users to press the 'Enter' key to make selections
  onEnter <- function(){
    focus <- as.character(tkfocus())
    if (length(grep('.2.2.3.1$', focus)))
      onApply()
    else
      tryCatch(tkinvoke(focus), error=function(er){})
  }
  tkbind(dlg, '<Return>', onEnter)
  
  ## configure the toplevel
  tkfocus(tableList)
  if (as.logical(tkwinfo('viewable', dlg)))
    tkgrab.set(dlg)
  tkwait.window(dlg)
  return(returnVal)
  
  invisible()
}

## User preferences GUI
ep <- function(dispPane=0){
  
  ##create main window
  tclCheck()
  dlg <- myToplevel('ep')
  if (is.null(dlg))
    return(invisible())
  tkwm.title(dlg, 'Preferences')
  tkfocus(dlg)
  tkwm.deiconify(dlg)
  newDef <- defaultSettings
  
  ##create paned notebook
  epBook <- ttknotebook(dlg, padding=3)
  tkgrid(epBook, column=1, row=1, sticky='nsew', padx=6, pady=6)
  
  ##create individual panes
  genFrame <- ttkframe(epBook, padding=6) 
  tkadd(epBook, genFrame, text=' General ')	
  grFrame <- ttkframe(epBook, padding=6)
  tkadd(epBook, grFrame, text=' Graphics ')
  coFrame <- ttkframe(epBook, padding=6)
  tkadd(epBook, coFrame, text=' Colors ')
  psFrame <- ttkframe(epBook, padding=6)
  tkadd(epBook, psFrame, text=' Plotting ')
  ppFrame <- ttkframe(epBook, padding=6)
  tkadd(epBook, ppFrame, text=' Peak Picking ')
  roiFrame <- ttkframe(epBook, padding=6)
  tkadd(epBook, roiFrame, text=' ROIs ')
  tkselect(epBook, dispPane)
  
  #####create widgets for genFrame
  ##create a label with instructions
  genLabel <- ttklabel(genFrame, wraplength=350, text=paste('The following', 
                                                            'settings will be applied when the digestR package is loaded.  Press', 
                                                            'the "?" button for more information on each setting.'))
  
  ##create window dimension settings frame
  sizeFrame <- ttklabelframe(genFrame, text='Window dimensions', padding=6)
  widthLabel <- ttklabel(sizeFrame, text='width:')
  heightLabel <- ttklabel(sizeFrame, text='height:')
  
  ##create main plot window widgets
  mainFrame <- ttkframe(sizeFrame)
  mainLabel <- ttklabel(mainFrame, text='Main plot')
  mainWdVar <- tclVar(newDef$size.main[1])
  mainWdEntry <- ttkentry(mainFrame, textvariable=mainWdVar, width=6,
                          justify='center')
  mainHtVar <- tclVar(newDef$size.main[2])
  mainHtEntry <- ttkentry(mainFrame, textvariable=mainHtVar, width=6,
                          justify='center')
  
  ##create subplot window widgets
  subFrame <- ttkframe(sizeFrame)
  subLabel <- ttklabel(subFrame, text='Subplot')
  subWdVar <- tclVar(newDef$size.sub[1])
  subWdEntry <- ttkentry(subFrame, textvariable=subWdVar, width=6,
                         justify='center')
  subHtVar <- tclVar(newDef$size.sub[2])
  subHtEntry <- ttkentry(subFrame, textvariable=subHtVar, width=6,
                         justify='center')
  
  ##create multiple file window widgets
  multiFrame <- ttkframe(sizeFrame)
  multiLabel <- ttklabel(multiFrame, text='Multi file')
  multiWdVar <- tclVar(newDef$size.multi[1])
  multiWdEntry <- ttkentry(multiFrame, textvariable=multiWdVar, width=6,
                           justify='center')
  multiHtVar <- tclVar(newDef$size.multi[2])
  multiHtEntry <- ttkentry(multiFrame, textvariable=multiHtVar, width=6,
                           justify='center')
  
  ##create SDI checkbox
  checkFrame <- ttkframe(genFrame)
  sdiVar <- tclVar(as.character(newDef$sdi))
  sdiButton <- ttkcheckbutton(checkFrame, onvalue='TRUE', offvalue='FALSE', 
                              variable=sdiVar, text=paste(' Run digestR using\n separate windows'))
  
  ##create update checkbox
  updateVar <- tclVar(as.character(newDef$update))
  updateButton <- ttkcheckbutton(checkFrame, onvalue='TRUE', offvalue='FALSE', 
                                 variable=updateVar, text=' Check for updates\n when digestR loads')
  
  ##create auto backup checkbox
  backupVar <- tclVar(as.character(newDef$autoBackup))
  backupButton <- ttkcheckbutton(checkFrame, onvalue='TRUE', offvalue='FALSE', 
                                 variable=backupVar, text=paste(' Enable automatic\n backup'))
  
  ##create working directory selection widgets
  dirFrame <- ttkframe(genFrame)
  dirLabel <- ttklabel(dirFrame, text='Default directory:')
  dirVar <- tclVar(newDef$wd)
  dirEntry <- ttkentry(dirFrame, textvariable=dirVar, justify='left', width=38, 
                       state='readonly')
  onAddLib <- function(){
    newWd <- myDir(parent=dlg)
    if (nzchar(newWd))
      tclvalue(dirVar) <- newWd
  }
  browseButton <- ttkbutton(dirFrame, text='Browse', width=10, command=onAddLib)
  
  ##add widgets to genFrame
  tkgrid(genLabel, column=1, columnspan=2, row=1, padx=c(6, 10), pady=c(10, 0), 
         sticky='w')
  
  tkgrid(sizeFrame, column=1, row=2, padx=4, pady=8, sticky='nsw')
  tkgrid(widthLabel, column=1, row=2, pady=4, sticky='sw')
  tkgrid(heightLabel, column=1, row=3, pady=4, sticky='sw')
  
  tkgrid(mainFrame, column=2, row=1, rowspan=3, padx=6)
  tkgrid(mainLabel, column=1, row=1, pady=c(0, 2), sticky='w')
  tkgrid(mainWdEntry, column=1, row=2, pady=4, sticky='e')
  tkgrid(mainHtEntry, column=1, row=3, pady=4, sticky='e')
  
  tkgrid(subFrame, column=3, row=1, rowspan=3, padx=6)
  tkgrid(subLabel, column=1, row=1, pady=c(0, 2), sticky='w')
  tkgrid(subWdEntry, column=1, row=2, pady=4, sticky='e')
  tkgrid(subHtEntry, column=1, row=3, pady=4, sticky='e')
  
  tkgrid(multiFrame, column=4, row=1, rowspan=3, padx=6)
  tkgrid(multiLabel, column=1, row=1, pady=c(0, 2), sticky='w')
  tkgrid(multiWdEntry, column=1, row=2, pady=4, sticky='e')
  tkgrid(multiHtEntry, column=1, row=3, pady=4, sticky='e')
  
  tkgrid(checkFrame, column=2, row=2, padx=4, pady=8, sticky='ns')
  tkgrid(sdiButton, column=1, row=1, pady=c(3, 0), sticky='w')
  tkgrid(updateButton, column=1, row=2, pady=6, sticky='w')
  tkgrid(backupButton, column=1, row=3, sticky='w')
  
  tkgrid(dirFrame, column=1, columnspan=2, row=3, padx=4, pady=4, sticky='w')
  tkgrid(dirLabel, column=1, row=1, pady=c(0, 3), sticky='w')
  tkgrid(dirEntry, column=1, row=2, padx=c(12, 4), pady=2)
  tkgrid(browseButton, column=2, row=2, pady=2, sticky='e')
  
  #####create widgets for grFrame
  ##create plot margin options
  marFrame <- ttklabelframe(grFrame, text='Plot margins', padding=4)
  topLabel <- ttklabel(marFrame, text='top:')
  bottomLabel <- ttklabel(marFrame, text='bottom:')
  leftLabel <- ttklabel(marFrame, text='left:')
  rightLabel <- ttklabel(marFrame, text='right:')
  
  marMainLabel <- ttklabel(marFrame, text='Main plot')
  topMarVar <- tclVar(newDef$mar[3])
  topMarEntry <- ttkentry(marFrame, textvariable=topMarVar, width=8,
                          justify='center')
  botMarVar <- tclVar(newDef$mar[1])
  botMarEntry <- ttkentry(marFrame, textvariable=botMarVar, width=8,
                          justify='center')
  leftMarVar <- tclVar(newDef$mar[2])
  leftMarEntry <-	ttkentry(marFrame, textvariable=leftMarVar, width=8,
                           justify='center')
  rightMarVar <- tclVar(newDef$mar[4])
  rightMarEntry <- ttkentry(marFrame, textvariable=rightMarVar, width=8,
                            justify='center')
  
  marSubLabel <- ttklabel(marFrame, text='Subplot')
  topSubMarVar <- tclVar(newDef$mar.sub[3])
  topSubMarEntry <-	ttkentry(marFrame, textvariable=topSubMarVar, width=8,
                             justify='center')
  botSubMarVar <- tclVar(newDef$mar.sub[1])
  botSubMarEntry <-	ttkentry(marFrame, textvariable=botSubMarVar, width=8,
                             justify='center')
  leftSubMarVar <- tclVar(newDef$mar.sub[2])
  leftSubMarEntry <-	ttkentry(marFrame, textvariable=leftSubMarVar, width=8,
                              justify='center')
  rightSubMarVar <- tclVar(newDef$mar.sub[4])
  rightSubMarEntry <-	ttkentry(marFrame, textvariable=rightSubMarVar, width=8,
                               justify='center')
  
  marMultiLabel <- ttklabel(marFrame, text='Multi file')
  topMultiMarVar <- tclVar(newDef$mar.multi[3])
  topMultiMarEntry <-	ttkentry(marFrame, textvariable=topMultiMarVar, width=8,
                               justify='center')
  botMultiMarVar <- tclVar(newDef$mar.multi[1])
  botMultiMarEntry <-	ttkentry(marFrame, textvariable=botMultiMarVar, width=8,
                               justify='center')
  leftMultiMarVar <- tclVar(newDef$mar.multi[2])
  leftMultiMarEntry <-	ttkentry(marFrame, textvariable=leftMultiMarVar, 
                                width=8, justify='center')
  rightMultiMarVar <- tclVar(newDef$mar.multi[4])
  rightMultiMarEntry <-	ttkentry(marFrame, textvariable=rightMultiMarVar, 
                                 width=8, justify='center')
  
  ##create x-axis gridlines checkbox
  gridFrame <- ttklabelframe(grFrame, text='Gridlines', padding=2)
  if (is.na(newDef$xtck))
    xtckVar <- tclVar(0)
  else
    xtckVar <- tclVar(1)
  xgridButton <- ttkcheckbutton(gridFrame, variable=xtckVar, text='x-axis')
  
  ##create y-axis gridlines checkbox
  if (is.na(newDef$ytck))
    ytckVar <- tclVar(0)
  else
    ytckVar <- tclVar(1)
  ygridButton <- ttkcheckbutton(gridFrame, variable=ytckVar, text='y-axis')
  
  ##create text magnification options
  cexFrame <- ttklabelframe(grFrame, text='Text magnification', padding=4)
  cexMainLabel <- ttklabel(cexFrame, text='Main plot')
  titleLabel <- ttklabel(cexFrame, text='title:')
  cexMainVar <- tclVar(newDef$cex.main)
  titleEntry <-	ttkentry(cexFrame, textvariable=cexMainVar, width=8,
                         justify='center')
  axesLabel <- ttklabel(cexFrame, text='axes:')
  cexAxesVar <- tclVar(newDef$cex.axis)
  axesEntry <-	ttkentry(cexFrame, textvariable=cexAxesVar, width=8,
                        justify='center')
  
  cexMultiLabel <- ttklabel(cexFrame, text='Multi file')
  roiMultiLabel <- ttklabel(cexFrame, text='ROI names:')
  cexRoiMultiVar <- tclVar(newDef$cex.roi.multi)
  roiMultiEntry <-	ttkentry(cexFrame, textvariable=cexRoiMultiVar, width=8,
                            justify='center')
  fileLabel <- ttklabel(cexFrame, text='file names:')
  cexFilesMultiVar <- tclVar(newDef$cex.files.multi)
  fileEntry <-	ttkentry(cexFrame, textvariable=cexFilesMultiVar, width=8,
                        justify='center')
  
  cexSubLabel <- ttklabel(cexFrame, text='Subplot')
  cexRoiSubVar <- tclVar(newDef$cex.roi.sub)
  roiSubEntry <-	ttkentry(cexFrame, textvariable=cexRoiSubVar, width=8,
                          justify='center')
  
  ##add widgets to grFrame	
  tkgrid(marFrame, row=1, column=1, padx=5, pady=3, sticky='we')
  tkgrid(topLabel, row=2, column=1, padx=3, pady=3, sticky='e')
  tkgrid(bottomLabel, row=3, column=1, padx=3, pady=3, sticky='e')
  tkgrid(leftLabel, row=4, column=1, padx=3, pady=3, sticky='e')
  tkgrid(rightLabel, row=5, column=1, padx=3, pady=3, sticky='e')
  
  tkgrid(marMainLabel, row=1, column=2, padx=4, pady=c(0, 1))
  tkgrid(topMarEntry, row=2, column=2, padx=4, pady=1)
  tkgrid(botMarEntry, row=3, column=2, padx=4, pady=1)
  tkgrid(leftMarEntry, row=4, column=2, padx=4, pady=1)
  tkgrid(rightMarEntry, row=5, column=2, padx=4, pady=1)
  
  tkgrid(marSubLabel, row=1, column=3, padx=4, pady=c(0, 1))
  tkgrid(topSubMarEntry, row=2, column=3, padx=4, pady=1)
  tkgrid(botSubMarEntry, row=3, column=3, padx=4, pady=1)
  tkgrid(leftSubMarEntry, row=4, column=3, padx=4, pady=1)
  tkgrid(rightSubMarEntry, row=5, column=3, padx=4, pady=1)
  
  tkgrid(marMultiLabel, row=1, column=4, padx=c(4, 7), pady=c(0, 1))
  tkgrid(topMultiMarEntry, row=2, column=4, padx=c(4, 7), pady=1)
  tkgrid(botMultiMarEntry, row=3, column=4, padx=c(4, 7), pady=1)
  tkgrid(leftMultiMarEntry, row=4, column=4, padx=c(4, 7), pady=1)
  tkgrid(rightMultiMarEntry, row=5, column=4, padx=c(4, 7), pady=1)
  
  tkgrid(gridFrame, row=1, column=2, padx=5, pady=3, sticky='nwe')
  tkgrid(xgridButton, row=1, column=1, padx=4, pady=5)
  tkgrid(ygridButton, row=2, column=1, padx=4, pady=5)
  
  tkgrid(cexFrame, row=2, column=1, columnspan=2, padx=5, pady=6, sticky='we')
  tkgrid(titleLabel, row=2, column=1, padx=c(8, 2), pady=1, sticky='e')
  tkgrid(axesLabel, row=3, column=1, padx=c(8, 2), pady=1, sticky='e')
  tkgrid(cexMainLabel, row=1, column=2, padx=3, pady=c(0, 1))
  tkgrid(titleEntry, row=2, column=2, padx=3, pady=1, sticky='e')
  tkgrid(axesEntry, row=3, column=2, padx=3, pady=1, sticky='e')
  
  tkgrid(roiMultiLabel, row=2, column=3, padx=c(8, 3), pady=1, sticky='e')
  tkgrid(fileLabel, row=3, column=3, padx=c(8, 3), pady=1, sticky='e')
  tkgrid(cexMultiLabel, row=1, column=4, pady=c(0, 1))
  tkgrid(roiMultiEntry, row=2, column=4, padx=3, pady=1, sticky='e')
  tkgrid(fileEntry, row=3, column=4, padx=3, pady=1, sticky='e')
  
  tkgrid(cexSubLabel, row=1, column=6, padx=c(12, 4), pady=c(0, 1))
  tkgrid(roiSubEntry, row=2, column=6, padx=c(12, 4), pady=1, sticky='e')
  
  #####create widgets for coFrame
  ##color change function
  colorChange <- function(colorVars, canvases){
    usrColor <- tclvalue(tcl("tk_chooseColor", parent=dlg, 
                             initialcolor=tclvalue(colorVars[[1]])))
    if (nzchar(usrColor)){
      for (i in seq_along(colorVars))
        tclObj(colorVars[[i]]) <- usrColor
      for (i in seq_along(canvases))
        tkconfigure(canvases[[i]], background=usrColor)
    }
  }
  
  ##create roi box color buttons
  roiColFrame <- ttklabelframe(coFrame, text='ROI colors')
  boxColLabel <- ttklabel(roiColFrame, text='Boxes:')
  aboxColVar <- tclVar(newDef$roi.bcolor[1])
  aboxColButton <- ttkbutton(roiColFrame, text='Active', width=9, 
                             command=function(...) colorChange(list(aboxColVar), list(aboxCanvas)))
  aboxCanvas <- tkcanvas(roiColFrame, width=35, height=20, 
                         background=newDef$roi.bcolor[1], relief='sunken', borderwidth=2)
  
  iboxColVar <- tclVar(newDef$roi.bcolor[2])
  iboxColButton <- ttkbutton(roiColFrame, text='Inactive', width=9, 
                             command=function(...) colorChange(list(iboxColVar), list(iboxCanvas)))
  iboxCanvas <- tkcanvas(roiColFrame, width=35, height=20, 
                         background=newDef$roi.bcolor[2], relief='sunken', borderwidth=2)
  
  ##create text color buttons
  textColLabel <- ttklabel(roiColFrame, text='Labels:')
  atextColVar <- tclVar(newDef$roi.tcolor[1])
  atextColButton <- ttkbutton(roiColFrame, text='Active', width=9, 
                              command=function(...) colorChange(list(atextColVar), list(atextCanvas)))
  atextCanvas <- tkcanvas(roiColFrame, width=35, height=20, 
                          background=newDef$roi.tcolor[1], relief='sunken', borderwidth=2)
  
  itextColVar <- tclVar(newDef$roi.tcolor[2])
  itextColButton <- ttkbutton(roiColFrame, text='Inactive', width=9, 
                              command=function(...) colorChange(list(itextColVar), list(itextCanvas)))
  itextCanvas <- tkcanvas(roiColFrame, width=35, height=20, 
                          background=newDef$roi.tcolor[2], relief='sunken', borderwidth=2)
  
  ##create axis color button
  fgColVar <- labColVar <- mainColVar <- subColVar <- colVar <- axColVar <- 
    tclVar(newDef$col.axis)
  axColButton <- ttkbutton(coFrame, text='Axes', width=11, 
                           command=function(...) colorChange(list(axColVar, fgColVar, labColVar,
                                                                  mainColVar, subColVar, colVar), list(axCanvas)))
  axCanvas <- tkcanvas(coFrame, width=40, height=20, background=newDef$col.axis,
                       relief='sunken', borderwidth=2)
  
  ##create peak color button
  peakColVar <- tclVar(newDef$peak.color)
  peakColButton <- ttkbutton(coFrame, text='Peak labels', width=11, 
                             command=function(...) colorChange(list(peakColVar), list(peakCanvas)))
  peakCanvas <- tkcanvas(coFrame, width=40, height=20, 
                         background=newDef$peak.color, relief='sunken', borderwidth=2)
  
  ##create background color button
  bgColVar <- tclVar(newDef$bg)
  bgColButton <- ttkbutton(coFrame, text='BG', width=11, 
                           command=function(...) colorChange(list(bgColVar), list(bgCanvas)))
  bgCanvas <- tkcanvas(coFrame, width=40, height=20, background=newDef$bg, 
                       relief='sunken', borderwidth=2)
  
  ##create projection color button
  projColVar <- tclVar(newDef$proj.color)
  projColButton <- ttkbutton(coFrame, text='1D', width=11, 
                             command=function(...) colorChange(list(projColVar), list(projCanvas)))
  projCanvas <- tkcanvas(coFrame, width=40, height=20, 
                         background=newDef$proj.color, relief='sunken', borderwidth=2)
  
  ##create positive contour color button
  posColVar <- tclVar(newDef$pos.color)
  posColButton <- ttkbutton(coFrame, text='+ Contour', width=11, 
                            command=function(...) colorChange(list(posColVar), list(posCanvas)))
  posCanvas <- tkcanvas(coFrame, width=40, height=20, 
                        background=newDef$pos.color, relief='sunken', borderwidth=2)
  
  ##create negative contour color button
  negColVar <- tclVar(newDef$neg.color)
  negColButton <- ttkbutton(coFrame, text='- Contour', width=11, 
                            command=function(...) colorChange(list(negColVar), list(negCanvas)))
  negCanvas <- tkcanvas(coFrame, width=40, height=20, 
                        background=newDef$neg.color, relief='sunken', borderwidth=2)
  
  ##create high contrast colors button
  onContrast <- function(){
    tclObj(aboxColVar) <- 'red'
    tkconfigure(aboxCanvas, background='red')
    tclObj(iboxColVar) <- 'black'
    tkconfigure(iboxCanvas, background='black')
    tclObj(atextColVar) <- 'red'
    tkconfigure(atextCanvas, background='red')
    tclObj(itextColVar) <- 'black'
    tkconfigure(itextCanvas, background='black')
    tclObj(axColVar) <- 'black'
    tclObj(fgColVar) <- 'black'
    tclObj(labColVar) <- 'black'
    tclObj(mainColVar) <- 'black'
    tclObj(subColVar) <- 'black'
    tclObj(colVar) <- 'black'	
    tkconfigure(axCanvas, background='black')
    tclObj(bgColVar) <- 'white'
    tkconfigure(bgCanvas, background='white')
    tclObj(posColVar) <- 'black'
    tkconfigure(posCanvas, background='black')
    tclObj(peakColVar) <- 'black'
    tkconfigure(peakCanvas, background='black')
    tclObj(projColVar) <- 'black'
    tkconfigure(projCanvas, background='black')
    tclObj(negColVar) <- 'black'
    tkconfigure(negCanvas, background='black')
  }
  hcColButton <- ttkbutton(coFrame, text='Hight Contrast', width=15, 
                           command=onContrast)
  
  ##create default colors button
  onDefCol <- function(){
    defSet <- createObj('defaultSettings', returnObj=TRUE)
    tclObj(aboxColVar) <- defSet$roi.bcolor[1]
    tkconfigure(aboxCanvas, background=defSet$roi.bcolor[1])
    tclObj(iboxColVar) <- defSet$roi.bcolor[2]
    tkconfigure(iboxCanvas, background=defSet$roi.bcolor[2])
    tclObj(atextColVar) <- defSet$roi.tcolor[1]
    tkconfigure(atextCanvas, background=defSet$roi.tcolor[1])
    tclObj(itextColVar) <- defSet$roi.tcolor[2]
    tkconfigure(itextCanvas, background=defSet$roi.tcolor[2])
    tclObj(axColVar) <- defSet$col.axis
    tclObj(fgColVar) <- defSet$fg
    tclObj(labColVar) <- defSet$col.lab
    tclObj(mainColVar) <- defSet$col.main
    tclObj(subColVar) <- defSet$col.sub
    tclObj(colVar) <- defSet$col
    tkconfigure(axCanvas, background='white')
    tclObj(bgColVar) <- defSet$bg
    tkconfigure(bgCanvas, background=defSet$bg)
    tclObj(posColVar) <- defSet$pos.color
    tkconfigure(posCanvas, background=defSet$pos.color)
    tclObj(peakColVar) <- defSet$peak.color
    tkconfigure(peakCanvas, background=defSet$peak.color)
    tclObj(projColVar) <- defSet$proj.color
    tkconfigure(projCanvas, background=defSet$proj.color)
    tclObj(negColVar) <- defSet$neg.color
    tkconfigure(negCanvas, background=defSet$neg.color)
  }
  defColButton <- ttkbutton(coFrame, text='Default Colors', width=15, 
                            command=onDefCol)
  
  ##add widgets to coFrame
  tkgrid(roiColFrame, column=1, columnspan=4, row=1, padx=10, pady=10)
  tkgrid(boxColLabel, column=1, row=1, padx=c(15, 3), pady=4)
  tkgrid(aboxColButton, column=2, row=1, padx=2, pady=4)
  tkgrid(aboxCanvas, column=3, row=1, padx=c(2, 20), pady=4)
  tkgrid(iboxColButton, column=4, row=1, padx=2, pady=4)
  tkgrid(iboxCanvas, column=5, row=1, padx=c(2, 15), pady=4)
  
  tkgrid(textColLabel, column=1, row=2, padx=c(15, 3), pady=4)
  tkgrid(atextColButton, column=2, row=2, padx=2, pady=4)
  tkgrid(atextCanvas, column=3, row=2, padx=c(2, 20), pady=4)
  tkgrid(itextColButton, column=4, row=2, padx=2, pady=4)
  tkgrid(itextCanvas, column=5, row=2, padx=c(2, 15), pady=4)
  
  tkgrid(axColButton, column=1, row=2, pady=c(10, 2), padx=c(35, 2))
  tkgrid(axCanvas, column=2, row=2, pady=c(10, 2), padx=c(2, 15))
  
  tkgrid(bgColButton, column=1, row=3, pady=2, padx=c(35, 2))
  tkgrid(bgCanvas, column=2, row=3, pady=2, padx=c(2, 15))
  
  tkgrid(posColButton, column=1, row=4, pady=2, padx=c(35, 2))
  tkgrid(posCanvas, column=2, row=4, pady=2, padx=c(2, 15))
  
  tkgrid(peakColButton, column=3, row=2, pady=c(10, 2), padx=c(15, 2))
  tkgrid(peakCanvas, column=4, row=2, pady=c(10, 2), padx=c(2, 35))
  
  tkgrid(projColButton, column=3, row=3, pady=2, padx=c(15, 2))
  tkgrid(projCanvas, column=4, row=3, pady=2, padx=c(2, 35))
  
  tkgrid(negColButton, column=3, row=4, pady=2, padx=c(15, 2))
  tkgrid(negCanvas, column=4, row=4, pady=2, padx=c(2, 35))
  
  tkgrid(hcColButton, column=1, columnspan=2, row=5, pady=c(15, 5), 
         padx=c(77, 0), sticky='w')
  tkgrid(defColButton, column=3, columnspan=2, row=5, pady=c(15, 5), 
         padx=c(4, 0), sticky='w')
  
  
  
  
  #####create widgets for psFrame
  ##create plot type frame
  typeFrame <- ttklabelframe(psFrame, text='Plot type', padding=4)
  typeVar <- tclVar(switch(newDef$type, 'auto'='auto', 'image'='image', 
                           'contour'='contour', 'filled'='filled contour', 'l'='line', 
                           'p'='points', 'b'='both'))
  typeBox <- ttkcombobox(typeFrame, textvariable=typeVar, values=c('auto', 
                                                                   'image', 'contour', 'filled contour', 'line', 'points', 'both'), 
                         exportselection=FALSE, width=11, state='readonly')
  
  ##create 1D settings frame
  onedFrame <- ttklabelframe(psFrame, text='1D settings', padding=4)
  pos1DLabel <- ttklabel(onedFrame, text='Baseline (0 - 99):')
  pos1DVar <- tclVar(newDef$position.1D)
  pos1DEntry <- ttkentry(onedFrame, textvariable=pos1DVar, width=12,
                         justify='center')
  
  offLabel <- ttklabel(onedFrame, text='Offset (-100 - 100):')
  offVar <- tclVar(newDef$offset)
  offEntry <- ttkentry(onedFrame, textvariable=offVar, width=12,
                       justify='center')
  
  ##create 2D settings frame
  twodFrame <- ttklabelframe(psFrame, text='2D settings', padding=4)
  conLabel <- ttklabel(twodFrame, text='Contour display:')
  if (all(newDef$conDisp))
    conVar <- tclVar('both')
  else if (newDef$conDisp[1])
    conVar <- tclVar('positive')
  else
    conVar <- tclVar('negative')
  conEntry <- ttkcombobox(twodFrame, textvariable=conVar, values=c('positive', 
                                                                   'negative', 'both'), width=9, exportselection=FALSE, state='readonly')
  
  clevelLabel <- ttklabel(twodFrame, text='Contour threshold:\n(positive)')
  clevelVar <- tclVar(newDef$clevel)
  clevelEntry <- ttkentry(twodFrame, textvariable=clevelVar, width=12,
                          justify='center')
  
  nlevelsLabel <- ttklabel(twodFrame, text='Contour levels:\n(0 - 1000)')
  nlevelsVar <- tclVar(newDef$nlevels)
  nlevelsEntry <- ttkentry(twodFrame, textvariable=nlevelsVar, width=12,
                           justify='center')
  
  ##create projection settings frame
  projFrame <- ttklabelframe(psFrame, text='Projections', padding=4)
  filterLabel <- ttklabel(projFrame, text='Type:')
  if (isTRUE(all.equal(newDef$filter, function(x){max(abs(x))})))
    filterVar <- tclVar('absolute max')
  else if (isTRUE(all.equal(newDef$filter, pseudo1D)))
    filterVar <- tclVar('pseudo1D')
  else if (isTRUE(all.equal(newDef$filter,  function(x){max(x)})))
    filterVar <- tclVar('max')
  else
    filterVar <- tclVar('min')
  
  filterBox <- ttkcombobox(projFrame, textvariable=filterVar, width=11, 
                           values=c('pseudo1D', 'max', 'min', 'absolute max'), exportselection=FALSE, 
                           state='readonly')
  
  dispLabel <- ttklabel(projFrame, text='Display:')
  dispVar <- tclVar(switch(newDef$proj.type, 'l'='line', 'p'='points', 
                           'b'='both'))
  dispBox <- ttkcombobox(projFrame, textvariable=dispVar, width=11, 
                         values=c('line', 'points', 'both'), exportselection=FALSE, 
                         state='readonly')
  
  dimLabel <- ttklabel(projFrame, text='Dimension:')
  if (newDef$proj.direct == 1)
    dimVar <- tclVar('direct')
  else
    dimVar <- tclVar('indirect')
  dimBox <- ttkcombobox(projFrame, textvariable=dimVar, width=11,	
                        values=c('direct', 'indirect'), exportselection=FALSE, state='readonly')
  
  ##add widgets to psFrame
  tkgrid(typeFrame, row=1, column=1, padx=4, pady=c(6, 0), sticky='ns')
  tkgrid(typeBox, pady=2, padx=5)
  
  tkgrid(onedFrame, row=2, column=1, padx=4, pady=c(4, 6), sticky='ns')
  tkgrid(pos1DLabel, row=1, column=1, pady=2, sticky='w')
  tkgrid(pos1DEntry, row=2, column=1, pady=3, padx=c(6, 0))
  tkgrid(offLabel, row=3, column=1, pady=c(6, 2), sticky='w')
  tkgrid(offEntry, row=4, column=1, pady=3, padx=c(6, 0))
  
  tkgrid(twodFrame, row=1, rowspan=2, column=2, padx=4, pady=6, sticky='ns')
  tkgrid(conLabel, row=1, column=1, columnspan=2, pady=c(0, 2), sticky='w')
  tkgrid(conEntry, row=2, column=1, columnspan=2,  pady=2, padx=c(6, 4))
  tkgrid(clevelLabel, row=3, column=1, columnspan=2,  pady=c(4, 2), sticky='w')
  tkgrid(clevelEntry, row=4, column=1, columnspan=2,  pady=2, padx=c(6, 4))
  tkgrid(nlevelsLabel, row=5, column=1, columnspan=2,  pady=c(4, 2), sticky='w')
  tkgrid(nlevelsEntry, row=6, column=1, columnspan=2,  pady=2, padx=c(6, 4))
  
  tkgrid(projFrame, row=1, rowspan=2, column=3, padx=4, pady=6, sticky='ns')
  tkgrid(filterLabel, row=1, column=1, pady=c(4, 6), sticky='w')
  tkgrid(filterBox, row=2, column=1, padx=6)
  tkgrid(dispLabel, row=3, column=1, pady=c(8, 6), sticky='w')
  tkgrid(dispBox, row=4, column=1, padx=6)
  tkgrid(dimLabel, row=5, column=1, pady=c(8, 6), sticky='w')
  tkgrid(dimBox, row=6, column=1, padx=6)
  
  #####create widgets for ppFrame
  ##create peak pch entry box
  markerFrame <- ttklabelframe(ppFrame, text='Peak markers', padding=4)
  pchLabel <- ttklabel(markerFrame, text='Symbol:')
  pchVar <- tclVar(newDef$peak.pch)
  pchEntry <- ttkentry(markerFrame, textvariable=pchVar, width=10, 
                       justify='center')
  
  ##create peak cex entry box
  peakCexLabel <- ttklabel(markerFrame, text='Magnification:')
  peakCexVar <- tclVar(newDef$peak.cex)
  peakCexEntry <- ttkentry(markerFrame, textvariable=peakCexVar, width=10, 
                           justify='center')
  
  ##create peak label position combo box
  peakPosLabel <- ttklabel(markerFrame, text='Label position:')
  peakPosVar <- tclVar(newDef$peak.labelPos)
  peakPosBox <- ttkcombobox(markerFrame, textvariable=peakPosVar, width=7, 
                            values=c('top', 'bottom', 'left', 'right', 'center'), 
                            exportselection=FALSE, state='readonly')
  
  ##create peak noiseFilt combo box
  pickSetFrame <- ttklabelframe(ppFrame, text='Pick settings', padding=4)
  peakFiltLabel <- ttklabel(pickSetFrame, text='Noise filter:')
  if (newDef$peak.noiseFilt == 2)
    peakFiltVar <- tclVar('strong')
  else if (newDef$peak.noiseFilt == 1)
    peakFiltVar <- tclVar('weak')
  else
    peakFiltVar <- tclVar('none')
  peakFiltBox <- ttkcombobox(pickSetFrame, textvariable=peakFiltVar, width=7, 
                             values=c('none', 'weak', 'strong'), exportselection=FALSE, 
                             state='readonly')
  
  ##create peak threshold entry box
  peakThreshLabel <- ttklabel(pickSetFrame, text='1D threshold:')
  peakThreshVar <- tclVar(newDef$thresh.1D)
  peakThreshEntry <- ttkentry(pickSetFrame, textvariable=peakThreshVar, 
                              width=10, justify='center')
  
  ##add widgets to ppFrame
  tkgrid(markerFrame, column=1, row=1, padx=10, pady=8)
  tkgrid(pchLabel, column=1, row=1, pady=c(0, 6), sticky='w')
  tkgrid(pchEntry, column=1, row=2, padx=c(15, 5))
  tkgrid(peakCexLabel, column=1, row=3, pady=c(10, 6), sticky='w')
  tkgrid(peakCexEntry, column=1, row=4, padx=c(15, 5), pady=c(0, 4))
  tkgrid(peakPosLabel, column=2, row=1, pady=c(0, 6), sticky='w')
  tkgrid(peakPosBox, column=2, row=2, padx=c(15, 8))
  
  tkgrid(pickSetFrame, column=2, row=1, pady=8)
  tkgrid(peakFiltLabel, column=1, row=1, pady=c(0, 6), sticky='w')
  tkgrid(peakFiltBox, column=1, row=2, padx=c(15, 6))
  tkgrid(peakThreshLabel, column=1, row=3, pady=c(10, 6), sticky='w')
  tkgrid(peakThreshEntry, column=1, row=4, padx=c(15, 5), pady=c(0, 4))
  
  #####create widgets for roiFrame
  ##create appearance frame and labels
  appFrame <- ttklabelframe(roiFrame, text='Appearance')
  activeLabel <- ttklabel(appFrame, text='Active')
  inactiveLabel <- ttklabel(appFrame, text='Inactive')
  
  ##create box type combo boxes
  ltyLabel <- ttklabel(appFrame, text='Box type:')
  altyVar <- tclVar(newDef$roi.lty[1])
  altyBox <- ttkcombobox(appFrame, textvariable=altyVar, width=8, 
                         values=c('solid', 'dashed', 'dotted', 'dotdash', 'longdash', 'twodash', 
                                  'blank'), justify='center', exportselection=FALSE, state='readonly')
  iltyVar <- tclVar(newDef$roi.lty[2])
  iltyBox <- ttkcombobox(appFrame, textvariable=iltyVar, width=8, 
                         values=c('solid', 'dashed', 'dotted', 'dotdash', 'longdash', 'twodash', 
                                  'blank'), justify='center', exportselection=FALSE, state='readonly')
  
  ##create line width entry boxes
  lwdLabel <- ttklabel(appFrame, text='Line width:')
  alwdVar <- tclVar(newDef$roi.lwd[1])
  alwdEntry <- ttkentry(appFrame, textvariable=alwdVar, width=11, 
                        justify='center')
  ilwdVar <- tclVar(newDef$roi.lwd[2])
  ilwdEntry <- ttkentry(appFrame, textvariable=ilwdVar, width=11, 
                        justify='center')
  
  ##create text magnification entry boxes
  roiCexLabel <- ttklabel(appFrame, text='Magnification:')
  acexVar <- tclVar(newDef$roi.cex[1])
  acexEntry <- ttkentry(appFrame, textvariable=acexVar, width=11, 
                        justify='center')
  icexVar <- tclVar(newDef$roi.cex[2])
  icexEntry <- ttkentry(appFrame, textvariable=icexVar, width=11, 
                        justify='center')
  
  ##create label horizontal adjustment entry box
  roiPosLabel <- ttklabel(appFrame, text='Label position:')
  roiPosVar <- tclVar(newDef$roi.labelPos[1])
  roiPosBox <- ttkcombobox(appFrame, textvariable=roiPosVar, width=9, 
                           values=c('top', 'bottom', 'left', 'right', 'center'), 
                           exportselection=FALSE, state='readonly')
  
  ##configure roi size widgets
  onSize <- function(){
    
    ##configure padding widgets
    if (tclvalue(fixedW1) != '0' && tclvalue(fixedW2) != '0')
      padState <- 'disabled'
    else
      padState <- 'normal'
    tkconfigure(roiPadLabel, state=padState)
    tkconfigure(roiPadEntry, state=padState)
    
    ##configure w1 size widgets
    if (tclvalue(fixedW1) == '0'){
      w1State <- 'disabled'
      tclObj(roiW1Var) <- 0
    }else
      w1State <- 'normal'
    tkconfigure(roiW1Entry, state=w1State)
    
    ##configure w2 size widgets
    if (tclvalue(fixedW2) == '0'){
      w2State <- 'disabled'
      tclObj(roiW2Var) <- 0
    }else
      w2State <- 'normal'
    tkconfigure(roiW2Entry, state=w2State)
  }
  
  ##create ROI fixed w1 size checkbutton
  autoFrame <- ttklabelframe(roiFrame, text='Auto generation settings', 
                             padding=4)
  roiSizeLabel <- ttklabel(autoFrame, text='ROI size:')
  if (newDef$roi.w1 == 0){
    fixedW1 <- tclVar(0)
    w1State <- 'disabled'
  }else{
    fixedW1 <- tclVar(1)
    w1State <- 'normal'
  }
  roiW1Button <- ttkcheckbutton(autoFrame, variable=fixedW1, text='Fixed W1', 
                                command=onSize)	
  
  ##create ROI fixed w2 size checkbutton
  if (newDef$roi.w2 == 0){
    fixedW2 <- tclVar(0)
    w2State <- 'disabled'
  }else{
    fixedW2 <- tclVar(1)
    w2State <- 'normal'
  }
  roiW2Button <- ttkcheckbutton(autoFrame, variable=fixedW2, text='Fixed W2', 
                                command=onSize)	
  
  ##create w1 size entry box
  roiW1Var <- tclVar(newDef$roi.w1)
  roiW1Entry <- ttkentry(autoFrame, textvariable=roiW1Var, width=7, 
                         justify='center', state=w1State)
  
  ##create w2 size entry box
  roiW2Var <- tclVar(newDef$roi.w2)
  roiW2Entry <- ttkentry(autoFrame, textvariable=roiW2Var, width=7, 
                         justify='center', state=w2State)
  
  ##create w2 size entry box
  if (newDef$roi.w1 && newDef$roi.w2)
    padState <- 'disabled'
  else
    padState <- 'normal'
  roiPadLabel <- ttklabel(autoFrame, text='Padding (%)', state=padState)
  roiPadVar <- tclVar(newDef$roi.pad)
  roiPadEntry <- ttkentry(autoFrame, textvariable=roiPadVar, width=7, 
                          justify='center', state=padState)
  
  ##create noise filter entry box
  roiFiltLabel <- ttklabel(autoFrame, text='Noise filter:')
  if (newDef$roi.noiseFilt == 2)
    roiFiltVar <- tclVar('strong')
  else if (newDef$roi.noiseFilt == 1)
    roiFiltVar <- tclVar('weak')
  else
    roiFiltVar <- tclVar('none')
  roiFiltBox <- ttkcombobox(autoFrame, textvariable=roiFiltVar, width=7, 
                            values=c('none', 'weak', 'strong'), exportselection=FALSE, 
                            state='readonly')
  
  ##add widgets to roiFrame
  tkgrid(appFrame, column=1, row=1, padx=6, pady=8, sticky='we')
  tkgrid(activeLabel, column=2, row=1, pady=1)
  tkgrid(inactiveLabel, column=3, row=1, pady=1)
  
  tkgrid(ltyLabel, column=1, row=2, padx=3, pady=1, sticky='e')
  tkgrid(altyBox, column=2, row=2, padx=4, pady=1)
  tkgrid(iltyBox, column=3, row=2, padx=4, pady=1)
  
  tkgrid(lwdLabel, column=1, row=3, padx=3, pady=1, sticky='e')
  tkgrid(alwdEntry, column=2, row=3, padx=4, pady=1)
  tkgrid(ilwdEntry, column=3, row=3, padx=4, pady=1)
  
  tkgrid(roiCexLabel, column=1, row=4, padx=3, pady=1, sticky='e')
  tkgrid(acexEntry, column=2, row=4, padx=4, pady=c(1, 4))
  tkgrid(icexEntry, column=3, row=4, padx=4, pady=c(1, 4))
  
  tkgrid(roiPosLabel, column=4, row=1, padx=8, pady=1)
  tkgrid(roiPosBox, column=4, row=2, padx=c(16, 6), pady=3)
  
  tkgrid(autoFrame, column=1, row=2, padx=6, pady=8, sticky='we')
  tkgrid(roiSizeLabel, column=1, row=1, pady=2, padx=4, sticky='w')
  tkgrid(roiW1Button, column=1, row=2, padx=6)
  tkgrid(roiW1Entry, column=1, row=3, padx=6, pady=c(3, 5))
  
  tkgrid(roiW2Button, column=2, row=2, padx=6)
  tkgrid(roiW2Entry, column=2, row=3, padx=6, pady=c(3, 5))
  
  tkgrid(roiPadLabel, column=3, row=2, padx=6)
  tkgrid(roiPadEntry, column=3, row=3, padx=6, pady=c(3, 5))
  
  tkgrid(roiFiltLabel, column=4, row=1, padx=6, pady=3, sticky='w')
  tkgrid(roiFiltBox, column=4, row=2, padx=15, pady=3)
  
  #####create widgets for bottomFrame
  ##create help button
  bottomFrame <- ttkframe(dlg)
  onHelp <- function(){
    myHelp('user_manual', TRUE)
  }
  helpButton <- ttkbutton(bottomFrame, text='?', width=3, command=onHelp)
  
  ##saves current set of preferences
  savePref <- function(){
    
    ##get values for the variables in the GUI
    varList <- list(list(mainWdVar, mainHtVar), list(subWdVar, subHtVar), 
                    list(multiWdVar, multiHtVar), sdiVar, updateVar, backupVar, dirVar, 
                    list(botMarVar, leftMarVar, topMarVar, rightMarVar), 
                    list(botSubMarVar, leftSubMarVar, topSubMarVar, rightSubMarVar), 
                    list(botMultiMarVar, leftMultiMarVar, topMultiMarVar, rightMultiMarVar),
                    xtckVar, ytckVar, cexMainVar, cexAxesVar, cexFilesMultiVar, 
                    cexRoiMultiVar, cexRoiSubVar,  list(aboxColVar, iboxColVar), 
                    list(atextColVar, itextColVar), axColVar, fgColVar, labColVar, 
                    mainColVar, subColVar, colVar, peakColVar, bgColVar, projColVar, 
                    posColVar, negColVar, typeVar, pos1DVar, offVar, conVar, clevelVar, 
                    nlevelsVar, filterVar, dispVar, dimVar, pchVar, peakCexVar, peakPosVar, 
                    peakFiltVar, peakThreshVar, list(alwdVar, ilwdVar), 
                    list(altyVar, iltyVar),	list(acexVar, icexVar), roiPosVar, roiW1Var, 
                    roiW2Var, roiPadVar, roiFiltVar)
    valList <- as.list(rep(NA, length(varList)))
    for (i in seq_along(varList)){
      if (length(varList[[i]]) > 1){
        vectorVals <- NULL
        for (j in varList[[i]])
          vectorVals <- c(vectorVals, tclvalue(j))
        valList[[i]] <- vectorVals
      }else
        valList[[i]] <- tclvalue(varList[[i]])
    }
    varNames <- c('size.main', 'size.sub', 'size.multi', 'sdi', 'update', 
                  'autoBackup', 'wd', 'mar', 'mar.sub', 'mar.multi', 'xtck', 'ytck', 
                  'cex.main', 'cex.axis', 'cex.files.multi', 'cex.roi.multi', 
                  'cex.roi.sub', 'roi.bcolor', 'roi.tcolor', 'col.axis', 'fg', 'col.lab', 
                  'col.main', 'col.sub', 'col', 'peak.color', 'bg', 'proj.color', 
                  'pos.color', 'neg.color', 'type', 'position.1D', 'offset', 'conDisp', 
                  'clevel', 'nlevels', 'filter', 'proj.type', 'proj.direct', 'peak.pch', 
                  'peak.cex', 'peak.labelPos', 'peak.noiseFilt', 'thresh.1D', 'roi.lwd', 
                  'roi.lty', 'roi.cex', 'roi.labelPos', 'roi.w1', 'roi.w2', 'roi.pad', 
                  'roi.noiseFilt')
    names(valList) <- varNames
    
    ##format values from psFrame
    valList$type <- switch(valList$type, 'auto'='auto', 'image'='image', 
                           'contour'='contour', 'filled contour'='filled', 'line'='l', 
                           'points'='p', 'both'='b')
    if (valList$conDisp == 'both')
      valList$conDisp <- c(TRUE, TRUE)
    else if (valList$conDisp == 'positive')
      valList$conDisp <- c(TRUE, FALSE)
    else
      valList$conDisp <- c(FALSE, TRUE)
    if (valList$filter == 'absolute max')
      valList$filter <- function(x){max(abs(x))}
    else if (valList$filter == 'pseudo1D')
      valList$filter <- pseudo1D
    else if (valList$filter == 'max')
      valList$filter <- function(x){max(x)}
    else if (valList$filter == 'min')
      valList$filter <- function(x){min(x)}
    else
      valList$filter <- 'custom'
    valList$proj.type <- unlist(strsplit(valList$proj.type, ''))[1]
    if (valList$proj.direct == 'direct')
      valList$proj.direct <- 1
    else
      valList$proj.direct <- 2
    
    ##format values from grFrame
    if (valList$xtck == '1')
      valList$xtck <- 1
    else
      valList$xtck <- NA_real_
    if (valList$ytck == '1')
      valList$ytck <- 1
    else
      valList$ytck <- NA_real_
    
    ##format values from coFrame
    colorNames <- c('roi.bcolor', 'roi.tcolor', 'col.axis', 'fg', 'col.lab', 
                    'col.main', 'col.sub', 'col', 'peak.color', 'bg', 'proj.color', 
                    'pos.color', 'neg.color')
    for (i in colorNames){
      colVal <- valList[[i]]
      for (j in seq_along(colVal)){
        if (length(grep('{', colVal[j], fixed=TRUE)))
          valList[[i]][j] <- strsplit(strsplit(colVal[j], '{', 
                                               fixed=TRUE)[[1]][2], '}', fixed=TRUE)[[1]][1]
      }
    }
    
    ##format values from ppFrame
    if (valList$peak.noiseFilt == 'strong')
      valList$peak.noiseFilt <- 2
    else if (valList$peak.noiseFilt == 'weak')
      valList$peak.noiseFilt <- 1
    else
      valList$peak.noiseFilt <- 0
    
    ##format values from roiFrame
    if (valList$roi.noiseFilt == 'strong')
      valList$roi.noiseFilt <- 2
    else if (valList$roi.noiseFilt == 'weak')
      valList$roi.noiseFilt <- 1
    else
      valList$roi.noiseFilt <- 0
    
    ##check, assign, and write out the new settings
    newDef[varNames] <- valList[varNames]
    newDef <- checkDef(newDef)
    prevSdi <- defaultSettings$sdi
    myAssign('defaultSettings', newDef, save.backup=FALSE)
    writeDef(defSet=newDef)
    
    ##apply changes to globalSettings
    newGlobal <- globalSettings
    globalPars <- c('offset', 'position.1D', 'filter', 'proj.direct', 'proj.mode', 
                    'proj.type', 'peak.disp', 'peak.noiseFilt', 'thresh.1D', 'peak.pch', 
                    'peak.cex', 'peak.labelPos', 'roiMain', 'roiMax', 'roi.bcolor', 
                    'roi.tcolor', 'roi.lwd', 'roi.lty', 'roi.cex', 'roi.labelPos', 
                    'roi.noiseFilt', 'roi.w1', 'roi.w2', 'roi.pad', 'cex.roi.multi', 
                    'cex.files.multi', 'cex.roi.sub', 'size.main', 'size.sub', 'size.multi', 
                    'mar', 'mar.sub', 'mar.multi', 'overlay.text', 'overlay.textSuppress', 'processSpeciesID', 
                    'processSingleFile', 'speciesList', 'speciesFiles', 'sd_noise_multiplier', 'geneDisp', 'vectorType', 'plotStyle')		
    
    for (i in globalPars)
      newGlobal[i] <- defaultSettings[i]
    myAssign('globalSettings', newGlobal)
    
    ##apply SDI/MDI setting
    if (.Platform$OS.type == 'windows'  && .Platform$GUI == 'Rgui' && 
        prevSdi != defaultSettings$sdi){
      
      ## Reads in the Rconsole file from the R user directory
      conPath <- file.path(Sys.getenv('R_USER'), 'Rconsole')
      if (file.exists(file.path(Sys.getenv('R_USER'), 'Rconsole'))){
        readCon <- file(conPath)
        conText <- readLines(readCon, warn=FALSE)
        close(readCon)
        file.remove(conPath)
        
        ## Reads in the Rconsole file from the R home directory
      }else if (file.access(file.path(R.home('etc'), 'Rconsole'), 2) == 0){
        conPath <- file.path(R.home('etc'), 'Rconsole')
        readCon <- file(conPath)
        conText <- readLines(readCon, warn=FALSE)
        close(readCon)
        file.remove(conPath)
        
        ## Copies Rconsole file from R home directory to R user directory
      }else{
        conPath <- file.path(R.home('etc'), 'Rconsole')
        readCon <- file(conPath)
        conText <- readLines(readCon, warn=FALSE)
        close(readCon)
        file.copy(conPath, file.path(Sys.getenv('R_USER'), 'Rconsole'))
      }
      
      ## Writes out a new Rconsole file
      outFile <- conText
      file.create(conPath)
      writeCon <- file(conPath, 'w')
      matches <- NULL
      if (prevSdi){
        prevMdi <- 'no'
        newMdi <- 'yes'
      }else{
        prevMdi <- 'yes'
        newMdi <- 'no'
      }
      sdiOptions <- paste(c('MDI = ', 'MDI= ', 'MDI =', 'MDI='), prevMdi, 
                          sep='')
      for (i in sdiOptions)
        matches <- c(matches, length(grep(i, outFile)) != 0)
      if (any(matches)){
        for (i in sdiOptions)
          outFile <- gsub(i, paste('MDI =', newMdi), outFile)
        writeLines(outFile, writeCon)
      }else
        writeLines(outFile, writeCon)
      close(writeCon)
    }
  }
  
  ##updates current set of preferences to match defaultSettings
  onDefault <- function(tab=NULL, defSet=defaultSettings){
    
    ##update general tab
    if (is.null(tab) || tab == 0){
      tclvalue(mainWdVar) <- defSet$size.main[1]
      tclvalue(mainHtVar) <- defSet$size.main[2]
      tclvalue(subWdVar) <- defSet$size.sub[1]
      tclvalue(subHtVar) <- defSet$size.sub[2]
      tclvalue(multiWdVar) <- defSet$size.multi[1]
      tclvalue(multiHtVar) <- defSet$size.multi[2]
      tclvalue(sdiVar) <- as.character(defSet$sdi)
      tclvalue(updateVar) <- as.character(defSet$update)
      tclvalue(backupVar) <- as.character(defSet$autoBackup)
      tclvalue(dirVar) <- defSet$wd
    }
    
    ##reset graphics tab
    if (is.null(tab) || tab == 1){
      tclvalue(topMarVar) <- defSet$mar[3]
      tclvalue(botMarVar) <- defSet$mar[1]
      tclvalue(leftMarVar) <- defSet$mar[2]
      tclvalue(rightMarVar) <- defSet$mar[4]
      tclvalue(topSubMarVar) <- defSet$mar.sub[3]
      tclvalue(botSubMarVar) <- defSet$mar.sub[1]
      tclvalue(leftSubMarVar) <- defSet$mar.sub[2]
      tclvalue(rightSubMarVar) <- defSet$mar.sub[4]
      tclvalue(topMultiMarVar) <- defSet$mar.multi[3]
      tclvalue(botMultiMarVar) <- defSet$mar.multi[1]
      tclvalue(leftMultiMarVar) <- defSet$mar.multi[2]
      tclvalue(rightMultiMarVar) <- defSet$mar.multi[4]
      if (is.na(defSet$xtck))
        tclvalue(xtckVar) <- 0
      else
        tclvalue(xtckVar) <- 1
      if (is.na(defSet$ytck))
        tclvalue(ytckVar) <- 0
      else
        tclvalue(ytckVar) <- 1
      tclvalue(cexMainVar) <- defSet$cex.main
      tclvalue(cexAxesVar) <- defSet$cex.axis
      tclvalue(cexRoiMultiVar) <- defSet$cex.roi.multi
      tclvalue(cexFilesMultiVar) <- defSet$cex.files.multi
      tclvalue(cexRoiSubVar) <- defSet$cex.roi.sub
    }
    
    ##reset colors tab
    if (is.null(tab) || tab == 2){
      tclvalue(aboxColVar) <- defSet$roi.bcolor[1]
      tkconfigure(aboxCanvas, background=defSet$roi.bcolor[1])
      tclvalue(iboxColVar) <- defSet$roi.bcolor[2]
      tkconfigure(iboxCanvas, background=defSet$roi.bcolor[2])
      tclvalue(atextColVar) <- defSet$roi.tcolor[1]
      tkconfigure(atextCanvas, background=defSet$roi.tcolor[1])
      tclvalue(itextColVar) <- defSet$roi.tcolor[2]
      tkconfigure(itextCanvas, background=defSet$roi.tcolor[2])
      tclvalue(fgColVar) <- tclvalue(labColVar) <- tclvalue(mainColVar) <- 
        tclvalue(subColVar) <- tclvalue(colVar) <- tclvalue(axColVar) <- 
        defSet$col.axis
      tkconfigure(axCanvas, background=defSet$col.axis)
      tclvalue(peakColVar) <- defSet$peak.color
      tkconfigure(peakCanvas, background=defSet$peak.color)
      tclvalue(bgColVar) <- defSet$bg
      tkconfigure(bgCanvas, background=defSet$bg)
      tclvalue(projColVar) <- defSet$proj.color
      tkconfigure(projCanvas, background=defSet$proj.color)
      tclvalue(posColVar) <- defSet$pos.color
      tkconfigure(posCanvas, background=defSet$pos.color)
      tclvalue(negColVar) <- defSet$neg.color
      tkconfigure(negCanvas, background=defSet$neg.color)
    }
    
    ##reset plotting tab
    if (is.null(tab) || tab == 3){
      tclvalue(typeVar) <- switch(defSet$type, 'auto'='auto', 'image'='image', 
                                  'contour'='contour', 'filled'='filled contour', 'l'='line', 
                                  'p'='points', 'b'='both')
      tclvalue(pos1DVar) <- defSet$position.1D
      tclvalue(offVar) <- defSet$offset
      if (all(defSet$conDisp))
        tclvalue(conVar) <- 'both'
      else if (defSet$conDisp[1])
        tclvalue(conVar) <- 'positive'
      else
        tclvalue(conVar) <- 'negative'
      tclvalue(clevelVar) <- defSet$clevel
      tclvalue(nlevelsVar) <- defSet$nlevels
      if (isTRUE(all.equal(defSet$filter, function(x){max(abs(x))})))
        tclvalue(filterVar) <- 'absolute max'
      else if (isTRUE(all.equal(defSet$filter, pseudo1D)))
        tclvalue(filterVar) <- 'pseudo1D'
      else if (isTRUE(all.equal(defSet$filter,  function(x){max(x)})))
        tclvalue(filterVar) <- 'max'
      else
        tclvalue(filterVar) <- 'min'
      tclvalue(dispVar) <- switch(defSet$proj.type, 'l'='line', 'p'='points', 
                                  'b'='both')
      if (defSet$proj.direct == 1)
        tclvalue(dimVar) <- 'direct'
      else
        tclvalue(dimVar) <- 'indirect'
    }
    
    ##reset peak picking tab
    if (is.null(tab) || tab == 4){
      tclvalue(pchVar) <- defSet$peak.pch
      tclvalue(peakCexVar) <- defSet$peak.cex
      tclvalue(peakPosVar) <- defSet$peak.labelPos
      if (defSet$peak.noiseFilt == 2)
        tclvalue(peakFiltVar) <- 'strong'
      else if (defSet$peak.noiseFilt == 1)
        tclvalue(peakFiltVar) <- 'weak'
      else
        tclvalue(peakFiltVar) <- 'none'
      tclvalue(peakThreshVar) <- defSet$thresh.1D
    }
    
    ##reset ROIs tab
    if (is.null(tab) || tab == 5){
      tclvalue(altyVar) <- defSet$roi.lty[1]
      tclvalue(iltyVar) <- defSet$roi.lty[2]
      tclvalue(alwdVar) <- defSet$roi.lwd[1]
      tclvalue(ilwdVar) <- defSet$roi.lwd[2]
      tclvalue(acexVar) <- defSet$roi.cex[1]
      tclvalue(icexVar) <- defSet$roi.cex[2]
      tclvalue(roiPosVar) <- defSet$roi.labelPos[1]
      if (defSet$roi.w1)
        tclvalue(fixedW1) <- 1
      else
        tclvalue(fixedW1) <- 0
      tclvalue(roiW1Var) <- defSet$roi.w1
      if (defSet$roi.w2)
        tclvalue(fixedW2) <- 1
      else
        tclvalue(fixedW2) <- 0
      tclvalue(roiW2Var) <- defSet$roi.w2
      tclvalue(roiPadVar) <- defSet$roi.pad
      if (defSet$roi.noiseFilt == 2)
        tclvalue(roiFiltVar) <- 'strong'
      else if (defSet$roi.noiseFilt == 1)
        tclvalue(roiFiltVar) <- 'weak'
      else
        tclvalue(roiFiltVar) <- 'none'
    }
  }
  
  ##create default button
  defaultButton <- ttkbutton(bottomFrame, text='Defaults', width=10, 
                             command=function(...) onDefault(as.numeric(tkindex(epBook, 'current')), 
                                                             createObj('defaultSettings', returnObj=TRUE)))
  
  ##create OK button
  onOk <- function(){
    savePref()
    tkdestroy(dlg)
  }
  okButton <- ttkbutton(bottomFrame, text='OK', width=10, command=onOk)
  
  ##create cancel button
  cancelButton <- ttkbutton(bottomFrame, text='Cancel', command=function(...)
    tkdestroy(dlg), width=10)
  
  ##create apply all button
  onApply <- function(){
    
    ##save preferences
    savePref()
    
    ##close and redisplay splash screen (if open)
    prevDev <- dev.list()
    if (is.null(fileFolder) || !length(fileFolder)){
      if (length(prevDev)){
        for (i in prevDev)
          dev.off(i)
        if (2 %in% prevDev){
          if (.Platform$OS == 'windows')
            dev.new(title='Main Plot Window', 
                    width=defaultSettings$size.main[1],
                    height=defaultSettings$size.main[2])
          else
            X11(title='Main Plot Window', width=defaultSettings$size.main[1],	
                height=defaultSettings$size.main[2])
          #splashScreen()
        }
      }
      return(invisible())
    }
    
    ##update fileFolder
    newFolder <- fileFolder
    for (i in seq_along(newFolder)){
      newFolder[[i]]$graphics.par <- defaultSettings
      newFolder[[i]]$graphics.par$usr <- fileFolder[[i]]$graphics.par$usr
    }
    myAssign('fileFolder', newFolder)
    
    ##close and reopen previously displayed devices
    if (length(prevDev)){
      for (i in prevDev)
        dev.off(i)
      if (2 %in% prevDev)
        dd()
      if (3 %in% prevDev)
        rvs()
      if (4 %in% prevDev)
        rvm()
    }
    
    ##update preferences GUI
    onDefault()
  }
  applyButton <- ttkbutton(bottomFrame, text='Apply All', width=10, 
                           command=onApply)
  
  ##add button to bottom of gui
  tkgrid(bottomFrame, column=1, row=2, pady=c(0, 8), sticky='we')
  tkgrid(helpButton, column=1, row=2, padx=c(15, 10))
  tkgrid(defaultButton, column=2, row=2, padx=c(0, 40))
  tkgrid(okButton, column=3, row=2, padx=c(0, 6))
  tkgrid(cancelButton, column=4, row=2, padx=c(0, 6))
  tkgrid(applyButton, column=5, row=2, padx=c(0, 15))
  
  ##make buttons a bit smaller on non-windows systems
  if (.Platform$OS.type != 'windows'){
    tkconfigure(browseButton, width=9)
    tkconfigure(defaultButton, width=9)
    tkconfigure(okButton, width=9)
    tkconfigure(cancelButton, width=9)
    tkconfigure(applyButton, width=9)
  }
  
  return(invisible())
}


################################################################################
##                                                                            ##
##     Internal functions that run when the digestR package is loaded         ##
##                                                                            ##
################################################################################

#' Update My Package
#'
#' This function updates the R package by installing an updated version from GitHub.
#'
#' If the package "digestR" is currently loaded, it will be detached before updating.
#'
#' @importFrom devtools install_github
#'
#' @return Invisible NULL.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' update_my_package()
#' }
#'
#' @keywords internal
updater <- function() {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  library(devtools)
  
  # Check if the package is loaded
  if ("digestR" %in% installed.packages()) {
    if ("digestR" %in% search()) {
      cat("Detaching package 'digestR'...\n")
      detach("package:digestR", unload = TRUE)
    }
    cat("Updating package 'digestR'...\n")
    devtools::install_github("LewisResearchGroup/digestR")
  } else {
    cat("Installing package 'digestR'...\n")
    devtools::install_github("LewisResearchGroup/digestR")
  }
  
  cat("Package updated successfully!\n")
  invisible(NULL)
}

## Updates digestR old
updater_old <- function(auto=FALSE){
  return(NULL)  # switch off updates for now

  ##display message
  if (auto){
    if (!defaultSettings$update)
      return(invisible())
    cat('\nChecking digestR for updates . . . ')
  }
  
  ##check for installed packages
  pkgs <- as.data.frame(installed.packages(), stringsAsFactors=FALSE)
  digestRpkgs <- pkgs[grep('digestR', rownames(pkgs)), ]
  digestRpaths <- digestRpkgs[, 'LibPath']
  digestRvers <- digestRpkgs[, 'Version']
  
  ##check for write permission
  writeDirs <- which(file.access(digestRpaths, mode=2) == 0)
  if (!length(writeDirs)){
    if (auto){
      cat('cancelled.\n  User does not have write permission.\n')
      return(invisible())
    }else
      err(paste('Could not update digestR.  You do not have write permission for ', 
                'the digestR package directory.\nTo update digestR you must run R as ', 
                'an administrator or change the write permissions for the ', 
                'directory below.\n\n', 'digestR library location:  ', 
                digestRpaths[writeDirs[1]], sep=''))
  }
  writeVers <- digestRvers[writeDirs]
  
  ##check for the newest version of digestR currently installed
  currVers <- writeVers[1]
  if (length(writeDirs) > 1){
    for (i in 2:length(writeVers)){
      cv <- compareVersion(currVers, writeVers[i])
      if (cv < 0)
        currVers <- writeVers[i]
    }
  }
  digestRpath <- digestRpaths[match(currVers, digestRvers)]	
  
  ##check if the package is up to date
  newVers <- suppressWarnings(
    available.packages(
      contriburl=contrib.url(repos='http://digestR.nmrfam.wisc.edu/R/', type='source')
      )
  )

  if (is.null(newVers) || !length(newVers) || !nzchar(newVers)){
    if (auto){
      cat('cancelled.\n  Could not access digestR repository.\n')
      return(invisible())
    }else
      err(paste('Could not access digestR repository.\n', 
                'Check your network settings and try again.', sep=''))
  }

  newVers <- newVers['digestR', 'Version']
  
  ##download and install the latest digestR package if update was selected manually
  if (!auto){
    detach(package:digestR)
    tryCatch(install.packages('digestR', lib=digestRpath, 
                              repos='http://digestR.nmrfam.wisc.edu/R/', type='source'),
             error=function(er)
               stop(paste('Could not update digestR.\nGo to the digestR homepage to ',
                          'download and manually install the latest version.', 
                          sep='')), call.=FALSE)
    
    ##check installation
    tryCatch(require(digestR, quietly=TRUE, warn.conflicts=FALSE), 	
             error=function(er)
               stop(paste('Could not update digestR.\nGo to the digestR homepage to ',
                          'download and manually install the latest version.', 
                          sep='')), call.=FALSE)
    currVers <- suppressWarnings(packageDescription('digestR', fields='Version', 
                                                    lib.loc=digestRpath))
    if (currVers == newVers){
      try(Sys.chmod(system.file('linux/digestR.sh', package='digestR'), mode = '555'), 
          silent=TRUE)
      myMsg('digestR update successful.  Restart R to apply the changes.   ', 
            icon='info')
    }
    return(invisible())
  }
  
  ##download and install digestR if the current package is not up-to-date
  if (compareVersion(newVers, currVers) < 1){
    if (auto)
      cat('done.\n')
    return(invisible())
  }
  usrSel <- myMsg(paste('digestR version ', newVers, ' is now available.\n', 
                        'Version ', pkgVar$version, ' is currently installed.\n', 
                        'Would you like to update digestR?', sep=''), type='yesno')
  if (usrSel == 'no'){
    if (auto)
      cat('cancelled.\n')
    return(invisible())
  }
  
  tryCatch(update.packages(repos='http://digestR.nmrfam.wisc.edu/R/', ask=FALSE,
                           type='source', lib.loc=digestRpath),
           error=function(er){ 
             require(digestR, quietly=TRUE, warn.conflicts=FALSE)
             myMsg(paste('Could not update digestR.\nGo to the ',
                         'digestR homepage to download and manually install the latest ', 
                         'version.', sep=''), icon='error')
           })
  
  ##check installation
  currVers <- suppressWarnings(packageDescription('digestR', fields='Version', 
                                                  lib.loc=digestRpath))
  if (currVers == newVers){
    try(Sys.chmod(system.file('linux/digestR.sh', package='digestR'), mode = '555'), 
        silent=TRUE)
    myMsg('digestR update successful.  Please restart R to apply changes.', 
          icon='info')
    q('no')
  }
  if (auto)
    cat('done.\n')
  
  return(invisible())
}

## Check for updates to BMRB standards library
updateLib <- function(auto=FALSE){
  # Remove update for testing
  return(invisible())

  ##exit if auto updates are turned off
  if (auto && !defaultSettings$libUpdate)
    return(invisible())
  
  ##display update message
  cat('\nChecking for updates to standards library . . .')
  flush.console()
  
  ##check for write permission
  libDir <- dirname(createObj('defaultSettings', returnObj=TRUE)$libLocs)
  if (file.access(libDir, mode=2)){
    if (auto){
      cat('cancelled.\n  User does not have write permission.\n')
      return(invisible())
    }else
      err(paste('Could not update library.  You do not have write permission ', 
                'for the digestR package directory.\nTo update the standards ',
                'library you must run R as an administrator or change the write ',
                'permissions for the directory below:', libDir, sep=''))
  }
  
  ##read remote library index file
  libUrl <- 'http://digestR.nmrfam.wisc.edX/pages/data/files/RSD_libraries'
  remoteIndex <- try(read.table(file.path(libUrl, 'index.txt'), head=TRUE, 
                                sep='\t', stringsAsFactors=FALSE), silent=TRUE)
  if (is.null(remoteIndex$Name)){
    if (auto){
      cat('cancelled.\n  Could not access digestR server.\n')
      return(invisible())
    }else
      err(paste('Could not access digestR server.\n', 
                'Check your network settings and try again.', sep=''))
  }
  remoteRsds <- file.path(remoteIndex$Library, remoteIndex$Name)
  
  ##read local library index file
  localIndex <- read.table(file.path(libDir, 'index.txt'), head=TRUE, 
                           sep='\t', stringsAsFactors=FALSE)
  localRsds <- file.path(localIndex$Library, localIndex$Name)
  
  ##check for new additions
  newFiles <- which(!remoteRsds %in% localRsds)
  if (length(newFiles))
    newFiles <- remoteRsds[newFiles]
  
  ##check for updated RSD files
  remoteIndex$Updated <- as.Date(remoteIndex$Updated)
  localIndex$Updated <- as.Date(localIndex$Updated)
  fileUpdates <- NULL
  for (i in seq_along(localRsds)){
    remoteMatch <- match(localRsds[i], remoteRsds)
    if (localIndex$Updated[i] < remoteIndex$Updated[remoteMatch])
      fileUpdates <- c(fileUpdates, remoteRsds[remoteMatch])
  }
  newFiles <- c(fileUpdates, newFiles)
  if (!length(newFiles)){
    cat('\nStandards library is up-to-date.\n')
    return(invisible())
  }
  
  ##ask user if they would like to update library
  if (auto){
    updateMsg <- paste('Updates are available for the digestR standards library.', 
                       'Would you like to download the updates now?\n', sep='\n')
    usrSel <- buttonDlg(updateMsg, c('Yes', 'No', 
                                     'Don\'t display this message again'), TRUE, default='No')
    
    ##edit defaultSettings if user doesn't want message to be displayed again
    if (as.logical(usrSel[2])){
      defaultSettings$libUpdate <- FALSE
      writeDef(defSet=defaultSettings)
      myAssign('defaultSettings', defaultSettings, FALSE)
    }
    if (usrSel[1] == 'No'){
      cat('cancelled.\n')
      return(invisible())
    }
  }
  
  ##create directories for new libraries
  updatedLibs <- unique(dirname(newFiles))
  for (i in updatedLibs){
    localLibPath <- file.path(libDir, i)
    if (!file.exists(localLibPath))
      dir.create(localLibPath, FALSE, TRUE)
    
    ##download library files
    cat('\n')
    download.file(paste(libUrl, i, 'metadata.txt', sep='/'), 
                  file.path(localLibPath, 'metadata.txt'))
    download.file(paste(libUrl, i, 'details.txt', sep='/'), 
                  file.path(localLibPath, 'details.txt'))
    download.file(paste(libUrl, i, 'roiRecord.roi', sep='/'), 
                  file.path(localLibPath, 'roiRecord.roi'))
    download.file(paste(libUrl, i, 'hash.ucsf', sep='/'), 
                  file.path(localLibPath, 'hash.ucsf'))
  }
  
  ##download RSD files
  for (i in newFiles)
    download.file(file.path(libUrl, i), file.path(libDir, i), mode='wb')
  
  ##download index file
  download.file(file.path(libUrl, 'index.txt'), file.path(libDir, 
                                                          'index.txt'))
  cat('Library update complete.\n')
  
  return(invisible())
}

## Checks R version (2.13.1) for presence of bug in image() function
checkImage <- function(){
  
  ##get R version
  ver <- R.Version()
  if (ver$major == '2' && ver$minor == '13.1'){
    
    ##create main window
    dlg <- tktoplevel()
    tcl('wm', 'attributes', dlg, topmost=TRUE)
    tkwm.resizable(dlg, FALSE, FALSE)
    tkwm.title(dlg, 'digestR - WARNING')
    tkfocus(dlg)
    tkwm.deiconify(dlg)
    
    ##create font for text
    fonts <- tcl('font', 'name')
    if (!'msgFont' %in% as.character(fonts)){
      msgFont <- as.character(tcl('font', 'configure', 'TkDefaultFont'))
      tcl('font', 'create', 'msgFont', msgFont[1], msgFont[2], msgFont[3], 
          '10', msgFont[5], msgFont[6], msgFont[7], msgFont[8], 
          msgFont[9], msgFont[10], msgFont[11], msgFont[12])
    }
    
    ##create text box
    textFrame <- ttkframe(dlg)
    textBox <- tktext(textFrame, wrap='word', height=13, width=60,
                      font='msgFont', cursor='arrow')
    
    ##add text to widget
    msg <- paste('R version 2.13.1 contains a bug in the image() function', 
                 'that may prevent certain spectra from being displayed correctly.  ',
                 'After a file is opened, the spectrum may be displayed with horizontal',
                 'or vertical lines running across it, or the spectrum may not be', 
                 'displayed at all.\n\nTo avoid this issue, choose a plot type option', 
                 'other than "auto" or "image" (select "Plot settings" from the', 
                 'Graphics menu).  This issue is not present in previous versions of',
                 'R and we expect that it will be resolved with the next R release. ', 
                 'Previous versions of R may be downloaded using the links below.\n\n')
    tkinsert(textBox, '1.0', msg)
    tkinsert(textBox, '10.0', '                 ')
    tkinsert(textBox, '10.end', 'Windows installer', 'winUrl')
    tkinsert(textBox, '10.end', '                 ')
    tkinsert(textBox, '10.end', 'Mac OS X installer', 'macUrl')
    tcl(textBox, 'tag', 'configure', 'winUrl', foreground='blue')
    tcl(textBox, 'tag', 'configure', 'macUrl', foreground='blue')
    tkconfigure(textBox, state='disabled')
    
    ##configure windows installer text to display webpage when clicked
    tcl(textBox, 'tag', 'bind', 'winUrl', '<Enter>', function(...){
      tcl(textBox, 'tag', 'configure', 'winUrl', underline=TRUE)
      tkconfigure(textBox, cursor='hand2')
    })
    tcl(textBox, 'tag', 'bind', 'winUrl', '<Leave>', function(...){
      tcl(textBox, 'tag', 'configure', 'winUrl', underline=FALSE)
      tkconfigure(textBox, cursor='arrow')
    })
    winPage <- 'http://cran.opensourceresources.org/bin/windows/base/old/2.13.0'
    tcl(textBox, 'tag', 'bind', 'winUrl', '<Button-1>', function(...)
      browseURL(winPage))
    
    ##configure mac installer text to display webpage when clicked
    tcl(textBox, 'tag', 'bind', 'macUrl', '<Enter>', function(...){
      tcl(textBox, 'tag', 'configure', 'macUrl', underline=TRUE)
      tkconfigure(textBox, cursor='hand2')
    })
    tcl(textBox, 'tag', 'bind', 'macUrl', '<Leave>', function(...){
      tcl(textBox, 'tag', 'configure', 'macUrl', underline=FALSE)
      tkconfigure(textBox, cursor='arrow')
    })
    macPage <- 'http://cran.opensourceresources.org/bin/macosx/old/R-2.13.0.pkg'
    tcl(textBox, 'tag', 'bind', 'macUrl', '<Button-1>', function(...)
      browseURL(macPage))
    
    ##create ok button
    okButton <- ttkbutton(dlg, text='OK', width=12, default='active', 
                          command=function(...) tkdestroy(dlg))
    
    ##add widgets to toplevel window
    tkgrid(textFrame, row=1, sticky='nswe', pady=8, padx=8)
    tkgrid(textBox)
    tkgrid(okButton, row=2, pady=c(4, 10))
  }
  
  return(invisible())
}

## Executes a set of tasks whenever the digestR package loads
.onLoad <- function(lib, pkg){

  #log_message('Loading digestR backend functions')

  ## Create or update necessary digestR objects
  digestR:::patch()
  
  ## Exit if digestR has not been installed
  if (length(grep('apple', Sys.getenv('R_PLATFORM')))){
    installDir <- tryCatch(find.package('digestR', lib.loc='~'), 
                           error=function(er) return(NULL))
    if (is.null(installDir))
      installDir <- tryCatch(dirname(find.package('digestR')), 
                             error=function(er) return(NULL))
    if (is.null(installDir))
      return(invisible())
      
    defPath <- paste(installDir, '/defaultSettings', sep='')
    if (file.access(defPath) == -1){
      digestR:::writeDef(defPath)
      return(invisible())
    }

  }
  
  ## Autoload functions in digestR namespace
  digestRfun <- c('aa', 'appendPeak', 'bringFocus', 'buttonDlg', 'ca', 'cl', 'cf', 
               'closeGui', 'co', 'ct', 'ctd', 'ctu', 'cw', 'da', 'dd', 'di', 'dp', 'dr', 
               'draw2D', 'drawPeptides', 'drf', 'ed', 'ep', 'err', 'export', 'fc', 'ff', 'fo', 
               'fs', 'getTitles', 'gui', 'hideGui', 'import', 'isNoise', 'loc', 
               'localMax', 'matchShift', 'maxShift', 'mmcd', 'myAssign', 'myDialog', 
               'myDir', 'myFile', 'myMsg', 'myOpen', 'mySave', 'mySelect', 'myToplevel', 
               'nf', 'ol', 'pa', 'paAll', 'pd', 'pDel', 'pDelAll', 'pe', 'peakPick',	
               'peakPick1D', 'peakPick2D', 'peakVolume', 'per', 'ph', 'pj', 'pjv', 'pl', 
               'plot1D', 'plot2D', 'pm', 'pp', 'pr', 'pReg', 'pu', 'pv', 'pw', 'pwAll', 
               'pz', 'ra', 'rb', 'rc', 'rcd', 'rci', 'rd', 'rdAll', 'rDel', 're', 
               'red', 'refresh', 'regionMax', 'rei', 'reset', 'rmd', 'rml', 'rmr', 'rmu', 
               'rn', 'roi', 'rotc', 'rotcc', 'rotd', 'rotu', 'rp', 'rpAll', 'rr', 'rs', 
               'rsAll', 'rsf', 'rSum', 'rv', 'rvm', 'rvs', 'se', 'setGraphics', 
               'setWindow', 'shiftToROI', 'showGui', 'spin', 'sr', 'ss', 'tableEdit', 
               'tclCheck', 'ucsf1D', 'ucsf2D', 'ud', 'vp', 'vpd', 'vpu', 'vs', 'wc', 
               'wl', 'writeUcsf', 'ws', 'zc', 'zf', 'zi', 'zm', 'zo', 'zp', 'zz', 'vd', 'up', 'csp', 'pd', 'cs')
  
  for (i in digestRfun)
    suppressPackageStartupMessages(autoload(i, 'digestR', warn.conflicts=FALSE))
  
  ## Turn on HTML help
  options(htmlhelp=TRUE, help_type='html', chmhelp=FALSE)
    
  ## Check if the Rgui is running in SDI (multiple windows) mode
  if (.Platform$GUI == 'Rgui')
    sdiCheck()
  
  # ## Add digestR to list of repositories and update package if applicable
  # if (is.na(match('digestR', names(getOption('repos')))))
  #   options(repos=c(getOption('repos'), digestR='http://digestR.nmrfam.wisc.edu/R'))
  # errMsg <- tryCatch(digestR:::updater(TRUE), error=function(er) 
  #   return(er$message))
  # if (!is.null(errMsg))
  #   packageStartupMessage('Non-fatal error occurred while checking for', 
  #                         ' updates:\n  "', errMsg, '"\n')
  
  # ## Check for standards library updates
  # errMsg <- tryCatch(digestR:::updateLib(TRUE), error=function(er) 
  #   return(er$message))
  # if (!is.null(errMsg))
  #   packageStartupMessage('Non-fatal error occurred while checking for ', 
  #                         'standards library updates:\n  "', errMsg, '"\n')
  
  ## Check for image() bug
  checkImage()
}

## Perform necessary actions from .onLoad when running digestR from source code
if (!'package:digestR' %in% search() && !exists('fileFolder')){
  
  ##assign digestR objects
  tclCheck()
  patch(FALSE)
  
  ## Set X11 options and display open file message
  if (.Platform$OS == 'windows'){
    dev.new(title='Main Plot Window', width=defaultSettings$size.main[1], 
            height=defaultSettings$size.main[2])
  }else{
    X11.options(type='Xlib')
    X11(title='Main Plot Window', width=defaultSettings$size.main[1],	
        height=defaultSettings$size.main[2])
  }
  
    tryCatch(splashScreen(), error=function(er){
      if (.Platform$OS != 'windows'){
        invisible(myMsg(paste('Your computer does not have the required ', 
                              'fonts to support fast X11 graphics in R.\n',
                              'To correct this issue you may need to download some or', 
                              ' all of the following X11 fonts:     \n\n', 
                              '                              xorg-x11-fonts-75dpi\n',
                              '                              xorg-x11-fonts-100dpi\n', 
                              '                              xorg-x11-fonts-truetype\n',
                              '                              xorg-x11-fonts-Type1\n\n', 
                              'Please refer to the R Installation and Administration',
                              ' Manual for more information:\n', 
                              'http://cran.r-project.org/doc/manuals/R-admin.html#X11-',
                              'issues', 
                              sep=''), 'ok', 'info'))
        dev.off()
        X11.options(type='cairo')
        X11(title='Main Plot Window', width=defaultSettings$size.main[1], 
            height=defaultSettings$size.main[2])
        splashScreen()
      }
    })
  
  ## Use a functional version of ::tk::dialog::file:: on older Linux systems
  tclVer <- as.character(tcl('info', 'patchlevel'))
  tclVer <- unlist(strsplit(tclVer, '.', fixed=TRUE))
  if (tclVer[1] < 8 || (tclVer[1] == 8 && tclVer[2] < 5) ||
      (tclVer[1] == 8 && tclVer[2] == 5 && tclVer[3] < 5)){
    filePath <- system.file('tcltk/tkfbox.tcl', package='digestR')
    tcl('source', filePath)
  }
  
  ## Add the tablelist package to the Tcl search path and load the package
  invisible(addTclPath(system.file('tcltk/tablelist', package='digestR')))
  invisible(tclRequire('tablelist_tile'))
  if (.Platform$OS == 'windows')
    tcl('option', 'add', '*Tablelist*selectBackground', 'SystemHighlight')
  tcl('option', 'add', '*Tablelist*stripeBackground', '#ececff')
  
  ## Correct problems with the "xpnative" theme for the treeview widget
  if (.Platform$OS == 'windows'){
    tcl('ttk::style', 'configure', 'Treeview', '-background', 'white')
    tcl('ttk::style', 'configure', 'Row', '-background', 'white')
    tcl('ttk::style', 'configure', 'Cell', '-background', 'white')
    tcl('ttk::style', 'map', 'Row', 
        '-background', c('selected', 'SystemHighlight'), 
        '-foreground', c('selected', 'SystemHighlightText'))
    tcl('ttk::style', 'map', 'Cell', 
        '-background', c('selected', 'SystemHighlight'), 
        '-foreground', c('selected', 'SystemHighlightText'))
    tcl('ttk::style', 'map', 'Item', 
        '-background', c('selected', 'SystemHighlight'), 
        '-foreground', c('selected', 'SystemHighlightText'))
  }
  
  ## Assign the digestR icon to GUIs
  #createTclImage('digestRIcon', DIGESTR_GIF_PATH)
  tt <- tktoplevel()
  #tcl('wm', 'iconphoto', tt, '-default', 'digestRIcon')
  
  ## Make sure Ttk widgets display the same color background as toplevels
  defBgColor <- as.character(tkcget(tt, '-background'))
  tkdestroy(tt)
  tcl('ttk::style', 'configure', 'TRadiobutton', '-background', defBgColor)
  tcl('ttk::style', 'map', 'TRadiobutton', '-background', c('disabled', 
                                                            defBgColor))
  tcl('ttk::style', 'configure', 'TCheckbutton', '-background', defBgColor)
  tcl('ttk::style', 'map', 'TCheckbutton', '-background', c('disabled', 
                                                            defBgColor))
  tcl('ttk::style', 'configure', 'TSizegrip', '-background', defBgColor)
  tcl('ttk::style', 'map', 'TSizegrip', '-background', c('disabled', 
                                                         defBgColor))
  tcl('ttk::style', 'configure', 'TLabel', '-background', defBgColor)
  tcl('ttk::style', 'map', 'TLabel', '-background', c('disabled', 
                                                      defBgColor))
  tcl('ttk::style', 'configure', 'TNotebook', '-background', defBgColor)
  tcl('ttk::style', 'map', 'TNotebook', '-background', c('disabled', 
                                                         defBgColor))
  tcl('ttk::style', 'configure', 'Treeview', '-background', 'white')
  tcl('ttk::style', 'map', 'Treeview', '-background', c('disabled', 
                                                        defBgColor))
  tcl('ttk::style', 'configure', 'TFrame', '-background', defBgColor)
  tcl('ttk::style', 'configure', 'TLabelframe', '-background', defBgColor)
  #gui()
  
  ## Turn on HTML help
  # options(htmlhelp=TRUE, help_type='html', chmhelp=FALSE)
  
  ## Print message on package load
  #cat("\n", "digestR version", pkgVar$version, "\n", 
  #    "Copyright (C) 2023, Dimitri Desmonts de Lamache, Raied Aburashed, Travis A. Bingemann, Sören Wacker and Ian A. Lewis\n",
  #    "DigestR is free software and comes with ABSOLUTELY NO WARRANTY.\n",
  #    "DigestR may be modified and redistributed under certain conditions.\n",
  #    "Go to http://www.r-project.org/Licenses/GPL-3 for more details.\n\n", 
  #    "Citation:\n",
  #    "Dimitri Desmonts de Lamache, Raied Aburashed, Travis A. Bingemann, Sören Wacker and Ian A. Lewis\n",
  #    "Magn. Reson. Chem. 47, S123-S126 (2023).\n\n")
}


.onAttach <- function(libname, pkgname) {
  ## Print message on package load
  packageStartupMessage("DigestR version 1.0.0", pkgVar$version, "\n", 
                        "Copyright (C) 2023, Dimitri Desmonts de Lamache, Raied Aburashed, Travis A. Bingemann, Sören Wacker and Ian A. Lewis\n",
 			"DigestR is free software and comes with ABSOLUTELY NO WARRANTY.\n",
  			"DigestR may be modified and redistributed under certain conditions.\n",
  		        "Go to http://www.r-project.org/Licenses/GPL-3 for more details.\n\n", 
  			"Citation:\n",
      			"Desmonts de Lamache, D., Aburashed A., Bingemann T. A., Wacker S., and Ian A. Lewis\n",
      			"Magn. Reson. Chem. 47, S123-S126 (2023).\n\n")
}

##################################################
# DigestR functions 

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
  log_message('span')
  log_message(span)
  
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

  if( is.null(lSpecies)) {
    log_message(paste0('Species file ', fileName, ' no special function defined, using generic loading function.'))
    lSpecies <- loadProteome(fileName)
  }

  return(lSpecies)
}

batchConvert <- function()
{

  ## load a species object based on the user choice made in processGUI
  lSpecies <- loadSpecies(globalSettings$speciesFiles[globalSettings$processSpeciesID])
  
  log_message('Loaded species file')

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
              log_message(paste0('Mascot file could not be read or found, ', saveFileName, " could not be generated."))
            else
              log_message(paste0('No matches found, ', saveFileName, " could not be generated."))	
          }else
          {
            sName <- writeDIANA(saveFileName, lSpecies, result) 
            log_message(paste0(sName, ' saved at ', format(Sys.time(), "%H:%M"), '.'))
          }
        }
      }
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
  #log_message(firstGeneIndex)
  #log_message(lastGeneIndex)	
  geneNames = lSpecies$genes$name[firstGeneIndex:lastGeneIndex]
  #log_message(geneNames)
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
      #		log_message("manually call plot1D")
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
      log_message('Insufficient data to generate noise estimate using FDR.')
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
  return(Rprof(NULL))
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
          log_message(w2Range)
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

#' Gene Name Threshold GUI
#'
#' This function opens a graphical user interface (GUI) that allows users to set threshold values for gene names.
#' Users can interact with the GUI to select files, set and apply threshold values, and reset to default values.
#'
#' @return NULL (The function primarily works with GUI elements.)
#'
#' @import tcltk2
#'
#' @export
gene_labeling <- function()
{
  ##creates main window
  tclCheck()
 
  # Destroy any existing window with the same ID to avoid conflicts
  if (as.logical(tcl('winfo', 'exists', '.gl'))) {
    tkdestroy('.gl')  # Destroy the previous window
  }

  dlg <- tktoplevel()
  tkwm.title(dlg, 'Gene Name Threshold')

  # Withdraw the window to prevent flickering during setup
  tkwm.withdraw(dlg)

  tkfocus(dlg)

  # dlg <- myToplevel('gl')
  # if (is.null(dlg))
  #   return(invisible())
  # tkwm.title(dlg, 'Gene Name Threshold')
  # tkfocus(dlg)
  # tkwm.deiconify(dlg)
  
  ##create file list box
  fileFrame <- ttklabelframe(dlg, text='Files')
  fileList <- tclVar()
  fileNames <- names(fileFolder)
  tclObj(fileList) <- getTitles(fileNames)
  fileBox <- tklistbox(fileFrame, height=10, width=34, listvariable=fileList,
                       selectmode='extended', active='dotbox',	exportselection=FALSE, bg='white', 
                       xscrollcommand=function(...) tkset(xscr, ...), yscrollcommand=function(...) tkset(yscr, ...),
		       font = "Helvetica 12 bold",  # Increased font size and made bold
                       bg = "white", fg = "black")
  
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

 ## Now that the window is set up, deiconify it to show the window properly
  tkwm.deiconify(dlg)
  tcl("update")  # Ensure the window is fully drawn and updated
  
  invisible()
}

look_up_species_files <- function(){
  #look up all files in subfolder 'species' and return a list of files
  species <- list.files(path=PROTEOMES_PATH, pattern = '.csv', full.names = TRUE)
  globalSettings$speciesFiles <<- species
  globalSettings$speciesList <<- species
  return(species)
}

#' Process Mascot Files
#'
#' This function provides a GUI for processing Mascot files. 
#' The user can select a species and choose to either convert a single file or all files in a folder 
#' and its subfolders.
#'
#' @return Invisible. The function internally updates global settings based on user inputs 
#' from the GUI and then calls a batch conversion function.
#'
#' @note The function uses the `tcltk` package to create the GUI.
#'
#' @examples
#' \dontrun{
#' pm()
#' }
#'
#' @import tcltk
#' @export
# process_mascot <- function()
# {

#   species <- look_up_species_files()

#   ##creates main window
#   tclCheck()
#   dlg <- myToplevel('pm')
#   if (is.null(dlg))
#     return(invisible())
  
#   tkwm.title(dlg, 'Process Mascot Files')
#   tkfocus(dlg)
#   tkwm.deiconify(dlg)
  
#   speciesLabelFrame <- ttklabelframe(dlg, text='Species')
  
#   speciesList <- tclVar()
#   tclObj(speciesList) <- species
  
#   tmp <- tclvalue(tkfont.actual(dlg))
  
#   tmp2 <- strsplit(tmp, ' ', fixed=TRUE)
#   fam <- tmp2[[1]][2]
#   siz <- tmp2[[1]][4]
#   wgt <- tmp2[[1]][6]
#   slt <- tmp2[[1]][8]
#   slt <- 'italic' ## change to italics
  
#   myFont <- tkfont.create(family = fam, size = siz, weight = wgt, slant = slt)
  
#   speciesBox <- tklistbox(speciesLabelFrame, height=10, width=34, listvariable=speciesList,
#                           selectmode='extended', active='dotbox',	exportselection=FALSE, bg='white', font = myFont,
#                           xscrollcommand=function(...) tkset(xscr, ...), yscrollcommand=function(...) tkset(yscr, ...))
  
#   xscr <- ttkscrollbar(speciesLabelFrame, orient='horizontal', command=function(...) tkxview(speciesBox, ...))
#   yscr <- ttkscrollbar(speciesLabelFrame, orient='vertical',	command=function(...) tkyview(speciesBox, ...))
  
#   if (length(species) > 2)
#   {
#     for (i in seq(0, length(species) - 1, 2))
#       tkitemconfigure(speciesBox, i, background='#ececff')
#   }
  
#   tkselection.set(speciesBox, 0)	
#   #################
#   ##create single vs. multiple file conversion radiobuttons
#   fileConvFrame <- ttklabelframe(dlg, text='Convert File(s):')
#   singleFileConvVal <- tclVar(TRUE)
  
#   singleFileOnButton <- ttkradiobutton(fileConvFrame, variable=singleFileConvVal, 
#                                        value=TRUE, text='Single File')
  
#   singleFileOffButton <- ttkradiobutton(fileConvFrame, variable=singleFileConvVal, 
#                                         value=FALSE, text='All Files in folder and sub folders')
#   #################
#   buttonFrame <- ttkframe(dlg)
  
#   onApply <- function()
#   {
#     idx <- 1 + as.integer(tkcurselection(speciesBox))
#     #selectedSpeciesName <- basename(species[idx]) # Added
#     globalSettings$processSpeciesID <<- idx
#     globalSettings$processSingleFile <<- (as.logical(tclObj(singleFileConvVal)) == TRUE)
#     tkdestroy(dlg)
#     #loadProteome(species[idx], selectedSpeciesName)  # Added Pass the selected file and name
#     batchConvert()
#   }
#   apply <- ttkbutton(buttonFrame, text='Apply', width=10, command=onApply)
  
#   ##add widgets to speciesFrame
#   tkgrid(speciesLabelFrame, column=1, row=1, sticky='nswe', pady=6, padx=6)
#   tkgrid(speciesBox, column=1, row=1, sticky='nswe')
#   tkgrid(yscr, column=2, row=1, sticky='ns')
#   tkgrid(xscr, column=1, row=2, sticky='we')
  
#   ##make fileFrame stretch when window is resized
#   tkgrid.columnconfigure(dlg, 1, weight=1)
#   tkgrid.rowconfigure(dlg, 1, weight=10)
#   tkgrid.columnconfigure(speciesLabelFrame, 1, weight=1)
#   tkgrid.rowconfigure(speciesLabelFrame, 1, weight=1)
  
#   tkgrid(fileConvFrame, column=1, row=2, sticky='we', pady=6, padx=6)
#   tkgrid(singleFileOnButton, column=1, row=1, sticky='nswe', pady=6, padx=6)
#   tkgrid(singleFileOffButton, column=2, row=1, sticky='nswe', pady=6, padx=6)
  
#   tkgrid(buttonFrame, column=1, row=3, sticky='nswe', pady=6, padx=6)
#   tkgrid(apply, column=1, row=1, sticky='w')
#   tkgrid.columnconfigure(buttonFrame, 1, weight=1)
#   tkgrid.rowconfigure(buttonFrame, 1, weight=1)
  
#   tkgrid(ttksizegrip(dlg), column=2, row=3, sticky='se')
  
#   invisible()
# }

####################################################################
# process mascot v3

process_mascot <- function() {
  # Fetch the species list from look_up_species_files
  species <- look_up_species_files()

  # Create or reinitialize the window
  tclCheck()
  
  # Destroy any previous window with the same ID to avoid conflicts
  if (as.logical(tcl('winfo', 'exists', '.pm'))) {
    tkdestroy('.pm')  # This will force the window to be destroyed and recreated
  }

  dlg <- tktoplevel()
  tkwm.title(dlg, 'Process Mascot Files')
  tkwm.geometry(dlg, "600x400+200+400")  # Position window at 200, 200 on screen
  
  # Withdraw the window initially to prevent flickering
  tkwm.withdraw(dlg)
	
  # Create custom font for the listbox items (increase font size)
  #listboxFont <- tkfont.create(family = "Helvetica", size = 12)
  
  ## Species Frame
  speciesLabelFrame <- ttklabelframe(dlg, text='Species')

  speciesList <- tclVar()
  tclObj(speciesList) <- species

  # Create listbox for species
  speciesBox <- tklistbox(speciesLabelFrame, height = 10, width = 34, listvariable = speciesList,
                          selectmode = 'extended', active = 'dotbox', exportselection = FALSE,
			  font = "Helvetica 12 bold",
			  bg = 'white', fg = "black")

  xscr <- ttkscrollbar(speciesLabelFrame, orient = 'horizontal', command = function(...) tkxview(speciesBox, ...))
  yscr <- ttkscrollbar(speciesLabelFrame, orient = 'vertical', command = function(...) tkyview(speciesBox, ...))

  # Alternate row colors
  if (length(species) > 2) {
    for (i in seq(0, length(species) - 1, 2)) {
      tkitemconfigure(speciesBox, i, background = '#ececff')
    }
  }

  tkselection.set(speciesBox, 0)

  ## Create file conversion radiobuttons
  fileConvFrame <- ttklabelframe(dlg, text = 'Convert File(s):')
  singleFileConvVal <- tclVar(TRUE)

  singleFileOnButton <- ttkradiobutton(fileConvFrame, variable = singleFileConvVal, value = TRUE, text = 'Single File')
  singleFileOffButton <- ttkradiobutton(fileConvFrame, variable = singleFileConvVal, value = FALSE, text = 'All Files in folder and sub folders')

  ## Apply button frame
  buttonFrame <- ttkframe(dlg)

  onApply <- function() {
    idx <- 1 + as.integer(tkcurselection(speciesBox))
    globalSettings$processSpeciesID <<- idx
    globalSettings$processSingleFile <<- (as.logical(tclObj(singleFileConvVal)) == TRUE)
    tkdestroy(dlg)
    batchConvert()
  }

  apply <- ttkbutton(buttonFrame, text = 'Apply', width = 10, command = onApply)

  ## Layout using grid
  tkgrid(speciesLabelFrame, column = 1, row = 1, sticky = 'nswe', pady = 6, padx = 6)
  tkgrid(speciesBox, column = 1, row = 1, sticky = 'nswe')
  tkgrid(yscr, column = 2, row = 1, sticky = 'ns')
  tkgrid(xscr, column = 1, row = 2, sticky = 'we')

  ## Make fileFrame stretch when window is resized
  tkgrid.columnconfigure(dlg, 1, weight = 1)
  tkgrid.rowconfigure(dlg, 1, weight = 10)
  tkgrid.columnconfigure(speciesLabelFrame, 1, weight = 1)
  tkgrid.rowconfigure(speciesLabelFrame, 1, weight = 1)

  tkgrid(fileConvFrame, column = 1, row = 2, sticky = 'we', pady = 6, padx = 6)
  tkgrid(singleFileOnButton, column = 1, row = 1, sticky = 'nswe', pady = 6, padx = 6)
  tkgrid(singleFileOffButton, column = 2, row = 1, sticky = 'nswe', pady = 6, padx = 6)

  tkgrid(buttonFrame, column = 1, row = 3, sticky = 'nswe', pady = 6, padx = 6)
  tkgrid(apply, column = 1, row = 1, sticky = 'w')
  tkgrid.columnconfigure(buttonFrame, 1, weight = 1)
  tkgrid.rowconfigure(buttonFrame, 1, weight = 1)

  tkgrid(ttksizegrip(dlg), column = 2, row = 3, sticky = 'se')

  # Now that the window is set up, deiconify it to show the window properly
  tkwm.deiconify(dlg)
  tkraise(dlg)
  tcl("update")  # Ensure the window is fully drawn and updated
  
  # Wait for the window to close
  tkwait.window(dlg)
  
  invisible()
}
		       
# loadProteome <- function(sFilename, selectedSpeciesName) {
#   log_message(sFilename)
#   df <- read.csv(sFilename, head = TRUE, stringsAsFactors = FALSE)
#   fileInfo <- file.info(sFilename)
#   ID <- as.integer(fileInfo$mtime)
#   df <- subset(df, select = c("GeneName", "seq", "chromosome", "start"))
#   names(df)[1] <- "name"
#   #return(prepareSpecies(selectedSpeciesName, ID, df))
#   return(prepareSpecies('Proteome:', ID, df))
# }

loadProteome <- function(sFilename, selectedSpeciesName) {
  log_message(sFilename)
  
  # Try reading the CSV file
  df <- tryCatch({
    read.csv(sFilename, head = TRUE, stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("Error reading the file:", e$message, "\n")
    return(NULL)
  })
  
  # If df is NULL, return early
  if (is.null(df)) {
    return(NULL)
  }
  
  # Define the required columns
  required_columns <- c("GeneName", "seq", "chromosome", "start")
  
  # Check if all required columns are present
  if (!all(required_columns %in% colnames(df))) {
    missing_cols <- setdiff(required_columns, colnames(df))
    cat("Error: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    cat('Check Column Names "GeneName", "seq", "chromosome", "start"\n')
    return(NULL)
  }
  
  # Proceed with subsetting the data if the columns exist
  df <- subset(df, select = c("GeneName", "seq", "chromosome", "start"))
  names(df)[1] <- "name"
  
  # Call prepareSpecies with the updated dataframe
  return(prepareSpecies('Proteome:', as.integer(file.info(sFilename)$mtime), df))
}

		       
#' Save File As DIANA Compressed File
#'
#' This function allows the user to save the current DIANA file in a compressed format.
#'
#' @param saveFileName Character string containing the desired save file name. If not provided,
#'                     a file dialog will prompt the user to choose a save location.
#' @return None (invisible return).
#'
#' @export
save <- function(saveFileName = '')
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
  
  #dlg <- myToplevel('es')
  # if (is.null(dlg))
  # return(invisible())

  # Destroy any existing window with the same ID to avoid conflicts # New
  if (as.logical(tcl('winfo', 'exists', '.es'))) {
    tkdestroy('.es')  # Destroy the previous window # New
  }
  dlg <- tktoplevel() # New
  tkwm.title(dlg, 'Edit Variables')
  tkwm.geometry(dlg, "550x500+200+200")

  # Withdraw the window to prevent flickering while setting up
  tkwm.withdraw(dlg)
	
  tkfocus(dlg)
  #tkwm.deiconify(dlg)
  
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
    log_message(length(points$x1))
  
  
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

#' Plot Settings GUI
ps <- function(dispPane='co'){
  
  ##create main window
  current <- wc()
  tclCheck()

 # Destroy any existing window with the same ID # New
  if (as.logical(tcl('winfo', 'exists', '.ps'))) {
    tkdestroy('.ps')  # Destroy the previous window
  }

  dlg <- tktoplevel() # New
  tkwm.title(dlg, 'Plot Settings') # New
  tkwm.geometry(dlg, "650x500+200+200")  # Increase window size (Width=600, Height=400)

  # Withdraw the window to prevent flickering during setup
  tkwm.withdraw(dlg) # New
  
  tkfocus(dlg) # New

  # dlg <- myToplevel('ps')
  # if (is.null(dlg))
  # {
  #   if (dispPane == 'co')
  #     tkselect('.ps.1', 0)
  #   else
  #     tkselect('.ps.1', 1)
    
  #   return(invisible())
  # }
  
  #tkwm.title(dlg, 'Plot Settings')
  #tkfocus(dlg)
  #tkwm.deiconify(dlg)
  
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
                         exportselection=FALSE, 	
			 xscrollcommand=function(...) tkset(coXscr, ...), 
                         yscrollcommand=function(...) tkset(coYscr, ...),
			 font = "Helvetica 12 bold",  # Increased font size and made bold
                         bg = "white", fg = "black")
  
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
  onedFileFrame <- ttklabelframe(onedFrame, text='Files - Double click to switch spectra')
  onedFileList <- tclVar()
  onedFileNames <- names(fileFolder)[which(sapply(fileFolder, 
                                                  function(x){x$file.par$number_dimensions}) == 1)]
  tclObj(onedFileList) <- getTitles(onedFileNames)
  onedFileBox <- tklistbox(onedFileFrame, width=30, listvariable=onedFileList, 
                           selectmode='extended', active='dotbox',	exportselection=FALSE,
                           xscrollcommand=function(...) tkset(onedXscr, ...), 
                           yscrollcommand=function(...) tkset(onedYscr, ...),
			   font = "Helvetica 12 bold",  # Increased font size and made bold
                           bg = "white", fg = "black")
			 
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
      log_message('Refreshing plot...')
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
        log_message(paste0(geneName, ' is not a valid gene of ', species$name))
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
      if (length(names(fileFolder)) &&
	  currentSpectrum %in% onedFileNames)
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

  ## Now that the window is set up, deiconify it to show the window properly
  tkwm.deiconify(dlg)
  tcl("update")  # Ensure the window is fully drawn and updated
  
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

#' Overlay Spectra and Generate Plot
#'
#' This function overlays spectra based on the provided list of spectrum indices and
#' generates a plot with the overlaid spectra. It supports different types of overlay
#' palettes and allows for customization of the overlay behavior.
#'
#' @param askUsr A logical value indicating whether to prompt the user for overlay list creation.
#'               Default is \code{TRUE}.
#' @param offset A numeric value indicating the offset between overlaid spectra. If not provided,
#'               the offset from global settings is used.
#' @param ... Additional parameters passed to internal plotting functions.
#'
#' @return None (invisible return).
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item It defines default overlay palettes for different numbers of colors.
#'   \item It determines the current spectrum and its dimensions.
#'   \item If \code{askUsr} is \code{TRUE} or \code{overlayList} is \code{NULL}, it opens
#'     the GUI for creating an overlay list and returns.
#'   \item If \code{askUsr} is \code{FALSE} and an overlay list is present, it overlays
#'     the spectra as specified by the overlay list.
#'   \item The function handles plot text overlay, removing the current spectrum, and plotting
#'     the overlaid spectra.
#' }
#'
#' @examples
#' \dontrun{
#'   # Overlay spectra based on the provided list
#'   ol(askUsr = FALSE, offset = 0.1)
#' }
#'
#' @importFrom graphics legend par plot
#' @export
overlay <- function(askUsr = TRUE, offset = NULL, ...)
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

   # Destroy any existing window with the same ID to avoid conflicts # New
  if (as.logical(tcl('winfo', 'exists', '.os'))) {
    tkdestroy('.os')  # Destroy the previous window # New
  } 
  dlg <- tktoplevel() # New
  tkwm.title(dlg, if (dispPane == 'ol') 'Overlays' else 'Other Title') # New
  tkwm.geometry(dlg, "800x500+200+200") # New
  
  # Withdraw the window to prevent flickering during setup
  tkwm.withdraw(dlg) # New

  tkfocus(dlg) #New
	     
  # dlg <- myToplevel('os')
  # if (is.null(dlg))
  # {
  #   if (dispPane == 'ol')
  #   {
  #     tkwm.title('.os', 'Overlays')
  #     tkselect('.os.1', 0)
  #   }
  #   return(invisible())
  # }
  
  # tkfocus(dlg)
  # tkwm.deiconify(dlg)
  # if (dispPane == 'ol')
  #   tkwm.title(dlg, 'Overlays')
  
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
                         exportselection=FALSE,  
                         xscrollcommand=function(...) tkset(olXscr, ...), 
                         yscrollcommand=function(...) tkset(olYscr, ...),
                         font = "Helvetica 12 bold",
			 bg = "white", fg = "black")
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
                          active='dotbox', 
                          xscrollcommand=function(...) tkset(overlayXscr, ...), 
                          yscrollcommand=function(...) tkset(overlayYscr, ...),
                          bg = "white", fg = "black")
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

  ## Now that the window is set up, deiconify it to show the window properly
  tkwm.deiconify(dlg) # New
  tcl("update")  # Ensure the window is fully drawn and updated # New
  
  invisible()
}

#' Manipulate Files using mf
#'
#' This function opens a graphical user interface (GUI) for manipulating files. Users can select multiple files,
#' perform various mathematical operations on them (such as addition, subtraction, multiplication, and division),
#' merge files, and calculate the mean of selected files.
#'
#' @return NULL (The function primarily works with GUI elements.)
#'
#' @import tcltk2
#'
#' @export
mf <- function()
{
  ##create main window
  current <- wc()
  tclCheck()

  # Destroy any existing window with the same ID to avoid conflicts # New
  if (as.logical(tcl('winfo', 'exists', '.mf'))) {
    tkdestroy('.mf')  # Destroy the previous window # New
  }
  dlg <- tktoplevel() # New
  tkwm.title(dlg, 'Manipulate Files') # New
  tkwm.geometry(dlg, "400x300+200+200") # new

  # Withdraw the window to prevent flickering during setup
  tkwm.withdraw(dlg) # New
  
  tkfocus(dlg) # New

  # dlg <- myToplevel('mf')
  
  # tkfocus(dlg)
  # tkwm.deiconify(dlg)
  
  #tkwm.title(dlg, 'Manipulate Files')
  
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
    
    #		log_message(paste0('Selected Files: ', selectedFiles))
    #		log_message(paste0('myFactorVal: ', myFactorVal))
    
    if (!length(usrSel))
      err(paste('You must select one or more files before',
                'pressing the execute button'))
    else
    {
      for(i in usrSel)
        log_message(names(fileFolder)[i])
      
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
                 log_message('Please enter the number of the file you wish to act as the negative operand, in the operand entry field.')
                 
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
  
  ##resets widgets whenever the mouse enters the GUI
  onMouse <- function()
  {
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
  tkwm.deiconify(dlg)
  tcl("update")  # Ensure the window is fully drawn and updated
  
  invisible()
}


checkColumnList <- function(sStr)
{
  return(!grepl("[^[:space:]\\,\\:0-9+]", sStr))
}

im <- function()
{
  tclCheck()
  
 # Destroy any existing window with the same ID to avoid conflicts
  if (as.logical(tcl('winfo', 'exists', '.im'))) {
    tkdestroy('.im')  # Destroy the previous window
  }

  dlg <- tktoplevel()
  tkwm.title(dlg, 'Import Maven File')
  tkwm.geometry(dlg, "300x400+200+200") 

  # Withdraw the window to prevent flickering during setup
  tkwm.withdraw(dlg)
  
  tkfocus(dlg)	
	
  #   dlg <- myToplevel('im')
  #   if (is.null(dlg))
  #   return(invisible())
  
  # tkwm.title(dlg, 'Import Maven File')
  # tkfocus(dlg)
  # tkwm.deiconify(dlg)

  listboxFont <- tkfont.create(family = "Helvetica", size = 12)
  
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
      log_message('Please use the \'browse\' button to select a Maven file to import. ')	
    }else if(!file.exists(tclvalue(fnOutput)))
    {
      log_message('The selected file does not exist. Please use the \'browse\' button to select a new file.')
    }else
    {
      if((nchar(tclvalue(ciEntry)) == 0) || (nchar(tclvalue(biEntry)) == 0) || (nchar(tclvalue(siEntry)) == 0))
      {
        log_message('Please ensure you have input the column numbers for compounds, blanks and samples.')
      }else
      {
        compoundsIdx <- tclvalue(ciEntry)
        blanksIdx <- tclvalue(biEntry)
        samplesIdx <- tclvalue(siEntry)
        thresholdMult <- tclvalue(tmEntry)
        
        if(checkColumnList(compoundsIdx) && checkColumnList(blanksIdx) && checkColumnList(samplesIdx) && is.numeric(thresholdMult) )
        {
          log_message(paste0(compoundsIdx, '_', blanksIdx, '_', samplesIdx, '_', thresholdMult))
          parseMaven(sMavenFileName = tclvalue(fnOutput), compoundsIdx = compoundsIdx, vBlanksIdx = blanksIdx, vSamplesIdx = samplesIdx, kThresholdMult = thresholdMult)
          tkdestroy(dlg)			
        }else
        {
          if(!checkColumnList(compoundsIdx))
            log_message('Please ensure that the compound column list contains only integer values separated by commas or a colon.')
          
          if(!checkColumnList(blanksIdx))
            log_message('Please ensure that the blanks column list contains only integer values separated by commas or a colon.')
          
          if(!checkColumnList(samplesIdx))
            log_message('Please ensure that the samples column list contains only integer values separated by commas or a colon.')
          
          if(!is.numeric(thresholdMult))
          {
            log_message('Please ensure that the threshold multiplier is a numeric value.')
            log_message(paste0('threshold: ', thresholdMult))
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

  tkwm.deiconify(dlg)
  tcl("update")  # Ensure the window is fully drawn and updated
  
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

## Open files without refreshing, creating an undo point, or console messages
## fileNames - character string or vector; full file paths for spectra being
##	opened
## Note - This function behaves differently than fo() and is meant for internal
##	use.  Attempting to open a file that already exists in fileFolder will 
##	create a duplicate entry.  3D spectra are not supported.
# silentOpen <- function(fileList){
  
#   ## Read selected files
#   fileList <- sort(fileList)
#   for (i in seq_along(fileList)){
    
#     ##Read Sparky Header and file info from binary
#     new.file <- try(dianaHead(file.name=fileList[i], print.info=FALSE), 
#                     silent=TRUE)
#     if (!is.list(new.file) || !length(new.file$file.par)){
#       fileList <- fileList[-i]
#       next
#     }	
    
#     ## Fetch the default graphics settings 
#     new.file$graphics.par <- defaultSettings
    
#     ## Set initial plotting range
#     if (new.file$file.par$number_dimensions == 1){
#       new.file$graphics.par$usr <- c(new.file$file.par$downfield_ppm[1],
#                                      new.file$file.par$upfield_ppm[1], 
#                                      new.file$file.par$zero_offset - 
#                                        (new.file$file.par$max_intensity - new.file$file.par$zero_offset) 
#                                      * globalSettings$position.1D,
#                                      new.file$file.par$max_intensity)
#     }else{         
#       new.file$graphics.par$usr <- c(new.file$file.par$downfield_ppm[2],
#                                      new.file$file.par$upfield_ppm[2], 
#                                      new.file$file.par$downfield_ppm[1],
#                                      new.file$file.par$upfield_ppm[1])
#     }    
    
#     ## Change name if spectrum already exists in fileFolder
#     fileNames <- as.vector(sapply(fileFolder, function(x) x$file.par$file.name))
#     matches <- which(fileNames == fileList[i])
#     if (length(matches)) {
#       newName <- paste(fileList[i], '(', length(matches) + 1, ')', sep='')
#       new.file$file.par$user_title <- basename(newName)
#       fileList[i] <- newName
#     }
    
#     ## Add spectrum to fileFolder
#     fileFolder[[(length(fileFolder) + 1)]] <- new.file
#     names(fileFolder)[[(length(fileFolder))]] <- fileList[i]
#   }
  
#   ## Assign the new objects to the global environment
#   myAssign('fileFolder', fileFolder, save.backup=FALSE)
  
#   invisible(fileList)
# }

## Create RSD file
## fileName - character string; full file path for input spectrum
## outPath - character string; full file path to save the output RSD to
## roiTable - character string; ROI table or full path to ROI table text file
createRsd <- function(fileName=currentSpectrum, outPath, rois=roiTable){
  
  if (missing(outPath)){
    outPath <- paste(unlist(strsplit(basename(fileName), '.', fixed=TRUE))[1], 
                     'rsd', sep='.')
    outPath <- file.path(dirname(fileName), outPath)
  }
  suppressWarnings(dir.create(dirname(outPath), recursive=TRUE))
  
  ## Get file data
  inFile <- ucsf2D(fileName)
  fp <- inFile$file.par
  
  ## Get ROI table
  if (is.null(rois))
    err('ROI table does not exist or no ROIs have been created.')
  else if (is.character(rois)){
    if (!file.access(rois, mode=4))
      rois <- read.table(rois,	head=TRUE, sep='\t', stringsAsFactors=FALSE)
    else
      err(paste('File does not exist or can\'t be accessed. ', 
                'You must provide a valid ROI table.'))
  }
  rois <- rois[order(rois$w1_upfield, rois$w2_upfield), ]
  
  ###### Write main header
  writeCon <- file(outPath, "w+b")
  
  ## format name
  writeBin('RSD_NMR', writeCon)
  writeBin(as.integer(rep(0, 2)), writeCon, size=1, endian='big')
  
  ## format version
  writeBin('1.2', writeCon)
  writeBin(as.integer(rep(0, 6)), writeCon, size=1, endian='big')
  
  ## number of dimensions
  writeBin(as.integer(fp$number_dimensions), writeCon, size=4, endian='big')
  
  ## number of blocks (ROIs)
  writeBin(as.integer(nrow(rois)), writeCon, size=4, endian='big')
  
  ## noise estimate from original spectrum
  writeBin(as.numeric(fp$noise_est), writeCon, size=4, endian='big')
  
  ## make header size equal to 100 bytes by filling remainder with 0s
  writeBin(as.integer(rep(0, 68)), writeCon, size=1, endian='big')
  
  ###### Write original spectrum header data for each dimension 
  for (i in 1:fp$number_dimensions){
    
    ## dimension index
    writeBin(as.integer(i - 1), writeCon, size=4, endian='big')
    
    ## number of points along this dimension
    writeBin(as.integer(fp$matrix_size[i]), writeCon, size=4, endian='big')
    
    ## downfield chemical shift in this dimension
    writeBin(as.numeric(fp$downfield_ppm[i]), writeCon, size=4, endian='big')
    
    ## upfield chemical shift in this dimension
    writeBin(as.numeric(fp$upfield_ppm[i]), writeCon, size=4, endian='big')
    
    ## spectral width in this dimension (Hz)
    writeBin(as.numeric(fp$spectrum_width_Hz[i]), writeCon, size=4, 
             endian='big')
    
    ## spectrometer frequency in this dimension (MHz)
    writeBin(as.numeric(fp$transmitter_MHz[i]), writeCon, size=4, endian='big')
    
    ## nucleus name in this dimension
    writeBin(as.character(fp$nucleus[i]), writeCon)
    writeBin(as.integer(rep(0, (10 - nchar(fp$nucleus[i]) - 1))), writeCon, 
             size=1, endian='big')
    
    ## make header size equal to 50 bytes by filling remainder with 0s
    writeBin(as.integer(rep(0, 16)), writeCon, size=1, endian='big')
  }
  
  ###### Write block headers (contains information on each ROI)
  for (i in 1:fp$number_dimensions){
    for (j in 1:nrow(rois)){
      
      ## dimension index
      writeBin(as.integer(i - 1), writeCon, size=4, endian='big')
      
      ## get upfield and downfield chemical shifts for current block
      if (fp$number_dimensions == 1)
        shifts <- matchShift(inFile, w2=c(rois[j, 'w2_upfield'], 
                                          rois[j, 'w2_downfield']), return.seq=TRUE)$w2
      else{
        if (i == 1)
          shifts <- matchShift(inFile, w1=c(rois[j, 'w1_upfield'], 
                                            rois[j, 'w1_downfield']), return.seq=TRUE)$w1
        else
          shifts <- matchShift(inFile, w2=c(rois[j, 'w2_upfield'], 
                                            rois[j, 'w2_downfield']), return.seq=TRUE)$w2
      }
      
      ## number of points in this dimension for current block
      np <- length(shifts)
      writeBin(as.integer(np), writeCon, size=4, endian='big')
      
      ## downfield chemical shift in this dimension for current block
      writeBin(as.numeric(shifts[np]), writeCon, size=4, endian='big')
      
      ## upfield chemical shift in this dimension for current block
      writeBin(as.numeric(shifts[1]), writeCon, size=4, endian='big')
      
      ## make header size equal to 30 bytes by filling remainder with 0s
      writeBin(as.integer(rep(0, 14)), writeCon, size=1, endian='big')
    }
  }
  
  ## Write data to file one block at a time
  for (i in 1:nrow(rois)){
    
    if (fp$number_dimensions == 1){
      bounds <- matchShift(inFile, w2=c(rois[i, 'w2_upfield'], 
                                        rois[i, 'w2_downfield']), return.inc=TRUE, invert=TRUE)
      writeBin(as.numeric(inFile$data[bounds$w2[1]:bounds$w2[2]]), writeCon, 
               size=4, endian='big')
    }else{
      bounds <- matchShift(inFile, w1=c(rois[i, 'w1_upfield'], 
                                        rois[i, 'w1_downfield']), w2=c(rois[i, 'w2_upfield'], 
                                                                       rois[i, 'w2_downfield']), return.inc=TRUE, invert=TRUE)
      writeBin(as.numeric(inFile$data[bounds$w2[1]:bounds$w2[2], 
                                      bounds$w1[1]:bounds$w1[2]]), writeCon, size=4, endian='big')
    }
  }
  close(writeCon)
  
  return(outPath)
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
  log_message(msg)
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
    #dfMascot <- read.csv( sMascotFileName, head = TRUE, stringsAsFactors = FALSE, skip = 0)
    dfMascot <- read.csv( sMascotFileName, head = TRUE, stringsAsFactors = FALSE, skip = 3)
  }else
    #dfMascot <- read.csv( sMascotFileName, head = TRUE, stringsAsFactors = FALSE, skip = 0)
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
  flush.console()
  print(paste('Parsing mascot file: ', fileName))

  resList <- parseMascot(fileName)
  
  if(is.null(resList))
  {
    #print('No records found.')
     print('Check file format, pep_seq column not found.')
    return(NULL)
  }else
  {
    print("Align and score peptides.")
    return(alignAndScorePeptides(resList, species))
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
  log_message(msg)
  flush.console()	
  lMatches <- findAlignments(dfPeptides, lSpecies$seq)
  ## map matches
  if(!is.null(lMatches)) 
  {
    msg <- paste0('Generating coincidence map.')
    log_message(msg)
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
  log_message(msg)
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
  
  if (is.null(lMatches) & is.null(geneMap))
  {
    return(-2)
  }


  if (exists('globalSettings')) {}
  else
  {
    globalSettings <- new.env()
    globalSettings$vectorType <- list('EndPoints' = 'EndPoints', 'PepPoints' = 'PepPoints')
  }

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
  log_message('Done')
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
                 log_message("Inappropriate number of files selected for addition operation")
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
                 log_message("Inappropriate number of files selected for subtraction operation")
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
                 log_message("Inappropriate number of files selected for multiplication operation")
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
                 log_message("Inappropriate number of files selected for division operation")
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
  log_message(paste0(sName, ' saved at ', format(Sys.time(), "%H:%M"), '.'))	
  
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
    log_message(paste0('Retrieving peptides from: ', par$file.name, '.'))
    flush.console()		
    if(i == 1)
    {
      pepListHolder <- retrievePeptidesList(par$file.name, par$endian, par$match_list_location, species)
      
    }else
    {
      pepListHolder <- append(pepListHolder, retrievePeptidesList(par$file.name, par$endian, par$match_list_location, species))
    }
  }
  
  log_message('Compiling master peptide list.')
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
  log_message(msg)
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
      log_message('fileFolder does not exist')
      
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
# Old Prepare species files
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
#' Venn Diagram Generator for Peptides
#'
#' This function provides a graphical user interface (GUI) for generating Venn diagrams 
#' from peptide sequences stored in CSV files. Users can select multiple CSV files from 
#' a directory, and the Venn diagram will display common and unique peptides across the 
#' selected files. The resulting Venn diagram can be saved as a PDF file.
#'
#' @return A graphical user interface that allows users to upload CSV files, select files 
#' for analysis, and generate a Venn diagram of peptide sequences.
#'
#' @details
#' The function offers the following main features:
#' \itemize{
#'   \item \strong{File Upload:} Users can select a directory containing CSV files and 
#'   choose specific files for Venn diagram generation.
#'   \item \strong{Venn Diagram Creation:} The function generates a Venn diagram based 
#'   on the unique peptide sequences found in each selected CSV file. The diagram displays 
#'   intersections between groups of peptides.
#'   \item \strong{Saving Output:} The resulting Venn diagram is saved as a PDF file. 
#'   The user can choose the file name and location for the saved output.
#' }
#'
#' @section CSV File Format:
#' Each CSV file must contain a column named \code{pep_seq}, which includes the peptide 
#' sequences. Only the unique peptide sequences in each file will be used for generating 
#' the Venn diagram.
#'
#' @section Instructions:
#' \itemize{
#'   \item Upload CSV files from a directory using the "Browse" button.
#'   \item If no files are selected, the function will use all the CSV files in the chosen directory.
#'   \item After selecting the files, click the "Generate Venn Diagram" button to create the plot.
#'   \item The Venn diagram will be saved as a PDF at the location specified by the user.
#' }
#'
#' @section Dependencies:
#' The function depends on the following R packages:
#' \itemize{
#'   \item \code{tcltk}: For creating the graphical user interface.
#'   \item \code{VennDiagram}: For generating Venn diagrams.
#'   \item \code{readr}: For reading CSV files.
#'   \item \code{grid}: For displaying the generated Venn diagram.
#'   \item \code{tkrplot}: For interactive plotting.
#' }
#'
#' @examples
#' \dontrun{
#' venn_diagram_generator()
#' }
#'
#' @export
# Function to launch the Tcl/Tk Venn Diagram Generator with file selection logic
venn_diagram_generator <- function() {
  
  # Load required libraries
  library(tcltk)
  library(VennDiagram)
  library(readr)
  library(grid)
  library(tkrplot)
  
  # Internal variables for the function
  csv_dir <- tclVar("")  # Initialize the directory variable
  selected_files <- character()  # Initialize a vector to hold selected files
  file_list <- tclVar("")  # Initialize file list variable for the listbox
  
  venn.plot <- NULL  # Store the generated Venn plot
  
  # Function to open the file dialog for selecting the CSV directory
  selectCSVDir1 <- function() {
    csv_dir_value <- tclvalue(tcl("tk_chooseDirectory"))  # Get the selected directory
    tclvalue(csv_dir) <- csv_dir_value  # Update the tclVar with the selected directory
    tkconfigure(csvDirEntry, text = csv_dir_value)  # Update the text of the CSV directory entry
    updateFileList(csv_dir_value)  # Update the file list in the listbox
  }
  
  # Function to update the file list in the listbox
  updateFileList <- function(dir) {
    csv_files <- list.files(path = dir, pattern = "*.csv", full.names = FALSE)  # List CSV files in the directory
    if (length(csv_files) == 0) {
      csv_files <- "No CSV files found"  # Show message if no files found
    }
    tclvalue(file_list) <- csv_files  # Update the file list variable
    tkdelete(listbox, 0, "end")  # Clear the listbox
    for (file in csv_files) {
      tkinsert(listbox, "end", file)  # Insert files into the listbox
    }
  }
  
  # Function to handle multiple file selection from the listbox
  selectFile <- function() {
    selected_indices <- tkcurselection(listbox)  # Get the selected indices from the listbox
    if (length(selected_indices) == 0) {
      selected_files <<- character()  # No files selected, clear the selection
    } else {
      selected_indices <- as.integer(selected_indices)
      selected_files <<- sapply(selected_indices, function(index) {
        tclvalue(tkget(listbox, index))
      })  # Store the selected files
      selected_files <<- file.path(tclvalue(csv_dir), selected_files)  # Get full paths for the selected files
    }
  }
  
  # Function to generate the Venn diagram and save it as a PDF
  generate_venn <- function() {
    if (length(selected_files) == 0) {
      # Use all files if no files are selected
      selected_files <<- list.files(path = tclvalue(csv_dir), pattern = "*.csv", full.names = TRUE)
    }
    
    if (length(selected_files) == 0) {
      tkmessageBox(title = "Error", message = "No files available for plotting!", icon = "error")
      return()
    }
    
    peptide_sets <- lapply(selected_files, function(file) {
      # Suppress messages and hide column types
      data <- suppressMessages(read_csv(file, skip = 3, show_col_types = FALSE))
      if (!"pep_seq" %in% colnames(data)) {
        tkmessageBox(title = "Error", message = paste("Column 'pep_seq' not found in", file), icon = "error")
        return(NULL)
      }
      unique(data[["pep_seq"]])
    })
    
    if (any(sapply(peptide_sets, is.null))) {
      return()  # Return if any file was invalid
    }
    
    category_names <- basename(selected_files)
    
    # Allow the user to choose the file name and location for saving the PDF
    saveFile <- tclvalue(tkgetSaveFile(defaultextension = ".pdf", filetypes = "{{PDF Files} {.pdf}}", initialfile = "venn_diagram.pdf"))
    
    if (saveFile != "") {
      pdf(file = saveFile, width = 8, height = 8)  # Create the PDF file
      venn.plot <<- venn.diagram(
        x = peptide_sets,
        category.names = category_names,
        fill = c("red", "blue", "green", "yellow")[1:length(selected_files)],
        alpha = 0.5,  # Keep transparency
        cex = 1.5,
        cat.cex = 1.5,
        cat.col = c("red", "blue", "green", "yellow")[1:length(selected_files)],
        filename = NULL
      )
      grid.draw(venn.plot)
      dev.off()  # Close the PDF device
      tkmessageBox(title = "Success", message = "Venn diagram saved as PDF")
      
      # Print custom message to the console
      print("Venn Diagram generated and saved")
    }
  }
     
  # Main Tk window
  win <- tktoplevel()
  tkwm.title(win, "Venn Diagram Generator for Peptides")
  
  # Create and configure the CSV directory selection frame with a border and label
  fileSelectionFrame <- ttklabelframe(win, text = "File Selection", padding = 10)
  tkgrid(fileSelectionFrame, padx = 20, pady = 10)
  
  # CSV Directory Selection section inside the bordered frame
  csvDirLabel <- ttklabel(fileSelectionFrame, text = "Choose CSV Directory:", font = "Helvetica 12 bold")
  csvDirEntry <- ttkentry(fileSelectionFrame, textvariable = csv_dir, width = 50)  # Increased the width to 60
  csvDirButton <- ttkbutton(fileSelectionFrame, text = "Browse", command = selectCSVDir1)
  
  tkgrid(csvDirLabel, row = 0, column = 0, sticky = "w", pady = 5)
  tkgrid(csvDirEntry, row = 1, column = 0, padx = 10)
  tkgrid(csvDirButton, row = 1, column = 1, padx = 5)
  
  # Listbox Frame setup inside a labelframe with a label
  listboxFrame <- ttklabelframe(win, text = "Loaded files", padding = 10)
  tkgrid(listboxFrame, padx = 20, pady = 10)
  
  # Adding scrollbars to the listbox
  yscroll <- tkscrollbar(listboxFrame, orient = "vertical", command = function(...) tkyview(listbox, ...))
  xscroll <- tkscrollbar(listboxFrame, orient = "horizontal", command = function(...) tkxview(listbox, ...))
  
  # Configure the listbox for file selection
  listbox <- tklistbox(listboxFrame, height = 10, width = 40, selectmode = "multiple", 
                       yscrollcommand = function(...) tkset(yscroll, ...), 
                       xscrollcommand = function(...) tkset(xscroll, ...), 
                       font = "Helvetica 12 bold",  # Increased font size and made bold
                       bg = "white", fg = "black")
  
  # Grid the listbox and scrollbars inside the labeled frame
  tkgrid(listbox, row = 0, column = 0)
  tkgrid(yscroll, row = 0, column = 1, sticky = "ns")
  tkgrid(xscroll, row = 1, column = 0, sticky = "ew")
  
  # Bind the listbox to the selectFile function to handle file selection
  tkbind(listbox, "<<ListboxSelect>>", selectFile)
  
  # Generate Venn Diagram button styled similarly to the browse button
  generateButton <- ttkbutton(win, text = "Generate Venn Diagram", command = generate_venn)
  tkgrid(generateButton, pady = 10)
  
  #tkmessageBox(title = "Instructions", message = "Upload at least two CSV files with a 'pep_seq' column.\n\nFor plotting purposes, if no files are selected, all files in the directory will be plotted.")
  response <- tclvalue(tkmessageBox(
    title = "Instructions", 
    message = paste(message = "Upload at least two CSV files with a 'pep_seq' column.\n\nFor plotting purposes, if no files are selected, all files in the directory will be plotted.")
  ))
  
  # Check the response
  if (response == "ok") {
    print("Roger, Roger")
  } else {
    print("Other action taken.")
  }
}

# Run the updated Venn Diagram Generator
# venn_diagram_generator()

################################################################################
## Change code so that it uses the current directory to export the file to.

#' Find and Save Unique Peptides from a given mascot file
#'
#' This function prompts the user to select directories and files, reads input CSV files,
#' and identifies unique peptides present in one file but not in another. It then saves
#' these unique peptides to a new CSV file.
#'
#' @return A character string indicating the path of the destination CSV file
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item The user is prompted to select a directory to work in.
#'   \item The user is prompted to select an input CSV file yo use as a Query (L1).
#'   \item The user is prompted to select another input CSV file as the experimental group (L2).
#'   \item The CSV files are read with headers assumed to be present.
#'   \item The function checks if the 'pep_seq' column exists in both L1 and L2 files.
#'     If not, an error is raised.
#'   \item Unique peptides are identified in L2 that are not present in L1.
#'   \item A new filename is created for the output CSV file by appending "_unique" to
#'     the L2 file name (without extension).
#'   \item The unique peptides are saved to the destination CSV file.
#' }
#'
#' @examples
#' \dontrun{
#'   # Call the function
#'   uniquePeptides()
#' }
#'
#' @importFrom utils read.csv write.csv
#' @export
unique_peptides <- function() {
  ## creates main window
  tclCheck()
  dlg <- tktoplevel()
  if (is.null(dlg))
    return(invisible())
  tkwm.title(dlg, 'Unique peptides')
  tkwm.geometry(dlg, "300x105")  # Set the window size
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
# unique_peptides()

################################################################################
## Add the mean and peptide numbers to be displayed on the plot
	 
#' Create Density Plot for Peptide Length Distribution
#'
#' This function generates density plots for the distribution of peptide lengths across
#' different groups based on input CSV files. The function provides options to create
#' different types of plots including overlay density plot, ridges plot, and colored
#' ridges plot.
#'
#' @param csv_dir A character string specifying the directory where the CSV files are located.
#' @param plot_type A character string specifying the type of density plot to generate.
#'                  Options are "overlay", "ridges", or "colored_ridges".
#'
#' @return A ggplot2 object representing the generated density plot.
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item The CSV files in the specified directory are grouped based on a second string
#'     in their filenames.
#'   \item The function calculates the number of characters in the 'pep_seq' column for
#'     each group and calculates the mean length.
#'   \item Depending on the \code{plot_type} argument, the function creates an overlay
#'     density plot, a ridges plot, or a colored ridges plot using the \pkg{ggplot2}
#'     and \pkg{ggridges} packages.
#'   \item The mean length is annotated on the plot.
#' }
#'
#' @examples
#' \dontrun{
#'   # Generate an overlay density plot
#'   overlay_plot <- create_density_plot(csv_dir = "path/to/csv/files", plot_type = "overlay")
#'
#'   # Generate a ridges plot
#'   ridges_plot <- create_density_plot(csv_dir = "path/to/csv/files", plot_type = "ridges")
#'
#'   # Generate a colored ridges plot
#'   colored_ridges_plot <- create_density_plot(csv_dir = "path/to/csv/files", plot_type = "colored_ridges")
#' }
#'
#' @importFrom ggplot2 ggplot geom_density labs theme_bw geom_vline geom_text scale_fill_viridis_c annotate
#' @importFrom ggridges geom_density_ridges2 geom_density_ridges_gradient
#' @importFrom dplyr group_by summarise
#' @export
# peptides_distribution <- function() {
#   # Define variables
#   csv_dir <- ""
#   plot_type <- ""
  
#   # Function to select the CSV directory
#   selectCSVDir1 <- function() {
#     csv_dir <<- tclvalue(tcl("tk_chooseDirectory"))
#     csv_dir <- tclVar(csv_dir)  # Update the csv_dir variable
#     tkconfigure(csvDirEntry, text = csv_dir)  # Update the csvDirEntry widget
#   }
  
#   # Function to set the selected plot type
#   setPlotType <- function(value) {
#     plot_type <<- value
#   }
  
#   group_csv_files <- function(csv_dir) {
#     # Get the list of CSV files in the directory
#     csv_files <- list.files(path = csv_dir, pattern = "*.csv", full.names = TRUE)
    
#     # Create an empty list to store the grouped dataframes
#     grouped_df <- list()
    
#     # Loop through the CSV files and group them based on the second string
#     for (i in 1:length(csv_files)) {
#       # Get the filename without the directory path
#       filename <- basename(csv_files[i])
#       # Get the second string in the filename by splitting on "_"
#       second_str <- strsplit(filename, "_")[[1]][2]
#       # Read the CSV file and extract the pep_seq column
#       csv_data <- read.csv(csv_files[i], skip = 0, header = TRUE)
#       # Check if the 'pep_seq' column exists
#       if (!("pep_seq" %in% colnames(csv_data))) {
#         stop("Error: 'pep_seq' column not found in the CSV file.")
#       }
#       pep_seq_col <- csv_data$pep_seq
#       # Create a new dataframe with the pep_seq column and the group name
#       pep_seq_df <- data.frame(group_name = second_str, pep_seq = pep_seq_col)
#       # If the group already exists in the list, append the dataframe
#       if (second_str %in% names(grouped_df)) {
#         grouped_df[[second_str]] <- rbind(grouped_df[[second_str]], pep_seq_df)
#       }
#       # If the group doesn't exist in the list, create a new list element
#       else {
#         grouped_df[[second_str]] <- pep_seq_df
#       }
#     }
    
#     # Return the list of grouped dataframes
#     return(grouped_df)
#   }
  
#   # Function to create the density plot
#   createDensityPlot <- function() {
#     grouped_df <- group_csv_files(csv_dir)
    
#     # Create an empty data frame to store the nchar values for all groups
#     nchar_df <- data.frame(group_name = character(), nchar = integer())
    
#     for (i in 1:length(grouped_df)) {
#       group_name <- names(grouped_df)[i]
#       pep_seq_col <- grouped_df[[i]]$pep_seq
#       group_nchar_df <- data.frame(group_name = group_name, nchar = nchar(pep_seq_col))
#       nchar_df <- rbind(nchar_df, group_nchar_df)
#     }
    
#     nchar_mean_df <- nchar_df %>%
#       group_by(group_name) %>%
#       summarise(mean_nchar = mean(nchar))
    
#      # Create the appropriate density plot based on the plot_type argument
#   if (plot_type == "overlay") {
#     plot <- ggplot(nchar_df, aes(x = nchar, group = group_name, fill = group_name, color = group_name)) +
#       geom_density(alpha = 0.5) +
#       labs(title = "Peptide Length Distribution", 
#            x = "Peptide length in AA", y = "Density", fill = "Groups") +
#       theme_bw() +
#       geom_vline(data = nchar_mean_df, aes(xintercept = mean_nchar, color = group_name),
#                  size = 1) +
#       labs(color = "Means") +
#       geom_text_repel(data = nchar_mean_df, aes(x = mean_nchar, y = 0, label = round(mean_nchar, 2)), 
#                       color = 'black', size = 4, nudge_x = 5)  
#   } else if (plot_type == "ridges") {
#     plot <- ggplot(nchar_df, aes(x = nchar, y = group_name, fill = group_name)) +
#       geom_density_ridges2(alpha = 0.5, rel_min_height = 0.01, scale = 7, quantile_lines = TRUE, quantile_fun = function(x,...)mean(x)) +
#       labs(title = "Peptide Length Distribution", 
#            x = "Peptide length in AA", y = "", fill = "Group") +
#       theme_bw()  +
#       geom_text(data = nchar_mean_df, aes(x = mean_nchar, y = group_name, label = round(mean_nchar, 2)),
#                 color = 'black', size = 3.5, hjust = 0.5, vjust = - 0.05)
#   } else if (plot_type == "colored_ridges") {
#     plot <- ggplot(nchar_df, aes(x = nchar, y = group_name, fill = stat(x))) +
#       geom_density_ridges_gradient(alpha = 0.5, rel_min_height = 0.02, scale = 3, quantile_lines = TRUE, quantile_fun = function(x,...)mean(x)) +
#       labs(title = "Peptide Length Distribution", 
#            x = "Peptide length in AA", y = "", fill = "Group") +
#       scale_fill_viridis_c(name = "Petide length", option = "H") +
#       geom_text(data = nchar_mean_df, aes(x = mean_nchar, y = group_name, label = round(mean_nchar, 2)),
#                 color = 'black', size = 3.5, hjust = - 0.5, vjust = -1) # Adjust hjust for label positioning
#     } else {
#       stop("Invalid plot_type argument. Please choose either 'Overlay', 'Ridges', or 'Colored Ridges'.")
#     }
    
#     print(plot)
#   }
  
#   # Create the main window
#   win <- tktoplevel()
#   tkwm.title(win, "Density Plot Creator")
  
#   # Create and configure the CSV directory selection frame
#   csvDirFrame <- ttkframe(win)
#   tkgrid(csvDirFrame, padx = 20, pady = 20)
  
#   # Create and configure the CSV directory label
#   csvDirLabel <- ttklabel(csvDirFrame, text = "CSV Directory:")
#   tkgrid(csvDirLabel, column = 1, row = 1, sticky = "w")
  
#   # Create and configure the CSV directory entry
#   csvDirEntry <- ttkentry(csvDirFrame, textvariable = csv_dir)
#   tkgrid(csvDirEntry, column = 2, row = 1)
  
#   # Create and configure the CSV directory browse button
#   csvDirButton <- ttkbutton(csvDirFrame, text = "Browse", command = selectCSVDir1)
#   tkgrid(csvDirButton, column = 3, row = 1, padx = 5)
  
#   # Create and configure the plot type selection frame
#   plotTypeFrame <- ttkframe(win)
#   tkgrid(plotTypeFrame, padx = 10, pady = 10)
  
#   # Create and configure the plot type label
#   plotTypeLabel <- ttklabel(plotTypeFrame, text = "Plot Type:")
#   tkgrid(plotTypeLabel, column = 1, row = 1, sticky = "w")
  
#   # Create and configure the plot type radiobuttons
#   plotTypeRadioFrame <- ttkframe(plotTypeFrame)
#   tkgrid(plotTypeRadioFrame, column = 2, row = 1)
  
#   # Create and configure the Overlay radiobutton
#   overlayRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Overlay", variable = plot_type, value = "overlay", command = function() setPlotType("overlay"))
#   tkgrid(overlayRadio, column = 1, row = 1, sticky = "w")
  
#   # Create and configure the Ridges radiobutton
#   ridgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Ridges", variable = plot_type, value = "ridges", command = function() setPlotType("ridges"))
#   tkgrid(ridgesRadio, column = 1, row = 2, sticky = "w")
  
#   # Create and configure theColored Ridges radiobutton
#   coloredRidgesRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Colored Ridges", variable = plot_type, value = "colored_ridges", command = function() setPlotType("colored_ridges"))
#   tkgrid(coloredRidgesRadio, column = 1, row = 3, sticky = "w")
  
#   # Create and configure the create plot button
#   createPlotButton <- ttkbutton(win, text = "Create Plot", command = createDensityPlot)
#   tkgrid(createPlotButton, pady = 10)
  
#   # Start the event loop
#   tkwait.visibility(win)
# }
# # call the function
# # peptides_distribution()
##########################################################################################################################################

#' Peptide Distribution Plot Generator
#'
#' This function creates a graphical user interface (GUI) for generating peptide length distribution plots 
#' from CSV files. The function allows users to upload multiple CSV files, select a plot type, and visualize 
#' peptide length distributions across different groups based on file naming conventions. The available plot 
#' types include overlay density plots, ridge plots, and colored ridge plots.
#'
#' @return A GUI that lets users upload CSV files, select files for analysis, and generate various peptide 
#' length distribution plots.
#' 
#' @details
#' The function provides the following main features:
#' \itemize{
#'   \item \strong{File Upload:} Users can select a directory containing CSV files and choose specific files for plotting.
#'   \item \strong{Plot Types:} Three types of plots are available: overlay density plots, ridge plots, and colored ridge plots.
#'   \item \strong{Grouping by File Name:} Files are grouped based on the second string in their filenames (e.g., for file \emph{group_sample.csv}, the group is \emph{sample}).
#' }
#'
#' @section CSV File Format:
#' Each CSV file should include a column named \code{pep_seq}, which contains peptide sequences. The function uses this column to calculate peptide lengths for plotting.
#'
#' @section Plot Types:
#' \itemize{
#'   \item \strong{Overlay Plot:} Displays density curves of peptide lengths for each group on the same axes.
#'   \item \strong{Ridge Plot:} Displays density curves of peptide lengths for each group as individual ridges.
#'   \item \strong{Colored Ridge Plot:} Displays ridges with colors representing peptide lengths.
#' }
#'
#' @section Dependencies:
#' This function depends on the following R packages:
#' \itemize{
#'   \item \code{tcltk}: For creating the graphical user interface.
#'   \item \code{ggplot2}, \code{ggridges}: For creating the plots.
#'   \item \code{dplyr}: For data manipulation.
#'   \item \code{ggrepel}: For adding labels to the plot.
#'   \item \code{viridis}: For applying color gradients in colored ridge plots.
#' }
#'
#' @examples
#' \dontrun{
#' peptides_distribution()
#' }
#'
#' @export

peptides_distribution <- function() {
  # Load required libraries
  library(tcltk)
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(ggrepel)
  library(viridis)
  
  # Define variables
  csv_dir <- tclVar("")  # Initialize csv_dir as tclVar
  plot_type <- tclVar("Overlay")  # Initialize plot_type as tclVar with a default value
  selected_files <- character()  # Initialize a vector to hold selected files
  file_list <- tclVar("")  # Initialize file list variable for the listbox
  
  # Function to select the CSV directory
  selectCSVDir1 <- function() {
    csv_dir_value <- tclvalue(tcl("tk_chooseDirectory"))  # Get the selected directory
    tclvalue(csv_dir) <- csv_dir_value  # Update the tclVar with the selected directory
    tkconfigure(csvDirEntry, text = csv_dir_value)  # Update the text of the CSV directory entry
    updateFileList(csv_dir_value)  # Update the file list in the listbox
  }
  
  # Function to update the file list in the listbox
  updateFileList <- function(dir) {
    csv_files <- list.files(path = dir, pattern = "*.csv", full.names = FALSE)  # List CSV files in the directory
    if (length(csv_files) == 0) {
      csv_files <- "No CSV files found"  # Show message if no files found
    }
    tclvalue(file_list) <- csv_files  # Update the file list variable
    tkdelete(listbox, 0, "end")  # Clear the listbox
    for (file in csv_files) {
      tkinsert(listbox, "end", file)  # Insert files into the listbox
    }
  }
  
  # Function to handle multiple file selection from the listbox
  selectFile <- function() {
    selected_indices <- tkcurselection(listbox)  # Get the selected indices from the listbox
    if (length(selected_indices) == 0) {
      selected_files <<- character()  # No files selected
    } else {
      selected_indices <- as.integer(selected_indices)
      selected_files <<- sapply(selected_indices, function(index) {
        tclvalue(tkget(listbox, index))
      })  # Store the selected files
      selected_files <<- file.path(tclvalue(csv_dir), selected_files)  # Get full paths for the selected files
    }
  }
  
  # Function to group the selected CSV file(s)
  group_csv_files <- function() {
    # Use all files if no files are selected
    if (length(selected_files) == 0) {
      selected_files <<- list.files(path = tclvalue(csv_dir), pattern = "*.csv", full.names = TRUE)
    }
    
    if (length(selected_files) == 0) {
      stop("No files available for plotting!")
    }
    
    # Create an empty list to store the grouped dataframes
    grouped_df <- list()
    
    # Loop through the selected files and group them based on the second string
    for (i in 1:length(selected_files)) {
      # Get the filename without the directory path
      filename <- basename(selected_files[i])
      # Get the second string in the filename by splitting on "_"
      second_str <- strsplit(filename, "_")[[1]][2]
      # Read the CSV file and extract the pep_seq column
      csv_data <- read.csv(selected_files[i], skip = 3, header = TRUE)
      
      # Check if the 'pep_seq' column exists
      if (!("pep_seq" %in% colnames(csv_data))) {
        stop("Error: 'pep_seq' column not found in the CSV file.")
      }
      
      # Remove duplicates from the 'pep_seq' column
      csv_data <- csv_data[!duplicated(csv_data$pep_seq), ]
      
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
    grouped_df <- group_csv_files()  # Call the function to group selected files
    
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
    
    plot_type_value <- tclvalue(plot_type)  # Get the selected plot type
    
    # Convert to lowercase for consistency and ensure valid input
    plot_type_value <- tolower(plot_type_value)
    
    # Create the appropriate density plot based on the plot_type argument
    if (plot_type_value == "overlay") {
      plot <- ggplot(nchar_df, aes(x = nchar, group = group_name, fill = group_name, color = group_name)) +
        geom_density(alpha = 0.5) +
        labs(title = "Peptide Length Distribution", 
             x = "Peptide length in AA", y = "Density", fill = "Groups") +
        theme_bw() +
        geom_vline(data = nchar_mean_df, aes(xintercept = mean_nchar, color = group_name),
                   size = 1) +
        labs(color = "Means") +
        geom_text_repel(data = nchar_mean_df, aes(x = mean_nchar, y = 0, label = round(mean_nchar, 2)), 
                        color = 'black', size = 4, nudge_x = 5, box.padding = 0.5, 
                        segment.color = 'grey50')  # Prevents overlap
    } else if (plot_type_value == "ridges") {
      plot <- ggplot(nchar_df, aes(x = nchar, y = group_name, fill = group_name)) +
        geom_density_ridges2(alpha = 0.5, rel_min_height = 0.01, scale = 7, quantile_lines = TRUE, quantile_fun = function(x,...)mean(x)) +
        labs(title = "Peptide Length Distribution", 
             x = "Peptide length in AA", y = "", fill = "Group") +
        theme_bw() +
        geom_text_repel(data = nchar_mean_df, aes(x = mean_nchar, y = group_name, label = round(mean_nchar, 2)),
                        color = 'black', size = 3.5, hjust = 0.5, vjust = -0.05, 
                        box.padding = 0.5, segment.color = 'grey50')  # Prevents overlap
    } else if (plot_type_value == "colored ridges") {
      plot <- ggplot(nchar_df, aes(x = nchar, y = group_name, fill = stat(x))) +
        geom_density_ridges_gradient(alpha = 0.5, rel_min_height = 0.02, scale = 3, quantile_lines = TRUE, quantile_fun = function(x,...)mean(x)) +
        labs(title = "Peptide Length Distribution", 
             x = "Peptide length in AA", y = "", fill = "Group") +
        scale_fill_viridis_c(name = "Peptide length", option = "H") +
        geom_text_repel(data = nchar_mean_df, aes(x = mean_nchar, y = group_name, label = round(mean_nchar, 2)),
                        color = 'black', size = 3.5, hjust = -0.5, vjust = -1, 
                        box.padding = 0.5, segment.color = 'grey50')  # Prevents overlap
    } else {
      stop("Invalid plot_type argument. Please choose either 'Overlay', 'Ridges', or 'Colored Ridges'.")
    }
    
    print(plot)
  }
  
  # Function to handle window closing event
  on_closing <- function() {
    tkdestroy(win)  # Destroy the window to avoid errors
  }
  
  # Create the main window
  win <- tktoplevel()
  tkwm.title(win, "Density Plot Generator")
  
  # Bind the close event handler to the window
  tkwm.protocol(win, "WM_DELETE_WINDOW", on_closing)  # Handle window close event
  
  # Create and configure the CSV directory selection frame with a border and label
  fileSelectionFrame <- ttklabelframe(win, text = "File Selection", padding = 10)
  tkgrid(fileSelectionFrame, padx = 20, pady = 10)
  
  # Create and configure the CSV directory label and entry
  csvDirLabel <- ttklabel(fileSelectionFrame, text = "Choose CSV Directory:", font = "Helvetica 12 bold")
  csvDirEntry <- ttkentry(fileSelectionFrame, textvariable = csv_dir, width = 50)
  csvDirButton <- ttkbutton(fileSelectionFrame, text = "Browse", command = selectCSVDir1)
  
  tkgrid(csvDirLabel, row = 0, column = 0, sticky = "w", pady = 5)
  tkgrid(csvDirEntry, row = 1, column = 0, padx = 10)
  tkgrid(csvDirButton, row = 1, column = 1, padx = 5)
  
  # Listbox Frame setup inside a labelframe with a label
  listboxFrame <- ttklabelframe(win, text = "Loaded files", padding = 10)
  tkgrid(listboxFrame, padx = 20, pady = 10)
  
  # Adding Scrollbars
  yscroll <- tkscrollbar(listboxFrame, orient = "vertical", command = function(...) tkyview(listbox, ...))
  xscroll <- tkscrollbar(listboxFrame, orient = "horizontal", command = function(...) tkxview(listbox, ...))
  
  # Configure the listbox
  listbox <- tklistbox(listboxFrame, height = 10, width = 40, selectmode = "multiple", 
                       yscrollcommand = function(...) tkset(yscroll, ...), 
                       xscrollcommand = function(...) tkset(xscroll, ...), 
                       font = "Helvetica 12 bold", bg = "white", fg = "black")  # Increased font size and made bold
  
  # Grid the listbox and scrollbars
  tkgrid(listbox, row = 0, column = 0)
  tkgrid(yscroll, row = 0, column = 1, sticky = "ns")
  tkgrid(xscroll, row = 1, column = 0, sticky = "ew")
  
  # Bind the listbox to the selectFile function (to handle file selection)
  tkbind(listbox, "<<ListboxSelect>>", selectFile)
  
  # Create and configure the plot type selection frame
  plotTypeFrame <- ttkframe(win)
  tkgrid(plotTypeFrame, padx = 20, pady = 10)
  
  # Create and configure the plot type label and combobox
  plotTypeLabel <- ttklabel(plotTypeFrame, text = "Select Plot Type:", font = "Helvetica 12 bold")
  plotTypeCombo <- ttkcombobox(plotTypeFrame, values = c("Overlay", "Ridges", "Colored Ridges"), textvariable = plot_type, state = "readonly")
  
  tkgrid(plotTypeLabel, row = 0, column = 0, pady = 5, sticky = "w")
  tkgrid(plotTypeCombo, row = 1, column = 0, padx = 10)
  
  # Create and configure the create plot button
  createPlotButton <- ttkbutton(win, text = "Generate Plot", command = createDensityPlot)
  tkgrid(createPlotButton, pady = 10)

  # Create and configure the instruction box
  response <- tclvalue(tkmessageBox(title = "Instructions", message = "Upload CSV files with a 'pep_seq' column.\n\nFor plotting purposes, if no files are selected, all files in the directory will be plotted."))
    # Check the response
  if (response == "ok") {
    print("Roger, Roger")
  } else {
    print("Other action taken.")
  }
}

# Call the function
# peptides_distribution()

################################################################################
# Old Cleavage site distribution function
# Updated with a seq logo distribution function
################################################################################

#'  Cleavage Sites Distribution 
#'
#' This function generates various plots to visualize cleavage site distributions from
#' input CSV files. The plots include percentages of mean first and last letter counts
#' for different groups of data.
#'
#' @param csv_dir A character string specifying the directory containing the input CSV files.
#' @param plot_type A character string specifying the type of plot to generate. Options are
#'                  "Nter" (N-terminal cleavage sites), "Cter" (C-terminal cleavage sites),
#'                  or "Combined" (combined plots for mean first and last letter counts).
#'
#' @return A ggplot2 object representing the generated plot.
#'
#' @details The function performs the following steps:
#' \itemize{
#'   \item The CSV files in the specified directory are grouped based on a second string
#'     in their filenames.
#'   \item First and last letters are extracted from each sequence in the CSV files.
#'   \item Percentages of mean first and last letter counts are calculated for each group.
#'   \item The function creates plots based on the specified plot type: "Nter" (N-terminal),
#'     "Cter" (C-terminal), or "Combined" (both mean first and last letter counts).
#' }
#'
#' @examples
#' \dontrun{
#'   # Generate a plot for N-terminal cleavage sites
#'   nter_plot <- plotCutSites(csv_dir = "path/to/csv/files", plot_type = "Nter")
#'
#'   # Generate a plot for C-terminal cleavage sites
#'   cter_plot <- plotCutSites(csv_dir = "path/to/csv/files", plot_type = "Cter")
#'
#'   # Generate a combined plot for mean first and last letter counts
#'   combined_plot <- plotCutSites(csv_dir = "path/to/csv/files", plot_type = "Combined")
#' }
#'
#' @importFrom dplyr group_by summarise
#' @importFrom ggplot2 ggplot geom_bar geom_errorbar ggtitle xlab ylab theme_bw
#' @export

# cut_sites_distribution <- function() {
#   # Load required libraries
#   library(dplyr)
#   library(ggplot2)
#   library(tcltk)

#   # Define variables
#   csv_directory <- ""
#   PlotType <- ""
#   protein_seq <- ""
  
#   # Function to select the CSV directory
#   selectCSVDir <- function() {
#     csv_directory <<- tclvalue(tcl("tk_chooseDirectory"))
#     csv_directory <- tclVar(csv_directory)  # Update the csv_directory variable
#     tkconfigure(csvDirEntry, text = csv_directory)  # Update the csvDirEntry widget
#   }
  
#   # Function to set the selected plot type
#   setPlotType <- function(value) {
#     PlotType <<- value
#   }

  
#   calculateLetterPercentages <- function(protein_seq) {

#     #protein_sequence <- tclvalue(protein_seq)
#     protein_seq <- as.character(protein_seq)
#     print(protein_seq)
#     letter_percentages <<- table(strsplit(protein_seq, "")[[1]]) / nchar(protein_seq) * 100
#   }
  
  
#   # Function to create the bar plot
#   createBarPlot <- function(protein_seq, PlotType) {
    
#     library(dplyr)
#     library(ggplot2)

#     valid_plot_types <- c("Nter", "Cter", "Cter_normalized", "Nter_normalized", "Combined_normalized")
#     if (!(PlotType %in% valid_plot_types)) {
#       stop("Error: Invalid PlotType argument. Valid options are 'Nter', 'Cter', 'Cter_normalized', 'Nter_normalized' and 'Combined_normalized'.")
#     }
    
#     extract_first_last_letters_from_csv <- function(file_path) {
#       # Read the CSV file skipping the first 3 lines (assuming headers are in line 4)
#       data <- read.csv(file_path, skip = 0, header = TRUE)
      
#       # Check if the 'pep_seq' column exists
#       if (!("pep_seq" %in% colnames(data))) {
#         stop("Error: 'pep_seq' column not found in the CSV file.")
#       }
      
#       # Extract the first and last letters from each sequence
#       data$first_letter <- substr(data$pep_seq, 1, 1)
#       data$last_letter <- substr(data$pep_seq, nchar(data$pep_seq), nchar(data$pep_seq))
      
#       # Count the occurrences of each first letter
#       first_letter_counts <- table(factor(data$first_letter, levels = c("A", "C", "D", "E", "F", "G","H", "I" ,"K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")))
      
#       # Count the occurrences of each last letter
#       last_letter_counts <- table(factor(data$last_letter, levels = c("A", "C", "D", "E", "F", "G","H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")))
      
#       # Combine the counts for first and last letters
#       combined_counts <- merge(
#         data.frame(letter = names(first_letter_counts), first_letter_count = as.numeric(first_letter_counts)),
#         data.frame(letter = names(last_letter_counts), last_letter_count = as.numeric(last_letter_counts)),
#         by = "letter", all = TRUE
#       )
#       combined_counts$combined_count <- rowSums(combined_counts[, c("first_letter_count", "last_letter_count")], na.rm = TRUE)
      
#       # Replace NA values with 0
#       combined_counts[is.na(combined_counts)] <- 0
      
#       # Return the combined counts
#       return(combined_counts)
#     }
    
#     group_csv_files <- function(csv_directory) {
#       # Get the list of CSV files in the directory
#       csv_files <- list.files(path = csv_directory, pattern = "*.csv", full.names = TRUE)
      
#       # Create an empty list to store the grouped CSV files
#       grouped_csv <- list()
      
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
      
#       # Return the list of grouped CSV files
#       return(grouped_csv)
#     }
    
#     # Load the necessary libraries
#     library(dplyr)
#     library(ggplot2)
    
#     # Group the CSV files
#     csv <- group_csv_files(csv_directory)
    
#     # Initialize an empty dataframe to store the results
#     dat <- data.frame()
    
#     # Iterate over each group of CSV files
#     for (groups in csv) {
      
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
#         Letter_Percentage = double(),
#         stringsAsFactors = FALSE
#       )
      
#       # Iterate over each CSV file in the group
#       for (i in groups) {
#         print(i)
        
#         # Extract the first and last letters from the CSV file
#         file_results <- extract_first_last_letters_from_csv(i)
        
#         # Get the filename without the directory path
#         filename <- basename(i)
        
#         # Get the second string in the filename by splitting on "_"
#         grp_str <- rep(strsplit(filename, "_")[[1]][2], nrow(file_results))
        
#         # Store the label for the group
#         file_results$group <- data.frame(grp_str)
        
#         # Append the results to the dataframe
#         dat <- rbind(dat, file_results)
#       }
#     }
    
    
#     # Add a call to calculateLetterPercentages before the dplyr::summarize function
#     calculateLetterPercentages(protein_seq)
    
#     # Then, use letter_percentages in your dplyr::summarize function
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
#         sum_last_letter_count = sum(last_letter_count),
#         mean = mean(combined_count),
#         sd = sd(combined_count),
#         se = sd / sqrt(n()),
#         sum = sum(combined_count),
#         Letter_Percentage = letter_percentages[letter]
#       )
    
    
#     # Calculate percentage of mean_Last_letter_count and mean_first_letter_count with SE
#     agg_df_sum <- agg_df %>%
#       group_by(group) %>%
#       dplyr::summarize(
#         total_sum_combined_count = sum(sum_combined_count),
#         total_sum_first_letter_count = sum(sum_first_letter_count),
#         total_sum_last_letter_count = sum(sum_last_letter_count),
#         total_sum = sum(sum)
#       )
    
#     agg_df <- agg_df %>%
#       left_join(agg_df_sum, by = "group") %>%
#       dplyr::mutate(
#         perc_mean_first_letter_count = mean_first_letter_count / total_sum_first_letter_count * 100,
#         perc_mean_last_letter_count = mean_last_letter_count / total_sum_last_letter_count * 100,
#         SE_perc_mean_first_letter_count = sqrt(perc_mean_first_letter_count / 100 * (100 - perc_mean_first_letter_count) / total_sum_first_letter_count),
#         SE_perc_mean_last_letter_count = sqrt(perc_mean_last_letter_count / 100 * (100 - perc_mean_last_letter_count) / total_sum_last_letter_count),
#         percentage = sum / total_sum * 100,
#         SE_percentage = sqrt(percentage/100 * (100 - percentage)/sum),
#         Normalized_Percentage_combined = percentage / Letter_Percentage,
#         Normalized_SE_combined = sqrt(Normalized_Percentage_combined/100 * (100 - Normalized_Percentage_combined)/sum),
#         Normalized_Percentage_Nter = perc_mean_first_letter_count / Letter_Percentage,
#         Normalized_SE_Nter = sqrt(Normalized_Percentage_Nter/100 * (100 - Normalized_Percentage_Nter)/sum),
#         Normalized_Percentage_Cter = perc_mean_last_letter_count / Letter_Percentage,
#         Normalized_SE_Cter = sqrt(Normalized_Percentage_Cter/100 * (100 - Normalized_Percentage_Cter)/sum)
#       )
    
#     # Convert the grouped dataframe to a regular dataframe
#     agg_df <- data.frame(agg_df)
    
#     # Convert the group variable to a factor
#     grp <- agg_df$group[[1]]
#     grp <- as.factor(grp)
#     agg_df$group <- grp
    
    
#     # Create the appropriate density plot based on the PlotType argument
    
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
    
#     # Cter Normalized Plot
#     else if (PlotType == "Cter_normalized") {
#       cutsiteplots <- ggplot(agg_df, aes(x = letter, y = Normalized_Percentage_Cter, fill = group)) +
#         geom_bar(position = "dodge", stat = "identity") +
#         geom_errorbar(
#           aes(ymin = Normalized_Percentage_Cter - Normalized_SE_Cter,
#               ymax = Normalized_Percentage_Cter + Normalized_SE_Cter),
#           position = position_dodge(width = 0.9), width = 0.25
#         ) +
#         labs(title = "C-terminal cleavage sites normalized",
#              x = "", y = "Percentage normalized to AA distribution in protein sequence") +
#         theme_grey(base_size = 12) +
#         theme(
#           axis.text.x = element_text(size = 14, face = "bold"),
#           axis.text.y = element_text(size = 12, face = "bold"),
#           title = element_text(size = 14, face = "bold")
#         )
#     }
    
#     # nter Normalized Plot
#     else if (PlotType == "Nter_normalized") {
#       cutsiteplots <- ggplot(agg_df, aes(x = letter, y = Normalized_Percentage_Nter, fill = group)) +
#         geom_bar(position = "dodge", stat = "identity") +
#         geom_errorbar(
#           aes(ymin = Normalized_Percentage_Nter - Normalized_SE_Nter,
#               ymax = Normalized_Percentage_Nter + Normalized_SE_Nter),
#           position = position_dodge(width = 0.9), width = 0.25
#         ) +
#         labs(title = "N-terminal cleavage sites normalized",
#              x = "", y = "Percentage normalized to AA distribution in protein sequence") +
#         theme_grey(base_size = 12) +
#         theme(
#           axis.text.x = element_text(size = 14, face = "bold"),
#           axis.text.y = element_text(size = 12, face = "bold"),
#           title = element_text(size = 14, face = "bold")
#         )
#     }
    
#     # Combined Normalized Plot
#     else if (PlotType == "Combined_normalized") {
#       cutsiteplots <- ggplot(agg_df, aes(x = letter, y = Normalized_Percentage_combined, fill = group)) +
#         geom_bar(position = "dodge",stat = "identity") +
#         geom_errorbar(aes(ymin = Normalized_Percentage_combined - Normalized_SE_combined,
#                           ymax = Normalized_Percentage_combined + Normalized_SE_combined),
#                       position = position_dodge(width = 0.9), width = 0.25) +
#         labs(title = "Combined Normalized",
#              x = "", y = "Percentage normalized to AA distribution in protein sequence") +
#         theme_grey(base_size = 12) +
#         theme(axis.text.x = element_text(size = 14, face = "bold"),
#               axis.text.y = element_text(size = 12, face = "bold"),
#               title = element_text( face = "bold")
#         )
#     } 
    
#     # If the PlotType argument is not valid
#     else {
#       stop("Error: Invalid PlotType argument. Valid options are 'Nter', 'Cter', 'Cter_normalized', 'Nter_normalized' and 'Combined_normalized'.")
#     }
    
#     # Return the plot
#     print(cutsiteplots)
#   }
  
#   # Create the main window
#   win <- tktoplevel()
#   tkwm.title(win, "Cut sites distribution")
  
#   # Create and configure the CSV directory selection frame
#   csvDirFrame <- ttkframe(win)
#   tkgrid(csvDirFrame, padx = 20, pady = 20)
  
#   # Create and configure the CSV directory label
#   csvDirLabel <- ttklabel(csvDirFrame, text = "CSV Directory:")
#   tkgrid(csvDirLabel, column = 1, row = 1, sticky = "e")  # Center text to the right
  
#   # Create and configure the CSV directory entry
#   csvDirEntry <- ttkentry(csvDirFrame, textvariable = csv_directory)
#   tkgrid(csvDirEntry, column = 2, row = 1, sticky = "we")  # Center input box within column
  
#   # Create and configure the CSV directory browse button
#   csvDirButton <- ttkbutton(csvDirFrame, text = "Browse", command = selectCSVDir)
#   tkgrid(csvDirButton, column = 3, row = 1, padx = 5, sticky = "w")  # Center button to the left
  
#   # Create a frame for the protein sequence label, entry widget, and import sequence button
#   seqFrame <- ttkframe(win)
#   tkgrid(seqFrame, padx = 20, pady = 10, columnspan = 3, row = 2)
  
#   # Create and configure the protein sequence label
#   seqLabel <- ttklabel(seqFrame, text = "Protein sequence:")
#   tkgrid(seqLabel, column = 1, row = 1, sticky = "e")  # Center text to the right
  
#   # Create and configure the sequence entry widget
#   seqEntry <- ttkentry(seqFrame)
#   tkgrid(seqEntry, column = 2, row = 1, sticky = "we")  # Center input box within column
  
#   # Create and configure the calculate button
#   calculateButton <- ttkbutton(seqFrame, text = "Import sequence", command = function(){
#     input_text <-tkget(seqEntry)
#     calculateLetterPercentages(input_text)
#   })
#   tkgrid(calculateButton, column = 3, row = 1, padx = 5, sticky = "w")  # Center button to the left
  
#   # Create and configure the plot type selection frame
#   plotTypeFrame <- ttkframe(win)
#   tkgrid(plotTypeFrame, padx = 10, pady = 10)
  
#   # Create and configure the plot type label
#   plotTypeLabel <- ttklabel(plotTypeFrame, text = "Plot Type:")
#   tkgrid(plotTypeLabel, column = 1, row = 1, sticky = "w")
  
#   # Create and configure the plot type radiobuttons
#   plotTypeRadioFrame <- ttkframe(plotTypeFrame)
#   tkgrid(plotTypeRadioFrame, column = 2, row = 1)
  
#   # Create and configure the Cter radiobutton
#   CterRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Cter", variable = PlotType, value = "Cter")
#   tkgrid(CterRadio, column = 1, row = 1, sticky = "w")
  
#   # Create and configure the Nter radiobutton
#   NterRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Nter", variable = PlotType, value = "Nter")
#   tkgrid(NterRadio, column = 1, row = 2, sticky = "w")
  
#   # Create and configure the Cter normalized radiobutton
#   CternormalRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Cter normalized", variable = PlotType, value = "Cter_normalized")
#   tkgrid(CternormalRadio, column = 1, row = 3, sticky = "w")
  
#   # Create and configure the Nter normalized radiobutton
#   NternormalRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Nter normalized", variable = PlotType, value = "Nter_normalized")
#   tkgrid(NternormalRadio, column = 1, row = 4, sticky = "w")
  
#   # Create and configure the Combined normalized radiobutton
#   CombinedNormalRadio <- ttkradiobutton(plotTypeRadioFrame, text = "Combined normalized", variable = PlotType, value = "Combined_normalized")
#   tkgrid(CombinedNormalRadio, column = 1, row = 5, sticky = "w")
  
#   # Create and configure the create plot button
#   createPlotButton <- ttkbutton(win, text = "Create Plot", command = function(){
#     input_text <-tkget(seqEntry)
#     print(input_text)
#     createBarPlot(input_text, tclvalue(PlotType))
#   }) 
#   tkgrid(createPlotButton, pady = 10)
#   tkwait.visibility(win)
  
# }

# Call the function
#cut_sites_distribution()

################################################################################
# Updated csd function with ggseqlogo
################################################################################

#' Cut Sites Distribution GUI
#'
#' This function creates a graphical user interface (GUI) for generating sequence logo plots based on specific positions within sequences from a CSV file. The GUI allows users to select a CSV file containing sequences and optionally provide a protein sequence for normalization. Users can generate either normalized or non-normalized sequence logo plots.
#'
#' @details The GUI includes the following features:
#' \itemize{
#'   \item A file input section to browse and select a CSV file containing sequences.
#'   \item An input box for entering a protein sequence used for normalization.
#'   \item Buttons to generate normalized and non-normalized sequence logo plots.
#'   \item The function uses the \code{ggseqlogo} package to create sequence logo plots for specified positions within the sequences.
#' }
#'
#' @return This function does not return a value. It opens a GUI window where users can interactively generate sequence logo plots.
#'
#' @importFrom tcltk tktoplevel tkwm.title tkconfigure tkframe tklabel tkgrid tkentry tclVar tkgetOpenFile tkbutton tkmessageBox tkbind
#' @importFrom ggseqlogo ggseqlogo
#' @importFrom dplyr table
#' @export
#'
#' @examples
#' \dontrun{
#'   # Launch the GUI
#'   cut_sites_distribution()
#' }

# cut_sites_distribution <- function() {
#   # Load required libraries
#   library(tcltk)
#   library(ggseqlogo)
#   library(dplyr)
#   library(ggplot2)  # Added for plotting the vertical dashed line
  
#   # Function to generate the GUI
#   generateGUI <- function() {
    
#     # Create main window
#     tt <- tktoplevel()
#     tkwm.title(tt, "Sequence Logo Plot for Specific Positions")
#     tkconfigure(tt, background = "#f5f5f5")
    
#     # Variables to store file path and protein sequence
#     file_var <- tclVar()
#     protein_seq_var <- tclVar("VDPENVFRKLL")
    
#     # Frame for file input
#     file_frame <- tkframe(tt, background = "#f5f5f5")
#     tkgrid(tklabel(file_frame, text = "Choose CSV File", background = "#f5f5f5", font = "Helvetica 10 bold"), pady = 5)
#     file_entry <- tkentry(file_frame, textvariable = file_var, width = 50)
#     tkgrid(file_entry, pady = 5)
    
#     # Define the command to be executed when the "Browse" button is clicked
#     file_button <- tkbutton(file_frame, text = "Browse", command = function() {
#       file_path <- tclvalue(tkgetOpenFile(filetypes = "{{CSV Files} {.csv}} {{All files} *}"))
#       if (file_path != "") {
#         tclvalue(file_var) <- file_path
#       }
#     })
#     tkgrid(file_button, pady = 5)
#     tkgrid(file_frame, padx = 20, pady = 10)
    
#     # Frame for protein sequence input (only for normalization)
#     protein_frame <- tkframe(tt, background = "#f5f5f5")
    
#     # Protein sequence label centered above the input box
#     protein_label <- tklabel(protein_frame, text = "Enter Protein Sequence for Normalization:", 
#                              background = "#f5f5f5", font = "Helvetica 10 bold")
#     tkgrid(protein_label, row = 0, column = 0, columnspan = 2, pady = 5)
    
#     protein_entry <- tkentry(protein_frame, textvariable = protein_seq_var, width = 50)
#     tkgrid(protein_entry, row = 1, column = 0, columnspan = 2, pady = 5)
    
#     # Define the command to display the sequence in the console
#     submit_sequence <- function() {
#       sequence <- tclvalue(protein_seq_var)
#       cat("Protein sequence:", sequence, "\n")
#     }
    
#     # Bind the Enter key to the protein_entry widget
#     tkbind(protein_entry, "<Return>", function() submit_sequence())
    
#     enter_button <- tkbutton(protein_frame, text = "Enter", command = submit_sequence)
#     tkgrid(enter_button, row = 1, column = 2, padx = 5, pady = 5)
    
#     tkgrid(protein_frame, padx = 20, pady = 10)
    
#     # Function to process data and generate the sequence logo plot
#     generatePlot <- function(normalized) {
      
#       file_path <- tclvalue(file_var)
#       protein_seq <- toupper(tclvalue(protein_seq_var))
      
#       # Validate file selection
#       if (file_path == "") {
#         tkmessageBox(message = "Please choose a CSV file.", icon = "error")
#         return()
#       }
      
#       # Read the CSV file
#       data <- tryCatch({
#         read.csv(file_path, stringsAsFactors = FALSE)
#       }, error = function(e) {
#         tkmessageBox(message = "Error reading the CSV file.", icon = "error")
#         return(NULL)
#       })
      
#       if (is.null(data)) return()
      
#       # Ensure 'pep_seq' column is present
#       if (!"pep_seq" %in% colnames(data)) {
#         tkmessageBox(message = "The uploaded file does not contain a 'pep_seq' column.", icon = "error")
#         return()
#       }
      
#       # Extract sequences from the CSV file
#       sequences <- data$pep_seq
      
#       # Adjust position extraction to handle sequences shorter than 4 amino acids
#       sequence_matrix <- sapply(sequences, function(seq) {
#         n <- nchar(seq)  # Get the length of the sequence
#         if (n >= 4) {
#           # Extract the last 4 positions (P4 to P1) and the first 4 positions (P1' to P4')
#           p4_p1 <- c(substr(seq, n - 3, n - 3), 
#                      substr(seq, n - 2, n - 2), 
#                      substr(seq, n - 1, n - 1), 
#                      substr(seq, n, n))
#           p1p_p4p <- c(substr(seq, 1, 1), 
#                        substr(seq, 2, 2), 
#                        substr(seq, 3, 3), 
#                        substr(seq, 4, 4))
#           return(c(p4_p1, p1p_p4p))  # Combine both parts
#         } else {
#           # Handle sequences shorter than 4 AA by filling the missing positions with NA
#           available_part <- c(strsplit(seq, "")[[1]], rep(NA, 8 - n))
#           return(available_part)
#         }
#       }, USE.NAMES = FALSE)
      
#       # Transpose the matrix to have positions as columns
#       sequence_matrix <- t(sequence_matrix)
      
#       if (normalized) {
#         # Calculate frequencies for each letter at each position, ignoring NAs
#         position_freq <- apply(sequence_matrix, 2, function(pos) {
#           table(factor(pos, levels = LETTERS), useNA = "no") / sum(!is.na(pos))
#         })
        
#         # Ensure protein sequence covers the required positions
#         if (nchar(protein_seq) < 10) {
#           tkmessageBox(message = "The protein sequence must be long enough to cover the required positions.", icon = "error")
#           return()
#         }
        
#         # Calculate the overall frequency of each letter in the protein sequence
#         protein_freq <- table(factor(unlist(strsplit(protein_seq, "")), levels = LETTERS)) / nchar(protein_seq)
        
#         # Normalize the position frequencies by the protein sequence frequencies
#         normalized_freq <- sweep(position_freq, 1, protein_freq, "/")
        
#         # Replace any non-finite or negative values with 0
#         normalized_freq[!is.finite(normalized_freq) | normalized_freq <= 0] <- 0
        
#         # Generate the sequence logo plot using ggseqlogo
#         logo_plot <- ggseqlogo(normalized_freq, method = "prob")
#       } else {
#         # For non-normalized sequences
#         selected_positions <- sapply(sequences, function(seq) {
#           n <- nchar(seq)
#           if (n >= 4) {
#             paste0(substr(seq, n - 3, n - 3),  # P4
#                    substr(seq, n - 2, n - 2),  # P3
#                    substr(seq, n - 1, n - 1),  # P2
#                    substr(seq, n, n),          # P1
#                    substr(seq, 1, 1),          # P1'
#                    substr(seq, 2, 2),          # P2'
#                    substr(seq, 3, 3),          # P3'
#                    substr(seq, 4, 4))          # P4'
#           } else {
#             paste0(rep("NA", 8), collapse = "")
#           }
#         }, USE.NAMES = FALSE)
        
#         # Remove any positions that contain "NA"
#         selected_positions <- selected_positions[!grepl("NA", selected_positions)]
        
#         # Create the sequence logo plot using ggseqlogo
#         logo_plot <- ggseqlogo(selected_positions, method = "prob")
#       }
      
#       # Customize x-axis labels to include "cs" and other positions
#       logo_plot <- logo_plot + 
#         scale_x_continuous(breaks = c(1:8, 9), 
#                            labels = c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'", "cs")) +
#         geom_vline(xintercept = 4.5, linetype = "dashed", color = "black", size = 1)  # Cleavage site at "cs"
      
#       # Plot the sequence logo
#       print(logo_plot)
#     }
    
#     # Buttons to generate plots
#     button_frame <- tkframe(tt, background = "#f5f5f5")
#     plot_normalized_button <- tkbutton(button_frame, text = "Generate Normalized Plot", command = function() generatePlot(TRUE))
#     tkgrid(plot_normalized_button, padx = 10, pady = 10)
    
#     plot_non_normalized_button <- tkbutton(button_frame, text = "Generate Non-normalized Plot", command = function() generatePlot(FALSE))
#     tkgrid(plot_non_normalized_button, padx = 10, pady = 10)
    
#     tkgrid(button_frame, padx = 20, pady = 10)
#   }
  
#   # Run the GUI
#   generateGUI()
# }

# # Call the function
# #cut_sites_distribution()

################################################################################
#' Peptide Distribution Plot Generator
#'
#' This function creates a graphical user interface (GUI) for generating peptide length distribution plots 
#' from CSV files. The function allows users to upload multiple CSV files, select a plot type, and visualize 
#' peptide length distributions across different groups based on file naming conventions. The available plot 
#' types include overlay density plots, ridge plots, and colored ridge plots.
#'
#' @return A GUI that lets users upload CSV files, select files for analysis, and generate various peptide 
#' length distribution plots.
#' 
#' @details
#' The function provides the following main features:
#' \itemize{
#'   \item \strong{File Upload:} Users can select a directory containing CSV files and choose specific files for plotting.
#'   \item \strong{Plot Types:} Three types of plots are available: overlay density plots, ridge plots, and colored ridge plots.
#'   \item \strong{Grouping by File Name:} Files are grouped based on the second string in their filenames (e.g., for file \emph{group_sample.csv}, the group is \emph{sample}).
#' }
#'
#' @section CSV File Format:
#' Each CSV file should include a column named \code{pep_seq}, which contains peptide sequences. The function uses this column to calculate peptide lengths for plotting.
#'
#' @section Plot Types:
#' \itemize{
#'   \item \strong{Overlay Plot:} Displays density curves of peptide lengths for each group on the same axes.
#'   \item \strong{Ridge Plot:} Displays density curves of peptide lengths for each group as individual ridges.
#'   \item \strong{Colored Ridge Plot:} Displays ridges with colors representing peptide lengths.
#' }
#'
#' @section Dependencies:
#' This function depends on the following R packages:
#' \itemize{
#'   \item \code{tcltk}: For creating the graphical user interface.
#'   \item \code{ggplot2}, \code{ggridges}: For creating the plots.
#'   \item \code{dplyr}: For data manipulation.
#'   \item \code{ggrepel}: For adding labels to the plot.
#'   \item \code{viridis}: For applying color gradients in colored ridge plots.
#' }
#'
#' @examples
#' \dontrun{
#' peptides_distribution()
#' }
#'
#' @export
cut_sites_distribution <- function() {
  # Load required libraries
  library(tcltk)
  library(ggseqlogo)
  library(dplyr)
  library(ggplot2)
  
  # Function to generate the GUI
  generateGUI <- function() {
    
    # Create main window
    tt <- tktoplevel()
    tkwm.title(tt, "P4-P4' Analysis")
    
    # Variables to store file path, protein sequence, plot type, normalization option, and position for bar plot
    file_var <- tclVar()
    protein_seq_var <- tclVar("VDPENVFRKLL")
    plot_type_var <- tclVar("Logo Plot")  # Default selection for the plot type
    normalization_var <- tclVar("Non-Normalized")  # Default selection for logo plot normalization
    position_var <- tclVar("P1")  # Default selection for bar plot (P1)
    
    # Frame for file input and protein sequence input (above plot type selection)
    input_frame <- ttkframe(tt, padding = 10)
    tkgrid(input_frame, padx = 20, pady = 10)
    
    # File Selection Section
    file_frame <- ttklabelframe(input_frame, text = "File Selection", padding = 10)
    tkgrid(file_frame, row = 0, column = 0, padx = 10, pady = 5, sticky = "ew")
    
    file_label <- ttklabel(file_frame, text = "Choose CSV File:", font = "Helvetica 12 bold")
    file_entry <- ttkentry(file_frame, textvariable = file_var, width = 50)
    file_button <- ttkbutton(file_frame, text = "Browse", command = function() {
      file_path <- tclvalue(tkgetOpenFile(filetypes = "{{CSV Files} {.csv}} {{All files} *}"))
      if (file_path != "") {
        tclvalue(file_var) <- file_path
      }
    })
    
    tkgrid(file_label, row = 0, column = 0, pady = 5, sticky = "w")
    tkgrid(file_entry, row = 1, column = 0, padx = 10)
    tkgrid(file_button, row = 1, column = 1, padx = 5)
    
    # Protein Sequence Input Section
    protein_frame <- ttklabelframe(input_frame, text = "Normalization", padding = 10)
    tkgrid(protein_frame, row = 1, column = 0, padx = 10, pady = 5, sticky = "ew")
    
    protein_label <- ttklabel(protein_frame, text = "Enter Protein Sequence for Normalization:", font = "Helvetica 12 bold")
    protein_entry <- ttkentry(protein_frame, textvariable = protein_seq_var, width = 50)
    enter_button <- ttkbutton(protein_frame, text = "Enter", command = function() {
      sequence <- tclvalue(protein_seq_var)
      cat("Protein sequence:", sequence, "\n")
    })
    
    tkgrid(protein_label, row = 0, column = 0, pady = 5, sticky = "w")
    tkgrid(protein_entry, row = 1, column = 0, padx = 10, pady = 5)
    tkgrid(enter_button, row = 1, column = 1, padx = 5)
    
    # Frame for plot type, position, and normalization selection (below file and protein input)
    selection_frame <- ttkframe(tt, padding = 10)
    tkgrid(selection_frame, padx = 20, pady = 10)
    
    # Plot Type Selection Section
    plot_type_label <- ttklabel(selection_frame, text = "Select Plot Type:", font = "Helvetica 12 bold")
    plot_type_combo <- ttkcombobox(selection_frame, values = c("Logo Plot", "Bar Plot"), textvariable = plot_type_var, state = "readonly")
    
    tkgrid(plot_type_label, row = 0, column = 0, pady = 5, sticky = "w")
    tkgrid(plot_type_combo, row = 1, column = 0, padx = 10, pady = 5)
    
    # Position Selection Section (disabled initially)
    position_label <- ttklabel(selection_frame, text = "Select Position for Bar Plot:", font = "Helvetica 12 bold")
    position_menu <- ttkcombobox(selection_frame, values = c("P1", "P1'"), textvariable = position_var, state = "disabled")
    
    tkgrid(position_label, row = 2, column = 0, pady = 5, sticky = "w")
    tkgrid(position_menu, row = 3, column = 0, padx = 10, pady = 5)
    
    # Normalization Selection Section
    normalization_label <- ttklabel(selection_frame, text = "Select Normalization Option:", font = "Helvetica 12 bold")
    normalization_combo <- ttkcombobox(selection_frame, values = c("Non-Normalized", "Normalized"), textvariable = normalization_var, state = "readonly")
    
    tkgrid(normalization_label, row = 4, column = 0, pady = 5, sticky = "w")
    tkgrid(normalization_combo, row = 5, column = 0, padx = 10, pady = 5)
    
    # Generate Plot Button
    generate_plot_button <- ttkbutton(tt, text = "Generate Plot", command = function() {
      plot_type <- tclvalue(plot_type_var)
      if (plot_type == "Logo Plot") {
        normalized <- ifelse(tclvalue(normalization_var) == "Normalized", TRUE, FALSE)
        generateLogoPlot(normalized, tclvalue(file_var), tclvalue(protein_seq_var))
      } else {
        generateBarPlot(tclvalue(position_var), tclvalue(file_var))
      }
    })
    tkgrid(generate_plot_button, padx = 20, pady = 10)
    
    # Function to toggle the state of the position drop-down (enabled/disabled based on plot type)
    togglePositionMenu <- function() {
      if (tclvalue(plot_type_var) == "Bar Plot") {
        tkconfigure(position_menu, state = "readonly")  # Enable position menu for bar plot
      } else {
        tkconfigure(position_menu, state = "disabled")  # Disable position menu for logo plot
      }
    }
    
    # Bind the drop-down menu selection to toggle position menu state
    tkbind(plot_type_combo, "<<ComboboxSelected>>", function() togglePositionMenu())
  }
  
  # Function to generate the bar plot for P1 or P1'
  generateBarPlot <- function(position, file_path) {
    # Validate file selection
    if (file_path == "") {
      tkmessageBox(message = "Please choose a CSV file.", icon = "error")
      return()
    }
    
    # Read the CSV file
    data <- tryCatch({
      read.csv(file_path, skip = 3, stringsAsFactors = FALSE)
    }, error = function(e) {
      tkmessageBox(message = "Error reading the CSV file.", icon = "error")
      return(NULL)
    })
    
    if (is.null(data)) return()
    
    # Ensure 'pep_seq' column is present
    if (!"pep_seq" %in% colnames(data)) {
      tkmessageBox(message = "The uploaded file does not contain a 'pep_seq' column.", icon = "error")
      return()
    }
    
    # Remove duplicates from 'pep_seq'
    data <- data[!duplicated(data$pep_seq), ]
    
    # Extract P1 or P1' positions
    sequences <- data$pep_seq
    position_letters <- sapply(sequences, function(seq) {
      n <- nchar(seq)
      if (n >= 4) {
        if (position == "P1") {
          return(substr(seq, n, n))  # Last position (P1)
        } else {
          return(substr(seq, 1, 1))  # First position (P1')
        }
      } else {
        return(NA)  # Handle sequences shorter than 4 AA
      }
    })
    
    position_letters <- na.omit(position_letters)
    position_freq <- table(position_letters) / length(position_letters)
    
    # Generate bar plot
    bar_plot <- ggplot(data = as.data.frame(position_freq), aes(x = position_letters, y = Freq)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      labs(title = paste0("Bar Plot for ", position), x = position, y = "Frequency") +
      theme_minimal()
    
    print(bar_plot)
  }
  
  # Function to generate the sequence logo plot (normalized or non-normalized)
  generateLogoPlot <- function(normalized, file_path, protein_seq) {
    # Validate file selection
    if (file_path == "") {
      tkmessageBox(message = "Please choose a CSV file.", icon = "error")
      return()
    }
    
    # Read the CSV file
    data <- tryCatch({
      read.csv(file_path, skip = 3, stringsAsFactors = FALSE)
    }, error = function(e) {
      tkmessageBox(message = "Error reading the CSV file.", icon = "error")
      return(NULL)
    })
    
    if (is.null(data)) return()
    
    # Ensure 'pep_seq' column is present
    if (!"pep_seq" %in% colnames(data)) {
      tkmessageBox(message = "The uploaded file does not contain a 'pep_seq' column.", icon = "error")
      return()
    }
    
    # Remove duplicates from 'pep_seq'
    data <- data[!duplicated(data$pep_seq), ]
    
    # Extract sequences from the CSV file
    sequences <- data$pep_seq
    
    # Process sequences and generate sequence logo plot
    sequence_matrix <- sapply(sequences, function(seq) {
      n <- nchar(seq)
      if (n >= 4) {
        p4_p1 <- c(substr(seq, n - 3, n - 3), 
                   substr(seq, n - 2, n - 2), 
                   substr(seq, n - 1, n - 1), 
                   substr(seq, n, n))
        p1p_p4p <- c(substr(seq, 1, 1), 
                     substr(seq, 2, 2), 
                     substr(seq, 3, 3), 
                     substr(seq, 4, 4))
        return(c(p4_p1, p1p_p4p))
      } else {
        available_part <- c(strsplit(seq, "")[[1]], rep(NA, 8 - n))
        return(available_part)
      }
    }, USE.NAMES = FALSE)
    
    sequence_matrix <- t(sequence_matrix)
    
    if (normalized) {
      # Calculate frequencies for each letter at each position, ignoring NAs
      position_freq <- apply(sequence_matrix, 2, function(pos) {
        table(factor(pos, levels = LETTERS), useNA = "no") / sum(!is.na(pos))
      })
      
      # Ensure protein sequence covers the required positions
      if (nchar(protein_seq) < 10) {
        tkmessageBox(message = "The protein sequence must be long enough to cover the required positions.", icon = "error")
        return()
      }
      
      # Calculate the overall frequency of each letter in the protein sequence
      protein_freq <- table(factor(unlist(strsplit(protein_seq, "")), levels = LETTERS)) / nchar(protein_seq)
      
      # Normalize the position frequencies by the protein sequence frequencies
      normalized_freq <- sweep(position_freq, 1, protein_freq, "/")
      
      # Replace any non-finite or negative values with 0
      normalized_freq[!is.finite(normalized_freq) | normalized_freq <= 0] <- 0
      
      # Generate the sequence logo plot using ggseqlogo
      logo_plot <- ggseqlogo(normalized_freq, method = "prob")
    } else {
      # For non-normalized sequences
      selected_positions <- sapply(sequences, function(seq) {
        n <- nchar(seq)
        if (n >= 4) {
          paste0(substr(seq, n - 3, n - 3),  # P4
                 substr(seq, n - 2, n - 2),  # P3
                 substr(seq, n - 1, n - 1),  # P2
                 substr(seq, n, n),          # P1
                 substr(seq, 1, 1),          # P1'
                 substr(seq, 2, 2),          # P2'
                 substr(seq, 3, 3),          # P3'
                 substr(seq, 4, 4))          # P4'
        } else {
          paste0(rep("NA", 8), collapse = "")
        }
      }, USE.NAMES = FALSE)
      
      # Remove any positions that contain "NA"
      selected_positions <- selected_positions[!grepl("NA", selected_positions)]
      
      # Create the sequence logo plot using ggseqlogo
      logo_plot <- ggseqlogo(selected_positions, method = "prob")
    }
    
    # Customize x-axis labels to include "cs" and other positions
    logo_plot <- logo_plot + 
      scale_x_continuous(breaks = c(1:8, 9), 
                         labels = c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'", "cs")) +
      geom_vline(xintercept = 4.5, linetype = "dashed", color = "black", size = 1)  # Cleavage site at "cs"
    
    # Plot the sequence logo
    print(logo_plot)
  }
  
  # Run the GUI
  generateGUI()
}

#cut_sites_distribution()


################################################################################
#   display_protease_cut_sites
#   Plot the cut sites of a given protease on a protein sequence.
# 
#   Parameters:
#   - prot (character): Name of the protease.
#   - colour (character): Color for the lines on the plot.
# 
#   Returns:
#   - Invisible: Returns the protease name invisibly.
# 
#   This function reads protease cut site data from a CSV file and visualizes
#   the cut sites of the specified protease on a protein sequence. It retrieves
#   the current gene information and finds the corresponding sequence based on
#   the gene index. Then, it identifies the cut sites of the specified protease
#   within the sequence and plots vertical lines at those positions using the
#   provided color.
# 
#   Note: The function assumes that the necessary data structures like
#   'globalSettings', 'species', and 'protCutSites' are available in the
#   environment.
# 
#   Example usage:
#   plotCutSite("Trypsin", "black")
#   """
# display_protease_cut_sites <- function() {
#   tclCheck()
#   dlg <- myToplevel('cs')
#   if (is.null(dlg))
#     return(invisible())
  
#   # Open a file dialog to choose the protease data CSV file
#   file_path <- tclvalue(tcltk::tkgetOpenFile())
  
#   if (identical(file_path, "")) {
#     tkmessageBox(
#       title = "Error",
#       message = "No file selected. Please choose a CSV file.",
#       icon = "error"
#     )
#     return(invisible())
#   }
  
#   # Read protease data from the chosen CSV file
#   protCutSites <- read.csv(file_path)
  
#   # Create the main window
#   mainWindow <- tktoplevel()
#   tkwm.title(mainWindow, "Plot Cut Site")

#   # Create protease label and dropdown menu
#   proteaseLabel <- tklabel(mainWindow, text = "Select Protease:")
#   tkgrid(proteaseLabel, padx = 10, pady = 10)
  
#   proteaseDropdown <- tklistbox(mainWindow)
#   tkgrid(proteaseDropdown, padx = 10, pady = 10)
  
#   # Populate protease dropdown menu with protease names
#   proteaseList <- unique(protCutSites$protease)
#   for (prot in proteaseList) {
#     tkinsert(proteaseDropdown, "end", prot)
#   }
  
#   # Create color label and entry field
#   colorLabel <- tklabel(mainWindow, text = "Enter Color for Lines:")
#   tkgrid(colorLabel, padx = 10, pady = 10)
  
#   colorEntry <- tkentry(mainWindow)
#   tkgrid(colorEntry, padx = 10, pady = 10)
  
#   # Create plot button
#   plotButton <- tkbutton(mainWindow, text = "Plot", command = function() {
#     # Get the selected protease and color
#     selectedProteaseIndex <- tclvalue(tkcurselection(proteaseDropdown))
#     selectedProtease <- tkget(proteaseDropdown, selectedProteaseIndex)
#     selectedColor <- tkget(colorEntry)
    
#     # Call the plotCutSite function with the selected protease and color
#     plotCutSite(selectedProtease, selectedColor, protCutSites, colorEntry)
    
#     # Close the GUI window
#     tkdestroy(mainWindow)
#   })
#   tkgrid(plotButton, padx = 10, pady = 10)
# }
# Function to plot cut sites
# plotCutSite <- function(prot, colour, protCutSites, colorEntry) {
#   # Retrieve current gene information
#   currGene <- globalSettings$geneDisp
  
#   # Find the index of the current gene in the species data
#   idx <- which(species$genes$name == currGene)
  
#   if (length(idx) > 0) {
#     # Retrieve the corresponding sequence based on the gene index
#     start <- species$genes$seqStartIdx[idx]
#     end <- start + species$genes$seqLength[idx]
#     currSeq <- substr(species$seq, start, end)
#   } else {
#     stop("No current gene information found. Function doesn't work with proteome-wide view. Call a specific protein.")
#   }
  
#   # Set gene sequence and colour variables
#   gene <- globalSettings$geneDisp
#   seq <- currSeq
#   colour <- tclvalue(tkget(colorEntry))
  
#   # Find the index of the protease in the protCutSites data
#   protIndex <- which(as.character(protCutSites$protease) == as.character(prot))
#   col_num <- ncol(protCutSites)
  
#   # Get the cut sites for the protease
#   protase_cutsite <- protCutSites[protIndex, 3:col_num]
#   loc <- 0
  
#   # Iterate through the cut sites
#   for (j in protase_cutsite) {
#     if (j == "") {
#       break
#     }
#     log_message('################################################################')
#     log_message(j)
#     log_message('################################################################')
    
#     prev_index <- 0
#     protcuts <- strsplit(seq, j)
    
#     # Iterate through the sequence segments
#     for (k in protcuts[[1]]) {
#       k <- paste(k, j, sep = '')
#       log_message(k)
      
#       # Find the location of the segment in the sequence
#       loc <- gregexpr(k, seq, fixed = TRUE)
#       start <- loc[[1]][1]
#       len <- nchar(k)
      
#       # Check if the position is greater than 0 and the segment length is greater than 0
#       if (start > 0 && len > 0) {
#         # Calculate the plot site
#         plot_site <- start + len - 1
        
#         if (k == "K") {
#           # Check if it is the last character in the sequence
#           if (plot_site == nchar(seq)) {
#             plot_site <- prev_index
#           } else {
#             plot_site <- prev_index + 1
#           }
#           log_message('IS IT HERE?')
#         }
        
#         log_message(plot_site)
        
#         # Assuming you have initialized a plot or set up the plotting device
#         abline(v = plot_site, col = colour, lty = 2)
        
        
#         prev_index <- plot_site
#       }
#     }
#   }
#   return(invisible(prot))
# }

#display_protease_cut_sites()	
##############################################################################################
#' Protease Cut Site Plotter
#'
#' This function provides a graphical user interface (GUI) for visualizing protease cut sites 
#' within a gene sequence. Users can select a CSV file containing protease cut site data, 
#' choose specific proteases, input a color for plotting, and generate cleavage site plots for 
#' selected proteases.
#'
#' @return A GUI that allows users to upload a CSV file, select proteases, specify a color, and plot 
#' protease cleavage sites within the gene sequence.
#'
#' @details
#' The function provides the following features:
#' \itemize{
#'   \item \strong{File Upload:} Users can select a CSV file that contains protease cut site data.
#'   \item \strong{Protease Selection:} A listbox is populated with proteases from the uploaded CSV file. 
#'   Users can select one or more proteases to plot.
#'   \item \strong{Color Selection:} Users can input a color for the plotted cleavage sites.
#'   \item \strong{Cleavage Site Plotting:} The function generates cleavage site plots for the selected 
#'   proteases, displaying the locations of the cuts on the gene sequence.
#' }
#'
#' @section CSV File Format:
#' The CSV file must contain a column named \code{protease} and columns that represent the cut site positions. 
#' Each row should correspond to a protease, and the columns from the third column onward should represent 
#' cleavage site locations.
#'
#' @section Instructions:
#' \itemize{
#'   \item Upload a CSV file using the "Browse" button.
#'   \item Select one or more proteases from the populated listbox.
#'   \item Input a color for the cleavage site lines in the plot.
#'   \item Click the "Plot cleavage sites" button to generate the cleavage site plots.
#' }
#'
#' @section Dependencies:
#' The function depends on the following R packages:
#' \itemize{
#'   \item \code{tcltk}: For creating the graphical user interface.
#'   \item \code{tkrplot}: For handling interactive plotting.
#'   \item \code{utils}: For reading the CSV file containing protease cut sites.
#' }
#'
#' @examples
#' \dontrun{
#' display_protease_cut_sites()
#' }
#'
#' @export
display_protease_cut_sites <- function() {
  # Load required libraries
  library(tcltk)
  library(tkrplot)
  
  # Create an environment to store global variables like protCutSites
  data_env <- new.env()
  
  # Main Tk window
  mainWindow <- tktoplevel()
  tkwm.title(mainWindow, "Protease Cut Site Plotter")
  
  # Create a labeled frame for file selection
  fileSelectionFrame <- ttklabelframe(mainWindow, text = "File Selection", padding = 10)
  tkgrid(fileSelectionFrame, padx = 20, pady = 10)
  
  # CSV File Selection label, entry, and button
  fileEntry <- ttkentry(fileSelectionFrame, width = 50)
  browseButton <- ttkbutton(fileSelectionFrame, text = "Browse", command = function() {
    file_path <- tclvalue(tkgetOpenFile())
    
    # Check if a file was selected
    if (file_path != "") {
      tkdelete(fileEntry, 0, "end")
      tkinsert(fileEntry, 0, file_path)
      
      # Read the CSV file and store it in the environment
      data_env$protCutSites <- read.csv(file_path)
      
      # Populate the protease listbox
      tkdelete(proteaseDropdown, 0, "end")
      proteaseList <- unique(data_env$protCutSites$protease)
      for (prot in proteaseList) {
        tkinsert(proteaseDropdown, "end", prot)
      }
    }
  })
  
  tkgrid(fileEntry, row = 1, column = 0, padx = 10)
  tkgrid(browseButton, row = 1, column = 1, padx = 5)
  
  # Protease Selection Frame
  proteaseFrame <- ttklabelframe(mainWindow, text = "Protease Selection", padding = 10)
  tkgrid(proteaseFrame, padx = 20, pady = 10)
  
  # Add a listbox with both scrollbars and set the background to white
  yscroll <- tkscrollbar(proteaseFrame, orient = "vertical", command = function(...) tkyview(proteaseDropdown, ...))
  xscroll <- tkscrollbar(proteaseFrame, orient = "horizontal", command = function(...) tkxview(proteaseDropdown, ...))
  
  proteaseDropdown <- tklistbox(proteaseFrame, height = 10, width = 40, selectmode = "multiple", 
                                yscrollcommand = function(...) tkset(yscroll, ...),
                                xscrollcommand = function(...) tkset(xscroll, ...),
                                background = "white", font = "Helvetica 12 bold", fg = "black")
  
  tkgrid(proteaseDropdown, row = 1, column = 0, pady = 10)
  tkgrid(yscroll, row = 1, column = 1, sticky = "ns")
  tkgrid(xscroll, row = 2, column = 0, sticky = "ew")
  
  # Create color label and entry field
  colorFrame <- ttklabelframe(mainWindow, text = "Color Selection", padding = 10)
  tkgrid(colorFrame, padx = 20, pady = 10)
  
  colorEntry <- ttkentry(colorFrame, width = 20, textvariable = tclVar("black"))
  tkgrid(colorEntry, row = 1, column = 0, padx = 10)
  
  # Plot button styled similarly to the other buttons
  plotButton <- ttkbutton(mainWindow, text = "Plot cleavage sites", command = function() {
    file_path <- tclvalue(tkget(fileEntry))
    
    if (identical(file_path, "")) {
      tkmessageBox(
        title = "Error",
        message = "No file selected. Please choose a CSV file.",
        icon = "error"
      )
      return(invisible())
    }
    
    # Get all selected protease indices
    selectedProteaseIndices <- as.integer(tkcurselection(proteaseDropdown))
    
    if (length(selectedProteaseIndices) == 0) {
      tkmessageBox(
        title = "Error",
        message = "No protease selected.",
        icon = "error"
      )
      return(invisible())
    }
    
    # Get all selected proteases
    selectedProteases <- sapply(selectedProteaseIndices, function(i) tclvalue(tkget(proteaseDropdown, i)))
    selectedColor <- tclvalue(tkget(colorEntry))
    
    # Call the plotCutSite function for each selected protease
    for (protease in selectedProteases) {
      plotCutSite(protease, selectedColor, data_env$protCutSites, colorEntry)
    }
  })
  tkgrid(plotButton, pady = 20)
  
  # Start the main window loop without blocking the R console
  tkfocus(mainWindow)
  
  response <- tclvalue(tkmessageBox(
    title = "Instructions", 
    message = paste(
      "1- Upload a CSV file containing proteases cleavage sites", 
      "An example file can be downloaded at:",
      "https://github.com/LewisResearchGroup/digestR",
      "",
      "2- Select one or multiple from the listbox.",
      "",
      "3- Type a color in the entry box",
      sep = "\n"
    )
  ))
  
  # Check the response
  if (response == "ok") {
    print("Roger, Roger")
  } else {
    print("Other action taken.")
  }
}

# Updated plotCutSite function
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
    
    prev_index <- 0
    protcuts <- strsplit(seq, j)
    
    # Iterate through the sequence segments
    for (k in protcuts[[1]]) {
      k <- paste(k, j, sep = '')
      
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
        }
        
        # Assuming you have initialized a plot or set up the plotting device
        abline(v = plot_site, col = colour, lty = 2)
        
        prev_index <- plot_site
      }
    }
  }
  
  # Print message to the console
  #cat("Cleavage sites plotted for", prot, "\n")
  cat(prot,"cleavage sites plotted", "\n")
  
  return(invisible(prot))
}

# Call the display_protease_cut_sites function to display the GUI
#display_protease_cut_sites()
			       
##############################################################################################

#' Generate a graphical user interface for the DigestR BioMart Downloader.
#'
#' This function sets up a Tkinter GUI for searching and downloading protein data using the BioMart API.
generate_proteome <- function() {
  library(tcltk2)
  library(biomaRt)
  
  # Initialize main window
  main_window <- tktoplevel()
  tkwm.title(main_window, "DigestR BioMart Downloader")
  
  # Set initial values for tcl variables
  biomart_var <- tclVar("genes")
  dataset_var <- tclVar("none")
  chromosomes_var <- tclVar("1, 2")
  search_pattern_var <- tclVar("taurus")
  
  #' Update the dataset tcl variable based on user selection in the listbox.
  update_dataset_var <- function(...) {
    selected_index <- tclvalue(tkcurselection(listbox))
    if (length(selected_index) > 0) {
      selected_entry <- tkget(listbox, selected_index)
      selected_dataset <- unlist(strsplit(as.character(selected_entry), "\\|"))[1]
      tclvalue(dataset_var) <- selected_dataset
      log_message('Selected dataset:', selected_dataset)
    }
  }
  
  #' Populate the listbox widget with search results.
  populate_listbox_with_search_results <- function(search_results) {
    tkdelete(listbox, 0, "end")
    for (i in 1:nrow(search_results)) {
      entry <- with(search_results[i, ], {
        sprintf("%-25s | %-30s | %-20s", dataset, description, version)
      })
      tkinsert(listbox, "end", entry)
    }
  }
  
  #' Search datasets using the user-provided search pattern.
  execute_dataset_search <- function() {
    search_pattern <- tclvalue(search_pattern_var)
    if (search_pattern != "") {
      ensembl <- useEnsembl(biomart = "genes")
      search_results <- searchDatasets(mart = ensembl, pattern = search_pattern)
      populate_listbox_with_search_results(search_results)
    }
  }
  
  #' Download protein data based on user settings and selections.
  download_proteome_data <- function() {
    mart <- tclvalue(biomart_var)
    dataset <- tclvalue(dataset_var)
    chromosomes <- tclvalue(chromosomes_var)

    log_message(chromosomes)

    if (chromosomes == "") {
      chromosome_list <- NULL
    } else {
      chromosome_list <- strsplit(chromosomes, ", ?")[[1]]
    }  
    
    biomart_instance <- BioMartData$new(biomart = mart, dataset = dataset)
    biomart_instance$get_data(chromosomes = chromosome_list)
  }
  
  # Create and configure widgets
  tkgrid(tklabel(main_window, text = "Biomart"),
         tk2combobox(main_window, values = c("genes", "ensembl"), textvariable = biomart_var, state = "readonly"),
         padx = 10, pady = 10, sticky = "nsew")

  tkgrid(tklabel(main_window, text = "Search Pattern"),
         tkentry(main_window, textvariable = search_pattern_var, width = 30, bg = 'white'),
         tkbutton(main_window, text = "Search Datasets", command = execute_dataset_search),
         padx = 10, pady = 10, sticky = "nsew")
  
  tkgrid(tklabel(main_window, text = "Dataset Results"),
         listbox <- tklistbox(main_window, height = 20, width = 100, font = 'Consolas 10', bg = 'white'),
         padx = 10, pady = 10, sticky = "nsew")
  
  tkbind(listbox, "<Double-1>", update_dataset_var)
  
  tkgrid(tklabel(main_window, text = "Chromosomes"),
         tkentry(main_window, textvariable = chromosomes_var, width = 30, bg = 'white'),
         tkbutton(main_window, text = "Download proteome", command = download_proteome_data),
         padx = 10, pady = 10, sticky = "nsew")
  
  # Configure grid layout to make it expandable
  tkgrid.columnconfigure(main_window, 0, weight = 1)
  tkgrid.rowconfigure(main_window, 3, weight = 1)

 response <- tclvalue(tkmessageBox(
    title = "Instructions and Warnings", 
    message = paste(
      "Instructions:",
      "1- Enter a search pattern and double-click on a data set to select it.",
      "2- Enter specific chromosomes in the Chromosome box, or leave the field blank to download the full proteome.",
      "",
      "Warnings:",
      "This step is long and may take several minutes to several hours.",
      "Proteomes including H. sapiens, M. musculus, B. taurus, and D. melanogaster, E.coli, P. aeruginosa, P.falciparum can be found on our GitHub page.",
      
      sep = "\n"
    )
  ))
  
  # Check the response
  if (response == "ok") {
    print("Roger, Roger")
  } else {
    print("Other action taken.")
  }
}


