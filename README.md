# #############################################################################
# Run the following lines in the R console
# Author: Dimitri Desmonts de Lamache 
###############################################################################


library(Rcpp)
library(tcltk)
setwd('C:/Users/dimit/OneDrive/Bureau/DigestR')
source('2017_09_16_rNMR_Travis_DD_Edits_Final.R')
source('digestR_Code_Dimitri_final.R')

################################################################

pairedColors <- 	c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
					  "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")

accentColors <- 	c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666")
	



##############################
> newX <- c(myX, rev(myX), myX[1])
> newX
[1]  1  2  3  4  5  6  7  8  9 10 10  9  8  7  6  5  4  3  2  1  1
> newY <- c(myAvg - tmp, rev(myAvg + tmp), (myAvg - tmp)[1])
> newY
[1] 2.308034 4.020336 4.244265 1.824281 5.607191 3.735477 3.844965 4.856742
[9] 5.715681 5.364692 7.709997 8.296436 6.068136 6.535266 5.178285 8.113273
[17] 3.734374 5.818067 5.982386 3.422849 2.308034
> polygon(newX, newY, density = 1, col='blue')
> polygon(newX, newY, density = 1, col='blue')
Error in segments(lx1, ly1, lx2, ly2, ...) : 
		plot.new has not been called yet
> plot(myX, myAvg, type='l')
> polygon(newX, newY, col='blue')
> plot(myX, myAvg, type='l')
> polygon(newX, newY, col='blue')


for(i in 1:length(temp))
{
	temp2[i] <- length(which(dat$compound == temp[i]))
}

##########################
in.folder = fileFolder[[wc()]]
w1Range=in.folder$graphics.par$usr[3:4] 
w2Range=in.folder$graphics.par$usr[1:2] 
col = in.folder$graphics.par$proj.color 
type = in.folder$graphics.par$type	
xlab = NULL
ylab = NULL 
main = in.folder$file.par$user_title 
roiMax = globalSettings$roiMax 
add = FALSE
axes = TRUE
offset = 0
geneNames = globalSettings$geneDisp


########################### Parsing Mascot Daemon Error Log #############
## add headers: Message, Time, Date, Task_ID, Task_Label, File
## Delete final column
## save as 'normal' csv file... by default is saved in unicode.
tmp <- myOpen()
myFile <- read.csv(tmp, stringsAsFactors=FALSE)
ofInterest <- which(myFile$Message == "Error")
myFile[ofInterest, c(5,6)]

###########################
#system.time(res <- readDIANA(loadFileName))
#summaryRprof("d:/Project Files/prof.out")
#loadFileName <- "R:/Travis/Mascot Output/01022013Llinas_Ian_Katie/3B-1_14_450_1800_1+.dcf"
###############################################
lSpecies <- loadSpecies(globalSettings$speciesFiles[globalSettings$processSpeciesID])
dir <- choose.dir()
lFilesCSV <- myOpen(initialdir = dir, multiple = TRUE)
resList <- parseMascot(lFilesCSV)
dfPeptides <- resList$pep
sGenotype <- resList$genotype
fileName <- resList$fileName
pep <- dfPeptides
targetSeq <- lSpecies$seq

totPeps <- length(pep$ID)
indices <- vector(mode = "integer")
lMatch <- list(indices = indices, pepLength = 0, score = 0)

lMatches <- vector(mode = "list", length = totPeps)
nonMatchIndices <- vector(mode = "integer", length=0)

numMatches <- 0
totPeps <- 250




text(in.folder$file.par$gOI$genomeAxis, jitter(in.folder$data[in.folder$file.par$gOI$geneAxis]), 
		pos = 4, col = 'blue', labels = in.folder$file.par$gOI$name, cex = in.folder$graphics.par$cex.axis * 0.8)


###################
#Old code and test code below
###################
tmp <- 'abcdefghijklmnopq'
tmp1 <- unlist(strsplit(tmp, "", fixed = TRUE)) ## convert string to vector
tmp2 <- sample(tmp1) ## randomly scramble vector of characters 
tmp3 <- paste(tmp2, collapse = '') ## concatonate vector of characters back into a string

#library(Rcpp)
#sourceCpp('C:/Users/tsbingem/cpp/test/src/test.cpp')
#> testCpp()
#[1] 3
# 
# 
### combine all currently open files in fileFolder


#####
idx <- which(aaMap > 0)


####

i <- 1
myPath <- extractPathAndName(fileFolder[[1]]$file.par$file.name)$path
newName <- fileFolder[[1]]$file.par$genotype
sFileNames <- fileFolder[[1]]$file.par$file.name
while(i < length(fileFolder))
{
	i <- i + 1
	sFileNames <- c(sFileNames, fileFolder[[i]]$file.par$file.name)
}
coAdd(sFileNames, paste(myPath, newName, '_Summary.dcf', sep='', collapse=''))




### go through fileFolder and combine by pairs
i <- 1
while(i < length(fileFolder))
{
	temp <- c(fileFolder[[i]]$file.par$file.name, fileFolder[[i+1]]$file.par$file.name)
	result <- extractPathAndName(temp[1])
	newName <- paste(result$name, '&', extractPathAndName(temp[2])$name, sep='', collapse='')
	## add back path and extension before submitting name to cAdd
	coAdd(temp, paste(result$path, newName, '.dcf', sep='', collapse=''))
	i <- i + 2
}

for(i in 1:length(fileFolder))
{
	fileFolder[[i]]$graphics.par$plotAA <- TRUE
	myAssign("fileFolder", fileFolder, FALSE)
}

hbbStart <- 5874328
goldberg = c(33, 34, 43, 72, 83, 131)
goldbergLines = goldberg + hbbStart
abline(v=goldbergLines, lty=2, col='grey', lwd=2)

#initDIANA()
#initHuman()
#########################################

tmp <- ucsf2D()
aaMap <- tmp$mapData$aminoAcidMap
t <- table(aaMap)
plot(t, type="l")


x <- 0:length(unique(aaMap))
px <- dpois(x, mean(aaMap))
plot(x, px, type="l", main="Poisson")

t






#########################################
set.seed(1)
x.poi<-rpois(n=numSamp,lambda=2.5) # a vector of random variables from the Poisson distr.

hist(x.poi,main="Poisson distribution")

lambda.est <- mean(x.poi) ## estimate of parameter lambda
#(tab.os<-table(aaMap)) ## table with empirical frequencies


freq.os<-aaMap
#for(i in 1: length(tab.os)) freq.os[i]<-tab.os[[i]]  ## vector of emprical frequencies

freq.ex<-(dpois(0:max(x.poi),lambda=lambda.est)*numSamp) ## vector of fitted (expected) frequencies

acc <- mean(abs(freq.os-trunc(freq.ex))) ## absolute goodness of fit index acc
acc/mean(freq.os)*100 ## relative (percent) goodness of fit index

h <- hist(x.poi ,breaks=length(tab.os))
xhist <- c(min(h$breaks),h$breaks)
yhist <- c(0,h$density,0)
xfit <- min(x.poi):max(x.poi)
yfit <- dpois(xfit,lambda=lambda.est)
plot(xhist,yhist,type="s",ylim=c(0,max(yhist,yfit)), main="Poisson density and histogram")
lines(xfit,yfit, col="red")

#Perform the chi-square goodness of fit test 
#In case of count data we can use goodfit() included in vcd package
library(vcd) ## loading vcd package
gf <- goodfit(x.poi,type= "poisson",method= "MinChisq")
summary(gf)
plot(gf,main="Count data vs Poisson distribution")



loadMenus <- function(top)
{
	tclCheck()
	
	if (is.null(top))
	{
		if (.Platform$OS.type == 'windows')
			top <- myToplevel('menu', width=800, height=600)
		else
			top <- myToplevel('menu', width=800, height=600)
		
		if (is.null(top))
			return(invisible())
		
		tkwm.title(top, 'DIANA')
		tcl('wm', 'attributes', top, topmost=TRUE)
		assign("winMain", top, inherits=FALSE, envir=.GlobalEnv)
	}
	topMenu <- tkmenu(top)
	tkconfigure(top, menu=topMenu)
	
	fileMenu <- tkmenu(topMenu, tearoff=FALSE)
	tkadd(fileMenu, 'command', label='Open files',  
			command=function() file_open())
	
	tkadd(fileMenu, 'command', label='Exit',  
			command=function() file_close())		
	
	tkadd(topMenu, 'cascade', label='File', menu=fileMenu)
	
	editMenu <- tkmenu(topMenu, tearoff=FALSE)
	tkadd(editMenu, 'command', label='Undo', accelerator='ud()', 
			command=function() ud())
	tkadd(editMenu, 'command', label='Redo', accelerator='rd()', 
			command=function() rd())
	
	tkadd(topMenu, 'cascade', label='Edit', menu=editMenu)
	
	graphicsMenu <- tkmenu(topMenu, tearoff=FALSE)
	tkadd(graphicsMenu, 'command', label='Plot colors', 
			accelerator='co()', command=function() co())
	tkadd(graphicsMenu, 'command', label='Plot settings', 
			accelerator='ct()',	command=function() ct())	
	tkadd(graphicsMenu, 'command', label='Perspective', accelerator='per()', 
			command=function() per())
	tkadd(topMenu, 'cascade', label='Graphics', menu=graphicsMenu)
	
	helpMenu <- tkmenu(topMenu, tearoff=FALSE)
	tkadd(helpMenu, 'command', label='Help topics',	
			command=function(...) rNMR:::myHelp('rNMR-package'))
	tkadd(helpMenu, 'command', label='List functions', 
			command=function(...) rNMR:::myHelp('more'))
	tkadd(helpMenu, 'command', label='Update rNMR', 
			command=function() rNMR:::updater())
	tkadd(helpMenu, 'command', label='User manual', 
			command=function(...) rNMR:::myHelp('user_manual', TRUE))
	tkadd(helpMenu, 'command', label='Developer\'s guide', command=function(...) 
				rNMR:::myHelp('developers_guide/developers_guide', TRUE))
	tkadd(helpMenu, 'command', label='Homepage', 
			command=function(...) browseURL('http://rnmr.nmrfam.wisc.edu'))
	tkadd(helpMenu, 'command', label='About rNMR',
			command=function() rNMR:::about())
	tkadd(topMenu, 'cascade', label='Help', menu=helpMenu)
	
	return(top)
}

loadGui2 <- function(top = NULL)
{
	library("tkrplot")
	
	top <- loadMenus(top)
	
	geneEntry <- tclVar("")
	statusText <- tclVar("Status: ")
	assign("status", statusText, inherits=FALSE, envir=.GlobalEnv)
	
	fMain <- tkframe(top)
	fControls <- tkframe(fMain, relief="groove", borderwidth=2)
	fPlot <- tkframe(fMain, borderwidth=0, background = "white")
	fStatus <- tkframe(fMain, relief="groove", borderwidth=2)
	
	geneLabel <- tklabel(fControls, text = "Enter Gene of Interest", justify = "left")
	txtGeneEntry <- tkentry(fControls, width="25", textvariable = geneEntry)
	
	genotypeLabel <- tklabel(fControls, text = "Select Genotype(s)", justify = "left")
	lbGenotype <- tklistbox(fControls, height=3, selectmode="multiple", 
			yscrollcommand=function(...)tkset(scr,...),
			background="white")
	
	butLoadBatchFile <- tkbutton(fControls,text="Load Batch List",command=batchProcess)
	butTestDCF <- tkbutton(fControls,text="Test DCF",command=testDCF)
	butZoom	<- tkbutton(fControls, text="Zoom", command=Zoom)
	
	scr <- tkscrollbar(fControls, repeatinterval=5, command = function(...)tkyview(top$env$gtCombobox, ...))
	img <- tkrplot(fPlot, fun=plotCurrent, hscale=1.5, vscale=1.5)
	
	statusLabel <- tklabel(fStatus, width="80", text=tclvalue(statusText), justify = "left")
	tkconfigure(statusLabel, textvariable=statusText)
	
	tkgrid(geneLabel, txtGeneEntry, genotypeLabel, lbGenotype, scr)
	tkgrid(butLoadBatchFile, butTestDCF, butZoom)
	tkgrid(fControls)
	tkgrid(img)
	tkgrid(fPlot, sticky="nsew")
	tkgrid(statusLabel, sticky="sew")
	tkgrid(fStatus, sticky="sew")
	tkgrid(fMain, sticky="nsew")
	
	tkgrid.columnconfigure(fControls, 0, weight=1)
	tkgrid.rowconfigure(fControls, 0, weight=0)	
	tkgrid.columnconfigure(top, 0, weight=1)
	tkgrid.rowconfigure(top, 0, weight=1)
	tkgrid.columnconfigure(fMain, 0, weight=1)
	tkgrid.rowconfigure(fMain, 0, weight=1)
	tkgrid.columnconfigure(fPlot, 0, weight=1)
	tkgrid.rowconfigure(fPlot, 0, weight=1)
	tkgrid.columnconfigure(img, 0, weight=1)
	tkgrid.rowconfigure(img, 0, weight=1)
	
	tkbind(img, "<Button-1>",buttonDwn)
	tkbind(img, "<ButtonRelease-1>", buttonRls)
	tkbind(top, "<Configure>", resize)
#	tkfocus(top)
#	tkwm.deiconify(top)
}




loadGui <- function(top = NULL)
{
	library("tkrplot")
	
	top <- loadMenus(top)
	
	geneEntry <- tclVar("")
	statusText <- tclVar("Status: ")
	assign("status", statusText, inherits=FALSE, envir=.GlobalEnv)
	
	fControls <- tkframe(top, relief="groove", borderwidth=2)
	fPlot <- tkframe(top, relief="groove", borderwidth=2, background = "white")
	fStatus <- tkframe(top, relief="groove", borderwidth=2)
	
	geneLabel <- tklabel(fControls, text = "Enter Gene of Interest", justify = "left")
	txtGeneEntry <- tkentry(fControls, width="25", textvariable = geneEntry)
	
	genotypeLabel <- tklabel(fControls, text = "Select Genotype(s)", justify = "left")
	lbGenotype <- tklistbox(fControls, width="90", height=5, selectmode="multiple", 
			yscrollcommand=function(...)tkset(scr,...),
			background="white")
	
	assign("selectedFiles", lbGenotype, envir = .GlobalEnv)
	
	butLoadBatchFile <- tkbutton(fControls,text="Load Batch List",command=batchProcess)
#	butTestDCF <- tkbutton(fControls,text="Test DCF",command=testDCF)
	butPlotSelected	<- tkbutton(fControls, text="Plot Selected", command=plotSelected)
	
	scr <- tkscrollbar(fControls, repeatinterval=5, command = function(...)tkyview(lbGenotype, ...))
	img <- tkrplot(fPlot, fun=plotCurrent, hscale=2.6, vscale=1.5)
	
	statusLabel <- tklabel(fStatus, width="140", text=tclvalue(statusText), justify = "left")
	tkconfigure(statusLabel, textvariable=statusText)
	
	tkgrid(geneLabel, txtGeneEntry, genotypeLabel, lbGenotype, scr, sticky="NS")
	tkgrid(butLoadBatchFile, butPlotSelected)
	tkgrid(fControls, sticky="nsew")
	tkgrid(img)
	tkgrid(fPlot, sticky="nsew")
	tkgrid(statusLabel, sticky="nsew")
	tkgrid(fStatus, sticky="nsew")
	
	tkgrid.configure(fControls, column=0, row=0)
	tkgrid.configure(fPlot, column=0, row=1)
	tkgrid.configure(fStatus, column=0, row=2)
	
	tkgrid.columnconfigure(top, 0, weight=1)
	tkgrid.rowconfigure(top, 1, weight=1)	
	
	tkgrid.rowconfigure(fPlot, 1, weight=1)
	
	tkbind(img, "<Button-1>",buttonDwn)
	tkbind(img, "<ButtonRelease-1>", buttonRls)
	tkbind(top, "<Configure>", resize)
	
	assign("currImg", img, env = .GlobalEnv)
	
	populateGenotypeSelector()
	
}

fileSelector <- function()
{
	
	if (.Platform$OS.type == 'windows')
		selector <- myToplevel('menu', width=800, height=600)
	else
		selector <- myToplevel('menu', width=800, height=600)
	
	if (is.null(selector))
		return(invisible())
	
	tkwm.title(selector, 'File Selector')
	tcl('wm', 'attributes', selector, topmost=TRUE)
	
	fFileSelector <- tkframe(selector, relief="groove", borderwidth=2)
	
	lbGenotype <- tklistbox(fFileSelector, width="90", height=5, selectmode="multiple", 
			yscrollcommand=function(...)tkset(scr,...),
			background="white")
	
	scr <- tkscrollbar(fFileSelector, repeatinterval=5, command = function(...)tkyview(lbGenotype, ...))
	
	butOk	<- tkbutton(fFileSelector, text="Ok", command=filesSelected)
	
	dir <- tkchooseDirectory()
	lFiles <- list.files(path = file.path(dir), pattern = '*.DCF', all.files = FALSE, full.names = FALSE, recursive = TRUE, ignore.case = TRUE)
	
	for (i in 1:length(lFiles))
	{
		tkinsert(lbGenotype, "end", lFiles[i])
	}	
	
	tkgrid(lbGenotype, scr)
	tkgrid(butOk)
	tkgrid(fFileSelector)
	
	assign("selectedFiles", lbGenotype, envir = .GlobalEnv)
}

filesSelected <- function()
{
	selectedFiles <- get("selectedFiles", envir = .GlobalEnv)
	tmp <- tkcurselection(selectedFiles)
}

buttonDwn <- function(x,y)
{
	tmp <- get("status", envir=.GlobalEnv)
	assign("xMin", x, envir = .GlobalEnv)
	assign("yMin", y, envir = .GlobalEnv)
	tclvalue(tmp) <- paste("Coords", x, ",", y, sep="", collapse="")	
}

buttonRls <- function(x,y)
{
	tmp <- get("status", envir=.GlobalEnv)
	gX <- get("xMin", envir = .GlobalEnv)
	gY <- get("yMin", envir = .GlobalEnv)
	
	assign("xMin", min(x, gX), envir = .GlobalEnv)
	assign("xMax", max(x, gX), envir = .GlobalEnv)
	
	assign("yMin", min(y, gY), envir = .GlobalEnv)
	assign("yMax", max(y, gY), envir = .GlobalEnv)
	
	tclvalue(tmp) <- paste( " Min x:", get("xMin", envir = .GlobalEnv), 
			" Min y:", get("yMin", envir = .GlobalEnv), 
			" Max x:", get("xMax", envir = .GlobalEnv), 
			" Max y:", get("yMax", envir = .GlobalEnv), sep="", collapse="")
	
	zoom()
}

plotSelected <- function()
{
	xMin <- get("xMin", envir = .GlobalEnv)
	yMin <- get("yMin", envir = .GlobalEnv)
	xMax <- get("xMax", envir = .GlobalEnv)
	yMax <- get("yMax", envir = .GlobalEnv)
	bZoom <- get("bZoom", envir = .GlobalEnv)
	
	xMin <- NULL
	yMin <- NULL
	xMax <- NULL
	yMax <- NULL
	bZoom <- FALSE
	
	stat <- get("status", envir=.GlobalEnv)
	tclvalue(stat) <- "Loading selected files..."
	
	resultList <- readCurrent()
	assign("resultList", resultList, envir = .GlobalEnv)
	
	tclvalue(stat) <- "Files loaded."
	
	img <- get("currImg", envir = .GlobalEnv)
	tkrreplot(img)	
}

resize <- function(x, y)
{
	tmp <- get("status", envir=.GlobalEnv)
	tclvalue(tmp) <- paste("Coords", x, ",", y, sep="", collapse="")		
}


file_open <- function()
{
	fileName <- myOpen("csv", list(csv = "Comma Separated Values Files, xls = Excel Files"), multiple = TRUE)
#	resultList <- vector(mode='list', length=length(fileNames))
#	winMain$env$gtCombobox$listvariable <- fileNames
#	tkconfigure(winMain$env$gtCombobox, values = unlist(resultList))
	#tkconfigure(winMain$env$gtCombobox, values = c('a','b','c')) <- this works
	for (i in 1:length(fileName))
	{
#		tkinsert(winMain$env$gtCombobox, "end", as.character(i))
	}
}

plotCurrent <- function()
{
	params <- par(bg="white")
	resultList <- get("resultList", envir = .GlobalEnv)
	lSpecies <- get("lSpecies", envir = .GlobalEnv)
	if(is.null(resultList))
	{
		x <- -10:10
		y <- 1:21
		plot(x, y, main="test")
	}else
	{
		plotAllByChrom(resultList, lSpecies)
	}
	
	par(params)
}

populateGenotypeSelector <- function()
{
	lb <- get("selectedFiles", envir = .GlobalEnv)
	
	dir <- tkchooseDirectory()
	lFiles <- list.files(path = file.path(dir), pattern = '*.DCF', all.files = FALSE, full.names = FALSE, recursive = TRUE, ignore.case = TRUE)
	assign("fileList", lFiles, envir = .GlobalEnv)
	assign("fileDirectory", dir, envir = .GlobalEnv)
	
	for (i in 1:length(lFiles))
	{
		tkinsert(lb, "end", lFiles[i])
	}
}


plotAllByGene <- function(resultList, lSpecies, genes, yLabel)
{
	maxY <- 0
	seqLen <- nchar(lSpecies$genes$seq)
	cumLen <- cumsum(seqLen)
	vSeq <- NULL
	sLegendLabels <- vector(mode = "character", length(resultList))
	vGeneIdx <- vector(mode = "integer", length = length(genes))
	vLen	<- vector(mode = "integer", length = length(genes))	
	vCum	<- vector(mode = "integer", length = length(genes))
	
	for (i in 1:length(genes))
	{
		vGeneIdx[i] <- which(lSpecies$genes$name == genes[i])
		vLen[i] <- seqLen[vGeneIdx[i]]
		
		# cumLen holds the cumulative length of sequences including the current ID - 
		# therefore to find the starting index for the current ID subtract the 
		# current gene sequence length minus 1
		# this looks after the case of it being the first gene which is not taken care of by
		# using the last index of the previous gene +1
		vCum[i] <- cumLen[vGeneIdx[i]] - (seqLen[vGeneIdx[i]] - 1)	
		
		vSeq <- paste(vSeq, lSpecies$genes$seq[vGeneIdx[i]], sep="")
	}
	vSeq <- unlist(strsplit(vSeq, ""))
	
	# vIndices hold the indices of genes of interest
	# allocate vector to hold indices for genes of interest
	seqIndices <- vector(mode = "integer", length = sum(seqLen[vGeneIdx]))
	
	lIndices <- list()
	for (i in 1:length(genes))
	{
		lIndices[[i]] <- seq(from = vCum[i], to = (vCum[i] + vLen[i] - 1))
	}
	# seqIndices are in amino acid space
	seqIndices <- unlist(lIndices)
	
	idxGeneLabelPos <- vLen / 2 # horizontal position for chromosome labels
	
	for(i in 1:length(resultList))
		maxY <- max(maxY, resultList[[indices[i]]]$aaMap[seqIndices])
	
	maxY <- maxY * 1.05 # maximum y value + 5% to provide padding in plot
	yLim <- c(0, maxY)
	
	plot((1:length(seqIndices)), resultList[[i]]$aaMap[seqIndices], xlab = genes, xaxt = "n", ylab = yLabel, type="l", ylim = yLim)
	sLegendLabels[1] <- resultList[[1]]$Genotype
	drawChromeBounds(cumsum(vLen))
	axis(side = 1, at = (1:length(seqIndices)), labels = vSeq, tick = FALSE, cex.axis=0.5)
	
	for(i in 2:length(resultList))
	{
		lines(1:length(seqIndices), resultList[[i]]$aaMap[seqIndices], type="l", col= i)
		sLegendLabels[i] <- resultList[[i]]$Genotype
	}
	legend("topright", sLegendLabels, lty=c(1,1), col=1:length(seqIndices))
}

plotAllByChrom <- function(resultList, lSpecies, yLabel = "")
{
	maxY <- 0
	
	if (lSpecies$name == "human")
	{
		tar <- c(as.character(1:22), "x", "y")
	}else{
		chrom <- lSpecies$genes$chrom
		tar <- unique(lSpecies$genes$chrom)
	}
	
	idxL <- vector(mode = "integer", length = length(tar))
	idxU <- vector(mode = "integer", length = length(tar))
	
	# set up vectors for indices of upper and lower chromosome boundaries along x-axis
	for (i in 1:length(tar))
	{
		idxL[i] <- min(which(lSpecies$genes$chrom == tar[i]))
		idxU[i] <- max(which(lSpecies$genes$chrom == tar[i]))
	}	
	idxChromLabelPos <- (idxL + idxU) / 2 # horizontal position for chromosome labels
	
	for(i in 1:length(resultList))
		maxY <- max(maxY, resultList[[i]]$geneMap)
	
	maxY <- maxY * 1.05 # maximum y value + 5% to provide padding in plot
	yLim <- c(0, maxY)
	
	plot(resultList[[1]]$geneMap, xlab = "Chromosome", xaxt = "n", ylab = yLabel, type="l", ylim = yLim)
	drawChromeBounds(idxU)
	axis(side = 1, at = idxChromLabelPos, labels = tar, tick = FALSE)
	
	if (length(resultList) > 1)
	{
		for(i in 2:length(resultList))
		{
			lines(resultList[[i]]$geneMap, type="l", col= i)
		}
	}
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
#plot1D <- function(	in.folder = fileFolder[[wc()]],
#		w1Range=in.folder$graphics.par$usr[3:4], 
#		w2Range=in.folder$graphics.par$usr[1:2], 
#		col = in.folder$graphics.par$proj.color, 
#		type = in.folder$graphics.par$type,	
#		xlab = NULL, ylab = NULL, 
#		main = in.folder$file.par$user_title, 
#		roiMax = globalSettings$roiMax, 
#		add = FALSE, axes = TRUE, offset = 0, ... )
#{
#	## Remove any non line types
#	if(!any(type == c('l', 'b', 'p')))
#		type <- 'l'
#
#	## Redirect 2D NMR data to slice/projection function
#	if(in.folder$file.par$number_dimensions > 1 )
#		plot2D(in.folder = in.folder, w1Range = w1Range, w2Range = w2Range,
#				col = col, main = main, add = add, axes = axes, type = type, ...)
#	else
#	{	
#		## Read in dataset if we don't have data
#		if( is.null(c(in.folder$data, in.folder$w2)) )
#		{
#			new.folder <- ucsf2D(file.name = in.folder$file.par$file.name,	w2Range = w2Range, file.par = in.folder$file.par, getNames = TRUE, ...)
#			in.folder$file.par <- new.folder$file.par
#			in.folder$w1 <- new.folder$w1
#			in.folder$w2 <- new.folder$w2
#			in.folder$data <- new.folder$data
#			in.folder$chromDetails <- new.folder$chromDetails
#			in.folder$geneDetails <- new.folder$geneDetails
#			genome <- new.folder$genome
#		}
#		## Erase old plot if add is false ## re-order this code so we only check for !add once ... do after loading new data if necessary
#		if(!add)
#		{
#			w2Range <- sort(w2Range) 
#			
#			if(is.null(xlab))
#				xlab <- "Chromosome"
#			
#			if(is.null(ylab))
#				ylab <- 'Intensity'
#			
#			plot(0,0, axes=FALSE, type= 'n', xlab = "", ylab="", main=main, xaxs="r", yaxs="r")
#			
#			title(main = main, sub=NULL, xlab = xlab, ylab=ylab, line=2)
#			box(col=in.folder$graphics.par$fg, bty='o')
#
#			par(usr = c(w2Range[1], w2Range[2], in.folder$file.par$min_intensity, (in.folder$file.par$max_intensity * 1.05) + offset))			
#	
#			if(axes)
#			{
#				if (!is.na(in.folder$graphics.par$xtck) &&	in.folder$graphics.par$xtck == 1)
#					xlty <- 2
#				else
#					xlty <- 1
#				
#				if (!is.na(in.folder$graphics.par$ytck) && 	in.folder$graphics.par$ytck == 1)
#					ylty <- 2
#				else
#					ylty <- 1
#				
#				drawAndLabelAxis(in.folder, in.folder$chromDetails$chromEnd, in.folder$chromDetails$labelPos, in.folder$chromDetails$chromName, xlty)		
#				
#				axis(side=2, lty=ylty, tck=in.folder$graphics.par$ytck, cex.axis=in.folder$graphics.par$cex.axis)				
#			}			
#		}
#		
#		lines(x=(in.folder$w2 + w2Range[1] - 1), y=in.folder$data + offset, col=col, type=type)
#
#		text((in.folder$file.par$gOI$idx + w2Range[1] - 1), jitter(in.folder$data[(in.folder$file.par$gOI$idx)] + offset), 
#				pos = 4, col = col, labels = in.folder$file.par$gOI$name, cex = in.folder$graphics.par$cex.axis * 0.8)
#		
#		points((in.folder$file.par$gOI$idx + w2Range[1] - 1), in.folder$data[in.folder$file.par$gOI$idx] + offset, col = col, pch = 16, type="p")
#		
#		if( roiMax &&  dev.cur() != 2 )
#		{
#			mShift <- maxShift( in.folder )
#			points( mShift$w2, mShift$Height)
#		}
#	}
#}

#
#writeDIANA_Test <- function(saveFileName, lSpecies, result)
#{
#	endian <- "little"
#	if(saveFileName == "")
#	{
#		saveFileName <- mySave(defaultextension = '.dcf', filetypes = list('dcf' = 'DIANA Compressed File'),	initialfile = '', title= 'Save')
#		result[["File"]] <- saveFileName
#	}
#	
#	writeCon <- file(saveFileName, "w+b")
#	
#	writeHeader(writeCon, lSpecies, result, endian)				
#	writeMatches_Test(writeCon, result, endian)
#	writeSpecies(writeCon, lSpecies, endian)
#	
#	close(writeCon)
#	
#	return(saveFileName)
#} 
#
#readDIANA_Test <- function(loadFileName, lSpecies)
#{
#	endian <- "little"
#	if(loadFileName == "")
#		loadFileName <- myOpen(defaultextension = 'dcf', filetypes = list('dcf' = 'DIANA Compressed File'),	initialfile = "", title= 'Load')
#	
#	readCon <- file(loadFileName, "rb")
#	
#	header <- readHeader(readCon, endian)
#	
#	mat <- readMatches_Test(readCon, endian)
#	
#	## only read species if species defined in header hasn't been loaded already
#	if(is.null(lSpecies) || lSpecies$name != header$species)
#	{
#		### this is a problem.. that we need to handle..  ###
#		### for now, overwrite current species with species stored in current file
#		nSpecies <- readSpecies(readCon, header$species, endian)
#		lSpecies <- get("lSpecies", envir=.GlobalEnv)
#		lSpecies <- nSpecies
#	}
#	
#	close(readCon)	
#	
#	result <- vector(mode="list", length=5)
#	names(result) <- c("Genotype", "File", "Matches", "Genome", "Gene")
#	
#	result[["Genotype"]] <- getGenotype(loadFileName)	
#	
#	result[["File"]] <- header$loadFileName	
#	#result[["Matches"]] <- lMatches
#	
#	mapData <- mapToGene(lSpecies, mat)
#	
#	result[["Genome"]] <- mapData$aminoAcidMap
#	result[["Gene"]] <- mapData$geneMap
#	
#	return(result)
#}



#readHeader <- function(readCon, endi)
#{
#	nameSpace <- 128L
#	startBinary <- readBin(readCon, what = integer(), endian = endi)	
#	versionDIANA <- readBin(readCon, what = numeric(), endian = endi)
#	
#	g_noiseEst <- readBin(readCon, what = numeric(), endian = endi)
#	g_min_intensity <- readBin(readCon, what = numeric(), endian = endi)
#	g_max_intensity <- readBin(readCon, what = numeric(), endian = endi)
#	g_zero_offset <- readBin(readCon, what = numeric(), endian = endi)
#	g_downfield_ppm <- readBin(readCon, what = integer(), endian = endi)
#	g_upfield_ppm <- readBin(readCon, what = integer(), endian = endi)
#	g_center_ppm <- readBin(readCon, what = integer(), endian = endi)
#
#	aa_noiseEst <- readBin(readCon, what = numeric(), endian = endi)
#	aa_min_intensity <- readBin(readCon, what = numeric(), endian = endi)
#	aa_max_intensity <- readBin(readCon, what = numeric(), endian = endi)
#	aa_zero_offset <- readBin(readCon, what = numeric(), endian = endi)
#	aa_downfield_ppm <- readBin(readCon, what = integer(), endian = endi)
#	aa_upfield_ppm <- readBin(readCon, what = integer(), endian = endi)
#	aa_center_ppm <- readBin(readCon, what = integer(), endian = endi)	
#	
#	speciesID <- readBin(readCon, what = integer(), endian = endi)
#	speciesNameSize <- readBin(readCon, what = integer(), endian = endi)
#	speciesName <- readChar(readCon, speciesNameSize + 1) ## not sure if 1 is required, needs to be tested
#	fileNameSize <- readBin(readCon, what = integer(), endian = endi)
#	fileName <- readChar(readCon, fileNameSize + 1)## not sure if 1 is required, needs to be tested
#	
#	padding <- nameSpace - fileNameSize - speciesNameSize + 1
#	
#	pad <- readChar(readCon, padding)# should take to end of fileName padding
#	
#	header <- data.frame(	
#			startBinary = startBinary,
#			version = versionDIANA,
#			speciesID = speciesID,
#			species = speciesName,
#			fileName = fileName,
#			g_noise = g_noiseEst,
#			g_min_int = g_min_intensity,
#			g_max_int = g_max_intensity,
#			g_zero = g_zero_offset,
#			g_downfield = g_downfield_ppm,
#			g_upfield = g_upfield_ppm,
#			g_center = g_center_ppm,
#			aa_noise = aa_noiseEst,
#			aa_min_int = aa_min_intensity,
#			aa_max_int = aa_max_intensity,
#			aa_zero = aa_zero_offset,
#			aa_downfield = aa_downfield_ppm,
#			aa_upfield = aa_upfield_ppm,
#			aa_center = aa_center_ppm,			
#			stringsAsFactors = FALSE)
#	return(header)
#}



#readHeader <- function(readCon, endi)
#{
#	startBinary <- readBin(readCon, what = integer(), endian = endi)	
#	versionDIANA <- readBin(readCon, what = "numeric", endian = endi)
#	speciesID <- readBin(readCon, what = integer(), endian = endi)
#	species <- readBin(readCon, what = "character", endian = endi)
#	fileNameSpace <- readBin(readCon, what = integer(), endian = endi)
#	fileNameSize <- readBin(readCon, what = integer(), endian = endi)
#	fileName <- readChar(readCon, fileNameSize + 1)
#	
#	padding <- readChar(readCon, fileNameSpace - fileNameSize + 1)# should take to end of fileName padding
# 
#	pad <- readChar(readCon, 256 - nchar(species) + 1) # should take to end of header padding
#	
#	header <- data.frame(	
#							startBinary = startBinary,
#							version = versionDIANA,
#							speciesID = speciesID,
#							species = species,
#							fileName = fileName,
#							stringsAsFactors = FALSE)
#	return(header)
#}

##  this version 3.78 seconds
#readSpecies <- function(readCon, endi)
#{
#	lSpecies <- vector(mode="list", length=4)
#	names(lSpecies) <- c("name", "ID", "genes", "seq")	
#	
#	len <-	readBin(readCon, what= integer(), endian=endi)
#	lSpecies$name <- readChar(readCon, len+1)
#	lSpecies$ID	<- readBin(readCon, what= integer(), endian=endi)
#	
#	totGenes <- readBin(readCon, what= integer(), endian=endi)
#	
#	lSpecies$genes <- data.frame("name" = character(totGenes), "chrom" = character(totGenes), 
#							"start" = integer(totGenes), "seqStartIdx" = integer(totGenes), 
#							"seqLength" = integer(totGenes), stringsAsFactors = FALSE)
#					
#	## this loop takes 4 seconds... 7 reads per gene
#	for (i in 1: totGenes)
#	{
#		len <-	readBin(readCon, what = integer(), endian=endi)
#		lSpecies$genes$name[i] <- readChar(readCon, len+1)
#		len <-	readBin(readCon, what = integer(), endian=endi)
#		lSpecies$genes$chrom[i] <- readChar(readCon, len+1)
#		lSpecies$genes$start[i] <- readBin(readCon, what = integer(), endian=endi)
#		lSpecies$genes$seqStartIdx[i] <- readBin(readCon, what = integer(), endian=endi)
#		lSpecies$genes$seqLength[i] <- readBin(readCon, what = integer(), endian=endi)
#	}
#	len <-	readBin(readCon, what = integer(), endian=endi)
#	lSpecies$seq <- readChar(readCon, len+1)
#		
#	return(lSpecies)
#}


#writeSpecies <- function(writeCon, lSpecies, endi)
#{
#	# write the number of genes
#	writeBin(length(lSpecies$genes$name), writeCon, endian=endi)
#	# write each of the gene names
#	for(i in 1:length(lSpecies$genes$name))
#	{
#		writeBin(lSpecies$genes$name[i], writeCon, endian=endi)		
#	}
##	some sequences are too long for max length of 10000 for readBin
##	before we write the sequences out, write the length of each sequence
#	writeBin(nchar(lSpecies$genes$seq), writeCon, endian=endi)
#
#	# write each sequence
#	for(i in 1:length(lSpecies$genes$name))
#	{
#		writeBin(lSpecies$genes$seq[i], writeCon, endian=endi)
#	}
#	# write which chromosome each gene is in
#	for(i in 1:length(lSpecies$genes$name))
#	{
#		writeBin(lSpecies$genes$chrom[i], writeCon, endian=endi)		
#	}
#	# write start location of each gene
#	for(i in 1:length(lSpecies$genes$name))
#	{
#		writeBin(lSpecies$genes$start[i], writeCon, endian=endi)		
#	}	
#	# write end location of each gene
#	for(i in 1:length(lSpecies$genes$name))
#	{
#		writeBin(lSpecies$genes$end[i], writeCon, endian=endi)		
#	}
#}

#readSpecies2 <- function(readCon, name, endi)
#{
#	# read number of genes
#	numGenes <-	readBin(readCon, what= integer(), endian=endi)
#		
#	# allocate vectors
#	names <- vector(mode = "character", length = numGenes)
#	seq <- vector(mode = "character", length = numGenes)	
#	chrom <- vector(mode = "character", length = numGenes)	
#	start <- vector(mode = "integer", length = numGenes)	
#	end <- vector(mode = "integer", length = numGenes)		
#	seqLen <- vector(mode = "integer", length = numGenes)	
#	
#	for(i in 1:numGenes)
#	{
#		names[i] <- readBin(readCon, what = "character", endian=endi)
#	}
#
#	for(i in 1:numGenes)
#	{	
#		seqLen[i] <- as.integer(readBin(readCon, what = "integer", endian=endi))	
#	}
#	
#	for(i in 1:numGenes)
#	{
#		seq[i] <- readChar(readCon, seqLen[i]+1) # +1 to get the string terminator
#	}
#
#	for(i in 1:numGenes)
#	{
#		chrom[i] <- readBin(readCon, what = "character", endian=endi)		
#	}
#	
#	for(i in 1:numGenes)
#	{
#		start[i] <- readBin(readCon, what = integer(), endian=endi)		
#	}
#	
#	for(i in 1:numGenes)
#	{
#		end[i] <- readBin(readCon, what = integer(), endian=endi)		
#	}
#	
#	genes <- data.frame(name = names,
#						seq = seq,
#						chrom = chrom,
#						start = start,
#						end = end,
#						stringsAsFactors = FALSE)
#				
#	lSpecies <- list(name = name,
#					genes = genes)
#			
#	lSpecies$seq <- paste(genes$seq, collapse = '', sep = '')
#	
#	return(lSpecies)
#}

#writeHeader <- function(writeCon, lSpecies, result, endi)
#{
#	versionDIANA <- 0.1 #get("DIANA_VERSION", envir = .GlobalEnv)
#	startBinary <- 243L
#	nameSpace <- 128L	
#	fileName <- result$File
#	fileNameSize <- as.integer(nchar(fileName))
#	speciesNameSize <- as.integer(nchar(lSpecies$name))
#	
#	file.range <- fivenum(result$Gene)
#	noiseEst <- as.numeric((diff(file.range[c(2, 4)]) / 2) / .674)
#	min_intensity = as.numeric(file.range[1])
#	max_intensity = as.numeric(file.range[5])
#	zero_offset = as.numeric(file.range[3])
#	downfield_ppm <- as.integer(1)
#	upfield_ppm <- as.integer(length(lSpecies$genes$name))
#	center_ppm <- as.integer(round(((downfield_ppm + upfield_ppm) / 2)))
#	
#	writeBin(startBinary, writeCon, endian = endi) # 4 bytes
#	writeBin(versionDIANA, writeCon, endian = endi)# 8 bytes = 12
#	
#	writeBin(noiseEst, writeCon, endian = endi) # 8
#	writeBin(min_intensity, writeCon, endian = endi) # 8
#	writeBin(max_intensity, writeCon, endian = endi) # 8
#	writeBin(zero_offset, writeCon, endian = endi)# 8
#	writeBin(downfield_ppm, writeCon, endian = endi) # 4
#	writeBin(upfield_ppm, writeCon, endian = endi) # 4
#	writeBin(center_ppm, writeCon, endian = endi) # 4 = 56
#	### new variable name to hold amino acid level data
#	file.range <- fivenum(result$Genome)
#	noiseEst <- as.numeric((diff(file.range[c(2, 4)]) / 2) / .674)
#	min_intensity = as.numeric(file.range[1])
#	max_intensity = as.numeric(file.range[5])
#	zero_offset = as.numeric(file.range[3])
#	downfield_ppm <- as.integer(1)
#	upfield_ppm <- as.integer(nchar(lSpecies$seq))
#	center_ppm <- as.integer( round( (downfield_ppm + upfield_ppm) / 2) )	
#	
#	writeBin(noiseEst, writeCon, endian = endi)
#	writeBin(min_intensity, writeCon, endian = endi)
#	writeBin(max_intensity, writeCon, endian = endi) # = 80
#	writeBin(zero_offset, writeCon, endian = endi)
#	writeBin(downfield_ppm, writeCon, endian = endi)
#	writeBin(upfield_ppm, writeCon, endian = endi)
#	writeBin(center_ppm, writeCon, endian = endi)	 # = 100
#	
#	writeBin(lSpecies$ID, writeCon, endian = endi) # = 104
#	writeBin(speciesNameSize, writeCon, endian = endi) # 4
#	writeBin(lSpecies$name, writeCon, endian = endi) # 5 = 114 ?
#	writeBin(fileNameSize, writeCon, endian = endi) # 118
#	writeBin(fileName, writeCon, endian = endi)
#	
#	padding <- nameSpace - fileNameSize - speciesNameSize 
#	writeBin(paste(as.character(rep(0, padding)), sep="", collapse=""), writeCon, endian = endi)	
#}

#plotByGene <- function(lSpecies, genes, map, yLabel = "")
#{
#	seqLen <- lSpecies$genes$seqLength
#	cumLen <- cumsum(seqLen)
#	vSeq <- NULL
#	
#	vGeneIdx <- vector(mode = "integer", length = length(genes))
#	vLen	<- vector(mode = "integer", length = length(genes))	
#	vCum	<- vector(mode = "integer", length = length(genes))
#	
#	for (i in 1:length(genes))
#	{
#		tmp <- which(lSpecies$genes$name == genes[i])
#		if (length(tmp > 1))# plot first matching gene if multiple matches
#			vGeneIdx[i] <- tmp[1]
#		else
#			vGeneIdx[i] <- tmp
#		
#		vLen[i] <- seqLen[vGeneIdx[i]]
#		
#		# cumLen holds the cumulative length of sequences including the current ID - 
#		# therefore to find the starting index for the current ID subtract the 
#		# current gene sequence length minus 1
#		# this looks after the case of it being the first gene which is not taken care of by
#		# using the last index of the previous gene +1
#		vCum[i] <- cumLen[vGeneIdx[i]] - (seqLen[vGeneIdx[i]] - 1)	
#
#		vSeq <- paste(vSeq, lSpecies$genes$seq[vGeneIdx[i]], sep="")
#	}
#	vSeq <- unlist(strsplit(vSeq, ""))
#	
#	# vIndices hold the indices of genes of interest
#	# allocate vector to hold indices for genes of interest
#	seqIndices <- vector(mode = "integer", length = sum(seqLen[vGeneIdx]))
#	
#	lIndices <- list()
#	for (i in 1:length(genes))
#	{
#		lIndices[[i]] <- seq(from = vCum[i], to = (vCum[i] + vLen[i] - 1))
#	}
#	# seqIndices are in amino acid space
#	seqIndices <- unlist(lIndices)
#
#	idxGeneLabelPos <- vLen / 2 # horizontal position for chromosome labels
#	
#	maxY <- max(map[seqIndices]) * 1.05 # maximum y value + 5% to provide padding in plot
#	yLim <- c(0, maxY)
#	
#	plot((1:length(seqIndices)), map[seqIndices], xlab = genes, xaxt = "n", ylab = yLabel, type="l", ylim = yLim)
#	drawChromeBounds(cumsum(vLen))
#	axis(side = 1, at = (1:length(seqIndices)), labels = vSeq, tick = FALSE, cex.axis=0.4)
#}

### to do:
### compare currently loaded files to currently selected files
### only read in files not currently loaded.
### unload files not currently needed.
readCurrent <- function()
{
	lFiles <- get("fileList", envir = .GlobalEnv)
	root <- get("fileDirectory", envir = .GlobalEnv)
	lb <- get("selectedFiles", envir = .GlobalEnv)
	stat <- get("status", envir = .GlobalEnv)
	lSpecies <- get("lSpecies", envir = .GlobalEnv)
	
	fileIndices <- as.numeric(tkcurselection(lb)) + 1	# tkcurselection returns a 0-based index, must convert to 1-based
	
	curFiles <- as.character(lFiles[fileIndices])
	numFiles <- length(curFiles) 
	
	resultList <- vector(mode='list', length=numFiles)
	
	for (i in 1:numFiles)
	{
		path <- paste(root, curFiles[i], sep="/", collapse="")
		resultList[[i]] <- readDIANA(path, lSpecies)
	}
	
	return(resultList)
}