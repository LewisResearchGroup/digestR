# DigestR

**Author**: Dimitri Desmonts de Lamache

## Overview
DigestR is an open-source software developed for the R statistical language environment. The DigestR package is based on rNMR and contains a collection of tools for visualizing and analyzing protein digestion. 
Users can interact interact with DigestR in two major ways: via point and click graphical user interfaces (GUIs) or by entering command directly in the R console. 
This guide is intended to give an overview of DigestR's functions.


## Prerequisites

	install.package('biomaRt')
	library(biomarRt)



## Run the following lines in the R console

	library(Rcpp)
	library(tcltk)
	setwd('C:/Users/dimit/OneDrive/Bureau/DigestR')
 
	source('2017_09_16_rNMR_Travis_DD_Edits_Final.R')
	source('digestR_Code_Dimitri_final.R')
	
	

	rm(list = ls())
	devtools::load_all("C:\\Users\\soere\\workspace\\digestR")


## Supported file formats 

DigestR supports Mascot (.csv) and MaxQuant (.mzXML) generated files. These files may be converted to .dcf files using DigestR file conversion functions pm() and im() respectively. 

## GUI Functions Documentation

This document provides details of several Graphical User Interface (GUI) functions implemented in the R programming language using the Tcl/Tk toolkit.

### Converting files. 

#### 1. Prepare Mascot files: rd()
Users can natively import in DigestR files generated from mascot (.csv) and from MaxQuant (.mzXML). However, these files need to be pre-processed to be analyzed by DigestR. To prepare the mascot files for digestR, call the function rd() or click on ‘Prepare Mascot file’ under the ‘Manipulate csv’ section. Select the file(s) to be processed. The software will automatically prepare and save the new file.

#### 2. Unique peptides: up()
DigestR allows users to compare two mascot files and identify peptides that are unique to a specific experimental group. The up() function opens a main GUI window with options to select two CSV files (Query and Experimental) and find unique peptides in the Experimental file that are not present in the Query file. The function then writes the peptides unique to the experimental group to a CSV file named with the Experimental filename appended with "_Unique".

#### 3. Process mascot files: pm() 
To create "digestion" maps, peptides identified by Mascot or MaxQuant need to be mapped to their proteomic location. First, the user needs to select a proteome to align peptides against. DigestR natively includes five proteomes, such as Human, Bos Taurus, and Plasmodium falciparum. Users can import new proteomes using the function xxx(). After proteome selection, users can align Mascot or MaxQuant identified peptides along the proteome by selecting single or multiple Mascot files. The alignment generates a "coincidence" or "digestion" map that users can interact with.

### 1. Remove Duplicates Function
The `rd()` function creates a main GUI window, which provides an option to remove duplicate entries from a CSV file. This function first creates a "Remove Duplicates" button. When clicked, the button opens a file dialog for the user to choose a CSV file. After the file is selected, the removeDuplicates function is called to remove duplicate entries based on the 'pep_seq' column. The resulting file is then saved to a location specified by the user.

### 2. Unique Peptides Function
The `up()` function creates a main GUI window, which provides options to select two CSV files (Query and Experimental) and find unique peptides in the Experimental file that are not present in the Query file. The function then writes the unique peptides to a CSV file named with the Experimental filename appended with "_Unique".

### 3. Density Plot Function
The `pd()` function creates a GUI that allows users to create density plots from grouped CSV files located in a specified directory. Users can choose from three types of plots: Overlay, Ridges, and Colored Ridges.

### 4. Cut Sites Distribution
The `csp()` function creates a GUI that allows users to create bar plots from grouped CSV files located in a specified directory. Users can choose from three types of plots: Nter, Cter, and Combined.

### 5. Plotting Cut Sites Function
The `cs()` function creates a GUI window that allows users to select a protease from a dropdown menu and input a color for the lines in the plot. Upon clicking the "Plot" button, the function calls the plotCutSite function, passing the selected protease and color as arguments, to generate a plot of cut sites.

### 6. Testing the Functions
For each function, to test the GUI, call the function without any arguments. For example, to test the up function, simply call up().
