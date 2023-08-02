# DigestR

**Author**: Dimitri Desmonts de Lamache

## Overview
DigestR is an open-source software developed for the R statistical language environment. The DigestR package is based on rNMR and contains a collection of tools for visualizing and analyzing protein digestion. 
Users can interact with DigestR in two major ways: via point and click graphical user interfaces (GUIs) or by entering command directly in the R console. 
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

## GUI Functions Documentation

This document provides details of several Graphical User Interface (GUI) functions implemented in the R programming language using the Tcl/Tk toolkit.

### Supported file formats 

DigestR supports Mascot (.csv) and MaxQuant (.mzXML) generated files. These files may be converted to .dcf files using DigestR file conversion functions pm() and im() respectively. 

### Converting files. 

Required files: .csv and .mzXML

#### 1. Prepare Mascot files: rd()
Users can natively import in DigestR files generated from mascot (.csv) and from MaxQuant (.mzXML). However, these files need to be pre-processed to be analyzed by DigestR. To prepare the mascot files for digestR, call the function rd() or click on ‘Prepare Mascot file’ under the ‘Manipulate csv’ section. Select the file(s) to be preprocessed. The software will automatically prepare and save the new file.

<img width="964" alt="Remove Duplicate GUI" src="https://github.com/LewisResearchGroup/DigestR/assets/139395028/e19ac9bb-5bc5-4f22-96f7-73cae0cfed61">


#### 2. Unique peptides: up()
DigestR allows users to compare two mascot files and identify peptides that are unique to a specific experimental group. The up() function opens a main GUI window with options to select two CSV files (Query and Experimental) and find unique peptides in the Experimental file that are not present in the Query file. The function then writes the peptides unique to the experimental group to a CSV file named with the Experimental filename appended with "_Unique".

<img width="959" alt="Unique Peptides GUI" src="https://github.com/LewisResearchGroup/DigestR/assets/139395028/afbd9554-1ac2-4492-9434-edf70c846e41">

#### 3. Process mascot files: pm() 
To create "digestion" maps, peptides identified by Mascot or MaxQuant need to be mapped to their proteomic location. First, the user needs to select a proteome to align peptides against. DigestR natively includes five proteomes, such as Human, Bos Taurus, and Plasmodium falciparum. Users can import new proteomes using the function xxx(). After proteome selection, users can align Mascot or MaxQuant identified peptides along the proteome by selecting single or multiple Mascot files. The alignment generates a "coincidence" or "digestion" map that users can interact with.

<img width="958" alt="Process Mascot GUI" src="https://github.com/LewisResearchGroup/DigestR/assets/139395028/67bf5d8f-cfbb-4d18-aac4-b7e511838809">


### Viewing and interacting with Digestion maps. 

Required files: .dcf

#### 1. Opening digestion maps: fo() / fs()
To open a "digestion" map in DigestR, either select "Open/Close Files" from the File menu or use the commands fo() or fs() in the R console. If multiple files have been opened, only the most recently opened spectrum will appear in the main plot window. To switch to another spectrum, double-click on a file name within the GUI. To close one or more files, select the desired files from the table and then press the "Close file" button.

#### 2. Manipulate dcf files: mf()
The mf() function in DigestR allow users to perform various mathematical operations with dcf files, facilitating comprehensive data manipulation and analysis. With mf(), users can add, substract, multiply, merge and divide, the data contained in multiple dcf files. This functionality allows users to perform mathematical operations tailored to their specific research needs, streamlining data processing and enhancing the overall analytical capabilities of DigestR.

<img width="958" alt="Manipulate files GUI" src="https://github.com/LewisResearchGroup/DigestR/assets/139395028/54e99813-9c85-400e-9eb8-aefaa27b0ee3">

#### 3. Plot settings: ct()
The ct() function allows users to interact with the "digestion" map directly through the graphical interface. Users can display the "digestion" map either at a proteome or protein level. By default, the full proteome view is displayed. 

<img width="960" alt="CT GUI" src="https://github.com/LewisResearchGroup/DigestR/assets/139395028/ca1cbb63-0a06-428b-801d-d7775805ba9c">



To display the digestion map of a specific protein, click on "Display Single Gene" and enter the protein name in capital letters (e.g., HBA1, BSA). It is also possible to display individual peptides and peptide endpoints. Additionally, variance can be displayed if the selected dataset is the mean of multiple datasets (see mf() function).

<img width="959" alt="CT_GENE GUI" src="https://github.com/LewisResearchGroup/DigestR/assets/139395028/8892b896-2c31-4be4-9cb9-4fdf9b44266d">

#### 4. Plot colors: co()
The plot color function allows users to easily manipulate the plot colors. To open the plot color GUI, enter the command co(). Color preferences can be applied to multiple spectra simultaneously by selecting names from the files list. Plot color options for the selected files may be configured individually using the buttons provided on the right side of the GUI. The "Axes" button changes the color of the x and y axes, "BG" changes the background color, and "Peak labels" changes the label color of identified peaks.

<img width="959" alt="Color GUI" src="https://github.com/LewisResearchGroup/DigestR/assets/139395028/9487351a-0567-488d-a477-0be519267493">

#### 5. Overlays: ol()
DigestR allows multiple "digestion" maps to be displayed concurrently on a single plot through the command ol(). To add or remove loaded files, select the digestion maps to overlay and click the "add" or "remove" buttons. The order of overlaid maps in the main plot window is taken directly from the order of digestion maps appearing in the overlays list box. Individual files can be assigned their own colors. The plot legend will be automatically generated, but it can be suppressed by unchecking the "Display names of the overlay spectrum on the plot" option. Similarly, the path of "digestion" maps can be suppressed by checking the corresponding checkbox.

#### 6. Display protease cut sites: cs()
Users can overlay known protease cut sites onto the "digestion" map(s) using the cs() command. The GUI will prompt the user to select the protease and choose the color of the lines representing cleavage sites. Note that this function requires a CSV file containing the names of the proteases and their cleavage sites. An example can be found here: XXXX

<img width="959" alt="cs GUI" src="https://github.com/LewisResearchGroup/DigestR/assets/139395028/2baa3d3e-2f37-49c7-a07a-615f34ebed33">

#### 7. Zoom: zm()
DigestR includes various zooming and scrolling commands, accessible through the zoom GUI by selecting "Zoom" from the View menu or using the command zm(). Digestion maps can be navigated using the arrow pad provided in the zoom GUI or by using the five distinct zoom functions called by the buttons provided on the right side of the zoom GUI. Many of these functions are iterative and must be exited by right-clicking in the main plot window.

<img width="960" alt="Zm GUI" src="https://github.com/LewisResearchGroup/DigestR/assets/139395028/a81055a2-8e06-46a0-bc1c-c5287e86a2df">

#### 8. Gene labelling: gl()
The gl() function allows users to override the threshold at which proteins are labeled when viewing data on the proteome-wide level. By lowering the default value, more peptides will be labeled.

<img width="959" alt="gl GUI" src="https://github.com/LewisResearchGroup/DigestR/assets/139395028/f639c95d-8dfb-4494-9c9e-117c90d90bb9">

### Peptides analyses 

Required files: .csv

#### 1. Peptide length distribution: pd()
Defect in proteolytic activity might have an impact on digested peptide length. Therefore, DigestR was developed to calculate and plot peptide length distributions in amino acids using the pd() command. Users can select a folder or subfolder and process all CSV files in that directory. Files will be grouped depending on the second string of the filename. Three types of density plots from grouped CSV files can be chosen by the user: Overlay, Ridges, and Colored Ridges.

<img width="959" alt="Density Plot GUI" src="https://github.com/LewisResearchGroup/DigestR/assets/139395028/25dff9dc-71da-43f4-87d9-882834331ac7">

#### 2. Cleavage site specificity: csp()
The csp() function allows users to plot amino acid distributions at C-terminus or N-terminus to track changes in cut site representation/specificity between groups. Similar to pd(), this function allows users to select a folder or subfolder and process all CSV files in that directory. Files will be grouped depending on the second string of the filename. The user can select to plot either C-terminus or N-terminus or both distributions.

<img width="959" alt="CSPGUI final" src="https://github.com/LewisResearchGroup/DigestR/assets/139395028/5576a77a-569b-4943-81e8-4bb757ddb2dd">




