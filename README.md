[![Install package](https://github.com/LewisResearchGroup/digestR/actions/workflows/R-package-windows.yml/badge.svg)](https://github.com/LewisResearchGroup/digestR/actions/workflows/R-package-windows.yml)
[![Install via GitHub](https://github.com/LewisResearchGroup/digestR/actions/workflows/R%20package%20on%20Windows%20-%20from%20GitHub.yml/badge.svg?branch=main)](https://github.com/LewisResearchGroup/digestR/actions/workflows/R%20package%20on%20Windows%20-%20from%20GitHub.yml)

# DigestR
 
**Author**: Dimitri Desmonts de Lamache

## Overview
DigestR is an open-source software developed for the R statistical language environment. The DigestR package is based on rNMR and contains a collection of tools for visualizing and analyzing protein digestion. 
Users can interact with DigestR in two major ways: via point and click graphical user interfaces (GUIs) or by entering command directly in the R console. 
This guide is intended to give an overview of DigestR's functions.

## How to install digestR package from GitHub
DigestR was tested with `R v4.3.1`, and `R v4.4.1`.

We recommend using R '4.3' or later versions. In case of dependencies issues, use the renv package to reproduce the exact environment used.
See section 'Reproducing R environment'

### Prerequisites

Before installing DigestR, ensure that `BiocManager` is installed to manage Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")
```
   
### Step 1: Install and load the devtools package
To install a R package, start by installing the devtools package. The best way to do this is from CRAN, by typing:

    install.packages("devtools")
    library(devtools)

### Step 2: Install the package of interest from GitHub

Install the the digestR package directly from GitHub
    
    install_github("LewisResearchGroup/digestR")

### Step 3: Load the package

    library(digestR)

## Reproducing the R Environment

This project uses `renv` to manage package dependencies. To reproduce the exact environment used:
### Step 1: Install and Load renv
   ```sh
install.packages("renv")
library(renv)
   ```
### Step 2: Clone the Repository
You can clone the repository using your system's terminal. Run the following command:
```sh
git clone https://github.com/LewisResearchGroup/DigestR.git
cd digestR
   ```
Alternatively, if you want to clone the repository from within R, you can do so using system commands:
```sh
system("git clone https://github.com/LewisResearchGroup/DigestR.git")
   ```
### Step 3: Navigate to the Project Directory
Change the working directory to where you cloned the repository:
```sh
setwd("path/to/DigestR")
   ```
Replace "path/to/DigestR" with the actual path where the repository was cloned.

### Setp 4: Initialize renv
```sh
renv::init()
   ```
### Step 5: Restore the Project Environment
```sh
renv::restore()
   ```
### Step 6: Install the Package
```sh
devtools::install()
   ```
### Step 7: load the Package
```sh
library(DigestR)
   ```

## GUI Functions Documentation

This document provides details of several Graphical User Interface (GUI) functions implemented in the R programming language using the Tcl/Tk toolkit.

### Supported file formats 

DigestR supports Mascot (.csv) generated files. These files may be converted to .dcf files using DigestR file conversion functions pm(). 

### Converting files. 

Required files: .csv

To read a Mascot file, DigestR requires the following information: protein accession number (prot_acc), protein description (prot_desc), and the peptide sequence (pep_seq).
When exporting data from Mascot you should uncheck the "Include search information" box to prevent the inclusion of the search information header. This header will prevent digestR from reading the CSV file. 
When importing the CSV file in DigestR, it should resemble the one below:

![](https://github.com/LewisResearchGroup/digestR/blob/main/Images/Mascot%20file%20format.png)

By default, Mascot creates a header, this 3 line header is required for the .csv preprocessing (see Prepare Mascot files). 

An example data file can be found here: https://github.com/LewisResearchGroup/digestR/blob/main/Example%20Files/Data_Example.csv

### Generating Proteome: gp() 
The generate_proteome function streamlines the process of accessing and downloading protein data from Ensembl BioMart, facilitating the creation of proteomes for comparison against experimental peptides. To generate a new proteome, users begin by selecting their desired Biomart library, using the dropdown menu – options include "genes" or "ensembl," with "genes" being the default value.

Following this, users input a search pattern to explore datasets within the BiomaRt database (e.g., "sapiens" or "taurus"). Upon clicking the "Search Datasets" button, the function connects to the BiomaRt servers and retrieves datasets matching the provided pattern. The outcomes are displayed in the "Dataset Results" listbox, showing the dataset names, descriptions, and versions. Double-clicking on a result selects the dataset for further processing.

Subsequently, users have the option to download data from specific chromosomes by entering a list of chromosomes separated by commas (e.g., "1, 2" or "10, X"). Alternatively, leaving this field empty will result in data retrieval for all chromosomes. Finally, users initiate the data retrieval process by clicking the "Download proteome" button. The function then acquires protein data from the selected dataset and chromosomes, saving the findings as CSV files within a dedicated "data/proteomes" subfolder. 

![](https://github.com/LewisResearchGroup/DigestR/blob/main/Images/generate%20proteome.png)

Warning: Generating proteomes can take several minutes to several hours depending on the size of the data to be retrieved from Biomart. 
For convernience, some proteomes have already been generated and can be found here: 
https://github.com/LewisResearchGroup/digestR/blob/main/Example%20Files/Human_proteome.csv
https://github.com/LewisResearchGroup/digestR/blob/main/Example%20Files/Pfal_3D7_proteome.csv

### Manipulating .csv files

#### 1. Process mascot files: pm() 
To create "digestion" maps, peptides identified by Mascot or MaxQuant need to be mapped to their proteomic location. First, the user needs to select a proteome to align peptides against (see Generate Proteome). DigestR will automatically detect and utilize all proteomes located within the "data/proteomes" subfolder. Users can also import their own proteomes into this subfolder. After proteome selection, users can align Mascot identified peptides along the selected proteome from a single or multiple files. These alignments generate "coincidence" or "digestion" maps that users can interact with.

![](https://github.com/LewisResearchGroup/DigestR/blob/main/Images/Process%20Mascot%20GUI.png)

#### 2. Cleavage site specificity: csd()
The csd() function allows users to plot amino acid distributions at C-terminus or N-terminus to track changes in cut site representation/specificity between groups. This function allows users to select a file to generate either a logo plots of the P4-P4' positions or bar plots at the P1 (Nterminus) or P1' (Cterminus) position. To identify cleavage sites of biological significance, it is possible to normalize the distribution with a specific amino acid sequence. Users can directly import a protein sequence in the appropriate box. The function then calculates the representation frequency for amino acid within the protein sequence to normalize the experimental amino acid cut-site distributions. 

![](https://github.com/LewisResearchGroup/digestR/blob/main/Images/LogoPlot_function.png)

#### 3. Peptide length distribution: pd()
Defect in proteolytic activity might have an impact on digested peptide length. Therefore, DigestR was developed to calculate and plot peptide length distributions in amino acids using the pd() command. Users can select a folder or subfolder and process all CSV files in that directory.  Files can be selected directly in the loaded files box. If no files are selected, all files will be used to generate the density plots. At least two files need to be imported in order to generate Venn diagrams Files will be grouped depending on the second string of the filename. Three types of density plots from grouped CSV files can be chosen by the user: Overlay, Ridges, and Colored Ridges.

![](https://github.com/LewisResearchGroup/digestR/blob/main/Images/peptide_distribution_function.png)

#### 4. Venn Diagrams: vd()
DigetR also allows for the creation of Venn diagrams in order to analyze peptide overlaps between groups. The vd function allows for users to import files contained in a specific folder and generate Venn diagram. Files can be selected directly in the loaded files box. If no files are selected, all files will be used to generate the Venn diagram . At least two files need to be imported in order to generate Venn diagrams

![](https://github.com/LewisResearchGroup/digestR/blob/main/Images/VennDiagram_function.png)

### Viewing and interacting with Digestion maps. 

Required files: .dcf

#### 1. Opening digestion maps: fo() / fs()
To open a "digestion" map in DigestR, either select "Open/Close Files" from the File menu or use the commands fo() or fs() in the R console. If multiple files have been opened, only the most recently opened spectrum will appear in the main plot window. To switch to another spectrum, double-click on a file name within the GUI. To close one or more files, select the desired files from the table and then press the "Close file" button.

#### 2. Manipulate dcf files: mf()
The mf() function in DigestR allow users to perform various mathematical operations with dcf files, facilitating comprehensive data manipulation and analysis. With mf(), users can add, substract, multiply, merge and divide, the data contained in multiple dcf files. This functionality allows users to perform mathematical operations tailored to their specific research needs, streamlining data processing and enhancing the overall analytical capabilities of DigestR.

![](https://github.com/LewisResearchGroup/DigestR/blob/main/Images/Manipulate%20files%20GUI.png)

#### 3. Plot settings: ct()
The ct() function allows users to interact with the "digestion" map directly through the graphical interface. Users can display the "digestion" map either at a proteome or protein level. By default, the full proteome view is displayed. 

![](https://github.com/LewisResearchGroup/DigestR/blob/main/Images/CT%20GUI.png)


To display the digestion map of a specific protein, click on "Display Single Gene" and enter the protein name in capital letters (e.g., HBA1, BSA). It is also possible to display individual peptides and peptide endpoints. Additionally, variance can be displayed if the selected dataset is the mean of multiple datasets (see mf() function).

![](https://github.com/LewisResearchGroup/DigestR/blob/main/Images/CT_GENE%20GUI.png)

#### 4. Plot colors: co()
The plot color function allows users to easily manipulate the plot colors. To open the plot color GUI, enter the command co(). Color preferences can be applied to multiple spectra simultaneously by selecting names from the files list. Plot color options for the selected files may be configured individually using the buttons provided on the right side of the GUI. The "Axes" button changes the color of the x and y axes, "BG" changes the background color, and "Peak labels" changes the label color of identified peaks.

![](https://github.com/LewisResearchGroup/DigestR/blob/main/Images/Color%20GUI.png)

#### 5. Overlays: ol()
DigestR allows multiple "digestion" maps to be displayed concurrently on a single plot through the command ol(). To add or remove loaded files, select the digestion maps to overlay and click the "add" or "remove" buttons. The order of overlaid maps in the main plot window is taken directly from the order of digestion maps appearing in the overlays list box. Individual files can be assigned their own colors. The plot legend will be automatically generated, but it can be suppressed by unchecking the "Display names of the overlay spectrum on the plot" option. Similarly, the path of "digestion" maps can be suppressed by checking the corresponding checkbox.

![](https://github.com/LewisResearchGroup/digestR/blob/main/Images/Overlay_function.png)

#### 6. Display protease cut sites: cs()
Users can overlay known protease cut sites onto the "digestion" map(s) using the cs() command. It is important to note that this function requires a CSV file containing the names of the proteases and their respective cleavage sites. An example CSV file can be found here: https://github.com/LewisResearchGroup/digestR/blob/main/tests/Proteasecutsiteslist.csv

| protease | abv | cutsites | X | X.1 | X.2 | X.3 | X.4 |
|---|---|---|---|---|---|---|--- |
| Trypsin | Trp | K | R |  |  |  |  |
| ChymotrypsinHA | ChymH | F | Y |  |  |  |  |
| ChymotrypsinLA | ChymL | F | Y | L | M | H |  |
| Pepsin1.3 | Pep1.3 | F | L |  |  |  |  |
| Pepsin2 | Pep2 | F | L | W | Y |  |  |
| Thermolysin  | Thermo | A | F | I | L | M | V |
| GluC | GluC | E |  |  |  |  |  |
| Proteinase K | ProtK | H | K | R | P |  |  |
| LysC | LysC | K |  |  |  |  |  |

Once the CSV file is loaded, users have to option to select a specific protease and choose the color for the lines representing the cleavage sites on the map. 

![](https://github.com/LewisResearchGroup/DigestR/blob/main/Images/cs%20GUI%20-%20Copy.png)

#### 7. Zoom: zm()
DigestR includes various zooming and scrolling commands, accessible through the zoom GUI by selecting "Zoom" from the View menu or using the command zm(). Digestion maps can be navigated using the arrow pad provided in the zoom GUI or by using the five distinct zoom functions called by the buttons provided on the right side of the zoom GUI. Many of these functions are iterative and must be exited by right-clicking in the main plot window.

![](https://github.com/LewisResearchGroup/DigestR/blob/main/Images/Zm%20GUI.png)

#### 8. Gene labelling: gl()
The gl() function allows users to override the threshold at which proteins are labeled when viewing data on the proteome-wide level. By lowering the default value, more peptides will be labeled.

![](https://github.com/LewisResearchGroup/DigestR/blob/main/Images/gl%20GUI.png)







