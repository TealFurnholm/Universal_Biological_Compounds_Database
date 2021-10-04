## Universal_Biological_Compounds_And_Reactions_Database
This script compiles and links biological reactions and compound information from across multiple public repositories. 
This is #2 of a series of pipelines that comprise this Universal Reference for various 'omics data analysis.
    <br>1. Universal Taxonomy Database: [found here](https://github.com/TealFurnholm/Universal-Taxonomy-Database)
    <br>2. Universal Compounds Database: *this repository*
    <br>3. Universal Reactions Database: *this repository*
    <br>4. Universal Protein Alignment Database: [found here](https://github.com/TealFurnholm/Universal_Microbiomics_Alignment_Database)
    <br>5. Universal ncRNA Alignment Database: [found here](https://github.com/TealFurnholm/Fix_RNACentral_Taxonomy)

## How Universal?
These databases span all kingdoms of life. The databases allow the simultaneous identification of microbial community phylogeny and functions. All the biological molecules/compounds and their data are linked to their enzymes and transporters to map the flow of metabolites in microbe-microbe or microbe-host interactions. 
The databases are used for (meta)transcriptomics, (meta)proteomics, metabolomics, metagenomics, and for novel binning and MAG quality control software I've created. This way many types of data can be combined into a secondary analysis, with taxonomy and functions being directly linked. It also covers both the protein and the non-coding fraction of sequencing. 

## One Manual Input: BioCyc
For the most part, this script is automated - except the BioCyc collection. You'll need to request a license to download their flat files. They will give you a username, password, and link to the download site - you'll be asked to asked to provide these for downloading the files. 

## Automation
This script automatically downloads data it needs off the various public repositories.
 - Upshot - you shouldn't have to do anything but run it
 - Downside - those repositories may change their links, so you may have to fix the links in the code
<p></p>I put as many "wget" commands as I could near the top of the script so you could easily mod them if need be.

## How to Use
There are several primary and secondary microbiome analysis pipelines on my main GitHub page that use these databases:
* RNAseq/Metatranscriptome Analysis [here](https://github.com/TealFurnholm/Strain-Level_Metatranscriptome_Analysis)
* Metagenome Primary Analysis [here](https://github.com/TealFurnholm/Strain-Level_Metagenome_Analysis)
* Metabolomics Analysis - TBD; for now you can directly link the metabolomics output to compounds->reactions->proteins/oranisms using the Functional and Protein databases
* Strain-level Metagenome Binning [here](https://github.com/TealFurnholm/Community-Based-Metagenome-Binning)
* Secondary Functional Analysis and Visualization [here](https://github.com/TealFurnholm/Meta-omics_Functional_Analysis)
    
# Get started: [wiki](https://github.com/TealFurnholm/Universal_Biological_Functions_Database/wiki)
