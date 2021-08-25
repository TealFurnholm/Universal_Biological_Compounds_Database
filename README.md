# Universal_Biological_Functions_Database
This script compiles and links biological reactions and compound information from across multiple public repositories. 
This is the second of a series of scripts that create/curate vital reference databases for microbiome studies.
    1. Universal Taxonomy: https://github.com/TealFurnholm/Universal-Taxonomy-Database
    2. Universal Functions: this
    3. Universal Protein Alignment Database: https://github.com/TealFurnholm/Universal_Microbiomics_Alignment_Database
    4. Non-coding RNA Alignment Database: https://github.com/TealFurnholm/Fix_RNACentral_Taxonomy

## How Universal?
These databases span all kingdoms of life. The databases allow the simultaneous identification of microbial community phylogeny and functions. All the biological molecules/compounds and their data are linked to their enzymes and transporters to map the flow of metabolites in microbe-microbe or microbe-host interactions. 
The databases are used for (meta)transcriptomics, (meta)proteomics, metabolomics, metagenomics, and for novel binning and MAG quality control software I've created. This way many types of data can be combined into a secondary analysis, with taxonomy and functions being directly linked. It also covers both the protein and the non-coding fraction of sequencing. 

## 1 Manual Input: BioCyc
For the most part, this script is automated - except the BioCyc collection.
First you'll need to request a license to download their flat files.
Then you'll need to download their biocyc - all tiers archive file: (you have to enter a user name and password)
<br> wget -N --no-check-certificate --user= --password=  http://brg-files.ai.sri.com/public/dist/biocyc-flatfiles.tar.gz
* this bloody file is huge (takes days to download) and has repeatedly been borked (BioCyc archiving issue)
* if the server gives out or your download is disrupted, add a -c to the wget command  

## Automation
This script automatically pulls in all the remaining data it needs off the various public repositories.
Upshot, you shouldn't have to do anything but run it with the options you want.
Downside - those repositories may change their links/structure, so you may have to fix the links in the code.
I put as many "wget" commands as I could near the top of the script so you could mod them if need be.

## How to Use
The wiki (https://github.com/TealFurnholm/Universal_Biological_Functions_Database/wiki) explains how these scripts create:
1. A compilation of function ids and their names
2. A "reactions database" which 
    * maps reaction ids across kegg, uniprot, EC, rhea, and Biocyc
    * has each databases respective reactants, products, and direction
3. A "compounds database" which has
    * each reactant/product compound id in the reactions database
    * the molecular mass, formula, and charge

These database files may or maynot be directly useful - they are used to annotate all the proteins in the UniProt database.
HOPEFULLY, I will have access to some metabolomics data soon so I can build upon the compound database, creating scripts to integrate the metabolome and NGS data.




