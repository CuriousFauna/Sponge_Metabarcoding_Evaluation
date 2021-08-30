# DNA metabarcoding evaluation of Sponges
This is a repository for code and data used in the analysis of the manuscript entitled: "Sponging up diversity: Evaluating metabarcoding performance for a taxonomically challenging phylum within a complex cryptobenthic community" - https://onlinelibrary.wiley.com/doi/full/10.1002/edn3.163. 

## Fastq Sequences
Fastq files were submited to NCBI through Genome -https://geome-db.org/workbench/project-overview. The NCBI Sequence Read Archive (SRA) numbers for this project are SRS7105074-SRS7105095 and can be obtained here: https://www.ncbi.nlm.nih.gov/sra.
The demultiplexed fastq files are named ArmsMeso1 thru 22. There are two files, one forward and one reverse for each treatment. For example: ArmsMeso1.1 = forward and ArmsMeso1.2 = reverse. 

Sequences were processed through the JAMP pipeline - https://github.com/VascoElbrecht/JAMP. See manuscript supplement for details.
The JAMP pipeline output is the Raw_MOTU_Table.csv file located in the Data Folder.

## Data and R Code
All data and R code used to obtain the results are presented in the Data and R Code folders. R packages needed to process the data are provided in the respective heading of each R file. See methods for details on sequence annotations and conducted statistical analyses. Figures and Tables presented in the manuscript are highighted within the R file SpongeComparison.R.



