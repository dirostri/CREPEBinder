# CREPE: A New ShinyR package for Transcription Factor Cataloguing 
This the CREPE GitHub repository. 

### Launch on Binder
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dirostri/CREPE/tree/main/HEAD)

### Running CREPE Locally
1. Open the App.R script in RStudio
2. Click on the Run App icon located at the top of the app
3. An html window will pop up with CREPE running

### Tutorial: Minimal Reproducible Example (MRE)
We provide tutorial files under the MinimalReproducibleExample/ directory.

#### TF Cataloguing
This function runs a domain analysis to find protein sequences that posess TF DNA-binding Domains (DBDs), and parses these results to provide a TF catalogue. 

1. Run CREPE (see section above)
2. Open the sidebar menu (Three Parallel Bars Button in the top left)
3. Under 'Choose Analysis To Run' click on TF Cataloguing. 
4. Under 'Protein Fasta' click on the 'Browse' button to search for the MRE fasta file
```
CREPE/MinimalReproducibleExample/TF_Cataloguing/Human_MRE.fasta
```
5. In the 'Species Name' text box write Human, or any other desired identifier. 
6. An indicator will pop-up in the bottom right of the screen indicating progress; may take several minutes to complete. 
7. A summary figure of the TF catalogue will be displayed after the analysis is complete. 
8. Click on the 'Download TF Catalogue' to download cataloguing results. 

#### TF Annotation
This function parses gene trees (in Newick format) to assign putative TFs the name of its nearest neighbor. The intention is to provide a name to TFs originating from non-model organisms. 

1. Run CREPE (see section above)
2. Open the sidebar menu (Three Parallel Bars Button in the top left)
3. Under 'Choose Analysis To Run' click on TF Annotation.
4. Under 'Phylogenetic Trees' click on the 'Select Folder' button to search for the MRE trees
```
CREPE/MinimalReproducibleExample/TF_Annotation/Trees/
```
5. Under Metadata click on the 'Browse' button to search for the MRE metadata file.
```
CREPE/MinimalReproducibleExample/TF_Annotation/metadata.csv
```
6. An indicator will pop-up in the bottom right of the screen indicating progress; may take several minutes to complete. 
7.  A table displaying the results will appear after analysis is complete. 
8.  To download the results, choose the species to map to under the 'Pick Mapping to Download' menu; the current options are human, fly or nearest species. Afterwards, click on the 'Download TF Annotation' button.

#### Custom Analysis (Advanced Users)
This function is for those who would like to perform the domain analysis separately using other software, or using different DBD models. 

1. Run CREPE (see section above)
2. Open the sidebar menu (Three Parallel Bars Button in the top left)
3. Under 'Choose Analysis To Run' click on Custom.
4. Under 'Protein Fasta' click on the 'Browse' button to search for the MRE fasta file
```
CREPE/MinimalReproducibleExample/TF_Cataloguing/Human_MRE.fasta
```
6. Under 'Tabular Domain Outfile' click on the 'Browse' button to search for the tabular domain analysis run separately by the user. In this case we added the 'Myc_N' domain as an example. 
**PLEASE NOTE: This is not a sequence specific TF DBD, it is only shown here as an example. **
```
CREPE/MinimalReproducibleExample/Custom/TabDomainOutfile.tsv
```
7. In the 'Custom Accension' text box write the accession for the additional DBDs used to make the Tabular Domain Outfile. In this case the accession is
```
PF01056
```
8. In the 'Species Name' text box write Human, or any other desired identifier. 
9. An indicator will pop-up in the bottom right of the screen indicating progress; may take several minutes to complete. 
10. A summary figure of the TF catalogue will be displayed after the analysis is complete.
11. Click on the 'Download TF Catalogue' to download cataloguing results.

##### Notes on making the Tabular Domain Outfile
There are several tools available to perform this task. However, to be compatible with CREPE the Tabular Domain Outfile must contain the following information in order: 

| DBD_Name        | Accession     | Query_Name   | Domain_Start   | Domain_Start   |
| ---------------:|:-------------:|-------------:|---------------:|---------------:|

Additionally, the Accensions of the DBD models referenced in CREPE originate from the PFAM database. If the analysis is performed using DBD models from other databases they will not be identified by CREPE. 







