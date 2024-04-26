This repository contains the flow data and scripts used to generate Fig. 4, Supplementary Fig. S7 and Table S5 in 'Promoter selection reduces phenotypic heterogeneity and improves bioproduction' (Patel et al.)

Each folder contains the .csv files generated after data cleanup in FlowJo - the analysis of which is below.

There are 5 biological replicates named either: 
- 140623_r1
- 140623_r2
- 210623_r1
- 210623_r2
- 300623_r1

Use each folder with the R script containing the same replicate name e.g. 140623_r1 folders for use with '140623_rep1.R'

The media condition is also in each folder name: 
- ypd_exp = YPD exponential phase
- ypd_sta = YPD stationary phase
- ypdx_exp = YPD10 exponential phase
- ypdx_sta = YPD10 stationary phase
- ynb_exp = YNB exponential phase
- ynb_sta = YNB stationary phase 


--------------------------------------------------------------------------------------------------------------
The R scripts:
The R scripts contain the Dip test, mixed modelling analyses (for subpopulation identification) and plots generated following flow cytometry of the model validation strains. Run the functions and then start from the 'mixedsort' function under each media subsection to generate the analyses and plots from the paper.

Plate layouts and .fcs files can be made available upon request. Specifics about how the scripts work are outlined within the comments.

-------------------------------------------------------------------------------------------------------------
FlowJo analysis to generate .csv files in the folders

Control strains:

- GFP-/RFP+ = BP-HK = A9
- GFP-/RFP- = BP- = A10
- GFP+/RFP- = 897- = A11
- negative sample control = BP with PI = A8

For all sample groups in the FlowJo workspace, wells 8-11 are the controls for each group. E.g. for another media condition, wells B8-B11 or wells C8-C11 may be the control wells. This will be specified within the scripts.

FlowJo Workspace Analysis Workflow:

1) singlets from BP- (no PI) (A10)
 - change axes to log
 - apply to all in group

2) comp from A9-A11
 - apply to all in group

3) Live/dead gating from BP-HK (A9) by comparing to BP with PI (A8)
 - Live gate made as an offshoot of the singlets gate so 
   it is indented underneath it. The Live gate is applied to all samples in the same media condition and growth phase group. 

4) sfGFP intensity files are exported as csvs showing BL1-H intensity for each cell

