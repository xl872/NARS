# NARS
23Na fMRI
README: Reproducible Analysis Workflows for NARS-fMRI

This repository contains the MATLAB and AFNI scripts used to generate the primary results reported in the NARS-fMRI study (Yu et al. Biorxiv, 2026, https://doi.org/10.64898/2026.02.09.704765)
The workflows are organized to enable reproduction of the major figures, statistical analyses, and validation experiments  described in the manuscript. Each analysis script is explicitly linked to manuscript figure panels to facilitate transparent reproduction and editorial assessment.
Manuscript Linkage Table
Script	Primary Output	Linked Manuscript Panels
QC_phase_$date.m	Phase stability and SNR metrics	Fig. 1
GLM_NARS_$date.m	ß / t-statistic maps	Fig. 23; Supp Fig. 45
Randomized_$date.m	Resampling stability analysis	Supp Fig. 78
NARS_Glu_$date.m	Peak-response scatter analysis	Fig. 4
T2_all.m / T1_all.m	Relaxometry maps	Supp Fig. 1

All scripts were developed and tested using MATLAB (2025a) and AFNI (Version 24.0.09). Minor path adjustments may be required depending on the local computing environment.

1. Software and Dependencies
- MATLAB (Version 2025a)
  - Signal Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - Image Processing Toolbox
- AFNI (Version 24.0.09)
Custom helper functions are located in /helper/ and must be added to the MATLAB path before execution.

2. Repository Structure
/NARS_pipeline/
/Glu_photometry/
/T1T2_mapping/
/helper/ (custom functions)

3. Example datasets:
- Mouse_0220(directory)
- Rat_0227(directory)
-Mouse_0522_Na_Glu(directory)
-T1T2_mapping(directory)

4. Analysis Workflow Overview (Recommended Execution Order)
?Quality Control and Trial Assessment
?GLM Functional Mapping
?Resampling-Based Stability Analysis
?NARSGlu Cross-Modal Analysis
?23Na T1/T2* Characterization

5. NARS-fMRI Functional Analysis: /NARS_pipeline/
a.  Trial Quality Control
Script: QC_phase_$date.m
Function:
- Reads individual trials
- Computes mean phase offsets and variability
- Estimates SNR as a function of trial averaging
Manuscript linkage:
- Generates the phase-stability and trial-averaging SNR analyses underlying Fig. 1, and provides pre-GLM quality-control metrics used throughout the functional mapping workflow.

b. GLM-Based Functional Mapping
Script: GLM_NARS_$date.m
Function:
- Performs voxelwise GLM analysis using AFNI-derived regressors
- Supports ideal function regression analysis 
- Generates beta-coefficient and t-statistic maps
Notes for reviewers:
Voxelwise FDR is not emphasized due to sparse focal activation patterns; reproducibility is instead evaluated using cross-species consistency and resampling-based stability.
Manuscript linkage:
- Produces voxelwise ß-coefficient and t-statistic maps shown in Fig. 23 and Supplementary Fig. 45. Implements the GLM framework used for cross-species NARS-fMRI functional mapping.

c. Resampling-Based Reliability Analysis
Script: Randomized_$date.m
Function:
- Performs randomized selection of trial subsets
- Executes GLM refitting over ~1000 repetitions
- Quantifies spatial and statistical stability
Manuscript linkage:
- Reproduces the resampling-based stability analyses shown in Supplementary Fig. 78 and quantifies robustness of localized activation patterns across randomized trial subsets.

6. Simultaneous NARS-fMRI and Glutamate Fiber Photometry: /Glu_photometry/
Script: NARS_Glu_$date.m
Function:
- Aligns trial-resolved NARS-fMRI responses with glutamate photometry signals
- Sorts trials based on Glu amplitude
- Generates scatter plots of peak sodium vs. Glu responses
Manuscript linkage:
- Generates the cross-modal peak-response analysis shown in Fig. 4 and supports interpretation of NARS-fMRI signals relative to glutamatergic activity dynamics.

7. 23Na T1 and T2* mapping: /T1T2_mapping/ »
Script: T2_all.m  and T1_all.m
Function:
- Performs T2* biexponential fitting using fixed-ratio and free-ratio models
- Fits multi-TR acquisitions to estimate T1 relaxation maps
Outputs:
- Produces T1 and biexponential T2* fitting results shown in Supplementary Fig. 1 and provides relaxation parameters used to interpret NARS-fMRI signal mechanisms.

Data Availability and Reproducibility Notes
- Example datasets are included to demonstrate full pipeline functionality.
- AFNI preprocessing scripts for BOLD can be assessed through XXXX Gib link by Xiaochen
- Randomized resampling is implemented to evaluate robustness against pseudo-positive activation patterns.
- Scripts were used without blinding, as analyses were fully automated.

Contact
Xin Yu, Xiaochen Liu
Email: xyu9@mgh.harvard.edu
Translational Neuroimaging & Neural Control Laboratory
Massachusetts General Hospital / Harvard Medical School
