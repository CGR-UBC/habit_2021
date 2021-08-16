# habit_2021

### Description

This repository contains the scripts to reproduce the analyses described in ‘Behavioral analysis of habit formation in modern slot machine gambling’. Code for the six models presented in the manuscript is contained in habit_analysis_script.Rmd, and supporting functions are in habit_functions.R.

### Instructions

Within the R Markdown setup chunk, set the working directory to the folder containing these two files and the data file. If all required R packages are installed, the analysis script can be run in full. Alternatively, code chunks for the six models can be run independently, provided the setup and data preparation chunks have been run. Chunks that produce figures and results tables will require the corresponding models to be run in advance. The majority of results, tables, and figures will be saved in a folder called 'results' within the working directory. Some results will be printed inline.