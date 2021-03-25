# Analysis of Data-Driven, Neuron-Type Specific CA3 SNNs
This repository includes information as to how to analyze a full-scale spiking neural network (SNN) model of hippocampal subregion CA3. The following instructions assume that the user has MATLAB installed and that they have installed CARLsim according to instructions listed in this [README](https://github.com/UCI-CARL/CARLsim4/tree/feat/meansdSTPPost_hc).

## Choosing a network to analyze
There are three directories reflecting the type of CA3 SNNs that can be analyzed: [baseline_analysis](https://github.com/Hippocampome-Org/snn_analysis/tree/main/baseline_analysis), which provides code to analyze the full-scale baseline SNN that maintains neuron and connection-type specificity; [class_analysis](https://github.com/Hippocampome-Org/snn_analysis/tree/main/class_analysis) provides code to analyze the full-scale class SNN that maintains neuron-type specificity while removing connection-type specificity; [archetype_analysis](https://github.com/Hippocampome-Org/snn_analysis/tree/main/archetype_analysis) provides code to analyze the full-scale archetype SNN that maintains connection-type specificity while removing neuron-type specificity.

## Network-Independent Analysis Instructions
Before any analysis can be performed, a directory that will contain the simulation output must be created. The main analysis function assumes a directory structure where the simulation output folders are in their own separate directory, and this directory should reside where the analysis code for the network type of choice is being analyzed (e.g. in the baseline_analysis directory).

1. Move to the directory where the CARLsim MATLAB Analysis Toolbox (MAT) has been installed and run the initOAT command which provides functions that are necessary for performing analysis.

  ```
 cd MAT_directory
 initOAT
  ```
 
 2. Move to the directory of the SNN type you wish to analyze, and open up the file evaluateSim.m. Update the fileLoc variable name with a directory of your choice, which is where the plots from the analysis code will be output. Then run the main function evaluteSim, which will perform all analysis necessary to measure average frequencies for each neuron type of the network, power spectral density analysis, and spike-to-phase (SPC) relationships.

  ```
 cd analysis_directory_of_choice
 fileLoc = 'yourChoiceDrive:\your_directory_choice'
 evaluateSim
  ```

3. Resultant plots then can be viewed in the output directory chosen. Additionally, the matrices that contain simulation statistics, average frequencies, power spectral density peaks and powers, and SPC relationships can be viewed.
