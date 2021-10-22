# NeurIPS-2021-Paper-9402


Matlab files are used to generate the theoretical 3D plots while the python notebook generates the images from the experiments with an implementation of the RF model using Tensorflow.


### MATLAB ###


* The main entry points are the "config*.m" files where * is a number. A config file generates theoretical data and writes them in the "output/config*" directory. 
* the "MSEplot.m" file uses the former output to generate the theoretical heat maps (or 3D landscapes)
* The former script requires the multigradient.m file from https://github.com/lrkrol/multigradient to define the color gradient for the heatmaps (not mandatory)
* "plotIntro.m" generates the images in the introduction


### PYTHON ###

* All the images were generated using a standard (free) instance of Google Colaboratory
* The 4 CSV files are theoretical data generated from the Matlab scripts for comparison with the experimental curves (from the scripts config0.m, htest_profile.m, trainingDoubleDescentEvidence.m)

