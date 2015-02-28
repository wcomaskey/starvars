
StarVars Repository: Classifying variability throughout the Universe
========


Files:
----------

code/:

classify_linear.R -  Run classifier (trained on Ogle+Hipparcos+ASAS) to classify LINEAR
color_convert_ugriz_BVR.R -  Function to compute colors from LINEAR
utils_classify.R - Set of routines used by classify_linear.R

data/:

training_set_features.dat - Training set features for Ogle+Hipparcos+ASAS labeled set
training_set_id_class.dat - Traiing set classes
200k_final_combo.csv - LINEAR features (for ~200k data set)
masterMain.dat.txt - metadata for LINEAR set, including ra/dec, colors, etc.
LINEARattributesFinalApr2013.dat - set of manually labeled LINEAR objects
QSO-Linear.dat - set of LINEAR QSOs (from SDSS)

plots/ and tables/:

These hold all of the output figures and tables from running classify_linear.R