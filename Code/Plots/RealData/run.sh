# KKI 3D Plots
R CMD BATCH kki_3dplots.R

# KKI Connectivity Graphs
R CMD BATCH kkiResultAnalyze.R

rm *.Rout 
rm *.pdf

# HPC Prediction Plots
R CMD BATCH hcpMakePredictions.R

rm *.Rout
