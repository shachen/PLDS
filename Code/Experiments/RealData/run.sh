# KKI Data Analysis
# select latent dimension
R CMD BATCH kkiDimensionSelection.R
# generate intial value
R CMD BATCH kki_group_pca_init.R
rm *.Rout
# analyze data
matlab -nodesktop -nosplash -r kkiMotorAnalyze\(\)\;exit\;
# analyze results
R CMD BATCH kkiResultAnalyze.R
rm *.Rout

# HCP Data Analysis
# generate inital value
R CMD BATCH hcpInitialValue.R
rm *.Rout
# analyze HCP data
matlab -nodesktop -nosplash -r hcpMotorAnalyze\(\)\;exit\;
