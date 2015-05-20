# generate simulation data
matlab -nodesktop -nosplash -r simdata_gen\(\'300\'\,\'10\'\,\'100\'\,\'20\'\,\'1\'\)\;exit\;
matlab -nodesktop -nosplash -r simdata_gen\(\'300\'\,\'10\'\,\'100\'\,\'20\'\,\'2\'\)\;exit\;
matlab -nodesktop -nosplash -r simdata_gen\(\'300\'\,\'10\'\,\'100\'\,\'20\'\,\'3\'\)\;exit\;

# run simulation 0 - 3
matlab -nodesktop -nosplash -r simulation0\;exit\;
matlab -nodesktop -nosplash -r simulation1\;exit\;
matlab -nodesktop -nosplash -r simulation2\;exit\;
matlab -nodesktop -nosplash -r simulation3\;exit\;
