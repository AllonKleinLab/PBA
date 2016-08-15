import numpy as np, sys, os, getopt
#========================================================================================#
def printHelp():
	print '''Argument flags:	
	-X <path_to_expression_matrix>            (required; .npy or .csv)
	-E <minimum_mean_expression>              (default = -10000; used to filter genes)
	-V <minimum_CV>                           (default = 0; used to filter genes)
	-N <Normalize>                            (default = False; used to normalize expression data for knn graph)
	-p <PCA dimension>                        (default = 50; used to compute distance matrix)
	-k <number of nearest neighbors>          (default = 10; used to compute edge list)
    -e <path_to_edge_list>                    (required if no expression matrix is supplied)
    -R <path_to_sources_sinks_vector>         (required; .npy or .csv)
    -S <path_to_lineage_specific_sink_matrix> (optional, needed to compute fate probabilities; .csv or .npy)
    -D <diffusion_constant>                   (default = 1.0; controls the level of stochasticity in the model)\n'''


#========================================================================================#
def main(argv):
	try:
		opts,args = getopt.getopt(argv, 'X:E:V:p:k:e:R:S:D:N:')
	except:
		print '\nInputs formatted incorrectly'
		printHelp(); sys.exit(2)

	#get the arguments and turn them into variables
	path_to_expression_matrix = None
	minimum_mean_expression = -10000
	minimum_CV = 0
	normalize = False
	k = 10
	p = 60	
	path_to_edge_list = None
	path_to_R = None
	path_to_S = None
	D = 1.0
	
	for o,a in opts:
		if o == '-X': path_to_expression_matrix = a
		if o == '-E': minimum_mean_expression = float(a)
		if o == '-V': minimum_CV = float(a)
		if o == '-N': normalize = (a == 'True')
		if o == '-p': p = int(a)
		if o == '-k': k = int(a)
		if o == '-e': path_to_edge_list = a
		if o == '-R': path_to_R = a
		if o == '-S': path_to_S = a
		if o == '-D': D = float(a)

	#====================================================================================#
	
	for path in [path_to_expression_matrix, path_to_edge_list, path_to_R, path_to_S]:
		if path != None and not os.path.exists(path): print 'Error: The file '+path+' does not exist'; sys.exit(2)
		
	if path_to_expression_matrix == None and path_to_edge_list == None: 
		print 'Error: You must input either an expression matrix (-X) or knn edge list (-e)'; sys.exit(2)
	elif path_to_edge_list == None:
		print '\n## Running compute_knn_graph.py'
		os.system('python compute_knn_graph.py -E '+repr(minimum_mean_expression)+' -V '+repr(minimum_CV)+' -k '+repr(k)+' -p '+repr(p)+' -X '+ path_to_expression_matrix + ' -N '+repr(normalize))
		path_to_edge_list = '/'.join(path_to_expression_matrix.split('/')[:-1] + ['edge_list.csv'])
	
	print '\n## Running compute_Linv.py'
	os.system('python compute_Linv.py -e '+path_to_edge_list)
	
	print '\n## Running compute_potential.py'
	os.system('python compute_potential.py -R '+path_to_R)
	
	print '\n## Running compute_fate_probabilities.py'
	os.system('python compute_fate_probabilities.py  -S '+path_to_S+' -e '+path_to_edge_list+' -D '+repr(D))
	
if __name__ == '__main__':
	main(sys.argv[1:])
	
	
	
	