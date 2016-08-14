import numpy as np, sys, os
#========================================================================================#
def printHelp():
	print '''Argument flags:	
    -R <path_to_sources_sinks_vector>    (required; .npy .csv)
    -L <path_to_pseudoinverse_laplacian> (default: "Linv.npy" in the same directory as R; .npy)'''
	
#========================================================================================#
def main(argv):
	try:
		opts,args = getopt.getopt(argv, 'L:R:')
	except:
		print 'Inputs formatted incorrectly'
		printHelp()

	#get the arguments and turn them into variables
	path_to_Linv = None
	path_to_R = None
	
	for o,a in opts:
		if o == '-L': path_to_Linv = a
		if o == '-R': path_to_R = a

	#====================================================================================#

	if path_to_R == None: print 'You must input a source/sink vector using the -R flag'; sys.exit(2)
	if not os.path.exists(path_to_R):  print 'The file '+path_to_R+' does not exist'; sys.exit(2)
	if   '.csv' in path_to_R: R = np.loadtxt(path_to_R, delimiter=',')
	elif '.npy' in path_to_R: R = np.load(path_to_R)
	else: print 'The source/sink vector R must end in ".npy" or ".csv"'; sys.exit(2)
	
	if path_to_Linv == None:
		tmp_path_to_Linv = '/'.join(path_to_R.split('/')[:-1] + ['Linv.npy'])
		if os.path.exists(tmp_path_to_Linv): Linv = np.load(tmp_path_to_Linv)
		else: print 'You must create the file '+tmp_path_to_Linv+'. Use "compute_Linv.py"'; sys.exit(2)
	else:
		if not os.path.exists(path_to_Linv): print 'The file '+path_to_Linv+' does not exist'; sys.exit(2)
		else: Linv = np.load(path_to_Linv)
	
	# Compute the potential
	print 'Computing the potential V'
	V = np.dot(Linv, R)
	
	outpath = '/'.join(path_to_R.split('/')[:-1] + ['V.npy'])
	np.save(outpath, V)