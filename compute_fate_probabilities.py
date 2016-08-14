import numpy as np, sys, os
#========================================================================================#
def printHelp():
	print '''Argument flags:	
    -S <path_to_lineage_specific_sink_matrix> (required; .csv or .npy)
    -V <path_to_potential_vector>             (default: "V.npy" in same directory as S; .npy or .csv)
    -e <path_to_edge_list>                    (default: "edge_list.csv" in same directory as S)'''

def load_edge_list(path):
	edge_text = open(path).read()
	edge_text = edge_text.replace('\r','\n').replace('\n\n','\n').strip('\n')
	out = []
	for l in edge_text.split('\n'):
		l = l.split(',')
		i = int(l[0])
		j - int(l[1])
		out += [(i,j)]
	return out

def make_adjacency_matrix(edges):
	N = np.max([i for i,j in edges]+[j for i,j in edges])+1
	A = np.zeros((N,N))
	for i,j in edges: 
		A[i,j] = 1.
		A[j,i] = 1.
	return A
	
def row_sum_normalize(A):
	d = np.sum(A,axis=1)
	A = A / (np.tile(d[:,None],(1,A.shape[1])) + .01)
	return A


#========================================================================================#
def main(argv):
	try:
		opts,args = getopt.getopt(argv, 'S:V:e:')
	except:
		print 'Inputs formatted incorrectly'
		printHelp()

	#get the arguments and turn them into variables
	path_to_S = None
	path_to_V = None
	path_to_edge_list = None
	
	for o,a in opts:
		if o == '-S': path_to_S = a
		if o == '-V': path_to_V = a
		if o == '-e': path_to_edge_list = a

	#====================================================================================#

	if path_to_S == None: print 'You must input a lineage-specific sink matrix using the -S flag'; sys.exit(2)
	if not os.path.exists(path_to_S):  print 'The file '+path_to_S+' does not exist'; sys.exit(2)
	if   '.csv' in path_to_S: S = np.loadtxt(path_to_S, delimiter=',')
	elif '.npy' in path_to_S: S = np.load(path_to_S)
	else: print 'The lineage-specific sink matrix S must end in ".npy" or ".csv"'; sys.exit(2)
	
	if path_to_V == None:
		tmp_path_to_V = '/'.join(path_to_S.split('/')[:-1] + ['V.npy'])
		if os.path.exists(tmp_path_to_V): V = np.load(tmp_path_to_V)
		else: print 'You must create the file '+tmp_path_to_V+'. Use "compute_potential.py"'; sys.exit(2)
	else:
		if not os.path.exists(path_to_V): print 'The file '+path_to_V+' does not exist'; sys.exit(2)
		else: V = np.load(path_to_V)

	if path_to_edge_list == None:
		tmp_path_to_edge_list = '/'.join(path_to_S.split('/')[:-1] + ['edge_list.csv'])
		if os.path.exists(tmp_path_to_edge_list): edges = load_edge_list(tmp_path_to_edge_list)
		else: print 'You must create the file '+tmp_path_to_edge_list+'. Use "compute_knn_graph.py"'; sys.exit(2)
	else:
		if not os.path.exists(path_to_edge_list): print 'The file '+path_to_edge_list+' does not exist'; sys.exit(2)
		else: edges = load_edge_list(path_to_edge_list)
		
	# Make Markov chain transition matrix
	print 'Making Markov chain transition matrix'	
	A = make_adjacency_matrix(edges)
	Vx,Vy = np.meshgrid(V,V)
	T = np.exp(A * np.minimum(Vy - Vx, 50))
	bigT = np.hstack((T,S))
	bigT = np.vstack((bigT,np.hstack((np.zeros((T.shape[1],S.shape[1])),np.identity(S.shape[1])))))
	bigT = row_sum_normalize(bigT)
	
	# compute fundamental matrix
	print 'Computing fundamental matrix'
	Q = bigT[:T.shape[0],:T.shape[0]]
	N = np.linalg.pinv(np.identity(Q.shape[0]) - Q)
	
	# compute fate probabilities
	print 'Computing fate probabilities'
	B = np.dot(N,S)
	
	outpath = '/'.join(path_to_S.split('/')[:-1] + ['B.npy'])
	np.save(outpath,B)
	
	
	
	
	
