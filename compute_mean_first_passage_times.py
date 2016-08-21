import numpy as np, sys, os, getopt
#========================================================================================#
def printHelp():
	print '''Argument flags:	
    -R <path_to_sources_sinks_vector>    (required; .npy .csv)
    -V <path_to_potential_vector>        (default; "V.npy" in same directory as R)
    -e <path_to_edge_list>               (default: "edge_list.csv" in same directory as R)
    -D <diffusion_constant>              (default = 1.0; controls the level of stochasticity in the model)\n\n'''

def load_edge_list(path):
	edge_text = open(path).read()
	edge_text = edge_text.replace('\r','\n').replace('\n\n','\n').strip('\n')
	out = []
	for l in edge_text.split('\n'):
		l = l.split(',')
		i = int(l[0])
		j = int(l[1])
		out += [(i,j)]
	return out

def make_adjacency_matrix(edges):
	N = np.max([i for i,j in edges]+[j for i,j in edges])+1
	A = np.zeros((N,N))
	for i,j in edges: 
		A[i,j] = 1.
		A[j,i] = 1.
		A[i,i] = 1.
		A[j,j] = 1.
	return A
	
def row_sum_normalize(A):
	d = np.sum(A,axis=1)
	A = A / (np.tile(d[:,None],(1,A.shape[1])))
	return A


#========================================================================================#
def main(argv):
	try:
		opts,args = getopt.getopt(argv, 'R:V:e:D:')
	except:
		print '\nInputs formatted incorrectly'
		printHelp(); sys.exit(2)

	#get the arguments and turn them into variables
	path_to_R = None
	path_to_V = None
	path_to_edge_list = None
	D = 1.0
	
	for o,a in opts:
		if o == '-R': path_to_R = a
		if o == '-V': path_to_V = a
		if o == '-e': path_to_edge_list = a
		if o == '-D': D = float(a)

	#====================================================================================#

	if path_to_R == None: print 'Error: You must input a source/sink vector using the -R flag'; sys.exit(2)
	if not os.path.exists(path_to_R):  print 'Error: The file '+path_to_R+' does not exist'; sys.exit(2)
	if   '.csv' in path_to_R: R = np.loadtxt(path_to_R, delimiter=',')
	elif '.npy' in path_to_R: R = np.load(path_to_R)
	else: print 'Error: The source/sink vector R must end in ".npy" or ".csv"'; sys.exit(2)
	
	if path_to_V == None:
		tmp_path_to_V = '/'.join(path_to_R.split('/')[:-1] + ['V.npy'])
		if os.path.exists(tmp_path_to_V): V = np.load(tmp_path_to_V)
		else: print 'Error: You must create the file '+tmp_path_to_V+'. Use "compute_potential.py"'; sys.exit(2)
	else:
		if not os.path.exists(path_to_V): print 'The file '+path_to_V+' does not exist'; sys.exit(2)
		else: V = np.load(path_to_V)

	if path_to_edge_list == None:
		tmp_path_to_edge_list = '/'.join(path_to_R.split('/')[:-1] + ['edge_list.csv'])
		if os.path.exists(tmp_path_to_edge_list): edges = load_edge_list(tmp_path_to_edge_list)
		else: print 'Error: You must create the file '+tmp_path_to_edge_list+'. Use "compute_knn_graph.py"'; sys.exit(2)
	else:
		if not os.path.exists(path_to_edge_list): print 'The file '+path_to_edge_list+' does not exist'; sys.exit(2)
		else: edges = load_edge_list(path_to_edge_list)
		
	# Make Markov chain transition matrix
	T = np.zeros((len(R),len(R)))
	A = make_adjacency_matrix(edges)
	V = V / D
	Vx,Vy = np.meshgrid(V,V)
	P = A * np.exp(np.minimum(Vy - Vx, 400))

	for i in range(len(R)):
		print 'Calculating MFPTs from all nodes to node',i	
		bigP = np.zeros((len(R)+1,len(R)+1))
		bigP[:len(R),:len(R)] = P
		bigP[:len(R),-1] = -np.minimum(R,0)
		bigP[-1,-1] = 1
		new_index = range(A.shape[0]+1)
		new_index.remove(i)
		new_index.append(i)
		new_index = np.array(new_index)
		bigP = bigP[new_index,:][:,new_index]
		bigP[-1,:] = 0
		bigP[-1,-1] = 1
		bigP = row_sum_normalize(bigP)
		Q  = bigP[:len(R)-1,:len(R)-1]
		RR = bigP[:len(R)-1,len(R)-1:]
		N = np.linalg.inv(np.identity(Q.shape[0])-Q)
		B = np.dot(N,RR)
		d = np.diag(B[:,-1])
		dinv = np.diag(1./B[:,-1])
		T[np.arange(len(R))!=i,i] = np.dot(np.dot(np.dot(dinv,N),d), np.ones(d.shape[0]))
		
	outpath = '/'.join(path_to_R.split('/')[:-1] + ['T.npy'])
	np.save(outpath,T)

if __name__ == '__main__':
	main(sys.argv[1:])
		
