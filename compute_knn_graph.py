import numpy as np, sys, os
#========================================================================================#
def printHelp():
	print '''Argument flags:	
	-X <path_to_expression_matrix>   (required; .npy or .csv)
	-E <minimum_mean_expression>     (default = 0.05; used to filter genes)
	-V <minimum_CV>                  (default = 2; used to filter genes)
	-p <PCA dimension>               (default = 50; used to compute distance matrix)
	-k <number of nearest neighbors> (default = 10; used to compute edge list)'''

def get_distance_matrix(M):
	D = np.zeros((M.shape[0],M.shape[0]))
	for i in range(M.shape[0]):
		Mtiled = np.tile(M[i,:][None,:],(M.shape[0],1))
		D[i,:] = np.sqrt(np.sum((Mtiled - M)**2, axis=1))
	return D

def row_normalize(E):
	total_counts = np.sum(E,axis=1)
	return E / np.tile(np.sum(E,axis=1)[:,None],(1,E.shape[1])) * np.mean(total_counts)

def Zscore(E):
	means = np.tile(np.mean(E,axis=0)[None,:],(E.shape[0],1))
	stds = np.tile(np.std(E,axis=0)[None,:],(E.shape[0],1))
	return (E - means) / (stds + .0001)

def filter_genes(E, Ecutoff, Vcutoff):
	mean_filter = np.mean(E,axis=0)> Ecutoff
	var_filter = np.var(E,axis=0) / (np.mean(E,axis=0)+.0001) > Vcutoff
	gene_filter = np.nonzero(np.all([mean_filter,var_filter],axis=0))[0]
	return E[:,gene_filter], gene_filter
	
def get_knn_edges(dmat, k):
	edges = set([])
	for i in range(dmat.shape[0]):
		for j in np.nonzero(dmat[i,:] <= sorted(dmat[i,:])[k])[0]:
			if i != j:
				ii,jj = tuple(sorted([i,j]))
				edges.add((ii,jj))
	return list(edges)

def get_PCA(A,numpc):
	try:
		from sklearn.decomposition import PCA
		pca = PCA(n_components=numpc)
		return pca.fit_transform(A)
	except:
		print 'WARNING: PCA is going to be very slow. Install the "scikit-learn" module for a big speedup'
		A -= A.mean(axis=0)
		R = np.cov(A, rowvar=False)
		evals, evecs = np.linalg.eigh(R)
		idx = np.argsort(evals)[::-1]
		evecs = evecs[:,idx]
		evals = evals[idx]
		evecs = evecs[:, :numpc]
		return np.dot(evecs.T, A.T).T

def save_edge_list(edges,path):
	out = []
	for i,j in edges:
		out += [repr(i)+','+repr(j)]
	open(path,'w').write('\n'.join(out))
	
#========================================================================================#
def main(argv):
	try:
		opts,args = getopt.getopt(argv, 'X:E:V:p:k:')
	except:
		print 'Inputs formatted incorrectly'
		printHelp()

	#get the arguments and turn them into variables
	path_to_expression_matrix = None
	minimum_mean_expression = 0.05
	minimum_CV = 2
	k = 10
	p = 50	
	
	for o,a in opts:
		if o == '-X': path_to_expression_matrix = a
		if o == '-E': minimum_mean_expression = float(a)
		if o == '-V': minimum_CV = float(a)
		if o == '-p': p = int(a)
		if o == '-k': k = int(a)

	#====================================================================================#
	
	# load input data
	print 'Loading expression data
	if path_to_expression_matrix == None: print 'You must input an expression matrix using the -X flag'; sys.exit(2)
	if not os.path.exists(path_to_expression_matrix): print 'The file '+path_to_expression_matrix+' does not exist'; sys.exit(2)
	if   '.csv' in path_to_expression_matrix: E = np.loadtxt(path_to_expression_matrix, delimiter=',')
	elif '.npy' in path_to_expression_matrix: E = np.load(path_to_expression_matrix)
	else: print 'The expression file must end in ".npy" or ".csv"'; sys.exit(2)
	
	# Filtering and normalization
	print 'Filtering and normalization'
	EE = row_normalize(E)
	EE,gene_filter = filter_genes(EE, minimum_mean_expression, minimum_CV)
	if np.sum(gene_filter) == 0: print 'No genes survived filtering, consider lowering the minimum mean expression (-E) or minimum CV (-V)'
	EZ = Zscore(EE)
	
	# PCA and distance matrix
	print 'Performing PCA'
	if p < EZ.shape[1]: Epca = get_PCA(EZ,p)
	else: Epca = EZ
	
	print 'Computing distance matrix'
	D = get_distance_matrix(Epca)
	
	# get knn edges
	print 'Getting knn edges'
	edges = get_knn_edges(D,k)
	
	outpath = '/'.join(path_to_expression_matrix.split('/')[:-1]+['edge_list.csv'])
	save_edge_list(edges, outpath)
	

if __name__ == '__main__':
	main(sys.argv[1:])
	
