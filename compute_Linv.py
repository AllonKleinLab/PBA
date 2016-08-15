import numpy as np, sys, os, getopt
#========================================================================================#
def printHelp():
	print '''Argument flags:	
	-e <path_to_edge_list>   (required)\n'''

def row_sum_normalize(A):
	s = np.sum(A,axis=1)
	X,Y = np.meshgrid(s,s)
	A = A / Y
	return A

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
	return A
	
#========================================================================================#
def main(argv):
	try:
		opts,args = getopt.getopt(argv, 'e:')
	except:
		print '\nInputs formatted incorrectly'
		printHelp(); sys.exit(2)

	#get the arguments and turn them into variables
	path_to_edge_list = None
	
	for o,a in opts:
		if o == '-e': path_to_edge_list = a

	#====================================================================================#
	
	if path_to_edge_list == None: print 'Error: You must input an edge list using the -e flag'; sys.exit(2)
	if not os.path.exists(path_to_edge_list):  print 'Error: The file '+path_to_edge_list+' does not exist'; sys.exit(2)
	
	# Make Laplacian matrix
	print 'Making Laplancian matrix'
	edges = load_edge_list(path_to_edge_list)
	A = make_adjacency_matrix(edges)
	L = np.identity(A.shape[0]) - row_sum_normalize(A)

	# Invert graph Laplacian
	print 'Inverting the Laplacian'
	Linv = np.linalg.pinv(L)	
	
	outpath = '/'.join(path_to_edge_list.split('/')[:-1] + ['Linv.npy'])
	np.save(outpath, Linv)


if __name__ == '__main__':
	main(sys.argv[1:])

		