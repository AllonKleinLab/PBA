import numpy as np, os, matplotlib.pyplot as plt

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

def color_by_fate(B):
	if B.shape[1] == 7: benchmarks = [1.5/7, 4./7, 5./7, 1., .5/7, 2.5/7, 6./7]
	if B.shape[1] == 2: benchmarks = [.3,.8]
	colors = [np.array(plt.cm.jet(benchmarks[i])[:3]) for i in range(B.shape[1])]
	
	c = np.zeros((B.shape[0],3))
	for i in range(B.shape[1]):
		ctile = 1 - np.tile(colors[i][None,:],(B.shape[0],1))
		ztile = np.tile(B[:,i][:,None],(1,3))
		c += ctile * ztile
	c = np.minimum(c,1)
	c = 1 - c
	c = c * .9	
	return c
	

#========================================================================================#

# first try to show results for the bifurcation example
if not os.path.exists('example_datasets/bifurcation/B.npy'):
	print '\nThere are no testing results for "example_datasets/bifurcation". To generate them, run the following:'
	print '\n    python PBA_pipeline.py -X example_datasets/bifurcation/X.npy -R example_datasets/bifurcation/R.npy -S example_datasets/bifurcation/S.npy\n'
else:
	coords = np.load('aux_files/bifurcation_coords.npy')
	xx = coords[:,0]; yy = coords[:,1]
	V = np.load('example_datasets/bifurcation/V.npy')
	B = np.load('example_datasets/bifurcation/B.npy')
	c = color_by_fate(B)
	fig,axs = plt.subplots(1,2)
	axs[0].scatter(xx,yy,c=V,edgecolor='',s=60)
	axs[0].set_title('Potnetial (V)')
	axs[1].scatter(xx,yy,c=c,edgecolor='',s=60)
	axs[1].set_title('Fate probabilities (B.npy)')
	for i in range(2):
		axs[i].set_xticks([])
		axs[i].set_yticks([])
	plt.suptitle('Testing results for example_datasets/bifurcation')
	plt.show()


# first try to show results for the hematopoiesis example
if not os.path.exists('example_datasets/hematopoiesis/B.npy'):
	print '\nThere are no testing results for "example_datasets/hematopoiesis". To generate them, run the following:'
	print '\n    python PBA_pipeline.py -X example_datasets/hematopoiesis/X.npy -R example_datasets/hematopoiesis/R.npy -S example_datasets/hematopoiesis/S.npy\n'
else:
	coords = np.load('aux_files/hematopoiesis_coords.npy')
	xx = coords[:,0]; yy = coords[:,1]
	V = np.load('example_datasets/hematopoiesis/V.npy')
	V = np.argsort(np.argsort(V))
	B = np.load('example_datasets/hematopoiesis/B.npy')
	c = color_by_fate(B)
	fig,axs = plt.subplots(1,2)
	axs[0].scatter(xx,yy,c=V,edgecolor='',s=20)
	axs[0].set_title('Potnetial (V.npy)')
	axs[1].scatter(xx,yy,c=c,edgecolor='k',s=20, linewidth=.1)
	axs[1].set_title('Fate probabilities (B.npy)')
	for i in range(2):
		axs[i].set_xticks([])
		axs[i].set_yticks([])
	plt.suptitle('Testing results for example_datasets/bifurcation')
	plt.show()
