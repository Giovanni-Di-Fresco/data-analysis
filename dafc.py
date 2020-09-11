import os
import numpy as np
import scipy
import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.sparse.csgraph import minimum_spanning_tree
from collections import defaultdict,Counter
from scipy.spatial.distance import squareform


import igraph as ig



                                                                            

def cmat(fpath, filename, se=',', nparray=False):
	'''
	Returns the correlation matrix. 
	
	Parameters:
		fpath (str): the path of the folder in which the file is stored.
		filename (str): the name of the file.
		se (str): the separator of the data in the file.
		nparray (Boolean): set as True for numpy array, False for pandas dataframe.
	
	Returns:
		arr (numpy array): the correlation matrix as a nxn numpy array.
		corrmatrix (pandas dataframe): the correlation matrix as a pandas dataframe. 		

	'''

	os.chdir(fpath)
	data=pd.read_csv(filename,sep=se) 
	corrmatrix=data.corr()
	if(nparray):
		arr=corrmatrix.to_numpy()
		return(arr)

	return(corrmatrix)  


def heatmap(a, fgsx=15, fgsy=15, mticks=False, title='Correlation Matrix', 
	ftitle=15, colbar=True, colsize=15,
	colormap='Spectral'):
	'''
	Returns the heatmap of a nxn pandas dataframe.
	
	Parameters:
		a (pandas dataframe): the nxn pandas dataframe which is to be plotted.
		fgsx (int): plot x size.
		fgsy (int): plot y size.
 		mticks (Boolean): set as True in order to plot ticks in the graph.
 		title (str): title of the graph.
 		ftitle (int): fontsize of the title.
 		colbar (Boolean): set as True in order to plot the color legend.
 		colormap (ste): name of the colormap theme, 
 		check https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html
 		for all the possibilities.
 	Returns:
 		none.
	'''

	f = plt.figure(figsize=(fgsx, fgsy))
	plt.imshow(a, cmap=colormap)
	if(mticks):
		plt.xticks(range(a.shape[1]), a.columns, fontsize=3, rotation=45)
		plt.yticks(range(a.shape[1]), a.columns, fontsize=3)
	if(colbar):
		cb = plt.colorbar()
		cb.ax.tick_params(labelsize=colsize)
	plt.title('Correlation Matrix', fontsize=ftitle);
	
	plt.show()


def hierarchical_order(corr_array, inplace=False, met='complete',
	dendogram=False, order=True):
	'''
	Perform the hierarchical clustering of a nxn matrix, and 
	return the ordered correlation matrix.
	
	Parameters:
		corr_array (numpy array or pandas dataframe): correlation matrix.
		inplace (Boolean): make a copy of corr_array.
		met (str): linkage method, chek scipy.cluster.hierarchy.linkage for alternatives.
		dendogram (Boolean): set as True in order to hide the dendogram plot.
		order (Boolean): set as True in order to return an ordered correlation matrix.
	
	Return:
		corr_array (numpy array or pandas dataframe): returns the ordered correlation 
		matrix.
	
 
	'''
	
 
	d=np.sqrt(2*(1-corr_array))
#	y=np.triu(d).flatten().T[np.triu(d).flatten().T !=0]
	y=squareform(d)
	pairwise_distances = y
	linkage = sch.linkage(pairwise_distances, method=met)
	fig = plt.figure(figsize=(25, 10))
	dn = sch.dendrogram(linkage,no_plot=dendogram)
	if(order):
		idx_to_cluster_array=sch.leaves_list(linkage)
		idx = np.argsort(idx_to_cluster_array) 
		if not inplace:
			corr_array = corr_array.copy()
    
		if isinstance(corr_array, pd.DataFrame):
			return corr_array.iloc[idx, :].T.iloc[idx, :]
		return corr_array[idx, :][:, idx]


def mst(corr, name=False, pandas= False):
	'''
	Returns the minimum spanning tree. 
	
	Parameters:
		corr(numpy array): nxn correlation matrix.
		name (list): a list with the name of the vertex. 
 		pandas (str): True if the imput is a pandas datadrame.
 		
 	Returns:
 		g (igraph): the minimum spanning tree. 
	''' 
	if(pandas):
		corr=corr.to_numpy()
	
	d=np.sqrt(2*(1-corr))
	l=np.int(len(corr)*len(corr)-len(corr)*(len(corr)-1)/2)
	a=squareform(d)
	b,c=np.triu_indices(len(d),k=1)
	d2=np.array(list(zip(a,c,b)))
	d2= d2[np.argsort(d2[:, 0])]
	ls=np.zeros((len(corr), 2))
	ls[:, 0]=np.arange(0, len(corr))
	ls[:, 1]=np.arange(0, len(corr))
	g=ig.Graph()
	g.add_vertices(len(corr))
	i=0
	h=0
	ls=np.zeros((len(corr), 2))
	ls[:, 0]=np.arange(0, len(corr))
	ls[:, 1]=np.arange(0, len(corr))
	g=ig.Graph()
	g.add_vertices(len(corr))
	if(name):
		g.vs["label"] = name
	counter=0
	while(counter<len(corr)-1):
		a1=int(d2[i, 1])
		a2=int(d2[i, 2])
 
		if(ls[a2, 1] != ls[a1, 1]):
			h=h+1
			g.add_edges([(a1, a2)])
			ck=ls[a2, 1]
			ck1=ls[a1, 1]
			ls[ls[:, 1]==ck, 1]=ck1
			counter+=1
		i=i+1
	return g
 
 
def boo(data, bootstrap):
	'''
	Returns a weighted graph containing all the edges of n minimun spanning tree
	obtained Bootstrapping the correlation matrix.
	
	Parameters:
		data (numpy array): nxn correlation matrix.
		bootstrap (int): number of shuffling
		
	Returns:
		g (igraph): weighted graph.
	
	'''
	
	d=defaultdict(int)
	for i in range(bootstrap):
		w=np.random.choice(range(len(data)),len(data))
		corr1=data.loc[w].corr()
		h=mst(corr1)
		k=h.get_edgelist()
		for j in range(len(k)):
			d[k[j]]+=1
	g=ig.Graph()
	g.add_vertices(len(d))
	g.add_edges(d.keys())
	g.es['weight']=list(d.values())	
	return(g)


def boo_t(data, bootstrap, tolerance):
	'''
	Returns a filtered weighted graph containing the edges of n minimun spanning tree
	obtained Bootstrapping the corralation matrix.
	
	Parameters:
		data (numpy array): nxn correlation matrix.
		bootstrap (int): number of shuffling
		tolerance (int): weight tolerance of the graph.
		
	Returns:
		g (igraph): weighted graph.

	'''
	
	g=boo(data, bootstrap)
	k1=g.get_edgelist()
	edges=[k1[i] for i in range(len(k1)) if (g.es(i)['weight'][0]<tolerance)]
	g.delete_edges(edges)
	return(g)




 
 
 