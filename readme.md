# data-analysis

A library for data analysis that includes functions that return: correlation matrix 
heatmap, hierarchical clustering, minimum spanning tree. There is also implemented a 
function that use a bootstrap method as a filter for correlation matrices.


## Data

The function cmat takes as input a data file (.dat or .csv), with its format options,
filled with time series such as the returns of different stocks.
cmat returns the correlation matrix of this time serie as a nxn numpy array or a 
pandas dataframe, that can be used as input for the other functions. 

## Setup

The code is written in Python 3.

### Install pandas
https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html

### Install Igraph
https://igraph.org/python/

## Documentation

The reference for Correlation, hierarchies, and networks in financial markets is:
https://doi.org/10.1016/j.jebo.2010.01.004

The reference for the Bootstrap filtering method is:
https://doi.org/10.1016/j.physa.2018.08.020
