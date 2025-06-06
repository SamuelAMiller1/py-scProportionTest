import matplotlib.pyplot as plt
import seaborn as sns

def cluster_map(results, 
                index='group2', 
                columns='cell_type', 
                data_values='observed_diff', 
                pvalue_values='adj_p_value',
                cluster_on='both', 
                asterisk_size=10,
                asterisk_color='black',
                alpha_threshold=0.05, 
                cluster_method='average', 
                cluster_metric='euclidean'):
    """
    Create a clustered heatmap using Seaborn's clustermap with observed differences and asterisks for significance.
    
    Parameters:
    - results: DataFrame containing the results of the permutation test.
    - index (str, optional): Index column for the pivot operation. Default is 'group2'.
    - columns (str, optional): Columns for the pivot operation. Default is 'cell_type'.
    - data_values (str, optional): Values column for the data_matrix pivot operation. Default is 'observed_diff'.
    - pvalue_values (str, optional): Values column for the pvalue_matrix pivot operation. Default is 'adj_p_value'.
    - cluster_on: Axis on which to perform clustering. Options are 'both', 'y', or 'x'.
    - asterisk_size: Size of the asterisk font.
    - asterisk_color: Color of the asterisk.
    - alpha_threshold: Threshold for significance for the adjusted p-value.
    - cluster_method (str, optional): Linkage algorithm to use for clustering. 
        Options['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'].
        Default is 'average'.
        
    - cluster_metric (str, optional): Distance metric to use for clustering. 
        Options['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'dice', 'euclidean', 
        'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 
        'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule'].
        Default is 'euclidean'.
    """
    
    if cluster_on not in ['both', 'y', 'x']:
        raise ValueError("'cluster_on' must be one of 'both', 'y', or 'x'")
    
    # Convert results into a matrix format for observed differences
    data_matrix = results.pivot(index=index, columns=columns, values=data_values)

    # Convert results into a matrix format for adjusted p-values
    pvalue_matrix = results.pivot(index=index, columns=columns, values=pvalue_values)

    # Define the custom colormap
    cmap = sns.diverging_palette(250, 15, s=75, l=40, n=9, center="light", as_cmap=True)

    # Determine clustering preferences based on 'cluster_on' argument
    row_cluster = True if cluster_on in ['both', 'y'] else False
    col_cluster = True if cluster_on in ['both', 'x'] else False

    # Generate the clustermap
    g = sns.clustermap(data_matrix, 
                       cmap=cmap, 
                       center=0,
                       method=cluster_method, # specify the clustering method
                       metric=cluster_metric, # specify the clustering metric
                       linewidths=0.75, 
                       figsize=(10, 8),
                       tree_kws=dict(linewidths=1.5),
                       annot=pvalue_matrix.applymap(lambda x: '*' if x < alpha_threshold else ''),
                       fmt='',
                       cbar_pos=(0.05, 0.8, 0.05, 0.18),
                       row_cluster=row_cluster,
                       col_cluster=col_cluster,
                       annot_kws={"size": asterisk_size, "color": asterisk_color})
    
    # Rotate the y-axis labels to print horizontally
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    
    # Add title to the colorbar
    g.cax.set_title("Observed log2(FD)", pad=15)

    #return object
    return g

    plt.show()
