import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import TwoSlopeNorm
import seaborn as sns

def point_range_plot(results, alpha=0.05, fold_difference=1, figsize=(6,6), dot_size=6, significance_color='darkred', 
                     ci_line_width=1.2, plot_title='Permutation Test Results: Observed Log2 Fold Differences', 
                     x_label='Log2 Fold Difference', y_label='Cell Type'):
    """
    Generate a point-range plot of the permutation test results, including the observed log2 fold differences
    and bootstrapped confidence intervals for each cell type. Points and confidence intervals are colored 
    based on adjusted p-value significance and an absolute log2 fold difference threshold.

    Parameters
    ----------
    results : DataFrame
        A pandas DataFrame with the permutation p-values, observed proportional differences,
        and bootstrapped confidence intervals for each cell type.
    alpha : float, optional
        The significance level for the p-value threshold. Default is 0.05.
    fold_difference : float, optional
        The absolute log2 fold difference threshold for coloring. Default is 1.
    figsize : tuple, optional
        The size of the figure for the plot. Default is (6, 6).
    dot_size : float, optional
        Size of the dot in the plot. Default is 6.
    significance_color : str, optional
        Color for significant dots and confidence intervals. Default is 'darkred'.
        For a full list of named colors supported in Matplotlib, see:
        https://matplotlib.org/stable/gallery/color/named_colors.html
    ci_line_width : float, optional
        Width of the confidence interval lines. Default is 1.2.
    plot_title : str, optional
        Title for the plot. Default is 'Permutation Test Results: Observed Log2 Fold Differences'.
    x_label : str, optional
        Label for the x-axis. Default is 'Log2 Fold Difference'.
    y_label : str, optional
        Label for the y-axis. Default is 'Cell Type'.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object of the plot.
    """

    # Error handling for argument types and values
    if not isinstance(results, pd.DataFrame):
        raise ValueError("`results` should be a pandas DataFrame.")

    if not all(col in results.columns for col in ['adj_p_value', 'observed_diff', 'lower_ci', 'upper_ci', 'cell_type']):
        raise ValueError("`results` DataFrame must contain 'adj_p_value', 'observed_diff', 'lower_ci', 'upper_ci', and 'cell_type' columns.")

    if not (0 <= alpha <= 1):
        raise ValueError("`alpha` should be between 0 and 1.")

    if not isinstance(figsize, tuple) or len(figsize) != 2:
        raise ValueError("`figsize` should be a tuple of two values representing the width and height of the figure.")
    
    if not isinstance(dot_size, (int, float)):
        raise ValueError("`dot_size` should be an integer or float.")
        
    if not isinstance(significance_color, str):
        raise ValueError("`significance_color` should be a string representing a color.")
    
    if not isinstance(ci_line_width, (int, float)):
        raise ValueError("`ci_line_width` should be an integer or float.")
    
    if not isinstance(plot_title, str):
        raise ValueError("`plot_title` should be a string.")
    
    if not isinstance(x_label, str):
        raise ValueError("`x_label` should be a string.")
    
    if not isinstance(y_label, str):
        raise ValueError("`y_label` should be a string.")
    
    # Begin plotting
    try:
        fig, ax = plt.subplots(figsize=figsize)
        
        # Add gridlines
        ax.grid(color='grey', linestyle='-', linewidth=0.5, which='both', axis='x')

        # Add dotted black lines at pos fold_difference and neg fold_difference
        ax.axvline(x=fold_difference, color='black', linestyle='--', linewidth=ci_line_width)
        ax.axvline(x=-fold_difference, color='black', linestyle='--', linewidth=ci_line_width)

        # Plot the observed log2 fold differences with confidence intervals
        for index, row in results.iterrows():
            significant = (row['adj_p_value'] < alpha) and (abs(row['observed_diff']) > fold_difference)
            color = significance_color if significant else 'grey'
            label = f'p-adj < {alpha}, |log2FD| > {fold_difference}' if significant else 'Not Significant'
            ax.plot([row['lower_ci'], row['upper_ci']], [row['cell_type'], row['cell_type']], color=color, linewidth=ci_line_width)  # Range line
            ax.plot(row['observed_diff'], row['cell_type'], 'o', color=color, markersize=dot_size, label=label)  # Point

        # Label the plot
        ax.set_title(plot_title)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

        # Add legend without duplicate labels
        handles, labels = ax.get_legend_handles_labels()
        unique_labels = dict(zip(labels, handles))
        ax.legend(unique_labels.values(), unique_labels.keys(), loc='center left', bbox_to_anchor=(1, 0.5))

        return fig
        
    except Exception as e:
        raise RuntimeError(f"An error occurred during plotting: {e}")



def multiple_permutation_dot_plot(results, 
                                  alpha_threshold=0.05, 
                                  observed_diff_threshold=0.58, 
                                  size_scale=500, 
                                  min_size=50,
                                  edgecolor='orange',
                                  linewidth=3,
                                  plot_title='Observed Differences and Adjusted P-value Sizes for Each Cell Type',
                                  xlabel='Cell Types',
                                  ylabel='Comparison Groups',
                                  rotation=45):
    """
    Plot the results of a permutation test with observed differences and adjusted p-values.

    Parameters:
    - results (DataFrame): DataFrame containing the results of the permutation test. 
                           It must contain the following columns:
                           - 'cell_type': Cell types being compared.
                           - 'group2': Comparison groups.
                           - 'observed_diff': Observed differences between the groups.
                           - 'adj_p_value': Adjusted p-values for the comparisons.
    - alpha_threshold (float, optional): Threshold for significance for the adjusted p-value. 
                                         Points with adjusted p-values below this threshold are considered significant. 
                                         Default is 0.05.
    - observed_diff_threshold (float, optional): Threshold for the absolute value of the observed log2 fold difference. 
                                                 Points with an absolute observed difference greater than this threshold are considered significant. 
                                                 Default is 0.58.
    - size_scale (float, optional): Scaling factor for the size of the points based on the negative logarithm of the adjusted p-value. 
                                    Higher values will result in larger points. Default is 500.
    - min_size (float, optional): Minimum size for the points, to ensure very small p-values are visible. Default is 50.
    - plot_title (str, optional): Title of the plot. Default is 'Observed Differences and Adjusted P-value Sizes for Each Cell Type'.
    - xlabel (str, optional): Label for the x-axis. Default is 'Cell Types'.
    - ylabel (str, optional): Label for the y-axis. Default is 'Comparison Groups'.
    - rotation (int, optional): Rotation for x-axis labels. Default is 45.

    Returns:
    - None: The function creates a plot and displays it using plt.show().

    The plot uses a diverging color palette to represent the observed differences and scales the size of the points based on the negative logarithm of the adjusted p-values. 
    Points that are considered significant based on the alpha_threshold and observed_diff_threshold are outlined in orange.

    Example:
    >>> multiple_permutation_dot_plot(results)
    """
    
    # Ensure 'results' is a DataFrame and contains required columns
    if not isinstance(results, pd.DataFrame):
        raise ValueError("The 'results' parameter must be a pandas DataFrame.")
    
    required_columns = ['cell_type', 'group2', 'observed_diff', 'adj_p_value']
    missing_columns = [col for col in required_columns if col not in results.columns]
    if missing_columns:
        raise ValueError(f"The 'results' DataFrame is missing the following required columns: {', '.join(missing_columns)}")
    
    try:
        # Add a new column for -log10(adj_p_value) and apply the size scaling
        results['size'] = np.maximum(-np.log10(results['adj_p_value']) * size_scale, min_size)

        # Determine the number of unique cell types and comparison groups
        num_cell_types = results['cell_type'].nunique()
        num_groups = results['group2'].nunique()

        # Set the figure size
        fig_width = max(8, num_cell_types * 0.1)
        fig_height = max(6, num_groups * 1)

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        # Create a diverging color palette
        cmap = sns.diverging_palette(250, 15, s=75, l=40, n=9, center="light", as_cmap=True)

        # Create a custom color normalization object
        divnorm = TwoSlopeNorm(vmin=results['observed_diff'].min(), vcenter=0, vmax=results['observed_diff'].max())

        # Manually plot points using plt.scatter
        for _, row in results.iterrows():
            is_significant = np.isclose(row['adj_p_value'], alpha_threshold, atol=1e-8) or row['adj_p_value'] < alpha_threshold
            is_significant &= abs(row['observed_diff']) > observed_diff_threshold
            plt.scatter(x=[row['cell_type']], y=[row['group2']], 
                        s=[row['size']], 
                        c=[cmap(divnorm(row['observed_diff']))], 
                        edgecolor=edgecolor if is_significant else None,
                        linewidth=linewidth if is_significant else 0)

        # Create color bar
        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=divnorm)
        sm.set_array([])
        cbar = plt.colorbar(sm, orientation="vertical", label="Observed log2(FD)", ax=ax)

        # Create a manual legend for sizes
        for s in [1, .1, .01]:
            plt.scatter([], [], c='gray', s=s * size_scale, label=str(s))
        plt.legend(scatterpoints=1, title='scaled -log10(p-adj)', bbox_to_anchor=(1.25, 0.0), loc="lower left", borderaxespad=0.1, borderpad=0.2, labelspacing=1.5, handletextpad=0.4)

        plt.title(plot_title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.xticks(rotation=rotation)

        # Add gridlines
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.gca().set_axisbelow(True)
        
        plt.show()
    except Exception as e:
        raise RuntimeError(f"An error occurred while plotting: {str(e)}")
