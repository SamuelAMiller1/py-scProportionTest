import matplotlib.pyplot as plt



def point_range_plot(results, alpha=0.05, fold_difference=1, figsize=(6,6)):
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

    Returns
    -------
    matplotlib.figure.Figure
        The figure object of the plot.
    """

    fig, ax = plt.subplots(figsize=figsize)

    # Add gridlines
    ax.grid(color='grey', linestyle='-', linewidth=0.5, which='both', axis='x')

    # Add dotted black lines at pos fold_difference and neg fold_difference
    ax.axvline(x=fold_difference, color='black', linestyle='--', linewidth=1.2)
    ax.axvline(x=-fold_difference, color='black', linestyle='--', linewidth=1.2)

    # Plot the observed log2 fold differences with confidence intervals
    for index, row in results.iterrows():
        significant = (row['adj_p_value'] < alpha) and (abs(row['observed_diff']) > fold_difference)
        color = 'darkred' if significant else 'grey'
        label = f'p-adj < {alpha}, |log2FD| > {fold_difference}' if significant else 'Not Significant'
        ax.plot([row['lower_ci'], row['upper_ci']], [row['cell_type'], row['cell_type']], color=color)  # Range line
        ax.plot(row['observed_diff'], row['cell_type'], 'o', color=color, label=label)  # Point

    # Label the plot
    ax.set_title('Permutation Test Results: Observed Log2 Fold Differences')
    ax.set_xlabel('Log2 Fold Difference')
    ax.set_ylabel('Cell Type')

    # Add legend without duplicate labels
    handles, labels = ax.get_legend_handles_labels()
    unique_labels = dict(zip(labels, handles))
    ax.legend(unique_labels.values(), unique_labels.keys(), loc='center left', bbox_to_anchor=(1, 0.5))

    return fig