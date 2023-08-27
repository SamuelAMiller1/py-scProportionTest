import matplotlib.pyplot as plt

def multiple_point_range_plot(results, cell_type, reference_treatment,
                                alpha=0.05, fold_difference=1, figsize=(6,6), dot_size=6, significance_color='darkred', 
                                ci_line_width=1.2, plot_title='Log2 Fold Differences Across Treatments', 
                                x_label='Log2 Fold Difference', y_label='Treatments', ascending=True):
    """
    Generate a point-range plot of the permutation test results, including the observed log2 fold differences
    and bootstrapped confidence intervals for a single cell type across different treatments.

    Parameters
    ----------
    results : DataFrame
        A pandas DataFrame with the permutation p-values, observed proportional differences,
        and bootstrapped confidence intervals for each treatment compared to a reference.
    cell_type : str
        The specific cell type for which the comparisons are to be plotted.
    reference_treatment : str
        The treatment that serves as a reference for the comparisons.
    ascending : bool, optional
        Sort treatments based on observed_diff in ascending order if True, else in descending order. Default is True.
    Other parameters are similar to the `point_range_plot` function.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object of the plot.
    """
    
    # Filter results for the specific cell type and where the reference treatment is in group1
    filtered_results = results[(results['cell_type'] == cell_type) & (results['group1'] == reference_treatment)]
    
    # Sort the results based on observed_diff in the specified order
    filtered_results = filtered_results.sort_values(by='observed_diff', ascending=ascending)
    
    try:
        fig, ax = plt.subplots(figsize=figsize)
        
        # Add gridlines
        ax.grid(color='grey', linestyle='-', linewidth=0.5, which='both', axis='x')

        # Add dotted black lines at pos fold_difference and neg fold_difference
        ax.axvline(x=fold_difference, color='black', linestyle='--', linewidth=ci_line_width)
        ax.axvline(x=-fold_difference, color='black', linestyle='--', linewidth=ci_line_width)

        # Plot the observed log2 fold differences with confidence intervals
        for index, row in filtered_results.iterrows():
            significant = (row['adj_p_value'] < alpha) and (abs(row['observed_diff']) > fold_difference)
            color = significance_color if significant else 'grey'
            label = f'p-adj < {alpha}, |log2FD| > {fold_difference}' if significant else 'Not Significant'
            ax.plot([row['lower_ci'], row['upper_ci']], [row['group2'], row['group2']], color=color, linewidth=ci_line_width)  # Range line
            ax.plot(row['observed_diff'], row['group2'], 'o', color=color, markersize=dot_size, label=label)  # Point

        # Plot labels
        ax.set_title(plot_title)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

        # Add legend without duplicate labels
        handles, labels = ax.get_legend_handles_labels()
        unique_labels = dict(zip(labels, handles))
        ax.legend(unique_labels.values(), unique_labels.keys(), loc='center left', bbox_to_anchor=(1, 0.5))

        
    except Exception as e:
        raise RuntimeError(f"An error occurred during plotting: {e}")
