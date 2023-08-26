import matplotlib.pyplot as plt
import pandas as pd

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
