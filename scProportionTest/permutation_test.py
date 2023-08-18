import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm



def bootstrap_ci(adata, group1_cells, group2_cells, cell_type_col, cell_type, n_bootstrap=10000, alpha=0.05):
    """
    Calculate the bootstrapped confidence intervals for the observed proportional difference 
    for a specific cell type between two groups.

    Parameters:
        adata (AnnData): Annotated data matrix containing cell data.
        group1_cells (list): List of cell names or indices for group 1.
        group2_cells (list): List of cell names or indices for group 2.
        cell_type_col (str): Column name in adata.obs containing cell type labels.
        cell_type (str): Specific cell type for which to calculate the confidence interval.
        n_bootstrap (int, optional): Number of bootstrap iterations. Default is 10000.
        alpha (float, optional): Significance level for the confidence interval. Default is 0.05.

    Returns:
        tuple: Lower and upper bounds of the bootstrapped confidence interval.
    """
    # Bootstrapped confidence interval calculation
    bootstrapped_diffs = []
    for _ in range(n_bootstrap):
        bootstrap_group1_cells = np.random.choice(group1_cells, size=len(group1_cells), replace=True)
        bootstrap_group2_cells = np.random.choice(group2_cells, size=len(group2_cells), replace=True)

        prop_group1 = (adata.obs.loc[bootstrap_group1_cells, cell_type_col] == cell_type).mean()
        prop_group2 = (adata.obs.loc[bootstrap_group2_cells, cell_type_col] == cell_type).mean()

        bootstrapped_diff = np.log2(prop_group1 / prop_group2) if prop_group2 > 0 else np.nan
        bootstrapped_diffs.append(bootstrapped_diff)

    return np.nanpercentile(bootstrapped_diffs, [alpha / 2 * 100, (1 - alpha / 2) * 100])



def permutation_test(adata, group1, group2, group_col='group', cell_type_col='cell_type', nperm=10000, alpha=0.05, n_bootstrap=10000, verbose=True):
    """
    Perform a permutation test to evaluate the differences in cell type proportions between two groups.
    Calculates p-values, adjusted p-values, and bootstrapped confidence intervals for the observed
    log2 fold differences in proportions for each cell type. The log2 fold difference is calculated
    as log2(prop_group1 / prop_group2), where group1 serves as the reference group.
    
    Parameters:
        adata (AnnData): Annotated data matrix containing cell data.
        group1 (str): Name of the reference group (numerator in the log2 fold difference calculation).
        group2 (str): Name of the group to be compared against the reference group (denominator in the log2 fold difference calculation).
        group_col (str, optional): Column name in adata.obs containing group labels. Default is 'group'.
        cell_type_col (str, optional): Column name in adata.obs containing cell type labels. Default is 'cell_type'.
        nperm (int, optional): Number of permutation iterations. Default is 10000.
        alpha (float, optional): Significance level for the confidence interval. Default is 0.05.
        n_bootstrap (int, optional): Number of bootstrap iterations. Default is 10000.
        verbose (bool, optional): If True, displays a progress bar. Default is True.
    
    Returns:
        pd.DataFrame: DataFrame containing the results, including cell type, p-value, adjusted p-value,
                      observed log2 fold difference, and lower and upper bounds of the bootstrapped confidence interval.
    """
    
    cells_group1 = adata.obs_names[adata.obs[group_col] == group1].tolist()
    cells_group2 = adata.obs_names[adata.obs[group_col] == group2].tolist()
    all_cells = cells_group1 + cells_group2

    cell_types = adata.obs[cell_type_col].unique()
    observed_diffs = {}
    bootstrapped_cis = {}
    null_diffs_collection = {cell_type: [] for cell_type in cell_types}
    
    # Wrap the loop over cell types with tqdm for progress bar
    for cell_type in tqdm(cell_types, desc="Processing cell types", disable=not verbose):
        prop_group1 = (adata.obs.loc[cells_group1, cell_type_col] == cell_type).mean()
        prop_group2 = (adata.obs.loc[cells_group2, cell_type_col] == cell_type).mean()
        observed_diffs[cell_type] = np.log2(prop_group1 / prop_group2) if prop_group2 > 0 else np.nan

        # Calculate bootstrapped confidence interval
        lower_ci, upper_ci = bootstrap_ci(adata, cells_group1, cells_group2, cell_type_col, cell_type, n_bootstrap, alpha)
        bootstrapped_cis[cell_type] = (lower_ci, upper_ci)

        # Perform the permutations
        null_diffs = []
        for _ in range(nperm):
            np.random.shuffle(all_cells)
            perm_group1 = all_cells[:len(cells_group1)]
            perm_group2 = all_cells[len(cells_group1):]

            perm_prop_group1 = (adata.obs.loc[perm_group1, cell_type_col] == cell_type).mean()
            perm_prop_group2 = (adata.obs.loc[perm_group2, cell_type_col] == cell_type).mean()
            perm_diff = np.log2(perm_prop_group1 / perm_prop_group2) if perm_prop_group2 > 0 else np.nan
            null_diffs.append(perm_diff)

        # Store null differences for this cell type
        null_diffs_collection[cell_type] = null_diffs

    # Calculate p-values based on null distributions
    p_values = {}
    for cell_type, observed_diff in observed_diffs.items():
        null_diffs = null_diffs_collection[cell_type]
        p_value = (np.sum(np.abs(null_diffs) >= np.abs(observed_diff)) + 1) / (nperm + 1)
        p_values[cell_type] = p_value

    # Create DataFrame to hold results
    results = pd.DataFrame({
        'cell_type': list(p_values.keys()),
        'p_value': list(p_values.values()),
        'observed_diff': [observed_diffs[cell_type] for cell_type in p_values.keys()],
        'lower_ci': [bootstrapped_cis[cell_type][0] for cell_type in p_values.keys()],
        'upper_ci': [bootstrapped_cis[cell_type][1] for cell_type in p_values.keys()],
    })

    # Adjust the p-values using Benjamini-Hochberg procedure
    results['adj_p_value'] = multipletests(results['p_value'], method='fdr_bh')[1]

    return results
