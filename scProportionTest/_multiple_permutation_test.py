import numpy as np
import pandas as pd
import os
from statsmodels.stats.multitest import multipletests
from concurrent.futures import ThreadPoolExecutor, as_completed

def _bootstrap_ci(adata, group1_cells, group2_cells, cell_type_col, cell_type, n_bootstrap=10000, alpha=0.05):
    """
    Compute the bootstrapped confidence interval for the log2 proportion difference between two groups of cells.

    Parameters:
    - adata: An AnnData object containing cell-wise observations and annotations.
    - group1_cells: List of cell identifiers belonging to group 1.
    - group2_cells: List of cell identifiers belonging to group 2.
    - cell_type_col: Column name in `adata.obs` containing cell type information.
    - cell_type: Specific cell type for which confidence interval should be computed.
    - n_bootstrap: Number of bootstrap iterations (default is 10000).
    - alpha: Significance level for confidence interval computation (default is 0.05).

    Returns:
    - tuple: Lower and upper bounds of the bootstrapped confidence interval.
    """
    
    # Initializing a list to store bootstrapped differences
    bootstrapped_diffs = []
    for _ in range(n_bootstrap):
        # Resampling the cells for each group with replacement
        bootstrap_group1_cells = np.random.choice(group1_cells, size=len(group1_cells), replace=True)
        bootstrap_group2_cells = np.random.choice(group2_cells, size=len(group2_cells), replace=True)
        
        # Calculating proportions of the specific cell type for each bootstrapped group
        prop_group1 = (adata.obs.loc[bootstrap_group1_cells, cell_type_col] == cell_type).mean()
        prop_group2 = (adata.obs.loc[bootstrap_group2_cells, cell_type_col] == cell_type).mean()

        # Handling cases where proportions are zero to avoid NaN values
        if prop_group1 == 0 and prop_group2 == 0:
            bootstrapped_diff = 0
        elif prop_group1 == 0:
            bootstrapped_diff = np.inf
        elif prop_group2 == 0:
            bootstrapped_diff = -np.inf
        else:
            bootstrapped_diff = np.log2(prop_group2 / prop_group1)
        
        bootstrapped_diffs.append(bootstrapped_diff)
        
    # Calculating the confidence interval based on the specified significance level
    return np.percentile(bootstrapped_diffs, [alpha / 2 * 100, (1 - alpha / 2) * 100])

def _single_comparison(adata, group1, group2, group_col, cell_type_col, nperm, alpha, n_bootstrap):
    """
    Perform a single comparison between two groups to compute observed differences and permutation-based p-values.

    Parameters:
    - adata: An AnnData object containing cell-wise observations and annotations.
    - group1: Name of the first group.
    - group2: Name of the second group.
    - group_col: Column name in `adata.obs` containing group information.
    - cell_type_col: Column name in `adata.obs` containing cell type information.
    - nperm: Number of permutations for p-value computation.
    - alpha: Significance level for confidence interval computation.
    - n_bootstrap: Number of bootstrap iterations for confidence interval computation.

    Returns:
    - DataFrame: Contains observed differences, p-values, and confidence intervals for each cell type between the two groups.
    """
    
    # Extracting cells associated with each group
    cells_group1 = adata.obs_names[adata.obs[group_col] == group1].tolist()
    cells_group2 = adata.obs_names[adata.obs[group_col] == group2].tolist()

    # Combining both groups for permutation
    all_cells = cells_group1 + cells_group2

    # Getting unique cell types
    cell_types = adata.obs[cell_type_col].unique()
    
    observed_diffs = {}
    bootstrapped_cis = {}
    null_diffs_collection = {cell_type: [] for cell_type in cell_types}
    
    for cell_type in cell_types:
        # Calculating proportions of the specific cell type for each group
        prop_group1 = (adata.obs.loc[cells_group1, cell_type_col] == cell_type).mean()
        prop_group2 = (adata.obs.loc[cells_group2, cell_type_col] == cell_type).mean()

        # Handling cases where proportions are zero to avoid NaN values
        if prop_group1 == 0 and prop_group2 == 0:
            observed_diff = 0
        elif prop_group1 == 0:
            observed_diff = np.inf
        elif prop_group2 == 0:
            observed_diff = -np.inf
        else:
            observed_diff = np.log2(prop_group2 / prop_group1)

        observed_diffs[cell_type] = observed_diff

        # Computing bootstrapped confidence intervals
        lower_ci, upper_ci = _bootstrap_ci(adata, cells_group1, cells_group2, cell_type_col, cell_type, n_bootstrap, alpha)
        bootstrapped_cis[cell_type] = (lower_ci, upper_ci)

        null_diffs = []
        for _ in range(nperm):
            # Permuting the cells and splitting into two pseudo-groups
            np.random.shuffle(all_cells)
            perm_group1 = all_cells[:len(cells_group1)]
            perm_group2 = all_cells[len(cells_group1):]
            
            # Calculating proportions of the specific cell type for each pseudo-group
            perm_prop_group1 = (adata.obs.loc[perm_group1, cell_type_col] == cell_type).mean()
            perm_prop_group2 = (adata.obs.loc[perm_group2, cell_type_col] == cell_type).mean()

            # Handling cases where proportions are zero to avoid NaN values
            if perm_prop_group1 == 0 and perm_prop_group2 == 0:
                perm_diff = 0
            elif perm_prop_group1 == 0:
                perm_diff = np.inf
            elif perm_prop_group2 == 0:
                perm_diff = -np.inf
            else:
                perm_diff = np.log2(perm_prop_group2 / perm_prop_group1)

            null_diffs.append(perm_diff)
        null_diffs_collection[cell_type] = null_diffs

    # Computing permutation-based p-values
    p_values = {}
    for cell_type, observed_diff in observed_diffs.items():
        null_diffs = null_diffs_collection[cell_type]
        p_value = (np.sum(np.abs(null_diffs) >= np.abs(observed_diff)) + 1) / (nperm + 1)
        p_values[cell_type] = p_value

    # Consolidating the results into a DataFrame
    results = pd.DataFrame({
        'group1': group1,
        'group2': group2,
        'cell_type': list(p_values.keys()),
        'p_value': list(p_values.values()),
        'observed_diff': [observed_diffs[cell_type] for cell_type in p_values.keys()],
        'lower_ci': [bootstrapped_cis[cell_type][0] for cell_type in p_values.keys()],
        'upper_ci': [bootstrapped_cis[cell_type][1] for cell_type in p_values.keys()],
    })
    
    return results

def multiple_permutation_test(adata, group1, group2_list, group_col='group', cell_type_col='cell_type', nperm=10000, alpha=0.05, n_bootstrap=10000, p_adjust_method='fdr_bh', seed=None):
    """
    Conduct permutation tests for multiple group comparisons in parallel.

    Parameters:
    - adata: An AnnData object containing cell-wise observations and annotations.
    - group1: Name of the primary group against which all groups in `group2_list` will be compared.
    - group2_list: List of group names to be compared against `group1`.
    - group_col: Column name in `adata.obs` containing group information (default is 'group').
    - cell_type_col: Column name in `adata.obs` containing cell type information (default is 'cell_type').
    - nperm: Number of permutations for p-value computation (default is 10000).
    - alpha: Significance level for confidence interval computation (default is 0.05).
    - n_bootstrap: Number of bootstrap iterations for confidence interval computation (default is 10000).
    - p_adjust_method: Method for multiple testing correction (default is 'fdr_bh').
    - seed: Seed for reproducibility.

    Returns:
    - DataFrame: Consolidated results containing observed differences, p-values, adjusted p-values, and confidence intervals for each cell type across all comparisons.
    """
    np.random.seed(seed)
    
    # Checking for the presence of required columns in the AnnData object
    if group_col not in adata.obs.columns:
        raise ValueError(f"'{group_col}' not found in 'adata.obs'.")
    if cell_type_col not in adata.obs.columns:
        raise ValueError(f"'{cell_type_col}' not found in 'adata.obs'.")
    if not isinstance(group2_list, list):
        raise ValueError("'group2_list' should be of type list.")

    # Determining the number of available CPUs for parallel processing
    num_cpus = os.cpu_count()
    print(f"Using {num_cpus} CPUs for parallel processing.")

    all_results = []
    # Using multi-threading to perform multiple group comparisons in parallel
    with ThreadPoolExecutor(max_workers=num_cpus) as executor:
        futures = [executor.submit(_single_comparison, adata, group1, group2, group_col, cell_type_col, nperm, alpha, n_bootstrap) for group2 in group2_list]

        # Collecting results from all threads
        for future in as_completed(futures):
            if future.exception() is not None:
                raise future.exception()
            result = future.result()
            all_results.append(result)

    # Concatenating all results into a single DataFrame
    final_results = pd.concat(all_results, ignore_index=True)

    # Applying multiple testing correction to the p-values
    final_results['adj_p_value'] = multipletests(final_results['p_value'], method=p_adjust_method)[1]

    return final_results
