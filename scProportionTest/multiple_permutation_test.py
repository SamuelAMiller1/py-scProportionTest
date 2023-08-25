import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import current_thread

def bootstrap_ci(adata, group1_cells, group2_cells, cell_type_col, cell_type, n_bootstrap=10000, alpha=0.05):
    """
    Compute the bootstrapped confidence interval for the log2 proportion difference between two groups of cells.

    Parameters:
    - adata (AnnData): An AnnData object containing cell-wise observations and annotations.
    - group1_cells (list): List of cell identifiers belonging to group 1.
    - group2_cells (list): List of cell identifiers belonging to group 2.
    - cell_type_col (str): Column name in `adata.obs` containing cell type information.
    - cell_type (str): Specific cell type for which confidence interval should be computed.
    - n_bootstrap (int, optional): Number of bootstrap iterations. Default is 10000.
    - alpha (float, optional): Significance level for confidence interval computation. Default is 0.05.

    Returns:
    - tuple: Lower and upper bounds of the bootstrapped confidence interval.

    This function resamples the cells of each group with replacement for `n_bootstrap` times. 
    For each iteration, it computes the log2 proportion difference for the given cell type between the two groups.
    At the end, the confidence interval is computed based on the desired significance level.
    """
    if not cell_type_col in adata.obs.columns:
        raise ValueError(f"'{cell_type_col}' not found in 'adata.obs'.")
    
    bootstrapped_diffs = []
    for _ in range(n_bootstrap):
        bootstrap_group1_cells = np.random.choice(group1_cells, size=len(group1_cells), replace=True)
        bootstrap_group2_cells = np.random.choice(group2_cells, size=len(group2_cells), replace=True)
        prop_group1 = (adata.obs.loc[bootstrap_group1_cells, cell_type_col] == cell_type).mean()
        prop_group2 = (adata.obs.loc[bootstrap_group2_cells, cell_type_col] == cell_type).mean()
        bootstrapped_diff = np.log2(prop_group1 / prop_group2) if prop_group2 > 0 else np.nan
        bootstrapped_diffs.append(bootstrapped_diff)
    return np.nanpercentile(bootstrapped_diffs, [alpha / 2 * 100, (1 - alpha / 2) * 100])

def single_comparison(adata, group1, group2, group_col, cell_type_col, nperm, alpha, n_bootstrap):
    """
    Perform single comparison between two groups to compute observed differences and permutation-based p-values.

    Parameters:
    - adata (AnnData): An AnnData object containing cell-wise observations and annotations.
    - group1 (str): Name of the first group.
    - group2 (str): Name of the second group.
    - group_col (str): Column name in `adata.obs` containing group information.
    - cell_type_col (str): Column name in `adata.obs` containing cell type information.
    - nperm (int): Number of permutations to be performed for p-value computation.
    - alpha (float): Significance level for confidence interval computation.
    - n_bootstrap (int): Number of bootstrap iterations for confidence interval computation.

    Returns:
    - DataFrame: Contains observed differences, p-values, and confidence intervals for each cell type between the two groups.

    This function first computes the observed log2 proportion differences for each cell type between the two groups. 
    It then uses bootstrapping to estimate confidence intervals for these observed differences.
    Additionally, it performs permutation tests to compute p-values for the observed differences.
    """
    if not group_col in adata.obs.columns:
        raise ValueError(f"'{group_col}' not found in 'adata.obs'.")
    
    if not cell_type_col in adata.obs.columns:
        raise ValueError(f"'{cell_type_col}' not found in 'adata.obs'.")

    print(f"Running for group2: {group2} on {current_thread().name}")
    
    cells_group1 = adata.obs_names[adata.obs[group_col] == group1].tolist()
    cells_group2 = adata.obs_names[adata.obs[group_col] == group2].tolist()
    all_cells = cells_group1 + cells_group2
    cell_types = adata.obs[cell_type_col].unique()
    observed_diffs = {}
    bootstrapped_cis = {}
    null_diffs_collection = {cell_type: [] for cell_type in cell_types}
    for cell_type in cell_types:
        prop_group1 = (adata.obs.loc[cells_group1, cell_type_col] == cell_type).mean()
        prop_group2 = (adata.obs.loc[cells_group2, cell_type_col] == cell_type).mean()
        observed_diffs[cell_type] = np.log2(prop_group2 / prop_group1) if prop_group1 > 0 else np.nan
        lower_ci, upper_ci = bootstrap_ci(adata, cells_group1, cells_group2, cell_type_col, cell_type, n_bootstrap, alpha)
        bootstrapped_cis[cell_type] = (lower_ci, upper_ci)
        null_diffs = []
        for _ in range(nperm):
            np.random.shuffle(all_cells)
            perm_group1 = all_cells[:len(cells_group1)]
            perm_group2 = all_cells[len(cells_group1):]
            perm_prop_group1 = (adata.obs.loc[perm_group1, cell_type_col] == cell_type).mean()
            perm_prop_group2 = (adata.obs.loc[perm_group2, cell_type_col] == cell_type).mean()
            perm_diff = np.log2(perm_prop_group1 / perm_prop_group2) if perm_prop_group2 > 0 else np.nan
            null_diffs.append(perm_diff)
        null_diffs_collection[cell_type] = null_diffs
    p_values = {}
    for cell_type, observed_diff in observed_diffs.items():
        null_diffs = null_diffs_collection[cell_type]
        p_value = (np.sum(np.abs(null_diffs) >= np.abs(observed_diff)) + 1) / (nperm + 1)
        p_values[cell_type] = p_value
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


def multiple_permutation_test(adata, group1, group2_list, group_col='group', cell_type_col='cell_type', nperm=10000, alpha=0.05, n_bootstrap=10000, p_adjust_method='fdr_bh'):
    """
    Conduct permutation tests for multiple group comparisons in parallel.

    Parameters:
    - adata (AnnData): An AnnData object containing cell-wise observations and annotations.
    - group1 (str): Name of the primary group against which all groups in `group2_list` will be compared.
    - group2_list (list): List of group names to be compared against `group1`.
    - group_col (str, optional): Column name in `adata.obs` containing group information. Default is 'group'.
    - cell_type_col (str, optional): Column name in `adata.obs` containing cell type information. Default is 'cell_type'.
    - nperm (int, optional): Number of permutations to be performed for p-value computation. Default is 10000.
    - alpha (float, optional): Significance level for confidence interval computation. Default is 0.05.
    - n_bootstrap (int, optional): Number of bootstrap iterations for confidence interval computation. Default is 10000.
    - p_adjust_method (str, optional): Method for multiple testing correction. Default is 'fdr_bh'.

    Returns:
    - DataFrame: Consolidated results containing observed differences, p-values, adjusted p-values, and confidence intervals for each cell type across all comparisons.

    This function performs the group comparisons in parallel using multi-threading. For each group in `group2_list`, 
    it compares against `group1` to compute observed differences, confidence intervals, and permutation-based p-values.
    Finally, it applies multiple testing correction to the p-values.
    """
    if not group_col in adata.obs.columns:
        raise ValueError(f"'{group_col}' not found in 'adata.obs'.")
    
    if not cell_type_col in adata.obs.columns:
        raise ValueError(f"'{cell_type_col}' not found in 'adata.obs'.")
    
    if not isinstance(group2_list, list):
        raise ValueError("'group2_list' should be of type list.")
    num_cpus = os.cpu_count()
    print(f"Using {num_cpus} CPUs for parallel processing.")

    all_results = []
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = [executor.submit(single_comparison, adata, group1, group2, group_col, cell_type_col, nperm, alpha, n_bootstrap) for group2 in group2_list]
        
        for future in as_completed(futures):
            if future.exception() is not None:
                raise future.exception()
            result = future.result()
            all_results.append(result)
    
    final_results = pd.concat(all_results, ignore_index=True)
    final_results['adj_p_value'] = multipletests(final_results['p_value'], method=p_adjust_method)[1]
    return final_results
