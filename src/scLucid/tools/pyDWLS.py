import numpy as np
import pandas as pd
from scipy import optimize
from scipy.stats import trim_mean
from sklearn.linear_model import Ridge
import warnings
warnings.filterwarnings('ignore')


class DWLS:
    """
    Dampened Weighted Least Squares (DWLS) for gene expression deconvolution
    
    This class implements the DWLS algorithm for estimating cell-type composition
    from bulk RNA-seq data using single-cell RNA-seq reference data.
    """
    
    def __init__(self, signature_matrix=None, bulk_data=None):
        """
        Initialize DWLS deconvolution object
        
        Parameters:
        -----------
        signature_matrix : pd.DataFrame, optional
            Gene expression signature matrix (genes x cell types)
        bulk_data : pd.DataFrame, optional
            Bulk RNA-seq data (genes x samples)
        """
        self.signature_matrix = signature_matrix
        self.bulk_data = bulk_data
        self.results = None
        
    
    def build_signature_matrix(self, sc_data, cell_type_labels, genes_to_use=None, 
                               trim_percent=0.0, min_cells=10):
        """
        Build signature matrix from single-cell data
        
        Parameters:
        -----------
        sc_data : pd.DataFrame
            Single-cell expression data (genes x cells)
        cell_type_labels : pd.Series or array-like
            Cell type labels for each cell
        genes_to_use : list, optional
            Specific genes to use. If None, uses all genes
        trim_percent : float, default=0.0
            Percentage to trim from both ends when calculating mean (0-0.5)
        min_cells : int, default=10
            Minimum number of cells required per cell type
        
        Returns:
        --------
        pd.DataFrame : Signature matrix (genes x cell types)
        """
        print("Building signature matrix from single-cell data...")
        
        # Convert to DataFrame if necessary
        if not isinstance(sc_data, pd.DataFrame):
            sc_data = pd.DataFrame(sc_data)
        
        if not isinstance(cell_type_labels, pd.Series):
            cell_type_labels = pd.Series(cell_type_labels, index=sc_data.columns)
        
        # Filter genes if specified
        if genes_to_use is not None:
            sc_data = sc_data.loc[sc_data.index.isin(genes_to_use)]
        
        # Get unique cell types
        cell_types = cell_type_labels.unique()
        
        # Filter cell types with insufficient cells
        cell_type_counts = cell_type_labels.value_counts()
        valid_cell_types = cell_type_counts[cell_type_counts >= min_cells].index
        
        print(f"Found {len(valid_cell_types)} cell types with >= {min_cells} cells")
        
        # Build signature matrix
        signature_list = []
        
        for ct in valid_cell_types:
            # Get cells of this type
            ct_cells = cell_type_labels[cell_type_labels == ct].index
            ct_data = sc_data[ct_cells]
            
            # Calculate trimmed mean for each gene
            if trim_percent > 0:
                ct_signature = ct_data.apply(
                    lambda x: trim_mean(x, proportiontocut=trim_percent), 
                    axis=1
                )
            else:
                ct_signature = ct_data.mean(axis=1)
            
            signature_list.append(ct_signature)
        
        # Combine into signature matrix
        signature_matrix = pd.concat(signature_list, axis=1)
        signature_matrix.columns = valid_cell_types
        
        # Remove genes with zero variance
        gene_vars = signature_matrix.var(axis=1)
        signature_matrix = signature_matrix[gene_vars > 0]
        
        print(f"Signature matrix shape: {signature_matrix.shape}")
        
        self.signature_matrix = signature_matrix
        return signature_matrix
    
    
    def select_genes(self, sc_data, cell_type_labels, n_genes=50, 
                     method='ratio', log_transform=True):
        """
        Select marker genes for each cell type
        
        Parameters:
        -----------
        sc_data : pd.DataFrame
            Single-cell expression data (genes x cells)
        cell_type_labels : pd.Series or array-like
            Cell type labels for each cell
        n_genes : int, default=50
            Number of genes to select per cell type
        method : str, default='ratio'
            Method for gene selection ('ratio', 'difference', 'fold_change')
        log_transform : bool, default=True
            Whether to log-transform the data before selection
        
        Returns:
        --------
        list : Selected gene names
        """
        print(f"Selecting marker genes using {method} method...")
        
        if not isinstance(cell_type_labels, pd.Series):
            cell_type_labels = pd.Series(cell_type_labels, index=sc_data.columns)
        
        # Log transform if needed
        if log_transform:
            sc_data_trans = np.log2(sc_data + 1)
        else:
            sc_data_trans = sc_data.copy()
        
        cell_types = cell_type_labels.unique()
        selected_genes = set()
        
        for ct in cell_types:
            ct_cells = cell_type_labels[cell_type_labels == ct].index
            other_cells = cell_type_labels[cell_type_labels != ct].index
            
            ct_mean = sc_data_trans[ct_cells].mean(axis=1)
            other_mean = sc_data_trans[other_cells].mean(axis=1)
            
            if method == 'ratio':
                # Use ratio of means
                scores = ct_mean / (other_mean + 1e-10)
            elif method == 'difference':
                # Use difference of means
                scores = ct_mean - other_mean
            elif method == 'fold_change':
                # Use log fold change
                scores = np.log2((ct_mean + 1) / (other_mean + 1))
            else:
                raise ValueError(f"Unknown method: {method}")
            
            # Select top genes
            top_genes = scores.nlargest(n_genes).index.tolist()
            selected_genes.update(top_genes)
        
        selected_genes = list(selected_genes)
        print(f"Selected {len(selected_genes)} unique marker genes")
        
        return selected_genes
    
    
    def solve_dampened_wls(self, S, b, params=None):
        """
        Solve dampened weighted least squares problem
        
        Parameters:
        -----------
        S : np.ndarray
            Signature matrix (genes x cell types)
        b : np.ndarray
            Bulk expression vector (genes,)
        params : dict, optional
            Additional parameters for optimization
        
        Returns:
        --------
        np.ndarray : Estimated cell-type proportions
        """
        n_genes, n_celltypes = S.shape
        
        if params is None:
            params = {}
        
        # Set default parameters
        dampen_factor = params.get('dampen_factor', 1.0)
        use_nonneg = params.get('use_nonneg', True)
        normalize = params.get('normalize', True)
        
        # Calculate gene-specific weights
        # Weight inversely proportional to expression level to dampen highly expressed genes
        gene_means = S.mean(axis=1)
        weights = 1.0 / (gene_means ** dampen_factor + 1e-10)
        
        # Normalize weights
        weights = weights / weights.sum() * n_genes
        
        # Create weighted matrices
        W = np.diag(np.sqrt(weights))
        S_weighted = W @ S
        b_weighted = W @ b
        
        # Solve weighted least squares with non-negativity constraint
        if use_nonneg:
            def objective(x):
                return np.sum((S_weighted @ x - b_weighted) ** 2)
            
            # Non-negative constraint
            bounds = [(0, None) for _ in range(n_celltypes)]
            
            # Initial guess
            x0 = np.ones(n_celltypes) / n_celltypes
            
            # Optimize
            result = optimize.minimize(
                objective, 
                x0, 
                method='L-BFGS-B',
                bounds=bounds,
                options={'maxiter': 1000}
            )
            
            proportions = result.x
        else:
            # Standard weighted least squares
            proportions = np.linalg.lstsq(S_weighted, b_weighted, rcond=None)[0]
            proportions = np.maximum(proportions, 0)  # Ensure non-negativity
        
        # Normalize to sum to 1 if requested
        if normalize and proportions.sum() > 0:
            proportions = proportions / proportions.sum()
        
        return proportions
    
    
    def solve_dwls_single(self, bulk_sample, params=None):
        """
        Solve DWLS for a single bulk sample
        
        Parameters:
        -----------
        bulk_sample : pd.Series or np.ndarray
            Bulk expression data for one sample
        params : dict, optional
            Additional parameters
        
        Returns:
        --------
        pd.Series : Estimated cell-type proportions
        """
        if self.signature_matrix is None:
            raise ValueError("Signature matrix not set. Please build or provide one.")
        
        # Get common genes
        if isinstance(bulk_sample, pd.Series):
            common_genes = self.signature_matrix.index.intersection(bulk_sample.index)
        else:
            common_genes = self.signature_matrix.index
        
        if len(common_genes) == 0:
            raise ValueError("No common genes between signature and bulk data")
        
        # Extract data for common genes
        S = self.signature_matrix.loc[common_genes].values
        if isinstance(bulk_sample, pd.Series):
            b = bulk_sample.loc[common_genes].values
        else:
            b = bulk_sample
        
        # Solve DWLS
        proportions = self.solve_dampened_wls(S, b, params)
        
        # Return as Series
        result = pd.Series(proportions, index=self.signature_matrix.columns)
        
        return result
    
    
    def deconvolve(self, bulk_data=None, params=None, verbose=True):
        """
        Deconvolve bulk RNA-seq data
        
        Parameters:
        -----------
        bulk_data : pd.DataFrame, optional
            Bulk RNA-seq data (genes x samples). If None, uses self.bulk_data
        params : dict, optional
            Parameters for DWLS solver
        verbose : bool, default=True
            Whether to print progress
        
        Returns:
        --------
        pd.DataFrame : Cell-type proportions (samples x cell types)
        """
        if bulk_data is None:
            bulk_data = self.bulk_data
        
        if bulk_data is None:
            raise ValueError("No bulk data provided")
        
        if self.signature_matrix is None:
            raise ValueError("Signature matrix not set")
        
        if verbose:
            print(f"Deconvolving {bulk_data.shape[1]} bulk samples...")
        
        # Deconvolve each sample
        results_list = []
        
        for i, sample_id in enumerate(bulk_data.columns):
            if verbose and (i + 1) % 10 == 0:
                print(f"  Processed {i + 1}/{bulk_data.shape[1]} samples")
            
            bulk_sample = bulk_data[sample_id]
            proportions = self.solve_dwls_single(bulk_sample, params)
            results_list.append(proportions)
        
        # Combine results
        results = pd.concat(results_list, axis=1).T
        results.index = bulk_data.columns
        
        self.results = results
        
        if verbose:
            print("Deconvolution complete!")
        
        return results
    
    
    def cross_validate(self, sc_data, cell_type_labels, n_folds=5, 
                       n_cells_per_bulk=100, random_state=42):
        """
        Perform cross-validation by creating pseudo-bulk samples
        
        Parameters:
        -----------
        sc_data : pd.DataFrame
            Single-cell expression data (genes x cells)
        cell_type_labels : pd.Series
            Cell type labels
        n_folds : int, default=5
            Number of cross-validation folds
        n_cells_per_bulk : int, default=100
            Number of cells to combine per pseudo-bulk sample
        random_state : int, default=42
            Random seed
        
        Returns:
        --------
        dict : Cross-validation results including correlations and errors
        """
        np.random.seed(random_state)
        
        print(f"Performing {n_folds}-fold cross-validation...")
        
        # Convert labels to Series if needed
        if not isinstance(cell_type_labels, pd.Series):
            cell_type_labels = pd.Series(cell_type_labels, index=sc_data.columns)
        
        all_cells = sc_data.columns.tolist()
        n_cells = len(all_cells)
        
        # Shuffle cells
        shuffled_cells = np.random.permutation(all_cells)
        fold_size = n_cells // n_folds
        
        correlations = []
        rmse_values = []
        
        for fold in range(n_folds):
            print(f"\nFold {fold + 1}/{n_folds}")
            
            # Split into training and test
            test_start = fold * fold_size
            test_end = (fold + 1) * fold_size if fold < n_folds - 1 else n_cells
            
            test_cells = shuffled_cells[test_start:test_end]
            train_cells = np.setdiff1d(shuffled_cells, test_cells)
            
            # Build signature from training data
            train_data = sc_data[train_cells]
            train_labels = cell_type_labels[train_cells]
            
            self.build_signature_matrix(train_data, train_labels)
            
            # Create pseudo-bulk from test data
            test_data = sc_data[test_cells]
            test_labels = cell_type_labels[test_cells]
            
            # Randomly sample cells to create pseudo-bulk
            n_pseudo_bulks = len(test_cells) // n_cells_per_bulk
            
            fold_corrs = []
            fold_rmses = []
            
            for i in range(n_pseudo_bulks):
                # Sample cells
                sampled_indices = np.random.choice(
                    len(test_cells), 
                    size=n_cells_per_bulk, 
                    replace=False
                )
                sampled_cells = test_cells[sampled_indices]
                
                # Create pseudo-bulk
                pseudo_bulk = test_data[sampled_cells].sum(axis=1)
                
                # True proportions
                true_props = test_labels[sampled_cells].value_counts()
                true_props = true_props / true_props.sum()
                
                # Deconvolve
                estimated_props = self.solve_dwls_single(pseudo_bulk)
                
                # Align cell types
                all_celltypes = self.signature_matrix.columns
                true_aligned = pd.Series(0.0, index=all_celltypes)
                true_aligned[true_props.index] = true_props
                
                # Calculate metrics
                corr = np.corrcoef(true_aligned.values, estimated_props.values)[0, 1]
                rmse = np.sqrt(np.mean((true_aligned.values - estimated_props.values) ** 2))
                
                fold_corrs.append(corr)
                fold_rmses.append(rmse)
            
            correlations.append(np.mean(fold_corrs))
            rmse_values.append(np.mean(fold_rmses))
            
            print(f"  Mean correlation: {correlations[-1]:.4f}")
            print(f"  Mean RMSE: {rmse_values[-1]:.4f}")
        
        results = {
            'mean_correlation': np.mean(correlations),
            'std_correlation': np.std(correlations),
            'mean_rmse': np.mean(rmse_values),
            'std_rmse': np.std(rmse_values),
            'fold_correlations': correlations,
            'fold_rmses': rmse_values
        }
        
        print(f"\nOverall CV Results:")
        print(f"  Correlation: {results['mean_correlation']:.4f} ± {results['std_correlation']:.4f}")
        print(f"  RMSE: {results['mean_rmse']:.4f} ± {results['std_rmse']:.4f}")
        
        return results


# Utility functions

def normalize_data(data, method='cpm'):
    """
    Normalize gene expression data
    
    Parameters:
    -----------
    data : pd.DataFrame
        Expression data (genes x samples)
    method : str, default='cpm'
        Normalization method ('cpm', 'tpm', 'log_cpm')
    
    Returns:
    --------
    pd.DataFrame : Normalized data
    """
    if method == 'cpm':
        # Counts per million
        normalized = data / data.sum(axis=0) * 1e6
    elif method == 'tpm':
        # Transcripts per million (assuming counts)
        normalized = data / data.sum(axis=0) * 1e6
    elif method == 'log_cpm':
        # Log CPM
        normalized = np.log2(data / data.sum(axis=0) * 1e6 + 1)
    else:
        raise ValueError(f"Unknown normalization method: {method}")
    
    return normalized


def filter_genes(data, min_cells=10, min_expression=1):
    """
    Filter genes based on detection criteria
    
    Parameters:
    -----------
    data : pd.DataFrame
        Expression data (genes x samples)
    min_cells : int
        Minimum number of cells/samples expressing the gene
    min_expression : float
        Minimum expression threshold
    
    Returns:
    --------
    pd.DataFrame : Filtered data
    """
    # Count cells/samples with expression above threshold
    n_detected = (data > min_expression).sum(axis=1)
    
    # Keep genes detected in sufficient cells
    keep_genes = n_detected >= min_cells
    
    filtered_data = data[keep_genes]
    
    print(f"Kept {filtered_data.shape[0]}/{data.shape[0]} genes")
    
    return filtered_data


def create_pseudo_bulk(sc_data, cell_labels, n_cells=100, random_state=None):
    """
    Create pseudo-bulk sample by aggregating single cells
    
    Parameters:
    -----------
    sc_data : pd.DataFrame
        Single-cell data (genes x cells)
    cell_labels : pd.Series
        Cell type labels
    n_cells : int
        Number of cells to aggregate
    random_state : int, optional
        Random seed
    
    Returns:
    --------
    tuple : (pseudo_bulk expression, true proportions)
    """
    if random_state is not None:
        np.random.seed(random_state)
    
    # Sample cells
    sampled_cells = np.random.choice(sc_data.columns, size=n_cells, replace=False)
    
    # Create pseudo-bulk
    pseudo_bulk = sc_data[sampled_cells].sum(axis=1)
    
    # Calculate true proportions
    true_props = cell_labels[sampled_cells].value_counts()
    true_props = true_props / true_props.sum()
    
    return pseudo_bulk, true_props


# Example usage
def example_usage():
    """
    Example demonstrating how to use the DWLS class
    """
    print("=" * 60)
    print("DWLS Python Implementation - Example Usage")
    print("=" * 60)
    
    # Generate synthetic data for demonstration
    np.random.seed(42)
    
    # Create synthetic single-cell data
    n_genes = 2000
    n_cells = 500
    n_celltypes = 5
    
    print("\n1. Generating synthetic single-cell data...")
    
    # Gene names
    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    
    # Cell type labels
    cell_types = [f"CellType_{i}" for i in range(n_celltypes)]
    cells_per_type = n_cells // n_celltypes
    cell_labels = np.repeat(cell_types, cells_per_type)
    
    # Generate cell-type-specific expression
    sc_data_list = []
    for i in range(n_celltypes):
        # Base expression
        base_expr = np.random.negative_binomial(5, 0.3, size=(n_genes, cells_per_type))
        
        # Add cell-type-specific genes (marker genes)
        marker_start = i * 50
        marker_end = (i + 1) * 50
        base_expr[marker_start:marker_end, :] *= 5
        
        sc_data_list.append(base_expr)
    
    sc_data = np.concatenate(sc_data_list, axis=1)
    sc_data = pd.DataFrame(sc_data, index=gene_names)
    cell_labels = pd.Series(cell_labels, index=sc_data.columns)
    
    print(f"   Single-cell data shape: {sc_data.shape}")
    print(f"   Cell types: {cell_types}")
    
    # Initialize DWLS
    print("\n2. Initializing DWLS...")
    dwls = DWLS()
    
    # Select marker genes
    print("\n3. Selecting marker genes...")
    marker_genes = dwls.select_genes(
        sc_data, 
        cell_labels, 
        n_genes=50,
        method='ratio'
    )
    
    # Build signature matrix
    print("\n4. Building signature matrix...")
    signature = dwls.build_signature_matrix(
        sc_data, 
        cell_labels,
        genes_to_use=marker_genes
    )
    
    print(f"   Signature matrix shape: {signature.shape}")
    print(f"   Cell types in signature: {list(signature.columns)}")
    
    # Create pseudo-bulk samples
    print("\n5. Creating pseudo-bulk samples...")
    n_bulk_samples = 20
    bulk_data_list = []
    true_props_list = []
    
    for i in range(n_bulk_samples):
        pseudo_bulk, true_props = create_pseudo_bulk(
            sc_data, 
            cell_labels, 
            n_cells=100,
            random_state=i
        )
        bulk_data_list.append(pseudo_bulk)
        true_props_list.append(true_props)
    
    bulk_data = pd.concat(bulk_data_list, axis=1)
    bulk_data.columns = [f"Sample_{i}" for i in range(n_bulk_samples)]
    
    print(f"   Bulk data shape: {bulk_data.shape}")
    
    # Deconvolve
    print("\n6. Performing deconvolution...")
    results = dwls.deconvolve(bulk_data, verbose=True)
    
    print(f"\n   Results shape: {results.shape}")
    print(f"\n   First few samples:")
    print(results.head())
    
    # Calculate accuracy
    print("\n7. Evaluating accuracy...")
    correlations = []
    for i, sample_id in enumerate(results.index):
        true_props = true_props_list[i]
        estimated_props = results.loc[sample_id]
        
        # Align cell types
        all_types = results.columns
        true_aligned = pd.Series(0.0, index=all_types)
        true_aligned[true_props.index] = true_props
        
        corr = np.corrcoef(true_aligned.values, estimated_props.values)[0, 1]
        correlations.append(corr)
    
    print(f"   Mean correlation: {np.mean(correlations):.4f} ± {np.std(correlations):.4f}")
    
    # Cross-validation
    print("\n8. Performing cross-validation...")
    cv_results = dwls.cross_validate(
        sc_data, 
        cell_labels, 
        n_folds=3,
        n_cells_per_bulk=100
    )
    
    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)
    
    return dwls, results, cv_results


if __name__ == "__main__":
    # Run example
    dwls_model, deconv_results, cv_results = example_usage()