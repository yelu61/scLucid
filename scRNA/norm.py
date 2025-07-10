import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns

def normalize_data(
    adata: sc.AnnData,
    layer: str = "counts",
    target_sum: float = 1e4,
    exclude_highly_expressed: bool = False,
    max_fraction: float = 0.05,
    plot: bool = True,
):
    """
    Normalize and log-transform the single-cell data.

    The normalization step scales the count data to account for differences in sequencing depth between cells.
    Normalize each cell by total counts over all genes, so that every cell has the same total count after normalization.
    The log-transformation reduces the effect of high expression values on downstream analyses.

    This function also plots the distributions of total counts and log-transformed data before and after normalization.

    Args:
        adata (AnnData): AnnData object containing single-cell data.
        layer (str, optional): Name of the layer in adata.layers containing the raw count data. Defaults to 'counts'.
        target_sum (float, optional): Total count to which each cell should be normalized. Defaults to 1e4.
            If None, after normalization, each observation (cell) has a total count equal to the median of total counts for observations (cells) before normalization.
            If choosing target_sum=1e6, this is CPM normalization.
        exclude_highly_expressed (bool, optional): Whether to exclude very highly expressed genes for the computation of the normalization factor (size factor) for each cell. Defaults to False.
        max_fraction (float, optional): If exclude_highly_expressed=True, consider cells as highly expressed that have more counts than max_fraction of the original total counts in at least one cell. Defaults to 0.05.
        plot (bool, optional): Whether to plot the distributions before and after normalization. Defaults to True.

    Returns:
        adata (AnnData): AnnData object with normalized and log-transformed data.

    Raises:
        ValueError: If `target_sum` is negative or `max_fraction` is not between 0 and 1.
    """
    
    # Check input parameter validity
    if target_sum is not None and target_sum <= 0:
        raise ValueError("target_sum must be a positive number.")
    if not 0 < max_fraction < 1:
        raise ValueError("max_fraction must be between 0 and 1 (exclusive).")
    
    # Check if the specified layer exists, if not, create it
    if layer not in adata.layers:
        print(f"Layer '{layer}' not found in adata.layers. Creating a new layer.")
        adata.layers[layer] = adata.X.copy()

    print("Begin data normalization and log-transformation.")

    # Normalize counts per cell
    X_norm = sc.pp.normalize_total(
        adata, 
        layer=layer,
        target_sum=target_sum, 
        exclude_highly_expressed=exclude_highly_expressed,
        max_fraction=max_fraction,
        inplace=False, 
    )['X']
    adata.layers['log1p_norm'] = X_norm

    # Log-transform the normalized data
    sc.pp.log1p(adata, layer = 'log1p_norm')
    
    print("Done.")

    if plot:
        # Visualize the distributions
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 5))
        fig.suptitle("Data Distributions Before and After Normalization", fontsize=16)

        # Plot total counts before normalization
        sns.histplot(
            adata.layers[layer].sum(1),
            bins=100,
            kde=True,
            ax=axes[0],
            color="navy",
        )
        axes[0].set_title("Total Counts (Before Normalization)", fontsize=14)
        axes[0].set_xlabel("Total Counts", fontsize=12)

        # Plot log-transformed values after normalization
        sns.histplot(
            adata.layers["log1p_norm"].sum(1),
            bins=100,
            kde=True,
            ax=axes[1],
            color="crimson",
        )
        axes[1].set_title("Shifted Logarithm (After Normalization)", fontsize=14)
        axes[1].set_xlabel("Log-Transformed Values", fontsize=12)

        plt.tight_layout()
        plt.show()

    return adata