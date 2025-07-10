import matplotlib.pyplot as plt
import scanpy as sc
from scipy import sparse
from typing import Optional, List, Union, Literal

def score_cell_cycle(
    adata: sc.AnnData,
    species: Literal["human", "mouse", "rat"] = "human",
    s_genes: List[str] = None,
    g2m_genes: List[str] = None,
    copy: bool = False,
    plot: bool = True,
    save_dir: Optional[str] = None,
    regress_out: bool = False,
    layer: Optional[str] = "log1p_norm",
) -> sc.AnnData:
    """
    Score cell cycle phases and optionally regress out cell cycle effects.
    
    Args:
        adata: AnnData object
        species: Species of the dataset. Options: "human", "mouse", "rat"
        s_genes: List of S phase marker genes. If provided, overrides the species-specific list.
        g2m_genes: List of G2M phase marker genes. If provided, overrides the species-specific list.
        copy: Whether to return a copy of the AnnData object
        plot: Whether to plot cell cycle scores
        save_dir: Directory to save plots. If None, plots are not saved to disk.
        regress_out: Whether to regress out cell cycle effects
        layer: Layer to use for regression if regress_out=True
        
    Returns:
        AnnData object with cell cycle scores added
        
    Example:
        >>> # For human data
        >>> adata = pp.score_cell_cycle(adata, species="human")
        >>> # For mouse data
        >>> adata = pp.score_cell_cycle(adata, species="mouse")
        >>> # To regress out cell cycle effects
        >>> adata = pp.score_cell_cycle(adata, species="mouse", regress_out=True)
    """
    # Species-specific gene lists
    species_genes = {
        "human": {
            "s_genes": [
                'MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2',
                'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2',
                'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7',
                'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN',
                'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8'
            ],
            "g2m_genes": [
                'HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80',
                'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A',
                'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E',
                'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK',
                'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8',
                'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5',
                'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA'
            ]
        },
        "mouse": {
            "s_genes": [
                'Mcm5', 'Pcna', 'Tyms', 'Fen1', 'Mcm2', 'Mcm4', 'Rrm1', 'Ung', 'Gins2',
                'Mcm6', 'Cdca7', 'Dtl', 'Prim1', 'Uhrf1', 'Cenpu', 'Hells', 'Rfc2',
                'Rpa2', 'Nasp', 'Rad51ap1', 'Gmnn', 'Wdr76', 'Slbp', 'Ccne2', 'Ubr7',
                'Pold3', 'Msh2', 'Atad2', 'Rad51', 'Rrm2', 'Cdc45', 'Cdc6', 'Exo1', 'Tipin',
                'Dscc1', 'Blm', 'Casp8ap2', 'Usp1', 'Clspn', 'Pola1', 'Chaf1b', 'Brip1', 'E2f8'
            ],
            "g2m_genes": [
                'Hmgb2', 'Cdk1', 'Nusap1', 'Ube2c', 'Birc5', 'Tpx2', 'Top2a', 'Ndc80',
                'Cks2', 'Nuf2', 'Cks1b', 'Mki67', 'Tmpo', 'Cenpf', 'Tacc3', 'Fam64a',
                'Smc4', 'Ccnb2', 'Ckap2l', 'Ckap2', 'Aurkb', 'Bub1', 'Kif11', 'Anp32e',
                'Tubb4b', 'Gtse1', 'Kif20b', 'Hjurp', 'Cdca3', 'Hn1', 'Cdc20', 'Ttk',
                'Cdc25c', 'Kif2c', 'Rangap1', 'Ncapd2', 'Dlgap5', 'Cdca2', 'Cdca8',
                'Ect2', 'Kif23', 'Hmmr', 'Aurka', 'Psrc1', 'Anln', 'Lbr', 'Ckap5',
                'Cenpe', 'Ctcf', 'Nek2', 'G2e3', 'Gas2l3', 'Cbx5', 'Cenpa'
            ]
        },
        "rat": {
            "s_genes": [
                'Mcm5', 'Pcna', 'Tyms', 'Fen1', 'Mcm2', 'Mcm4', 'Rrm1', 'Ung', 'Gins2',
                'Mcm6', 'Cdca7', 'Dtl', 'Prim1', 'Uhrf1', 'Cenpu', 'Hells', 'Rfc2',
                'Rpa2', 'Nasp', 'Rad51ap1', 'Gmnn', 'Wdr76', 'Slbp', 'Ccne2', 'Ubr7',
                'Pold3', 'Msh2', 'Atad2', 'Rad51', 'Rrm2', 'Cdc45', 'Cdc6', 'Exo1', 'Tipin',
                'Dscc1', 'Blm', 'Casp8ap2', 'Usp1', 'Clspn', 'Pola1', 'Chaf1b', 'Brip1', 'E2f8'
            ],
            "g2m_genes": [
                'Hmgb2', 'Cdk1', 'Nusap1', 'Ube2c', 'Birc5', 'Tpx2', 'Top2a', 'Ndc80',
                'Cks2', 'Nuf2', 'Cks1b', 'Mki67', 'Tmpo', 'Cenpf', 'Tacc3', 'Fam64a',
                'Smc4', 'Ccnb2', 'Ckap2l', 'Ckap2', 'Aurkb', 'Bub1', 'Kif11', 'Anp32e',
                'Tubb4b', 'Gtse1', 'Kif20b', 'Hjurp', 'Cdca3', 'Hn1', 'Cdc20', 'Ttk',
                'Cdc25c', 'Kif2c', 'Rangap1', 'Ncapd2', 'Dlgap5', 'Cdca2', 'Cdca8',
                'Ect2', 'Kif23', 'Hmmr', 'Aurka', 'Psrc1', 'Anln', 'Lbr', 'Ckap5',
                'Cenpe', 'Ctcf', 'Nek2', 'G2e3', 'Gas2l3', 'Cbx5', 'Cenpa'
            ]
        }
    }
    
    # Check if species is valid
    if species not in species_genes:
        raise ValueError(f"Unknown species: {species}. Valid options are: {', '.join(species_genes.keys())}")
    
    # Use provided gene lists or species-specific defaults
    s_genes = s_genes if s_genes is not None else species_genes[species]["s_genes"]
    g2m_genes = g2m_genes if g2m_genes is not None else species_genes[species]["g2m_genes"]
    
    # Try to detect species if not specified
    if species == "human" and not any(gene in adata.var_names for gene in s_genes[:10]):
        # Check if mouse genes might be present
        mouse_genes = species_genes["mouse"]["s_genes"]
        if any(gene in adata.var_names for gene in mouse_genes[:10]):
            print("Human genes not found, but mouse genes detected. Switching to mouse gene list.")
            s_genes = mouse_genes
            g2m_genes = species_genes["mouse"]["g2m_genes"]
        else:
            # Check if rat genes might be present
            rat_genes = species_genes["rat"]["s_genes"]
            if any(gene in adata.var_names for gene in rat_genes[:10]):
                print("Human genes not found, but rat genes detected. Switching to rat gene list.")
                s_genes = rat_genes
                g2m_genes = species_genes["rat"]["g2m_genes"]
    
    # Filter gene lists to include only genes present in the dataset
    s_genes_in_data = [gene for gene in s_genes if gene in adata.var_names]
    g2m_genes_in_data = [gene for gene in g2m_genes if gene in adata.var_names]
    
    if len(s_genes_in_data) < 5 or len(g2m_genes_in_data) < 5:
        print(f"Warning: Few cell cycle genes found in data (S: {len(s_genes_in_data)}, G2M: {len(g2m_genes_in_data)})")
        print(f"This may indicate that your data is from a different species than '{species}'.")
        print("Consider specifying the correct species or providing custom gene lists.")
        
        # Show some of the genes that weren't found
        missing_s = [gene for gene in s_genes[:10] if gene not in adata.var_names]
        missing_g2m = [gene for gene in g2m_genes[:10] if gene not in adata.var_names]
        if missing_s:
            print(f"Example missing S genes: {', '.join(missing_s[:5])}")
        if missing_g2m:
            print(f"Example missing G2M genes: {', '.join(missing_g2m[:5])}")
            
        # Suggest solutions
        print("\nPotential solutions:")
        print("1. Try a different species: pp.score_cell_cycle(adata, species='mouse')")
        print("2. Provide custom gene lists matching your dataset")
        print("3. Check your gene symbols for correct capitalization or naming conventions")
        
        if len(s_genes_in_data) < 3 or len(g2m_genes_in_data) < 3:
            raise ValueError("Insufficient cell cycle genes found. Cannot proceed with scoring.")
    
    print(f"Scoring cell cycle using {len(s_genes_in_data)} S-phase genes and {len(g2m_genes_in_data)} G2M-phase genes")
    if len(s_genes_in_data) < 10:
        print(f"S-phase genes used: {', '.join(s_genes_in_data)}")
    if len(g2m_genes_in_data) < 10:
        print(f"G2M-phase genes used: {', '.join(g2m_genes_in_data)}")
    
    # Create a copy if requested
    if copy:
        adata = adata.copy()
    
    # Score cell cycle
    sc.tl.score_genes_cell_cycle(
        adata, 
        s_genes=s_genes_in_data, 
        g2m_genes=g2m_genes_in_data
    )
    
    # Add information about species used
    adata.uns["cell_cycle"] = {
        "species": species,
        "s_genes_used": s_genes_in_data,
        "g2m_genes_used": g2m_genes_in_data
    }
    
    # Plot cell cycle scores
    if plot:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # Scatter plot
        sc.pl.scatter(
            adata, 
            x='S_score', 
            y='G2M_score', 
            color='phase',
            title=f'Cell Cycle Scores ({species.capitalize()})',
            show=False,
            ax=axes[0]
        )
        
        # Distribution by phase
        phase_counts = adata.obs['phase'].value_counts()
        axes[1].bar(phase_counts.index, phase_counts.values, color=['#1f77b4', '#ff7f0e', '#2ca02c'])
        for i, v in enumerate(phase_counts.values):
            axes[1].text(i, v + 0.5, str(v), ha='center')
        axes[1].set_title('Cell Cycle Phase Distribution')
        axes[1].set_ylabel('Number of Cells')
        
        plt.tight_layout()
        plt.show()
        
        if save_dir:
            import os
            os.makedirs(save_dir, exist_ok=True)
            plt.savefig(os.path.join(save_dir, f"cell_cycle_scores_{species}.png"), dpi=300)
            plt.close()
    
    # Regress out cell cycle effects if requested
    if regress_out:
        print("Regressing out cell cycle effects...")
        from .normalize import regress_out as reg_out
        reg_out(adata, keys=['S_score', 'G2M_score'], layer=layer, output_layer="cell_cycle_regressed")
        print("Cell cycle effects regressed out and stored in 'cell_cycle_regressed' layer.")
    
    return adata