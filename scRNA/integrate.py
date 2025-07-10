import anndata as ad
import numpy as np
import scanorama
import scvi
import scanpy.external as sce

def integration_scanorama(adata, batch_key= "study", dims=50):
    
    var_select = adata.var.highly_variable_nbatches > 2
    var_genes = var_select.index[var_select]
    print("Genes that are variable in at least 2 datasets and use for remaining analysis: ", len(var_genes))
    
    # split per batch into new objects.
    batches = adata.obs[batch_key].cat.categories.tolist()
    alldata = {}
    for batch in batches:
        alldata[batch] = adata[adata.obs[batch_key] == batch,]

    #subset the individual dataset to the variable genes we defined at the beginning
    alldata2 = dict()
    for ds in alldata.keys():
        print(ds)
        alldata2[ds] = alldata[ds][:,var_genes]

    #convert to list of AnnData objects
    adatas = list(alldata2.values())

    #run scanorama.integrate
    scanorama.integrate_scanpy(adatas, dimred = dims) 
    
     #scanorama adds the corrected matrix to adata.obsm in each of the datasets in adatas.
    adatas[0].obsm['X_scanorama'].shape
    
    # Get all the integrated matrices.
    scanorama_int = [ad_.obsm['X_scanorama'] for ad_ in adatas]

    # make into one matrix.
    all_s = np.concatenate(scanorama_int)
    print(all_s.shape)

    # add to the AnnData object
    adata.obsm["X_Scanorama"] = all_s
    
    return adata


def integration_scvi(
    adata, 
    layer="counts", 
    batch_key="sampleID", 
    batch_size=256,
    max_epochs=500,
    n_layers=2, 
    n_latent=30, 
    #n_hidden=128,
    ):
    scvi.model.SCVI.setup_anndata(adata, layer=layer, batch_key=batch_key)
    
    model = scvi.model.SCVI(
        adata, 
        n_layers=n_layers, 
        n_latent=n_latent,
        #n_hideen=n_hidden, 
        gene_likelihood="nb")
    
    model.train(
        batch_size=batch_size,
        max_epochs=max_epochs,
        early_stopping=True,
    )
    
    adata.obsm["X_scVI"] = model.get_latent_representation()
 
    return adata


def integration_harmony(adata, batch_key = "sampleID", basis='X_pca'):
    sce.pp.harmony_integrate(adata, 
                             key = batch_key,
                             basis = basis)
    return adata