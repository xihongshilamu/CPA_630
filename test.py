# some standard packages to assist this tutorial
import cpa
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

 
def test(data_path, save_path):
    adata = sc.read(data_path) # 14811*4999
    cpa_api = cpa.api.API(
        adata, 
        pretrained=save_path
    )
    cpa_api.compute_comb_emb(thrh=0)
    uncer = cpa_api.compute_uncertainty(
                        cov={'cell_type': 'A549'}, 
                        pert='Nutlin', 
                        dose='1.0'
                    )
    print(uncer)

def train(adata_path, model_save_dir):
    adata = sc.read(adata_path)
    # adata.obs = adata.obs.astype('category')
    # adata.obs['split'] = adata.obs['split'].astype('category')
    # adata.obs['knock_gene'] = adata.obs['knock_gene'].astype('category')
    # adata.obs['knock_or_not'] = adata.obs['knock_or_not'].astype('category')
    # adata.obs['control'] = adata.obs['control'].astype('category')

    adata.obs['knock_gene'].cat.categories
    adata.obs['knock_or_not'].cat.categories


    cpa_api = cpa.api.API(
        adata, 
        perturbation_key='knock_gene',
        doser_type='logisim',
        dose_key='batch',
        split_key='split',
        covariate_keys=['knock_or_not'],
        save_dir=model_save_dir,
        only_parameters=False,
        hparams={}, 
    )  
    
    cpa_api.train(
    max_epochs=1000, 
    run_eval=True, 
    checkpoint_freq=20,
    filename='bulk_cpa.pt', 
    max_minutes=2*60
)
    
def main():
    data_path = '/data/share/cnic02/projects/perturb_cwt/bulk_3072_train_test_CPA.h5ad'
    model_save_dir = '/data/share/cnic02/models/cpa_models/0424_bulk/'
    model_save_path ='pretrained_models/gsm.pt'
    train(data_path, model_save_dir)
    
    
if __name__ == '__main__':
    main()