# imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import pyranges as pr

import pickle
import dill
import os
import sys
import warnings
warnings.filterwarnings("ignore")

_stderr = sys.stderr   
null = open(os.devnull,'wb')

# set directories
outDir = '/home/jupyter/AM_barcoding_analysis_final/edit/PDAC_scenic_plus_output/'
tmpDir = '/home/jupyter/AM_barcoding_analysis_final/edit/PDAC_scenic_plus_scratch/'

# get data
adata = sc.read_h5ad(os.path.join(outDir, 'deepika_modified_gem.h5ad'))
cistopic_obj = dill.load(open(os.path.join(outDir, 'cisTopic/cistopic_obj.pkl'), 'rb'))
menr = dill.load(open(os.path.join(outDir, 'pycisTarget/motifs/menr.pkl'), 'rb'))

# create scenic plus object
print('creating scenic plus object')
from scenicplus.scenicplus_class import create_SCENICPLUS_object

scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata.raw.to_adata(),
    cisTopic_obj = cistopic_obj,
    menr = menr,
    bc_transform_func = lambda x: f'{x}___cisTopic' #function to convert scATAC-seq barcodes to scRNA-seq ones
)
scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())
print('finished creating scenic plus object')

biomart_host = "http://sep2019.archive.ensembl.org/"

scplus_obj.dr_cell['GEX_X_pca'] = scplus_obj.dr_cell['GEX_X_pca'].iloc[:, 0:2]

# run scenic plus
print('running scenic plus')
from scenicplus.wrappers.run_scenicplus import run_scenicplus
try:
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['GEX_states'],
        species = 'hsapiens',
        assembly = 'hg38',
        tf_file = os.path.join(outDir, 'utoronto_human_tfs_v_1.01.txt'),
        save_path = os.path.join(outDir, 'scenicplus'),
        biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = True,
        export_to_loom_file = True,
        export_to_UCSC_file = True,
        path_bedToBigBed = '/home/jupyter/AM_barcoding_analysis_final/edit/PDAC_scenic_plus_output',
        n_cpu = 24,
        _temp_dir = tmpDir)
except Exception as e:
    #in case of failure, still save the object
    dill.dump(scplus_obj, open(os.path.join(outDir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)
    raise(e)
print('finished running')