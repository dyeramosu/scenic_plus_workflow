# imports
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from pycisTopic.cistopic_class import *

import pickle

# set dirs
outDir = '/home/jupyter/AM_barcoding_analysis_final/edit/PDAC_scenic_plus_output/'
tmpDir = '/home/jupyter/AM_barcoding_analysis_final/edit/PDAC_scenic_plus_scratch/'

# get object
cisTopic_obj = pickle.load(open(os.path.join(outDir, 'cisTopic/cistopic_obj.pkl'), 'rb'))
print(cisTopic_obj)

# run LDA
print("running models now")

models=run_cgs_models_mallet('/home/jupyter/Mallet/bin/mallet',
                    cisTopic_obj,
                    n_topics=list(range(3, 75, 3)),
                    n_cpu=24,
                    n_iter=500,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    tmp_path=tmpDir, #Use SCRATCH if many models or big data set
                    save_path=None)

# save models
pickle.dump(models,
            open(os.path.join(outDir, 'cisTopic/models/PDAC_500_iter_LDA_3_75.pkl'), 'wb'))

print("finished saving models")