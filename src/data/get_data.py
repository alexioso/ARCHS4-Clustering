import os
import requests
import shutil
import h5py
import numpy as np
import pandas as pd
import time



def main(input_filepath, output_filepath):
    """ Runs data processing scripts to turn raw data from (../../data/raw) into
        cleaned data ready to be analyzed (saved in ../../data/processed).
    """
   
    out_file = output_filepath

    if not os.path.isdir("../../data"):
        os.mkdir("../../data")
    if not os.path.isdir("../../data/raw"):
        os.mkdir("../../data/raw")
    if not os.path.isdir("../../data/processed"):
        os.mkdir("../../data/processed")
    if not os.path.isdir("../../data/interim"):
        os.mkdir("../../data/interim")
    if not os.path.isdir("../../data/external"):
        os.mkdir("../../data/external")
    if not os.path.isdir("../../data/processed/tables"):
        os.mkdir("../../data/processed/tables")
    if not os.path.isdir("../../data/processed/figures"):
        os.mkdir("../../data/processed/figures")
    
    # Check if gene expression file was already downloaded, if not in data/raw directory then download file from S3
    if(not os.path.exists(out_file)):
        print("Downloading compressed gene expression matrix.")
        url = input_filepath
        r = requests.get(url, stream=True)
        if r.status_code == 200:
            with open(out_file, 'wb') as f:
                shutil.copyfileobj(r.raw, f)
        del r 
    else:
        print("Local file already exists.")




if __name__ == '__main__':
    
    input_url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix_v9.h5"
    output_filepath = "../../data/raw/human_matrix_v9.h5"
    main(input_url, output_filepath)
