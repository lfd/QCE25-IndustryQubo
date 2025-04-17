# Copyright 2025 OTH - Laboratory for Digitalisation (LfD)
# Written by Lukas Schmidbauer
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pandas as pd
import numpy as np
import pprint as pp
import glob

#Subsumes in each subfolder: local_parameters.csv and sampleset.csv and then joins them

folder_Batches = ["results/simulated/raw/",
                  "results/simulated/rounded/",
                  "results/simulated/scaled/"
                  ]

for folder_Batch in folder_Batches:
    if "quantum" in folder_Batch:
        quantum = True
    else:
        quantum = False
    
    qubo_type = folder_Batch.split("/")[-2]

    filenames_local_parameters = glob.glob(folder_Batch + "*/local_parameters.csv")
    filenames_sampleset = glob.glob(folder_Batch + "*/sampleset.csv")

    lpli = []
    for fn in filenames_local_parameters:
        df = pd.read_csv(fn, index_col=0)
        lpli.append(df)

    dflp = pd.concat(lpli, axis=0, ignore_index=True)
    
    sali = []
    for fn in filenames_sampleset:
        df = pd.read_csv(fn, index_col=0) 

        names = []
        for n in df.columns:
            if n not in ["chain_break_fraction", "energy", "num_occurrences" , "JobId"]:
                names.append(n)

        df["Bitstr"] = df[names].apply(
            lambda x: ''.join(x.dropna().astype(str)),
            axis = 1
        )
        sali.append(df)

    dfsa = pd.concat(sali, axis=0, ignore_index=True)

    dfres = pd.merge(dfsa, dflp, how="left", on="JobId")
    dfres["Quantum"] = quantum

    dfres.to_csv(folder_Batch + "Results_" + qubo_type + ".csv")