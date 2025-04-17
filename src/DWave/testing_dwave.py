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

from DWave_params import *
import numpy as np
import dimod
from dwave.system import EmbeddingComposite, DWaveSampler
from dwave.samplers import SimulatedAnnealingSampler
import json
import os
import pandas as pd
import shutil
import dill

# dwave-ocean-sdk==7.1.0

def genParameterList():
    """
        returns list of parameters for experiments
    """
    params = []

    number_shots = [500]
    folder_root = "../QUBOS/"
    rep = 1 
    
    # For Simulated Annealing: False; For Quantum Annealing: True
    for quantum in [False]:
        if quantum:
            res_dir_path = "results/quantum/"
        else:
            res_dir_path = "results/simulated/"

        draw = {}
        draw["quantum"] = quantum
        draw["res_dir_path"] = res_dir_path + "raw/"
        draw["number_shots"] = number_shots
        draw["folder_root"] = folder_root
        draw["rep"] = rep
        draw["n_toolkits"] = [3,9,13,16,18,19]
        draw["n_press_machines"] = 2
        draw["penalties_1"] = [3,4,5]
        draw["penalties_2"] = [7,8,9]
        draw["qubo_type"] = "raw"
        params.append(draw)

        drounded = {}
        drounded["quantum"] = quantum
        drounded["res_dir_path"] = res_dir_path + "rounded/"
        drounded["number_shots"] = number_shots
        drounded["folder_root"] = folder_root
        drounded["rep"] = rep
        drounded["n_toolkits"] = [3,9,13,16,18,19]
        drounded["n_press_machines"] = 2
        drounded["penalties_1"] = [0]
        drounded["penalties_2"] = [0]
        drounded["qubo_type"] = "rounded"
        params.append(drounded)

        dscaled = {}
        dscaled["quantum"] = quantum
        dscaled["res_dir_path"] = res_dir_path + "scaled/"
        dscaled["number_shots"] = number_shots
        dscaled["folder_root"] = folder_root
        dscaled["rep"] = rep
        dscaled["n_toolkits"] = [3,9,13,16,18,19]
        dscaled["n_press_machines"] = 2
        dscaled["penalties_1"] = [0]
        dscaled["penalties_2"] = [1,10]
        dscaled["qubo_type"] = "scaled"
        params.append(dscaled)

    return params
    
params = genParameterList()
annealing_time = [10, 20, 40, 80, 160, 320, 640, 1280]
auto_scale = [True]
flux_drift_compensation = [True]
programming_thermalization = [1000.0]
readout_thermalization = [100.0]
reduce_intersample_correlation = [False]

i = 0
global_params = []

dwaveParamsDict = {}
for dic in params:
    qubo_type = dic["qubo_type"]
    rep = dic["rep"]
    n_toolkits = dic["n_toolkits"]
    penalties_1 = dic["penalties_1"]
    penalties_2 = dic["penalties_2"]
    n_press_machines = dic["n_press_machines"]
    folder_root = dic["folder_root"]
    number_shots = dic["number_shots"]
    res_dir_path = dic["res_dir_path"]
    quantum = dic["quantum"]

    if quantum:
        anneal_Schedule = [
            "hc_bowunder",
            "hc_bowover",
            "hc_steepflatsteep",
            "linear"
        ]
    else:
        anneal_Schedule = ["linear"]

    for _rep in range(rep):
        for k in range(len(n_toolkits)):
            for penalty1 in penalties_1:
                for penalty2 in penalties_2:
                    if qubo_type == "raw":
                        qubo_name = qubo_type + "_qubo_capacity_penalty_1e+0" + str(penalty1) + "_assignment_penalty_1e+0" + str(penalty2) + "_" + str(n_press_machines) + "-machines_" + str(n_toolkits[k]) + "-diesets.npy"
                    elif qubo_type == "scaled":
                        qubo_name = qubo_type + "_qubo_ratio_capacity_to_assignment_" + str(penalty2) + "_" + str(n_press_machines) + "-machines_" + str(n_toolkits[k]) + "-diesets.npy"
                    else: #rounded
                        qubo_name = qubo_type + "_qubo_" + str(n_press_machines) + "-machines_" + str(n_toolkits[k]) + "-diesets.npy"

                    qubofile = folder_root + qubo_name

                    if not os.path.isfile(qubofile):
                        print("Could not open qubo file")
                    qubo = np.load(qubofile)

                    Q = dimod.BinaryQuadraticModel(qubo, "BINARY").to_qubo()[0]
                    for n_shots in number_shots:
                        for asched in anneal_Schedule:
                            for atime in annealing_time:
                                for autsc in auto_scale:
                                    for fdriftcomp in flux_drift_compensation:
                                        for progtherm in programming_thermalization:
                                            for readtherm in readout_thermalization:
                                                for (
                                                    interscorr
                                                ) in reduce_intersample_correlation:
                                                    i += 1
                                                    print(i)
                                                    results_dir_path = (
                                                        res_dir_path + str(i) + "/"
                                                    )
                                                    if not os.path.exists(results_dir_path):
                                                        os.makedirs(results_dir_path)

                                                    hgsched = [
                                                        [
                                                            [0, 1],
                                                            [atime, 1],
                                                        ]
                                                    ]

                                                    local_params = pd.DataFrame(
                                                        {
                                                            "JobId": i,
                                                            "_Rep": _rep,
                                                            "Rep": rep,
                                                            "Schedule": asched,
                                                            "Annealing_time": atime,
                                                            "Auto_scale": autsc,
                                                            "Flux_drift_compensation": fdriftcomp,
                                                            "H_gain_schedule": str(hgsched),
                                                            "Programming_Thermalisation": progtherm,
                                                            "Readout_thermalisation": readtherm,
                                                            "Reduce_intersample_correlation": interscorr,
                                                            "Toolkits": n_toolkits[k],
                                                            "PressMachines": n_press_machines,
                                                            "NumQubits": len(qubo),
                                                            "Penalty1": penalty1,
                                                            "Penalty2": penalty2,
                                                            "NumShots": number_shots,
                                                            "QuboName": qubo_name,
                                                            "Quantum": str(quantum),
                                                            "Variant": qubo_type,
                                                        },
                                                        index = [0]
                                                    )
                                                    if i >= 0:                                                    
                                                        local_params.to_csv(
                                                            results_dir_path
                                                            + "local_parameters.csv"
                                                        )

                                                        file_name = (
                                                            "dw_"
                                                            + str(quantum)
                                                            + "_results_"
                                                            + "QUBO_"
                                                            + str(n_press_machines)
                                                            + "_machines_"
                                                            + str(n_toolkits[k])
                                                            + "_toolkits_1e"
                                                            + str(penalty1)
                                                            + "_1e"
                                                            + str(penalty2)
                                                            + "_penalty_"
                                                            + str(len(qubo))
                                                            + "_qubits_"
                                                            + str(n_shots)
                                                            + "_shots_"
                                                            + str(rep)
                                                            + ".json"
                                                        )

                                                        if quantum:
                                                            sampler = EmbeddingComposite(
                                                                DWaveSampler()
                                                            )

                                                            print("Submitting...")
                                                            sampleset = sampler.sample_qubo(
                                                                Q,
                                                                answer_mode="raw",
                                                                num_reads=n_shots,
                                                                anneal_schedule=scale_schedule(
                                                                    get_schedule(name=asched),
                                                                    atime,
                                                                ),
                                                                auto_scale=autsc,
                                                                flux_drift_compensation=fdriftcomp,
                                                                programming_thermalization=progtherm,
                                                                readout_thermalization=readtherm,
                                                                reduce_intersample_correlation=interscorr,
                                                                return_embedding = True,
                                                            )

                                                            with open(results_dir_path + "SamplerProperties.json", "w") as json_file:
                                                                json.dump(
                                                                    sampler.properties,
                                                                    json_file,
                                                                    ensure_ascii=False,
                                                                    indent=2,
                                                                )
                                                            with open(results_dir_path + "SamplerParameters.json", "w") as json_file:
                                                                json.dump(
                                                                    sampler.parameters,
                                                                    json_file,
                                                                    ensure_ascii=False,
                                                                    indent=2,
                                                                )

                                                        else:
                                                            sampler = SimulatedAnnealingSampler()
                                                            
                                                            sampleset = sampler.sample_qubo(
                                                                Q,
                                                                seed=1234,
                                                                answer_mode="raw",
                                                                num_reads=n_shots,
                                                                num_sweeps= atime
                                                            )

                                                        pdsampleset = sampleset.to_pandas_dataframe()
                                                        pdsampleset["JobId"] = i
                                                        if quantum:
                                                            embedding = sampleset.to_serializable()["info"]["embedding_context"]["embedding"]
                                                            embedding_size = 0
                                                            for var in embedding.keys():
                                                                embedding_size += len(embedding[var])
                                                            pdsampleset["EmbeddingSize"] = embedding_size
                                                            pdsampleset["Embedding"] = str(embedding)
                                                        
                                                        pdsampleset.to_csv(results_dir_path + "sampleset.csv")

                                                        shutil.copy("testing_dwave.py", results_dir_path + "testing_dwave.py")
                                                        shutil.copy("DWave_params.py", results_dir_path + "DWave_params.py")
                                                        shutil.copy(qubofile, results_dir_path + qubo_name)

                                                        try:
                                                            with open(results_dir_path + "sampleset.pkl", 'wb') as fb:
                                                                dill.dump(sampleset.to_serializable(), fb)
                                                        except Exception as error:
                                                            print("Dill error: " + type(error).__name__, error)

                                                        file_path = os.path.join(results_dir_path, file_name)
                                                        with open(file_path, "w") as json_file:
                                                            json.dump(
                                                                sampleset.to_serializable(),
                                                                json_file,
                                                                ensure_ascii=False,
                                                                indent=2,
                                                            )
