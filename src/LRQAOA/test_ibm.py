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

import numpy as np
import pandas as pd
import json
import os
from csv import writer
from pathlib import Path
import shutil

from qiskit import qasm3, qpy

from qiskit_aer import AerSimulator

from qiskit_aer.primitives import SamplerV2 as Sampler

from qiskit_ibm_runtime import QiskitRuntimeService, Session
from qiskit import QuantumCircuit

from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit.quantum_info import SparsePauliOp

from LRQAOA import *
from util import *


def genParameterList():
    """
    returns list of parameters for experiments
    """
    params = []

    number_shots = [1000]
    folder_root = "../QUBOS/"
    rep = 1

    # For Simulated LRQAOA: False; For LRQAOA on IBM hardware: True
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
        draw["n_toolkits"] = [3, 9, 13, 16, 18, 19]
        draw["n_press_machines"] = 2
        draw["penalties_1"] = [3, 4, 5]
        draw["penalties_2"] = [7, 8, 9]
        draw["qubo_type"] = "raw"
        params.append(draw)

        drounded = {}
        drounded["quantum"] = quantum
        drounded["res_dir_path"] = res_dir_path + "rounded/"
        drounded["number_shots"] = number_shots
        drounded["folder_root"] = folder_root
        drounded["rep"] = rep
        drounded["n_toolkits"] = [3, 9, 13, 16, 18, 19]
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
        dscaled["n_toolkits"] = [3, 9, 13, 16, 18, 19]
        dscaled["n_press_machines"] = 2
        dscaled["penalties_1"] = [0]
        dscaled["penalties_2"] = [1, 10]
        dscaled["qubo_type"] = "scaled"
        params.append(dscaled)

    return params


params = genParameterList()

service_type = "free"
dbetaL = [0.6]
dgammaL = [0.9]
p = [1, 2, 5, 10]
opt_levels = [0, 3]

error_mitigation = False

instance = {
    "free": "ibm-q/open/main",
}
with open("credentials.json", "r") as file:
    tokens = json.load(file)

List = []
List_rows = [
    "JobId",
    "Toolkits",
    "PressMachines",
    "Penalty1",
    "Penalty2",
    "NumQubits",
    "NumShots",
    "p",
    "dbeta",
    "dgamma",
    "ServiceType",
    "QiskitOptLevel",
    "QuboName",
    "BackendName",
    "SessionID",
    "PrimitiveID",
    "Version",
    "CreationDate",
    "ErrorMitigation",
    "LCNumQubits",
    "LCOperations",
    "LCNumInstructions",
    "LCDepth",
    "LCNumNonLocalGates",
    "TCNumQubits",
    "TCOperations",
    "TCNumInstructions",
    "TCDepth",
    "TCNumNonLocalGates",
]
i = 0
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
    for _rep in range(rep):
        for opt_level in opt_levels:
            for k in range(len(n_toolkits)):
                for penalty1 in penalties_1:
                    for penalty2 in penalties_2:
                        for n_shots in number_shots:
                            file_base = (
                                "ibm_"
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
                                + str(n_shots)
                                + "_shots_"
                                + str(rep)
                            )
                            file_name = file_base + ".json"

                            # Load qubo matrix
                            if qubo_type == "raw":
                                qubo_name = (
                                    qubo_type
                                    + "_qubo_capacity_penalty_1e+0"
                                    + str(penalty1)
                                    + "_assignment_penalty_1e+0"
                                    + str(penalty2)
                                    + "_"
                                    + str(n_press_machines)
                                    + "-machines_"
                                    + str(n_toolkits[k])
                                    + "-diesets.npy"
                                )
                            elif qubo_type == "scaled":
                                qubo_name = (
                                    qubo_type
                                    + "_qubo_ratio_capacity_to_assignment_"
                                    + str(penalty2)
                                    + "_"
                                    + str(n_press_machines)
                                    + "-machines_"
                                    + str(n_toolkits[k])
                                    + "-diesets.npy"
                                )
                            else:  # rounded
                                qubo_name = (
                                    qubo_type
                                    + "_qubo_"
                                    + str(n_press_machines)
                                    + "-machines_"
                                    + str(n_toolkits[k])
                                    + "-diesets.npy"
                                )

                            qubofile = folder_root + qubo_name

                            if not os.path.isfile(qubofile):
                                print("Could not open qubo file")
                            qubo = np.load(qubofile)

                            num_qubits = len(qubo)

                            for dbeta in dbetaL:
                                for dgamma in dgammaL:
                                    if quantum:
                                        service = QiskitRuntimeService(
                                            channel="ibm_quantum",
                                            instance=instance[service_type],
                                            token=tokens[service_type],
                                        )
                                        if service_type == "free":
                                            backend = service.backend(name="ibm_kyiv")
                                        elif service_type == "paid":
                                            backend = service.backend(name="ibm_marrakesh")

                                        sampler = Sampler(mode=backend)
                                        sampler.options.default_shots = n_shots

                                        if error_mitigation:
                                            sampler.options.dynamical_decoupling.enable = (
                                                True
                                            )
                                            sampler.options.dynamical_decoupling.sequence_type = (
                                                "XY4"
                                            )
                                            sampler.options.twirling.enable_gates = True
                                            sampler.options.twirling.num_randomizations = (
                                                "auto"
                                            )

                                        for p_ in p:
                                            print(i)
                                            i += 1
                                            qc = create_LR_QAOA_Circuit(
                                                qubo, dbeta, dgamma, p_
                                            )
                                            qc.measure_all()

                                            pm = generate_preset_pass_manager(
                                                backend=backend,
                                                optimization_level=opt_level,
                                            )
                                            ic = pm.run(qc)

                                            job = sampler.run([ic])
                                            print(f"job id: {job.job_id()}")

                                            file_baseC = (
                                                res_dir_path
                                                + "job-"
                                                + job.job_id()
                                                + "/"
                                                + file_base
                                                + str(p_)
                                                + "_"
                                                + str(dbeta)
                                                + "_"
                                                + str(dgamma)
                                            )

                                            Path(
                                                res_dir_path + "job-" + job.job_id()
                                            ).mkdir(parents=True, exist_ok=True)

                                            with open(
                                                file_baseC + "Qasm_logical_", "w"
                                            ) as f:
                                                qasm3.dump(qc, f)
                                            with open(
                                                file_baseC + "Qpy_logical_", "wb"
                                            ) as fb:
                                                qpy.dump(qc, fb)
                                            with open(
                                                file_baseC + "Qasm_transpiled_", "w"
                                            ) as f:
                                                qasm3.dump(ic, f)
                                            with open(
                                                file_baseC + "Qpy_transpiled_", "wb"
                                            ) as fb:
                                                qpy.dump(ic, fb)

                                            job_folder = (
                                                res_dir_path
                                                + "job-"
                                                + job.job_id()
                                                + "/"
                                            )
                                            shutil.copy(
                                                folder_root + qubo_name,
                                                job_folder + qubo_name,
                                            )
                                            shutil.copy(
                                                "test_ibm.py",
                                                job_folder + "test_ibm.py",
                                            )
                                            shutil.copy(
                                                "LRQAOA.py", job_folder + "LRQAOA.py"
                                            )
                                            shutil.copy(
                                                "util.py", job_folder + "util.py"
                                            )

                                            List.append(
                                                [
                                                    job.job_id(),
                                                    n_toolkits[k],
                                                    n_press_machines,
                                                    penalty1,
                                                    penalty2,
                                                    num_qubits,
                                                    n_shots,
                                                    p_,
                                                    dbeta,
                                                    dgamma,
                                                    service_type,
                                                    opt_level,
                                                    qubo_name,
                                                    job.backend().name,
                                                    job.session_id,
                                                    job.primitive_id,
                                                    job._version,
                                                    str(job.creation_date),
                                                    str(error_mitigation),
                                                    qc.num_qubits,
                                                    str(qc.count_ops()),
                                                    qc.size(),
                                                    qc.depth(),
                                                    qc.num_nonlocal_gates(),
                                                    ic.num_qubits,
                                                    str(ic.count_ops()),
                                                    ic.size(),
                                                    ic.depth(),
                                                    ic.num_nonlocal_gates(),
                                                ]
                                            )

                                            with open(
                                                file_baseC + "_experiment.csv", "w"
                                            ) as f:
                                                writer_obj = writer(f)
                                                writer_obj.writerow(List_rows)
                                                writer_obj.writerow(List[-1])
                                                f.close()
                                    else:
                                        for p_ in [1,2,5,10,20,50,100]:
                                            if n_toolkits[k] <= 3:
                                                i+=1
                                                print(i)
                                    
                                                qc = create_LR_QAOA_Circuit(
                                                    qubo, dbeta, dgamma, p_
                                                )
                                                qc.measure_all()

                                                sampler = Sampler()
                                                ressamp = sampler.run([qc], shots=n_shots).result()[0].data.meas.get_counts()

                                                dfDict = dict()
                                                dfDict["Bitstr"] = []
                                                dfDict["Counts"] = []
                                                dfDict["Energy"] = []
                                                dfDict["p"] = p_
                                                dfDict["Toolkits"] = n_toolkits[k]
                                                dfDict["PressMachines"] = n_press_machines
                                                dfDict["Penalty1"] = penalty1
                                                dfDict["Penalty2"] = penalty2
                                                dfDict["NumShots"] = n_shots
                                                dfDict["dbeta"] = dbeta
                                                dfDict["dgamma"] = dgamma
                                                dfDict["QuboName"] = qubo_name
                                                dfDict["Variant"] = qubo_type
                                                dfDict["LCNumQubits"] = qc.num_qubits
                                                dfDict["LCOperations"] = str(qc.count_ops())
                                                dfDict["LCNumInstructions"] = qc.size()
                                                dfDict["LCDepth"] = qc.depth()
                                                dfDict["LCNumNonLocalGates"] = qc.num_nonlocal_gates()
                                                for bitstr in ressamp.keys():
                                                    bitstr2 = [int(x) for x in invertBitStr(bitstr[::-1])]
                                                    dfDict["Energy"].append(np.array(bitstr2).T @ qubo @ np.array(bitstr2))
                                                    dfDict["Bitstr"].append(invertBitStr(bitstr[::-1]))
                                                    dfDict["Counts"].append(ressamp[bitstr])
                                                
                                                df = pd.DataFrame(dfDict, columns=["Bitstr", "Counts", "Energy", "p", "Toolkits", "PressMachines", "Penalty1", "Penalty2","NumShots","dbeta","dgamma","QuboName", "Variant", "LCNumQubits","LCOperations","LCNumInstructions","LCDepth","LCNumNonLocalGates"])
                                                
                                                Path(
                                                res_dir_path + "/" + str(i)
                                                    ).mkdir(parents=True, exist_ok=True)
                                                df.to_csv(res_dir_path + "/" + str(i) + "/" + "res.csv")

exists = Path(res_dir_path + "experiments.csv").exists()
with open(res_dir_path + "experiments.csv", "a", newline="\n") as f_object:
    writer_object = writer(f_object)
    if not exists:
        writer_object.writerow(List_rows)
    for L in List:
        writer_object.writerow(L)
    f_object.close()
