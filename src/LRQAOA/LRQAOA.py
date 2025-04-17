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

from qiskit import QuantumCircuit
import numpy as np
import networkx as nx
import warnings
from edge_coloring import *

from util import *

def create_LR_QAOA_Circuit(qubo, dbeta, dgamma, p):
    """
    IMPORTANT NOTICE: The optimal bitstring will be inverted (1->0, 0->1) and flipped (::-1) (e.g. if x = 10010 is optimal in the qubo, s = 10110 is the optimal solution for the Ising Model)
    Creates a depth optimised LR-QAOA circuit with p layers based on the symmetric qubo provided as a matrix.
    This circuit needs to undergo a hardware transpile step to solve the qubit assignment and qubit routing problem.
    Input:
        qubo: numpy 2D Array
        dbeta: float
        dgamma: float
        p: int
    """
    qubo = np.array(qubo, dtype=np.float64)
    #Do a symmetry check on the incoming qubo
    for r in range(len(qubo)):
        for c in range(r, len(qubo[r])):
            if qubo[r][c] != qubo[c][r]:
                warnings.warn("Input Qubo not symmetric! Continuing operation...")

    #Make qubo symmetric, to allow for arbitrary matrix inputs
    for r in range(len(qubo)):
        for c in range(r+1, len(qubo[r])):
            val = (qubo[r][c] + qubo[c][r]) / 2.0
            qubo[r][c] = val
            qubo[c][r] = val

    # Create Ising model
    J = np.zeros(shape=(len(qubo), len(qubo)))
    h = np.zeros(shape=(len(qubo)))
    O = 0

    for r in range(len(qubo)):
        sum = 0
        for c in range(len(qubo[r])):
            # Note that the diagonal is added to the constant term
            if r != c:
                J[r][c] = qubo[r][c] / 4.0
            sum += qubo[r][c]
            O += qubo[r][c] / 4.0
        h[r] = sum / 2.0

    # Add the diagonal of qubo to the constant term
    for r in range(len(qubo)):
        O += qubo[r][r] / 4.0

    # Make J upper triangular
    for r in range(len(J)):
        for c in range(r + 1, len(J[r])):
            J[r][c] = J[r][c] + J[c][r]
            J[c][r] = 0

    # evaled = evalIsing(J, h, O)
    # print("Values: ", evaled)
    # print("Minimum: ", getMinFromEval(evaled))

    # Compute maximum values
    maxJ = np.max(np.abs(J))
    maxh = np.max(np.abs(h))

    # print("Pre Normalization J:\n", J)
    # print("Pre Normalization h:\n", h)
    # print("Pre Normalization O:\n", O)

    # From now on we use the ising model
    # Normalize J,h,O (Be aware of the Ising transformation factors we applied before)
    Normalizationfactor = np.max([maxJ, maxh])

    #print("Normalization factor: ", Normalizationfactor)

    for r in range(len(J)):
        for c in range(len(J[r])):
            J[r][c] = J[r][c] / Normalizationfactor
        h[r] = h[r] / Normalizationfactor
    O = O / Normalizationfactor

    # evaled = evalIsing(J, h, O)
    # print("Values: ", evaled)
    # print("Minimum: ", getMinFromEval(evaled))

    # print("Post Normalization J:\n", J)
    # print("Post Normalization h:\n", h)
    # print("Post Normalization O:\n", O)

    # Compute edge coloring to map J
    G = nx.Graph()
    for r in range(len(J)):
        for c in range(r + 1, len(J[r])):
            if J[r][c] != 0:
                G.add_edge(r, c, weight=J[r][c])
    # {<edge>: <color>, ...}
    colored_edges = edge_coloring(G)

    # Transform into Layer form: {<color1>: [<edges>], }
    color_dict = dict()
    for key in colored_edges.keys():
        # This implicit sorting is important for the later introduction of gates
        if key[0] < key[1]:
            if color_dict.__contains__(colored_edges[key]):
                color_dict[colored_edges[key]].append(key)
            else:
                color_dict[colored_edges[key]] = [key]
    qc = QuantumCircuit(len(J))
    # Add all hadamard gates
    for i in range(len(J)):
        qc.h(i)

    # p of these layers
    for layer in range(p):
        betai = (1 - (layer / p)) * dbeta
        gammai = ((1 + layer) / p) * dgamma

        # Problem Hamiltonian: Simple Mapper
        ##h
        for r in range(len(h)):
            qc.rz(phi=2 * h[r] * gammai, qubit=r)

        #qc.barrier()
        ##J
        # Simple ordering
        # for r in range(len(J)):
        #     for c in range(r + 1, len(J[r])):
        #         qc.rzz(theta=2 * J[r][c] * gammai, qubit1=r, qubit2=c)

        # Edge colored ordering
        for col in color_dict.keys():
            for gate in color_dict[col]:
                r = gate[0]
                c = gate[1]
                qc.rzz(theta=2 * J[r][c] * gammai, qubit1=r, qubit2=c)

        #qc.barrier()
        # Mixer
        for r in range(len(J)):
            qc.rx(theta=-2 * betai, qubit=r)
        #qc.barrier()

    return qc
