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

def evalQubo(qubo, bitstr=None):
    """
    Input:  qubo as matrix
            (optional) bitstr: Bitstring at which to evaluate the given qubo
    Output: dict of {<bitvector>: <value>}
    """
    out = dict()
    
    if bitstr == None:
        for bitstr in binaryRep(len(qubo)):
            val = 0
            for r in range(len(qubo)):
                for c in range(len(qubo)):
                    if bitstr[r] == "1" and bitstr[c] == "1":
                        val += qubo[r][c]
            out[bitstr] = round(val, 1)
    else:
        val = 0
        for r in range(len(qubo)):
            for c in range(len(qubo)):
                if bitstr[r] == "1" and bitstr[c] == "1":
                    val += qubo[r][c]
        out[bitstr] = round(val, 1)

    return out                

def evalIsing(J, h, O):
    """
    Assumes the use of x_i = (1+s_i)*0.5 => s_i = 1 <=> x_i = 1; s_i = -1 <=> x_i = 0
    Input:  J: matrix
            h: array
            O: float
    Output: dict of {<bitvector>: <value>}
    """

    out = dict()
    for bitstr in binaryRep(len(J)):
        val = 0
        for r in range(len(J)):
            for c in range(len(J)):
                if r!=c:
                    if bitstr[r] == bitstr[c]: 
                        val += J[r][c]
                    else:
                        val -= J[r][c]

            if bitstr[r] == "1":
                val += h[r]
            else:
                val -= h[r]
        val += O

        out[bitstr] = round(val, 1)

    return out 


def binaryRep(n):
      if n == 0:
            return ['']
      else:
            return [i + '0' for i in binaryRep(n-1)] + [i + '1' for i in binaryRep(n-1)]

def getMinFromEval(evalDict):
     """
     Input: Full value dict from evalQubo
     Output: Minimum
     """

     key = min(evalDict, key=evalDict.get)
     return {key: evalDict[key]}

def getTestQubo(idx, vt = 1e8):
     """
     returns a test qubo, prints its value distribution and its minimum
     """
     qubos = [
            #3
            [[0,-1, vt],
            [-1, 0,vt],
            [vt,vt,0]],
            #4
            [[0,-1,vt,vt],
            [-1,0,vt,vt],
            [vt,vt,0,0],
            [vt,vt,0,0]],
            #5
            [[0,-1,vt,vt,vt],
             [-1,0,vt,vt,vt],
             [vt,vt,0,0,0],
             [vt,vt,0,0,0],
             [vt,vt,0,0,0]],
            #6
            [[0,-1,vt,vt,vt,vt],
             [-1,0,vt,vt,vt,vt],
             [vt,vt,0,0,0,0],
             [vt,vt,0,0,0,0],
             [vt,vt,0,0,0,0],
             [vt,vt,0,0,0,0]],
            #7
             [[0,-1,vt,vt,vt,vt,vt],
             [-1,0,vt,vt,vt,vt,vt],
             [vt,vt,0,0,0,0,0],
             [vt,vt,0,0,0,0,0],
             [vt,vt,0,0,0,0,0],
             [vt,vt,0,0,0,0,0],
             [vt,vt,0,0,0,0,0]],
            #8
            [[0,-1,vt,vt,vt,vt,vt,vt],
             [-1,0,vt,vt,vt,vt,vt,vt],
             [vt,vt,0,0,0,0,0,0],
             [vt,vt,0,0,0,0,0,0],
             [vt,vt,0,0,0,0,0,0],
             [vt,vt,0,0,0,0,0,0],
             [vt,vt,0,0,0,0,0,0],
             [vt,vt,0,0,0,0,0,0]],
             #9
            [[0,-1,vt,vt,vt,vt,vt,vt,vt],
             [-1,0,vt,vt,vt,vt,vt,vt,vt],
             [vt,vt,0,0,0,0,0,0,0],
             [vt,vt,0,0,0,0,0,0,0],
             [vt,vt,0,0,0,0,0,0,0],
             [vt,vt,0,0,0,0,0,0,0],
             [vt,vt,0,0,0,0,0,0,0],
             [vt,vt,0,0,0,0,0,0,0],
             [vt,vt,0,0,0,0,0,0,0]],
            
            [[0, 1, 2],
             [1, 0, 4],
             [2, 4, 0],],

            [[0,  1,  2,  4],
             [1,  0,  8, 16],
             [2,  8,  0, 32],
             [4, 16, 32,  0],
             ],

            [[0,  1,  2,  4,  8],
             [1,  0, 16, 32, 64],
             [2, 16,  0,128,256],
             [4, 32,128,  0,512],
             [8, 64,256,512,  0],
             ],

             [[0 ,   1,   2,   4,   8,  16],
              [1 ,   0,  32,  64, 128, 256],
              [2 ,  32,   0, 512,1024,2048],
              [4 ,  64, 512,   0,4096,8192],
              [8 , 128,1024,4096,   0,16384],
              [16, 256,2048,8192,16384, 0],],

            np.add(np.array(
            [[0, 1, 2],
             [1, 0, 4],
             [2, 4, 0],]), -2),
            np.add(np.array(
            [[0,  1,  2,  4],
             [1,  0,  8, 16],
             [2,  8,  0, 32],
             [4, 16, 32,  0],
             ]), -2),
            np.add(np.array(
            [[0,  1,  2,  4,  8],
             [1,  0, 16, 32, 64],
             [2, 16,  0,128,256],
             [4, 32,128,  0,512],
             [8, 64,256,512,  0],
             ]), -2),
            np.add(np.array(
            [[0 ,   1,   2,   4,   8,  16],
              [1 ,   0,  32,  64, 128, 256],
              [2 ,  32,   0, 512,1024,2048],
              [4 ,  64, 512,   0,4096,8192],
              [8 , 128,1024,4096,   0,16384],
              [16, 256,2048,8192,16384, 0],]), -2),

          np.ones(shape=(4,4)),
          [[1,0,0,4],
           [0,1,0,0],
           [0,0,1,0],
           [4,0,0,1]],
           
           [[1,0,0,-4],
            [0,1,0,0],
            [0,0,1,0],
            [-4,0,0,1]],

            [[1,-2,0,-4],
            [-2,1,3,0],
            [0,3,1,0],
            [-4,0,0,1]],

            [[1,0,-1,0,0,-3],
             [0,0,0,0,2,0],
             [-1,0,0,1,0,0],
             [0,0,1,0,0,5],
             [0,2,0,0,0,0],
             [-3,0,0,5,0,1],
             ],

            [[0,0,2,0,0,-3,0,0],
             [0,0,0,0,0,0,0,0],
             [2,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0],
             [-3,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0],
             [0,0,0,0,0,0,0,0],
            ]
     ]

     print("Qubo: ", qubos[idx])
     evaled = evalQubo(qubos[idx])
     print("Values: ", evaled)
     print("Minimum: ", getMinFromEval(evaled))

     return qubos[idx]
     

def invertBitStr(bitstr):
    return ''.join('1' if x == '0' else '0' for x in bitstr)

