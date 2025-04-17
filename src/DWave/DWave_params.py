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

import random
import math
import matplotlib.pyplot as plt


def scale_schedule(schedule, time):
    """
    Takes a normalized schedule and scales it according to total annealing time in [us]
    """
    out = []
    for i in range(len(schedule)):
        out.append((schedule[i][0] * time, schedule[i][1]))

    return out

def plot_schedule(schedule, filename):
    """
    Plots a given schedule
    """
    plt.plot(*zip(*schedule))
    plt.savefig(filename)
    plt.show()

def get_schedule(name, param=None):
    """
    return a normalized schedule identified by name
    """
    match name:
        case "hc_bowunder":
            return [
                (0.0, 0.0),
                (0.1, 0.01),
                (0.2, 0.03),
                (0.3, 0.08),
                (0.4, 0.15),
                (0.5, 0.23),
                (0.6, 0.33),
                (0.7, 0.45),
                (0.8, 0.61),
                (0.9, 0.85),
                (1.0, 1.0),
            ]
        case "hc_bowover":
            return [
                (0.0, 0.0),
                (0.1, 0.20),
                (0.2, 0.35),
                (0.3, 0.50),
                (0.4, 0.60),
                (0.5, 0.67),
                (0.6, 0.75),
                (0.7, 0.82),
                (0.8, 0.90),
                (0.9, 0.95),
                (1.0, 1.0),
            ]
        case "hc_steepflatsteep":
            return [
                (0.0, 0.0),
                (0.1, 0.30),
                (0.2, 0.50),
                (0.3, 0.53),
                (0.4, 0.56),
                (0.5, 0.59),
                (0.6, 0.62),
                (0.7, 0.65),
                (0.8, 0.75),
                (0.9, 0.90),
                (1.0, 1.0),
            ]
        case "linear":
            return [(0.0, 0.0), (1.0, 1.0)]