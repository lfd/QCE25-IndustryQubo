###
Create and run the docker file:
```
docker build --no-cache -t ext .

docker run -it ext
```

# Tree structure
## /Additional_Figures
- **/\<k>_Toolkits**: Plots for energy and valid solutions (with and without mitigation) for LR-QAOA, Quantum and Simulated Annealing
## /src
- **/DWave**: Code for hardware experiments on D-Wave and Simulated Annealing
- **/LRQAOA**: Code for hardware experiments on IBM via LRQAOA and simulated LRQAOA
- **/img-gen**: Figures that are shown in the paper
- **/QUBOS**: QUBO formulations for toolkits and penalty weights
- **/results**: .csv files that contain the experimental results

Create Plots:

```
cd src
make
```

Run (Simulated) Annealing:
```
cd src/DWave/
python testing_dwave.py
```

Run (Simulated) LRQAOA:
```
cd src/LRQAOA/
python test_ibm.py
```