import numpy as np
import os
import json

results = {}
sizes = [5 * 2**x for x in range(1,7)]

for size in sizes:
    max_list = []
    sim_list = []

    for i in range(20):
        gomea_command = f" ./RV-GOMEA -b -f -8 -s -r 19 {size} -115 -100 45 0.35 10 25 0.9 1 0 1e-1 100 0.0 1"

        gomea_result = os.popen(gomea_command).readlines()
        max_list.append(float(gomea_result[-3].split(" ")[2]))
        sim_list.append(float(gomea_result[-2].split(" ")[4]))
        print(f"size: {size}, run: {i}")

    results[size] = {"max":np.mean(max_list), "similarity": np.mean(sim_list)}

print(results)
file = open(f"{os.getcwd()}/data-8.txt", "w")
json.dump(results, file)

similarities = []
for key, value in results.items():
    similarities.append(value["similarity"])

import pylab
pylab.loglog(sizes, similarities, marker='o')
pylab.savefig(f"plots_sim-8")
pylab.show()