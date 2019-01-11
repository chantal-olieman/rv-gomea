import pylab
import numpy as np
import os
import json

sizes = [6] # [2, 4, 6, 8]

for size in sizes:
    path_file =f"{os.getcwd()}/data-{size}.txt"
    file = open(path_file)
    results = json.load(file)
    similarities = []
    x_values = []
    for key, value in results.items():
        similarities.append(float(value["similarity"]))
        x_values.append(float(key))
    print(x_values)
    print(similarities)
    pylab.loglog(x_values, similarities, marker='o')

pylab.title("Similarity between biggest subsets")
pylab.xlabel("problem size")
pylab.ylabel("similarity")
pylab.savefig(f"plots_sim_all")
pylab.show()

