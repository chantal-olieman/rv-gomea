import os
import numpy
import json

runs = 10

linkage_options = [1, -8]

result_times = [[] for i  in range(len(linkage_options)+1)]
result_evals = [[] for i  in range(len(linkage_options)+1)]
circles = [4 ,9, 16, 25, 36, 49]
vtr_file = open(f"{os.getcwd()}/../plot-cias/optima.txt")
vtr = json.load(vtr_file)
blackbox = ""
rotation = 0
problem = 16

for circle in circles:
    times = []
    evals = []
    for linkage in linkage_options:
        learned_linkage = f"./RV-GOMEA -s -r {blackbox} -f {linkage} {problem} {circle*2} 0 1 {rotation} " \
                          f"0.35 100 1 0.9 1 0 {1.0001*float(vtr[str(circle)])} 100 0.0 300"
        print(f"settings: {learned_linkage}")
        linkage_evals = []
        linkage_time = []
        max_time = 0
        for i in range(runs):
            results = os.popen(learned_linkage).readlines()
            if i % 5 == 0:
                print(f"run {i}")
            try:
                linkage_evals.append(float(results[0].split(" ")[1]))
                linkage_time.append(float(results[0].split(" ")[5]))
                if linkage_time[-1] > 600:
                    max_time += 1
                    if max_time >= 3:
                        break
            except Exception as e:
                print(f"error: {results}")
        evaluations = numpy.median(linkage_evals)
        time = numpy.median(linkage_time)
        times.append(time)
        evals.append(evaluations)
    result_times[0].append(circle)
    result_evals[0].append(circle)
    for i in range(len(times)):
        result_times[i+1].append(times[i])
    for i in range(len(evals)):
        result_evals[i+1].append(evals[i])
    print(result_evals)
    print(result_times)




print("done")