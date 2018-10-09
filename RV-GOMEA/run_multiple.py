import os
import numpy

result_times = [[],[],[],[]]
result_evals = [[],[],[],[]]
runs = 10

linkage_options = [1, -7]
populations = [1215]
blackbox = ""
rotation = 0
problem = 0

for population in populations:
    times = []
    evals = []
    for linkage in linkage_options:
        learned_linkage = f"./RV-GOMEA -s -r {blackbox} -f {linkage} {problem} {population} -115 -100 {rotation} " \
                          f"0.35 10 25 0.9 1 0 1e-10 100 0.0 300"
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
    result_times[0].append(population)
    result_evals[0].append(population)
    for i in range(len(times)):
        result_times[i+1].append(times[i])
    for i in range(len(evals)):
        result_evals[i+1].append(evals[i])
    file = open("output_time_sphere.txt", "w")
    file.write(str(result_times))
    file = open("output_evals_sphere.txt", "w")
    file.write(str(result_evals))



file = open("output_time_sphere.txt", "w")
file.write(str(result_times))
file = open("output_evals_sphere.txt", "w")
file.write(str(result_evals))
print("done")