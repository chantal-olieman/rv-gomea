import os
import numpy

linkage_tree_time = []
linkage_tree_obj = []
dependency_matrix_time = []
dependency_matrix_obj = []
runs = 20


for i in range(runs):
    results = os.popen('./RV-GOMEA -b -s -r -f -2 7 100 -115 -100 0 0.35 10 25 0.9 1 0 1e-10 100 0.0 30').readlines()
    try:
        linkage_tree_obj.append(float(results[0].split(" ")[3]))
        linkage_tree_time.append(float(results[0].split(" ")[5]))
    except IndexError:
        print(f"error: {results}")

for i in range(runs):
    results = os.popen('./RV-GOMEA -b -s -r -f -3 7 100 -115 -100 0 0.35 10 25 0.9 1 0 1e-10 100 0.0 30').readlines()
    try:
        dependency_matrix_obj.append(float(results[0].split(" ")[3]))
        dependency_matrix_time.append(float(results[0].split(" ")[5]))
    except IndexError:
        print(f"error: {results}")


print(f"the linkage tree runs for {numpy.mean(linkage_tree_time)} "
      f"on average with an objective of {numpy.mean(linkage_tree_obj)}")
print(f"the dependency tree runs for {numpy.mean(dependency_matrix_time)} "
      f"on average with an objective of {numpy.mean(dependency_matrix_obj)}")

result = "time-linkage:"+str(linkage_tree_time) +"\ntime-dependency:"+ str(dependency_matrix_time) +"\nobjective-linkage:"+ str(linkage_tree_obj) +"\nobjective-dependency:"+  str(dependency_matrix_obj)
file = open("output.txt", "w")
file.write(str(result))
print("done")