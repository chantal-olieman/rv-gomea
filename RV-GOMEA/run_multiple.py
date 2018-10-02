    import os
import numpy

linkage_tree_time = []
linkage_tree_obj = []
dependency_matrix_time = []
dependency_matrix_obj = []
runs = 20

population = 10
blackbox = "-b"
rotation = 45
problem = 13

learned_linkage = f"./RV-GOMEA -s -r {blackbox} -f -7 {problem} {population} -115 -100 {rotation} " \
                  f"0.35 10 25 0.9 1 0 1e-10 100 0.0 30"

dependency = f"./RV-GOMEA -s -r {blackbox} -f -6 {problem} {population} -115 -100 {rotation} " \
                  f"0.35 10 25 0.9 1 0 1e-10 100 0.0 30"

print(f"settings: {learned_linkage}")

for i in range(runs):
    results = os.popen(learned_linkage).readlines()
    if i % 5 == 0:
        print(f"run {i}")
    try:
        linkage_tree_obj.append(float(results[0].split(" ")[3]))
        linkage_tree_time.append(float(results[0].split(" ")[5]))
    except IndexError:
        print(f"error: {results}")

for i in range(runs):
    results = os.popen(dependency).readlines()
    if i % 5 == 0:
        print(f"run {i}")
    try:
        dependency_matrix_obj.append(float(results[0].split(" ")[3]))
        dependency_matrix_time.append(float(results[0].split(" ")[5]))
    except IndexError:
        print(f"error: {results}")


print(f"the linkage tree runs for {numpy.mean(linkage_tree_time)} "
      f"on average with an objective of {numpy.mean(linkage_tree_obj)}")
print(f"the dependency tree runs for {numpy.mean(dependency_matrix_time)} "
      f"on average with an objective of {numpy.mean(dependency_matrix_obj)}")

result = "time-linkage:"+str(linkage_tree_time) +"\ntime-dependency:"+ str(dependency_matrix_time) +"\nobjective-linkage:"+ str(linkage_tree_obj) +"\nobjective-dependency:"+  str(dependency_matrix_obj) +"\n"+learned_linkage+ "\n" + dependency
file = open("output.txt", "w")
file.write(str(result))
print("done")