import pylab
import os
import json
import numpy as np

runs = 1
dim = 10
vtr_file = open(f"{os.getcwd()}/optimal/optima.txt")
vtr = json.load(vtr_file)

original_vtr_file = open(f"{os.getcwd()}/optimal/optima_original.txt")
vtr = json.load(vtr_file)

gomea_command = f"./RV-GOMEA -s -r -b -f -8 16 {dim*2} 0 1 0 0.35 10 25 0.9 1 0 {float(vtr[str(dim)])*1.0001} 100 0.0 30"
dim = dim*2
heatmap = np.ones((dim, dim))
print(gomea_command)
directory = f"{os.getcwd()}/circles/{dim}"
if not os.path.exists(directory):
    os.makedirs(directory)

result_vec_count = {}
heatmaps = {}

for i in range(runs):
    gomea_result = os.popen(gomea_command).readlines()
    individual = [float(element) for element in gomea_result[1].split(",")[:-1]]
    print(gomea_result[0])
    result_vec_count[tuple(individual)] = result_vec_count.get(tuple(individual), 0) + 1
    # plot solution
    print(individual)
    x = individual[:int(len(individual) / 2)]
    y = individual[int(len(individual) / 2):]
    pylab.scatter(x, y)
    pylab.title(f"Cias points for n = {len(x)}")
    pylab.ylabel(f"y")
    pylab.xlabel("x")
    pylab.savefig(f"{directory}/{dim}_{i+1}.png")
    pylab.show()
    # plot heatmap x / y

    matrix = []
    for i in range(len(individual)):
        row = [float(element) for element in gomea_result[3 + i].split(",")[:-1]]
        matrix.append(row)
    matrix = np.array(matrix)
    heatmap += matrix * 1 / runs

    if result_vec_count[tuple(individual)] > 1:
        heatmaps[tuple(individual)] = heatmaps[tuple(individual)] + matrix
    else:
        heatmaps[tuple(individual)] = matrix

pylab.imshow(heatmap)
pylab.title(f"Dependency of points for n = {dim/2}\n ")
pylab.savefig(f"heatmaps/{dim}_heatmap.png")
pylab.show()
#
# combination_map = {(0.0, 0.0): 0, (0.0, 1.0): 1, (0.5, 0.5): 2, (1.0, 0.0): 3, (1.0, 1.0): 4}
#
# print(len(heatmaps))
# for key, value in heatmaps.items():
#     if result_vec_count[key] > 3 or key == tuple([0.0, 1.0, 0.5, 0.0, 1.0, 0.0, 0.0, 0.5, 1.0, 1.0]):
#         print(result_vec_count[key])
#         print(value)
#         pylab.imshow(value)
#         pylab.title(f"Dependency of points for n = {dim/2}\n {key}")
#         pylab.savefig(f"heatmaps/{dim}_heatmap.png")
#         pylab.show()
