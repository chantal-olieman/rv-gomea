import pylab
import os
import numpy as np

runs = 10
dim = 72
gomea_command = f"./RV-GOMEA -s -r -f -8 16 {dim} 0 1 0 0.35 10 25 0.9 1 0 50600 100 0.0 20"
heatmap = np.ones((dim,dim))
print(gomea_command)
directory = f"{os.getcwd()}/circles/{dim}"
if not os.path.exists(directory):
    os.makedirs(directory)

for i in range(runs):
    gomea_result = os.popen(gomea_command).readlines()
    individual = [float(element) for element in gomea_result[1].split(",")[:-1]]

    print(gomea_result[0])

    # plot solution
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
    heatmap += matrix*1/runs


pylab.imshow(heatmap)
pylab.title(f"Dependency of points for n = {len(x)}")
pylab.savefig(f"heatmaps/{dim}_heatmap.png")
pylab.show()

