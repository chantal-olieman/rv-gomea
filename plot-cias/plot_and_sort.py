import pylab
import os
import numpy as np

runs = 200
dim = 10
gomea_command = f"./RV-GOMEA-check-at-end -s -r -f -8 16 {dim} 0 1 0 0.35 10 25 0.9 1 0 20.5 100 0.0 20"
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
    # x = individual[:int(len(individual) / 2)]
    # y = individual[int(len(individual) / 2):]
    # pylab.scatter(x, y)
    # pylab.title(f"Cias points for n = {len(x)}")
    # pylab.ylabel(f"y")
    # pylab.xlabel("x")
    # pylab.savefig(f"{directory}/{dim}_{i+1}.png")
    # pylab.show()
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

combination_map = {(0.0, 0.0): 0, (0.0, 1.0): 1, (0.5, 0.5): 2, (1.0, 0.0): 3, (1.0, 1.0): 4}
total_heatmap = np.ones((dim, dim))
for key, value in heatmaps.items():
    print(key)
    for i in range(int(dim/2)):
        for j in range(i, int(dim/2)):
            x_tuple = combination_map[tuple([key[i], key[int(dim/2)+i]])]
            y_tuple = combination_map[tuple([key[j], key[int(dim/2)+j]])]
            total_heatmap[2*x_tuple][2*y_tuple] += value[2*i][2*j]
            total_heatmap[(2*x_tuple)+1][(2*y_tuple)] += value[(2*i)+1][(2*j)]
            total_heatmap[(2*x_tuple)][(2*y_tuple)+1] += value[(2*i)][(2*j)+1]
            total_heatmap[(2*x_tuple)+1][(2*y_tuple)+1] += value[(2*i)+1][(2*j)+1]

pylab.imshow(total_heatmap)
pylab.title(f"Dependency of total heatmap")
pylab.savefig(f"heatmaps/{dim}_heatmap.png")
pylab.show()

print(len(heatmaps))
for key, value in heatmaps.items():
    if result_vec_count[key] > 3 or key == tuple([0.0, 1.0, 0.5, 0.0, 1.0, 0.0, 0.0, 0.5, 1.0, 1.0]):
        print(result_vec_count[key])
        print(value)
        pylab.imshow(value)
        pylab.title(f"Dependency of points for n = {dim/2}\n {key}")
        pylab.savefig(f"heatmaps/{dim}_heatmap.png")
        pylab.show()
