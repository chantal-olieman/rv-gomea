import pylab
import os
import json
import numpy as np

runs = 2000
dim = 5
vtr_file = open(f"{os.getcwd()}/optimal/optima.txt")
vtr = json.load(vtr_file)

original_vtr_file = open(f"{os.getcwd()}/optimal/optima_original.txt")
original_vtr = json.load(original_vtr_file)

gomea_command = f"./RV-GOMEA -s -r -b -f -8 17 {dim*2} 0 1 0 0.35 10 25 0.9 1 0 {float(original_vtr[str(dim)])*-0.9999} 100 0.0 30"
# gomea_command = f"./RV-GOMEA -s -w -r -b -f -8 16 {dim*2} 0 1 0 0.35 10 25 0.9 1 0 {float(vtr[str(dim)])*1.001} 100 0.0 30"
# gomea_command = f"./RV-GOMEA -s -r -b -f -8 14 {dim*2} 0 1 0 0.35 10 25 0.9 1 0 -10e3 100 0.0 30"
# gomea_command = "./RV-GOMEA -b -f -14 -s -r 14 10 0 1 0 0.35 100 1 0.9 1 0 -0.7070360705084289 100 0.0 3600 &"
dim = dim*2
heatmap = np.ones((dim, dim))
print(gomea_command)
directory = f"{os.getcwd()}/circles/{dim}"
if not os.path.exists(directory):
    os.makedirs(directory)

result_vec_count = {}
heatmaps = {}


def calculate_original_fitness(solution):
    min_distance = 1
    num_circles = int(len(solution)/2)
    for i in range(num_circles):
        for j in range(i, num_circles):
            if i == j:
                continue
            x_diff = solution[i] - solution[j]
            y_diff = solution[i + num_circles] - solution[j + num_circles]
            distance = pow((x_diff * x_diff + y_diff * y_diff), 0.5)
            min_distance = min(min_distance, distance)
    print(f"\n\nfound min distance: {min_distance}, optimal min distance: {original_vtr[str(num_circles)]} "
          f"\nwhich is {min_distance/float(original_vtr[str(num_circles)]) } ")
    return min_distance


for i in range(runs):
    gomea_result = os.popen(gomea_command).readlines()
    individual = [float(element) for element in gomea_result[1].split(",")[:-1]]
    print(gomea_result[0])
    # plot solution
    individual = [round(i,1) for i in individual]
    print(individual)
    result_vec_count[tuple(individual)] = result_vec_count.get(tuple(individual), 0) + 1
    # calculate_original_fitness(individual)
    # x = individual[:int(len(individual) / 2)]
    # y = individual[int(len(individual) / 2):]
    # pylab.scatter(x, y)
    # pylab.title(f"Cias points for n = {len(x)}")
    # pylab.ylabel(f"y")
    # pylab.xlabel("x")
    # pylab.savefig(f"{directory}/{dim}_{i+1}.png")
    # pylab.show()
    # # plot heatmap x / y

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

print(len(heatmaps))
for key, value in heatmaps.items():
    if result_vec_count[key] > 20 or key == tuple([0.0, 1.0, 0.5, 0.0, 1.0, 0.0, 0.0, 0.5, 1.0, 1.0]):
        print(result_vec_count[key])
        pylab.imshow(value)
        key_show = [tuple([key[i], key[int(i+dim/2)]]) for i in range(int(dim/2))]
        pylab.title(f"Dependency of points for n = {dim/2} seen {result_vec_count[key]}\n {key_show}")
        pylab.savefig(f"heatmaps/{dim}_heatmap.png")
        pylab.show()
