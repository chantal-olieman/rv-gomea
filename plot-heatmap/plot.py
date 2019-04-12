import pylab
import subprocess
import os
import re
from subprocess import call
from subprocess import check_output
import json
import numpy as np

runs = 1
dim = 146
scale = dim
problem = 13
from_memory = 0



gomea_command = f"./RV-GOMEA -f -10  -s -r -b {21+problem} {dim} -100 100 0 0.35 50 25 0.9 1 3000000.0 0.1 100 0.0 1"
# gomea_command = f"./RV-GOMEA -f -8  -s -r -b {13} 10 -100 100 45 0.35 50 25 0.9 1 3000000.0 0.1 100 0.0 1"


gomea_command = " ./RV-GOMEA -f -10  -s -r  0 100 -115 -100 0 0.35 10 25 0.9 1 1500000.0 0.1 100 0.0 1"
# directory = f"{os.getcwd()}/heatmaps/{dim}"
# if not os.path.exists(directory):
#     os.makedirs(directory)

filename = f"{os.getcwd()}/heatmaps/src/matrix_{problem}.txt"
heatmaps = {}
print(gomea_command)
for i in range(runs):
    if from_memory:
        break
    path_command = "export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH"
    gomea_result = subprocess.Popen(path_command + " ; " + gomea_command, cwd=os.getcwd(), shell=True, stdout=subprocess.PIPE)
    out, err = gomea_result.communicate()
    print(out)
    correct = json.dumps(out.decode('utf8'))
    print(correct)
    split_list = correct.replace("n", "").split(",")[-((dim*dim)+1):-1]
    split_list[0] = split_list[0].split(" ")[-1]


    # plot heatmap x / y

    matrix = []
    for i in range(dim):
        # for t in range(dim):
        #     # print(split_list[t])
        row = []
        for j in range(dim):
            res = split_list[j+(i*dim)].replace("\\", "")

            if res:
                if float(res) == 1.0:
                    row.append(0.0)
                else:
                    row.append(float(res))

        matrix.append(row)
    matrix = np.array(matrix)
    print(matrix)
    # heatmap += matrix * 1 / runs

    # if result_vec_count[tuple(individual)] > 1:
    #     heatmaps[tuple(individual)] = heatmaps[tuple(individual)] + matrix
    # else:
    #     heatmaps[tuple(individual)] = matrix

    file = open(filename, "w")
    json.dump(matrix.tolist(), file)

if from_memory:
    file = open(filename)
    matrix = json.load(file)
#
# sorted_matrix = np.zeros((dim,dim))
# order_index = 1
# dependency_index = 0
# order_list = [0]
# order_set = set()
# order_set.add(0)
# for i in range(dim):
#     if i >= order_index:
#         for p in range(dim):
#             if p not in order_set:
#                 order_set.add(p)
#                 order_list.append(p)
#                 order_index += 1
#                 break
#     current_i = order_list[i]
#
#     for j in range(dim):
#         if j not in order_set:
#             if matrix[current_i][j] != 0.000:
#                 order_set.add(j)
#                 order_list.append(j)
#                 order_index += 1
#
#     for j in range(order_index):
#         sorted_matrix[i][j] = matrix[order_list[i]][order_list[j]]
# print(order_list[:100])

# scaled_matrix = [row[700:850] for row in matrix[700:850]]
pylab.imshow(matrix)
# pylab.imshow(scaled_matrix)
pylab.title(f"Heatmap of CEC 2013 f{problem}\n ")
pylab.savefig(f"heatmaps/{problem}_heatmap_{dim}.png")
pylab.show()

# print(len(heatmaps))
# for key, value in heatmaps.items():
#     if result_vec_count[key] > 3 or key == tuple([0.0, 1.0, 0.5, 0.0, 1.0, 0.0, 0.0, 0.5, 1.0, 1.0]):
#         print(result_vec_count[key])
#         print(value)
#         pylab.imshow(value)
#         pylab.title(f"Dependency of points for n = {dim/2}\n {key}")
#         pylab.savefig(f"heatmaps/{dim}_heatmap.png")
#         pylab.show()
