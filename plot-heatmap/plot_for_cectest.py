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
problem = 34
max_eval = "10e5"
from_memory = 0

flag = -8
method = "FB"

# gomea_command = f"./RV-GOMEA-CECs -f {flag} -s -r {problem} 1000 -100 -100 0 0.35 50 25 0.9 1 3000000 0.1 100 0.0 1"
# gomea_command = f"./RV-GOMEA -f -8  -s -r -b {13} 10 -100 100 45 0.35 50 25 0.9 1 3000000.0 0.1 100 0.0 1"

gomea_command = f"./RV-GOMEA-CECs -f -8  -s -r  {problem} 1000 -100 100 0 0.35 10 25 0.9 1 1500000.0 0.1 100 0.0 1"
# gomea_command = " ./RV-GOMEA -f -10  -s -r  0 100 -115 -100 0 0.35 10 25 0.9 1 1500000.0 0.1 100 0.0 1"
# directory = f"{os.getcwd()}/heatmaps/{dim}"
# if not os.path.exists(directory):
#     os.makedirs(directory)
# gomea_command = "./RV-GOMEA-CECs -f -8  -s -r  34 1000 -100 100 0 0.35 10 25 0.9 1 1500000.0 0.1 100 0.0 1"

filename = f"{os.getcwd()}/CEC/src/benchmark_{problem}.txt"
heatmap = np.zeros((dim,dim))
print(gomea_command)
for i in range(runs):
    if from_memory:
        break
    path_command = "export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH"
    gomea_result = subprocess.Popen(path_command + " ; " + gomea_command, cwd=os.getcwd(), shell=True, stdout=subprocess.PIPE)
    out, err = gomea_result.communicate()
    # print(out)
    correct = json.dumps(out.decode('utf8'))
    print(correct)
    # print(correct.replace("n", "").split(","))
    split_list = correct.replace("n", "").split(",")[-((dim*dim)+1):-1]
    split_list[0] = split_list[0].split(" ")[-1]
    # print(correct[-4])

    # plot heatmap x / y
    # print("splitting and working")
    matrix = []
    print(split_list)
    # print(len(split_list))
    # print(len(split_list[400]))
    for i in range(dim):
        # for t in range(dim):
        #     # print(split_list[t])
        row = []
        for j in range(dim):
            res = split_list[j+(i*dim)].replace("\\", "")

            if res:
                if float(res) == 1.0:
                    row.append(1.0)
                else:
                    row.append(float(res))

        matrix.append(row)
    matrix = np.array(matrix)
    # print(matrix)
    # print("done?")
    heatmap += matrix * 1 / runs
    print(heatmap)
    pylab.imshow(heatmap)
    # pylab.title(f"Heatmap of 50 dimensional {problem} problem\n ")
    pylab.savefig(f"CEC/{problem}.png")
    # pylab.plt.colorbar(orientation="horizontal")
    pylab.show()
    # if result_vec_count[tuple(individual)] > 1:
    #     heatmaps[tuple(individual)] = heatmaps[tuple(individual)] + matrix
    # else:
    #     heatmaps[tuple(individual)] = matrix
    file = open(filename, "w")
    json.dump(heatmap.tolist(), file)

if from_memory:
    file = open(filename)
    heatmap = json.load(file)

# heatmap[0,0] = 1.0

pylab.imshow(heatmap)
# pylab.title(f"Heatmap of 50 dimensional {problem} problem\n ")
pylab.savefig(f"CEC/{problem}.png")
# pylab.plt.colorbar(orientation="horizontal")
pylab.show()

scaled_matrix = [row[:50] for row in heatmap[:50]]
# pylab.imshow(heatmap)
pylab.imshow(scaled_matrix)
# pylab.title(f"Heatmap of 50 dimensional {problem} problem\n ")
pylab.savefig(f"CEC/{problem}_50.png")
pylab.show()

