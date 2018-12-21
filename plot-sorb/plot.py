import pylab
import os
import json
import numpy as np

runs = 1
dim = 25
scale = 6

gomea_command = f"./RV-GOMEA  -b -f -8 -s -r 19 {scale} -115 -100 45 0.35 10 25 0.9 1 0 1e-10 100 0.0 1"

heatmap = np.ones((dim, dim))
print(gomea_command)
# directory = f"{os.getcwd()}/heatmaps/{dim}"
# if not os.path.exists(directory):
#     os.makedirs(directory)

heatmaps = {}

for i in range(runs):
    gomea_result = os.popen(gomea_command).readlines()

    # plot heatmap x / y

    matrix = []
    for i in range(dim):
        row = [float(element) for element in gomea_result[1 + i].split(",")[:-1]]
        matrix.append(row)
    matrix = np.array(matrix)
    heatmap += matrix * 1 / runs

    # if result_vec_count[tuple(individual)] > 1:
    #     heatmaps[tuple(individual)] = heatmaps[tuple(individual)] + matrix
    # else:
    #     heatmaps[tuple(individual)] = matrix

pylab.imshow(heatmap)
pylab.title(f"Heatmap of soreb for 10^{scale}\n ")
pylab.savefig(f"heatmaps/{scale}_heatmap.png")
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
