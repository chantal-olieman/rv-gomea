import json

from matplotlib import pylab
import os

results = {}

def read_in_solution(num_circles):
    solution = []
    file = open(f"/home/chantal/Documents/coordinates/csq{num_circles}.txt")
    content = file.readlines()
    for i in range(num_circles):
        circle = content[i].split()
        solution.append(float(circle[1]))
        solution.append(float(circle[2]))
    # normalize x
    max_x = -1
    min_x = 1
    for i in range(num_circles):
        max_x = max(max_x, solution[i*2])
        min_x = min(min_x, solution[i*2])
    diff = abs(max_x - min_x)
    for i in range(num_circles):
        solution[i*2] = (solution[2*i] - min_x)/diff
    # normalize y
    max_y = 0-1
    min_y = 1
    for i in range(num_circles):
        max_y = max(max_y, solution[(i * 2)+1])
        min_y = min(min_y, solution[(i * 2)+1])
    diff = max_y - min_y
    for i in range(num_circles):
        solution[(i * 2)+1] = (solution[(2 * i)+1] - min_y) / diff
    return solution


def calculate_relaxes_fitness(solution):
    total = 0
    for i in range(int(len(solution) / 2)):
        for j in range(i, int(len(solution) / 2)):
            if i == j:
                continue
            x_diff = solution[i * 2] - solution[j * 2]
            y_diff = solution[(i * 2) + 1] - solution[(j * 2) + 1]
            distance = pow((x_diff * x_diff + y_diff * y_diff), 0.5)
            total += pow(distance, -4)
    return total


for i in range(2, 1000):
    results[str(i)] = str(calculate_relaxes_fitness(read_in_solution(i)))

file = open(f"{os.getcwd()}/optima.txt", "w")
file.write(json.dumps(results))

# dim = 10
# solution = read_in_solution(dim)
# print(solution)
#
# # print(f"The total results = {calculate_relaxes_fitness(solution)}")
#
# for i in range(int(len(solution) / 2)):
#     pylab.scatter(solution[(i * 2)], solution[(i * 2) + 1])
# pylab.title(f"Cias points for n = {dim}")
# pylab.ylabel(f"y")
# pylab.xlabel("x")
# pylab.show()
