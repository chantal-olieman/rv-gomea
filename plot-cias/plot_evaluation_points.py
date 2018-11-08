import pylab
import os
import numpy as np

dim = 10
first_individual = [0] * dim
second_individual = [0] * dim
elitist = [0] * dim
k = 0
j = 1


file = open(f"{os.getcwd()}/input.txt")
content = file.readlines()
for i in range(dim):
    line = content[i]
    splitted = line.split()
    first_individual[i] = float(splitted[1])
    second_individual[i] = float(splitted[3])
    elitist[i] = float(splitted[5])


print(first_individual)
print(second_individual)
print(elitist)

x = []
y = []

individuals = [first_individual, second_individual, elitist]
colors = ['g', 'c', 'r', 'b', 'y']
for j in range(int(dim/2)):
    for i in range(3):
        pylab.scatter(individuals[i][j*2], individuals[i][(j*2)+1], c=colors[j])

# pylab.scatter(x, y)
# pylab.scatter(x[0], y[0], c='g')
#
# individual_to_compare[0] = second_individual[0]
# x = []
# y = []
# for i in range(int(dim/2)):
#     x.append(individual_to_compare[i])
#     y.append(individual_to_compare[i+1])
# pylab.scatter(x[0], y[0], c='r')
#
# individual_to_compare[0] = first_individual[0]
# individual_to_compare[1] = second_individual[1]
# x = []
# y = []
# for i in range(int(dim/2)):
#     x.append(individual_to_compare[i])
#     y.append(individual_to_compare[i+1])
# pylab.scatter(x[0], y[0], c='r')
#
# individual_to_compare[0] = second_individual[0]
# x = []
# y = []
# for i in range(int(dim/2)):
#     x.append(individual_to_compare[i])
#     y.append(individual_to_compare[i+1])
# pylab.scatter(x[0], y[0], c='y')


pylab.title(f"Cias points for n = {len(x)}")
pylab.ylabel(f"y")
pylab.xlabel("x")
# pylab.savefig(f"{directory}/{dim}_{i+1}.png")
pylab.show()

    #individual = [float(element) for element in gomea_result[1].split(",")[:-1]]
    # plot solution
    # x = individual[:int(len(individual) / 2)]
    # y = individual[int(len(individual) / 2):]
    #
