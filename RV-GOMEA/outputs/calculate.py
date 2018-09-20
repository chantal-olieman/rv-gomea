time_linkage = [7.856, 7.098, 12.214, 9.222, 8.163, 30.002, 6.285, 8.144, 7.333, 6.55, 4.704, 8.839, 7.572, 8.687, 11.053, 20.431, 9.492]
time_dependency= [3.045, 3.562, 8.713, 3.075, 3.702, 3.551, 2.354, 10.543, 4.307, 2.419, 3.776, 2.649, 6.498, 2.175, 4.291, 2.806, 3.398, 2.987, 3.089, 2.127]
objective_linkage= [9.56e-11, 9.98e-11, 9.58e-11, 9.97e-11, 1e-10, 3.99, 8.99e-11, 9.91e-11, 9.97e-11, 9.95e-11, 9.97e-11, 9.77e-11, 9.99e-11, 8.72e-11, 9.8e-11, 9.81e-11, 1e-10]
objective_dependency=[9.55e-11, 9.52e-11, 8.32e-11, 9.97e-11, 9.97e-11, 9.99e-11, 9.91e-11, 9.83e-11, 9.99e-11, 9.83e-11, 9.77e-11, 9.97e-11, 9.02e-11, 1e-10, 1e-10, 9.79e-11, 6.94e-11, 9.43e-11, 9.95e-11, 8.67e-11]

print(len(time_linkage))
print(len(objective_linkage))
print(len(time_dependency))
print(len(objective_dependency))

import numpy

del[time_linkage[5]]
del[objective_linkage[5]]

print(numpy.mean(time_linkage))
print(numpy.mean(time_dependency))
print(numpy.mean(objective_linkage))
print(numpy.mean(objective_dependency))