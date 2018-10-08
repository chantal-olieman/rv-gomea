import json
import pylab

mesure_type = "time"
problem = "rosenbrock"

file = open(f"output_{mesure_type}_{problem}.txt")
results = json.load(file)
population = results[0]
univariate = results[1]
linkage_learning = results[2]
dependency_learning = results[3]

pylab.loglog(population[:], univariate, '-g', label='Uni-variate model',  marker='o')
pylab.loglog(population, linkage_learning, '-b', label='Mutial information',  marker='o')
pylab.loglog(population, dependency_learning, '-r', label='Differential grouping',  marker='o')
# pylab.title(f"{mesure_type} for the {problem} problem")
pylab.title(f"Number of evaluations on the {problem} problem")
# pylab.title(f"Runtime on the {problem} problem")
pylab.xlabel("Population size")
# pylab.ylabel(f"{mesure_type}")
pylab.ylabel(f"Number of evaluations")
# pylab.ylabel(f"time (s)")
pylab.legend(loc='upper left')
pylab.savefig(f"{problem}_{mesure_type}")
pylab.show()