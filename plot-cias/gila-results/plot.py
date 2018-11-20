import json
import pylab


def read_in_and_write(measure_type):
    file = open(f"{measure_type}_2.txt")
    results = json.load(file)
    index = 0
    population = results[index]
    univariate = results[index+1]
    index += 1
    univariate2 = results[index+1]
    differnetial_grouping = results[index+2]
    original_diff_grouping = results[index+3]
    epsi_05 = results[index+4]
    epsi_01 = results[index+5]
    epsi_005 = results[index+6]
    # no_pruning = results[index+7]

    pylab.loglog(population[:len(univariate)], univariate, label='Uni-variate model',  marker='o')
    pylab.loglog(population[:len(univariate2)], univariate2, label='Uni-variate 2 model',  marker='o')
    pylab.loglog(population[:len(differnetial_grouping)], differnetial_grouping, label='Differential grouping model',  marker='o')
    # pylab.loglog(population[:len(mutial_information) ], mutial_information[:], label='Mutual Information Linkage',
                # marker='o')
    # pylab.loglog(population[:len(evolve)], evolve, label='Evolving factor 1.0', marker='o')
    pylab.loglog(population[:len(original_diff_grouping)], original_diff_grouping, label='Original diff grouping', marker='o')
    pylab.loglog(population[:len(epsi_05)], epsi_05, label='evolve pruning 0.5', marker='o')
    pylab.loglog(population[:len(epsi_01)], epsi_01, label='evolve pruning 0.1', marker='o')
    pylab.loglog(population[:len(epsi_005)], epsi_005, label='evolve pruning 0.05', marker='o')
    # pylab.loglog(population[:len(no_pruning)], no_pruning, label='no pruning', marker='o')
    # pylab.loglog(population[:len(evolve_adapt)], evolve_adapt,  label='Adaptive evolving with factor 1.0',  marker='o')
    # pylab.loglog(population[:len(evolve_adapt_p2)], evolve_adapt_p2,  label='Adaptive evolving with pruning, sets of 1',  marker='o')
    # pylab.loglog(population[:len(evolve_adapt_p1)], evolve_adapt_p1,  label='Adaptive evolving with pruning, no minimal set size',  marker='o')
    blackbox_string = ""
    problem = "cias"

    if measure_type == "times":
        pylab.title(f"Runtime on the {problem} problem, based on parameters {blackbox_string}")
        pylab.ylabel(f"time (s)")
    elif measure_type == "evals":
        pylab.title(f"Number of evaluations on the {problem} problem, based on parameters {blackbox_string}")
        pylab.ylabel(f"Number of evaluations")

    pylab.xlabel("Problem size")

    pylab.legend(loc='upper left')
    pylab.savefig(f"cias_{measure_type}")
    pylab.show()


# problem = input("problem? \n")  # Python 3
# blackbox = input("_blackbox or empty string? \n")  # Python 3

# problem = "sphere"
# blackbox = "_blackbox"
# # blackbox = ""
read_in_and_write("evals")
read_in_and_write("times")
