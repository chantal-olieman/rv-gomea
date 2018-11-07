import os
import numpy



def run_tests(problem, rotation, blackbox):
    correct = "EVERYTHING GOOD"
    runs = 10
    timeout = 180

    linkage_options = [-8]
    populations = [10]
    while populations[-1] < 80:
        populations.append(populations[-1] * 2)

    # result_times = [[]] * le(linkage_options)
    # result_evals = [[]] * len(linkage_options)

    result_times = [[] for i in range(len(linkage_options))]
    result_evals = [[] for i in range(len(linkage_options))]

    normal_evals = {0:[[10, 20, 40, 80], [622.3, 673.9, 732.5, 787.3]],
                        7:[[10, 20, 40, 80], [176748.2, 366344.1, 153146.5, 137631.5]],
                         13:[[10, 20, 40, 80], [822.6, 862.7, 948.1, 1020.5]]}
    normal_time = {0: [[10, 20, 40, 80], [0.0098, 0.01219, 0.01452, 0.01781]],
                       7:[[10, 20, 40, 80], [0.18324, 0.5345, 0.61005, 1.08594]],
                   13:[[10, 20, 40, 80], [0.01265, 0.01685, 0.02214, 0.03186]]}

    result_times.insert(0, [])
    result_evals.insert(0, [])


    problem_names = {0: "sphere", 7: "rosenbrock", 13: "soREB", 11: "michalewicz", 12:"rastrigin"}
    print(problem_names[problem])
    linkage_options_timed_out = set()



    for population in populations:
        result_times[0].append(population)
        result_evals[0].append(population)
        for i in range(len(linkage_options)):
            linkage = linkage_options[i]
            if linkage in linkage_options_timed_out:
                continue
            run_times = []
            run_evals = []
            learned_linkage = f"./RV-GOMEA -s -r {blackbox} -f {linkage} {problem} {population+1} -115 -100 {rotation} " \
                              f"0.35 10 25 0.9 1 0 1e-10 100 0.0 600000"
            # print(f"settings: {learned_linkage}")
            j = 0
            while j < runs:
                results = os.popen(learned_linkage).readlines()

                try:
                    j += 1
                    run_evals.append(float(results[0].split(" ")[1]))
                    run_times.append(float(results[0].split(" ")[5]))
                except Exception as e:
                    linkage_options_timed_out.add(linkage)
                    print(f"error: {results}, {e}, {learned_linkage}")
            evaluations = numpy.mean(run_evals)
            time = numpy.mean(run_times)
            if time > timeout:
                linkage_options_timed_out.add(linkage)

            result_times[i + 1].append(round(time, 5))
            result_evals[i + 1].append(round(evaluations,1))

    correct_run = 1
    for i in range(len(populations)):
        if result_evals[1][i]<normal_evals[problem][1][i]*1.3:
            pass
        else:
            correct = "there is a mistake"
            correct_run = 0
            print(f"Not ok {problem}, {populations[i]}")

    for i in range(len(populations)):
        if result_times[1][i]<normal_time[problem][1][i]*1.3:
            pass
        else:
            correct_run = 0
            correct = "there is a mistake"
            print(f"Not ok {problem}, {populations[i]}")

    if not correct_run:
        print(f"ours:   {result_evals}")
        print(f"normal: {normal_evals[problem]}")


        print(f"ours:   {result_times}")
        print(f"normal: {normal_time[problem]}")

    print(correct)



# run_tests(13, 45, "")
run_tests(13, 0, "")
# run_tests(11, 0, "")
run_tests(7, 0, "")
run_tests(0, 0, "")
