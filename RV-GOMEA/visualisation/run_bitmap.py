import os
import numpy as np


runs = 1


population = 10
rotation = 45
problem = 13
bitmap = np.zeros((population, population))


learned_linkage = f"./RV-GOMEA -s -r -f -7 {problem} {population} -115 -100 {rotation} " \
                  f"0.35 10 25 0.9 1 0 1e-10 100 0.0 300"
print(f"settings: {learned_linkage}")


for i in range(runs):
    results = os.popen(learned_linkage).readlines()
    if i % 5 == 0:
        print(f"run {i}")
    try:
        for i in range(population):
            row = results[i+1].split(", ")[:-1]
            for j in range(population):
                bitmap[i][j] = bitmap[i][j] + float(row[j])*float(1/runs)


    except Exception as e:
        print(e)
        print(f"error: {results}")

print(bitmap)
for i in range(population):
    bitmap[i][i] = bitmap.max()

from matplotlib import pyplot as plt
img = plt.imshow(bitmap)
plt.savefig("bitmap_soREB_10.png")
plt.show()

# 
# file = open("bitmap_sphere_10.txt", "w")
# file.write(str(bitmap))
print("done")