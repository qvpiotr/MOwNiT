import numpy as np
import matplotlib.pyplot as plt
import math
import os

def monteCarlo(N):
    x = np.random.uniform(low = 0, high = 1, size = [N, 1])
    y = np.random.uniform(low = 0, high = 1, size = [N, 1])

    smaller_bool = x**2 > y
    # smaller_bool = 1/math.sqrt(x) > y

    approx = np.sum(smaller_bool)/N
    # print('{} {}'.format(N, abs(1/3 - approx)/(1/3)*100))
    print('{}'.format(approx))

I = [10 ** i for i in range(1,7)]
copyI = []
for i in I:
    copyI.append(i)
for i in copyI:
    I.append(int(i/2))
I.sort()
for i in I:
    monteCarlo(i)

dir_path = os.path.dirname(os.path.realpath(__file__))
f = open(dir_path + '/' + 'data.txt', 'r')
data = np.loadtxt(f)


x = data[:, 0]
y = data[:, 1]

print(x)
print(y)
plt.plot(x, y, 'ro')
plt.xscale('log')
plt.show()


