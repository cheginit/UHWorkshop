import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/residual", dtype=np.float)

c = 10
plt.plot(data[c:, 0], data[c:, 1], label="tot")
plt.plot(data[c:, 0], data[c:, 2], label="u")
plt.plot(data[c:, 0], data[c:, 3], label="v")
plt.plot(data[c:, 0], data[c:, 4], label="p")
plt.plot(data[c:, 0], data[c:, 4], label="div")

plt.legend(loc="best")
plt.show()
