import matplotlib.pyplot as plt
import numpy as np

## Create data
dt = 0.01
t = np.arange(0.0, 5.0, 0.01)
r = np.exp(-t)
x = 0.02*np.random.randn(len(t))

## Plot
plt.plot(t, r+x)
plt.xlabel("Time (s)")
plt.ylabel("Signal")
plt.title("Noisy decay")

# plt.savefig("test.png")
# print "Figure saved as test.png"

plt.show()