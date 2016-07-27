import matplotlib.pyplot as plt
import numpy as np


## 1D histogram
x = np.random.rand(1000)
plt.hist(x, bins=40)

plt.show()

## 2D histogram
x = np.random.randn(1000)
y = np.random.randn(1000) + 5

plt.hist2d(x, y, bins=40)
plt.show()

