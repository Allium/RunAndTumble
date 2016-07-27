import numpy as np
import matplotlib.pyplot as plt

## Data
t = np.arange(0.01, 5.0, 0.01)
s1 = np.sin(2*np.pi*t)
s2 = np.exp(-t)
s3 = np.sin(6*np.pi*t)

## First subplot
ax1 = plt.subplot(211)
plt.plot(t, s1, c="b")
plt.plot(t, s3, c="r", ls=":", lw=2.0)

## Second subplot
ax2 = plt.subplot(212, sharex=ax1)
plt.plot(t, s2)

plt.show()