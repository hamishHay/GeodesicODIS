import numpy as np
import matplotlib.pyplot as plt
import time


# data = np.loadtxt("OUTPUT.txt", skiprows=77, usecols=(3, 6), delimiter=' ').T
success = False
skip = 90
# while not success:
try:
    success = True
    data = np.genfromtxt("OUTPUT.txt", skip_header=skip, skip_footer=6, usecols=(3,6)).T
    skip += 1
except ValueError:
    success = False
        # continue

t = data[0,:]
diss = data[1,:]*1e9
print(diss)
plt.ion()
fig, ax = plt.subplots(figsize=(16,4))

line1, = ax.semilogy(t, diss)

ax.set_xlabel("Time [$t/T$]")
ax.set_ylabel("Dissipated Energy")


success = False
while True:
    data = np.genfromtxt("OUTPUT.txt", skip_header=skip, skip_footer=6, usecols=(3,6)).T

    try:
        t = data[0,:]
        diss = data[1,:]*1e9

        line1.set_xdata(t)
        line1.set_ydata(diss)

        ax.set_xlim([0, t[-1]])
        ax.set_ylim([np.max(diss)*1.1, np.min(diss[10:])*0.9][::-1])

        fig.canvas.draw()
        fig.canvas.flush_events()


        time.sleep(0.1)
    except IndexError:
        time.sleep(2)

plt.show()
