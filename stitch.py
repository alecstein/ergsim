import matplotlib.pyplot as plt
import matplotlib.animation as animation
import signal
import argparse

signal.signal(signal.SIGINT, signal.SIG_DFL)

# Parse args from command line

parser = argparse.ArgumentParser(description='Stitch histograms together')
parser.add_argument('n', type=int, nargs=1, help='Number of particles')
parser.add_argument('t', type=int, nargs=1, help='Number of timesteps')
args = parser.parse_args()
n = args.n[0]
t = args.t[0]


filename = f"N_{n}_Time_{t}_x.txt"
f = open(filename, 'r')
raw_data = []
for line in f:
    raw_data.append(float(line.strip()))

# break the raw data into subslices
data = [raw_data[i:i + n] for i in range(0, len(raw_data), n)]

fig, ax = plt.subplots(figsize=(8, 5))

bin_ct = 100
bin_arr = [0.01 * x for x in range(0, bin_ct)]


def animate(i):
    ax.clear()
    ax.hist(data[i], density=True, histtype='step', bins=bin_arr, label = f"t = {i}")
    ax.set_title(f"Position distribution, ${n}$ particles\n $v_0 = const$")
    ax.set_xlabel("$x$ (m)")
    ax.set_ylabel("$p(x)$")
    ax.set_ylim(0, 5)
    ax.legend()
    ax.set_xlim(0, 1)


# If the time is t, and there are 100 timesteps, then the animation will be
# 1000 frames long. The time per frame is 1/100.
    

ani = animation.FuncAnimation(fig,
                              animate,
                              frames=len(data),
                              repeat=True,
                              interval=50)
ani.save(f'animation_{n}_{t}_positions.mp4')

plt.show()
