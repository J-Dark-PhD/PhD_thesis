import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append("../")

component_to_test = "everything"

max_size_range = []
inital_value = 0.005
max_size_range.append(inital_value)
# for i in range(1, 35):
for i in range(1, 28):
    max_size_range.append(inital_value * 0.9**i)

max_sizes = []
mesh_sizes = []
sim_times = []

factor = component_to_test + "_max_size"

for size in max_size_range:
    max_sizes.append(size)
    data_folder = "../h_transport_results/{}_{:.1e}/".format(factor, size)
    data_sim = np.genfromtxt(data_folder + "sim_data.csv", delimiter=",", names=True)

    mesh_sizes.append(data_sim["cells"])
    sim_times.append(data_sim["sim_time"])

max_cell_sizes = np.array(max_sizes)
mesh_sizes = np.array(mesh_sizes)
sim_times = np.array(sim_times)

# ##### Plotting ##### #

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

plt.figure()
plt.plot(mesh_sizes, sim_times, marker="x", color="black")
plt.xlabel("Mesh size (cells)")
plt.ylabel("Time (s)")
plt.title("Mesh size vs Simulation time")

plt.show()
