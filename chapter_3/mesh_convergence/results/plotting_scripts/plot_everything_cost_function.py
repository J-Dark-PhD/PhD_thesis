import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append("../")

component_to_test = "everything"

max_size_range = []
inital_value = 0.005
max_size_range.append(inital_value)
# for i in range(1, 35):
for i in range(1, 25):
    max_size_range.append(inital_value * 0.9**i)

inventory_lipb = []
inventory_structure = []
inventory_baffle = []
inventory_tungsten = []
inventory_pipes = []
max_sizes = []
mesh_sizes = []
sim_times = []

factor = component_to_test + "_max_size"

for size in max_size_range:
    max_sizes.append(size)
    data_folder = "../h_transport_results/{}_{:.1e}/".format(factor, size)
    data = np.genfromtxt(
        data_folder + "derived_quantities.csv", delimiter=",", names=True
    )
    data_sim = np.genfromtxt(data_folder + "sim_data.csv", delimiter=",", names=True)

    inventory_lipb.append(data["Total_solute_volume_6"])
    inventory_structure.append(data["Total_retention_volume_7"])
    inventory_baffle.append(data["Total_retention_volume_8"])
    inventory_tungsten.append(data["Total_retention_volume_9"])
    inventory_pipes.append(
        sum(data["Total_retention_volume_{}".format(i)] for i in range(10, 21))
    )

    mesh_sizes.append(data_sim["cells"])
    sim_times.append(data_sim["sim_time"])

inventory_lipb = np.array(inventory_lipb)
inventory_structure = np.array(inventory_structure)
inventory_baffle = np.array(inventory_baffle)
inventory_tungsten = np.array(inventory_tungsten)
inventory_pipes = np.array(inventory_pipes)
max_cell_sizes = np.array(max_sizes)
mesh_sizes = np.array(mesh_sizes)
sim_times = np.array(sim_times)
# y = np.abs(np.diff(inventories) / np.diff(mesh_sizes))
# standard_value = inventories.min()
reference_value_lipb = inventory_lipb[-1]
reference_value_structure = inventory_structure[-1]
reference_value_baffle = inventory_baffle[-1]
reference_value_tungsten = inventory_tungsten[-1]
reference_value_pipes = inventory_pipes[-1]

y_lipb = np.abs(np.diff(inventory_lipb) / np.diff(mesh_sizes))
y_structure = np.abs(np.diff(inventory_structure) / np.diff(mesh_sizes))
y_baffle = np.abs(np.diff(inventory_baffle) / np.diff(mesh_sizes))
y_tungsten = np.abs(np.diff(inventory_tungsten) / np.diff(mesh_sizes))
y_pipes = np.abs(np.diff(inventory_pipes) / np.diff(mesh_sizes))

# ##### Plotting ##### #

# plt.rc("text", usetex=True)
# plt.rc("font", family="serif", size=12)

plt.figure()
plt.plot(mesh_sizes[:-1], y_lipb * sim_times[:-1], marker="x", color="black")
ax = plt.gca()
ax.invert_xaxis()
ax.set_xscale("log")
ax.set_xlabel("Mesh size (cells)")
plt.ylabel("cost x time")
plt.title("Lipb")

plt.figure()
plt.plot(mesh_sizes[:-1], y_structure * sim_times[:-1], marker="x", color="black")
ax = plt.gca()
ax.invert_xaxis()
ax.set_xscale("log")
ax.set_xlabel("Mesh size (cells)")
plt.ylabel("cost x time")
plt.title("Structure")

plt.figure()
plt.plot(mesh_sizes[:-1], y_baffle * sim_times[:-1], marker="x", color="black")
ax = plt.gca()
ax.invert_xaxis()
ax.set_xscale("log")
ax.set_xlabel("Mesh size (cells)")
plt.ylabel("cost x time")
plt.title("Baffle Plate")

plt.figure()
plt.plot(mesh_sizes[:-1], y_tungsten * sim_times[:-1], marker="x", color="black")
ax = plt.gca()
ax.invert_xaxis()
ax.set_xscale("log")
ax.set_xlabel("Mesh size (cells)")
plt.ylabel("cost x time")
plt.title("First Wall")

plt.figure()
plt.plot(mesh_sizes[:-1], y_pipes * sim_times[:-1], marker="x", color="black")
ax = plt.gca()
ax.invert_xaxis()
ax.set_xscale("log")
ax.set_xlabel("Mesh size (cells)")
plt.ylabel("cost x time")
plt.title("BZ Pipes")

plt.show()
