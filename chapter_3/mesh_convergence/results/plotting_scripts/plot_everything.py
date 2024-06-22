import matplotlib.pyplot as plt
import numpy as np

component_to_test = "everything"

max_size_range = []
inital_value = 0.005
max_size_range.append(inital_value)
# for i in range(1, 35):
for i in range(1, 28):
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
reference_value_lipb = inventory_lipb[-1]
reference_value_structure = inventory_structure[-1]
reference_value_baffle = inventory_baffle[-1]
reference_value_tungsten = inventory_tungsten[-1]
reference_value_pipes = inventory_pipes[-1]

# ##### Plotting ##### #

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

plt.figure()
plt.plot(mesh_sizes, inventory_lipb, marker="x", color="black")
plt.xlabel("Mesh size (cells)")
plt.ylabel("Inventory (H m$^{-1}$)")
plt.xscale("log")
plt.ylim(0, 6e21)
plt.title("Lipb")
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)


plt.figure()
plt.plot(mesh_sizes, inventory_structure, marker="x", color="black")
plt.xlabel("Mesh size (cells)")
plt.ylabel("Inventory (H m$^{-1}$)")
plt.xscale("log")
plt.ylim(0, 3e20)
plt.title("Structure")
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)


plt.figure()
plt.plot(mesh_sizes, inventory_baffle, marker="x", color="black")
plt.xlabel("Mesh size (cells)")
plt.ylabel("Inventory (H m$^{-1}$)")
plt.xscale("log")
plt.ylim(0, 2e18)
plt.title("Baffle plate")
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)


plt.figure()
plt.plot(mesh_sizes, inventory_tungsten, marker="x", color="black")
plt.xlabel("Mesh size (cells)")
plt.ylabel("Inventory (H m$^{-1}$)")
plt.xscale("log")
plt.ylim(0, 5e17)
plt.title("First wall")
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)


plt.figure()
plt.plot(mesh_sizes, inventory_pipes, marker="x", color="black")
plt.xlabel("Mesh size (cells)")
plt.ylabel("Inventory (H m$^{-1}$)")
plt.xscale("log")
plt.ylim(0, 1.5e19)
plt.title("BZ Pipes")
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)


plt.show()
