import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append("../../")

from parametric_values import (
    max_size_range,
)

# for value in max_size_range:
#     print("{:.2e}".format(value))
# quit()

component_to_test = "lipb"

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
    data_folder = "../h_transport_results/{}_{:.2e}/".format(factor, size)
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

ref_lipb = inventory_lipb[0]
ref_structure = inventory_structure[0]
ref_baffle = inventory_baffle[0]
ref_tungsten = inventory_tungsten[0]
ref_pipes = inventory_pipes[0]

inventory_lipb_normalised = np.abs(inventory_lipb / ref_lipb) - 1
inventory_structure_normalised = np.abs(inventory_structure / ref_structure) - 1
inventory_baffle_normalised = np.abs(inventory_baffle / ref_baffle) - 1
inventory_tungsten_normalised = np.abs(inventory_tungsten / ref_tungsten) - 1
inventory_pipes_normalised = np.abs(inventory_pipes / ref_pipes) - 1

# ##### Plotting ##### #

# plt.rc("text", usetex=True)
# plt.rc("font", family="serif", size=12)

plt.figure(figsize=[8, 4.8])
plt.plot(max_size_range, inventory_lipb_normalised, marker="x", label="Lipb")
plt.plot(max_size_range, inventory_structure_normalised, marker="x", label="Structure")
plt.plot(max_size_range, inventory_baffle_normalised, marker="x", label="Baffle Plate")
plt.plot(max_size_range, inventory_tungsten_normalised, marker="x", label="Firstwall")
plt.plot(max_size_range, inventory_pipes_normalised, marker="x", label="BZ Pipes")
alpha = 0.25
plt.axhspan(
    0.02,
    -0.02,
    color="tab:blue",
    alpha=alpha,
)
plt.axhspan(
    0.01,
    -0.01,
    color="tab:blue",
    alpha=alpha,
)
plt.axhspan(
    0.005,
    -0.005,
    color="tab:blue",
    alpha=alpha,
)
plt.axhspan(
    0.001,
    -0.001,
    color="tab:blue",
    alpha=alpha,
)
leg = plt.legend(loc="upper left", bbox_to_anchor=(1.02, 1.05))
# plt.ylim(-0.01, 0.01)
plt.xlabel("Maxiumum cell size (m)")
plt.ylabel("Inventory (H m$^{-3}$)")
plt.title("Lipb increasing")
plt.tight_layout()


plt.show()
