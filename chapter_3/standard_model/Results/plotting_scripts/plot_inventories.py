import matplotlib.pyplot as plt
import numpy as np

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=12)

red_W = (171 / 255, 15 / 255, 26 / 255)
dark_factor = 0.8
grey_eurofer = (
    dark_factor * 153 / 255,
    dark_factor * 153 / 255,
    dark_factor * 153 / 255,
)
green_lipb = (dark_factor * 146 / 255, dark_factor * 196 / 255, dark_factor * 125 / 255)


# ## Read data
volumes_bz_pipes = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

folder = "../transient_15d"

# traps
data_traps = np.genfromtxt(
    folder + "/derived_quantities.csv", delimiter=",", names=True
)
t_traps = data_traps["ts"] / 86400
lipb_traps = data_traps["Total_solute_volume_6"]
strucure_traps = (
    data_traps["Total_retention_volume_8"] + data_traps["Total_retention_volume_9"]
)
first_wall_traps = data_traps["Total_retention_volume_7"]
bz_pipe_traps = data_traps["Total_retention_volume_10"]

bz_pipes_traps = sum(
    [data_traps["Total_retention_volume_{}".format(vol)] for vol in volumes_bz_pipes]
)

# ## Plot


fig, axs = plt.subplots(2, 2, sharex=True)

plt.sca(axs[0, 0])
plt.plot(t_traps, lipb_traps, color=green_lipb)
plt.annotate(
    "LiPb",
    xy=(3, 2e21),
    color=green_lipb,
)
plt.ylabel(r"Inventory (T m$^{-1}$)")
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.sca(axs[0, 1])
plt.plot(t_traps, strucure_traps, color=grey_eurofer)
plt.annotate(
    "Structure",
    xy=(3, 1e20),
    color=grey_eurofer,
)
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.sca(axs[1, 0])
plt.plot(t_traps, bz_pipes_traps, color=grey_eurofer)
plt.annotate(
    "BZ pipes",
    xy=(3, 5e18),
    color=grey_eurofer,
)
plt.ylabel(r"Inventory (T m$^{-1}$)")
plt.xlabel(r"Time (days)")
plt.xlim(0, 6)
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)


plt.sca(axs[1, 1])
plt.plot(t_traps, first_wall_traps, color=red_W)
plt.annotate(
    "First wall",
    xy=(3, 4.2e17),
    color=red_W,
)
plt.xlabel(r"Time (days)")
plt.xlim(0, 6)
ax = plt.gca()
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# for ax in axs:
#     ax.spines["top"].set_visable("none")
#     ax.spines["right"].set_visable("none")


plt.show()
quit()


# fig = plt.figure(figsize=(4.5, 7))
fig = plt.figure()
linestyle_no_trap = "dashed"

# create one big plot to have a common y label
ax = fig.add_subplot(111)
ax.spines["top"].set_color("none")
ax.spines["bottom"].set_color("none")
ax.spines["left"].set_color("none")
ax.spines["right"].set_color("none")
ax.tick_params(labelcolor="w", top=False, bottom=False, left=False, right=False)
ax.set_ylabel(r"Inventory (T m$^{-1}$)")

# LiPb
ax1 = fig.add_subplot(2211)
plt.sca(ax1)

colour_lipb = green_lipb
plt.plot(t_traps, lipb_traps, color=colour_lipb)
plt.annotate(
    "LiPb",
    (t_traps[-1], lipb_traps[-1]),
    xytext=(t_traps[-1] * 1.1, lipb_traps[-1]),
    color=colour_lipb,
)


# Structure
ax2 = fig.add_subplot(2212)
plt.sca(ax2)
# colour_structure = 'grey'
colour_structure = grey_eurofer
plt.plot(t_traps, strucure_traps, color=colour_structure)
plt.annotate(
    "Structure",
    (t_traps[-1], strucure_traps[-1]),
    xytext=(t_traps[-1] * 1.1, strucure_traps[-1]),
    color=colour_structure,
)

# Pipes
ax3 = fig.add_subplot(2221)
plt.sca(ax3)
# colour_pipe11 = 'darkred'
colour_pipe11 = grey_eurofer
plt.plot(t_traps, bz_pipes_traps, color=colour_pipe11)
plt.annotate(
    "BZ pipes",
    (t_traps[-1], bz_pipes_traps[-1]),
    xytext=(t_traps[-1] * 1.1, bz_pipes_traps[-1]),
    color=colour_pipe11,
)


# FW
ax4 = fig.add_subplot(2222)
plt.sca(ax4)
# colour_first_wall = 'darkgoldenrod'
colour_first_wall = red_W
plt.plot(t_traps, first_wall_traps, color=colour_first_wall)
plt.annotate(
    "First wall",
    (t_traps[-1], first_wall_traps[-1]),
    xytext=(t_traps[-1] * 1.1, first_wall_traps[-1]),
    color=colour_first_wall,
)

plt.xlabel(r"Time (days)")

for ax in [ax1, ax2, ax3, ax4]:
    plt.sca(ax)
    ax.get_shared_x_axes().join(ax, ax4)
    # plt.ylabel(r"Inventory (T m$^{-1}$)")  # "Inventory \n" + r"(T m$^{-1}$)"
    plt.ylim(bottom=0)
    plt.xlim(0, 8)
    # plt.xlim(left=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# remove the xticks for top plots
ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax3.set_xticklabels([])

plt.tight_layout()

plt.show()
