import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib import cm, rcParams
from matplotlib.colors import LogNorm
import shutil
from mpl_toolkits.axes_grid1 import make_axes_locatable
import generate_data

rcParams["text.usetex"] = True if shutil.which("latex") else False


def plot_annealed_trap_fitting():
    """
    orginal data from A.Zaloznik et al, available at https://doi.org/10.1088/0031-8949/t167/1/014031
    """

    test_temperatures = [370, 400, 500, 600, 800]

    # read fitting data
    annealed_defect_type_1_densities = np.genfromtxt(
        "data/annealed_defect_1_densities.txt"
    )
    annealed_defect_type_2_densities = np.genfromtxt(
        "data/annealed_defect_2_densities.txt"
    )
    annealed_defect_type_3_densities = np.genfromtxt(
        "data/annealed_defect_3_densities.txt"
    )

    T_values = np.linspace(1, 900, num=1000)

    # ##### Plotting ##### #
    green_ryb = (117 / 255, 184 / 255, 42 / 255)
    firebrick = (181 / 255, 24 / 255, 32 / 255)
    electric_blue = (83 / 255, 244 / 255, 255 / 255)

    plt.figure()

    # defect type I
    plt.scatter(
        test_temperatures,
        generate_data.defect_type_1_densities,
        marker="x",
        color=green_ryb,
    )
    plt.plot(
        T_values,
        annealed_defect_type_1_densities,
        color=green_ryb,
        label=r"Defect type I",
    )

    # defect type II
    plt.scatter(
        test_temperatures,
        generate_data.defect_type_2_densities,
        marker="x",
        color=firebrick,
    )
    plt.plot(
        T_values,
        annealed_defect_type_2_densities,
        color=firebrick,
        label=r"Defect type II",
    )

    # defect type III
    plt.scatter(
        test_temperatures,
        generate_data.defect_type_3_densities,
        marker="x",
        color=electric_blue,
    )
    plt.plot(
        T_values,
        annealed_defect_type_3_densities,
        color=electric_blue,
        label=r"Defect type III",
    )

    plt.ylim(bottom=0)
    plt.xlim(350, 850)
    plt.ylabel(r"Trap density, n$_{\mathrm{t}}$ (at. \%)")
    plt.legend()
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.xlabel(r"Annealing temperature (K)")
    plt.tight_layout()

    plt.savefig("annealing_trap_fitting.pdf")

def plot_trap_density_distribution():
    x_0 = 2.3e-06
    dx_0 = 1.0e-07
    x_values = np.linspace(0, 1e-5, num=1000)

    # distrubution
    traps = 1 / (1 + np.exp((x_values - x_0) / dx_0))

    # ##### plotting ##### #
    plt.figure()
    plt.plot(x_values, traps, color="black")
    plt.xlim(0, 1e-05)
    plt.ylim(0, 1.01)
    plt.ylabel("Trap density ratio")
    plt.xlabel("x (m)")
    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.tight_layout()
    
    plt.savefig("trap_distribution.pdf")


def plot_TDS_data_fitted():
    """
    TDS data from T.Swartz-Selinger, currently unpublished
    """
    implantation_time = 72 * 3600
    resting_time = 0.5 * 24 * 3600
    sample_area = 12e-03 * 15e-03
    dpa_values = [2.5, 0.5, 0.23, 0.1, 0.023, 0.005, 0.001, 0]

    # ##### get tds and fitting data #### #
    tds_T_and_flux = []
    fitting_data = []

    for dpa in dpa_values:
        tds_data_file = "data/tds_data_schwartz_selinger/{}_dpa.csv".format(dpa)
        tds_data = np.genfromtxt(tds_data_file, delimiter=",", dtype=float)
        tds_data_T = tds_data[:, 0]
        tds_data_flux = tds_data[:, 1] / sample_area
        tds_T_and_flux.append([tds_data_T, tds_data_flux])

        results_folder = "data/damaged_sample_tds_fittings/"
        fitting_file = results_folder + "dpa_{}/last.csv".format(dpa)
        T_sim = []
        flux_1 = []
        flux_2 = []
        with open(fitting_file, "r") as csvfile:
            plots = csv.reader(csvfile, delimiter=",")
            for row in plots:
                if "t(s)" not in row:
                    if float(row[0]) >= implantation_time + resting_time * 0.75:
                        T_sim.append(float(row[1]))
                        flux_1.append(float(row[2]))
                        flux_2.append(float(row[3]))
        flux = -np.asarray(flux_1) - np.asarray(flux_2)
        fitting_data.append([T_sim, flux])

    # ##### plotting ##### #

    # colourbar
    plot_dpa_values = dpa_values[:-1]
    norm = LogNorm(vmin=min(plot_dpa_values), vmax=max(plot_dpa_values))
    colorbar = cm.viridis
    sm = plt.cm.ScalarMappable(cmap=colorbar, norm=norm)
    colours = [colorbar(norm(dpa)) for dpa in dpa_values]

    plt.figure(figsize=(6.4, 5.5))

    plt.plot(
        tds_T_and_flux[-1][0],
        tds_T_and_flux[-1][1],
        linewidth=3,
        color="black",
    )

    for case, colour in zip(tds_T_and_flux[:-1], colours):
        T_values = case[0]
        flux_values = case[1]
        plt.plot(T_values, flux_values, color=colour, linewidth=3)

    for case in fitting_data[:-1]:
        T_values = case[0]
        flux_values = case[1]
        plt.plot(
            T_values,
            flux_values,
            color="grey",
            linewidth=2,
            linestyle="dashed",
            alpha=0.5,
        )
    plt.plot(
        fitting_data[-1][0],
        fitting_data[-1][1],
        color="grey",
        linewidth=2,
        linestyle="dashed",
        alpha=0.5,
        label=r"Simulation",
    )

    plt.xlim(300, 1000)
    plt.ylim(0, 1e17)
    plt.xlabel(r"Temperature (K)")
    plt.ylabel(r"Desorption flux (D m$ ^{-2}$ s$ ^{-1}$)")
    plt.legend()
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.annotate("undamaged", xy=(320, 3e16), color="black")
    ax.arrow(370, 2.8e16, 130, -2e16, color="black", zorder=10, width=0.1)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    plt.colorbar(sm, cax=cax, label=r"Damage (dpa)")
    
    plt.savefig("tds_data.pdf")
    

def plot_0_1_dpa_TDS_data_fitted():
    green_ryb = (117 / 255, 184 / 255, 42 / 255)
    firebrick = (181 / 255, 24 / 255, 32 / 255)
    pewter_blue = (113 / 255, 162 / 255, 182 / 255)
    blue_jeans = (96 / 255, 178 / 255, 229 / 255)
    electric_blue = (83 / 255, 244 / 255, 255 / 255)

    implantation_time = 72 * 3600
    resting_time = 0.5 * 24 * 3600

    tds_data_file = "data/tds_data_schwartz_selinger/0.1_dpa.csv"
    tds_data = np.genfromtxt(tds_data_file, delimiter=",")

    temperatures = tds_data[:, 0]

    area = 12e-03 * 15e-03
    flux = tds_data[:, 1] / area

    plt.figure()
    plt.scatter(temperatures, flux, label="Exp", color="black")

    T_sim = []
    flux1, flux2 = [], []
    solute = []
    ret = []
    t = []
    trap1 = []
    trap2 = []
    trap3 = []
    trap4 = []
    trap5 = []
    trap6 = []

    with open("data/damaged_sample_tds_fittings/dpa_0.1/last.csv", "r") as csvfile:
        plots = csv.reader(csvfile, delimiter=",")
        for row in plots:
            if "t(s)" not in row:
                if float(row[0]) >= implantation_time + resting_time * 0.75:
                    t.append(float(row[0]))
                    T_sim.append(float(row[1]))
                    flux1.append(float(row[2]))
                    flux2.append(float(row[3]))
                    solute.append(float(row[4]))
                    ret.append(float(row[5]))
                    trap1.append(float(row[6]))
                    trap2.append(float(row[7]))
                    trap3.append(float(row[8]))
                    trap4.append(float(row[9]))
                    trap5.append(float(row[10]))
                    trap6.append(float(row[11]))


    trap_1_contrib = (np.diff(trap1) / np.diff(t)) * -1
    trap_2_contrib = (np.diff(trap2) / np.diff(t)) * -1
    trap_3_contrib = (np.diff(trap3) / np.diff(t)) * -1
    trap_4_contrib = (np.diff(trap4) / np.diff(t)) * -1
    trap_5_contrib = (np.diff(trap5) / np.diff(t)) * -1
    trap_6_contrib = (np.diff(trap6) / np.diff(t)) * -1
    
    plt.plot(
        T_sim,
        -np.asarray(flux1) - np.asarray(flux2),
        label="Simulation",
        color="orange",
        linewidth=2,
    )

    plt.plot(
        T_sim[1:],
        trap_1_contrib,
        linestyle="dashed",
        color=green_ryb,
        label=r"Trap 1",
    )
    plt.plot(
        T_sim[1:],
        trap_2_contrib,
        linestyle="dashed",
        color=firebrick,
        label=r"Trap D1",
    )
    plt.plot(
        T_sim[1:],
        trap_3_contrib,
        linestyle="dashed",
        color=blue_jeans,
        label=r"Trap D2",
    )
    plt.plot(
        T_sim[1:],
        trap_4_contrib,
        linestyle="dashed",
        color=electric_blue,
        label=r"Trap D3",
    )
    plt.plot(
        T_sim[1:],
        trap_5_contrib,
        linestyle="dashed",
        color=pewter_blue,
        label=r"Trap D4",
    )
    plt.plot(
        T_sim[1:],
        trap_6_contrib,
        linestyle="dashed",
        color="black",
        label=r"Trap D5",
    )
    plt.fill_between(T_sim[1:], 0, trap_1_contrib, color="grey", alpha=0.1)
    plt.fill_between(T_sim[1:], 0, trap_2_contrib, color="grey", alpha=0.1)
    plt.fill_between(T_sim[1:], 0, trap_3_contrib, color="grey", alpha=0.1)
    plt.fill_between(T_sim[1:], 0, trap_4_contrib, color="grey", alpha=0.1)
    plt.fill_between(T_sim[1:], 0, trap_5_contrib, color="grey", alpha=0.1)
    plt.fill_between(T_sim[1:], 0, trap_6_contrib, color="grey", alpha=0.1)



    plt.xlim(300, 1000)
    plt.ylim(bottom=0)
    plt.xlabel(r"Temperature (K)")
    plt.ylabel(r"Desorption flux (D m$ ^{-2}$ s$ ^{-1}$)")
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.legend()
    plt.tight_layout()
    
    plt.savefig("tds_fitted_0.1_dpa.pdf")
  

def plot_damaged_trap_fitting():
    """
    TDS data from T.Swartz-Selinger, currently unpublished
    """
    dpa_values = [0, 0.001, 0.005, 0.023, 0.1, 0.23, 0.5, 2.5]

    # read fitting data
    trap_D1_fitting = np.genfromtxt("data/damage_trap_D1_fitting.txt")
    trap_D2_fitting = np.genfromtxt("data/damage_trap_D2_fitting.txt")
    trap_D3_fitting = np.genfromtxt("data/damage_trap_D3_fitting.txt")
    trap_D4_fitting = np.genfromtxt("data/damage_trap_D4_fitting.txt")
    trap_D5_fitting = np.genfromtxt("data/damage_trap_D5_fitting.txt")
    dpa_x_values = np.linspace(0, 3, num=1000)

    # ##### plotting ##### #

    green_ryb = (117 / 255, 184 / 255, 42 / 255)
    firebrick = (181 / 255, 24 / 255, 32 / 255)
    pewter_blue = (113 / 255, 162 / 255, 182 / 255)
    electric_blue = (83 / 255, 244 / 255, 255 / 255)

    fig = plt.figure()

    plt.scatter(
        dpa_values, generate_data.trap_D1_densities, marker="x", color=firebrick
    )
    plt.plot(dpa_x_values, trap_D1_fitting, color=firebrick, label=r"Trap D1")
    plt.scatter(
        dpa_values, generate_data.trap_D2_densities, color=pewter_blue, marker="x"
    )
    plt.plot(dpa_x_values, trap_D2_fitting, color=pewter_blue, label=r"Trap D2")
    plt.scatter(
        dpa_values, generate_data.trap_D3_densities, color=electric_blue, marker="x"
    )
    plt.plot(dpa_x_values, trap_D3_fitting, color=electric_blue, label=r"Trap D3")
    plt.scatter(
        dpa_values, generate_data.trap_D4_densities, color=green_ryb, marker="x"
    )
    plt.plot(dpa_x_values, trap_D4_fitting, color=green_ryb, label=r"Trap D4")
    plt.scatter(dpa_values, generate_data.trap_D5_densities, color="black", marker="x")
    plt.plot(dpa_x_values, trap_D5_fitting, color="black", label=r"Trap D5")
    
    plt.legend(loc="lower right")
    plt.xlabel(r"Damage (dpa)")
    plt.ylabel(r"Trap density (m$^{-3}$)")
    plt.ylim(0, 7e25)
    plt.xlim(left=0)

    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    
    plt.savefig("damaged_trap_fitting.pdf")

if __name__ == "__main__":
    plot_annealed_trap_fitting()
    plot_trap_density_distribution()
    plot_TDS_data_fitted()
    plot_0_1_dpa_TDS_data_fitted()
    plot_damaged_trap_fitting()

    # plt.show()
