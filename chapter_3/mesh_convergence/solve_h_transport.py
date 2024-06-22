from fenics import *
from parameters_2D import my_model, id_lipb
from properties import S_0_lipb, E_S_lipb
import FESTIM as F
import csv


def h_transport_sim_FESTIM(
    mesh_domains_filename,
    mesh_boundaries_filename,
    temperature_field,
    velocity_field,
    case,
):
    """ """
    # define mesh to use
    my_model.mesh = F.MeshFromXDMF(
        volume_file=mesh_domains_filename,
        boundary_file=mesh_boundaries_filename,
    )

    # project temperature field on mesh
    my_model.T = F.TemperatureFromXDMF(filename=temperature_field, label="T")

    # idea to use function rather than filename
    # V_CG1 = FunctionSpace(my_model.mesh.mesh, "CG", 1)
    # my_model.T = F.Temperature()
    # my_model.T.T = project(temperature_field, V_CG1, solver_type="mumps")
    # my_model.T.create_functions(my_model.mesh)

    # detail results filename and location
    results_folder = "results/h_transport_results/{}/".format(case)
    for export in my_model.exports.exports:
        if isinstance(export, F.DerivedQuantities):
            export.filename = results_folder + "derived_quantities.csv"
        elif isinstance(export, F.XDMFExport):
            export.folder = results_folder
            export.append = False
            export.define_xdmf_file()

    my_model.initialise()

    # project velocity field
    u = velocity_field

    # modify the form F
    id_flow = id_lipb
    test_function_solute = my_model.h_transport_problem.mobile.test_function
    solute = my_model.h_transport_problem.mobile.solution

    dx = my_model.mesh.dx
    S_lipb = S_0_lipb * exp(-E_S_lipb / F.k_B / my_model.T.T)
    my_model.h_transport_problem.F += inner(
        dot(grad(S_lipb * solute), u), test_function_solute
    ) * dx(id_flow)

    my_model.run()

    number_of_cells = len(my_model.mesh.mesh.cells())
    simulation_time = my_model.timer.elapsed()[0]
    data = [number_of_cells, simulation_time]
    header = ["cells", "sim_time"]
    with open(results_folder + "sim_data.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(data)


if __name__ == "__main__":
    temperature_field_file = "results/standard_case/T_sl.xdmf"
    velocity_field = "results/standard_case/u_2D.xdmf"
    # project velocity field
    mesh_folder = "meshes/standard_case/"
    mesh_domains_file = mesh_folder + "mesh_domains_2D.xdmf"
    mesh_boundaries_file = mesh_folder + "mesh_boundaries_2D.xdmf"

    mesh = Mesh()
    XDMFFile(mesh_domains_file).read(mesh)
    V_ele = VectorElement("CG", mesh.ufl_cell(), 2)
    V_u = FunctionSpace(mesh, V_ele)
    u_full = Function(V_u, name="velocity")
    XDMFFile(velocity_field).read_checkpoint(u_full, "u", -1)

    h_transport_sim_FESTIM(
        mesh_domains_filename=mesh_domains_file,
        mesh_boundaries_filename=mesh_boundaries_file,
        temperature_field=temperature_field_file,
        velocity_field=u_full,
        case="standard_case",
    )
