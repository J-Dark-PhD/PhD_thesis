from fenics import *
import FESTIM as F


def slicer(
    temperature_field_3D, mesh_domains_filename, mesh_boundaries_filename, case, exports
):
    mesh_folder = "meshes/standard_case/"
    mesh_3D = F.MeshFromXDMF(
        volume_file=mesh_folder + "mesh_domains_3D.xdmf",
        boundary_file=mesh_folder + "mesh_boundaries_3D.xdmf",
    )

    V_CG1 = FunctionSpace(mesh_3D.mesh, "CG", 1)
    temp_field = temperature_field_3D
    T = Function(V_CG1)
    XDMFFile(temp_field).read_checkpoint(T, "T", -1)

    class u_slice(UserExpression):
        def eval(self, value, x):

            value[0] = T(x[0], x[1], 0.116 / 2)

        def value_shape(self):
            return ()

    mesh_2D = F.MeshFromXDMF(
        volume_file=mesh_domains_filename,
        boundary_file=mesh_boundaries_filename,
    )
    V_2D = FunctionSpace(mesh_2D.mesh, "CG", 1)

    T.set_allow_extrapolation(True)
    V_2D = FunctionSpace(mesh_2D.mesh, "CG", 1)
    print("Projecting onto 2D mesh")
    T_sl = interpolate(u_slice(), V_2D)
    if exports:
        XDMFFile("results/temperature_fields/T_{}.xdmf".format(case)).write_checkpoint(
            T_sl, "T", 0, XDMFFile.Encoding.HDF5, append=False
        )
    return T_sl


if __name__ == "__main__":
    mesh_folder = "meshes/standard_case/"
    mesh_domains_filename = mesh_folder + "mesh_domains_2D.xdmf"
    mesh_boundaries_filename = mesh_folder + "mesh_boundaries_2D.xdmf"
    temperature_field_3D = "results/3D_results/T.xdmf"
    T_sl = slicer(
        temperature_field_3D,
        mesh_domains_filename,
        mesh_boundaries_filename,
        case="sl",
        exports=False,
    )
    XDMFFile("results/standard_case/T_sl.xdmf").write_checkpoint(
        T_sl, "T", 0, XDMFFile.Encoding.HDF5, append=False
    )
