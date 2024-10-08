from fenics import *
import festim as F


def slicer(filename):
    mesh_folder = "meshes/"
    mesh_3D = F.MeshFromXDMF(
        volume_file=mesh_folder + "mesh_domains_3D.xdmf",
        boundary_file=mesh_folder + "mesh_boundaries_3D.xdmf",
    )

    V_CG1 = FunctionSpace(mesh_3D.mesh, "CG", 1)
    temperature_field = filename
    T = Function(V_CG1)
    XDMFFile(temperature_field).read_checkpoint(T, "T", -1)

    class u_slice(UserExpression):
        def eval(self, value, x):

            value[0] = T(x[0], 0.255 / 2, x[2])

        def value_shape(self):
            return ()

    mesh_2D = F.MeshFromXDMF(
        volume_file=mesh_folder + "mesh_domains_2D.xdmf",
        boundary_file=mesh_folder + "mesh_boundaries_2D.xdmf",
    )

    V_2D = FunctionSpace(mesh_2D.mesh, "CG", 1)

    T.set_allow_extrapolation(True)
    V_2D = FunctionSpace(mesh_2D.mesh, "CG", 1)
    print("Projecting onto 2D mesh")
    T_sl = interpolate(u_slice(), V_2D)

    return T_sl


if __name__ == "__main__":
    temperature_field_3D = "Results/3D_results/T_3D.xdmf"

    T_sl = slicer(temperature_field_3D)

    results_folder = "Results/3D_results/"
    XDMFFile(results_folder + "T_slice.xdmf").write_checkpoint(
        T_sl, "T", 0, XDMFFile.Encoding.HDF5, append=False
    )
