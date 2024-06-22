import meshio

mesh = meshio.read("WCLL_2D.med")

# In order to use MeshFunction of FEniCS
# The tag must be a positive number (size_t)

# print(mesh.cell_data_dict["cell_tags"].keys())

# this doesn't work cause cell_data_dict doesn't have a setter
mesh.cell_data_dict["cell_tags"]["triangle"] *= -1
mesh.cell_data_dict["cell_tags"]["line"] *= -1

# print the correspondance
print("This is the correspondance dict")
print(mesh.cell_tags)
# mesh.cell_tags = {-6: ['Down'], -7: ['Top'], -8: ['Lying on Top']}


# Export mesh that contains only triangular faces
# along with tags
meshio.write_points_cells(
    "mesh_domains_2D.xdmf",
    mesh.points,
    [mesh.cells[1]],
    cell_data={"f": [-1 * mesh.cell_data["cell_tags"][1]]},
)


# Export mesh that contains only lines
# along with tags

meshio.write_points_cells(
    "mesh_boundaries_2D.xdmf",
    mesh.points,
    [mesh.cells[0]],
    cell_data={"f": [-1 * mesh.cell_data["cell_tags"][0]]},
)
