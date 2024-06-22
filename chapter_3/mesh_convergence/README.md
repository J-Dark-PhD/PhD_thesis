# Mesh Refinement Studies

## 1. Mesh

To create a mesh, run:
```
C:\SALOME-9.8.0\run_salome.bat -t mesh_with_salome.py
```
You may need to replace the path to `run_salome.bat`.

This will create a .med file.

Then to convert the .med file in .xdmf files (readable by FESTIM), run:

```
python convert_mesh.py
```
If needed, install [meshio](https://github.com/nschloe/meshio) with:

```
pip install meshio[all]
```

## 2. Run FESTIM

Run a FEniCS container:

```
docker run -ti -v ${PWD}:/home/fenics/shared quay.io/fenicsproject/stable:latest
```

Install FESTIM v0.10.0:

```
pip install git+https://github.com/RemDelaporteMathurin/FESTIM@v0.10.0
```

To run the FESTIM simulation:

```
python main.py
```

This will produce several .xdmf files
