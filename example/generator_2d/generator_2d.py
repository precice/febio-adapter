import numpy as np
import time
import precice

# preCICE setup
interface = precice.Interface("Generator", "../precice-config.xml", 0, 1)
mesh_id = interface.get_mesh_id("Generator-Mesh")
data_id = interface.get_data_id("P_ext", mesh_id)

np_axis = 10  # Number of points along one axis
x_min, x_max = -0.5, 0.5
y_min, y_max = -0.5, 0.5
x_coords, y_coords = np.meshgrid(np.linspace(x_min, x_max, np_axis), np.linspace(y_min, y_max, np_axis))

dims = 2  # Dimensions of the generator mesh, even though the coupled problem is itself 3D
nv = np_axis ** dims  # Total number of vertices

coords = np.zeros((nv, dims + 1))  # Coordinates of vertices, to set the coupling mesh in preCICE (1 layer in Z axis)
write_data = np.zeros(nv)  # Data to write to preCICE in each time step

# Use same number of points in x and y axis because we want a square mesh
for y in range(np_axis):
    for x in range(np_axis):
        n = x + y * np_axis
        coords[n, 0] = x_coords[x, y]
        coords[n, 1] = y_coords[x, y]
        coords[n, 2] = 0.0

vertex_ids = interface.set_mesh_vertices(mesh_id, coords)

precice_dt = interface.initialize()

t = 0
dt = 0.01

while interface.is_coupling_ongoing():
    dt = np.minimum(dt, precice_dt)

    print("Generating data")
    write_data.fill(t*100)  # Change this value to change which constant data is written to preCICE

    interface.write_block_scalar_data(data_id, vertex_ids, write_data)

    precice_dt = interface.advance(dt)

    t = t + dt

interface.finalize()
