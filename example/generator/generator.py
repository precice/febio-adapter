import numpy as np
import time
import precice


n = 20
dn = 1 / n

# generate mesh
y = np.linspace(0, 1, n + 1)

# preCICE setup
participant_name = "Generator"
config_file_name = "../precice-config.xml"
solver_process_index = 0
solver_process_size = 1
interface = precice.Interface(participant_name, config_file_name, solver_process_index, solver_process_size)

interface.initialize()

t = 0

while True:

    dt = 0.01

    print("Generating data")
    time.sleep(0.2)
    u = 1 - 2 * np.random.rand(n)

    t = t + dt
    if t > 0.1:
        break


interface.finalize()
