import numpy as np
import lather

spots = np.loadtxt('g2013.txt', usecols=[2, 10, 11, 12], unpack=True)

spots[:, -1] /= 1e6  # Correct millionths of a stellar hemisphere to hemispheres

sim = lather.Simulation('config.cfg')
sim.clear_spots()

for spot, next_spot in zip(spots, spots[1:]):
    sim.add_spot(spot[4], spot[3], spot[2], False)

    if next_spot[0] != spot[0]:
        print(sim)
        exit()
        sim.observe(np.array([]))




