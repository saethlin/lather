import numpy as np
import lather

spots = []
with open('g2013.txt') as spots_file:
    for line in spots_file:
        day = float(line[6:8])
        size = float(line[40:44])/1e6
        longitude = float(line[57:62])
        latitude = float(line[63:68])

        spots.append((day, (latitude, longitude, size, False)))

sim = lather.Simulation('config.cfg')
sim.clear_spots()

results = []

for spot, next_spot in zip(spots, spots[1:]):
    sim.add_spot(*spot[1])
    if spot[0] != next_spot[0]:
        results.append(sim.observe(np.array([spot[0]])))
        sim.clear_spots()