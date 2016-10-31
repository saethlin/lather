import time
import numpy as np
from matplotlib import pyplot as plt
import lather

obstime = np.linspace(0, 25.05, 50, endpoint=False)

sim = lather.Simulation('config.cfg')
results = sim.observe(obstime, True)

flux_observed = results['flux']
sim.clear_spots()
sim.add_spot(29., 180., 0.1, False)
sim.add_spot(-29., 0., 0.1, False)

print(time.time()-start)
results = sim.observe(obstime, observe_rv=True)
print(time.time()-start)

soap_flux = np.loadtxt('soap_flux.txt')
soap_rv = np.loadtxt('soap_rv.txt')

soap_rv -= soap_rv[0]
results['flux'] /= results['flux'][0]

f, axarr = plt.subplots(3, sharex=True)

axarr[0].set_title('LATHER vs SOAP-2')

axarr[0].plot(obstime, soap_rv*1000)
axarr[0].plot(obstime, results['rv']*1000)
axarr[0].set_ylabel('RV (m/s)')
axarr[0].set_xlim(0, 25.05)

axarr[1].plot(obstime, soap_rv*1000 - results['rv']*1000)
axarr[1].set_ylabel('Difference in RV (m/s)')
axarr[1].set_xlim(0, 25.05)

axarr[2].plot(obstime, soap_flux)
axarr[2].plot(obstime, results['flux'])
axarr[2].set_ylabel('Difference in relative flux')
axarr[2].set_xlabel('Time (days)')

fig = plt.gcf()
fig.set_size_inches(19*2/3, 10*2/3)
fig.tight_layout()
fig.savefig('comparison_plot.pdf')

plt.show()
