import numpy as np
from matplotlib import pyplot as plt
import scipy.signal
import lather
import matplotlib

def anomaly(time, period, ecc, argument_of_periastron):
    """
    Loosely based on Wright and Howard's calcnu function
    """
    ecc %= 1
    phase = (time - argument_of_periastron) / period
    mean_anomaly = 2 * np.pi * (phase - np.floor(phase))

    # Newton's method
    ecc_anomaly = mean_anomaly.copy()
    while True:
        diff = ecc_anomaly - ecc * np.sin(ecc_anomaly) - mean_anomaly
        ecc_anomaly -= diff / (1 - ecc * np.cos(ecc_anomaly))
        if np.all(abs(diff) < 1e-8):
            break

    return 2 * np.arctan(np.sqrt((1 + ecc) / (1 - ecc)) * np.tan(ecc_anomaly / 2))


def rv_model(time, period, mass, eccentricity, w, t_p):
    """
    Expects orbit parameters Period, semiamplitude, eccentricity,
    argument of periastron, t_p
    """
    model = np.zeros(time.size)
    if len(args) > 5:
        y = args[-1]
        orbit_parameters = np.array(args[:-1]).reshape((-1, 5))
        model += y
    else:
        orbit_parameters = np.array(args).reshape((-1, 5))

    for P, M, e, w, t_p in orbit_parameters:
        K = (2*np.pi*6.67e-11/P)**(1/3) * (M*6e24/(2e30*0.5 + 6e24*M)**(2/3)) * 1/(np.sqrt(1-e**2))
        f = anomaly(time, P, e, t_p)
        model += K * np.cos(w + f) + e * np.cos(w)

    return model


def periodogram(time, data, start, stop, N=1e4):
    """
    Use scipy's lombscargle algorithm, faster on small data sets and produces
    same results as matlab's lomb()
    """

    freq = 1/(10**np.linspace(np.log10(start), np.log10(stop), N))

    nfactor = 2/(data.size * np.var(data))

    power = scipy.signal.lombscargle(time, data-np.mean(data), 2*np.pi*freq) * nfactor

    return 1/freq, power


sim = lather.Simulation('config.cfg')
obstime = np.arange(100).astype(np.float64)
obstime += 1/24*np.random.rand(obstime.size)
spot_signal = sim.observe(obstime, observe_rv=True)['rv'] * 1000

total = spot_signal.copy()

# Add planet and observation noise
planet_signal = rv_model(obstime, 35, 1.3, 0.0, 0.0, 0.0)
total += planet_signal
total += np.random.rand(total.size)

periods, power = periodogram(obstime, total, 2.0, 150)


fig, axarr = plt.subplots(2)

axarr[0].set_title('Individual signals')

axarr[0].plot(obstime, total, 'ko')
axarr[0].plot(obstime, planet_signal)
axarr[0].plot(obstime, spot_signal)
axarr[0].set_ylabel('RV (m/s)')
axarr[0].set_xlabel('Time (d)')

axarr[1].set_title('Total periodogram')
axarr[1].semilogx(periods, power)
axarr[1].set_xlim(periods.min(), periods.max())
axarr[1].set_ylabel('Periodogram power')
axarr[1].set_xlabel('Period (d)')

plt.show()

fig.set_size_inches(10, 8)
fig.tight_layout()
fig.savefig('detection_test.pdf')

