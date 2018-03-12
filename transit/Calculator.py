from .imports import *
from craftroom.units import *

def K( planet_mass=None, stellar_mass=None, period=None):
    '''Calculate K, in cm/s, given:
            planet_mass in Earths
            stellar_mass in Suns
            period in days
    '''


    # convert units
    stellar_mass = self.stellar_mass*Msun
    period = self.period*day

    # everything below here is in cgs units
    q = (planet_mass*Mearth)/(stellar_mass*Msun)
    return q*(2*np.pi*G*(stellar_mass*Msun)/(period*day))**(1.0/3.0)

def aoverrs( period=None, stellar_mass=None, stellar_radius=None):
    '''Calculate aoverrs (unitless), given:
            period in days
            stellar_mass in Earths
            stellar_radius in Earths
    '''

    return (G*(period*day)**2*(stellar_mass*Msun)/4/np.pi**2/(stellar_radius*Rsun)**3)**(1./3.)
