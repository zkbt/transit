import numpy as np

import matplotlib.pyplot as plt
import astropy.io
import zachopy.color, zachopy.utils
from zachopy.Talker import Talker
import scipy.interpolate
import george
import pemcee as emcee
import zachopy.borrowed.mpfit.mpfit as mpfit
import zachopy.oned
import time, datetime
ppm = 1e6
