'''
from .Instrument import Instrument
from .Parameter import Parameter
from .Parameters import Parameters
from .Planet import Planet
from .Star import Star
from .TLC import TLC
from .TM import TM
from .Synthesizer import Synthesizer, LM, MCMC
import Calculator
from .PhasedPlots import PhasedPlots
from .IndividualPlots import IndividualPlots
from .GroupedPlots import GroupedPlots
from .DiagnosticsPlots import DiagnosticsPlots
from .SmooshedPlot import SmooshedPlot
from .RVPhasedPlot import RVPhasedPlot
from .RVUnphasedPlot import RVUnphasedPlot
from .QuicklookTransitPlot import QuicklookTransitPlot
from .MultiplexPlot import MultiplexPlot
'''

__version__ = '0.0.0'

# specify whether we're calling this from .within setup.py
try:
    __TRANSITSETUP__
except NameError:
    __TRANSITSETUP__ = False

if not __TRANSITSETUP__:
    # (run this stuff if it's not form within setup.py)
    pass
