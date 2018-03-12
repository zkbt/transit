from .Plots import *

class MultiplexPlot(Talker):
    def __init__(self, tlcs, plotter, figsize=None, dpi=None, **kwargs):
        Talker.__init__(self)

        self.tlcs = tlcs
        # what are we going to populate each plot with?
        self.plotter = plotter

        # set up an overall grid, to match all the light curves
        self.figure = plt.figure('multiplexed', figsize=figsize, dpi=dpi)
        self.gs_overall = plt.matplotlib.gridspec.GridSpec(1,len(tlcs))
        for i, tlc in enumerate(tlcs):
            self.plotter(tlc=tlc, gs_destination=self.gs_overall[i], notfirst=i > 0, **kwargs)
            self.speak('after multiplex plot for {}'.format(tlc))
