import emcee
import matplotlib.pyplot as plt, numpy as np
import datetime
from zachopy.Talker import Talker

class EnsembleSampler(emcee.EnsembleSampler, Talker):
    """Exactly the same as Dan F-M's ensemble sampler, but with some built-in plots."""
    def __init__(self, *args, **kwargs):
        emcee.EnsembleSampler.__init__(self, *args, **kwargs)
        Talker.__init__(self)
        self.nametag = 'emcee'
        self.labels = None

    def run_mcmc_with_progress(self, pos0, N, updates=15, **kwargs):
        '''Run the emcee sampler, but print an update every [updates] seconds.'''

        # calculate how long it take for one step
        before = datetime.datetime.now()
        self.run_mcmc(pos0,1)
        after = datetime.datetime.now()
        dt = after - before
        timeperstep = dt.microseconds/1.0e6
        self.speak('each step of {nwalkers} walkers takes {timeperstep:0.6f}s'.format(nwalkers=self.k, timeperstep=timeperstep))

        # calculate of
        nperupdate = np.int(max(np.float(updates)/timeperstep, 1))
        nchunks = np.int(np.ceil(np.float(N)/nperupdate))

        pos = pos0
        for i in range(nchunks):
            pos, prob, state = self.run_mcmc(pos, nperupdate)
            self.speak('{0:10.0f}/{1} steps, {2:2.2f} minutes remaining'.format(i*nperupdate, N, (nchunks-i-1)*nperupdate*timeperstep/60.0))

        return pos, prob, state

    def addLabels(self, labels):
        self.labels = labels

    def setupHistoryPlot(self):

        self.figure = plt.figure('emcee parameter history')
        hr = np.ones(self.dim+1)
        hr[-1] *=3
        self.gs = plt.matplotlib.gridspec.GridSpec(self.dim+1, 1, hspace=0.1, wspace=0, height_ratios=hr, left=0.3)
        self.ax_history = {}
        ax = plt.subplot(self.gs[-1])
        ax.set_ylabel('lnp')
        self.ax_history['lnp'] = ax

        # loop over the individual parameters
        for i in range(self.dim):
            key = self.labels[i]
            ax = plt.subplot(self.gs[i], sharex=self.ax_history['lnp'])
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.setp(ax.get_yticklabels(), visible=False)
            ax.set_ylabel(key, rotation=0, ha='right')
            self.ax_history[key] = ax


    def HistoryPlot(self, limits):
        try:
            self.ax_history
        except:
            self.setupHistoryPlot()

        limits[1] = min(limits[1],self.chain.shape[1]-1)
        kw = dict(color='black', alpha=0.1, linewidth=1)
        x = np.arange(limits[0], limits[1],1)
        for walker in range(self.k):
            self.ax_history['lnp'].plot(x,self.lnprobability[walker,limits[0]:limits[1]], **kw)
            for p in range(self.dim):
                key = self.labels[p]
                self.ax_history[key].plot(x, self.chain[walker,limits[0]:limits[1],p], **kw)
                if walker == 0:
                    self.ax_history[key].set_ylim(*np.percentile(self.flatchain[:,p], [1,99]))
        self.ax_history['lnp'].set_xlim(0, limits[1])
        self.ax_history['lnp'].set_ylim(np.percentile(self.flatlnprobability,25), np.max(self.flatlnprobability))
        plt.draw()
