import emcee
import matplotlib.pyplot as plt, numpy as np
import datetime
from craftroom.Talker import Talker
from tqdm import tqdm

class EnsembleSampler(emcee.EnsembleSampler, Talker):
    """Exactly the same as Dan F-M's ensemble sampler, but with some built-in plots."""
    def __init__(self, *args, **kwargs):
        emcee.EnsembleSampler.__init__(self, *args, **kwargs)
        Talker.__init__(self)
        self.nametag = 'emcee'
        self.labels = None

    def run_mcmc_with_progress(self, pos0, N, updates=10, **kwargs):
        '''Run the emcee sampler, but print an update every [updates] iterations.'''

        # calculate how long it take for one step
        before = datetime.datetime.now()
        self.speak('running MCMC, and will provide updates every {updates} steps'.format(**locals()))


        nchunks = np.int(np.ceil(np.float(N)/updates))
        pos = pos0
        count = 0
        for i in tqdm(range(nchunks)):

            pos, prob, state = self.run_mcmc(pos, updates)
            after = datetime.datetime.now()
            count += updates

            span = np.maximum((after - before).seconds, 1)/60.0
            rate = count/span
            remaining = (N-count)/rate
            #self.speak('completed {0}/{1} steps, in {2} minutes, with {3} minutes remaining'.format(count, N, span, remaining))

        return pos, prob, state

    def addLabels(self, labels):
        self.labels = labels

    def setupHistoryPlot(self, keys=None):

        if keys is None:
            self.toplot = self.labels
        else:
            self.toplot = keys
        self.indices = [self.labels.index(k) for k in self.toplot]
        self.nrows = len(self.toplot)

        # set up figure
        self.figure = plt.figure('emcee parameter history', figsize=(7,1*self.nrows+3))
        hr = np.ones(self.nrows+1)
        hr[-1] *=3
        self.gs = plt.matplotlib.gridspec.GridSpec(self.nrows+1, 1, hspace=0.1, wspace=0, height_ratios=hr, left=0.4, top=0.95, bottom=0.05, right=0.95)
        self.ax_history = {}
        ax = plt.subplot(self.gs[-1])
        ax.set_ylabel('lnp')
        self.ax_history['lnp'] = ax

        # loop over the individual parameters
        for i in range(self.nrows):
            key = self.toplot[i]
            ax = plt.subplot(self.gs[i], sharex=self.ax_history['lnp'])
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.setp(ax.get_yticklabels(), visible=False)
            ax.set_ylabel(key, rotation=0, ha='right')
            self.ax_history[key] = ax


    def HistoryPlot(self, limits, keys=None, maxwalkers=50):
        try:
            plt.scf(self.figure)
            self.ax_history
        except AttributeError:
            self.setupHistoryPlot(keys=keys)

        limits[1] = min(limits[1],self.chain.shape[1]-1)

        kw = dict(color='black', alpha=0.1, linewidth=1)
        x = np.arange(limits[0], limits[1]+1,1)
        for walker in range(np.minimum(self.k, maxwalkers)):
            self.ax_history['lnp'].plot(x,self.lnprobability[walker,limits[0]:limits[1]+1], **kw)
            for p in range(self.nrows):
                key = self.toplot[p]
                index = self.indices[p]
                self.ax_history[key].plot(x, self.chain[walker,limits[0]:limits[1]+1,index], **kw)
                if walker == 0:
                    self.ax_history[key].set_ylim(*np.percentile(self.flatchain[:,index], [1,99]))
        self.ax_history['lnp'].set_xlim(limits[0], limits[1])
        self.ax_history['lnp'].set_ylim(np.percentile(self.flatlnprobability,25), np.max(self.flatlnprobability))
        plt.draw()
