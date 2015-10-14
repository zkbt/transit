from Plots import Plot
from imports import *
from zachopy.painting import ink_errorbar
import zachopy.cmaps

scale = 1e3
class RVPhasedPlot(Plot):
    def __init__(self, **kwargs):
        self.label = 'phased'
        Plot.__init__(self, **kwargs)

    def setup(self, rvcs=None, **kwargs):

        inchinmm = 25.4#/2.0

        aspect = 0.6#2./3.
        # create a figure to hold all the phased light curves
        plt.figure(self.label, figsize=(89/inchinmm,89/inchinmm*aspect), dpi=100)
        gs = plt.matplotlib.gridspec.GridSpec(1, 1, hspace=0.1, wspace=0, bottom=0.18, left=0.16, right=0.91, top=0.95)

        # set up empty dictionary of axes
        self.axes = {}
        self.axes[self.label] = plt.subplot(gs[-1])
        self.axes[self.label].set_ylabel('Radial velocity (m/s)', fontsize=8)
        self.axes[self.label].set_xlabel('Phased time from mid-transit (days)', fontsize=8)

    def x(self, rvc):
        return rvc.TM.planet.timefrommidtransit(rvc.bjd)

    def y(self, rvc):
        return scale*(rvc.rv - rvc.TM.star.gamma.value)

    def yerr(self, rvc):
        return scale*rvc.uncertainty



    def plot(self, rvcs, xlim=(-.1, 0.1), ylim=(-15,15), binsize=6.0/24.0/60.0, title='', **kwargs):

        for rvc in rvcs:
            x, y, yerr = self.x(rvc), self.y(rvc), self.yerr(rvc)

            cmap = zachopy.cmaps.one2another(bottom='white', top='orangered')

            kw = dict(linewidth=0, elinewidth=1.5, capthick=1.5, marker='o', markersize=4, markeredgecolor='none')

            plt.sca(self.axes[self.label])

            minimumuncertainty = np.min(rvc.effective_uncertainty)
            weights = np.minimum((minimumuncertainty/rvc.effective_uncertainty)**2, 1)
            ink_errorbar(x, y, yerr, colors=cmap(weights), alpha=1, zorder=weights, **kw)
            print weights
            print cmap(weights)

            ghosts = [-1,1]
            for ghost in ghosts:
                bwcmap = zachopy.cmaps.one2another(bottom='white', top='silver')
                ink_errorbar(x+ ghost*rvc.TM.planet.period.value, y, yerr, colors=bwcmap(weights), zorder=weights, **kw)#grayscale=True,

            tm = rvc.TM


            chunks = ghosts
            chunks.append(0)
            for ghost in ghosts:
                t = np.arange(tm.planet.period.value*(-0.5 + ghost), tm.planet.period.value*(0.5 + ghost), 0.01)
                x, y = t, scale*(tm.stellar_rv(t=t + tm.planet.t0.value) - tm.star.gamma.value)
                sorted = np.argsort(x)
                x = x[sorted]
                y = y[sorted]

                if ghost == 0:
                    self.axes[self.label].plot(x, y, linewidth=1.5, alpha=1, color=cmap(0.5), zorder=0.5)
                else:
                    self.axes[self.label].plot(x, y, linewidth=1.5, alpha=1, color='lightgray', zorder=0.5)


        kw = dict(linestyle='--', linewidth=2, alpha=0.3, color='gray', zorder=-102)
        self.axes[self.label].axvline(tm.planet.period.value/2.0, **kw)
        self.axes[self.label].axvline(-tm.planet.period.value/2.0, **kw)
        xlim = rvc.TM.planet.period.value*np.array([-0.75, 0.75])
        self.axes[self.label].set_xlim(*xlim)
        self.axes[self.label].set_ylim(*ylim)
        ax = self.axes[self.label]
        plt.setp(ax.get_xticklabels(), fontsize=8)
        plt.setp(ax.get_yticklabels(), fontsize=8)
