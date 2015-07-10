from Plots import *



class RVPhasedPlot(Plot):
    def __init__(self, **kwargs):
        self.label = 'phased'
        Plot.__init__(self, **kwargs)

    def setup(self, rvcs=None, **kwargs):

        # create a figure to hold all the phased light curves
        plt.figure(self.label, figsize=(8,5), dpi=100)
        gs = plt.matplotlib.gridspec.GridSpec(1, 1, hspace=0.1, wspace=0, bottom=0.15, left=0.2)

        # set up empty dictionary of axes
        self.axes = {}
        self.axes[self.label] = plt.subplot(gs[-1])
        self.axes[self.label].set_ylabel('Radial Velocity (m/s)')
        self.axes[self.label].set_xlabel('Time from Mid-transit (days)')

    def x(self, rvc):
        return rvc.TM.planet.timefrommidtransit(rvc.bjd)

    def y(self, rvc):
        return rvc.rv - rvc.TM.star.gamma.value

    def yerr(self, rvc):
        return rvc.uncertainty



    def plot(self, rvcs, xlim=(-.1, 0.1), ylim=(0.985, 1.01), binsize=6.0/24.0/60.0, title=''):

        for rvc in rvcs:
            x, y, yerr = self.x(rvc), self.y(rvc), self.yerr(rvc)


            kw = dict(linewidth=0, elinewidth=2, capthick=2, marker='o', markeredgecolor='none')

            plt.sca(self.axes[self.label])
            ink_errorbar(x, y, yerr, colors=rvc.colors['points'], alpha=1, **kw)


            ghosts = [-1,1]
            for ghost in ghosts:
                ink_errorbar(x+ ghost*rvc.TM.planet.period.value, y, yerr, colors=rvc.colors['points'], grayscale=True, alpha=0.5, **kw)

            tm = rvc.TM


            chunks = ghosts
            chunks.append(0)
            for ghost in ghosts:
                t = np.arange(tm.planet.period.value*(-0.5 + ghost), tm.planet.period.value*(0.5 + ghost), 0.01)
                x, y = t, tm.stellar_rv(t=t + tm.planet.t0.value) - tm.star.gamma.value
                sorted = np.argsort(x)
                x = x[sorted]
                y = y[sorted]

                if ghost == 0:
                    self.axes[self.label].plot(x, y, linewidth=3, alpha=0.5, color=rvc.colors['lines'])
                else:
                    self.axes[self.label].plot(x, y, linewidth=3, alpha=0.5, color='slategray')



        xlim = rvc.TM.planet.period.value*np.array([-0.75, 0.75])
        self.axes[self.label].set_xlim(*xlim)
