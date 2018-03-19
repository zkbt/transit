from .Plots import *

class QuicklookTransitPlot(Plot):
    '''Plot one TLC, with its currently set model (potentially sampling from .a GP).'''

    def __init__(self, **kwargs):
        Plot.__init__(self, **kwargs)

    def setup(self, tlc=None, **kwargs):

        # store the TLC
        self.tlc = tlc
        try:
            for k,v in self.tlc.quicklook.iteritems():
                self.__dict__[k] = v
            return
        except AttributeError:
            tlc.quicklook = self

        self.figure = plt.figure(self.tlc.name, figsize=(15,20))
        gs = plt.matplotlib.gridspec.GridSpec(5,2,
                    height_ratios=[1,1,1,1,1],width_ratios=[1, .5],
                    wspace=0,hspace=0.1,left=0.15)

        # create the axes
        self.ax = {}
        self.ax['raw'] = plt.subplot(gs[0,0])
        for i, k in enumerate(['cleaned',  'residuals', 'cleanedfromgp', 'residualsfromgp']):
            self.ax[k] = plt.subplot(gs[i+1,0], sharex=self.ax['raw'])

        self.ax['text'] = plt.subplot(gs[:,1])
        self.ax['text'].axis('off')

        for k, a in self.ax.iteritems():
            if k != 'residualsfromgp':
                plt.setp(a.get_xticklabels(), visible=False)

        self.ax['raw'].set_ylabel('Raw')
        self.ax['cleaned'].set_ylabel('Without\nInstrument')
        self.ax['cleanedfromgp'].set_ylabel('Without\nGP')
        self.ax['residuals'].set_ylabel('O-C')
        self.ax['residualsfromgp'].set_ylabel('O-C-GP')
        self.ax['residualsfromgp'].set_xlabel('Time from .Mid-transit (days)')

        t = tlc.TM.planet.timefrommidtransit(tlc.bjd)

        self.ax['raw'].set_xlim(np.min(t), np.max(t))
        self.smoothed = tlc.fake(np.linspace(np.min(tlc.bjd), np.max(tlc.bjd), 1000))

    def plot(self, **kwargs):
        tlc = self.tlc
        self.points, self.lines = {}, {}

        ok = tlc.ok
        t = tlc.TM.planet.timefrommidtransit(tlc.bjd)
        t_smooth = tlc.TM.planet.timefrommidtransit(self.smoothed.bjd)

        color = tlc.color
        pkw = dict(color=color, markeredgecolor=color, linewidth=0, marker='o', markersize=5, alpha=0.5, )
        lkw = dict(color='gray', linewidth=2,  alpha=0.5, )

        # have to do this to make sure we recomputed the gp
        lnp = tlc.gp_lnprob()

        mean_wiggle, cov = tlc.gp.predict(tlc.residuals()[tlc.ok], self.smoothed.bjd)
        i = np.arange(len(mean_wiggle))
        std = np.sqrt(cov[i,i])



        k = 'raw'
        self.points[k] = self.ax[k].plot(t[ok],tlc.flux[ok],**pkw)

        k = 'cleaned'
        self.points[k] = self.ax[k].plot(t[ok],tlc.corrected()[ok],**pkw)

        #self.lines[k] = self.ax[k].plot()
        k = 'residuals'
        self.points[k] = self.ax[k].plot(t[ok],tlc.residuals()[ok],**pkw)
        #self.lines[k] = self.ax[k].plot(t_smooth, mean+std,**lkw)
        #self.lines[k] = self.ax[k].plot(t_smooth, mean-std,**lkw)

        n = 10
        for i in range(n+1):
            wiggle = tlc.gp.sample_conditional(tlc.residuals()[tlc.ok], self.smoothed.bjd)
            if i == n:
                lkw['color'] = 'black'
                wiggle = mean_wiggle

            k = 'residuals'
            self.lines[k] = self.ax[k].plot(t_smooth, wiggle ,**lkw)
            k = 'cleaned'
            self.lines[k] = self.ax[k].plot(t_smooth, wiggle + tlc.TM.planet_model(tlc=self.smoothed) ,**lkw)
            k = 'raw'
            self.lines[k] = self.ax[k].plot(t_smooth, wiggle + tlc.TM.model(tlc=self.smoothed) ,**lkw)


        mean_wiggle, cov = tlc.gp.predict(tlc.residuals()[tlc.ok], tlc.bjd[tlc.ok])
        k = 'cleanedfromgp'
        self.points[k] = self.ax[k].plot(t[ok],tlc.corrected()[ok] - mean_wiggle,**pkw)
        self.lines[k] = self.ax[k].plot(t_smooth, tlc.TM.planet_model(tlc=self.smoothed) ,**lkw)

        k = 'residualsfromgp'
        self.points[k] = self.ax[k].plot(t[ok],tlc.residuals()[ok]-mean_wiggle,**pkw)
        self.lines[k] = self.ax[k].plot(t[ok], np.zeros_like(t[ok]),**lkw)

        self.printParameters()
        plt.draw()


    def printParameters(self):
        toprint = '\n'*5
        toprint += self.tlc.name + '\n'*2
        for p in self.tlc.TM.parameters:
            toprint += '{0} = {1}\n'.format(p.name, p.value)

        self.ax['text'].text(0.1, 0.9, toprint, va='top', fontsize=11)
'''# create smoothed TLC structures, so the modeling will work
self.TM.smooth_phased_tlc = self.fake( np.linspace(-self.TM.planet.period.value/2.0 + self.TM.planet.t0.value + 0.01, self.TM.planet.period.value/2.0 + self.TM.planet.t0.value-0.01, 10000))
self.TM.smooth_unphased_tlc = self.fake(np.linspace(np.min(self.TLC.bjd), np.max(self.TLC.bjd), 10000))



t_phased = self.TM.planet.timefrommidtransit(self.bjd)
t_unphased = self.bjd - self.TM.planet.t0.value

try:
assert(self.ready)
except:
self.points_phased = self.ax_phased.plot(t_phased, self.flux, **kw)
self.points_phased_zoom = self.ax_phased_zoom.plot(t_phased, self.flux, **kw)
self.points_unphased = self.ax_unphased.plot(t_unphased, self.flux, **kw)
self.points_unphased_zoom = self.ax_unphased_zoom.plot(t_unphased, self.flux, **kw)
self.ready = True

for phased in [self.points_phased[0], self.points_phased_zoom[0]]:
phased.set_data(t_phased, self.flux)
for unphased in [self.points_unphased[0], self.points_unphased_zoom[0]]:
unphased.set_data(t_unphased, self.flux)
self.TM.plot()

nsigma=5
self.ax_phased.set_ylim(np.min(self.flux)-nsigma*np.mean(self.uncertainty), np.max(self.flux)+nsigma*np.mean(self.uncertainty))
self.ax_unphased.set_xlim(np.min(t_unphased), np.max(t_unphased))
self.ax_phased.set_xlim(-self.TM.planet.period.value/2.0, self.TM.planet.period.value/2.0)
self.ax_phased_zoom.set_xlim(-self.TM.planet.duration, self.TM.planet.duration)
self.ax_unphased_zoom.set_xlim(-self.TM.planet.duration, self.TM.planet.duration)

plt.draw()'''
