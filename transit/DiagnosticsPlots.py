from .Plots import *

class DiagnosticsPlots(Plot):

    def printParameters(self):
        #NEEDS WORKS!
                    todisplay = self.tlc.TM.lastfit.pdf.listParameters()

                    todisplay += '\n'
                    todisplay += 'chisq = {chisq:.1f}/{dof:.0f}\nrescaling = {rescaling:.2f}'.format(**self.tlc.TM.lastfit.notes)

                    self.ax_text.text(-0.1, 1.1, todisplay, va='top', fontsize=11)

    def setup(self, tlc=None, everything=True, gs_destination=None, correlations='fitted', notfirst=False, **kwargs):
        '''Setup the axes needed for plotting a light curve.'''

        self.tlc = tlc
        assert(self.tlc is not None)

        # create a figure to populate
        if everything:
            label = 'diagnostics'
        else:
            label = 'summary'

        if gs_destination is None:
            self.figure_diagnostics = plt.figure('light curve {0}'.format(label), figsize=(20,12), dpi=50)

        # set up the light curve plots
        if everything:
            # create two columns, to populate with lightcurves on left and diagnostics on right
            if gs_destination is None:
                gs_overarching = plt.matplotlib.gridspec.GridSpec(1, 3, width_ratios=[1, 0.3, 1.0], wspace=0.25, left=0.09, right=0.98)
            else:
                gs_overarching = plt.matplotlib.gridspec.GridSpecFromSubplotSpec(1, 3, gs_destination, width_ratios=[1, 0.3, 1.0], wspace=0.25)

            # create plots for light curves
            gs_lightcurve = plt.matplotlib.gridspec.GridSpecFromSubplotSpec(5, 2, hspace=0.05, wspace=0, width_ratios = [1,.2], height_ratios=[1,0.5,0.1,1,0.5], subplot_spec = gs_overarching[0])

            # noise properties
            gs_diagnostics = plt.matplotlib.gridspec.GridSpecFromSubplotSpec(3,1, height_ratios = [3,1,1], hspace=0.3, subplot_spec = gs_overarching[1])

            # potentially correlated variables
            self.to_correlate = []
            if correlations == 'fitted':
                for k in self.tlc.TM.floating:
                    if k != 'C':
                        if k in self.tlc.TM.instrument.__dict__.keys() and k != 'rescaling':
                            self.to_correlate.append(k.split('_tothe')[0])
            elif correlations == 'all':
                for k in self.tlc.externalvariables.keys():
                  if (k == 'ok') or (k == 'bad'):
                      # # a = raw_input('line 53')
                      continue
                  if np.std(self.tlc.externalvariables[k]) > 0:
                      self.to_correlate.append(k.split('_tothe')[0])
            ncor = np.int(np.ceil(np.sqrt(len(self.to_correlate))))
            gs_external = plt.matplotlib.gridspec.GridSpecFromSubplotSpec(len(self.to_correlate), 2, subplot_spec = gs_overarching[2], width_ratios=[1,5], hspace=0.1, wspace=0)
        else:
            # or, if desired, only plot the uncorrected and corrected light curves
            if gs_destination is None:
                gs_lightcurve = plt.matplotlib.gridspec.GridSpec(5, 1, hspace=0.05, wspace=0, height_ratios=[1,0.5,0.1,1,0.5])
            else:
                gs_lightcurve = plt.matplotlib.gridspec.GridSpecFromSubplotSpec(5, 2, gs_destination, hspace=0.05, wspace=0, width_ratios = [1,.2],  height_ratios=[1,0.5,0.1,1,0.5])
            self.ax_correlations = {}


        # set up the light curve (and residual) panels (leaving a space between the uncorrected and the corrected
        self.ax_raw = plt.subplot(gs_lightcurve[0,0])
        self.ax_instrument = plt.subplot(gs_lightcurve[1,0], sharex=self.ax_raw)

        self.ax_corrected = plt.subplot(gs_lightcurve[-2,0], sharex=self.ax_raw)
        self.ax_residuals = plt.subplot(gs_lightcurve[-1,0], sharex=self.ax_raw)
        self.ax_residuals_histogram = plt.subplot(gs_lightcurve[-1,1], sharey=self.ax_residuals)
        self.ax_instrument_histogram = plt.subplot(gs_lightcurve[1,1], sharey=self.ax_residuals, sharex=self.ax_residuals_histogram)

        # hide the tick labels on most of the light curve panels
        for a in [self.ax_raw, self.ax_corrected, self.ax_instrument]:
            plt.setp(a.get_xticklabels(), visible=False)
        for a in [self.ax_residuals_histogram, self.ax_instrument_histogram]:
              plt.setp(a.get_xticklabels(), visible=False)
              plt.setp(a.get_yticklabels(), visible=False)

        # set up the labels for the light curve panels
        try:
            self.ax_raw.set_title('E={epoch} | {left:.0f}-{right:.0f} nm'.format(name=self.tlc.name, epoch=self.tlc.epoch, left=self.tlc.left/10, right=self.tlc.right/10))
        except AttributeError:
            pass
        if notfirst == False:
            self.ax_raw.set_ylabel('Basic Photometry')
            self.ax_instrument.set_ylabel('transit\nresiduals\n(ppm)')
            self.ax_corrected.set_ylabel('Corrected Photometry')
            self.ax_residuals.set_ylabel('final\nresiduals\n(ppm)')
        self.ax_residuals.set_xlabel('Time from .Syzygy (days)')

        # set up the y limits (is this necessary?)
        #self.ax_raw.set_ylim(np.min(self.tlc.flux), np.max(self.tlc.flux))
        #self.ax_corrected.set_ylim(np.min(self.tlc.flux), np.max(self.tlc.flux))



        if everything:

            self.ax_text = plt.subplot(gs_diagnostics[0])
            self.ax_text.axis('off')

            # create a plot to store an autocorrelation function
            self.ax_acf = plt.subplot(gs_diagnostics[1])
            self.ax_acf.set_xlabel('Lag (in datapoints)')
            self.ax_acf.set_ylabel('ACF')

            # create a plot to show the RMS as a function of binning
            self.ax_binnedrms = plt.subplot(gs_diagnostics[2])
            self.ax_binnedrms.set_xlabel('# of datapoints in a bin')
            self.ax_binnedrms.set_ylabel('Binned RMS (ppm)')

            self.ax_correlations, self.ax_timeseries = {}, {}

            for i in range(len(self.to_correlate)):
                k = self.to_correlate[i]
                self.ax_correlations[k] = plt.subplot(gs_external[i,0])
                a = self.ax_correlations[k]
                plt.setp(a.get_xticklabels(), visible=False)
                plt.setp(a.get_yticklabels(), visible=False)
                a.set_ylabel(k, rotation=45, ha='right', size=10)

                self.ax_timeseries[k] = plt.subplot(gs_external[i,1], sharex=self.ax_residuals, sharey=self.ax_correlations[k])
                a = self.ax_timeseries[k]
                a.get_yaxis().get_major_formatter().set_useOffset(False)

                plt.setp(a.get_xticklabels(), visible=False)
                plt.setp(a.get_yticklabels(), visible=False)
                if i == len(self.to_correlate) - 1:
                    self.ax_correlations[k].set_xlabel('transit\nresiduals')

    def plot(self, noiseassumedforplotting=0.001, directory=None, ylim=[0.981, 1.007], binsize=6.0/60.0/24.0, mintimespan=None, phaseofinterest=0.0, **kwargs):
        '''A quick tool to plot what the light curve (and external variables) looks like.'''
        ppm = 1e6

        self.speak('plotting light curves with diagnostics')
        # create smoothed TLC structures, so the modeling will work
        self.tlc.TM.smooth_phased_tlc = self.tlc.fake( np.linspace(-self.tlc.TM.planet.period.value/2.0 + self.tlc.TM.planet.t0.value + 0.01, self.tlc.TM.planet.period.value/2.0 + self.tlc.TM.planet.t0.value-0.01, 100000))
        self.tlc.TM.smooth_unphased_tlc = self.tlc.fake(np.linspace(np.min(self.tlc.TLC.bjd), np.max(self.tlc.TLC.bjd), 100000))

        goodkw = {'color':self.tlc.colors['points'], 'alpha':0.5, 'marker':'o', 'edgecolor':'none', 's':10}
        badkw = {'color':self.tlc.colors['points'], 'alpha':0.25, 'marker':'x', 'edgecolor':'none', 's':10}


        time = self.tlc.TM.planet.timefrommidtransit(self.tlc.bjd)
        if phaseofinterest == 0.5:
            time[time < 0] += self.tlc.TM.planet.period.value
            #KLUDGE?

        notok = self.tlc.bad
        # # a = raw_input('early in the plotting for {0}?'.format(self.tlc))

        for good in [False,True]:

            if good:
                ok = (self.tlc.bad == 0)#.nonzero()
                kw = goodkw
                if np.sum(ok) == 0:
                    pass#return
            else:
                ok = (self.tlc.bad)#.nonzero()
                kw = badkw
            if np.sum(ok) == 0:
                # # a = raw_input('line 167')
                continue

            # # a = raw_input('good is {0}?'.format(good))

            self.ax_raw.scatter(time[ok], self.tlc.flux[ok], **kw)
            self.ax_corrected.scatter(time[ok], self.tlc.flux[ok]/self.tlc.TM.instrument_model()[ok], **kw)
            self.ax_residuals.scatter(time[ok], ppm*self.tlc.residuals()[ok], **kw)
            self.ax_instrument.scatter(time[ok], ppm*self.tlc.instrumental()[ok], **kw)
            # # a = raw_input('some scatters')

            if False:
                zorder = 10
                bkw = dict(zorder=10, markersize=0, elinewidth=3, linewidth=0, capthick=0, alpha=1 )
                plt.sca(self.ax_raw)
                plotbinned(time[ok], self.tlc.flux[ok],
                            uncertainty=self.tlc.uncertainty[ok], bin=binsize,
                            **bkw)

                plt.sca(self.ax_corrected)
                plotbinned(time[ok], self.tlc.flux[ok]/self.tlc.TM.instrument_model()[ok],
                            uncertainty=self.tlc.uncertainty[ok], bin=binsize,
                            **bkw)
                plt.sca(self.ax_residuals)
                plotbinned(time[ok], ppm*self.tlc.residuals()[ok],
                            uncertainty=ppm*self.tlc.uncertainty[ok],
                            bin=binsize, **bkw)
                plt.sca(self.ax_instrument)
                plotbinned(time[ok], ppm*self.tlc.instrumental()[ok],
                            uncertainty=ppm*self.tlc.uncertainty[ok], bin=binsize, **bkw)
                # # a = raw_input('lots of binned')
            plt.draw()

            self.tlc.points_correlations = {}
            kw['s'] = 10
            for k in self.ax_correlations.keys():
                #self.ax_correlations[k].plot(ppm*self.tlc.instrumental()[ok], self.tlc.externalvariables[k][ok], **kw)[0]
                #self.ax_timeseries[k].plot( time[ok], self.tlc.externalvariables[k][ok], **kw)

                # plot the external variables
                x = self.tlc.externalvariables[k]
                res, inst = self.tlc.residualsexceptfor(k, modeltoo=True)
                self.ax_correlations[k].scatter(ppm*(res[ok]), x[ok], **kw)
                if good:
                    binwidth = (np.nanmax(x[ok]) - np.nanmin(x[ok]))/10.0
                    bx, by, be = zachopy.oned.binto(x[ok], ppm*(res[ok]), yuncertainty=self.tlc.uncertainty[ok], binwidth=binwidth, robust=False, sem=True)
                    self.ax_correlations[k].errorbar(by, bx, None, be, color='black', alpha=0.3, elinewidth=3, marker='o', linewidth=0, capthick=3)
                    sorted = np.argsort(x[ok])
                    self.ax_correlations[k].plot(ppm*inst[ok][sorted], x[ok][sorted], color='gray', linewidth=3, alpha=0.75 )
                self.ax_timeseries[k].scatter( time[ok], self.tlc.externalvariables[k][ok], **kw)

                # set limits of plot windows
                if good:
                    self.ax_correlations[k].set_xlim(ppm*np.min(self.tlc.instrumental()[ok]), ppm*np.max(self.tlc.instrumental()[ok]))
                    self.ax_timeseries[k].set_ylim(np.min(self.tlc.externalvariables[k][ok]), np.max(self.tlc.externalvariables[k][ok]))

        # a = raw_input('done with loop')
        ok = (self.tlc.bad == 0).nonzero()
        which = np.median(np.round((self.tlc.TM.smooth_unphased_tlc.bjd - self.tlc.TM.planet.t0.value)/self.tlc.TM.planet.period.value))
        modeltime = self.tlc.TM.smooth_unphased_tlc.bjd - self.tlc.TM.planet.t0.value - which*self.tlc.TM.planet.period.value
        assert(len(modeltime) == len(self.tlc.TM.model(self.tlc.TM.smooth_unphased_tlc)))
        kw = {'color':self.tlc.colors['lines'], 'linewidth':3, 'alpha':1.0}

        self.tlc.line_raw = self.ax_raw.plot(modeltime, self.tlc.TM.model(self.tlc.TM.smooth_unphased_tlc), **kw)[0]
        self.tlc.line_corrected = self.ax_corrected.plot(modeltime, self.tlc.TM.planet_model(self.tlc.TM.smooth_unphased_tlc), **kw)[0]
        self.tlc.line_residuals = self.ax_residuals.plot(modeltime, ppm*np.zeros_like(modeltime), **kw)[0]
        justinstrument = (self.tlc.TM.instrument_model(self.tlc.TM.smooth_unphased_tlc)/np.median(self.tlc.TM.instrument_model()) - 1)
        self.tlc.line_instrument = self.ax_instrument.plot(modeltime, ppm*justinstrument, **kw)[0]


        # plot histograms of the residuals
        kw = {'color':self.tlc.colors['points'], 'linewidth':3, 'alpha':1.0}

        expectation = [0, ppm*np.mean(self.tlc.uncertainty)]
        zachopy.oned.plothistogram(ppm*self.tlc.instrumental()[ok], nbins=100, ax=self.ax_instrument_histogram,expectation =expectation , **kw)
        zachopy.oned.plothistogram(ppm*self.tlc.residuals()[ok], nbins=100, ax=self.ax_residuals_histogram, expectation =expectation , **kw)

        try:
            # plot binned RMS

            zachopy.oned.plotbinnedrms(self.tlc.residuals()[ok], ax=self.ax_binnedrms, yunits=1e6,  **kw)

            # plot the ACF
            zachopy.oned.plotautocorrelation(self.tlc.residuals()[ok], ax =self.ax_acf,  **kw)
        except AttributeError:
            pass

        # print text about the fit


        self.tlc.line_correlations = {}
        kw['alpha'] = 0.5
        kw['linewidth'] = 1
        #for k in self.ax_correlations.keys():
            #self.tlc.line_correlations[k] = self.ax_correlations[k].plot(ppm*justinstrument, self.tlc.TM.smooth_unphased_tlc.externalvariables[k],  **kw)[0]
        #assert(np.std(ppm*justinstrument)>1)

        nsigma = 5
        if noiseassumedforplotting is None:
            self.tlc.noiseassumedforplotting = np.mean(self.tlc.uncertainty)
        else:
            self.tlc.noiseassumedforplotting = noiseassumedforplotting
        if self.tlc.noiseassumedforplotting is not None:
            scale = self.tlc.noiseassumedforplotting*nsigma
        else:
            scale = nsigma*np.mean(self.tlc.uncertainty[ok])
        ppm = 1e6
        self.ax_residuals.set_ylim(-scale*ppm, scale*ppm)
        self.ax_instrument.set_ylim(-scale*ppm, scale*ppm)

        if mintimespan is None:
            self.ax_residuals.set_xlim(np.min(time), np.max(time))
        else:
            self.ax_residuals.set_xlim(np.mean(time)-mintimespan/2.0, np.mean(time)+mintimespan/2.0)

        buffer = 0.0075

        if phaseofinterest == 0.5:
            self.ax_corrected.set_ylim(*ylim)
            self.ax_corrected.get_yaxis().get_major_formatter().set_useOffset(False)
            self.ax_raw.set_ylim(*np.percentile(self.tlc.flux[self.tlc.ok], [1,99]))
        else:
            self.ax_corrected.set_ylim(1.0 - self.tlc.TM.planet.depth - buffer, 1.0 + buffer)
            self.ax_raw.set_ylim(1.0 - self.tlc.TM.planet.depth - buffer*3, 1.0 + buffer*3)

        plt.draw()
        if directory is None:
            directory = self.tlc.TM.directory

        #filename = directory + self.tlc.name.translate(None, '!@#$%^&*()<>') + '_lightcurveDiagnostics.pdf'
        #self.speak('saving light curve diagnostic plot to {0}'.format(filename))
        filename = self.tlc.directory + 'diagnostics.png'
        plt.savefig(filename, dpi=100)
        # a = raw_input('okay with the plot for {0}?'.format(self.tlc))
