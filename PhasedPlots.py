from Plots import *

class PhasedPlots(Plot):


    def setup(self, tlcs=None, **kwargs):
        self.tlcs = tlcs
        self.telescopes = np.unique([tlc.telescope for tlc in tlcs])

        if len(self.telescopes) == 1:
            h = 1
            n = 1
        else:
            h = np.ones(len(self.telescopes)+1)
            h[-1] += 0#2
            n = len(self.telescopes) + 1
        # create a figure to hold all the phased light curves
        plt.figure('all light curves', figsize=(10,n*4), dpi=50)
        gs = plt.matplotlib.gridspec.GridSpec(n, 1, hspace=0.1, wspace=0, bottom=0.2)

        # set up empty dictionary of axes
        self.axes = {}

        # populate the axes that has the phased light curves
        if len(self.telescopes) == 1:
            self.axes['all'] = plt.subplot(111)
        else:
            self.axes['all'] = plt.subplot(gs[-1])
        count = 0
        for telescope in self.telescopes:
            if len(self.telescopes) == 1:
                self.axes[telescope] = self.axes['all']
            else:
                self.axes[telescope]= plt.subplot(gs[count], sharex=self.axes['all'], sharey=self.axes['all'])
                count += 1
                self.axes[telescope].yaxis.set_label_position("right")
                self.axes[telescope].set_ylabel('{0}'.format(tlc.telescope))
                self.axes[telescope].set_ylabel(telescope)


    def plot(self, tlcs, xlim=(-.1, 0.1), ylim=(0.985, 1.01), binsize=6.0/24.0/60.0, title=''):

        bjd, flux, uncertainty, whichtelescope, bigok = [],[],[],[], []

        count =0
        wholetel = {}
        for i in range(len(tlcs)):

            # pull out telescope light curve
            tlc = tlcs[i]
            ok = tlc.bad == False
            tm = tlc.TM
            tlc.plot(tm)
            bjd.extend(tlc.bjd)
            flux.extend(tlc.corrected())
            uncertainty.extend(tlc.uncertainty*tlc.rescaling)
            bigok.extend(ok)
            tels = np.empty(tlc.n).astype(np.str)
            tels[:] = tlc.name.replace(' ','')
            whichtelescope.extend(tels)

            telescope = tlc.telescope
            try:
                wholetel[telescope]
            except KeyError:
                wholetel[telescope] = dict(bjd=[],flux=[],uncertainty=[])


            wholetel[telescope]['bjd'].extend(tlc.bjd[ok])
            wholetel[telescope]['flux'].extend(tlc.corrected()[ok])
            wholetel[telescope]['uncertainty'].extend(tlc.uncertainty[ok]*tlc.rescaling)

            ax = self.axes[telescope]
            if len(self.telescopes) > 1:
                plt.sca(ax)
                tlc.plot(model=tm)
            plt.sca(self.axes['all'])
            tlc.plot(model=tm)

        if len(self.telescopes) > 1:
            for telescope in wholetel.keys():
                ax = self.axes[telescope]
                plt.sca(ax)

                for k in wholetel[telescope].keys():
                    wholetel[telescope][k] = np.array(wholetel[telescope][k])
                plotbinned(tm.planet.timefrommidtransit(wholetel[telescope]['bjd']), wholetel[telescope]['flux'], uncertainty=wholetel[telescope]['uncertainty'], bin=binsize)
                plt.setp(ax.get_xticklabels(), visible=False)

            #tlc.plot(model=tm, alpha=0.1)



        bjd = np.array(bjd)
        flux = np.array(flux)
        uncertainty = np.array(uncertainty)
        bigok = np.array(bigok)
        telescope = np.array(telescope)
        plt.sca(self.axes['all'])
        plotbinned(tm.planet.timefrommidtransit(bjd[bigok]), flux[bigok], uncertainty=uncertainty[bigok], alpha=0.75, bin=binsize)
        self.axes['all'].set_xlim(*xlim)
        self.axes['all'].set_ylim(*ylim)
        self.axes['all'].set_xlabel('Time from Mid-Transit (days)')
        self.axes['all'].set_ylabel('Relative Flux')
        self.axes['all'].set_title(title)
        table = astropy.table.Table(dict(bjd=bjd,
                                    flux=flux,
                                    uncertainty=uncertainty,
                                    telescope=whichtelescope))

        table.write('merged_lc.txt',
                    format='ascii.fixed_width',
                    bookend=False)
        tlc.setupSmooth()
        tm.plotPhased(linewidth=3,  alpha=0.5, color='black')
        for k in self.axes.keys():
            a = self.axes[k]
            plt.sca(a)
            tm.plotPhased(linewidth=3,  alpha=0.5, color='black')

        plt.savefig('combined_lightcurves.pdf')
