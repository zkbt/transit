from Plots import *

class IndividualPlots(Plot):


    @property
    def individuals(self):
        # individual telescopes
        return np.array([self.identifier(t) for t in self.tlcs])


    def identifier(self, tlc):
        return str(tlc)

    def setup(self, tlcs=None, stretch=2, epochs=[287,295,306,317], telescopes=['MEarth13','MEarth12','MEarth14', 'MEarth18', 'TRAPPIST', 'PISCOg', 'PISCOi'], **kwargs):
        # set up the axes for the plot
        self.tlcs = tlcs
        if telescopes is None:
            self.telescopes = np.unique([tlc.telescope for tlc in self.tlcs])
        else:
            self.telescopes= np.array(telescopes)
        if epochs is None:
            self.epochs = np.unique([tlc.epoch for tlc in self.tlcs])
        else:
            self.epochs = np.array(epochs)
        self.ncols = len(self.epochs)#np.ceil(np.sqrt(self.npanels)).astype(np.int)/2
        self.nrows = len(self.telescopes)#np.ceil(self.npanels/np.float(self.ncols)).astype(np.int)

        self.figure = plt.figure('light curves', figsize=(self.ncols*stretch*3, self.nrows*3), dpi=50)
        self.gs = plt.matplotlib.gridspec.GridSpec(self.nrows, self.ncols, hspace=0.1, wspace=0.1, bottom=0.1, left=0.1, right=0.9, top=0.9)

        # define all the axes
        self.axes = {}
        share = None
        for i in range(len(self.tlcs)):
            tlc = self.tlcs[i]
            ident = self.identifier(tlc)
            row = (self.telescopes == tlc.telescope).nonzero()[0]
            col = (self.epochs == tlc.epoch).nonzero()[0]

            if len(col) == 0:
                continue

            self.axes[ident] = plt.subplot(self.gs[row, col], sharex=share, sharey=share)
            share = self.axes[ident]

            if col != 0:
                plt.setp(self.axes[ident].get_yticklabels(), visible=False)
            else:
                self.axes[ident].set_ylabel('Relative Flux')
            if row != self.nrows-1:
                plt.setp(self.axes[ident].get_xticklabels(), visible=False)
            else:
                self.axes[ident].set_xlabel('Time from Mid-Transit (days)')




    def plot(self, tlcs, xlim=(-.075, 0.075), ylim=(0.99, 1.005), binsize=6.0/24.0/60.0, title=''):

        done = {}
        for tlc in self.tlcs:
            key = self.identifier(tlc)
            try:
                plt.sca(self.axes[key])
            except KeyError:
                continue
            ok = tlc.bad == False
            tm = tlc.TM
            tlc.plot(model=tm)
            try:
                done[key]
            except KeyError:
                tlc.setupSmooth()
                tm.plotPhased(linewidth=3,  alpha=0.5, color='black')
                done[key] = True

            plotbinned(tm.planet.timefrommidtransit(tlc.bjd)[ok], tlc.corrected()[ok], uncertainty=tlc.uncertainty[ok], bin=binsize,   alpha=0.5)

            if tlc.telescope == self.telescopes[0]:
                self.axes[key].xaxis.set_label_position("top")
                self.axes[key].set_xlabel('E={0}'.format(tlc.epoch))

            if tlc.epoch == self.epochs[-1]:
                self.axes[key].yaxis.set_label_position("right")
                self.axes[key].set_ylabel('{0}'.format(tlc.telescope))
        self.axes[key].set_xlim(*xlim)
        self.axes[key].set_ylim(*ylim)


        '''bjd.extend(tlc.bjd)
            flux.extend(tlc.corrected())
            uncertainty.extend(tlc.uncertainty*tlc.rescaling)
            bigok.extend(ok)
            tels = np.empty(tlc.n).astype(np.str)
            tels[:] = tlc.name.replace(' ','')
            whichtelescope.extend(tels)

            telescope = str(tlc)
            try:
                wholetel[telescope]
            except KeyError:
                wholetel[telescope] = dict(bjd=[],flux=[],uncertainty=[])


            wholetel[telescope]['bjd'].extend(tlc.bjd[ok])
            wholetel[telescope]['flux'].extend(tlc.corrected()[ok])
            wholetel[telescope]['uncertainty'].extend(tlc.uncertainty[ok]*tlc.rescaling)

            ax = self.axes[telescope]
            if len(self.getsownpanel) > 1:
                plt.sca(ax)
                tlc.plot(model=tm)
            plt.sca(self.axes['all'])
            tlc.plot(model=tm)'''

        '''if len(self.getsownpanel) > 1:
            for telescope in wholetel.keys():
                ax = self.axes[telescope]
                plt.sca(ax)

                for k in wholetel[telescope].keys():
                    wholetel[telescope][k] = np.array(wholetel[telescope][k])
                plotbinned(tm.planet.timefrommidtransit(wholetel[telescope]['bjd']), wholetel[telescope]['flux'], bin=binsize)
                plt.setp(ax.get_xticklabels(), visible=False)

            #tlc.plot(model=tm, alpha=0.1)



        bjd = np.array(bjd)
        flux = np.array(flux)
        uncertainty = np.array(uncertainty)
        bigok = np.array(bigok)
        telescope = np.array(telescope)
        plt.sca(self.axes['all'])
        plotbinned(tm.planet.timefrommidtransit(bjd[bigok]), flux[bigok], alpha=0.75, bin=binsize)
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
            tm.plotPhased(linewidth=3,  alpha=0.5, color='black')'''

        plt.savefig('individual_lightcurves.pdf')
