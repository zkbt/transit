from Plots import *

def modify(s):
    if 'MEarth' in s:
        return s.replace('1', ' #')
    else:
        return s

class IndividualPlots(Plot):

    @property
    def individuals(self):
        # individual telescopes
        return np.array([self.identifier(t) for t in self.tlcs])

    def identifier(self, tlc):
        return str(tlc)

    def setup(self, tlcs=None, stretch=2.0, epochs=[-19, -11,  0, 11], telescopes=['MEarth13','MEarth12','MEarth14', 'MEarth18', 'TRAPPIST', 'PISCOg', 'PISCOi'], gskw=dict(bottom=0.05, left=0.1, right=0.95, top=0.95), dpi=30, **kwargs):
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

        inchinmm = 25.4/2.0
        scale = 183/inchinmm
        self.figure = plt.figure('light curves', figsize=(scale, scale*self.nrows/(self.ncols*stretch)), dpi=dpi)
        self.gs = plt.matplotlib.gridspec.GridSpec(self.nrows, self.ncols, hspace=0.1, wspace=0.05, **gskw)

        # define all the axes
        self.axes = {}
        self.maxrow, self.maxcol, self.mincol = {}, {}, {}
        share = None
        for i in range(len(self.tlcs)):
            tlc = self.tlcs[i]
            ident = self.identifier(tlc)
            row = (self.telescopes == tlc.telescope).nonzero()[0]
            col = (self.epochs == tlc.epoch).nonzero()[0]

            if (len(col) == 0) or (len(row) == 0):
                continue
            else:
                row = row[0]
                col = col[0]


            self.axes[ident] = plt.subplot(self.gs[row, col], sharex=share, sharey=share)
            share = self.axes[ident]

            try:
                self.maxrow[tlc.epoch] = np.maximum(row, self.maxrow[tlc.epoch])
            except KeyError:
                self.maxrow[tlc.epoch] = row

            try:
                self.mincol[tlc.telescope] = np.minimum(col, self.mincol[tlc.telescope])
            except KeyError:
                self.mincol[tlc.telescope] = col


            try:
                self.maxcol[tlc.telescope] = np.maximum(col, self.maxcol[tlc.telescope])
            except KeyError:
                self.maxcol[tlc.telescope] = col


    def plot(self, tlcs, synthesizer=None, xlim=(-1.5/24, 1.5/24), ylim=(0.985, 1.015), binsize=5.0/24.0/60.0, title='',  epochs=None, telescopes=None, gskw=None, **kwargs):
        self.synthesizer=synthesizer
        done = {}
        self.synthesizer.fromPDF(option='best')
        bestdeterministic = {}
        # first plot the central value of the fit and the GP
        for tlc in self.tlcs:
            # select this particular light curve
            key = self.identifier(tlc)
            ident = key
            row = (self.telescopes == tlc.telescope).nonzero()[0]
            col = (self.epochs == tlc.epoch).nonzero()[0]
            if (len(col) == 0) or (len(row) == 0):
                continue
            else:
                row = row[0]
                col = col[0]



            # if a panel hasn't be set; skip it!
            try:
                plt.sca(self.axes[key])
            except KeyError:
                continue

            if self.synthesizer.gp:
                points = tlc.gp_points()
            else:
                points = tlc.points()

            if 'PISCOg' in tlc.telescope:
                nudge = -0.005
            else:
                nudge = 0.0

            plt.plot(points['t'], points['raw']+nudge, color='gray', marker='o', markersize=3, markeredgecolor='none', linewidth=0, alpha=0.5)

            bx, by, be = zachopy.oned.binto(points['t'], points['raw'], binsize,
                yuncertainty=tlc.effective_uncertainty[tlc.ok], robust=False, sem=True)
            plt.errorbar(bx, by+nudge, be, color='black', alpha=0.75, elinewidth=3, markersize=6, linewidth=0, capthick=0, zorder=40, **kwargs)


            #plotbinned(points['t'], points['raw'], color=tlc.color, marker='o', markersize=10, alpha=1)


            if self.synthesizer.gp:
                lines = tlc.gp_lines(mean=True)
            else:
                lines = tlc.lines()
            bestdeterministic[key] = lines['raw'] - lines['raw']
            plt.plot(lines['t'], lines['raw']+nudge, linewidth=1, alpha=1, color='#1f78b4', zorder=20)

            #print 'best!'
            hyperparameters = tlc.TM.instrument.gplna.value, tlc.TM.instrument.gplntau.value
            #print j, i, which, walker
            #print hyperparameters

            #print row, col, self.maxrow[tlc.epoch], self.mincol[tlc.telescope], self.maxcol[tlc.telescope]
            if tlc.telescope == self.telescopes[0]:
                a = self.axes[key].twiny()
                a.plot(None)
                a.xaxis.set_label_position("top")
                a.set_xlabel('E={0}'.format(tlc.epoch))
                plt.setp(a.get_xticklabels(), visible=False)
                #print 'adding title'


            if row == self.maxrow[tlc.epoch]:
                self.axes[key].xaxis.set_label_position("bottom")
                self.axes[ident].set_xlabel('Time from Mid-Transit (days)')
                plt.xticks([-0.05, 0, 0.05])

                #print 'adding xlabel'
            else:
                plt.setp(self.axes[ident].get_xticklabels(), visible=False)
                #print 'hiding xlabel'

            if col == self.mincol[tlc.telescope]:
                self.axes[key].yaxis.set_label_position("left")
                self.axes[ident].set_ylabel('Relative\nFlux')
                plt.yticks([0.995, 1.0, 1.005])

                #print 'adding ylabel'
            else:
                plt.setp(self.axes[ident].get_yticklabels(), visible=False)
                #print 'hiding ylabel'

            if col == self.maxcol[tlc.telescope]:
                a = self.axes[key].twinx()
                a.yaxis.set_label_position("right")
                a.plot(None)
                plt.setp(a.get_yticklabels(), visible=False)

                a.set_ylabel(modify(tlc.telescope), rotation=270, labelpad=15)
                #print 'adding ytitle'

        for i in range(5):
            # point at a random sample in the PDF
            self.synthesizer.fromPDF(option='random')#, verbose=True)

            # loop over the light curves and plot them
            for tlc in self.synthesizer.tlcs:
                #tlc.pithy=False
                # select this particular light curve
                key = self.identifier(tlc)

                # if a panel hasn't be set; skip it!
                try:
                    plt.sca(self.axes[key])
                except KeyError:
                    continue

                if self.synthesizer.gp:
                    extralines = tlc.gp_lines(mean=False)
                    #extralines = tlc.gp_lines(mean=True)
                else:
                    extralines = tlc.lines()
                    #assert(False)
                #self.speak('{0}'.format(tlc))

                if 'PISCOg' in tlc.telescope:
                    nudge = -0.005
                else:
                    nudge = 0.0

                self.axes[key].plot(extralines['t'], extralines['raw'] - bestdeterministic[key]+nudge, linewidth=1, alpha=0.75, color='#b2df8a', zorder=10)
                hyperparameters = tlc.TM.instrument.gplna.value, tlc.TM.instrument.gplntau.value
                #print j, i, which, walker
                #print hyperparameters
                #tlc.pithy=True
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

        '''if pdf:
            plt.savefig('individual_lightcurves.pdf')
        else:
            plt.draw()'''
