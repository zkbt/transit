from Plots import *

class IndividualPlots(Plot):

    @property
    def individuals(self):
        # individual telescopes
        return np.array([self.identifier(t) for t in self.tlcs])

    def identifier(self, tlc):
        return str(tlc)

    def setup(self, tlcs=None, stretch=2, epochs=[-19, -11,  0, 11], telescopes=['MEarth13','MEarth12','MEarth14', 'MEarth18', 'TRAPPIST', 'PISCOg', 'PISCOi'], dpi=30, **kwargs):
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

        self.figure = plt.figure('light curves', figsize=(self.ncols*stretch*3, self.nrows*3), dpi=dpi)
        self.gs = plt.matplotlib.gridspec.GridSpec(self.nrows, self.ncols, hspace=0.1, wspace=0.1, bottom=0.1, left=0.1, right=0.9, top=0.9)

        # define all the axes
        self.axes = {}

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

            if col != 0:
                plt.setp(self.axes[ident].get_yticklabels(), visible=False)
            else:
                self.axes[ident].set_ylabel('Relative Flux')
            if row != self.nrows-1:
                plt.setp(self.axes[ident].get_xticklabels(), visible=False)
            else:
                self.axes[ident].set_xlabel('Time from Mid-Transit (days)')


    def plot(self, tlcs, synthesizer=None, xlim=(-.075, 0.075), ylim=(0.99, 1.005), binsize=6.0/24.0/60.0, title='',  **kwargs):
        self.synthesizer=synthesizer
        done = {}
        self.chain = synthesizer.sampler.flatchain
        self.best = self.chain[np.argmax(synthesizer.sampler.flatlnprobability)]
        self.synthesizer.fromArray(self.best)
        bestdeterministic = {}
        # first plot the central value of the fit and the GP
        for tlc in self.tlcs:
            # select this particular light curve
            key = self.identifier(tlc)

            # if a panel hasn't be set; skip it!
            try:
                plt.sca(self.axes[key])
            except KeyError:
                continue

            if self.synthesizer.gp:
                points = tlc.gp_points()
            else:
                points = tlc.points()
            plt.plot(points['t'], points['cleaned'], color=tlc.color, marker='o', linewidth=0, alpha=0.25)

            if self.synthesizer.gp:
                lines = tlc.gp_lines(mean=True)
            else:
                lines = tlc.lines()
            bestdeterministic[key] = lines['raw'] - lines['cleaned']
            plt.plot(lines['t'], lines['cleaned'], linewidth=2, alpha=1, color='black')

            #print 'best!'
            hyperparameters = tlc.TM.instrument.gplna.value, tlc.TM.instrument.gplntau.value
            #print j, i, which, walker
            #print hyperparameters


            if tlc.telescope == self.telescopes[0]:
                self.axes[key].xaxis.set_label_position("top")
                self.axes[key].set_xlabel('E={0}'.format(tlc.epoch))

            if tlc.epoch == self.epochs[-1]:
                self.axes[key].yaxis.set_label_position("right")
                self.axes[key].set_ylabel('{0}'.format(tlc.telescope))

        # then, plot some samples
        colors = {2:'blue', 3:'gray'}
        chain = self.synthesizer.sampler.chain
        nwalkers, nsteps, ndim = chain.shape
        #self.speak('shape is {0}'.format(chain.shape))
        #for j in [2,3]:
        for j in [3]:
            quartile = nsteps*j/4
            for i in range(5):
                which = int(np.random.uniform(quartile,quartile+nsteps/4))
                walker = int(np.random.uniform(0, nwalkers))
                #new = np.array(self.best) + 0.0#
                #new[-3] = new[-3] + i*0.01
                new = chain[walker, which,:]

                self.synthesizer.fromArray(new)#, verbose=True)
                for tlc in self.synthesizer.tlcs:
                    tlc.pithy=False
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
                    self.axes[key].plot(extralines['t'], extralines['raw'] - bestdeterministic[key], linewidth=1, alpha=0.5, color=colors[j])
                    hyperparameters = tlc.TM.instrument.gplna.value, tlc.TM.instrument.gplntau.value
                    #print j, i, which, walker
                    #print hyperparameters
                    tlc.pithy=True
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
