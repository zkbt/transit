from Plots import *

class PhasedPlots(Plot):
    def plot(self, tlcs, xlim=(-.1, 0.1), ylim=(0.99, 1.005), binsize=6.0/24.0/60.0, title=''):



        # figure out height ratios
        telescopes = np.unique([tlc.telescope for tlc in tlcs])
        if len(telescopes) == 1:
            h = 1
            n = 1
        else:
            h = np.ones(len(telescopes)+1)
            h[-1] += 0#2
            n = len(telescopes) + 1

        # create a figure to hold all the phased light curves
        plt.figure('all light curves', figsize=(10,n*4), dpi=50)
        gs = plt.matplotlib.gridspec.GridSpec(n, 1, hspace=0.1, wspace=0, height_ratios=h, bottom=0.2)


        bjd, flux, uncertainty, telescope, bigok = [],[],[],[], []
        if len(telescopes) == 1:
            ax_all = plt.subplot(111)
        else:
            ax_all = plt.subplot(gs[-1])
        axes = {}
        count =0
        wholetel = {}
        for i in range(len(tlcs)):

            # pull out this light curve
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
            telescope.extend(tels)

            this = (telescopes == tlc.telescope).nonzero()[0][0]
            try:
                axes[this]
            except:
                if len(telescopes) == 1:
                    axes[this] = ax_all
                else:
                    axes[this]= plt.subplot(gs[count], sharex=ax_all, sharey=ax_all)
                    count += 1
                    axes[this].set_ylabel(telescopes[this])
                wholetel[this] = dict(bjd=[],flux=[],uncertainty=[])

            wholetel[this]['bjd'].extend(tlc.bjd[ok])
            wholetel[this]['flux'].extend(tlc.corrected()[ok])
            wholetel[this]['uncertainty'].extend(tlc.uncertainty[ok]*tlc.rescaling)

            ax = axes[this]
            if len(telescopes) > 1:
                plt.sca(ax)
                tlc.plot(model=tm)
            plt.sca(ax_all)
            tlc.plot(model=tm)

        if len(telescopes) > 1:
            for this in wholetel.keys():
                ax = axes[this]
                plt.sca(ax)

                for k in wholetel[this].keys():
                    wholetel[this][k] = np.array(wholetel[this][k])
                plotbinned(tm.planet.timefrommidtransit(wholetel[this]['bjd']), wholetel[this]['flux'], bin=binsize)
                plt.setp(ax.get_xticklabels(), visible=False)

            #tlc.plot(model=tm, alpha=0.1)



        bjd = np.array(bjd)
        flux = np.array(flux)
        uncertainty = np.array(uncertainty)
        bigok = np.array(bigok)
        telescope = np.array(telescope)
        plt.sca(ax_all)
        plotbinned(tm.planet.timefrommidtransit(bjd[bigok]), flux[bigok], alpha=0.75, bin=binsize)
        ax_all.set_xlim(*xlim)
        ax_all.set_ylim(*ylim)
        ax_all.set_xlabel('Time from Mid-Transit (days)')
        ax_all.set_ylabel('Relative Flux')
        ax_all.set_title(title)
        table = astropy.table.Table(dict(bjd=bjd,
                                    flux=flux,
                                    uncertainty=uncertainty,
                                    telescope=telescope))

        table.write('merged_lc.txt',
                    format='ascii.fixed_width',
                    bookend=False)
        tlc.setupSmooth()
        tm.plotPhased(linewidth=3,  alpha=0.5, color='black')
        for k in axes.keys():
            a = axes[k]
            plt.sca(a)
            tm.plotPhased(linewidth=3,  alpha=0.5, color='black')

        plt.savefig('combined_lightcurves.pdf')
