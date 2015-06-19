from imports import *

def plotbinned(x, flux, alpha=0.5, bin=6.0/24.0/60.0):
    bx, by, be = zachopy.oned.binto(x, flux, bin, robust=False, sem=True)
    plt.errorbar(bx, by, be, color='black', alpha=alpha, elinewidth=3, marker='o', linewidth=0, capthick=3)

class Plot(Talker):
    def __init__(self):
        Talker.__init__(self)


class Phased(Plot):
    def plot(self, tlcs):

        # create a figure to hold all the phased light curves
        plt.figure('all light curves', figsize=(10,15), dpi=30)

        # figure out height ratios
        telescopes = np.unique([tlc.telescope for tlc in tlcs])
        h = np.ones(len(telescopes)+1)
        h[-1] += 0#2

        gs = plt.matplotlib.gridspec.GridSpec(len(telescopes)+1, 1, hspace=0.1, wspace=0, height_ratios=h)


        bjd, flux, uncertainty, telescope, bigok = [],[],[],[], []
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
                axes[this]= plt.subplot(gs[count], sharex=ax_all, sharey=ax_all)
                count += 1
                axes[this].set_ylabel(telescopes[this])
                wholetel[this] = dict(bjd=[],flux=[],uncertainty=[])

            wholetel[this]['bjd'].extend(tlc.bjd[ok])
            wholetel[this]['flux'].extend(tlc.corrected()[ok])
            wholetel[this]['uncertainty'].extend(tlc.uncertainty[ok]*tlc.rescaling)

            ax = axes[this]
            plt.sca(ax)
            tlc.plot(model=tm)
            plt.sca(ax_all)
            tlc.plot(model=tm)

        for this in wholetel.keys():
            ax = axes[this]
            plt.sca(ax)

            for k in wholetel[this].keys():
                wholetel[this][k] = np.array(wholetel[this][k])
            plotbinned(tm.planet.timefrommidtransit(wholetel[this]['bjd']), wholetel[this]['flux'])
            plt.setp(ax.get_xticklabels(), visible=False)

            #tlc.plot(model=tm, alpha=0.1)



        bjd = np.array(bjd)
        flux = np.array(flux)
        uncertainty = np.array(uncertainty)
        bigok = np.array(bigok)
        telescope = np.array(telescope)
        plt.sca(ax_all)
        plotbinned(tm.planet.timefrommidtransit(bjd[bigok]), flux[bigok], alpha=0.75, bin=6.0/24.0/60.0)
        ax_all.set_xlim(-0.1, 0.1)
        ax_all.set_ylim(0.99, 1.006)
        ax_all.set_xlabel('Time from Mid-Transit (days)')
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
