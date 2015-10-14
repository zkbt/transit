from imports import *

def plotbinned(x, flux, uncertainty=None, alpha=0.5, bin=6.0/24.0/60.0, **kwargs):
    bx, by, be = zachopy.oned.binto(x, flux,  bin, yuncertainty=uncertainty, robust=False, sem=True)
    plt.errorbar(bx, by, be, color='black', alpha=alpha,  **kwargs)

def ink_errorbar(x, y, yerr=None, xerr=None, colors=None, grayscale=False, alpha=1.0, **kwargs):

    kw = {}
    yerrtoplot = yerr
    xerrtoplot = xerr
    for k in kwargs.keys():
        kw[k] = kwargs[k]
    for i in range(len(x)):
        color = colors[i] + 0.0
        assert(len(color) == 4)
        color[-1] *= alpha
        if grayscale:
            color[0:3] = np.mean(color[0:3])
        kw['color'] = color
        kw['ecolor'] = color
        kw['markeredgecolor'] = color
        kw['markerfacecolor'] = color
        #print color
        if yerr is not None:
            if len(yerr.shape) == 1:
                yerrtoplot = yerr[i]
            if len(yerr.shape) == 2:
                yerrtoplot = [yerr[:,i]]

        if xerr is not None:
            if len(xerr.shape) == 1:
                xerrtoplot = xerr[i]
            if len(xerr.shape) == 2:
                xerrtoplot = [xerr[:,i]]

        plt.errorbar(x[i], y[i], xerr=xerrtoplot, yerr=yerrtoplot, **kw)

class Plot(Talker):
    def __init__(self, **kwargs):
        Talker.__init__(self)
        self.setup(**kwargs)
        self.plot(**kwargs)
