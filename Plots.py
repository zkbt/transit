from imports import *

def plotbinned(x, flux, uncertainty=None, alpha=0.5, bin=6.0/24.0/60.0, **kwargs):
    bx, by, be = zachopy.oned.binto(x, flux,  bin, yuncertainty=uncertainty, robust=False, sem=True)
    plt.errorbar(bx, by, be, color='black', alpha=alpha, elinewidth=3, marker='o', linewidth=0, capthick=3, **kwargs)

def ink_errorbar(x, y, yerr, colors, grayscale=False, alpha=1.0, **kwargs):
    for i in range(len(x)):
        color = colors[i] + 0.0
        assert(len(color) == 4)
        color[-1] *= alpha
        if grayscale:
            color[0:2] = np.mean(color[0:2])
        plt.errorbar(x[i], y[i], yerr[i], color=color, **kwargs)

class Plot(Talker):
    def __init__(self, **kwargs):
        Talker.__init__(self)
        self.setup(**kwargs)
        self.plot(**kwargs)
