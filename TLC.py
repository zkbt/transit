from imports import *
from Planet import Planet
from Star import Star
from TM import TM
from PDF import PDF

ppm = 1e6

class TLC(Talker):
    '''Transit Light Curve class,
        to store both light curve and auxiliary variables.'''

    def fromFile(self, filename):
        '''specific light curve formats inheriting from TLC
            might need to replace this default reader'''

        assert(filename is not None)
        # read in the ascii file
        data = astropy.io.ascii.read(filename)

        # create bjd, flux, and uncertainty arrays
        #  (converting magnitudes to relative flux)
        try:
            bjd = data['bjd']
            flux = data['flux']
            uncertainty = data['uncertainty']
        except KeyError:
            raise KeyError("your TLC file {0} doesn't have 'bjd', 'flux', 'uncertainty'".format(filename))
        self.speak(
            'added {n} times, fluxes, and uncertainties'.format(n=len(bjd)))
        kwargs = {}
        for key in data.colnames:
            kwargs[key] = data[key]
            self.speak('added {0} elements of {1}'.format(len(kwargs[key]), key))


        medflux = np.median(data['flux'])
        kwargs['uncertainty'] /= medflux
        kwargs['flux'] /= medflux

        # populate the light curve from these arrays
        self.fromArrays(**kwargs)

    def __init__(self, bjd=None, flux=None, uncertainty=None,
                        inputfilename=None, directory=None,
                        path=None,
                        left=None, right=None,
                        name=None,
                        remake=False,
                        color='slategray',
                        isfake=False,
                        telescope=None, epoch=None,
                        **kwargs):

        # initialize the Talker object
        Talker.__init__(self)
        self.pithy=True
        self.color = color
        self.speak('creating an empty TLC')

        self.rescaling = 1.0
        self.beta = 1.0

        # define a dictionary of flags that can be used for bad data
        self.flags = dict(outlier=1, saturation=2, custom=4)

        # specify the left and right wavelengths of this bandpass
        self.left = left
        self.right = right
        #print '!!!!!!!!', self.left, self.right
        # keep track of the telescope (and epoch)

        self.telescope=telescope
        self.epoch=epoch

        # is this a real light curve, or a fake one?
        #  (e.g., one at high resolution for plotting)
        self.isfake = isfake

        # specify the original filename
        self.inputfilename = inputfilename

        # specify a directory in which to store saved versions of this TLC,
        #     as well as other outputs (plots, etc...)
        if path is not None:
            zachopy.utils.mkdir(path)
            directory = path + 'T={telescope}_E={epoch}/'.format(**self.__dict__)

        self.directory = directory
        if self.directory is None:
            try:
                self.directory = '/'.join(self.inputfilename.replace('data/', 'fits/').split('/')[:-1]) + '/'
                zachopy.utils.mkdir(self.directory)
            except AttributeError:
                pass
        self.speak('assigning it the directory {0}'.format(self.directory))

        # initialize the TLC by filling it with data
        self.initialize(bjd=bjd, flux=flux, uncertainty=uncertainty,
                        remake=remake,
                        telescope=telescope, epoch=epoch, **kwargs)



        # pick the central wavelength of this light curve's bandpass
        if self.left is None or self.right is None:
            self.wavelength = None
        else:
            self.wavelength = (self.left + self.right)/2.0


        # assign a name to this lightcurve
        if name is None:
            if self.telescope is None:
                name = '???'
            else:
                name = self.telescope
            if self.epoch is None:
                pass
            else:
                name += ',E={0}'.format(self.epoch)
        self.name = name

        # assign the colors for this light curve
        self.setupColors(color=color)


    def initialize(self, bjd=None, flux=None, uncertainty=None,
                        remake=False, telescope=None, epoch=None,
                        **kwargs):
        '''If possible, try to load the light curve from its directory
            otherwise, create it from raw input file.'''

        try:
            assert(remake == False)
            self.load(self.directory)
            assert(self.bad.shape == self.flux.shape)
            self.speak('initialized TLC from pre-saved file in {0}'.format(self.directory))
            haschanged = False
        except (IOError,AssertionError):
            self.speak("failed to load!")
            # if possible, initialize from arrays; if not, load from scratch
            if bjd is not None and flux is not None:
                self.fromArrays(bjd, flux, uncertainty, **kwargs)
                self.speak('initialized TLC from {0}-element arrays'.format(self.n))
            else:
                self.fromFile(self.inputfilename, **kwargs)
                self.speak('initialized TLC from {0}'.format(self.inputfilename))
            haschanged = True

        self.telescope=telescope
        self.epoch = epoch
        try:
            self.flux
        except AttributeError:
            self.bad = np.array([])
            return None

        # make sure an array of "bad" values is defined
        try:
            self.bad
        except AttributeError:
            haschanged = True
            self.speak("$$$$$$$$$ bad wasn't defined!")
            try:
                self.bad = kwargs['bad']
            except KeyError:
                self.bad = np.isfinite(self.flux) == False
        assert(self.bad.shape == self.flux.shape)

        if haschanged & (self.isfake == False):
            zachopy.utils.mkdir(self.directory)
            self.save(self.directory)
    @property
    def effective_uncertainty(self):
        return self.uncertainty*self.rescaling*self.beta

    def setupColors(self, color='eye', minimumuncertainty=0.001):
        '''Method to set the line and point colots for this light curve.
            color = 'eye': make colors as they would appear to the human eye
            [add other options (e.g. specify a color specific color)]'''

        # set up the appropriate colors to use
        try:
            self.colors
        except AttributeError:
            self.colors = {}
        self.color = color

        if type(self.color) == str:
            self.colors = {}
            if color=='eye':
                self.colors['points'] = zachopy.color.nm2rgb([self.left/10, self.right/10], intensity=1.0)
                self.colors['lines'] = zachopy.color.nm2rgb([self.left/10, self.right/10], intensity=2.0)
            else:
                self.colors['lines'] = color
                r, g, b = zachopy.color.name2color(color.lower())
                rgba = np.zeros((self.n, 4))
                rgba[:,0] = r
                rgba[:,1] = g
                rgba[:,2] = b
                eu = self.effective_uncertainty
                eu[np.isfinite(eu) == False] = np.inf
                weights = np.minimum((minimumuncertainty/eu)**2, 1)
                over = weights > 1
                weights[over] = 1
                rgba[:,3] = 1.0*weights
                self.colors['points'] = rgba
        else:
            self.colors['points'] = self.color.color(self.wavelength/10)[0:3]
            self.colors['lines'] = 'gray'
    @property
    def ok(self):
        return self.bad == False

    def plot(self, model=None, alpha=1):
        ok = self.bad == False
        if model == None:
            x = self.bjd
        else:
            x = model.planet.timefrommidtransit(self.bjd)
        colors = self.colors['points']
        colors[:,3]*=alpha
        plt.scatter(x[ok], self.corrected()[ok], color=colors[ok], edgecolor='none', s=50, linewidth=0)
        colors[self.bad,3]*0.2
        plt.scatter(x[self.bad], self.corrected()[self.bad], color=colors[self.bad], marker='x', s=50)



    def restrictToNight(self, night=None):
        '''Trim to data from only a particular night.'''

        # if a night has been selected, restrict to it
        if night is not None:
            ok = np.abs(bjd - night) < 0.5
        else:
            ok = bjd > 0
        self.trimTo(ok)


    def trimTo(self, ok):
        '''Trim the TLC to data points where ok == True.'''

        self.speak('trimming from ')
        self.bjd


    def chisq(self):
        ok = self.bad == False
        return np.sum((self.residuals()/self.uncertainty)[ok]**2)

    def gp_compute(self, hyperparameters):
        a, tau = np.exp(hyperparameters[:2])

        #self.gp = george.GP(a*george.kernels.ExpSquaredKernel(tau))#, solver=george.HODLRSolver)
        ok = self.bad == False
        t = self.bjd[ok]
        if len(t) > 5000:
            solver = george.HODLRSolver
        else:
            solver = george.BasicSolver

        self.gp = george.GP(a*george.kernels.Matern32Kernel(tau), solver=solver)

        yerr = self.uncertainty[ok]*self.rescaling
        before = time.clock()
        self.gp.compute(t, yerr)
        after = time.clock()

        self.speak('used {3} to compute kernel {0} for {1} data points in {2} microseconds'.format(self.gp, len(t), 1e6*(after-before), solver))

    def gp_lnlikelihood(self):
        '''Return the GP calculated likelihood of *this* light curve, assuming the (not hyper-)parameters have been set elsewhere.'''

        #WHAT ER THE HYPERPAREMETES?
        hyperparameters = self.TM.instrument.gplna.value, self.TM.instrument.gplntau.value
        self.gp_compute(hyperparameters)
        ok = self.bad == False

        before = time.clock()
        lnp = self.gp.lnlikelihood(self.residuals()[ok], quiet=True)
        after = time.clock()

        self.speak('computed likelihood for {0} data points in {1} microseconds'.format(np.sum(ok), 1e6*(after-before)))
        return lnp

    def bestBeta(self, timescale=15.0/60.0/24.0):
        ok = self.ok
        if np.sum(ok) > 0:
            x = self.bjd[ok]
            y = self.residuals()[ok]
            err = self.uncertainty[ok]*self.rescaling

            bx, by, be = zachopy.oned.binto(x=x, y=y, yuncertainty=err, binwidth=timescale, sem=True, robust=False)
            fine = np.isfinite(by)
            if np.sum(fine) > 1:
                inflation = np.sqrt(np.mean(((by - 0)**2/be**2)[fine]))
                inflation = np.maximum(inflation, 1)
            else:
                inflation = 1.0
            assert(np.isfinite(inflation))
        else:
            inflation = 1.0
        self.speak('the best beta for {0} is {1} atop the rescaling of {2}'.format(self.name, inflation, self.rescaling))
        return inflation

    def fromArrays(self, bjd, flux, uncertainty=None, **kwargs):
        '''Populate a TLC from input arrays (used by transmission.py)'''

        # how many data points are in light curve?
        self.n = len(bjd)

        # define the times
        self.bjd = np.array(bjd)

        # define the flux array, and normalize it to its median
        self.flux = np.array(flux)

        # make sure the uncertainty is defined
        if uncertainty is None:
            uncertainty = np.ones_like(self.flux)
        else:
            self.uncertainty = np.array(uncertainty)

        self.uncertainty /= np.median(flux)
        self.flux /= np.median(flux)


        # populate the external variables, both as a list and as individual entries
        self.externalvariables = {}
        for key, value in kwargs.iteritems():
            if len(value) == len(self.bjd):
                if key != 'bjd' and key != 'flux' and key !='uncertainty' and key !='left' and key !='right' and key !='wavelength':#and key!= 'bad':
                    self.externalvariables[key] = value
                else:
                    self.speak( "   " +  key+  " was skipped")
            else:
                self.speak(key + ' has length '+  str(len(value)))

    def save(self, directory, verbose=True):
        '''Tool to save a light curve to a directory.'''
        self.directory = directory
        tosave = {}
        for k in self.__dict__.keys():
            if k != 'TM' and k!= 'TLC' and 'ax_' not in k and 'points_' not in k and 'line_' not in k and 'figure_' not in k:
                tosave[k] = self.__dict__[k]
                self.speak( " saved " + k)
                #if k == 'externalvariables':
                    #self.speak( "      " + k + ', including:')
                    #for evkey in tosave[k].keys():
                        #self.speak( "          " + evkey)

        zachopy.utils.mkdir(directory)
        filename = directory + str(self.__class__).split('.')[-1].split("'")[0] + '.npy'

        np.save(filename, tosave)
        self.speak("saving light curve to {directory}, including all its external variables".format(directory=directory))

    def load(self, directory):
        self.directory = directory
        filename = directory + str(self.__class__).split('.')[-1].split("'")[0] + '.npy'
        self.speak('trying to load TLC from {0}'.format(filename))
        loaded = np.load(filename)[()]
        for key in loaded.keys():
            self.speak( ' loaded ' + key)
            self.__dict__[key] = loaded[key]



    def setupDiagnostics(self, everything=True):

        '''Setup the axes needed for plotting a light curve.'''

        # create a figure to populate
        if everything:
            label = 'diagnostics'
        else:
            label = 'summary'
        self.figure_diagnostics = plt.figure('light curve {0}'.format(label), figsize=(20,12), dpi=50)
        try:
            # if the plot window is already set up, don't do anything!
            self.ax_raw
        except:



            # set up the light curve plots
            if everything:
                # create two columns, to populate with lightcurves on left and diagnostics on right
                gs_overarching = plt.matplotlib.gridspec.GridSpec(1, 3, width_ratios=[1, 0.3, 1.0], wspace=0.25, left=0.09, right=0.98)

                # create plots for light curves
                gs_lightcurve = plt.matplotlib.gridspec.GridSpecFromSubplotSpec(5, 2, hspace=0.05, wspace=0, width_ratios = [1,.2], height_ratios=[1,0.5,0.1,1,0.5], subplot_spec = gs_overarching[0])

                # noise properties
                gs_diagnostics = plt.matplotlib.gridspec.GridSpecFromSubplotSpec(3,1, height_ratios = [3,1,1], hspace=0.3, subplot_spec = gs_overarching[1])

                # potentially correlated variables
                self.to_correlate = []
                #for k in self.TM.floating:
                #    if k != 'C':
                #        if k in self.TM.instrument.__dict__.keys() and k != 'rescaling':
                #            self.to_correlate.append(k.split('_tothe')[0])
                for k in self.externalvariables.keys():
                    if np.std(self.externalvariables[k]) > 0:
                        self.to_correlate.append(k.split('_tothe')[0])
                ncor = np.int(np.ceil(np.sqrt(len(self.to_correlate))))
                gs_external = plt.matplotlib.gridspec.GridSpecFromSubplotSpec(len(self.to_correlate), 2, subplot_spec = gs_overarching[2], width_ratios=[1,5], hspace=0.1, wspace=0)
                assert(False)
            else:
                gs_lightcurve = plt.matplotlib.gridspec.GridSpec(5, 1, hspace=0.05, wspace=0, height_ratios=[1,0.5,0.1,1,0.5])

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
                #plt.setp(a.get_xticklines(), visible=False)
                #a.set_yticks([])
                #a.set_xticks([])
                plt.setp(a.get_xticklabels(), visible=False)
                plt.setp(a.get_yticklabels(), visible=False)

            # set up the labels for the light curve panels
            try:
                self.ax_raw.set_title('{name} | {left:.0f}-{right:.0f} angstroms'.format(name=self.name, left=self.left, right=self.right))
            except:
                pass
            self.ax_raw.set_ylabel('Basic Photometry')
            self.ax_instrument.set_ylabel('transit\nresiduals\n(ppm)')
            self.ax_corrected.set_ylabel('Corrected Photometry')
            self.ax_residuals.set_ylabel('final\nresiduals\n(ppm)')
            self.ax_residuals.set_xlabel('Time from Mid-Transit (days)')

            # set up the y limits (is this necessary?)
            self.ax_raw.set_ylim(np.min(self.flux), np.max(self.flux))
            self.ax_corrected.set_ylim(np.min(self.flux), np.max(self.flux))


            if everything:

                # create a plot to store an autocorrelation function
                #self.ax_text = plt.subplot(gs_diagnostics[0])

                try:
                    todisplay = self.TM.lastfit.pdf.listParameters()

                    todisplay += '\n'
                    todisplay += 'chisq = {chisq:.1f}/{dof:.0f}\nrescaling = {rescaling:.2f}'.format(**self.TM.lastfit.notes)

                    self.ax_text.text(-0.1, 1.1, todisplay, va='top', fontsize=11)
                    self.ax_text.axis('off')
                except AttributeError:
                    self.speak('no fit was found to print')
                # create a plot to store an autocorrelation function
                self.ax_acf = plt.subplot(gs_diagnostics[1])
                self.ax_acf.set_xlabel('Lag (in datapoints)')
                self.ax_acf.set_ylabel('ACF')

                # create a plot to show the RMS as a function of binning
                self.ax_binnedrms = plt.subplot(gs_diagnostics[2])
                self.ax_binnedrms.set_xlabel('# of datapoints in a bin')
                self.ax_binnedrms.set_ylabel('Binned RMS (ppm)')
                #self.ax_parameters = plt.subplot(gs_diagnostics[0], frameon=False)
                self.ax_correlations, self.ax_timeseries = {}, {}

                for i in range(len(self.to_correlate)):
                    k = self.to_correlate[i]
                    self.ax_correlations[k] = plt.subplot(gs_external[i,0])
                    a = self.ax_correlations[k]
                    plt.setp(a.get_xticklabels(), visible=False)
                    plt.setp(a.get_yticklabels(), visible=False)
                    a.set_ylabel(k, rotation=45, ha='right', size=9)

                    self.ax_timeseries[k] = plt.subplot(gs_external[i,1], sharex=self.ax_residuals, sharey=self.ax_correlations[k])
                    a = self.ax_timeseries[k]
                    plt.setp(a.get_xticklabels(), visible=False)
                    plt.setp(a.get_yticklabels(), visible=False)



    def setupLightcurvePlots(self, everything=True):

        '''Setup the axes needed for plotting a light curve, with both unphased and phased.'''
        try:
            # if the plot window is already set up, don't do anything!
            self.ax_unphased

        except:

            # set up the grid for plotting
            gs_lightcurve = plt.matplotlib.gridspec.GridSpec(2, 2, width_ratios=[2,1], wspace=0.0, hspace=0.25)

            # set up the light curve (and residual) panels (leaving a space between the uncorrected and the corrected
            self.ax_phased = plt.subplot(gs_lightcurve[0,0])
            self.ax_unphased = plt.subplot(gs_lightcurve[1,0], sharey=self.ax_phased)
            self.ax_phased_zoom = plt.subplot(gs_lightcurve[0,1], sharey=self.ax_phased)
            self.ax_unphased_zoom = plt.subplot(gs_lightcurve[1,1], sharey=self.ax_phased, sharex=self.ax_phased_zoom)

            # hide the tick labels on most of the light curve panels
            for a in [self.ax_phased_zoom]:
                plt.setp(a.get_xticklabels(), visible=False)

            for a in [self.ax_phased_zoom, self.ax_unphased_zoom]:
                plt.setp(a.get_yticklabels(), visible=False)

            # set up the labels for the light curve panels
            self.ax_unphased.set_xlabel('Time since {0:.3f}'.format(self.TM.planet.t0.value))
            self.ax_phased.set_xlabel('Phased Time from Mid-transit (days)')
            self.ax_unphased_zoom.set_xlabel('Time from Mid-transit (days)')
            self.ax_unphased.set_ylabel('Relative Flux')
            self.ax_phased.set_ylabel('Relative Flux')


    def fake(self, new_bjd):
        '''Create a fake transit light curve, using an input BJD array.'''

        # make an empty dictionary
        dict = {}

        # populate it with interpolated values
        dict['bjd'] = new_bjd
        dict['flux'] = np.interp(new_bjd, self.TLC.bjd, self.TLC.flux)
        dict['uncertainty'] = np.interp(new_bjd, self.TLC.bjd, self.TLC.uncertainty)

        # loop over the existing external variables, and populate them too
        for evkey in self.TLC.externalvariables.keys():
            #interpolator = scipy.interpolate.interp1d(self.TLC.bjd, #self.TLC.externalvariables[evkey])

            #dict[evkey] = interpolator(new_bjd)
            ok = self.bad == False
            dict[evkey] = np.interp(new_bjd, self.TLC.bjd[ok], self.TLC.externalvariables[evkey][ok])

        # create the fake TLC
        return TLC(left=self.left, right=self.right, directory=self.directory + 'fake/', isfake=True, **dict)

    def LightcurvePlots(self):
        '''A quick tool to plot what the light curve (and external variables) looks like.'''

        # set up the phase/unphased lightcurve plots
        self.setupLightcurvePlots()

        # create smoothed TLC structures, so the modeling will work
        self.TM.smooth_phased_tlc = self.fake( np.linspace(-self.TM.planet.period.value/2.0 + self.TM.planet.t0.value + 0.01, self.TM.planet.period.value/2.0 + self.TM.planet.t0.value-0.01, 10000))
        self.TM.smooth_unphased_tlc = self.fake(np.linspace(np.min(self.TLC.bjd), np.max(self.TLC.bjd), 1000))



        kw = {'marker':'.', 'color':self.colors['points'], 'alpha':0.5, 'linewidth':0, 'marker':'o', 'markeredgecolor':self.colors['points'], 'markersize':6}
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

        plt.draw()

    def residuals(self):
        return self.flux/self.TM.model() - 1

    def residualsexceptfor(self, exception=None, modeltoo=False):

        if exception is not None:
            original = {}
            instincluding = self.TM.instrument_model()
            for power in [1,2]:
                key = exception + '_tothe{0}'.format(power)
                try:
                    original[key] = self.TM.instrument.__dict__[key].value + 0.0
                    self.TM.instrument.__dict__[key].value = 0.0
                except KeyError:
                    pass
            ref = self.flux/self.TM.model() - 1
            instwithout = self.TM.instrument_model()
            for key in original.keys():
                self.TM.instrument.__dict__[key].value = original[key]

            if modeltoo:
                inst = instincluding - instwithout
                return ref, inst
            return ref

    def instrumental(self):
        about_instrument = self.flux/self.TM.planet_model()
        return  about_instrument/np.median(self.TM.instrument_model()) - 1

    def timefrommidtransit(self):
        return self.TM.planet.timefrommidtransit(self.bjd)

    def corrected(self):
        return self.flux/self.TM.instrument_model()


    def ExternalsMatrixPlots(self):
        new = dict(flux=self.instrumental(), time=self.timefrommidtransit(), residuals=self.residuals(), **self.externalvariables)
        keys = ['time']
        evkeys = self.externalvariables.keys()
        evkeys.remove('ok')
        evkeys.sort()
        keys.extend(evkeys)
        keys.append('flux')
        keys.append('residuals')
        self.speak('plotting instrumental residuals against a matrix of external variables')
        pdf = PDF(samples=new)
        pdf.plot(keys=keys)
        filename = self.TM.directory + 'matrixofexternals.png'
        self.speak('saving matrix of external variables to {0}'.format(filename))
        plt.savefig(filename)

    def setupSmooth(self):
        self.TM.smooth_phased_tlc = self.fake( np.linspace(-self.TM.planet.period.value/2.0 + self.TM.planet.t0.value + 0.01, self.TM.planet.period.value/2.0 + self.TM.planet.t0.value-0.01, 10000))


    def DiagnosticsPlots(self, noiseassumedforplotting=0.001, directory=None):
        '''A quick tool to plot what the light curve (and external variables) looks like.'''

        self.speak('plotting light curves with diagnostics')
        self.setupDiagnostics()
        ppm=1e6
        # create smoothed TLC structures, so the modeling will work
        self.TM.smooth_phased_tlc = self.fake( np.linspace(-self.TM.planet.period.value/2.0 + self.TM.planet.t0.value + 0.01, self.TM.planet.period.value/2.0 + self.TM.planet.t0.value-0.01, 10000))
        self.TM.smooth_unphased_tlc = self.fake(np.linspace(np.min(self.TLC.bjd), np.max(self.TLC.bjd), 10000))

        goodkw = {'color':self.colors['points'], 'alpha':1, 'linewidth':0, 'marker':'o', 'edgecolor':self.colors['points'], 's':20}
        badkw = {'color':self.colors['points'], 'alpha':0.25, 'linewidth':0, 'marker':'x', 'edgecolor':self.colors['points'], 's':20}
        time = self.TM.planet.timefrommidtransit(self.bjd)

        notok = self.bad
        for good in [False,True]:
            if good:
                ok = (self.bad == 0).nonzero()
                kw = goodkw
                if np.sum(ok) == 0:
                    return
            else:
                ok = (self.bad).nonzero()
                kw = badkw
            if np.sum(ok) == 0:
                continue
            self.ax_raw.scatter(time[ok], self.flux[ok], **kw)
            self.ax_corrected.scatter(time[ok], self.flux[ok]/self.TM.instrument_model()[ok], **kw)
            self.ax_residuals.scatter(time[ok], ppm*self.residuals()[ok], **kw)
            self.ax_instrument.scatter(time[ok], ppm*self.instrumental()[ok], **kw)
            self.points_correlations = {}
            kw['s'] = 10
            kw['alpha'] *= 0.5
            for k in self.ax_correlations.keys():
                #self.ax_correlations[k].plot(ppm*self.instrumental()[ok], self.externalvariables[k][ok], **kw)[0]
                #self.ax_timeseries[k].plot( time[ok], self.externalvariables[k][ok], **kw)

                # plot the external variables
                x = self.externalvariables[k]
                self.ax_correlations[k].scatter(ppm*(self.residualsexceptfor(k)[ok]), x[ok], **kw)
                if good:
                    binwidth = (np.nanmax(x[ok]) - np.nanmin(x[ok]))/10.0
                    bx, by, be = zachopy.oned.binto(x[ok], ppm*(self.residualsexceptfor(k)[ok]),  binwidth=binwidth, yuncertainty=self.uncertainty[ok], robust=False, sem=True)
                    self.ax_correlations[k].errorbar(by, bx, None, be, color='black', alpha=0.3, elinewidth=3, marker='o', linewidth=0, capthick=3)
                self.ax_timeseries[k].scatter( time[ok], self.externalvariables[k][ok], **kw)

                # set limits of plot windows
                if good:
                    self.ax_correlations[k].set_xlim(ppm*np.min(self.instrumental()[ok]), ppm*np.max(self.instrumental()[ok]))
                    self.ax_timeseries[k].set_ylim(np.min(self.externalvariables[k]), np.max(self.externalvariables[k]))

        ok = (self.bad == 0).nonzero()
        which = np.median(np.round((self.TM.smooth_unphased_tlc.bjd - self.TM.planet.t0.value)/self.TM.planet.period.value))
        modeltime = self.TM.smooth_unphased_tlc.bjd - self.TM.planet.t0.value - which*self.TM.planet.period.value
        assert(len(modeltime) == len(self.TM.model(self.TM.smooth_unphased_tlc)))
        kw = kw = {'color':self.colors['lines'], 'linewidth':3, 'alpha':1.0}

        self.line_raw = self.ax_raw.plot(modeltime, self.TM.model(self.TM.smooth_unphased_tlc), **kw)[0]
        self.line_corrected = self.ax_corrected.plot(modeltime, self.TM.planet_model(self.TM.smooth_unphased_tlc), **kw)[0]
        self.line_residuals = self.ax_residuals.plot(modeltime, ppm*np.zeros_like(modeltime), **kw)[0]
        justinstrument = (self.TM.instrument_model(self.TM.smooth_unphased_tlc)/np.median(self.TM.instrument_model()) - 1)
        self.line_instrument = self.ax_instrument.plot(modeltime, ppm*justinstrument, **kw)[0]


        # plot histograms of the residuals
        expectation = [0, ppm*np.mean(self.uncertainty)]
        zachopy.oned.plothistogram(ppm*self.instrumental()[ok], nbins=100, ax=self.ax_instrument_histogram,expectation =expectation , **kw)
        zachopy.oned.plothistogram(ppm*self.residuals()[ok], nbins=100, ax=self.ax_residuals_histogram, expectation =expectation , **kw)


        # plot binned RMS
        zachopy.oned.plotbinnedrms(self.residuals()[ok], ax=self.ax_binnedrms, yunits=1e6,  **kw)

        # plot the ACF
        zachopy.oned.plotautocorrelation(self.residuals()[ok], ax =self.ax_acf,  **kw)


        # print text about the fit


        self.line_correlations = {}
        kw['alpha'] = 0.5
        kw['linewidth'] = 1
        #for k in self.ax_correlations.keys():
            #self.line_correlations[k] = self.ax_correlations[k].plot(ppm*justinstrument, self.TM.smooth_unphased_tlc.externalvariables[k],  **kw)[0]
        #assert(np.std(ppm*justinstrument)>1)

        nsigma = 5
        if noiseassumedforplotting is None:
            self.noiseassumedforplotting = np.mean(self.uncertainty)
        else:
            self.noiseassumedforplotting = noiseassumedforplotting
        if self.noiseassumedforplotting is not None:
            scale = self.noiseassumedforplotting*nsigma
        else:
            scale = nsigma*np.mean(self.uncertainty[ok])
        ppm = 1e6
        self.ax_residuals.set_ylim(-scale*ppm, scale*ppm)
        self.ax_instrument.set_ylim(-scale*ppm, scale*ppm)

        self.ax_residuals.set_xlim(np.min(time), np.max(time))
        buffer = 0.0075
        self.ax_corrected.set_ylim(1.0 - self.TM.planet.depth - buffer, 1.0 + buffer)
        self.ax_raw.set_ylim(1.0 - self.TM.planet.depth - buffer*3, 1.0 + buffer*3)

        plt.draw()
        if directory is None:
            directory = self.TM.directory

        filename = directory + self.name.translate(None, '!@#$%^&*()<>') + '_lightcurveDiagnostics.pdf'
        self.speak('saving light curve diagnostic plot to {0}'.format(filename))
        plt.savefig(filename)

    def linkModel(self, tm):
        self.TM = tm
        self.TM.TLC = self
        self.TLC = self
        self.TM.TM = self.TM

    def gp_points(self):
        ok = self.ok

        # a dictionary of ok data points, in various stages of correction
        d = {}
        d['bjd'] = self.bjd[ok]
        d['t'] = self.TM.planet.timefrommidtransit(self.bjd)[ok]
        d['raw'] = self.flux[ok]
        d['cleaned'] = self.corrected()[ok]
        d['residuals'] = self.residuals()[ok]

        lnp = self.gp_lnlikelihood()
        mean_wiggle = self.gp.predict(d['residuals'], d['bjd'], mean_only=True)
        d['cleanedfromgp'] = d['cleaned'] - mean_wiggle
        d['residualsfromgp'] = d['residuals'] - mean_wiggle
        return d

    def points(self):
        ok = self.ok

        # a dictionary of ok data points, in various stages of correction
        d = {}
        d['bjd'] = self.bjd[ok]
        d['t'] = self.TM.planet.timefrommidtransit(self.bjd)[ok]
        d['raw'] = self.flux[ok]
        d['cleaned'] = self.corrected()[ok]
        d['residuals'] = self.residuals()[ok]
        return d

    def gp_lines(self, mean=True, resolution=300):
        '''if mean=True, will return the mean of the GP prediction; otherwise, will sample from it'''
        # make sure a smooth fake TLC is set up
        try:
            self.smoothed
        except AttributeError:
            self.smoothed = self.fake(np.linspace(np.min(self.bjd), np.max(self.bjd), resolution))

        # a dictionary of ok data points, in various stages of correction
        d = {}
        d['bjd'] = self.smoothed.bjd
        d['t'] = self.TM.planet.timefrommidtransit(d['bjd'])
        d['residualsfromgp'] = np.zeros_like(d['bjd'])
        d['cleanedfromgp'] = self.TM.planet_model(tlc=self.smoothed)

        ok = self.ok
        lnp = self.gp_lnlikelihood()
        if mean:
            wiggle, cov = self.gp.predict(self.residuals()[ok], d['bjd'])
            i = np.arange(len(wiggle))
            d['residualsstd'] = np.sqrt(cov[i,i])
        else:
            wiggle = self.gp.sample_conditional(self.residuals()[ok], d['bjd'])

        d['residuals'] = wiggle
        d['cleaned'] = d['cleanedfromgp'] + wiggle
        d['raw'] = self.TM.model(tlc=self.smoothed) + wiggle

        return d

    def lines(self):
        # make sure a smooth fake TLC is set up
        try:
            self.smoothed
        except AttributeError:
            self.smoothed = self.fake(np.linspace(np.min(self.bjd), np.max(self.bjd), 1000))

        # a dictionary of ok data points, in various stages of correction
        d = {}
        d['bjd'] = self.smoothed.bjd
        d['t'] = self.TM.planet.timefrommidtransit(d['bjd'])
        d['residualsfromgp'] = np.zeros_like(d['bjd'])
        d['cleanedfromgp'] = self.TM.planet_model(tlc=self.smoothed)

        wiggle = 0
        d['residuals'] = wiggle
        d['cleaned'] = d['cleanedfromgp'] + wiggle
        d['raw'] = self.TM.model(tlc=self.smoothed) + wiggle

        return d

    def splitIntoEpochs(self, planet=None, buffer=5, thresholdInTransit=1, thresholdOutOfTransit=1, newdirectory=None):
        assert(planet is not None)

        inTransit = np.abs(planet.timefrommidtransit(self.bjd)) < planet.duration/2.0
        epochNumbers = planet.thisepoch(self.bjd)
        uniqueEpochNumbers = np.unique(epochNumbers)

        newTLCs = []
        for u in uniqueEpochNumbers:
            inThisEpoch = epochNumbers == u
            nInThisTransit = np.sum(inTransit*inThisEpoch)
            nearThisTransit = inThisEpoch*(np.abs(planet.timefrommidtransit(self.bjd)) < (0.5 + buffer)*planet.duration)

            nNearThisTransit = np.sum(nearThisTransit)
            if (nInThisTransit > thresholdInTransit)&(nNearThisTransit > (thresholdOutOfTransit+thresholdOutOfTransit)):
                self.speak('creating a new TLC at epoch {0} with {1} data points'.format(u, nNearThisTransit))
                ok = nearThisTransit
                ev = {}
                for k in self.externalvariables.keys():
                    ev[k] = self.externalvariables[k][ok]
                ev['bad'] = self.bad[ok]
                if newdirectory is None:
                    newdirectory = self.directory
                newTLC = TLC(self.bjd[ok], self.flux[ok], self.uncertainty[ok],
                                left=self.left, right=self.right,
                                directory = newdirectory + '{0:.0f}/'.format(u), color=self.color,
                                epoch=u, telescope=self.telescope, name=None, remake=True,  **ev)
                newTLCs.append(newTLC)
            else:
                self.speak("there aren't enough data on epoch {0} to be worth while".format(u))
        return newTLCs

    def __repr__(self):
        return '<TLC|{name}|{telescope}|N={n} good data>'.format(name=self.name, telescope=self.telescope, epoch=self.epoch, n=np.sum(self.bad == False))

def demo():
    planet = Planet(J=0.00, rp_over_rs=0.1, rsum_over_a=1.0/20.0, cosi=0.000, q=0.000, period=1.58, t0=2456000.0, esinw=0.0, ecosw=0.0)
    star = Star(u1=0.3, u2=0.3, gd=0.32, albedo=0)
    tm = TM(planet=planet, star=star)
    (p1, p4, s1, s4) = planet.contacts()
    print planet.contacts()
    duration = (p4 - p1 + 1)*planet.period.value
    t = np.linspace(planet.t0.value - duration, planet.t0.value + duration, 1000)
    uncertainty = 0.003*np.ones_like(t)
    tlc = TLC(t, tm.model(t) + np.random.normal(len(t))*uncertainty, uncertainty)
    tlc.plot()
