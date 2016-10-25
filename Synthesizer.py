from imports import *
import transit.PDF as PDF
import transit
##@profile



def makeCode(d):
    '''convert a parameter dictionary into a string label'''
    return '{name}@{telescope}@{epoch}'.format(**d)


class Synthesizer(Talker):
    '''combine multiple TLCs (and maybe RVCs) together, a fit them simultaneously.'''

    def __init__(self, tlcs=[], rvcs=[], **kwargs):
        '''synthesizes a list of transit light curves (and their models),
            allowing them to be fit jointly.'''

        # setup talker
        Talker.__init__(self)

        # define the TLCs and the TMs, in the same order
        self.tlcs = tlcs
        self.rvcs = rvcs
        self.data = []

        # define self.data as the combination of all LCs and RVCs
        try:
            self.data.extend(self.tlcs)
        except TypeError:
            pass
        try:
            self.data.extend(self.rvcs)
        except TypeError:
            pass

        # keep track of the list of transit models, associated with each data
        self.tms = [thing.TM for thing in self.data]

        # determine the free parameters of the fit
        self.findParameters()

    @property
    def n(self):
        '''the total number of datapoints (both transit and RV)'''
        return self.nlc + self.nrv

    @property
    def nlc(self):
        '''the number of data-points in the light curves'''

        # loop over all TLCs
        tot = 0
        for tlc in self.tlcs:
            tot += len(tlc.bjd)
        return tot

    @property
    def nrv(self):
        '''the number of data-points in the radial velocity curves'''

        # loop over all RVCs
        tot = 0
        for rvc in self.rvcs:
            tot += len(rvc.bjd)
        return tot

    @property
    def m(self):
        '''the number of parameters that are going to be fit'''
        return len(self.parameters)

    def registerParameter(self, parameter, data):
        '''add a new parameter to the list this synthesizer cares about'''
        this = dict(name=parameter.name,
                    telescope=data.telescope,
                    epoch=data.epoch,
                    parameter=[parameter])
        this['code'] = makeCode(this)
        self.parameters.append(this)

    def popParameter(self, code):
        '''remove a parameter from the synthesizer's list'''
        codes = [par['code'] for par in self.parameters]
        return self.parameters.pop(codes.index(code))

    def findParameters(self):
        '''Search through all the TMs, and find all floating parameters. Merge
        them together into a big list with a fixed order, which can be used to
        set parameters when called by a fitter like MPFIT or emcee.'''

        # loop over data chunks (of either TLC or RVC type)
        self.parameters = []
        for i in range(len(self.data)):
            tm, data = self.tms[i], self.data[i]
            for p in tm.parameters:
                self.registerParameter(p, data)
                this = dict(name=p.name,                # general name of par.
                            telescope=data.telescope,    # the telescope
                            epoch=data.epoch,            # the epoch
                            parameter=[p])              # list of 1 or more
                this['code'] = makeCode(this)           # set the code

        self.speak('the list of all parameters is:')
        for i in range(len(self.parameters)):
            p = self.parameters[i]
            self.speak('   {0} -- {1}'.format(i, p['code']))

    def printParameters(self):
        '''Loop through and print all the parameters.'''
        for i in range(len(self.parameters)):
            p = self.parameters[i]
            n = len(p['parameter'])

            #if n == 1:
            self.speak("{i}: {code} {par}".format(i=i, code=p['code'], par=p['parameter'][0]))
            #else:
            #    self.speak("{i}: {code} (used {n} times) {par}".format(i=i, code=p['code'], par=oneParameter(p), n=n))

    def toArray(self):
        '''output parameters to (array of values, array of parinfos)'''
        par = np.array([p['parameter'][0].value for p in self.parameters])
        parinfo = [p['parameter'][0].parinfo for p in self.parameters]
        return par, parinfo

    def fromArray(self, a, verbose=False):
        '''set the parameters, from an array'''
        assert(len(a) == self.m)
        for i in range(self.m):
            for p in self.parameters[i]['parameter']:
                p.value = a[i]
                if verbose:
                    self.speak('{0} = {1}'.format(p.name, p.value))

    def tieAcrossEpochs(self, name):
        '''Merge together all parameters matching a particular name
        (and telescope and/or epoch), replacing their dictionary entry
        with one list of dictionary entries.'''

        # these parameters don't get touched at all
        original = []

        # how many parameters do we start with?
        n = len(self.parameters)

        # merge each telescope separately
        telescopes = np.unique([p['telescope'] for p in self.parameters if p['name'] == name])
        toMerge = {}
        for t in telescopes:
            toMerge[t] = []

        for i in range(n):
            p = self.parameters.pop()
            group = p['telescope']
            try:
                match = name == p['name']
            except KeyError:
                match = False

            if match:
                toMerge[group].append(p)
                #self.speak('adding {0} to "{1}" list'.format(p['code'], group))
            else:
                original.append(p)



        # keep the other parameters as they were
        self.parameters = original



        # add the merged parameters as a list of parameters
        epoch='global'
        for k in toMerge.keys():
            telescope=k

            # want to make sure we merge everything into one list, not nested dictionaries
            aslist = []
            for i in range(len(toMerge[k])):
                if type(toMerge[k][i]) == dict:
                    aslist.extend(toMerge[k][i]['parameter'])
                else:
                    aslist.append(toMerge[k][i])

            this = dict(name=name, telescope=telescope, epoch=epoch, parameter=aslist)
            this['code'] = makeCode(this)
            self.parameters.append(this)
            #self.speak('added {0} parameters under name {1}'.format(len(toMerge[k]), this['code']))

    def tieAcrossTelescopes(self, name):
        '''Merge together all parameters matching a particular name
        (and epoch), replacing their dictionary entry
        with one list of dictionary entries.'''

        # these parameters don't get touched at all
        original = []

        # how many parameters do we start with?
        n = len(self.parameters)

        # merge each telescope separately
        epochs = np.unique([p['epoch'] for p in self.parameters if p['name'] == name])
        toMerge = {}
        for e in epochs:
            toMerge[e] = []

        for i in range(n):
            p = self.parameters.pop()
            group = p['epoch']
            try:
                match = name == p['name']
            except KeyError:
                match = False

            if match:
                toMerge[group].append(p)
                #self.speak('adding {0} to "{1}" list'.format(p['code'], group))
            else:
                original.append(p)

        # keep the other parameters as they were
        self.parameters = original

        # add the merged parameters as a list of parameters
        telescope='global'
        for k in toMerge.keys():
            epoch=k
            # want to make sure we merge everything into one list, not nested dictionaries
            aslist = []
            for i in range(len(toMerge[k])):
                if type(toMerge[k][i]) == dict:
                    aslist.extend(toMerge[k][i]['parameter'])
                else:
                    aslist.append(toMerge[k][i])

            this = dict(name=name, telescope=telescope, epoch=epoch, parameter=aslist)
            this['code'] = makeCode(this)
            self.parameters.append(this)

    def deviatesTLC(self, p, fjac=None, plotting=False):
        '''Return the normalized deviates (an input for mpfit), collected from
            all the transit models.'''

        # update all the parameters inside the models, from the array
        self.fromArray(p)

        # set status to 0
        status = 0

        combined_deviates = []
        for tlc in self.tlcs:

            # if necessary, plot the light curve
            if plotting:
                tlc.LightcurvePlots()

            # weight only the good points
            ok = (tlc.bad == 0).nonzero()

            # calculate the deviates for this one TLC
            devs = (tlc.flux[ok] - tlc.TM.model()[ok])/tlc.effective_uncertainty[ok]

            # extend the list of deviates
            combined_deviates.extend(devs)

        # mpfit wants a list with the first element containing a status code
        return [status,combined_deviates]

    def deviatesRVC(self, p, fjac=None, plotting=False):
        '''Return the normalized deviates (an input for mpfit), collected from
            all the radial velocity curves.'''

        # update all the parameters inside the models, from the array
        self.fromArray(p)

        # set status to 0
        status = 0

        combined_deviates = []
        for rvc in self.rvcs:

            # weight only the good points
            ok = (rvc.bad == 0).nonzero()

            # calculate the deviates for this one TLC
            devs = (rvc.rv[ok] - rvc.TM.stellar_rv()[ok])/rvc.effective_uncertainty[ok]

            combined_deviates.extend(devs)

        # mpfit wants a list with the first element containing a status code
        return [status,combined_deviates]

    def deviates(self, p, fjac=None, plotting=False):
        '''return the (data) deviates, given a parameter array'''
        # create empty deviates
        dev = []

        # add the light curve deviates
        lcstatus, lightcurve = self.deviatesTLC(p, fjac=fjac, plotting=plotting)
        dev.extend(lightcurve)

        # add the radial velocity deviates
        rvstatus, rv = self.deviatesRVC(p, fjac=fjac, plotting=plotting)
        dev.extend(rv)

        # mpfit wants a list with the first element containing a status code
        return [(lcstatus | rvstatus), np.array(dev) ]

    def conjugatepriors_as_deviates(self):
        '''return conjugate priors, in the form of deviates'''
        status = 0
        effective_deviates = []

        # if specified, add a stellar density prior
        if (self.densityprior is not None):
            central, width = self.densityprior
            value = self.tms[0].planet.stellar_density
            prior = (value - central)/width
            effective_deviates.append(prior)

        # if specified, add period prior
        if (self.periodprior is not None):
            central, width = self.periodprior
            value = self.tms[0].planet.period.value
            prior = (value - central)/width
            effective_deviates.append(prior)

        # if specified, add t0 prior
        if (self.t0prior is not None):
            central, width = self.t0prior
            value = self.tms[0].planet.t0.value
            prior = (value - central)/width
            effective_deviates.append(prior)

        return [status, np.array(effective_deviates)]

    def deviates_including_priors(self, p, fjac=None, plotting=False):
        '''return the deviates, combining data and priors;
            (the equivalent of lnprob = lnlikelihood + lnprior)'''

        dev = []

        # get the data deviates
        datastatus, data = self.deviates(p, fjac=fjac, plotting=plotting)
        dev.extend(data)
        # get the conjugate prior deviates
        cpstatus, cp = self.conjugatepriors_as_deviates()
        dev.extend(cp)

        return [(datastatus | cpstatus), np.array(dev) ]

    def lnprob(self, p):
        """return the log posterior for a given parameter set"""

        try:
            self.budget['iteration'] += 1
        except (AttributeError,KeyError):
            self.budget = {}
            self.budget['iteration'] = 0

        lnp = 0.0


        # update all the parameters inside the models
        self.fromArray(p)

        # calculate the priors (including parameter region constraints)
        lnprior = self.lnprior()
        lnp += lnprior

        # if the priors are already zero, skip the likelihood calculation
        if np.isfinite(lnp) == False:
            return lnp

        lnlikelihood = self.lnlikelihood(p)
        lnp += lnlikelihood


        key = 'lnprob'
        self.budget[key] = lnp
        return lnp

    def lnlikelihood(self, p):
        '''return the lnlikelihood, using whichever method is set'''


        # calculate the likelihood, using the appropriate method
        if self.likelihoodtype == 'white':
            lnlikely = self.white_lnlikelihood
        elif self.likelihoodtype == 'red_beta':
            lnlikely = self.white_lnlikelihood
            # (same as above, but effective_uncertainty includes beta > 1)
        elif self.likelihoodtype == 'red_gp':
            lnlikely = self.gp_lnlikelihood
        lnlikelihood = lnlikely(p)

        # add additional terms
        lnlikelihood += self.additionalterms()

        return lnlikelihood

    def lnprior(self):
        '''return the prior probability of the current parameter set.'''

        # start at 0.0
        lnp = 0

        # return -np.inf if outside limits
        for p in self.parameters:
            par = p['parameter'][0]
            if par.value < par.limits[0]:
                key = '{0} out of limits'.format(par.name)
                self.budget[key] = -np.inf
                return -np.inf
            if par.value > par.limits[1]:
                key = '{0} out of limits'.format(par.name)
                self.budget[key] = -np.inf
                return -np.inf

        # return -np.inf if eccentricity is non-physical
        if (self.tms[0].planet.ecosw.value**2 + self.tms[0].planet.esinw.value**2) > 1.0:
            key = 'eccentricity out of limits'
            self.budget[key] = -np.inf
            return -np.inf

        # incorporate density prior, if need be
        try:
            assert(self.densityprior is not None)
            # pull out the prior parameters
            central, width = self.densityprior


            # calculate the density prior
            prior = -0.5*((self.tms[0].planet.stellar_density - central)/width)**2

            self.speak("{central}, {width}, {prior}".format(**locals()))
            # add it to the lnprior

            key = 'density prior'
            self.budget[key] = prior

            lnp +=  prior
        except (AttributeError,AssertionError):
            pass

        # incorporate a prior on the period, if need be
        try:
            assert(self.periodprior is not None)
            # pull out the prior parameters
            central, width = self.periodprior

            # calculate prior
            value = self.tms[0].planet.period.value
            prior = -0.5*((value - central)/width)**2

            key = 'period prior'
            self.budget[key] = prior

            # append it
            lnp += prior
        except (AttributeError,AssertionError):
            pass

        # incorporate a prior on t0, if need be
        try:
            assert(self.t0prior is not None)
            # pull out the prior parameters
            central, width = self.t0prior

            # calculate prior
            value = self.tms[0].planet.t0.value
            prior = -0.5*((value - central)/width)**2

            key = 't0 prior'
            self.budget[key] = prior

            # append it
            lnp += prior
        except (AttributeError,AssertionError):
            pass

        try:
            assert(self.eprior is not None)
            if (self.tms[0].planet.ecosw.fixed == False) or (self.tms[0].planet.esinw.fixed == False):
                # pull out the prior parameters
                a, b = self.eprior
                value = self.tms[0].planet.e
                #####

                key = 'eccentricity prior'
                self.budget[key] = prior

                prior = np.log(value**(a-1)*(1-value)**(b-1)/scipy.special.beta(a,b))
                lnp += prior
        except (AttributeError,AssertionError):
            pass

        assert(np.isfinite(lnp))
        return lnp

    def additionalterms(self):

        extra = 0.0
        # if using additative jitter for the RV's, need penalty
        jitterconstraint = 0.0
        for rvc in self.rvcs:
            # weight only the good points
            ok = (rvc.bad == 0).nonzero()
            thisjitter = np.sum(np.log(1.0/rvc.effective_uncertainty[ok]))
            key = 'non-chisq jitter term from {0}'.format(rvc.name)
            self.budget[key] = thisjitter
            jitterconstraint += thisjitter
        extra += jitterconstraint

        return extra

    def white_lnlikelihood(self,p):
        '''return the white lnlikelihood for a parameter set'''

        # calculate the chisq (note "effective_uncertainty" is used in deviates)
        chisq = np.sum(self.deviates(p)[-1]**2)

        # calculate the lnlikelihood from the chisq
        lnlikelihood = - chisq/2.0

        # don't let infinitely bad likelihoods break code
        if np.isfinite(lnlikelihood) == False:
            lnlikelihood = -np.inf

        key = 'white lnlikelihood (-chisq/2)'
        self.budget[key] = lnlikelihood

        return lnlikelihood

    def gp_lnlikelihood(self, p):
        '''return the GP lnlikelihood for a parameter set'''

        # loop over the tlcs
        lnlikelihood = 0
        for tlc in self.tlcs:
            # the TLC should already know what its hyperparameters are
            if np.sum(tlc.ok) > 0:
                lnlikelihood += tlc.gp_lnlikelihood()

        key = 'GP lnlikelihood'
        self.budget[key] = lnlikelihood

        return lnlikelihood

    @property
    def label(self):
        telescopes = np.unique([thing.telescope for thing in self.data])
        epochs = np.unique([thing.epoch for thing in self.data])

        return '{0}tlcs{1}telescopes{2}epochs{3}rvcs'.format(
                    len(self.tlcs),
                    len(telescopes),
                    len(epochs),
                    len(self.rvcs))

class Fit(Synthesizer):
    '''an object for fitting synthesized datasets (LM & MCMC inherit from it)'''
    def __init__(self,
                    tlcs=[],
                    rvcs=[],
                    directory=None,
                    likelihoodtype='white',
                    **kwargs):

        # keep track of which kinds of likelihood we're using
        self.likelihoodtype = likelihoodtype

        # create a basic synthesizer
        Synthesizer.__init__(self, tlcs=tlcs, rvcs=rvcs, **kwargs)

        # assign a directory for this fit
        if directory == None:
            base = 'synthesized/'
            zachopy.utils.mkdir(base)
            directory = base + self.label + '/'
        self.directory = directory
        zachopy.utils.mkdir(self.directory)

    def save(self):
        '''save this fit, so it be reloaded quickly next time.'''

        # create a directory, if need be
        zachopy.utils.mkdir(self.directory)

        # save the PDF, which contains all the information about the fit
        self.speak('saving LM fit to {0}'.format(self.directory))
        self.speak('  the PDF')
        self.pdf.save(self.directory + 'pdf.npy')

        # save a list dictionaries saying how the light curve has been modified
        modifications = []
        for tlc in self.tlcs:
            this = dict(name=tlc.name,
                        rescaling=tlc.rescaling,
                        chisq=tlc.chisq(),
                        ok=(tlc.bad == False))
            modifications.append(this)
        filename = self.directory + 'modifications.npy'
        np.save(filename, modifications)
        self.speak('saved modifications to {0}'.format(filename))

    def load(self):
        '''load a fit from a file, including parameters and LC modifcations'''

        self.speak('attempting to load from {0}'.format(self.directory))
        #self.speak('  the PDF')
        self.pdf = PDF.load(self.directory + 'pdf.npy')
        #self.speak('  the fitting notes')

        # pick up the covariance matrix
        self.covariance = self.pdf.covariance
        for i in range(self.m):
            match = np.array(self.pdf.names) == self.parameters[i]['code']
            for p in self.parameters[i]['parameter']:
                loaded = np.array(self.pdf.parameters)[match]
                assert(len(loaded) == 1)
                loaded = loaded[0]
                p.value = loaded.value
                p.uncertainty = loaded.uncertainty


        # incorporate all modifications that were made to the light curves
        # KLUDGE!
        if self.__class__.__name__ == 'MCMC':
            filename = self.directory + 'lm/modifications.npy'
        else:
            filename = self.directory + 'modifications.npy'
        modifications = np.load(filename)
        self.speak('loaded TLC modifications from {0}'.format(filename))
        assert(len(modifications) == len(self.tlcs))

        # loop over the modications, and apply them
        for i in range(len(self.tlcs)):
            this = modifications[i]
            tlc = self.tlcs[i]

            assert(tlc.name == this['name'])
            tlc.rescaling = this['rescaling']
            tlc.bad = this['ok'] == False
            assert(tlc.bad.shape == tlc.flux.shape)

        self.speak('applied the outlier clipping and rescaling')

    def setRescaling(self):
        '''after a fit has been performed (and applied!),
            rescale all light curves to have a chisq of 1'''

        # calculate a chisq rescaling for each tlc
        for tlc in self.tlcs:

            # select only the good points
            ok = tlc.ok

            # calculate the raw chisq of this individual light curve
            chisq = tlc.chisq()
            dof = len(ok) - self.m*(0.0 + len(ok))/self.n
            reduced_chisq = chisq/dof
            self.speak('{0} has a reduced chisq of {1} ({2} dof)'.format(
                        tlc.name, chisq, dof))


            tlc.rescaling = np.maximum(np.sqrt(reduced_chisq), 1)
            self.speak('rescaling original errors by {0} for next fit'.format(
                        tlc.rescaling))

            # if using the beta-approximation for red noise, calculate it
            if self.likelihoodtype == 'red_beta':
                tlc.beta = tlc.bestBeta()
                self.speak('set the beta to be {0} for next fit'.format(
                    tlc.beta))
            else:
                tlc.beta = 1.0

    def trainGP(self, tied=True, plot=False):
        '''after fit has been run, use residuals to train GP hyperparameters'''

        gplnarange = [-10,10]
        gplntaurange = [-10,10]
        # create parameters to describe the GP (time-only)
        for tlc in self.tlcs:
            # only add a GP for finite light curves
            if np.sum(tlc.ok) > 0:
                # pull out the instrument for this light curve
                instrument = tlc.TM.instrument
                # calculate the RMS, and use it as initial guess for "a"
                try:
                    instrument.gplna
                except AttributeError:
                    instrument.gplna = transit.Parameter(   'gplna', None)
                rms = np.std(tlc.residuals()[tlc.ok])
                instrument.gplna.float(0, gplnarange)

                # calculate cadence, and use it as initial guess for "tau"
                try:
                    instrument.gplntau
                except AttributeError:
                    instrument.gplntau = transit.Parameter( 'gplntau', None)
                instrument.gplntau.float(0,gplntaurange)

                self.registerParameter(instrument.gplna, tlc)
                self.registerParameter(instrument.gplntau, tlc)

        if tied:
            self.tieAcrossEpochs('gplna')
            self.tieAcrossEpochs('gplntau')

        # temporarily set the GP parameters to be the only ones we care about
        codes = [p['code'] for p in self.parameters]
        self.originalparameters = self.parameters
        self.parameters = []
        for p in self.originalparameters:
            if (p['name'] == 'gplna'):
                self.parameters.append(p)
            if (p['name'] == 'gplntau'):
                self.parameters.append(p)

        # play with the likelihood
        self.likelihoodtype = 'red_gp'
        pinitial = self.toArray()[0]
        self.budget = {}
        self.lnprob(pinitial)
        print self.budget

        self.gp_inputs = []
        self.gp_outputs = []
        def nlnprob(p):
            lnprob = self.lnprob(p)
            self.gp_inputs.append(p)
            self.gp_outputs.append(lnprob)
            return -lnprob

        d = self.directory+'gp/'
        zachopy.utils.mkdir(d)
        storedvalues = d + 'gplnagplntau.txt'
        try:
            gplna, gplntau = np.loadtxt(storedvalues)
            gplna_map, gplna_unc = gplna
            gplntau_map, gplntau_unc = gplntau

            self.speak('loaded gp parameters from file')
        except IOError:
            # use nelder-mead to optimize (often fails)
            #scipy.optimize.minimize(nlnprob, pinitial, method='nelder-mead')
            ndim = len(pinitial)
            nwalkers = 20

            initialwalkers = np.zeros((nwalkers, ndim))
            initialwalkers[:,0] = np.random.uniform(*gplnarange, size=nwalkers)
            initialwalkers[:,1] = np.random.uniform(*gplntaurange, size=nwalkers)
            p0 = initialwalkers
            sampler = emcee.EnsembleSampler(nwalkers, ndim, self.lnprob)
            plt.figure('training GP')
            plt.cla()
            def plotsampler(sampler, color=None, **kw):
                self.gp_inputs = sampler.flatchain
                self.gp_outputs = sampler.flatlnprobability

                array = np.array(self.gp_inputs)
                plt.figure('training GP')

                if color == None:
                    color = self.gp_outputs
                plt.scatter(array[:,0], array[:,1],c=color, edgecolor='none', **kw)
                plt.xlabel('gplna')
                plt.ylabel('gplntau')
                plt.xlim(gplnarange)
                plt.ylim(gplntaurange)
                plt.title(self.tlcs[0].telescope)
                plt.draw()

            print("Running first burn-in...")
            p0, _, _ = sampler.run_mcmc(p0, 100)
            plotsampler(sampler, color='gray', alpha=0.3)
            p = sampler.flatchain[np.argmax(sampler.flatlnprobability)]
            p0 = [p + 1e-8 * np.random.randn(ndim) for i in xrange(nwalkers)]
            sampler.reset()

            print("Running second burn-in...")
            p0, _, _ = sampler.run_mcmc(p0, 100)
            plotsampler(sampler, color='blue', alpha=0.3)
            sampler.reset()

            print("Running production...")
            sampler.run_mcmc(p0, 100)
            plotsampler(sampler)


            plt.savefig(d + 'optimizing.png')
            samples = {}
            samples['gplna'] = sampler.chain[:,:,0]
            samples['gplntau'] = sampler.chain[:,:,1]
            samples['lnprob'] = sampler.lnprobability[:,:]
            self.gppdf = PDF.Sampled(samples=samples, summarize=True)
            self.gppdf.storeforhuman(output=d + 'gpparameters')

            gplna_map, gplntau_map = sampler.flatchain[np.argmax(sampler.flatlnprobability)]
            gplna_unc = 1.48*zachopy.oned.mad(sampler.flatchain[:,0])
            gplntau_unc = 1.48*zachopy.oned.mad(sampler.flatchain[:,1])

            gplna = [gplna_map, gplna_unc]
            gplntau = [gplntau_map, gplntau_unc]
            np.savetxt(storedvalues, (gplna, gplntau))




        self.hyperparameters = [gplna_map, gplntau_map]

        # apply the hyperparameters
        self.fromArray([gplna_map, gplntau_map])
        # then, bring all the parameters back in
        self.parameters = self.originalparameters
        # fix the values for the hyperparameters
        for tm in self.tms:
            try:
                tm.instrument.gplna.fixed = True
                tm.instrument.gplntau.fixed = True
            except AttributeError:
                pass
        self.speak("gplna={0}\ngplntau={1}".format(gplna, gplntau))
        self.speak('just trained the GP')


class LM(Fit):
    '''perform a fit using LM optimization. this can only handle
        likelihoodtype = 'white', conjugate priors only, and no RV jitter'''

    def __init__(self, tlcs=[], rvcs=[], label='', **kwargs):

        # initialize the basic Fit
        Fit.__init__(self, tlcs=tlcs, rvcs=rvcs, **kwargs)

        # specify the directory for the LM
        self.directory = self.directory + 'lm/'
        zachopy.utils.mkdir(self.directory)

    def fit(self,   plot=False,
                    quiet=False,
                    firstpass=True,
                    secondpass=False,
                    remake=False,
                    densityprior=None,
                    periodprior=None,
                    t0prior=None,
                    eprior=None,
                    **kwargs):
        '''use LM (mpfit) to find the maximum probability parameters,
            and a covariance matrix. options for first, second, third pass
            are used to loop through several times to reject outliers and
            rescale uncertainties'''

        # either reload this fit if it's already been performed, or recreate it
        try:
            assert(remake == False)
            self.load()
            return
        except (AssertionError,IOError):
            self.speak('could not load this fit, making it from scratch!')

        # if priors on physical parameters are being used, keep track of them
        self.densityprior=densityprior
        self.periodprior=periodprior
        self.t0prior=t0prior
        self.eprior=eprior


        self.speak('performing a fast LM fit')

        # pull the parameter array out of the list of dictionaries
        p0, parinfo = self.toArray()

        # perform the LM fit, to get best fit parameters and covariance matrix
        self.speak('running mpfit minimization')
        self.mpfitted = mpfit.mpfit(self.deviates_including_priors, p0,
                                    parinfo=parinfo, quiet=quiet)

        # set the parameters to their fitted values
        self.fromArray(self.mpfitted.params)

        # on the first pass, reject outliers and rescale the errors
        # (roughly, just to set relative weights of datasets)
        if firstpass:

            # first, clip the outliers (so rescaling will be more meaningful)
            for tlc in self.tlcs:

                # select only the good points
                ok = tlc.ok

                # where are the residuals beyond 3 sigma?
                outlierthreshold = 3.0
                r = tlc.residuals()[ok]
                sigma = 1.48*zachopy.oned.mad(r)
                outlier = (np.abs(r) > outlierthreshold*sigma)
                tlc.bad[ok] = tlc.bad[ok] | (tlc.flags['outlier']*outlier)
                self.speak("identified {0} new points as bad"
                            "refitting without them".format(np.sum(outlier)))

            # rescaling the uncertainties
            self.setRescaling()

            # refit, after the outliers have been rejected, on second pass
            self.fit(plot=plot,
                    remake=remake, quiet=quiet,
                    firstpass=False, secondpass=True,
                    densityprior=self.densityprior,
                    periodprior=self.periodprior,
                    t0prior=self.t0prior,
                    **kwargs)
            return

        # on second pass, simply reset the rescalings
        # (one fit has already been done with ~correct weights)
        if secondpass:
            # rescaling the uncertainties
            self.setRescaling()

            # refit, after the outliers have been rejected
            self.fit(remake=remake, plot=plot, quiet=quiet,
                    firstpass=False, secondpass=False,
                    densityprior=self.densityprior,
                    periodprior=self.periodprior,
                    t0prior=self.t0prior, **kwargs)
            return

        # store the covariance matrix and 1D uncertainties on the parameters
        self.covariance = self.mpfitted.covar
        for i in range(self.m):
            for p in self.parameters[i]['parameter']:
                p.uncertainty = np.sqrt(self.covariance[i,i])

        # include all parameters in the covariance matrix, even with 0 unc.
        interesting = np.arange(self.m)

        # create a PDF structure out of this covariance matrix
        self.pdf = PDF.MVG(
                    names=np.array([p['code'] for p in self.parameters]),
                    parameters=np.array([p['parameter'][0] for p in np.array(self.parameters)[interesting]]),
                    covariance=self.covariance[interesting,:][:,interesting])

        # print the results
        self.pdf.printParameters()

        # save this fit
        self.save()


class MCMC(Fit):
    '''perform a fit using MCMC exploration. this object can handle
        likelihoodtype = 'white' | 'red_gp' | 'red_beta', any priors, and RV jitter'''

    def __init__(self, tlcs=[], rvcs=[], likelihoodtype='red_gp', **kwargs):
        Fit.__init__(self, tlcs=tlcs, rvcs=rvcs,
                        likelihoodtype=likelihoodtype, **kwargs)
        self.directory = self.directory + 'mcmc/'
        zachopy.utils.mkdir(self.directory)

    def fit(self,
        nburnin=500, ninference=1000, nwalkers=500, nleap=20, npreburnin=100,
        broad=True, ldpriors=True, fromcovariance=True,
        densityprior=None, periodprior=None, t0prior=None, eprior=None,
        plot=True, interactive=False, remake=False, updates=10, **kwargs):
        '''Use MCMC (with the emcee) to sample from the parameter probability distribution.'''

        # keep track of the priors
        self.densityprior=densityprior
        self.periodprior=periodprior
        self.t0prior=t0prior
        self.eprior=eprior

        # either load the MCMC fit, or run it from scratch
        self.speak('running an MCMC fit')
        try:
            assert(remake==False)
            self.load()
            return
        except (IOError,AssertionError):
            self.speak('could not load this MCMC fit, remaking it from scratch')

        # do an LM minization to clip outliers and rescale uncertainties
        if len(self.tlcs) > 0:    # create a LM fit to initialize the outlier rejection and rescaling

            # ignore the GP in the LM fit
            for tm in self.tms:
                try:
                    tm.instrument.gplna.fixed = True
                    tm.instrument.gplntau.fixed = True
                except AttributeError:
                    pass

            # set up the LM fit, using a white assumption
            self.lm = LM(self.tlcs, directory=self.directory, likelihoodtype='white')
            # set LM fit's (possibly linked) parameters to be same as this one
            self.lm.parameters = self.parameters
            # run (or load) the fit
            self.lm.fit(remake=remake, densityprior=self.densityprior, periodprior=self.periodprior, t0prior=self.t0prior, quiet=False)

            # take care of the kludge that left the GP out of the fit
            '''if self.likelihoodtype=='red_gp':
                self.lm.trainGP()

                gplna, gplntau = self.lm.hyperparameters
                for tm in self.tms:
                    try:
                        tm.instrument.gplna.fixed = True
                        tm.instrument.gplna.value = gplna
                        tm.instrument.gplntau.fixed = True
                        tm.instrument.gplntau.value = gplntau
                    except AttributeError:
                        pass'''


        # pull the parameter array out of the list of dictionaries
        p0, parinfo = self.toArray()
        nparameters = len(p0)

        # setup the initial walker positions
        self.speak('initializing {nwalkers} for each of the {nparameters}'.format(**locals()))
        initialwalkers = np.zeros((nwalkers, nparameters))

        # loop over the parameters
        for i in range(nparameters):
            # parameter = self.parameters[i]['parameter'][0]
            parameter = self.parameters[i]['parameter'][0]
            # take the initial uncertainties, and start around them
            if fromcovariance:
                nsigma=5
                #parameter = self.lm.parameters[i]['parameter'][0]
                bottom = parameter.value - nsigma*parameter.uncertainty
                top = parameter.value + nsigma*parameter.uncertainty
            else:

                bottom = np.min(parameter.limits)
                top = np.max(parameter.limits)

            initialwalkers[:,i] = np.random.uniform(bottom, top, nwalkers)
            self.speak('  {parameter.name} picked from uniform distribution spanning {bottom} to {top}'.format(**locals()))

        # set up the emcee sampler
        self.sampler = emcee.EnsembleSampler(nwalkers, nparameters, self.lnprob)
        self.names = [p['code'] for p in self.parameters]

        # add names to the sampler (for plotting progress, if desired)
        self.sampler.addLabels(self.names)

        # run a burn in step, and then reset
        done, saved = False, False
        burnt, inferred = False, False
        count = 0
        pos = initialwalkers
        index = 0

        # KLUDGE! to clean out the bad walker gunk
        self.speak("running {} preburn-in steps, with {} walkers".format(npreburnin, nwalkers))
        before = time.clock()
        pos, prob, state = self.sampler.run_mcmc(pos, npreburnin)
        after = time.clock()
        self.speak('it took {0} seconds'.format(after-before))
        self.speak()
        lastlnprob = self.sampler.lnprobability[:,-1]
        self.speak("the worst last walker was at lnprob = {}".format(np.min(lastlnprob)))
        limit = np.percentile(lastlnprob, 50)
        self.speak("reseting the walker positions to those above the last walkers' median ({})".format(limit))
        goodenough = self.sampler.flatlnprobability > limit
        options = self.sampler.flatchain[goodenough]
        pos = options[np.random.randint(0, len(options), nwalkers), :]


        # run the burn-in of the chain, creating plots along the way
        while (done and saved) == False:

            self.speak("{0}".format(datetime.datetime.now()))
            # run the chain for one "leap"



            if done == False:
                self.speak("running {0}-{1} of {2} burn-in steps, with {3} walkers".format(count, count+nleap, nburnin, nwalkers))
                #pos, prob, state = self.sampler.run_mcmc_with_progress(pos, nleap, updates=updates)
                before = time.clock()
                pos, prob, state = self.sampler.run_mcmc(pos, nleap)
                after = time.clock()
                self.speak('it took {0} seconds'.format(after-before))

            if done:
                plot=True
                self.speak('making plots and saving final chains')

            # output some summary plots
            if plot:
                # set up where plots should be saved
                if done:
                    d = self.directory + 'final/'
                    output = d + '{0:03.0f}_final'.format(index)
                else:
                    d = self.directory + 'status/'
                    if burnt:
                        output = d + '{0:03.0f}_inference'.format(index)
                    else:
                        output = d + '{0:03.0f}_burnin'.format(index)

                zachopy.utils.mkdir(d)
                plt.ioff()
                save = True

                nwalkers, nsteps, ndim = self.sampler.chain.shape

                if burnt:
                    which = nburnin + np.arange(nsteps - nburnin)
                    if done:
                        which = np.arange(ninference) + nburnin
                else:
                    which = nsteps/2 + np.arange(nsteps/2)

                # plot the trace of the parameters
                self.speak('creating plot of the parameter traces')
                before = time.clock()
                somekeys = [k for k in self.sampler.labels if 'global@global' in k]
                if len(somekeys) == 0:
                    somekeys = [k for k in self.sampler.labels if 'global' in k]
                self.sampler.HistoryPlot([0, nsteps],keys=somekeys)
                if save:
                    plt.savefig(output + '_parametertrace.png')
                after = time.clock()
                self.speak('it took {0} seconds'.format(after-before))

                self.speak('creating a PDF of the parameters')
                before = time.clock()
                samples = {}
                nc= len(self.sampler.flatlnprobability)
                for i in range(nparameters):
                    samples[self.names[i]] = self.sampler.flatchain[nc/2:,i]
                samples['lnprob'] = self.sampler.flatlnprobability[nc/2:]
                self.pdf = PDF.Sampled(samples=samples, summarize=done)
                after = time.clock()
                self.speak('it took {0} seconds'.format(after-before))

                # save summary
                self.speak('writing a text summary of the PDF')
                before = time.clock()
                s = ["{0}\n".format(datetime.datetime.now())]
                s.append("{nwalkers} walkers, {nsteps} steps, {ndim} parameters\n\n".format(**locals()))
                self.pdf.calculateUncertainties()
                s.extend(self.pdf.listParameters())
                txt = open(output + '_summary.txt', 'w')
                txt.write(''.join(s))
                txt.close()
                after = time.clock()
                self.speak('it took {0} seconds'.format(after-before))


                self.speak('saving the PDF to {0}'.format(self.directory))
                before = time.clock()
                self.pdf.save(self.directory + 'pdf.npy')
                self.pdf.export(self.directory + 'transitproperties.npy', keys=[k for k in self.pdf.names if 'global@global' in k])
                after = time.clock()
                self.speak('it took {0} seconds'.format(after-before))

                if done:
                    self.pdf.printParameters()
                    self.save()

                self.speak('creating a plot of the PDF')
                before = time.clock()
                keys = [n for n in self.pdf.names if 'global' in n]
                keys.append('lnprob')
                self.speak('{0}'.format(keys))
                figure = self.pdf.triangle(keys=keys, title=output,
                    plot_contours=True, plot_datapoints=False, plot_density=False,
                    alpha=0.5, show_titles=False)
                plt.savefig(output + '_parameterpdf.pdf')
                plt.close(figure)
                after = time.clock()
                self.speak('it took {0} seconds'.format(after-before))

                self.plotData(output=output)
                plt.close('all')
                #bla = self.input('bla')
                # plot
                '''for tlc in self.tlcs:
                    transit.QuicklookTransitPlot(tlc=tlc)
                    if save:
                        plt.savefig(output + '_lightcurves_{0}.pdf'.format(tlc.name.replace(',','_')))
                '''
                if done:
                    saved = True

            count += nleap
            index += 1
            burnt = count >= nburnin
            inferred = count >= (nburnin + ninference)
            done = burnt and inferred

        # set the parameter to their MAP values
        best = self.sampler.flatchain[np.argmax(self.sampler.flatlnprobability)]
        self.fromArray(best)

    def fromPDF(self, option='best'):
        self.speak('drawing {0} sample from PDF'.format(option))
        values = self.pdf.drawSample(keys=[p['code'] for p in self.parameters], option=option)
        self.fromArray(values)

    def plotData(self, output=None):
        if output is None:
            output = self.directory + 'data'
        #
        try:
            best = self.sampler.flatchain[np.argmax(self.sampler.flatlnprobability)]
        except AttributeError:
            orderednames = [p['code'] for p in self.parameters]
        self.fromPDF(option='best')

        if len(self.tlcs) > 0:
            self.speak('creating a plot of the light curves that have been fitted')
            before = time.clock()
            # is it the MEarth object?
            ismearth=False
            for tlc in self.tlcs:
                if 'MEarth' in tlc.telescope:
                    ismearth=True
            if ismearth:
                transit.IndividualPlots(tlcs=self.tlcs, synthesizer=self)
            else:
                xlim = self.tlcs[0].TM.planet.duration*2
                transit.IndividualPlots(tlcs=self.tlcs, synthesizer=self, telescopes=None, epochs=None, xlim=[-xlim, xlim], binsize=15./60/24, gskw=dict(top=0.9, bottom=0.2))
            plt.savefig(output + '_everything.pdf')
            after = time.clock()
            self.speak('it took {0} seconds'.format(after-before))

        if len(self.rvcs) > 0:
            self.fromPDF(option='best')

            self.speak('creating a plot of the RVs that have been fitted')
            before = time.clock()
            #ylim=[-15, 15]
            height = self.rvcs[0].TM.planet.semiamplitude.value*4*1e3
            ylim = [-height, height]
            kw = dict(ylim=ylim)
            p = transit.RVPhasedPlot(rvcs=self.rvcs, **kw)
            plt.savefig(output + '_rv.pdf')
            after = time.clock()
            self.speak('it took {0} seconds'.format(after-before))
