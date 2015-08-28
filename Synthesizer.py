from imports import *
import transit.PDF as PDF
import transit
##@profile



def makeCode(d):
    # convert a parameter dictionary into a string label
    return '{name}@{telescope}@{epoch}'.format(**d)


class Synthesizer(Talker):
    '''Combine multiple TLCs (and maybe RVCs) together, a fit them simultaneously.'''

    def __init__(self, tlcs=[], rvcs=[], **kwargs):
        '''Synthesizes a list of transit light curves (and their models),
            allowing them to be fit jointly.'''

        # setup talker
        Talker.__init__(self)

        # define the TLCs and the TMs, in the same order
        self.tlcs = tlcs
        self.rvcs = rvcs
        self.data = []
        try:
            self.data.extend(self.tlcs)
        except TypeError:
            pass
        try:
            self.data.extend(self.rvcs)
        except TypeError:
            pass

        self.tms = [thing.TM for thing in self.data]

        self.findParameters()

    @property
    def n(self):
        return self.nlc + self.nrv

    @property
    def nlc(self):
        # the number of data-points in the light curves
        tot = 0
        for tlc in self.tlcs:
            tot += len(tlc.bjd)
        return tot

    @property
    def nrv(self):
        # the number of data-points in the radial velocity curves
        tot = 0
        for rvc in self.rvcs:
            tot += len(rvc.bjd)
        return tot


    @property
    def m(self):
        # the number of parameters that are going to be fit
        return len(self.parameters)

    def findParameters(self):
        '''Search through all the TMs, and find all floating parameters. Merge
        them together into a big list with a fixed order, which can be used to
        set parameters when called by a fitter like MPFIT or emcee.'''

        self.parameters = []
        for i in range(len(self.data)):
            tm, data = self.tms[i], self.data[i]
            for p in tm.parameters:

                this = dict(name=p.name,                # general name of par.
                            telescope=data.telescope,    # the telescope
                            epoch=data.epoch,            # the epoch
                            parameter=[p])              # list of 1 or more
                this['code'] = makeCode(this)           # set the code
                self.parameters.append(this)

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

        par = np.array([p['parameter'][0].value for p in self.parameters])
        parinfo = [p['parameter'][0].parinfo for p in self.parameters]
        return par, parinfo

    def fromArray(self, a):
        assert(len(a) == self.m)
        for i in range(self.m):
            for p in self.parameters[i]['parameter']:
                p.value = a[i]


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

            combined_deviates.extend(devs)

        #assert(self.densityprior == False)
        #assert(self.densityprior is not None)
        if (self.densityprior is not None) & (len(combined_deviates) > 0):
            central, width = self.densityprior
            if len(self.tms) > 1:
                assert(self.tms[0].planet.stellar_density == self.tms[1].planet.stellar_density)
            density_prior = (self.tms[0].planet.stellar_density - central)/width
            combined_deviates.append(density_prior)
        #print density_prior
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
        dev = []
        lcstatus, lightcurve = self.deviatesTLC(p, fjac=fjac, plotting=plotting)
        dev.extend(lightcurve)
        rvstatus, rv = self.deviatesRVC(p, fjac=fjac, plotting=plotting)
        dev.extend(rv)


        if (self.periodprior is not None):
            central, width = self.periodprior
            value = self.tms[0].planet.period.value
            prior = (value - central)/width
            dev.append(prior)

        if (self.t0prior is not None):
            central, width = self.t0prior
            value = self.tms[0].planet.t0.value
            prior = (value - central)/width
            dev.append(prior)

        return [(lcstatus | rvstatus), np.array(dev) ]

    def gp_lnprob(self, p):

        # update all the parameters inside the models, from the array
        #   these include all hyperparameters too
        self.fromArray(p)

        for p in self.parameters:
            self.speak('{0} = {1}'.format(p['parameter'][0].name, p['parameter'][0].value))

        lnprior = self.lnprior()
        if np.isfinite(lnprior) == False:
            return -np.inf

        # loop over the tlcs
        lnlikelihood = 0
        for tlc in self.tlcs:
            # the TLC should already know what its hyperparameters are
            if np.sum(tlc.ok) > 0:
                lnlikelihood += tlc.gp_lnprob()

        lnprob = lnprior + lnlikelihood
        self.speak('lnlikelihood = {lnlikelihood:10.1f}, lnprior = {lnprior:10.1f}, lnprob = {lnprob:10.1f}'.format(**locals()))

        return lnprob

    def lnprior(self):

        lnp = 0

        for p in self.parameters:
            par = p['parameter'][0]
            if par.value < par.limits[0]:
                return -np.inf
            if par.value > par.limits[1]:
                return -np.inf

        # incorporate density prior, if need be
        if (self.densityprior is not None):
            # pull out the prior parameters
            central, width = self.densityprior

            # make sure multiple transit models agree on the calculated density
            if len(self.tms) > 1:
                assert(self.tms[0].planet.stellar_density == self.tms[1].planet.stellar_density)

            # calculate the density prior
            prior = -0.5*((self.tms[0].planet.stellar_density - central)/width)**2

            # add it to the lnprior
            lnp +=  prior

        # incorporate a prior on the period, if need be
        if (self.periodprior is not None):
            # pull out the prior parameters
            central, width = self.periodprior

            # calculate prior
            value = self.tms[0].planet.period.value
            prior = -0.5*((value - central)/width)**2

            # append it
            lnp += prior

        # incorporate a prior on t0, if need be
        if (self.t0prior is not None):
            # pull out the prior parameters
            central, width = self.t0prior

            # calculate prior
            value = self.tms[0].planet.t0.value
            prior = -0.5*((value - central)/width)**2

            # append it
            lnp += prior

        if (self.eprior is not None):
            # pull out the prior parameters
            a, b = self.eprior
            value = self.tms[0].planet.e
            #####

            prior = np.ln(value**(a-1)*(1-value)**(b-1)/scipy.special.beta(a,b))
            lnp += prior
        return lnp

    ##@profile
    def lnprob(self, p):
        """Return the log posterior, calculated from the deviates function (which may have included some conjugate Gaussian priors.)"""

        if self.gp:
            return self.gp_lnprob(p)

        chisq = np.sum(self.deviates(p)[-1]**2)
        #N = np.sum(self.TLC.bad == 0)

        # sum the deviates into a chisq-like thing
        #lnlikelihood = -N * np.log(self.instrument.rescaling.value) - chisq/self.instrument.rescaling.value**2



        lnlikelihood = - chisq/2.0#/self.instrument.rescaling.value**2


        # don't let infinitely bad likelihoods break code
        if np.isfinite(lnlikelihood) == False:
            lnlikelihood = -1e9

        # initialize an empty constraint, which could freak out if there's something bad about this fit
        constraints = 0.0



        # loop over the parameters
        for parameter in self.parameters:

            p = parameter['parameter'][0]
            # if a parameter is outside its allowed range, then make the constraint very strong!
            inside = (p.value < p.limits[1]) & (p.value > p.limits[0])
            if inside == False:
                constraints -= np.inf

        jitterconstraint = 0.0
        for rvc in self.rvcs:
            jitterconstraint += np.sum(np.log(1.0/rvc.effective_uncertainty))
            #self.speak('for jitter of {0}, penalty is {1}'.format(rvc.jitter, jitterconstraint))

        # return the constrained likelihood
        #print self.tms[0].planet.stellar_density, density_prior
        return lnlikelihood + constraints + jitterconstraint# + density_prior



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
    def __init__(self, tlcs=[], rvcs=[], directory=None, gp=False, **kwargs):
        self.gp = gp
        Synthesizer.__init__(self, tlcs=tlcs, rvcs=rvcs, **kwargs)

        if directory == None:
            base = 'synthesized/'
            zachopy.utils.mkdir(base)
            directory = base + self.label + '/'
        self.directory = directory
        zachopy.utils.mkdir(self.directory)

    def save(self):
        #Save this fit, so it be reloaded quickly next time.
        zachopy.utils.mkdir(self.directory)

        # save the PDF, which contains all the information about the fit
        self.speak('saving LM fit to {0}'.format(self.directory))
        self.speak('  the PDF')
        self.pdf.save(self.directory + 'pdf.npy')

        # save a list dictionaries saying how the light curve has been modified
        modifications = []
        for tlc in self.tlcs:
            print tlc
            this = dict(name=tlc.name,
                        rescaling=tlc.rescaling,
                        chisq=tlc.chisq(),
                        ok=tlc.bad == False)
            modifications.append(this)
        filename = self.directory + 'modifications.npy'
        np.save(filename, modifications)
        self.speak('saved modifications to {0}'.format(filename))

    def load(self):
        #Save this fit, so it be reloaded quickly next time.
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
                #self.speak('assigned {0} to {1}'.format(loaded, self.parameters[i]['code']))
                #self.speak('!!!!!!!!!!!!!!!!')

        # incorporate all modifications that were made to the light curves
        filename = self.directory + 'modifications.npy'
        modifications = np.load(filename)
        self.speak('loaded TLC modifications from {0}'.format(filename))
        assert(len(modifications) == len(self.tlcs))

        for i in range(len(self.tlcs)):
            this = modifications[i]
            tlc = self.tlcs[i]

            assert(tlc.name == this['name'])
            tlc.rescaling = this['rescaling']
            tlc.bad = this['ok'] == False
            assert(tlc.bad.shape == tlc.flux.shape)
            #assert(tlc.chisq() == this['chisq'])

        self.speak('applied the outlier clipping and rescaling')


    def setRescaling(self):
        # after a fit has been performed (and applied!), rescale all light curves to have a chisq of 1

        # calculate a chisq rescaling for each tlc
        for tlc in self.tlcs:

            # select only the good points
            ok = tlc.ok

            # calculate the raw chisq of this individual light curve
            chisq = tlc.chisq()
            dof = len(ok) - self.m*(0.0 + len(ok))/self.n
            reduced_chisq = chisq/dof
            self.speak('{0} has a reduced chisq of {1} ({2} dof)'.format(tlc.name, chisq, dof))


            tlc.rescaling = np.maximum(np.sqrt(reduced_chisq), 1)
            self.speak('rescaling original errors by {0} for next fit'.format(tlc.rescaling))

            if self.gp == False:
                tlc.beta = tlc.bestBeta()
                self.speak('set the beta to be {0} for next fit'.format(tlc.beta))


class LM(Fit):
    def __init__(self, tlcs=[], rvcs=[], label='', **kwargs):
        Fit.__init__(self, tlcs=tlcs, rvcs=rvcs, **kwargs)
        self.directory = self.directory + 'lm/'
        zachopy.utils.mkdir(self.directory)

    def fit(self, plot=False, quiet=False, firstpass=True, secondpass=False, remake=False, densityprior=None, periodprior=None, t0prior=None, eprior=None, **kwargs):
        '''Use LM (mpfit) to find the maximum probability parameters, and a covariance matrix.'''

        try:
            assert(remake == False)
            self.load()
            return
        except (AssertionError,IOError):
            self.speak('could not load this fit, making it from scratch!')

        self.densityprior=densityprior
        self.periodprior=periodprior
        self.t0prior=t0prior
        self.eprior=eprior


        self.speak('performing a fast LM fit')

        # pull the parameter array out of the list of dictionaries
        p0, parinfo = self.toArray()

        # perform the LM fit, to get best fit parameters and covariance matrix
        self.speak('running mpfit minimization')
        self.mpfitted = mpfit.mpfit(self.deviates, p0, parinfo=parinfo, quiet=quiet)

        # set the parameters to their fitted values
        self.fromArray(self.mpfitted.params)

        # on the first pass, reject outliers and rescale the errors (roughly, to set relative weight of datasets)
        if firstpass:

            # first, clip the outliers (so rescaling will be more meaningful)
            for tlc in self.tlcs:

                # select only the good points
                ok = tlc.ok

                # where are the residuals beyond 3 sigma?
                outlierthreshold = 3.0
                r = tlc.residuals()[ok]
                outlier = (np.abs(r) > outlierthreshold*1.48*zachopy.oned.mad(r))
                tlc.bad[ok] = tlc.bad[ok] | (tlc.flags['outlier']*outlier)
                self.speak("identified {0} new points as bad; refitting without them".format(np.sum(outlier)))

            # rescaling the uncertainties
            self.setRescaling()
            # refit, after the outliers have been rejected
            self.fit(plot=plot, remake=remake, quiet=quiet, firstpass=False, secondpass=True, densityprior=self.densityprior, periodprior=self.periodprior, t0prior=self.t0prior, **kwargs)
            return

        # on second pass, simply reset the rescalings, now that one fit has been done with ~correct weights
        if secondpass:
            # rescaling the uncertainties
            self.setRescaling()

            # refit, after the outliers have been rejected
            self.fit(remake=remake, plot=plot, quiet=quiet, firstpass=False, secondpass=False, densityprior=self.densityprior, periodprior=self.periodprior, t0prior=self.t0prior, **kwargs)
            return

        # store the covariance matrix of the fit, and the 1D uncertainties on the parameters
        self.covariance = self.mpfitted.covar
        for i in range(self.m):
            for p in self.parameters[i]['parameter']:
                p.uncertainty = np.sqrt(self.covariance[i,i])

        # pull out the parameters that actually varied and create a PDF object out of them
        # interesting = (self.covariance[range(self.m), range(self.m)] > 0).nonzero()[0]

        # include all parameters in the covariance matrix, even those with 0 uncertainty
        interesting = np.arange(self.m)

        # create a PDF structure out of this covariance matrix
        self.pdf = PDF.MVG(
                    names=np.array([p['code'] for p in self.parameters]),
                    parameters=np.array([p['parameter'][0] for p in np.array(self.parameters)[interesting]]),
                    covariance=self.covariance[interesting,:][:,interesting])

        # print the results
        self.pdf.printParameters()

        self.save()


class MCMC(Fit):
    def __init__(self, tlcs=[], rvcs=[], **kwargs):
        Fit.__init__(self, tlcs=tlcs, rvcs=rvcs, **kwargs)
        self.directory = self.directory + 'mcmc/'
        zachopy.utils.mkdir(self.directory)

    def fit(self,
        nburnin=1000, ninference=1000, nwalkers=500,
        broad=True, ldpriors=True, fromcovariance=True, densityprior=None, periodprior=None, t0prior=None, eprior=None,
        plot=True, interactive=False, remake=False, updates=10, **kwargs):
        '''Use MCMC (with the emcee) to sample from the parameter probability distribution.'''

        self.densityprior=densityprior
        self.periodprior=periodprior
        self.t0prior=t0prior
        self.eprior=eprior

        self.speak('running an MCMC fit')
        try:
            assert(remake==False)
            self.load()
            return
        except IOError:
            self.speak('could not load this MCMC fit, remaking it from scratch')

        #if fromcovariance:
        if len(self.tlcs) > 0:    # create a LM fit to initialize the outlier rejection and rescaling
            if self.gp:
                for tm in self.tms:
                    tm.instrument.gplna.fixed = True
                    tm.instrument.gplntau.fixed = True

            self.lm = LM(self.tlcs, directory=self.directory, gp=False)
            # set LM fit's (possibly linked) parameters to be same as this one
            self.lm.parameters = self.parameters
            # run (or load) the fit
            self.lm.fit(remake=remake, densityprior=self.densityprior, periodprior=self.periodprior, t0prior=self.t0prior, quiet=False)
            if self.gp:
                for tm in self.tms:
                    tm.instrument.gplna.fixed = False
                    tm.instrument.gplntau.fixed = False

            a = self.input('parameters okay?')

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
        burnt = False
        count = 0
        pos = initialwalkers
        while burnt == False:
            self.speak("running {0} burn-in steps, with {1} walkers.".format(nburnin, nwalkers))
            pos, prob, state = self.sampler.run_mcmc_with_progress(pos, nburnin, updates=updates)

            if interactive:
                plot = True
                output = None
                plt.ion()

            if plot:
                d = self.directory + 'status/'
                zachopy.utils.mkdir(d)
                output = d + '{0:03.0f}'.format(count/nburnin)
                plt.ioff()
                save = True

            if plot:
                # plot the chains

                self.sampler.HistoryPlot([count, count + nburnin],nmax=None)
                if save:
                    plt.savefig(output + '_parametertrace.pdf')

                # plot the PDF
                samples = {}
                for i in range(nparameters):
                    samples[self.names[i]] = self.sampler.flatchain[:,i]
                samples['lnprob'] = self.sampler.flatlnprobability[:]

                #if save:
                    #self.pdf = PDF.Sampled(samples=samples)
                    #self.pdf.printParameters()
                    #self.pdf.triangle(self.pdf.names)
                    #self.pdf.triangle(keys=self.pdf.names, title=output,
                            #plot_contours=True, plot_datapoints=False, plot_density=False,
                            #alpha=0.5, show_titles=False)
                    #plt.savefig(output + '_parameterpdf.pdf')

                # set the parameter to their MAP values
                best = self.sampler.flatchain[np.argmax(self.sampler.flatlnprobability)]
                self.fromArray(best)

                # plot
                for tlc in self.tlcs:
                    transit.QuicklookTransitPlot(tlc=tlc)
                    if save:
                        plt.savefig(output + '_lightcurves_{0}.pdf'.format(tlc.name.replace(',','_')))



            if interactive:
                answer = self.input('Do you think we have burned in?')
                if 'y' in answer:
                    burnt = True
            else:
                burnt = True

            count += nburnin

        # after the burn-in, reset the chain
        self.sampler.reset()

        # loop until satisfied with the inference samples
        self.speak('running for inference, using {0} steps and {1} walkers'.format(ninference, nwalkers))

        # start with the last set of walker positions from the burn in and run for realsies
        self.sampler.run_mcmc_with_progress(pos, ninference)

        # trim the chain to the okay values (should this be necessary?)
        # ok = self.sampler.flatlnprobability > (np.max(self.sampler.flatlnprobability) - 100)

        samples = {}
        for i in range(nparameters):
            samples[self.names[i]] = self.sampler.flatchain[:,i]
        samples['lnprob'] = self.sampler.flatlnprobability[:]

        # set the parameter to their MAP values
        best = self.sampler.flatchain[np.argmax(self.sampler.flatlnprobability)]
        self.fromArray(best)

        self.pdf = PDF.Sampled(samples=samples)
        self.pdf.printParameters()
        self.save()

            #if plot:
            #    self.model.TLC.DiagnosticsPlots(directory=self.directory)
