from zachopy.Talker import Talker
import matplotlib.pyplot as plt, numpy as np
import zachopy.borrowed.mpfit.mpfit as mpfit
import zachopy.oned
import pemcee as emcee
import transit.PDF as PDF
##@profile



def makeCode(d):
    return '{name}_{telescope}_{epoch}'.format(**d)



class Synthesizer(Talker):
    def __init__(self, tlcs, **kwargs):
        '''Synthesizes a list of transit light curves (and their models),
            allowing them to be fit jointly.'''

        # setup talker
        Talker.__init__(self)

        # define the TLCs and the TMs, in the same order
        self.tlcs = tlcs
        self.tms = [tlc.TM for tlc in self.tlcs]
        self.findParameters()

    @property
    def n(self):
        # the number of data-points
        tot = 0
        for tlc in self.tlcs:
            tot += len(tlc.bjd)
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
        for i in range(len(self.tlcs)):
            tm, tlc = self.tms[i], self.tlcs[i]
            for p in tm.parameters:

                this = dict(name=p.name,                # general name of par.
                            telescope=tlc.telescope,    # the telescope
                            epoch=tlc.epoch,            # the epoch
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
                match = name in p['name']
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
                match = name in p['name']
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

    def deviates(self, p, fjac=None, plotting=False):
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
            devs = (tlc.flux[ok] - tlc.TM.model()[ok])/tlc.uncertainty[ok]/tlc.rescaling

            combined_deviates.extend(devs)

        # mpfit wants a list with the first element containing a status code
        return [status,np.array(combined_deviates)]


class Fit(Synthesizer):
    def __init__(self, tlcs, **kwargs):
        Synthesizer.__init__(self, tlcs, **kwargs)

    def save(self):
        #Save this fit, so it be reloaded quickly next time.
        zachopy.utils.mkdir(self.directory)
        self.speak('saving LM fit to {0}'.format(self.directory))
        self.speak('  the PDF')
        self.pdf.save(self.directory + 'pdf.npy')
        self.speak('  the fitting notes')
        np.save(self.directory + 'fitting_notes.npy', self.notes)
        self.speak('  the best-fit model')
        self.model.save(self.directory)

    def load(self):
        #Save this fit, so it be reloaded quickly next time.
        self.speak('attempting to load from {0}'.format(self.directory))
        #self.speak('  the PDF')
        self.pdf = PDF.load(self.directory + 'pdf.npy')
        #self.speak('  the fitting notes')
        self.notes = np.load(self.directory + 'fitting_notes.npy')[()]
        #self.speak('  the best-fit model')
        self.model.load(self.directory)




class LM(Fit):
    def __init__(self, tlcs, **kwargs):
        Fit.__init__(self, tlcs, **kwargs)

    def fit(self, plot=False, quiet=False, firstpass=True, remake=False, **kwargs):
        '''Use LM (mpfit) to find the maximum probability parameters, and a covariance matrix.'''

        self.speak('performing a fast LM fit')

        p0, parinfo = self.toArray()

        # perform the LM fit, to get best fit parameters and covariance matrix
        self.speak('running mpfit minimization')
        self.mpfitted = mpfit.mpfit(self.deviates, p0, parinfo=parinfo, quiet=quiet)

        # set the parameters to their fitted values
        self.fromArray(self.mpfitted.params)

        # calculate a chisq rescaling for each tlc
        for tlc in self.tlcs:

            # select only the good points
            ok = tlc.ok

            chisq = np.sum((tlc.residuals()/tlc.uncertainty)[ok]**2)
            dof = len(ok) - self.m*(0.0 + len(ok))/self.n
            reduced_chisq = chisq/dof

            self.speak('{0} has a reduced chisq of {1} ({2} dof)'.format(tlc.name, chisq, dof))


            if firstpass:

                tlc.rescaling = np.maximum(np.sqrt(reduced_chisq), 1)
                self.speak('rescaling errors by {0} for next fit'.format(tlc.rescaling))
                # SHOULD PROBABLY SPLIT OUTLIER AND RESCALING TO SEPARATE STEPS

                # where are the residuals beyond 4 sigma?
                outlierthreshold = 3.0
                r = tlc.residuals()[ok]
                outlier = (np.abs(r) > outlierthreshold*1.48*zachopy.oned.mad(r))
                tlc.bad[ok] = tlc.bad[ok] | (tlc.flags['outlier']*outlier)
                self.speak("identified {0} new points as bad; refitting without them".format(np.sum(outlier)))

        if firstpass:
            # refit, after the outliers have been rejected
            self.fit(plot=plot, quiet=quiet, firstpass=False, **kwargs)
            return

        # store the covariance matrix of the fit, and the 1D uncertainties on the parameters
        self.covariance = self.mpfitted.covar
        for i in range(self.m):
            for p in self.parameters[i]['parameter']:
                p.uncertainty = np.sqrt(self.covariance[i,i])

        # pull out the parameters that actually varied and create a PDF object out of them
        interesting = (self.covariance[range(self.m), range(self.m)] > 0).nonzero()[0]

        # create a PDF structure out of this covariance matrix
        self.pdf = PDF.MVG(parameters=np.array([p['parameter'][0] for p in np.array(self.parameters)[interesting]]),
            covariance=self.covariance[interesting,:][:,interesting])
        self.pdf.printParameters()

            #self.save()
"""

class MCMC(Fit):
    def __init__(self, model, **kwargs):
        Fit.__init__(self, model)
        self.directory = self.model.directory + 'mcmc/'

    def fit(self, nburnin=500, ninference=500, nwalkers=100,
        broad=True, ldpriors=True,
        plot=True, interactive=False, remake=False, **kwargs):
        '''Use MCMC (with the emcee) to sample from the parameter probability distribution.'''

        self.speak('running an MCMC fit')
        try:
            assert(remake==False)
            self.load()
        except:


            self.model.instrument.rescaling.float(value=1.0,limits=[0.5, 2.0])

            self.model.defineParameterList()
            # populate an array with the parameters that are floating
            self.findFloating()

            # apply the limb darkening priors, if required
            if ldpriors:
                self.model.applyLDpriors()

            # pull out the parameters into an array for mpfit
            p0, parinfo = self.model.toArray()
            nparameters = len(self.floating)

            # setup the initial walker positions
            self.speak('initializing {nwalkers} for each of the {nparameters}'.format(**locals()))
            initialwalkers = np.zeros((nwalkers, nparameters))

            # loop over the parameters
            for i in range(nparameters):
                parameter = self.model.parameters[i]
                initialwalkers[:,i] = np.random.uniform(parameter.limits[0], parameter.limits[1], nwalkers)
                self.speak('  {parameter.name} picked from uniform distribution spanning {parameter.limits}'.format(**locals()))

            # set up the emcee sampler
            self.sampler = emcee.EnsembleSampler(nwalkers, nparameters, self.model.lnprob)
            self.names = [p.name for p in self.model.parameters]

            # add names to the sampler (for plotting progress, if desired)
            self.sampler.addLabels(self.names)

            # run a burn in step, and then reset
            burnt = False
            count = 0
            pos = initialwalkers
            while burnt == False:
                self.speak("running {0} burn-in steps, with {1} walkers.".format(nburnin, nwalkers))
                pos, prob, state = self.sampler.run_mcmc_with_progress(pos, nburnin)

                if plot:
                    self.sampler.HistoryPlot([count, count + nburnin])

                samples = {}
                for i in range(nparameters):
                    samples[self.floating[i]] = self.sampler.flatchain[:,i]

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
            ok = self.sampler.flatlnprobability > (np.max(self.sampler.flatlnprobability) - 100)
            samples = {}
            for i in range(nparameters):
            	samples[self.floating[i]] = self.sampler.flatchain[ok,i]

            # set the parameter to their MAP values
            best = self.sampler.flatchain[np.argmax(self.sampler.flatlnprobability)]
            self.model.fromArray(best)

            self.notes = {}
            self.notes['chisq'] = -2*self.model.lnprob(best)
            self.notes['dof'] = len(self.model.TLC.bjd) - len(self.floating)
            self.notes['reduced_chisq'] = self.notes['chisq']/self.notes['dof']
            self.notes['floating'] = self.floating
            self.pdf = PDF.Sampled(samples=samples)
            self.pdf.printParameters()
            self.notes['rescaling'] = np.array(self.pdf.values)[np.array(self.pdf.names) == 'rescaling'][0]
            self.save()

            if plot:
                self.model.TLC.DiagnosticsPlots(directory=self.directory)
"""
