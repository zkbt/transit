from zachopy.Talker import Talker
import matplotlib.pyplot as plt, numpy as np
import zachopy.borrowed.mpfit.mpfit as mpfit
import zachopy.oned
import pemcee as emcee
import transit.PDF as PDF
##@profile

class Fit(Talker):
    def __init__(self, model, **kwargs):
        Talker.__init__(self)
        self.model = model
        self.model.lastfit = self

    def findFloating(self):
        # determine which parameters are floating
        self.floating = []
        for x in (self.model.planet, self.model.star, self.model.instrument):
            d = x.__dict__
            for key in d.keys():
              try:
                if d[key].fixed == False:
                  self.floating.append(key)
                  self.speak('    '+key)
              except:
                pass


    def save(self):
        """Save this fit, so it be reloaded quickly next time."""
        zachopy.utils.mkdir(self.directory)
        self.speak('saving LM fit to {0}'.format(self.directory))
        self.speak('  the PDF')
        self.pdf.save(self.directory + 'pdf.npy')
        self.speak('  the fitting notes')
        np.save(self.directory + 'fitting_notes.npy', self.notes)
        self.speak('  the best-fit model')
        self.model.save(self.directory)


    def load(self):
        """Save this fit, so it be reloaded quickly next time."""
        self.speak('attempting to load from {0}'.format(self.directory))
        #self.speak('  the PDF')
        self.pdf = PDF.load(self.directory + 'pdf.npy')
        #self.speak('  the fitting notes')
        self.notes = np.load(self.directory + 'fitting_notes.npy')[()]
        #self.speak('  the best-fit model')
        self.model.load(self.directory)

class LM(Fit):
    def __init__(self, model, **kwargs):
        Fit.__init__(self, model)
        self.directory = self.model.directory + 'lm/'

    def fit(self, plot=False, quiet=True, ldpriors=True, identifyoutliers=True, remake=False, **kwargs):
        '''Use LM (mpfit) to find the maximum probability parameters, and a covariance matrix.'''

        self.speak('performing a fast LM fit')

        try:
            assert(remake==False)
            self.load()
        except:

            # populate an array with the parameters that are floating
            self.findFloating()
            # apply the limb darkening priors, if required
            if ldpriors:
                self.model.applyLDpriors()

            # pull out the parameters into an array for mpfit
            p0, parinfo = self.model.toArray()

            # perform the LM fit, to get best fit parameters and covariance matrix
            self.speak('running mpfit minimization')
            self.mpfitted = mpfit.mpfit(self.model.deviates, p0, parinfo=parinfo, quiet=quiet)

            # set the parameters to their fitted values
            for i in range(len(self.model.parameters)):
                self.model.parameters[i].value = self.mpfitted.params[i]

            # determine the uncertainties, including a rescaling term, by calculating the chisq of the good points
            ok = (self.model.TLC.bad == 0).nonzero()

            self.notes = {}
            self.notes['chisq'] = np.sum((self.model.TLC.residuals()[ok]/self.model.TLC.uncertainty[ok])**2)
            self.notes['dof'] =  self.mpfitted.dof
            self.notes['reduced_chisq'] = self.notes['chisq']/self.notes['dof']
            self.notes['rescaling'] = np.maximum(np.sqrt(self.notes['reduced_chisq']), 1)
            self.notes['floating'] = self.floating

            self.speak('acheived a chisq of {0:.2f}/{1} required a rescaling of {2:.2f}'.format(self.notes['chisq'] , self.notes['dof'], self.notes['rescaling']))

            # if we're tring to identify outliers, throw out the worst points and refit
            if identifyoutliers:

                # where are the residuals beyond 4 sigma?
                outlierthreshold = 4.0
                r = self.model.TLC.residuals()[ok]
                bad = (np.abs(r) > outlierthreshold*1.48*zachopy.oned.mad(r))

                # mark those points as bad
                self.model.TLC.bad[ok] = bad
                self.speak("identified {0} new points as bad; refitting without them".format(np.sum(bad)))

                # refit, after the outliers have been rejected
                self.fit(plot=plot, quiet=quiet, ldpriors=ldpriors, identifyoutliers=False,  **kwargs)

            # store the covariance matrix of the fit, and the 1D uncertainties on the parameters
            self.covariance = self.mpfitted.covar*self.notes['rescaling'] **2
            for i in range(len(self.model.parameters)):
                self.model.parameters[i].uncertainty = np.sqrt(self.covariance[i, i])

            # pull out the parameters that actually varied and create a PDF object out of them
            interesting = (self.covariance[range(len(self.model.parameters)), range(len(self.model.parameters))] > 0).nonzero()[0]

            # create a PDF structure out of this covariance matrix
            self.pdf = PDF.MVG(parameters=self.model.parameters[interesting],
                covariance=self.covariance[interesting,:][:,interesting])

            self.save()
            if plot:
                self.model.TLC.DiagnosticsPlots(directory=self.directory)


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
        except IOError:

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

            self.save()

            if plot:
                self.model.TLC.DiagnosticsPlots(directory=self.directory)
