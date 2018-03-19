import numpy as np, matplotlib.pyplot as plt
import eb
from .Planet import Planet
from .Star import Star
from .Instrument import Instrument
from .plots.SingleTransitPlot import SingleTransitPlot
import craftroom.color, craftroom.oned
import matplotlib.gridspec
import matplotlib.patches
import copy
from craftroom.Talker import Talker
from .Fits import LM, MCMC
ppm = 1e6

class TM(Talker):
	'''Transit Model object handles generation of model transit light curves.
		(relies heavily on Jonathan Irwin's "eb" code, which is an updated
		implementation of the classic EBOP and JKTEBOP, available at:
			https://github.com/mdwarfgeek'''

	def __init__(self, planet=None, star=None, instrument=None, directory=None,depthassumedforplotting=None, **kwargs):
		'''Initialize the parameters of the model.'''

		# setup the talker
		Talker.__init__(self)


		# keep track of a depth for plotting, if necessary
		self.depthassumedforplotting=depthassumedforplotting


		# load the model, if possible, and if a directory was given
		self.directory = directory
		if directory is not None and planet is None:
			self.load(directory)
		else:
			# define the subsets of the parameters, maybe by falling to defaults
			if planet is None:
				planet = Planet()
			if star is None:
				star = Star()
			if instrument is None:
				instrument = Instrument()
			self.planet = planet
			self.star = star
			self.instrument = instrument


	def load(self, directory):
		'''Load parameters from .directory.'''
		#self.directory = directory
		self.speak('loading TM from .{0}'.format(directory))
		self.planet = Planet(directory=directory)
		self.star = Star(directory=directory)
		self.instrument = Instrument(directory=directory)

	def save(self, directory):
		'''Save parameters to directory.'''
		#self.directory = directory
		for x in (self.planet, self.star, self.instrument):
			x.save(directory)

	def linkLightCurve(self, transitlightcurve):
		'''Attach a model to a light curve, defining all the TM and TLC attributes.'''
		self.TLC = transitlightcurve
		self.TM = self
		self.TLC.TM = self
		self.TLC.TLC = self.TLC

	def linkRV(self, rvcurve):
		'''Attach a model to a radial velocity curve, defining all the TM and RVC attributes.'''
		self.RVC = rvcurve
		self.TM = self
		self.RVC.TM  = self
		self.RVC.RVC = self.RVC


	##@profile
	def instrument_model(self, tlc=None):
		'''Model of the instrument.'''

		# by default, use the real TLC
		if tlc is None:
			tlc = self.TLC

		# use the instrument to spit out a corrective model
		return self.instrument.model(tlc)

	##@profile
	def model(self, tlc=None):
		'''Model including both instrument and planetary transit.'''

		# create the complete model, for either real or fake light curve
		return self.planet_model(tlc=tlc)*self.instrument_model(tlc=tlc)

	@property
	def parameters(self):
		'''Return a list of the parameters that are variable.'''
		try:
			assert(len(self._parameters) == len(self.floating))
			return self._parameters
		except:
			self.defineParameterList()
			return self._parameters

	@parameters.setter
	def parameters(self, **kwargs):
		pass

	def defineParameterList(self):
		'''Set up the parameter list, by pulling the variable parameters out of the subsets.'''

		# define a list containing the keys all the parameters that float
		self.floating = []
		list = []
		for x in (self.planet, self.star, self.instrument):
			d = x.__dict__
			for key in d.keys():
				try:
					if d[key].fixed == False:
						self.floating.append(key)
						list.append(d[key])
				except:
					pass
		self._parameters = np.array(list)

	def	fromArray(self, array):
		'''Use an input array to assign the internal parameter attributes.'''
		count = 0
		for parameter in self.parameters:
			parameter.value = array[count]
			count += 1
		#print count
		#print array
		assert(count > 0)

	def	fromlmfitParams(self, params):
		'''Use an input lmfit parameter list object to assign the internal parameter attributes.'''
		count = 0
		for parameter in self.parameters:
			parameter.value = params[parameter.name].value
			count += 1
		#print count
		#print array
		assert(count > 0)

	def toArray(self):
		'''Define an parameter array, by pulling them out of the internal parameter attributes.'''
		list = []
		parinfolist = []
		for parameter in self.parameters:
			list.append(parameter.value)
			parinfolist.append(parameter.parinfo)
		return list, parinfolist

	def plotPhased(self, smooth=None, **kw):
		'''Plot the light curve model, phased with the planetary period.'''
		t_phased = self.planet.timefrommidtransit(self.smooth_phased_tlc.bjd)
		assert(len(t_phased) == len(self.model(self.smooth_phased_tlc)))

		toplot = self.planet_model(self.smooth_phased_tlc)
		if smooth is not None:
			cadence = np.mean(self.smooth_phased_tlc.bjd[1:] - self.smooth_phased_tlc.bjd[:-1])
			n = np.round(smooth/cadence).astype(np.int)
			kernel = np.ones(n)/n
			toplot = np.convolve(toplot, kernel, 'same')
		else:
			toplot = self.planet_model(self.smooth_phased_tlc)
		#plt.plot(t_phased, toplot, **kw)
		# '''KLUDGED COMMENTED OUT!'''
		#try:
		#	for phased in [self.line_phased[0], self.line_phased_zoom[0]]:
		#		phased.set_data(t_phased, self.model(self.smooth_phased_tlc))
		#except AttributeError:
		self.line_phased = self.TLC.ax_phased.plot(t_phased,self.model(self.smooth_phased_tlc), **self.kw)
		self.line_phased_zoom = self.TLC.ax_phased_zoom.plot(t_phased, self.model(self.smooth_phased_tlc), **self.kw)

	def plotUnphased(self):
		'''Plot the light curve model, linear in time.'''
		t_unphased = self.smooth_unphased_tlc.bjd - self.planet.t0.value
		assert(len(t_unphased) == len(self.model(self.smooth_unphased_tlc)))
		#try:
		#	for unphased in [self.line_unphased[0], self.line_unphased_zoom[0]]:
		#		unphased.set_data(t_unphased, self.model(self.smooth_unphased_tlc))
		#except AttributeError:
		self.line_unphased = self.TLC.ax_unphased.plot(t_unphased, self.model(self.smooth_unphased_tlc), **self.kw)
		self.line_unphased_zoom = self.TLC.ax_unphased_zoom.plot(t_unphased, self.model(self.smooth_unphased_tlc), **self.kw)


	def plotDiagnostics(self):
		'''Plot the light curve model, linear in time.'''
		t_unphased = self.smooth_unphased_tlc.bjd - self.planet.t0.value
		assert(len(t_unphased) == len(self.model(self.smooth_unphased_tlc)))

		#try:
		#	self.line_raw.set_data(t_unphased, self.model(self.smooth_unphased_tlc))
		#	self.line_corrected.set_data(t_unphased, self.planet_model(self.smooth_unphased_tlc))
		#	self.line_residuals.set_data(t_unphased, ppm*np.zeros_like(t_unphased))
		#	self.line_instrument.set_data(t_unphased, ppm*self.instrument_model(self.smooth_unphased_tlc))
		#except:
		# NOT USING?!?
		self.line_raw = self.TLC.ax_raw.plot(t_unphased, self.model(self.smooth_unphased_tlc), **kw)[0]
		self.line_corrected = self.TLC.ax_corrected.plot(t_unphased, self.planet_model(self.smooth_unphased_tlc), **kw)[0]
		self.line_residuals = self.TLC.ax_residuals.plot(t_unphased, ppm*np.zeros_like(t_unphased), **kw)[0]
		self.line_instrument = self.TLC.ax_instrument.plot(t_unphased, ppm*self.instrument_model(self.smooth_unphased_tlc), **kw)[0]

	def plot(self):
		'''Plot the model lines over the existing light curve structures.'''
		self.kw = {'color':self.TLC.colors['lines'], 'linewidth':3, 'alpha':1.0}


		# create smoothed TLC structures, so the modeling will work
		self.smooth_phased_tlc = self.TLC.fake(np.linspace(-self.planet.period.value/2.0 + self.planet.t0.value + 0.01, self.planet.period.value/2.0 + self.planet.t0.value-0.01, 10000))
		self.smooth_unphased_tlc = self.TLC.fake(np.linspace(np.min(self.TLC.bjd), np.max(self.TLC.bjd), 1000))

		self.plotPhased()
		self.plotUnphased()

		if self.depthassumedforplotting is None:
			self.depthassumedforplotting = self.planet.rp_over_rs.value**2

		# need to work on displaying the parameters...
		try:
			xtext = self.planet.duration/2.0*1.1
			posttransit = self.planet.timefrommidtransit(bjd) > self.planet.duration/2.0
			ytext = np.mean(self.model()[posttransit]) - self.planet.rp_over_rs**2/2.0
			instrument_string = "Instrumental Parameters"
			self.TLC.ax_raw.text(xtext, ytext, instrument_string)
		except:
			pass


	##@profile
	def lmfitdeviates(self, params, fjac=None, plotting=False):
		'''Return the normalized deviates (an input for lmfit).'''

		# populate the parameter attributes, using the input array
		values = params
		self.fromlmfitParams(params)

		# if necessary, plot the light curve along with this step of the deviates calculation
		ok = (self.TLC.bad == 0).nonzero()
		devs = (self.TLC.flux[ok] - self.TM.model()[ok])/self.TLC.uncertainty[ok]

		if plotting:
			self.plotter = SingleTransitPlot(tlc=self.TLC)
			print('chisq = {}'.format(np.sum(devs**2)))
			plt.show()



		# add limb darkening priors (these are just independent right now, but that could be better)
		try:
			prioru1 = (self.star.u1.value - self.u1prior_value)/self.u1prior_uncertainty
			prioru2 = (self.star.u2.value - self.u2prior_value)/self.u2prior_uncertainty
			devs = np.append(devs, prioru1)
			devs = np.append(devs, prioru2)
			#print '==============================='
			#print 'u1: ({value} - {center})/{uncertainty}'.format(value = self.star.u1.value, center=self.u1prior_value, uncertainty =self.u1prior_uncertainty)
			#print 'u2: ({value} - {center})/{uncertainty}'.format(value = self.star.u2.value, center=self.u2prior_value, uncertainty =self.u2prior_uncertainty)

		except:
			pass

		# lmfit just wants the deviates
		return devs

	##@profile
	def deviates(self, p, fjac=None, plotting=False):
		'''Return the normalized deviates (an input for mpfit).'''

		# populate the parameter attributes, using the input array
		self.fromArray(p)

		# if necessary, plot the light curve along with this step of the deviates calculation
		status =0
		if plotting:
			self.TLC.LightcurvePlots()

		ok = (self.TLC.bad == 0).nonzero()

		devs = (self.TLC.flux[ok] - self.TM.model()[ok])/self.TLC.uncertainty[ok]

		# add limb darkening priors
		try:
			prioru1 = (self.star.u1.value - self.u1prior_value)/self.u1prior_uncertainty
			prioru2 = (self.star.u2.value - self.u2prior_value)/self.u2prior_uncertainty
			devs = np.append(devs, prioru1)
			devs = np.append(devs, prioru2)
			#print '==============================='
			#print 'u1: ({value} - {center})/{uncertainty}'.format(value = self.star.u1.value, center=self.u1prior_value, uncertainty =self.u1prior_uncertainty)
			#print 'u2: ({value} - {center})/{uncertainty}'.format(value = self.star.u2.value, center=self.u2prior_value, uncertainty =self.u2prior_uncertainty)

		except:
			pass

		# mpfit wants a list with the first element containing a status code
		return [status,devs]

	def applyLDpriors(self):
		self.speak('using atmosphere model as prior on LD coefficients')
		self.u1prior_value = self.star.u1.value + 0.0
		self.u2prior_value = self.star.u2.value + 0.0
		self.u1prior_uncertainty = self.star.u1.uncertainty + 0.0
		self.u2prior_uncertainty = self.star.u2.uncertainty + 0.0

	##@profile
	def lnprob(self, p):
		"""Return the log posterior, calculated from .the TM.deviates function (which may have included some conjugate Gaussian priors.)"""



		chisq = np.sum(self.deviates(p)[-1]**2)/2.0
		N = np.sum(self.TLC.bad == 0)

		# sum the deviates into a chisq-like thing
		lnlikelihood = -N * np.log(self.instrument.rescaling.value) - chisq/self.instrument.rescaling.value**2
		if np.isfinite(lnlikelihood) == False:
			lnlikelihood = -1e9

		# initialize an empty constraint, which could freak out if there's something bad about this fit
		constraints = 0.0

		# loop over the parameters


		for parameter in self.parameters:

			# if a parameter is outside its allowed range, then make the constraint very strong!
			inside = (parameter.value < parameter.limits[1]) & (parameter.value > parameter.limits[0])
			try:
				assert(inside)
			except AssertionError:
				constraints -= 1e6

		# return the constrained likelihood
		return lnlikelihood + constraints

	def fastfit(self, **kwargs):
		self.fitter = LM(self, **kwargs)
		self.fitter.fit(**kwargs)

	def slowfit(self, **kwargs):
		self.fitter = MCMC(self, **kwargs)
		self.fitter.fit(**kwargs)

	def __repr__(self):
		return self.TLC.__repr__().replace('TLC', 'TM')




class TMKreidberg(TM):
	#@profile
	def planet_model(self, tlc=None, t=None):
		'''Model of the planetary transit.'''
		self.set_ebparams()

		# by default, will use the linked TLC, but could use a custom TLC
		# 	 (e.g. a high-resolution one, for plotting)
		if tlc is None:
			try:
				tlc = self.TLC
			except AttributeError:
				self.speak('no TLC is connected to the TM')

		# if called with a "t=" keyword set, then will use those custom times
		if t is None:
			t = tlc.bjd

		# MAKE A BATMAN/SPIDERMAN MODEL!

		# work in relative flux (rather than magnitudes) -- why did I do this?
		return 10**(-0.4*eb.model(self.ebparams, t, typ))





class TMIrwin(TM):

	def __init__(self, *args, **kwargs):
		TM.__init__(self, *args, **kwargs)
		
		# create an empty array of parameters for eb
		self.ebparams = np.zeros(eb.NPAR, dtype=np.double)

	#@profile
	def set_ebparams(self):
		'''Set up the parameters required for eb. '''
		# These are the basic parameters of the model.
		self.ebparams[eb.PAR_J]      =  self.planet.surface_brightness_ratio  # J surface brightness ratio
		self.ebparams[eb.PAR_RASUM]  =  self.planet.rsum_over_a  # (R_1+R_2)/a
		self.ebparams[eb.PAR_RR]     =  self.planet.rp_over_rs  # R_2/R_1
		self.ebparams[eb.PAR_COSI]   =  self.planet.cosi  # cos i

		# Mass ratio is used only for computing ellipsoidal variation and
		# light travel time.  Set to zero to disable ellipsoidal.
		self.ebparams[eb.PAR_Q]      =  0#self.planet.q

		# Light travel time coefficient.
		#ktot = 55.602793  # K_1+K_2 in km/s
		#cltt = 1000*ktot / eb.LIGHT

		# Set to zero if you don't need light travel correction (it's fairly slow
		# and can often be neglected).
		self.ebparams[eb.PAR_CLTT]   =  0#cltt#*(self.planet.q.value == 0.0)      # ktot / c

		# Radiative properties of star 1.
		self.ebparams[eb.PAR_LDLIN1] = self.star.u1.value   # u1 star 1
		self.ebparams[eb.PAR_LDNON1] = self.star.u2.value  # u2 star 1
		self.ebparams[eb.PAR_GD1]    = self.star.gd.value     # gravity darkening, std. value
		self.ebparams[eb.PAR_REFL1]  = self.star.albedo.value      # albedo, std. value

		# Spot model.  Assumes spots on star 1 and not eclipsed.
		self.ebparams[eb.PAR_ROT1]   =  1.0# 0.636539  # rotation parameter (1 = sync.)
		self.ebparams[eb.PAR_FSPOT1] =  0.0       # fraction of spots eclipsed
		self.ebparams[eb.PAR_OOE1O]  =  0.0       # base spottedness out of eclipse
		self.ebparams[eb.PAR_OOE11A] =  0.0#0.006928  # *sin
		self.ebparams[eb.PAR_OOE11B] =  0.0 # *cos

		# PAR_OOE12* are sin(2*rot*omega) on star 1,
		# PAR_OOE2* are for spots on star 2.

		# Assume star 2 is the same as star 1 but without spots.
		self.ebparams[eb.PAR_LDLIN2] = self.ebparams[eb.PAR_LDLIN1]
		self.ebparams[eb.PAR_LDNON2] = self.ebparams[eb.PAR_LDNON1]
		self.ebparams[eb.PAR_GD2]    = self.ebparams[eb.PAR_GD1]
		self.ebparams[eb.PAR_REFL2]  = self.ebparams[eb.PAR_REFL1]

		# Orbital parameters.
		self.ebparams[eb.PAR_ECOSW]  =  self.planet.ecosw.value  # ecosw
		self.ebparams[eb.PAR_ESINW]  = self.planet.esinw.value  # esinw
		self.ebparams[eb.PAR_P]      = self.planet.period.value  # period
		self.ebparams[eb.PAR_T0]     = self.planet.t0.value + self.planet.dt.value # T0 (epoch of primary eclipse), with an offset of dt applied
		# OTHER NOTES:
		#
		# To do standard transit models (a'la Mandel & Agol),
		# set J=0, q=0, cltt=0, albedo=0.
		# This makes the secondary dark, and disables ellipsoidal and reflection.
		#
		# The strange parameterization of radial velocity is to retain the
		# flexibility to be able to model just light curves, SB1s, or SB2s.
		#
		# For improved precision, it's best to subtract most of the "DC offset"
		# from .T0 and the time array (e.g. take off the nominal value of T0 or
		# the midtime of the data array) and add it back on at the end when
		# printing self.ebparams[eb.PAR_T0] and vder[eb.PAR_TSEC].  Likewise the period
		# can cause scaling problems in minimization routines (because it has
		# to be so much more precise than the other parameters), and may need
		# similar treatment.

	def stellar_rv(self, rvc=None, t=None):
		self.set_ebparams()

		# by default, will use the linked RVC, but could use a custom one (e.g. high-resolution for plotting), or just times
		if rvc is None:
			rvc = self.RVC

		if t is None:
			t = rvc.bjd

		# make sure the types are okay for Jonathan's inputs
		typ = np.empty_like(t, dtype=np.uint8)
		typ.fill(eb.OBS_VRAD1)

		rv = self.planet.semiamplitude.value*eb.model(self.ebparams, t, typ) + self.star.gamma.value
		assert(np.isfinite(rv).all())
		return rv

	#@profile
	def planet_model(self, tlc=None, t=None):
		'''Model of the planetary transit.'''
		self.set_ebparams()

		# by default, will use the linked TLC, but could use a custom TLC
		# 	 (e.g. a high-resolution one, for plotting)
		if tlc is None:
			try:
				tlc = self.TLC
			except AttributeError:
				self.speak('no TLC is connected to the TM')

		# if called with a "t=" keyword set, then will use those custom times
		if t is None:
			t = tlc.bjd

		# make sure the types are okay for Jonathan's inputs
		typ = np.empty_like(t, dtype=np.uint8)
		typ.fill(eb.OBS_MAG)

		# work in relative flux (rather than magnitudes) -- why did I do this?
		return 10**(-0.4*eb.model(self.ebparams, t, typ))
