import numpy as np
import matplotlib.pyplot as plt
from Planet import Planet
from Star import Star
from TM import TM
from PDF import PDF
import matplotlib.gridspec
import astropy.io
import pickle
import zachopy.color
ppm = 1e6
from zachopy.Talker import Talker

class TLC(Talker):
	'''Transit Light Curve class, to store both light curve and auxiliary variables.'''
	def __init__(self, bjd=None, flux=None, uncertainty=None, mearth=None, directory=None, left=None, right=None, name=None, remake=False, threshold=None, noiseassumedforplotting=None, night=None, **kwargs):

		Talker.__init__(self)
		self.flags = dict(outlier=1, saturation=2, custom=4)
		# define wavelengths
		self.left = left
		self.right = right
		self.directory = directory

		try:
			assert(remake == False)
			self.load(directory)
		except:
			# if possible, initialize from arrays
			if bjd is not None and flux is not None:
				self.fromArrays(bjd, flux, uncertainty, **kwargs)
			# try to load from MEarth
			if mearth is not None:
				self.fromMEarth(mearth, night=night, threshold=threshold)
				self.left = 6000.0
				self.right = 7000.0

		try:
			self.bad
		except:
			self.bad = np.isfinite(self.flux) == False

		if self.left is None or self.right is None:
			self.wavelength = None
		else:
			self.wavelength = (self.left + self.right)/2.0


		self.name = name

		if noiseassumedforplotting is None:
			self.noiseassumedforplotting = np.mean(self.uncertainty)
		else:
			self.noiseassumedforplotting = noiseassumedforplotting


		self.setupColors()

	def setupColors(self):
		# set up the appropriate colors to use
		try:
			self.colors
		except:
			self.colors = {}
			try:
				self.colors['points'] = zachopy.color.nm2rgb([self.left/10, self.right/10], 0.25)
				self.colors['lines'] = zachopy.color.nm2rgb([self.left/10, self.right/10], intensity=3.0)
			except:
				self.colors['points'] = 'black'
				self.colors['lines'] = 'gray'


	def fromMEarth(self, filename, night=None, threshold=None):
		'''Populate a TLC from a MEarth file.'''
		hdus = astropy.io.fits.open(filename)
		data = hdus[1].data
		header = hdus[1].header
		self.data, self.header = data, header
		#print data.columns
		i_target = (data['class'] == 9).nonzero()[0]
		assert(len(i_target == 1))
		i_target = i_target[0]
		bjd = data['bjd'][i_target]
		flux = 10**(-0.4*(data['flux'][i_target]))
		uncertainty = (1-10**(-0.4*(data['fluxerr'][i_target])))*flux
		if night is not None:
			ok = np.abs(bjd - night) < 0.5
		else:
			ok = bjd > 0
		if threshold is not None:
			ok *= uncertainty/flux < (np.median(uncertainty/flux)*3)
		self.fromArrays(bjd[ok], flux[ok]/np.median(flux[ok]), uncertainty[ok]/np.median(flux[ok]))

		# populate the external variables, both as a list and as individual entries
		self.externalvariables = {}


		for key in ['tv', 'texp', 'off', 'rms', 'extc', 'see', 'ell', 'sky', 'nois', 'fang', 'iang'  ]:
			self.externalvariables[key] = np.zeros(len(bjd))

			for i in range(len(bjd)):
				self.externalvariables[key][i] = header['{0}{1:.0f}'.format(key,i+1)]

			self.externalvariables[key] = self.externalvariables[key][ok]

		for key in ['xlc','ylc','airmass','sky']:
			value=data[key][i_target][ok]
			if len(value) == len(self.bjd):
				self.externalvariables[key] = value

	def fromArrays(self, bjd, flux, uncertainty=None, **kwargs):
		'''Populate a TLC from input arrays (used by transmission.py'''
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
				if key != 'bjd' and key != 'flux' and key !='uncertainty' and key !='left' and key !='right' and key !='wavelength':
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
				#self.speak( "      " + k)
				#if k == 'externalvariables':
					#self.speak( "      " + k + ', including:')
					#for evkey in tosave[k].keys():
						#self.speak( "          " + evkey)

		np.save(directory + 'TLC.npy', tosave)
		self.speak("saving light curve to {directory}, including all its external variables".format(directory=directory))

	def load(self, directory):
		self.directory = directory
		filename = directory + str(self.__class__).split('.')[-1].split("'")[0] + '.npy'
		self.speak('trying to load TLC from {0}'.format(filename))
		loaded = np.load(filename)[()]
		for key in loaded.keys():
			#self.speak( '     ' + key)
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
			assert(False)
		except:



			# set up the light curve plots
			if everything:
				# create two columns, to populate with lightcurves on left and diagnostics on right
				gs_overarching = matplotlib.gridspec.GridSpec(1, 3, width_ratios=[1, 0.3, 1.0], wspace=0.25, left=0.09, right=0.98)

				# create plots for light curves
				gs_lightcurve = matplotlib.gridspec.GridSpecFromSubplotSpec(5, 2, hspace=0.05, wspace=0, width_ratios = [1,.2], height_ratios=[1,0.5,0.1,1,0.5], subplot_spec = gs_overarching[0])

				# noise properties
				gs_diagnostics = matplotlib.gridspec.GridSpecFromSubplotSpec(3,1, height_ratios = [3,1,1], hspace=0.3, subplot_spec = gs_overarching[1])

				# potentially correlated variables
				self.to_correlate = []
				for k in self.TM.floating:
					if k != 'C':
						if k in self.TM.instrument.__dict__.keys() and k != 'rescaling':
							self.to_correlate.append(k.split('_tothe')[0])
				ncor = np.int(np.ceil(np.sqrt(len(self.to_correlate))))
				gs_external = matplotlib.gridspec.GridSpecFromSubplotSpec(len(self.to_correlate), 2, subplot_spec = gs_overarching[2], width_ratios=[1,5], hspace=0.1, wspace=0)

			else:
				gs_lightcurve = matplotlib.gridspec.GridSpec(5, 1, hspace=0.05, wspace=0, height_ratios=[1,0.5,0.1,1,0.5])

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
				self.ax_text = plt.subplot(gs_diagnostics[0])

				todisplay = self.TM.lastfit.pdf.listParameters()

				todisplay += '\n'
				todisplay += 'chisq = {chisq:.1f}/{dof:.0f}\nrescaling = {rescaling:.2f}'.format(**self.TM.lastfit.notes)

				self.ax_text.text(-0.1, 1.1, todisplay, va='top', fontsize=11)
				self.ax_text.axis('off')

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
			gs_lightcurve = matplotlib.gridspec.GridSpec(2, 2, width_ratios=[2,1], wspace=0.0, hspace=0.25)

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
			dict[evkey] = np.interp(new_bjd, self.TLC.bjd, self.TLC.externalvariables[evkey])

		# create the fake TLC
		return TLC(left=self.left, right=self.right, **dict)

	def LightcurvePlots(self):
		'''A quick tool to plot what the light curve (and external variables) looks like.'''

		# set up the phase/unphased lightcurve plots
		self.setupLightcurvePlots()

		# create smoothed TLC structures, so the modeling will work
		self.TM.smooth_phased_tlc = self.fake( np.linspace(-self.TM.planet.period.value/2.0 + self.TM.planet.t0.value + 0.01, self.TM.planet.period.value/2.0 + self.TM.planet.t0.value-0.01, 100000))
		self.TM.smooth_unphased_tlc = self.fake(np.linspace(np.min(self.TLC.bjd), np.max(self.TLC.bjd), 100000))



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

	def residualsexceptfor(self, exception=None):

		if exception is not None:
			original = {}
			for power in [1,2]:
				key = exception + '_tothe{0}'.format(power)
				original[key] = self.TM.instrument.__dict__[key].value + 0.0
			 	self.TM.instrument.__dict__[key].value = 0.0
			ref = self.flux/self.TM.model() - 1
			for key in original.keys():
				self.TM.instrument.__dict__[key].value = original[key]
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

	def DiagnosticsPlots(self, noiseassumedforplotting=0.001, directory=None):
		'''A quick tool to plot what the light curve (and external variables) looks like.'''

		self.speak('plotting light curves with diagnostics')
		self.setupDiagnostics()
		ppm=1e6
		# create smoothed TLC structures, so the modeling will work
		self.TM.smooth_phased_tlc = self.fake( np.linspace(-self.TM.planet.period.value/2.0 + self.TM.planet.t0.value + 0.01, self.TM.planet.period.value/2.0 + self.TM.planet.t0.value-0.01, 100000))
		self.TM.smooth_unphased_tlc = self.fake(np.linspace(np.min(self.TLC.bjd), np.max(self.TLC.bjd), 100000))

		goodkw = {'color':self.colors['points'], 'alpha':0.5, 'linewidth':0, 'marker':'o', 'markeredgecolor':self.colors['points'], 'markersize':6}
		badkw = {'color':self.colors['points'], 'alpha':0.25, 'linewidth':0, 'marker':'x', 'markeredgecolor':self.colors['points'], 'markersize':6}
		time = self.TM.planet.timefrommidtransit(self.bjd)

		notok = self.bad
		for good in [False,True]:
			if good:
				ok = (self.bad == 0).nonzero()
				kw = goodkw
			else:
				ok = (self.bad).nonzero()
				kw = badkw
			if np.sum(ok) == 0:
				continue
			self.points_raw = self.ax_raw.plot(time[ok], self.flux[ok], **kw)[0]
			self.points_corrected = self.ax_corrected.plot(time[ok], self.flux[ok]/self.TM.instrument_model()[ok], **kw)[0]
			self.points_residuals = self.ax_residuals.plot(time[ok], ppm*self.residuals()[ok], **kw)[0]
			self.points_instrument = self.ax_instrument.plot(time[ok], ppm*self.instrumental()[ok], **kw)[0]
			self.points_correlations = {}
			kw['markersize'] = 2
			kw['alpha'] *= 0.5
			for k in self.ax_correlations.keys():
				#self.ax_correlations[k].plot(ppm*self.instrumental()[ok], self.externalvariables[k][ok], **kw)[0]
				#self.ax_timeseries[k].plot( time[ok], self.externalvariables[k][ok], **kw)

				# plot the external variables
				x = self.externalvariables[k]
				self.ax_correlations[k].plot(ppm*(self.residualsexceptfor(k)[ok]), x[ok], **kw)[0]
				self.ax_timeseries[k].plot( time[ok], self.externalvariables[k][ok], **kw)

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
		assert(np.std(ppm*justinstrument)>1)

		nsigma = 5
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

		filename = directory + 'lightcurveDiagnostics.pdf'
		self.speak('saving light curve diagnostic plot to {0}'.format(filename))
		plt.savefig(filename)

	def linkModel(self, tm):
		self.TM = tm
		self.TM.TLC = self
		self.TLC = self
		self.TM.TM = self.TM

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
