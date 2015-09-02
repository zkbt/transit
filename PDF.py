'''Keep track of a multidimensional PDF.'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec
from transit.Parameter import Parameter
from zachopy.Talker import Talker
import copy
import scipy.special
import triangle


# which attributes can be saved?
saveable = ['names', 'values', 'covariance', 'samples']

def load(filename):
	return Sampled(filename=filename)

class PDF(Talker):
	"""an instance of Probability Density Function can [visualize, print, marginalize] a distribution"""

	def __init__(self, filename=None, **kwargs):
		"""general initialization for PDF instances"""

		# initialize the talker
		Talker.__init__(self)

		# by default, PDF starts out uninitialized
		self._initialized = False

		# if possible, load the PDF from a file
		if filename is not None:
			self.load(filename)

	def listParameters(self):
		s = []
		for i in range(self.n):
			p = self.parameters[i]
			thisline = '{0} = {1}\n'.format(p.name, p.string())
			#if len(thisline) > 35:
			#	thisline = thisline.replace('=', '=\n')
			s.append(thisline)
		return s

	def printParameters(self):
		for p in self.listParameters():
			self.speak(p)

	def load(self, filename):
		"""load the PDF from a file"""

		# load the saved dictionary
		loaded = np.load(filename)[()]
		for k in saveable:
			self.__dict__[k] = loaded[k]

		# recreate the parameter objects
		self.parameters = []
		for i in range(len(self.names)):
			self.parameters.append(Parameter(self.names[i], self.values[i]))

		#self.populateUncertainties()

		# record the fact that this has been initialized
		self._initialized = True

		# provide an update
		self.speak('loaded PDF from {0}'.format(filename))

	def save(self, filename):
		"""save the PDF to a file"""

		# create a dictionary of the things to be saved
		tosave = {}
		for k in saveable:
			tosave[k] = self.__dict__[k]

		# save the dictionary
		np.save(filename, tosave)

		# provide an update
		self.speak('saved PDF to {0}'.format(filename))

	@property
	def n(self):
		return len(self.parameters)

	def populateUncertainties(self):
		for i in range(self.n):
			self.parameters[i].uncertainty = np.sqrt(self.covariance[i,i])

	@property
	def sigmas(self):
		i = np.arange(self.n)
		return np.sqrt(self.covariance[i, i])

	@property
	def correlation(self):
		return self.covariance/self.sigmas.reshape(self.n, 1)/self.sigmas.reshape(1, self.n)

	def export(self, filename, keys=None):

		if keys == None:
			keys = self.names

		self.speak('saving dictionary to {0} containing samples for'.format(filename))
		toexport = {}
		for k in keys:
			toexport[k] = self.samples[k]
			self.speak('   {0}'.format(k))
		np.save(filename, toexport)
		self.speak('done!')

	def triangle(self, keys=None, truths=None, quantiles=[0.16, 0.5, 0.84], title=None, dpi=70, **kwargs):
		data = []

		for k in keys:
			try:
				ok *= np.isfinite(self.samples[k])
			except UnboundLocalError:
				ok = np.isfinite(self.samples[k])

		for k in keys:
			data.append(self.samples[k][ok])

		data = np.vstack(data).transpose()

		# make the plot
		figure = triangle.corner(data,
								labels=keys,
								truths=truths,
								quantiles=quantiles,
								title_args={"fontsize": 12},
								levels = 1.0 - np.exp(-0.5 * np.array([1.0,2.0]) ** 2),
								**kwargs)
		figure.set_dpi(dpi)
		plt.draw()
		if title is not None:
			figure.gca().annotate(title, xy=(0.5, 1.0), xycoords="figure fraction",
                      xytext=(0, -5), textcoords="offset points",
                      ha="center", va="top")

		return figure


	def plot(self, keys=None, plotcovariance=False, onesigmalabels=False, subsample=10000, nbins=100, dye=None):
		'''Make a matrix plot of the PDF.'''

		self.plotcovariance = plotcovariance
		self.onesigmalabels = onesigmalabels
		self.subsample = subsample
		self.dye=dye

		# decide which elements to plot
		if keys is None:
			keys = self.names


		# set up the grid of subplots
		n = len(keys)
		scale = 2.0
		self.figure = plt.figure('matrix', figsize=(self.n*scale, self.n*scale), dpi=700.0/self.n/scale)

		gs = matplotlib.gridspec.GridSpec(n,n, hspace=0.05, wspace=0.05, left=0.1, right=0.95, bottom=0.1, top=0.95)
		self.pdfax = {}

		# loop over columns of the grid (i)
		for i in range(n):
			self.pdfax[keys[i]] = {}
			# loop over rows in the columns (j)
			for j in range(n):
				try:
					self.pdfax[keys[i]][keys[j]].cla()
				except:
					pass
				# setup axis sharing, always share the x axes
				try:
					assert(i != j)
					assert(i != (j-1))
					sharex = self.pdfax[keys[i]][keys[j-1]]
				except:
					sharex = None
				# don't share y axis with the histograms
				try:
					assert(i != j)
					assert((i-1) != j)
					sharey = self.pdfax[keys[i-1]][keys[j]]
				except:
					sharey = None

				# make a plot if in the lower left of the grid
				if i <= j:

					# create (and store) this subplot window
					ax = plt.subplot(gs[j,i], sharex=sharex, sharey=sharey)
					self.pdfax[keys[i]][keys[j]] = ax

					# set the axis labels, along the outside of the plot
					if i > 0:
						plt.setp(ax.get_yticklabels(), visible=False)
					else:
						ax.set_ylabel(keys[j])
					if j < (n-1):
						plt.setp(ax.get_xticklabels(), visible=False)
					else:
						ax.set_xlabel(keys[i])
						plt.sca(ax)
						locs, labels = plt.xticks()
						plt.setp(labels, rotation=90)

					# populate the subplots with data
					if i != j:
						self.plotpair(keys[i], keys[j])
					if i == j:
						self.plothist(keys[i], nbins=nbins)

					if self.onesigmalabels:
						self.fix_ticks('x',keys[i])
						if i != j:
							self.fix_ticks('y',keys[j])

	def fix_ticks(self, which, key, nsigma=2):
		'''Replace the nasty overlapping ticks with angles ones at +/- 2 sigma.'''

		# basics
		i = self.names.index(key)
		ax = plt.gca()

		# position ticks at -nsigma, mean, +nsigma
		ticks = [self.parameters[i].value - self.parameters[i].uncertainty*nsigma, self.parameters[i].value, self.parameters[i].value + self.parameters[i].uncertainty*nsigma]

		# format the tick labels
		def postdecimal(str):
			'''a helper function to figure out how many digits after the decimal point a string contains.'''
			return len(str.split('.')[-1])

		# keep two significant digits on the uncertainties
		uformat = '{0:+.3g}'.format

		# keep the same number of (total) digits for the central value
		vformatstring = '{0}'.format(postdecimal(uformat(self.parameters[i].uncertainty)))
		vformat = ('{0:.'+vformatstring+'f}').format

		# format the text for ticks at -nsigma, mean, +nsigma
		ticklabels = [uformat(-self.parameters[i].uncertainty*nsigma), vformat(self.parameters[i].value), uformat(self.parameters[i].uncertainty*nsigma)]

		# apply the ticks to the correct axis
		if which == 'x':
			lines = ax.set_xticks(ticks)
			labels = ax.set_xticklabels(ticklabels, rotation=45)
			ax.set_xlim(self.parameters[i].value - self.parameters[i].uncertainty*nsigma*2, self.parameters[i].value + self.parameters[i].uncertainty*nsigma*2)
		if which == 'y':
			lines = ax.set_yticks(ticks)
			labels = ax.set_yticklabels(ticklabels, rotation=45)
			ax.set_ylim(self.parameters[i].value - self.parameters[i].uncertainty*nsigma*2, self.parameters[i].value + self.parameters[i].uncertainty*nsigma*2)

		# nudge the tick label sizes
		labels[0].set_size('small')
		labels[1].set_weight('extra bold')
		labels[-1].set_size('small')

	def plothist(self, key, nbins=100):

		# get the plot to populate
		ax = self.pdfax[key][key]

		# plot the histogram of the samples
		ax.hist(self.samples[key], bins=nbins, linewidth=0, alpha=0.5, color=self.color, normed=True)

		# plot the Gaussian approximation
		if self.plotcovariance:
			x = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
			i = self.names.index(key)
			mean = self.parameters[i].value
			sigma = self.parameters[i].uncertainty
			plt.plot(x, 1.0/np.sqrt(2*np.pi)/sigma*np.exp(-0.5*(x-mean)**2/sigma**2), color=self.color, alpha=0.6, linewidth=3)

	def plotpair(self, keyx, keyy):

		stride = np.maximum(len(self.samples[keyx])/self.subsample, 1)
		theta = np.linspace(0, 2*np.pi, 100)

		i = self.names.index(keyx)
		j = self.names.index(keyy)


		ax = self.pdfax[keyx][keyy]
		ax.plot(self.samples[keyx][::stride], self.samples[keyy][::stride],  marker='o', alpha=0.1, linewidth=0, color=self.color, markeredgecolor='none')
		nsigma = 2

		if self.plotcovariance:
			bigcovar = self.covariance
			covar = bigcovar[[i,j],:][:,[i,j]]
			cholesky = np.linalg.cholesky(covar)
			theta = np.linspace(0, 2*np.pi, 100)
			x_ellipse, y_ellipse = cholesky.dot(np.array([np.sin(theta), np.cos(theta)]))
			for nsigma in [1,2]:
				kw = {"color":self.color, 'linewidth':5-nsigma*1, 'alpha':0.8 - 0.2*nsigma}
				plt.plot(nsigma*x_ellipse + self.parameters[i].value, nsigma*y_ellipse + self.parameters[j].value, **kw)


class Sampled(PDF):

	def __init__(self, samples=None, summarize=True, **kwargs):

		"""initialize a PDF from a dictionary of samples"""
		self.color = 'Black'
		PDF.__init__(self, **kwargs)
		if self._initialized:
			return

		# by default, the covariance matrix is empty (can calculate if need be)
		self.covariance = None

		# store the input samples
		self.samples = samples

		# pull out the parameter names
		self.names = self.samples.keys()
		if 'lnprob' in self.names:
			self.names.append(self.names.pop(self.names.index('lnprob')))

		# define the parameters, using the samples
		self.parameters = []
		for key in self.names:
			self.parameters.append(Parameter(key,None))

		# calculate the values
		self.recenter()

		if summarize:
			# calculate the covariance matrix
			self.calculateCovariance()

			# populate the parameter uncertainties
			self.populateUncertainties()


	def recenter(self, method=np.mean):
		"""recenter the parameter values from the samples, by default using the median"""
		for p in self.parameters:
			key = p.name
			p.value = method(self.samples[key])
		self.values = [parameter.value for parameter in self.parameters]

	def calculateCovariance(self):
		# figure out covariance from samples
		self.covariance = np.zeros((self.n, self.n))
		for i in range(self.n):
			for j in range(self.n):
				dx = self.samples[self.parameters[i].name][:] - self.parameters[i].value
				dy = self.samples[self.parameters[j].name][:] - self.parameters[j].value
				self.covariance[i,j] = np.mean(dx*dy)

	def calculateUncertainties(self, style='percentiles'):

		for i in range(self.n):

			if style == 'percentiles':
				edge = (1 - scipy.special.erf(1.0/np.sqrt(2)))/2.0
				limits = [100*edge, 100*(1-edge)]
				span = np.percentile(self.samples[self.parameters[i].name], limits)
				value = np.mean(span)
				uncertainty = (span[1] - span[0])/2.0
			elif style == 'std':
				value = np.mean(self.samples[self.parameters[i].name])
				uncertainty = np.std(self.samples[self.parameters[i].name])

			self.parameters[i].value = value
			self.parameters[i].uncertainty = uncertainty



class MVG(PDF):
	def __init__(self, names=None, parameters=None, covariance=None, **kwargs):
		"""initialize a PDF from a list of parameter objects and a covariance matrix"""
		self.color = 'SeaGreen'
		PDF.__init__(self, **kwargs)
		if self._initialized:
			return

		# by default, there are no samples (can make some if we need them)
		self.samples = None
		self.parameters = copy.deepcopy(parameters)
		if names is not None:
			assert(len(names) == len(parameters))
			for i in range(len(names)):
				self.parameters[i].name = names[i]



		self.names = [p.name for p in self.parameters]
		self.values = [p.value for p in self.parameters]
		self.covariance = covariance

		# populate the parameter uncertainties
		self.populateUncertainties()

	def simulateSamples(self, n=100):
		'''Use parameter values and covariances to generate samples
				(requires that these are both defined ahead of time.)'''
		parameters = self.parameters
		covariance = self.covariance
		means = [parameter.value for parameter in parameters]

		s = np.random.multivariate_normal(means, covariance, n)
		self.samples = {}
		for i in range(len(parameters)):
			self.samples[parameters[i].name] = s[:,i]
