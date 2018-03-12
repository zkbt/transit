from .imports import *
from .Data import Data
from .Planet import Planet
from .Star import Star
from .TM import TM
from .PDF import PDF

class RVC(Data):
	'''Radial Velocity Curve class,
		to store both radial velocities and auxiliary variables.'''

	def __init__(self, bjd=None, rv=None, uncertainty=None,
						inputfilename=None, directory=None,
						name=None,
						remake=False,
						color='slategray',
						isfake=False,
						telescope=None,
						**kwargs):

		# initialize the Talker object
		Data.__init__(self)
		self.speak('creating a RV curve')

		# initialize the base color
		self.color = color


		# define a dictionary of flags that can be used for bad data
		self.flags = dict(outlier=1, saturation=2, custom=4)


		# keep track of the telescope (and epoch)
		self.telescope=telescope
		self.epoch=None

		# is this a real light curve, or a fake one?
		#  (e.g., one at high resolution for plotting)
		self.isfake = isfake

		# specify the original filename
		self.inputfilename = inputfilename

		# specify a directory in which to store saved versions of this TLC,
		# 	as well as other outputs (plots, etc...)
		self.directory = directory
		if self.directory is None:
			try:
				self.directory = '/'.join(self.inputfilename.replace('data/', 'fits/').split('/')[:-1]) + '/'
				mkdir(self.directory)
			except AttributeError:
				pass
		self.speak('assigning it the directory {0}'.format(self.directory))

		# initialize the TLC by filling it with data
		self.initialize(bjd=bjd, rv=rv, uncertainty=uncertainty, remake=remake, **kwargs)

		# make sure an array of "bad" values is defined
		try:
			self.bad
		except AttributeError:
			try:
				self.bad = kwargs['bad']
			except KeyError:
				self.bad = np.isfinite(self.rv) == False

		# assign a name to this lightcurve
		if name is None:
			if self.telescope is None:
				name = '???'
			else:
				name = self.telescope
		self.name = name

		# assign the colors for this light curve
		self.setupColors(color=color, minimumuncertainty=None)

	@property
	def jitter(self):
		try:
			return self.TM.star.jitter.value
		except AttributeError:
			return 0.0

	@property
	def effective_uncertainty(self):
		return np.sqrt(self.uncertainty**2 + self.jitter**2)

	def initialize(self, bjd=None, rv=None, uncertainty=None,
						remake=False,
						**kwargs):
		'''If possible, try to load the RV curve from .its directory
			otherwise, create it from .raw input file.'''

		try:
			assert(remake == False)
			self.load(self.directory)
			self.speak('initialized RV from .pre-saved file in {0}'.format(self.directory))
		except:
			# if possible, initialize from .arrays; if not, load from .scratch
			if bjd is not None and rv is not None:
				self.fromArrays(bjd, rv, uncertainty, **kwargs)
				self.speak('initialized TLC from .{0}-element arrays'.format(self.n))
			else:
				self.fromFile(self.inputfilename)
				self.speak('initialized TLC from .{0}'.format(self.inputfilename))

			if self.isfake == False:
				self.save(self.directory)

	def residuals(self):
		return self.rv - self.TM.stellar_rv()

	def chisq(self):
		ok = self.bad == False
		return np.sum((self.residuals()/self.uncertainty)[ok]**2)


	def fromArrays(self, bjd, rv, uncertainty=None, **kwargs):
		'''Populate a TLC from .input arrays (used by transmission.py)'''

		# how many data points are in light curve?
		self.n = len(bjd)

		# define the times
		self.bjd = np.array(bjd)

		# define the rv array, and normalize it to its median
		self.rv = np.array(rv)

		# make sure the uncertainty is defined
		if uncertainty is None:
			uncertainty = np.ones_like(self.rv)
		else:
			self.uncertainty = np.array(uncertainty)


		# populate the external variables, both as a list and as individual entries
		self.cotrending = {}
		for key, value in kwargs.iteritems():
			if len(value) == len(self.bjd):
				if key != 'bjd' and key != 'rv' and key !='uncertainty' and key !='left' and key !='right' and key !='wavelength' and key!= 'bad':
					self.cotrending[key] = value
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
				#self.speak( "	  " + k)
				#if k == 'cotrending':
					#self.speak( "	  " + k + ', including:')
					#for evkey in tosave[k].keys():
						#self.speak( "		  " + evkey)

		mkdir(directory)
		filename = directory + str(self.__class__).split('.')[-1].split("'")[0] + '.npy'

		np.save(filename, tosave)
		self.speak("saving light curve to {directory}, including all its external variables".format(directory=directory))

	def load(self, directory):
		self.directory = directory
		filename = directory + str(self.__class__).split('.')[-1].split("'")[0] + '.npy'
		self.speak('trying to load TLC from .{0}'.format(filename))
		loaded = np.load(filename)[()]
		for key in loaded.keys():
			#self.speak( '	 ' + key)
			self.__dict__[key] = loaded[key]
