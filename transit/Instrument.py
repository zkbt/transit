from .Parameters import Parameters
import numpy as np

class Instrument(Parameters):
	'''An Instrument object to keep track of all the decorrelation nuisance parameters.'''

	def __init__(self, tlc=None, order=1, torder=0, directory=None, **kwargs):
		# initialize the Instrument, using the extra columns in the light curve

		dict = {}

		dict['rescaling'] = 1.0

		# a constant out of transit flux level
		try:
			dict['C'] = np.median(tlc.flux)
		except:
			dict['C'] = 1.0

		if tlc:
			# include terms of (ev)^power for all the external variables available
			for power in np.arange(order)+1:
				for evkey in tlc.cotrending.keys():
					dict[evkey + '_tothe{0:1d}'.format(power)] = 0.0

			# include (t)^power terms for time (CURRENTLY SET UP ONLY FOR SINGLE TRANSITS!)
			for power in np.arange(torder)+1:
				dict['t_tothe{0:1d}'.format(power)] = 0.0

			# include other custom instrument parameters
			for k in kwargs.keys():
				dict[k] = kwargs[k]
		#except AttributeError:
		#	pass

		print(dict)
		# include all the parameters explicitly defined here
		Parameters.__init__(self, directory=directory, **dict)


	def normalize(self, x, ok=None):
		if ok is None:
			s = np.std(x)
		else:
			s = np.std(x[ok])
		if s == 0:
			s = np.mean(x)
			if s == 0.0:
				s = 1.0
		return (x - np.mean(x))/s

	#@profile
	def model(self, tlc):

		keys = self.keys
		m = np.zeros_like(tlc.bjd)
		for k in keys:

			parameter = self.__dict__[k]
			if parameter.name == 'C':
				m += parameter.value

			if 'tothe' in parameter.name:
				# parse the parameter name string
				chunks = parameter.name.split('_tothe')

				# pull out the template vector appropriate
				name = chunks[0]

				try:
					# try to return a prenormalized timeseries
					ev = tlc._normalized[k]
					assert(len(ev) == len(tlc.flux))
				except (AttributeError,KeyError,AssertionError):
					# if that's impossible, calculate and store a normalized timeseries for this key

					# make the _normalized dictionary exists
					try:
						tlc._normalized
					except AttributeError:
						tlc._normalized = {}

					# check whether the external variable here is time
					if k == 't':
						x = tlc.bjd
					else:
						x = tlc.cotrending[name]
					tlc._normalized[k] = self.normalize(x, ok=(tlc.bad == False))
					ev = tlc._normalized[k]

				# figure out what power the template should be raised to
				power = np.int(chunks[-1])

				# add the template (to the power) into the model
				if parameter.value != 0.0:
					if (np.isfinite(ev)).all():
						m += parameter.value*ev**power

		return m
		'''if True:
			t = bjd - np.mean(bjd)
			return (1.0 + self.C.value + self.T1.value*t + self.T2.value*t**2)
		else:
			return np.ones_like(bjd)'''
