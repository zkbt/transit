from Parameters import Parameters
import numpy as np

class Instrument(Parameters):
	'''An Instrument object to keep track of all the decorrelation nuisance parameters.'''

	def __init__(self, tlc=None, order=1, torder=0, directory=None, **kwargs):
		# initialize the Instrument, using the extra columns in the light curve

		dict = {}

		# a constant out of transit flux level
		try:
			dict['C'] = np.median(tlc.flux)
		except:
			dict['C'] = 1.0

		try:
			# include terms of (ev)^power for all the external variables available
			for power in np.arange(order)+1:
				for evkey in tlc.externalvariables.keys():
					dict[evkey + '_tothe{0:1d}'.format(power)] = 0.0

			# include (t)^power terms for time (CURRENTLY SET UP ONLY FOR SINGLE TRANSITS!)
			for power in np.arange(torder)+1:
				dict['t_tothe{0:1d}'.format(power)] = 0.0
		except:
			pass

		# include all the parameters explicitly defined here
		Parameters.__init__(self, directory=directory, **dict)


	def normalize(self, x):
		return (x - np.mean(x))/np.std(x)

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

				except (AttributeError,KeyError):
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
						x = tlc.externalvariables[name]
					tlc._normalized[k] = self.normalize(x)
					ev = tlc._normalized[k]

				# figure out what power the template should be raised to
				power = np.int(chunks[-1])

				# add the template (to the power) into the model
				if parameter.value != 0:
					m += parameter.value*ev**power
		assert(np.isfinite(m).all())
		return m
		'''if True:
			t = bjd - np.mean(bjd)
			return (1.0 + self.C.value + self.T1.value*t + self.T2.value*t**2)
		else:
			return np.ones_like(bjd)'''
