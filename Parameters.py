import numpy as np
from Parameter import Parameter
from zachopy.Talker import Talker

class Parameters(Talker):
	'''An object to keep track of an array of parameter.'''
	def __init__(self, directory=None, **kwargs):
		Talker.__init__(self)
		if directory is not None:
			self.load(directory)
		else:
			for key, value in kwargs.iteritems():
				# set up as a parameter
				self.__dict__[key] = Parameter(key, value)

				# fix everything initially (probably unnecessary?)
				self.__dict__[key].fix()

	@property
	def keys(self):
		k = self.__dict__.keys()
		k.remove('mute')
		k.remove('pithy')
		k.remove('line')
		try:
			k.remove('prefix')
		except:
			pass
		return k

	def __str__(self):
		s = ''
		for parameter in self.__dict__.values():
			s += "{0:10}\n".format(parameter)
		return s

	def toArray(self):
		d = self.__dict__
		array = np.zeros(len(d))
		for k, count in zip(sorted(d.keys()), range(len(array))):
			array[count] = d[k].value
		return array

	def fromArray(self, array):
		d = self.__dict__
		assert(len(d) == len(array))
		for k, count in zip(sorted(d.keys()), range(len(array))):
			d[k].value = array[count]

	def load(self, directory):
		filename = directory + str(self.__class__).split('.')[-1].split("'")[0] + '.npy'
		loaded = np.load(filename)[()]
		for key in loaded.keys():
			self.__dict__[key] = loaded[key]
		self.speak('loaded from {0}'.format(filename))



	def save(self, directory):
		'''Tool to save a parameters to a directory.'''
		tosave = {}
		for k in self.keys:
			tosave[k] = self.__dict__[k]
			self.speak( '{0}'.format(self.__dict__[k]))
		filename = directory + str(self.__class__).split('.')[-1].split("'")[0] + '.npy'
		np.save(filename, tosave)
		self.speak('saved to {0}'.format(filename))
