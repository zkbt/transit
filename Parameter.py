import numpy as np

class Parameter(object):
	'''An object to handle everything related to a parameter, aiming to be useful either as an input to mpfit or to my MCMC codes.'''


	def __init__(self, name=None, value=None):
		"""
		'value' - the starting parameter value (but see the START_PARAMS
				 parameter for more information).
		"""
		self.value = value

		"""
		'fixed' - a boolean value, whether the parameter is to be held
				 fixed or not.  Fixed parameters are not varied by
				 MPFIT, but are passed on to MYFUNCT for evaluation.
		"""
		self.fixed = True

		"""
		'limited' - a two-element boolean array.  If the first/second
				   element is set, then the parameter is bounded on the
				   lower/upper side.  A parameter can be bounded on both
				   sides.  Both LIMITED and LIMITS must be given
				   together.
		"""
		self.limited = [False, False]

		"""
		'limits' - a two-element float array.  Gives the
				  parameter limits on the lower and upper sides,
				  respectively.  Zero, one or two of these values can be
				  set, depending on the values of LIMITED.  Both LIMITED
				  and LIMITS must be given together.
		"""
		self.limits = [-1.0e6, 1.0e6]

		"""
		'parname' - a string, giving the name of the parameter.  The
				   fitting code of MPFIT does not use this tag in any
				   way.  However, the default iterfunct will print the
				   parameter name if available.
		"""
		self.parname = name

		"""
		'step' - the step size to be used in calculating the numerical
				derivatives.  If set to zero, then the step size is
				computed automatically.  Ignored when AUTODERIVATIVE=0.
		"""
		self.step = 0.0


		"""
		'mpprint' - if set to 1, then the default iterfunct will print the
				   parameter value.  If set to 0, the parameter value
				   will not be printed.  This tag can be used to
				   selectively print only a few parameter values out of
				   many.  Default: 1 (all parameters printed)
		"""
		self.mpprint = True


		# an attribute to keep track of an uncertainty estimate (currently symmetric)
		self.uncertainty = 0.0
		self.name = self.parname


		self.independent = False

	def __repr__(self):
		return "<Parameter {name}|{value}\pm{uncertainty}>".format(**self.__dict__)

	@property
	def parinfo(self):
		return self.__dict__

	@parinfo.setter
	def parinfo(self, value):
		for key in value.keys():
			self.__dict__[key] = value[key]

	def float(self, value=None, limits=None, step=None, shrink = 100.0):
		'''Helper function to cause a parameter to be able to float.'''
		self.fixed = False
		if value is not None:
			self.value = value

		if limits is not None:
			self.limits = limits
			self.limited = [True,True]

		if self.value != 0:
			self.step = np.minimum(np.abs(self.limits[1] - self.limits[0]), self.value/shrink)
		if self.value == 0:
			self.step = np.abs(self.limits[1] - self.limits[0])/shrink


	def fix(self, value=None):
		'''Helper function to fix a parameter.'''
		self.fixed = True
		if value is not None:
			self.value=value


	def __float__(self):
		return self.value

	def exponent(self, number):
		s = '{0}'.format(number)
		if 'e' in s:
			return s.replace('e', r'\times 10^{') + '}'
		else:
			return s

	def string(self):
		ndigits = -np.round(np.log10(self.uncertainty)).astype(np.int)+1
		s = '${value} \pm {uncertainty}$'.format(value=self.exponent(np.round(self.value, decimals=ndigits)), uncertainty=self.exponent(np.round(self.uncertainty, decimals=ndigits)))
		return s
		#if ndigits < 0:
		#	form = '.{0:.0f}f'.format(np.round(-ndigits + 1))

			####### PICK UP HERE!
		"""
			if keyword_set(auto) then begin
				ndig =alog10(mean([(pos_err), (neg_err)])) > (-100)
				if ndig lt 0 then begin
					f_errors = '(D20.' + strcompress(round(-ndig+1), /remo) + ')'
				endif else begin
					f_errors = '(I20)'; + strcompress(round(-ndig+1), /remo) + ')'
				endelse
				f_central = f_errors
				threshold = 0.1
				if abs((pos_err - sym_err)/sym_err) lt threshold and abs((neg_err-sym_err)/sym_err) lt threshold then sym =1
			endif"""
