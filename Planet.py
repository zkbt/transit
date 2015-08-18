from Parameters import Parameters
import eb
import numpy as np
import zachopy.units

class Planet(Parameters):
	'''Parameters (both set and calculated) describing a planet.'''
	def __init__(self, **kwargs):

		# define defaults
		Parameters.__init__(self, 	J=0.0, \
									k=0.1, \
									rsovera=1.0/10.0, \
									b=0.0, \
									period=10.0, \
									t0=2456000.0, \
									dt=0.0, \
									semiamplitude=0.0, \
									esinw=0.0, \
									ecosw=0.0)

		# overwrite defaults with input keywords
		Parameters.__init__(self, **kwargs)

		# set up the parameter constraints
		self.J.parinfo['limited'] = [True, False]
		self.J.parinfo['limits'] = [0, 100]

		self.b.parinfo['limited'] = [True, True]
		self.b.parinfo['limits'] = [0.0, 1.0]

		#self.q.parinfo['limited'] = [True, False]
		#self.q.parinfo['limits'] = [0, 100]

		self.period.parinfo['limited'] = [True, False]
		self.period.parinfo['limits'] = [0, 100]

		self.t0.parinfo['limited'] = [True, False]
		self.t0.parinfo['limits'] = [0, 1000000000]



	@property
	def surface_brightness_ratio(self):
		return self.J


	@property
	def mass_ratio(self):
		return self.q

	@property
	def rp_over_rs(self):
		return self.k

	@property
	def depth(self):
		return self.rp_over_rs.value**2

	@property
	def a_over_rs(self):
		return 1.0/self.rsovera.value

	@property
	def rsum_over_a(self):
		#by definition of k = rp/rs
		k = self.k.value
		rsovera = self.rsovera.value
		return (1.0 + k)*rsovera

	@property
	def e(self):
		return np.sqrt(self.esinw.value**2 + self.ecosw.value**2)

	@property
	def cosi(self):
		# from Winn (2010)
		b = self.b.value
		rsovera = self.rsovera.value
		e = self.e
		esinw = self.esinw.value
		return b*rsovera*(1.0 + esinw)/(1-e*e)

	@property
	def sini(self):
		return np.sqrt(1 - self.cosi**2)

	@property
	def stellar_density(self):
		'''requires an estimate of the mass ratio!'''
		rsovera = self.rsovera.value
		period = self.period.value*zachopy.units.day
		q = 0.0#self.q.value
		return 3*np.pi/zachopy.units.G/period**2*(1.0/rsovera)**3/(1.0 + q)

	@property
	def duration(self):
		'''duration from 1st to 4th contact (in days)'''
		phasecontacts = self.contacts()
		return (phasecontacts[1] - (phasecontacts[0] - 1))*self.period.value

	def contacts(self):
		return eb.phicont(self.esinw.value, self.ecosw.value, self.cosi, self.rsum_over_a)

	def duration_total(self):
		period = self.period.value
		rsovera = self.rsovera.value
		k = self.k.value
		b = self.b.value
		sini = self.sini
		return period/np.pi*np.arcsin(rsovera/sini*np.sqrt((1.0 + k)**2 - b**2))

	def duration_full(self):
		period = self.period.value
		rsovera = self.rsovera.value
		k = self.k.value
		b = self.b.value
		sini = self.sini
		return period/np.pi*np.arcsin(rsovera/sini*np.sqrt((1.0 - k)**2 - b**2))



	def thisepoch(self, bjd):
		return np.round((bjd - self.t0.value)/self.period.value).astype(np.int)

	def thismidtransit(self, bjd):
		return np.round((bjd - self.t0.value)/self.period.value)*self.period.value + self.t0.value

	def timefrommidtransit(self, bjd):
		'''Calculate the time from the assigned mid-transit time.'''
		phasedtime = (bjd - self.thismidtransit(bjd))
		mask = phasedtime > 0.5*self.period.value
		phasedtime[mask] -= self.period.value
		return phasedtime
