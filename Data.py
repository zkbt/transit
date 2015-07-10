from imports import *
class Data(Talker):
	def __init__(self, **kwargs):
		Talker.__init__(self, **kwargs)


	def setupColors(self, color='eye', minimumuncertainty=0.001):
		'''Method to set the line and point colots for this light curve.
			color = 'eye': make colors as they would appear to the human eye
			[add other options (e.g. specify a color specific color)]'''

		# set up the appropriate colors to use
		#try:
		#	self.colors
		#except:

		if True:
			self.colors = {}
			#if color=='eye':
			#	self.colors['points'] = zachopy.color.nm2rgb([self.left/10, self.right/10], 0.25)
			#	self.colors['lines'] = zachopy.color.nm2rgb([self.left/10, self.right/10], intensity=3.0)
			#else:
			if True:
				self.colors['lines'] = color
				r, g, b = zachopy.color.name2color(color.lower())
				rgba = np.zeros((self.n, 4))
				rgba[:,0] = r
				rgba[:,1] = g
				rgba[:,2] = b
				if minimumuncertainty is None:
					minimumuncertainty = np.min(self.effective_uncertainty)
				weights = np.minimum((minimumuncertainty/self.effective_uncertainty)**2, 1)
				over = weights > 1
				weights[over] = 1
				rgba[:,3] = weights
				self.colors['points'] = rgba
