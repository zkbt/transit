from IndividualPlots import *

class GroupedPlots(IndividualPlots):
    def identifier(self, tlc):
        s = tlc.telescope
        if 'MEarth' in s:
            s = 'MEarth'
        return s
