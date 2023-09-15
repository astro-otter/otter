'''
Simple TDE class with information about an individual TDE
'''

class TDE:

    def __init__(self, name, ra, dec, z, sources):

        self.name = name
        self.ra = ra
        self.dec = dec
        self.z = z
        self.sources = sources
        if self.sources is not None:
            self.fancySources = self._fancySources()
        else:
            self.fancySources = ''

        self.path = self._getpath()
        
    def _fancySources(self):

        for s in self.sources:
            print(s)

    def _getpath(self):
        return '' # for now
            
    def __str__(self):
        return f'''
        TDE: {self.name}
        --------------------------------
        RA      : {self.ra}
        DEC     : {self.dec}
        Z       : {self.z}
        Sources : {self.fancySources}
        '''
        
