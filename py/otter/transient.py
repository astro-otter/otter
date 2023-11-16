'''
Class for a transient, 
basically just inherits the dict properties with some overwriting
'''

class Transient(dict):

    def __init__(self, d={}):
        '''
        Overwrite the dictionary init
        
        Args:
            d [dict]: A transient dictionary
        '''
        super(Transient, self).__init__(d)

        self.bibcodes = [ref['name'] for ref in self['reference_alias']]
        self.readable_sources = [ref['human_readable_name'] for ref in self['reference_alias']]
        
    def __repr__(self, html=False):
        if not html:
            return f'''
            Transient: {self['name']['default_name']}\n
            -----------------------------------------\n
            RA: {self['coordinate']['equitorial'][0]['ra']}\n
            Dec: {self['coordinate']['equitorial'][0]['dec']}\n
            Sources: {self.readable_sources}
            '''
        else:
            
