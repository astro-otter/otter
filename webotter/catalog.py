'''
The catalog pages
'''
from ast import literal_eval
from astropy.coordinates import SkyCoord
from otter import Otter
from flask import render_template, Blueprint, request
from . import plotgen
from . import db


bp = Blueprint('catalog', __name__)

def fix(key):

    if key in {'spectra', 'photometry'}:
        val = request.form.get(key)
        return val == 'on'
    else:
        val = request.form[key]
    
        trueType = {'minZ':float,
                    'maxZ':float,
                    'ra':str,
                    'dec':str,
                    'tdename':str,
                    'spectra':bool,
                    'photometry':bool,
                    'searchRadius':float
                    }
    
        if len(val) == 0:
            return None
    
        return trueType[key](val) 

# a simple page that says hello
@bp.route('/', methods=['GET', 'POST'])
def home():

    if request.method == 'POST':
        # here we can query the database differently
        tdename = fix('tdename')
        ra = fix('ra')
        dec = fix('dec')
        minZ = fix('minZ')
        maxZ = fix('maxZ')
        hasphoto = fix('photometry')
        hasspec = fix('spectra')
        searchRadius = fix('searchRadius')
        
        cat = Otter()
        
        if ra is not None and dec is not None:
            coord = SkyCoord(ra, dec, unit=('hourangle', 'deg'))
        else:
            coord=None

        if minZ is None:
            minZ = 0
            
        tdes = cat.query(names=tdename,
                         coords=coord,
                         minZ=minZ,
                         maxZ=maxZ,
                         radius=searchRadius,
                         hasSpec=hasspec,
                         hasPhot=hasphoto
                         )
        
    else:
        tdes = db.get_db().query() # get all data

    plotHTML = plotgen.plotCatalogSummary(tdes)
    return render_template('index.html', tdes=tdes, otherHTML=plotHTML)

@bp.route('/<tdename>')
def genTDEpages(tdename):
    '''
    Generate a tde page from a tde object
    '''
    tdes = db.get_db().query(names=tdename)[0] # should be safe to only take the first
    return render_template('tde.html', tde=tdes, plotPhotometry=plotgen.plotPhot,
                           plotSpectra=plotgen.plotSpec)
