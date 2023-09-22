'''
The catalog pages
'''
from ast import literal_eval
from otter import TDECatalog
from flask import render_template, Blueprint, request
from . import plotSummary
from . import db

bp = Blueprint('catalog', __name__)

def fix(key):

    val = request.form[key]
    
    trueType = {'z':float,
                'minZ':float,
                'maxZ':float,
                'ra':str,
                'dec':str,
                'tdename':str,
                'spectraType':str,
                'photoType':str
                }
    
    if len(val) == 0:
        return None
    
    return trueType[key](val) 

# a simple page that says hello
@bp.route('/', methods=['GET', 'POST'])
def home():

    if request.method == 'POST':

        print(request.form)
        
        # here we can query the database differently
        tdename = fix('tdename')
        ra = fix('ra')
        dec = fix('dec')
        z = fix('z')
        minZ = fix('minZ')
        maxZ = fix('maxZ')
        photoType = fix('photoType')
        spectraType = fix('spectraType')

        cat = TDECatalog()
        tdes = cat.query(names=tdename,
                         ra=ra,
                         dec=dec,
                         z=z,
                         minZ=minZ,
                         maxZ=maxZ,
                         photometryType=photoType,
                         spectraType=spectraType
                         )
        
    else:
        tdes = db.get_db().tdes

    plotHTML = plotSummary.plotAll(tdes)
    return render_template('index.html', tdes=tdes, otherHTML=plotHTML)

@bp.route('/<tdename>')
def genTDEpages(tdename):
    '''
    Generate a tde page from a tde object
    '''
    tdes = db.get_db().tdes
    return render_template('tde.html', tde=tdes[tdename])
