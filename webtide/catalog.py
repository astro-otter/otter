'''
The catalog pages
'''
from tidecat import TDECatalog
from flask import render_template, Blueprint, request
from . import plotSummary
from . import db

bp = Blueprint('catalog', __name__)

# a simple page that says hello
@bp.route('/', methods=['GET', 'POST'])
def home():

    if request.method == 'POST':
        # here we can query the database differently
        tdename = request.form['tdename']
        catalog = TDECatalog()
        tdes = {tdename: catalog.tdes[tdename]}
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
