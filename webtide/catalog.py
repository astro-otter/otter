
'''
The home page
'''

from flask import render_template, Blueprint
from . import db

bp = Blueprint('catalog', __name__)

# a simple page that says hello
@bp.route('/')
def home():
    database = db.get_db()
    return render_template('index.html', tdes=database.tdes)

@bp.route('/<tdename>')
def genTDEpages(tdename):
    '''
    Generate a tde page from a tde object
    '''
    tdes = db.get_db().tdes
    return render_template('tde.html', tde=tdes[tdename])
    
