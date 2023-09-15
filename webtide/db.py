from pyArango.connection import *

from flask import current_app, g

def get_db():
    if 'db' not in g:
        c = Connection(
            username = os.environ['arango_user'],
            password = os.environ['arango_pass']
        )
        g.db = c['tide']
        
    return g.db


def close_db(e=None):
    db = g.pop('db', None)

    if db is not None:
        db.close()

def init_db():
    db = get_db()

    tdes = db['tdes'].fetchAll(rawResults=True)

def init_app(app):
    app.teardown_appcontext(close_db)
