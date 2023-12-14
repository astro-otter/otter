const { db } = require("@arangodb");

// create the database
db._createDatabase("tide")
db._useDatabase("tide")

// create the tde collection
db._create("tdes")

// create our users
const users = require('@arangodb/users');
users.save('admin@tide', 'insecure'); // for full access
users.save('user@tide', 'password'); // for the api

// give admin full access
users.grantDatabase('admin@tide', 'tide', 'rw')
users.grantCollection('admin@tide', 'tide', 'tdes', 'rw')

// give users only read access
users.grantDatabase('user@tide', 'tide', 'ro')
users.grantCollection('user@tide', 'tide', 'tdes', 'ro')
