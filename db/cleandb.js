//db = require("@arangodb")
//db._dropDatabase("tide")

users = require("@arangodb/users")
users.remove("admin@tide")
users.remove("user@tide")
