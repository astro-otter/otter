Using the RESTful API Directly
==============================

If you prefer to use the RESTful API directly, here is an example to pass an AQL query to the database:

::

   curl -u "user-guest:test" \
        -X POST \
        --header 'accept: application/json' \
	--header 'Content-Type: application/json' \
	--data '{"query": "FOR t IN transients LIMIT 10 RETURN t", "count": true, "batchSize": 100}' \
	'https://otter.idies.jhu.edu/api/_db/otter/_api/cursor'

The POST endpoint `_db/otter/_api/cursor` is the most important because it is how to post AQL queries to
the OTTER database and then get out their results. If, instead, you want to do a cone search, you can use

::

   curl -u "user-guest:test" \
        -X POST \
        --header 'accept: application/json' \
	--header 'Content-Type: application/json' \
	--data '{"query": "FOR t IN transients FILTER (t._ra >= @ra - @sep AND t._ra <= @ra + @sep AND t._dec >= @dec - @sep AND t._dec <= @dec + @sep) FILTER ASTRO::CONE_SEARCH(t._ra, t._dec, @ra, @dec, @sep) RETURN t", "count": true, "batchSize": 2, "bindVars": {"ra": 185.0, "dec": 12.0, "sep": 1.0}}' \
	'https://otter.idies.jhu.edu/api/_db/otter/_api/cursor'

For more details on AQL see https://docs.arangodb.com/3.12/aql/ and for more details on the Arangodb HTTP REST API see
https://docs.arangodb.com/3.12/develop/http-api/. If you have questions about other less commonly used endpoints, please reach out to the developers. In the future, we will try to add all
of them here, but they are all well documented in the arangodb gui as a Swagger ui under the "support" option > "Rest API" tab. (NOTE: At some point we
should try to figure out how to port that documentation over to here automatically, not sure if there is a way though)
