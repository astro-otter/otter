Using the RESTful API Directly
==============================

If you prefer to use the RESTful API directly, here is an example to pass an AQL query to the database:

::

   curl -u "user-guest" -X POST --header 'accept: application/json' --data-binary @- --dump - 'http://localhost:8529/_db/otter/_api/cursor' <<'EOF'
   {
     "query": "FOR p IN transients LIMIT 2 RETURN p",
     "count": true,
     "batchSize": 2
   }
   EOF

The POST endpoint `_db/otter/_api/cursor` is the most important because it is how to post AQL queries to the OTTER database
and then get out their results. For more details on AQL see https://docs.arangodb.com/3.12/aql/ and for more details on the Arangodb HTTP REST API see
https://docs.arangodb.com/3.12/develop/http-api/.

If you have questions about other less commonly used endpoints, please reach out to the developers. In the future, we will try to add all
of them here, but they are all well documented in the arangodb gui as a Swagger ui under the "support" option > "Rest API" tab. (NOTE: At some point we
should try to figure out how to port that documentation over to here automatically, not sure if there is a way though)
