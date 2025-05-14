Installation Quickstart
-----------------------

User Installation
^^^^^^^^^^^^^^^^^

For non-developers, the :code:`astro-otter` API can be installed using :code:`pip`:

::

   python3 -m pip install astro-otter

Then, from a python terminal you can do things like

::

   import otter
   db = otter.Otter()
   db.query(names="18hyz")

Developer Installation
^^^^^^^^^^^^^^^^^^^^^^

For a local installation of everything, typically for developers, follow these steps instead:

1. First, install some other dependencies:
   a. ADS API: https://ads.readthedocs.io/en/latest/#getting-started
   b. Docker (install AND run):

2. Then, set the following environment variables in your .bashrc/.zshrc file:

::

   export OTTER_ROOT=/path/to/where/to/clone
   export ARANGO_URL=<insert url>
   export ARANGO_ROOT_PASSWORD=<insert password>
   export VETTING_PASSWORD=<insert password>
   export ARANGO_USER_PASSWORD=<insert password>
   export ADS_DEV_KEY=<insert api key>
   export TNS_BOT_ID=<bot id>
   export TNS_BOT_NAME=<bot name>
   export TNS_API_KEY=<api key>
   export OTTER_WEB_BASE_URL="/otter"

If you don't know what any of the values should be, reach out to an existing developer

3. Source your .bashrc/.zshrc

::

   source ~/.bashrc

4. Clone the relevant repositories

::

   git clone https://github.com/astro-otter/otter.git $OTTER_ROOT/otter
   git clone https://github.com/astro-otter/otterdb.git $OTTER_ROOT/otterdb
   git clone https://github.com/astro-otter/otter-web.git $OTTER_ROOT/otter-web
   git clone https://github.com/astro-otter/otter-docker.git $OTTER_ROOT/otter-docker

5. Install the otter python API from source and install our pre-commit
   hooks that help maintain code quality.

::

   cd $OTTER_ROOT/otter
   python -m pip install -e .
   pre-commit install


6. Get the database up and running:
   a. First, make sure docker is running
   b. If you have an arm64 cpu (like a Mac Silicon chip), you will need to build the docker container on
      that cpu architecture. From the otterdb directory run:

      ::
	 docker build -t otterdb:v0.3.6 .

   c. Then run the docker container from dockerhub (if you built it locally, make sure you change the
      name at the bottom of this command)

      ::
	 docker run \
	 -e ARANGO_ROOT_PASSWORD=$ARANGO_ROOT_PASSWORD \
	 -e VETTING_PASSWORD=$VETTING_PASSWORD \
	 -e VETTING_USER="vetting-user" \
	 -e ARANGO_USER_USERNAME="user-guest" \
	 -e ARANGO_USER_PASSWORD=$ARANGO_USER_PASSWORD \
	 -e DB_LINK_PORT_8529_TCP="http://127.0.0.1:8529" \
	 -p 8529:8529 \
	 noahfranz13/otterdb:v0.3.6

   d. Open a new terminal. Then, add a copy of the data to the database from the master copy at SciServer.
      From the otterdb directory, run

      ::
	 python3 import_from_sciserver.py

7. Test that arangodb is running correctly and the data is imported by navigating to http://localhost:8529
   in a browser and logging in as user "root" and the password set by $ARANGO_ROOT_PASSWORD

8. Install the web application by running

   ::
      cd $OTTER_DIR/otter-web
      python3 -m pip install -e .

9. Run the web application in developer mode. From the otter-web directory

   ::
      export ARANGO_URL="http://localhost:8529" && python3 start.py

10. Test that the website is running by navigating to http://localhost:8080/otter

11. Assuming everything is working as expected, you are now able to change things as you wish!
