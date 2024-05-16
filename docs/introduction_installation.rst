Installation Quickstart
-----------------------

User Installation
^^^^^^^^^^^^^^^^^

For non-developers, the :code:`astro-otter` API can be installed using :code:`pip`:

::

   python3 -m pip install astro-otter

Then the data repository must also be downloaded from GitHub. This is a temporary solution
until we get the data hosted in a more permanent location. To download the dataset make sure
git is installed and then run

::

   export OTTER_ROOT=/path/to/where/to/clone
   git clone https://github.com/astro-otter/otterdb.git $OTTER_ROOT/otterdb
   cd $OTTER_ROOT/otter/scripts/
   python3 gen_summary_table.py --otterroot $OTTER_ROOT

Developer Installation
^^^^^^^^^^^^^^^^^^^^^^

For a local installation, typically for developers, follow these steps instead:

1. First set an :code:`OTTER_ROOT` environment variable and download the source code and data

::

   export OTTER_ROOT=/path/to/where/to/clone
   git clone https://github.com/astro-otter/otter.git $OTTER_ROOT/otter
   git clone https://github.com/astro-otter/otterdb.git $OTTER_ROOT/otterdb

2. Install otter from the source code

::

   cd $OTTER_ROOT/otter
   python -m pip install -e .

3. Process the data to build the local "database" (although it is really just a directory). Then, you can build the "database" by running the following commands:

::

   cd $OTTER_ROOT/otter/scripts/
   python3 gen_summary_table.py --otterroot $OTTER_ROOT
