
=================================
Windows installation instructions
=================================

Since chemicalite is not available on the Windows platform, a workaround is needed to use HIPPO. In this case using a linux virtual machine within Docker provides an accessible route to using HIPPO. Follow the steps below

Setting up Docker
=================

Docker is the industry standard way to create and share application *containers*, which are standalone Linux emulations with customisable software environments. 

1. Create a Docker_ account
2. Install `Docker Desktop`_

Setting up the Jupyter extension
================================

HIPPO works best within a JupyterLab interface. This is provided on the Docker Extension Marketplace.

1. In Docker Desktop: go to *Add Extensions*
2. Search for `Jupyter Notebook Scientific Python Stack Extension` and install it
3. Open the newly installed extension

Installing HIPPO
================

1. Open a terminal in the JupyterLab interface
2. Follow the linux/mac installation procedure for HIPPO:

::

	pip install --upgrade hippo-db
	conda install -c conda-forge chemicalite

The HIPPO module will now be available to use in any Jupyter notebook running inside the Docker extension:

::

	import hippo
	animal = hippo.HIPPO('test', 'test.sqlite')

N.B. the Docker Extension is running inside a container with a persistent storage volume and software environment that will survive a Docker Desktop restart, but if you uninstall the extension the data (and HIPPO installation) will be lost.

Further resources
=================

See the other pages in this documentation for help with HIPPO, and further information on the Jupyter Docker Extension is documented here_.

.. _Docker: https://www.docker.com/get-started/
.. _Docker Desktop: https://docs.docker.com/desktop/install/windows-install/
.. _here: https://www.docker.com/blog/getting-started-with-jupyterlab-as-a-docker-extension/
