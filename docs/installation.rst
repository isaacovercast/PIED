.. _sec-installation:

============
Installation
============

``PIED`` requires Python >= 3.5. Installation is facilitated by the conda package
management system.

1. Download `miniconda <https://conda.io/miniconda.html>`_ and run the installer: ``bash Miniconda*``
2. Create a separate `conda environment <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_ to install PIED into:

.. code:: bash

    conda create -n PIED
    conda activate PIED

3. Install:

.. code:: bash

    conda install -c conda-forge -c PIED PIED

4. Test:

.. code:: bash

   PIED -v

Installation issues can be reported on the `PIED github <https://github.com/isaacovercast/PIED>`_.
