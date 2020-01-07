.. _tool-installation:

Tool installation
=================

Install the conda package manager
---------------------------------

We will use the package/tool managing system |conda| to install some programs
that we will use during the course.
It is not installed by default, thus we need to install it first to be able to use it. 
Let us download |conda| first:


.. code-block:: bash

    # download latest conda installer
    $ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh


.. Note::
   Should the conda installer download fail. Please find links to alternative locations on the
   :doc:`../downloads` page.


Now lets install |conda|:

.. code-block:: bash

    # run the installer
    $ bash Miniconda3-latest-Linux-x86_64.sh
    

After you accepted the license agreement |conda| will be installed.
At the end of the installation you will encounter the following:


.. code-block:: bash

    ...
    installation finished.
    Do you wish the installer to initialize Miniconda3
    by running conda init? [yes|no]
    [no] >>> 


Please type **"yes"** here. 
This will add some code to your `.bashrc` init file, which is important to work with |conda| correctly.


.. Attention::
    **Please close and reopen the terminal, to complete the installation.**


After closing and re-opening the shell/terminal, we should be able to use the |conda| command:


.. code-block:: bash

    $ conda update --yes conda


Installing conda channels to make tools available
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Different tools are packaged in what |conda| calls channels.
We need to add some channels to make the bioinformatics and genomics tools
available for installation:


.. code-block:: bash
    
    # Install some conda channels
    # A channel is where conda looks for packages
    $ conda config --add channels defaults
    $ conda config --add channels bioconda 
    $ conda config --add channels conda-forge     


.. Attention::
    The order of adding channels is important. Make sure you use the shown order of commands.

   
Creating environments
---------------------

We create a |conda| environment for some tools.
This is useful to work **reproducible** as we can easily re-create the tool-set with the same version numbers later on.


.. code-block:: bash

    $ conda create -n ngs python=3
    # activate the environment
    $ conda activate ngs

    
So what is happening when you type ``conda activate ngs`` in a shell.
The ``PATH`` variable of your shell gets temporarily manipulated and set to:


.. code-block:: bash

   $ echo $PATH
   /home/guest/miniconda3/bin:/home/guest/miniconda3/condabin:...
   $ conda activate ngs
   $ echo $PATH
   /home/guest/miniconda3/envs/ngs/bin:/home/guest/miniconda3/condabin: ...


Now it will look first in your environment's bin directory but afterwards in the general conda bin (``/home/guest/miniconda3/condabin``).
So basically everything you install generally with conda (without being in an environment) is also available to you but gets overshadowed if a similar program is in ``/home/guest/miniconda3/envs/ngs/bin`` and you are in the ``ngs`` environment.


Install software
----------------

To install software into the activated environment, one uses the command ``conda install``.

.. code-block:: bash
         
    # install more tools into the environment
    $ conda install package


.. note::
   To tell if you are in the correct conda environment, look at the command-prompt.
   Do you see the name of the environment in round brackets at the very beginning of the prompt, e.g. (ngs)?
   If not, activate the ``ngs`` environment with ``conda activate ngs`` before installing the tools.

    
General conda commands
----------------------

.. code-block:: bash

    # to search for packages
    $ conda search [package]
    
    # To update all packages
    $ conda update --all --yes

    # List all packages installed
    $ conda list [-n env]

    # conda list environments
    $ conda env list

    # create new env
    $ conda create -n [name] package [package] ...

    # activate env
    $ conda activate [name]

    # deactivate env
    $ conda deactivate
