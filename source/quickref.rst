Quick command reference
=======================

Shell commands
--------------

.. code:: bash

    # Where in the directory tree am I?
    pwd

    # List the documents and sub-directories in the current directory
    ls

    # a bit nicer listing with more information
    ls -laF

    # Change into your home directory
    cd ~

    # Change back into the last directory
    cd -

    # Change one directory up in the tree
    cd ..

    # Change explicitly into a directory "temp"
    cd temp

    # Quickly show content of a file "temp.txt"
    # exist the view with "q", navigate line up and down with "k" and "j"
    less temp.text

    # Show the beginning of a file "temp.txt"
    head temp.txt

    # Show the end of a file "temp.txt"
    tail temp.txt

General conda commands
----------------------

.. code:: bash

    # To update all packages
    conda update --all --yes

    # List all packages installed
    conda list [-n env]

    # conda list environments
    conda env list

    # create new env
    conda create -n [name] package [package] ...

    # activate env
    conda activate [name]

    # deavtivate env
    conda deactivate
