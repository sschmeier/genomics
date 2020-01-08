# Computational Genomics Tutorial pages

[![Build Status](https://travis-ci.org/sschmeier/genomics.svg?branch=master)](https://travis-ci.org/sschmeier/genomics) [![PyUp](https://pyup.io/repos/github/sschmeier/genomics/shield.svg)](https://pyup.io/repos/github/sschmeier/genomics/) [![ReadTheDocs](https://readthedocs.org/projects/genomics/badge/?version=latest)](https://genomics.readthedocs.io/en/latest/?badge=latest)

A computational genomics tutorial with continuous integration, using [Travis-CI](https://travis-ci.org/),
and deployment, using [Netlify](https://www.netlify.com/), and [sphinx](http://www.sphinx-doc.org/).

The production version of the deployed site: 

 - [https://genomics.netlify.com](https://genomics.netlify.com)
 - [https://genomics.readthedocs.io/](https://genomics.readthedocs.io/)
 - [https://genomics.sschmeier.com](https://genomics.sschmeier.com)


## Details

- Netlify runs Ubuntu Xenial 16.04 build, with a python3.5 installed.
- Thus, need to use py3.5

## Installing Locally

1. Set up a [python virtual environment](https://packaging.python.org/guides/installing-using-pip-and-virtualenv/)
   named `venv` => `pip install virtualenv; virtualenv --python=python3.5 venv`.
2. Activate the `venv` environment => `source venv/bin/activate`
3. Install the dependencies inside of it by running  `pip install -r requirements.txt`.
4. Run `make htmlwatch`.
5. Edit your rst-files.
6. Commit changes
7. Bump version using, e.g. `bump2version patch`
8. Push changes to remote
