Installation notes
==================

The `foldamers <https://github.com/shirtsgroup/foldamers>`_ package will eventually be available for installation via `Anaconda <https://www.anaconda.com/>`_.  This will make resolution of software conflicts much easier.  However, at present, because the package has not been made public, Anaconda installation is not yet possible, and software conflicts must be resolved by the individual.

Here we provide installation instructions which have been tested on multiple platforms.

Dependencies for the foldamers package
--------------------------------------

Due to conflicts among dependencies for the foldamers package with other Python versions, **we recommend using Python version 3.6**.

The following is a list of software dependencies for the foldamers package, with recommended version numbers in parentheses:

1) `cg_openmm <https://github.com/shirtsgroup/cg_openmm>`_ (version 0.0)

   **Dependencies for the 'cg_openmm' software package:**
   2) `OpenMM <http://openmm.org/>`_ (version 7.3.1)
   3) `Yank <http://getyank.org/latest/>`_ (version 0.24.1)

4) `pymbar <https://github.com/choderalab/pymbar>`_ (version 3.0.3)
5) `MDTraj <http://mdtraj.org/1.9.3/>`_ (version 1.9.3)
6) `MSMBuilder <http://msmbuilder.org/3.8.0/>`_ (version 3.8)
7) `scikit-learn <https://scikit-learn.org/stable/>`_ (version 0.18.1)
8) `sklearn-genetic <https://github.com/manuel-calzolari/sklearn-genetic>`_ (version 0.2)
9) `kHelios <https://pubs.acs.org/doi/abs/10.1021/acs.jcim.6b00721>`_

Recommended installation steps
------------------------------

We recommend installation of `Anaconda <https://www.anaconda.com/>`_ prior to installation of the 'foldamers' package, as this makes resolution of conflicts between dependencies much easier.

We direct users that have not installed `Anaconda <https://www.anaconda.com/>`_ to the `Download page <https://www.anaconda.com/distribution/>`_ in order to select the appropriate version for your platform (Windows, Linux, Mac).  (It shouldn't matter which version of Anaconda is installed.)

The following installation steps are recommended for users that have already installed `Anaconda <https://www.anaconda.com/>`_ on their system:

1) Install 'sklearn-genetic' (used for genetic algorithm optimization of model parameters)::

	pip install sklearn-genetic

2) Create an Anaconda environment for Python version 3.6 (the most stable Python version for 'foldamers')::

	conda create -n foldamers python=3.6 mdtraj=1.9.3 openmm=7.3.1 pymbar=3.0.3 msmbuilder=3.8 scikit-learn=0.18.1 yank=0.21.2

3) Activate this new environment::

	conda activate foldamers

4) Clone and install a fresh copy of 'cg_openmm'::

	git clone https://github.com/shirtsgroup/cg_openmm.git
	cd cg_openmm
	python setup.py install

5) Clone and install a fresh copy of 'foldamers' to your system::

	git clone https://github.com/shirtsgroup/foldamers.git
	cd foldamers
	python setup.py install


