Coarse grained model utilities
==============================

This page details the functions and classes in src/cg_model/cgmodel.py

The 'basic_cgmodel' function to build coarse grained oligomers
--------------------------------------------------------------

Shown below is the 'basic_cgmodel' function, which requires only a minimal set of input arguments to build a coarse grained model.  Given a set of input arguments this function creates a CGModel() class object, applying a set of default values for un-defined parameters.

.. automodule:: cg_model.cgmodel
    :members: basic_cgmodel

Full 'CGModel' class to build/model coarse grained oligomers
------------------------------------------------------------

Shown below is a detailed description of the full 'cgmodel' class object.

.. automodule:: cg_model.cgmodel
    :members: CGModel
    :inherited-members: CGModel
    :exclude-members: basic_cgmodel, get_parent_bead

Other coarse grained model utilities
------------------------------------

.. automodule:: cg_model.cgmodel
    :members: get_parent_bead
