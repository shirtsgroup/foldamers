Ensemble building tools
=======================

The foldamers package contains several tools for building conformational ensembles.  The `MDTraj <http://mdtraj.org>`_ and `MSMBuilder <http://msmbuilder.org/>`_ packages are leveraged to perform structural analyses in order to identify poses that are structurally similar.

Using MSMBuilder to generate conformational ensembles
-----------------------------------------------------

The foldamers package allows the user to apply K-means clustering tools from MSMBuilder in order to search for ensembles of poses that are structurally similar.  The centroid configurations for individual clusters are used as a reference, and ensembles are defined by including all structures that fall below an RMSD positions threshold (<2 Angstroms).

.. automodule:: ensembles.cluster
    :members: concatenate_trajectories, align_structures, get_cluster_centroid_positions

.. raw:: latex

    \newpage

Native structure-based ensemble generation tools
------------------------------------------------

The foldamers package allows the user to build "native" and "nonnative" structural ensembles, and to evaluate their energetic differences with the Z-score.  These tools require identification of a "native" structure.

.. automodule:: ensembles.ens_build
    :members: get_ensembles, get_ensembles_from_replica_positions, get_native_ensemble, get_nonnative_ensemble, get_native_structure, z_score

.. raw:: latex

    \newpage

Energy-based ensemble generation tools
--------------------------------------

The foldamers package allows the user to build structural ensembles that exhibit similar energies.  Shown below are tools that enable energy-based ensemble generation.

.. automodule:: ensembles.ens_build
    :members: get_ensemble, test_energy, improve_ensemble

.. raw:: latex

    \newpage

Writing and reading ensemble data from the 'foldamers' database
---------------------------------------------------------------

The foldamers package is designed to store the low-energy poses from simulation runs of new (previously un-modelled) coarse grained representations.  At present, the package does not enable storage of heteropolymers, in order to minimize the size of the database.  For homopolymers, the syntax for assigning directory names for coarse grained model data is as follows:

directory_name = str( "foldamers/ensembles/" + str(polymer_length) + "_" + str(backbone_length) + "_" + str(sidechain_length) "_" + str(sidechain_positions) + "_" + str(bb_bb_bond_length) + "_" + str(sc_bb_bond_length) + "_" + str(sc_sc_bond_length) )

For example, the directory name for a model with 20 monomers, all of which contain one backbone bead and one sidechain bead, and whose bond lengths are all 7.5 Angstroms, would be: "foldamers/ensembles/20_1_1_0_7.5_7.5_7.5".

The following functions are used to read and write ensemble data to the foldamers database (located in 'foldamers/ensembles').

.. automodule:: ensembles.ens_build
    :members: get_ensemble_directory, write_ensemble_pdb, get_pdb_list, get_ensemble_data

