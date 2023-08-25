**********************
Particle files and units
*************************

The particle directory is a place for files which describe a single particle, which could represent an atom, molecule or coarse-grained model.
These files are LAMMPS-inspired (https://lammps.sandia.gov/doc/read_data.html), but deviate significantly from LAMMPS as described in detail below.
Each file represents only a single particle (or molecule), and these files are used to define the types of particles (or molecules) that can exist in the simulation.
Each plugin many also contain a particle directory as additional examples and for use in tutorials.
The units used in these files (length, energy, etc) are not assumed (with one exception) and must be consistent with the units to initialize your simulations.

======
Units
======

The :doc:`charge plugin </plugin/charge/README>` is the only plugin that assumes units.
Otherwise, any other use of FEASST assumes that the user takes care to input values with a consistent set of units.
For example, the units of length for positions in the data file should match whatever is used for defining the Domain, and the :cpp:class:`FileXYZ <feasst::FileXYZ>` positions will also be consistent with that user choice.
The energy scale is given by epsilon-like parameters, as should be consistent with the beta given to :cpp:class:`ThermoParams <feasst::ThermoParams>`.
That is, beta has units of inverse energy.
Similarly, chemical_potential has units of energy.
Because pressure times volume is in units of energy, pressure is in units of energy per volume.

The :doc:`charge plugin </plugin/charge/README>` is an exception because the unit of charge is assumed to be elementary charge.
Because of an assumed conversion factor, the units of energy must be in kJ/mol, and therefore, beta must be given in units of mol/kJ.
Finally, the units of length must be in Angstroms.

===========
Data files
===========

Differences from LAMMPS
========================

LAMMPS data files are defined here: https://lammps.sandia.gov/doc/read_data.html

The data files used by FEASST have the following major differences:

- atoms in LAMMPS are analogous to sites in FEASST.
- molecules in LAMMPS are analogous to particles in FEASST.
- FEASST data files contain only one particle.
- all FEASST indices (types, sites, etc) begin with 0, not 1.
- characters are case sensitive,
- the number of spaces between characters does not matter,
- LMP Coeffs sections were replaced by Properties sections, with very different formatting,
- "2 dimensions" may be specified at the beginning of the file for a 2D simulation.
- FEASST data files contain no information about the :cpp:class:`Domain <feasst::Domain>` boundaries.
- The "Sites" section has the following three major differences: (1) these files describe only one particle so, unlike LAMMPS, molecule index is not included, (2) charge is a site-type property, and thus is not given for each site and (3) wrapping is not included.

Site Properties
================

Each site type may have a list of properties given by a label and value.

The format for this section is as follows:

[site type] [label_0] [value_0] ... [label_n] [value_n]

The labels are typically the ones described in :cpp:class:`ModelParams <feasst::ModelParams>`, such as sigma, epsilon, cutoff and charge.
The use of these parameters depends on the chosen models.
For example, :cpp:class:`LennardJones <feasst::LennardJones>` or :cpp:class:`SquareWell <feasst::SquareWell>`.

Other more complex models may utilize a custom defined ModelParam, such as :cpp:class:`patch_angle <feasst::PatchAngle>` as shown in /path/to/feasst/plugin/patch/particle/patch_one.fstprt.
In addition, the label "director" is described in :cpp:class:`VisitModelInnerPatch <feasst::VisitModelInnerPatch>` and utilized by :cpp:class:`Site <feasst::Site>`.

Sites
======

The format for this section is as follows:

[site index] [site type] [x-position] [y-position] [z-position]

Note that the site type matches the Site Properties.

Site Labels
=============

Site labels are used to attach a name to the site type, often used for visualization programs, etc.
Note: these are not currently implemented but serve as a place holder or for reference.

Bond, Angle and Dihedral Properties
======================================

The format for this section is as follows:

[bond type] [bond class name] [label_0] [value_0] ... [label_n] [value_n]

The labels are described for each bond class, such as :cpp:class:`RigidBond <feasst::RigidBond>`, :cpp:class:`AngleHarmonic <feasst::AngleHarmonic>` and :cpp:class:`DihedralTraPPE <feasst::DihedralTraPPE>`.

Bonds
=======

The format for this section is as follows:

[bond index] [bond type] [site index i] [site index j]

Angles
======

The format for this section is as follows:

[angle index] [angle type] [site index i] [site index j] [site index k]

Note that j is the vertex of the defined angle.
In 2D, angles are defined clockwise, such that angles ijk and kji are not the same.
For more information on the angle definition, see :cpp:func:`vertex_angle_radians <feasst::Position::vertex_angle_radians()>`

Dihedrals
==========

The format for this section is as follows:

[dihedral index] [dihedral type] [site index i] [site index j] [site index k] [site index l]

For more information on the dihedral definition, see :cpp:func:`torsion_angle_radians <feasst::Position::torsion_angle_radians()>`.

Comments
==========

For comments at the beginning of the file, begin each comment line with the "#" character.
Do not add comments anywhere else in the file.
Comments can only be one contiguous group of lines beginning with "#."
Always end comments with a blank line immediately after.

