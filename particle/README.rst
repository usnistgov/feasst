*************************
Particle files and units
*************************

The particle directory is a place for files which describe a single particle.
Particles could represent an atom, molecule, colloid or coarse-grained model.
These files are LAMMPS-inspired (https://docs.lammps.org/read_data.html), but deviate significantly from LAMMPS as described in detail below.
Each file represents only a single particle (or molecule), and these files are used to define the types of particles (or molecules) that can exist in the simulation.
Each plugin many also contain a particle directory as additional examples and for use in tutorials (list examples with the command "find /path/to/feasst/ -name '\*.fstprt'").

The units used in these files (length, energy, etc) are not assumed (with one exception) and must be consistent with the units to initialize your simulations.
The :doc:`charge plugin </plugin/charge/README>` is the only plugin that assumes units.
Otherwise, any other use of FEASST assumes that the user takes care to input values with a consistent set of units of their choosing.
For example, the units of length for positions in the data file should match whatever is used for defining the Domain, and the :cpp:class:`FileXYZ <feasst::FileXYZ>` positions will also be consistent with that user choice.
The energy scale is given by epsilon-like parameters, as should be consistent with the beta given to :cpp:class:`ThermoParams <feasst::ThermoParams>`.
That is, beta has units of inverse energy.
Similarly, chemical_potential has units of energy.
Because pressure times volume is in units of energy, pressure is in units of energy per volume.

The :cpp:class:`Potentials <feasst::Potential>` in :doc:`charge plugin </plugin/charge/README>` are exceptions because the unit of charge is assumed to be elementary charge.
Because of an assumed conversion factor, the units of energy must be in kJ/mol, and therefore, beta must be given in units of mol/kJ.
Finally, the units of length must be in Angstroms.

Site Properties
================

This section and all following sections begins with the Section name (in this case, "Site Properties"), followed by an empty line, then one line for each site type, and finally ends with another empty line.
Each site type line may have a list of properties given by a name and value.
The format for each line (e.g., each site type) in this section is as follows:

[site type name] [name_0]=[value_0] ... [name_n]=[value_n]

All names are strings and all values are floating point numbers.
After the site type name, the following names are typically the ones described in :cpp:class:`ModelParams <feasst::ModelParams>`, such as sigma, epsilon, cutoff and charge.
The use of these parameters depends on the chosen models.
For example, :cpp:class:`LennardJones <feasst::LennardJones>` or :cpp:class:`SquareWell <feasst::SquareWell>`.

Other more complex models may utilize a custom defined ModelParam, such as :cpp:class:`patch_angle <feasst::PatchAngle>` as shown in /path/to/feasst/plugin/patch/particle/patch_one.fstprt.
In addition, the name "director" is described in :cpp:class:`VisitModelInnerPatch <feasst::VisitModelInnerPatch>` and utilized by :cpp:class:`Site <feasst::Site>`.

Sites
======

As for all sections, this begins with "Sites", an empty line, a number of lines equal to the number of sites, and an empty line (or end of file).

The format for each line in this section is as follows:

[site name] [site type name] [x-position] [y-position] [z-position]

The site type must match one of those listed in Site Properties.
The z-position is not required if two-dimensional.

Bond, Angle and Dihedral Properties
======================================

As for all sections, these begin with "Bond Properties", "Angle Properties" or "Dihedral Properties", an empty line, a number of lines equal to the number of bond/angle/dihedral types, and an empty line (or end of file).

The format for each line in this section is as follows:

[bond type name] [bond class name] [name_0]=[value_0] ... [name_n]=[value_n]

All names are strings and all values are floating point numbers.
The names are described for each bond class, such as :cpp:class:`RigidBond <feasst::RigidBond>`, :cpp:class:`AngleHarmonic <feasst::AngleHarmonic>` and :cpp:class:`DihedralTraPPE <feasst::DihedralTraPPE>`.

Bonds
=======

As for all sections, this section begins with "Bonds", an empty line, a number of lines equal to the number of bonds, and an empty line (or end of file).
The format for each line (i.e., each bond) in this section is as follows:

[bond name] [bond type name] [site name i] [site name j]

Angles
======

As for all sections, this section begins with "Angles", an empty line, a number of lines equal to the number of angles, and an empty line (or end of file).
The format for each line (i.e., each angle) in this section is as follows:

[angle name] [angle type name] [site name i] [site name j] [site name k]

Note that site j is the vertex of the defined angle.
In 2D, angles are defined clockwise, such that angles ijk and kji are not the same.
For more information on the angle definition, see :cpp:func:`vertex_angle_radians <feasst::Position::vertex_angle_radians()>`

Dihedrals
==========

As for all sections, this section begins with "Dihedrals", an empty line, a number of lines equal to the number of dihedrals, and an empty line (or end of file).
The format for each line (i.e., each dihedral) in this section is as follows:

[dihedral name] [dihedral type name] [site name i] [site name j] [site name k] [site name l]

For more information on the dihedral definition, see :cpp:func:`torsion_angle_radians <feasst::Position::torsion_angle_radians()>`.

Comments
==========

For comments at the beginning of the file, begin each comment line with the "#" character.
Otherwise, commends can only be added between sections.
Use caution when adding comments anywhere else in the file.
Comments can only be one contiguous group of lines beginning with "#."
Use caution if you add comments anywhere else in the file.
End comments with a blank line immediately afterward.

Dimensions
================

A two dimensional particle can be initialized by providing the following line: "2 dimensions"

Differences from LAMMPS
========================

LAMMPS data files are defined here: https://docs.lammps.org/read_data.html

The data files used by FEASST have the following major differences:

- atoms in LAMMPS are analogous to sites in FEASST.
- molecules in LAMMPS are analogous to particles in FEASST.
- the numbers of sites, bonds, etc are not given explicitly but instead determined by the number of entries in the corresponding sections described below.
- FEASST data files contain only one particle.
- characters are case sensitive.
- the number of spaces between characters does not matter.
- LMP Coeffs sections were replaced by Properties sections, with very different formatting as described below.
- "2 dimensions" may be specified in the file for a 2D particle.
- FEASST data files contain no information about the :cpp:class:`Domain <feasst::Domain>` boundaries.
- The "Sites" section has the following three major differences: (1) these files describe only one particle so, unlike LAMMPS, molecule index is not included, (2) charge is a site-type property, and thus is not given for each site and (3) wrapping is not included.

