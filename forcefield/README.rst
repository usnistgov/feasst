**********************
Forcefield
**********************

The forcefield directory is a place for files which describe atoms, molecules and coarse-grained models.

===========
Data files
===========

Files which begin with the name "data" are LAMMPS-inspired data files.

https://lammps.sandia.gov/doc/read_data.html

They deviate from LAMMPS data files as follows:

Nomenclature and style
=======================

* atoms in LAMMPS are analogous to sites in FEASST.
* molecules in LAMMPS are analogous to particles in FEASST.
* FEASST data files may only contain one particle.
* all FEASST indices (types, sites, etc) begin with 0, not 1.
* characters are case sensitive
* the number of spaces between characters does not matter
* LMP Coeffs sections were replaced by Properties sections, with different formatting
* "2 dimensions" may be specified at the beginning of the file for a 2D simulation.

Site Properties
================

This is where the biggest difference from the LAMMPS data format is present.
Each site type may have a list of properties given by a label and value.

The format for this section is as follows:

[site index] [label_0] [value_0] ... [label_n] [value_n]

Sites
======

Charge is a site-type property, and thus is not included for each site.
In addition, wrapping is also not included.

The format for this section is as follows:

[site index] [site type] [x-position] [y-position] [z-position]

Site Labels
=============

Site labels are used to attach a name to the site type, often used for visualization programs, etc.
