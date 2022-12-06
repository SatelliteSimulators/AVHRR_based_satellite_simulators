# CLARA
# satellite product simulator for the CLARA-A2 and CLARA-A3 CDRs
#
# 16/12-2019

Contact:
Salomon.Eliasson@smhi.se

ABOUT:

This code simulates the CLARA-A2 or CLARA-A3 CDR cloud products
(https://wui.cmsaf.eu (last access: 30 July 2019)) based on a model atmosphere.
The simulator is described in (https://doi.org/10.5194/gmd-2019-174)

NOTE: 

*) The code is written in fortran 90

*) The simulator currently expects input from model netcdf files
(daily or monthly files have been tested) and writes daily
output files

*) The CLARA product simulator is aimed to be included in the COSP
in the near future. See https://doi.org/10.1175/2011BAMS2856.1 and
https://doi.org/10.5194/gmd-2017-148. Several subroutines from the
original COSP repository are used for compatibility.

*) The executable expects a namelist as input. An example
namelist containing mostly containing paths, dates and other
options can be used.

*) So far the simulator was tested on input from the EC Earth Earth system model.

*) Any bug discoveries or suggestions can be sent to salomon.eliasson@smhi.se

REQUIREMENTS:

- Needs access to netcdf libraries
- Need access to a netcdf land sea mask. For example, make a symlink in the top level data/ directory and call it land_sea_mask.nc

INSTRUCTIONS:


*) The code can be compiled using ifort or gfortan using 'make clara'

*) Run the simulator by '.clara.x namelist_clara

*) See simulator paper for details on input and output

