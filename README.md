# Cloud_cci
# satellite product simulator for the Cloud_cci CDR

Contact:
Salomon.Eliasson@smhi.se

ABOUT:

This code is the simulates the Cloud_cci CDR cloud products
(http://www.esa-cloud-cci.org/) based on a model atmosphere.
The simulator is described in (https://doi.org/10.5194/gmd-2018-212)

NOTE: 

*) The code is written in fortran 90

*) The simulator currently expects input from model netcdf files
(daily or monthly files have been tested) and writes daily
output files

*) The Cloud_cci product simulator is aimed to be included in the COSP
in the future. See https://doi.org/10.1175/2011BAMS2856.1 and
https://doi.org/10.5194/gmd-2017-148. Several subroutines from the
original COSP repository are used for compatibility.

*) The executable expects a namelist as input. An example
namelist containing mostly containing paths, dates and other
options can be used.

*) So far the simulator has been used on EC Earth and RACMO model input.

*) Any bug discoveries or suggestions can be sent to salomon.eliasson@smhi.se


INSTRUCTIONS:

*) The code can be compiled using ifort or gfortan using 'make cloud_cci'
*) Run the simulator by '.cloud_cci.x namelist_cloud_cci
*) See simulator paper for details on input and output

