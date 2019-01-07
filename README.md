# Cloud_cci
# satellite product simulator for the Cloud_cci CDR
#
#
#
# 
# Salomon.Eliasson@smhi.se

ABOUT:

This code is the simulates the Cloud_cci CDR cloud products
(http://www.esa-cloud-cci.org/) based on a model atmosphere.
The simulator is described in (https://doi.org/10.5194/gmd-2018-212)

NOTE: 

      - The code is written in fortran 90

      - The simulator currently expects input from model netcdf files
        (daily or monthly files have been tested) and writes daily
        output files

      - The Cloud_cci product simulator is aimed to be included in the
        COSP in the future. See https://doi.org/10.1175/2011BAMS2856.1
        and https://doi.org/10.5194/gmd-2017-148. Several low level
        codes are copied from original COSP repository. 

      - The executable expects a namelist as input. An example
        namelist containing mostly containing paths, dates and other
        options can be used.

      - So far the simulator has been used on EC Earth and RACMO model input.
      
      - Any bug discoveries or suggestions can be sent to salomon.eliasson@smhi.se


INSTRUCTIONS:
	- The code can be compiled using ifort or gfortan using 'make cloud_cci'
	- Run the simulator by '.cloud_cci.x namelist_cloud_cci


INPUT: 

'CC'   = Horizontal fraction of the grid box covered by cloud [0-1]. code=248
'CI'   = Sea Ice Cover [0-1]. code=31
'CIWC' = Specific cloud ice water content. Grid-box
         mean specific cloud ice water content (mass of condensate
         / mass of moist air).  [kg kg**-1]. code=247
'CLWC' = Specific cloud liquid water content. Grid-box mean
         specific cloud liquid water content (mass of condensate /
         mass of moist air) [kg kg**-1]. code=246
'PSURF'= Surface pressure [(Pa)]. code=1
'Q'    = Specific humidity. Grid box mean (mass of water vapour /
         mass of moist air). [kg kg**-1]. code=133
'SKT'  = Skin temperature. Temperature of the surface skin
         (radiative surface temperature). Before 01/10/2008, the skin
         temperature was equal to the bulk SST over the
         ocean. [K]. code=235
'T'    = Temperature [K]. code=130
'T2M'  = 2 metre temperature [K]. code=167
'TCC'  = Total cloud cover. Total cloud cover derived from model
         levels using the model's overlap assumption [0-1]. code=164
'TCWV' = Total column water vapour (Vertically integrated water
         vapour) [kg m**-2]. code=137

OUTPUT:

       cer_ice = simulated ice effective radius
       cer_liq = simulated liquid effective radiusd
           cfc = simulated total cloud cover
       cfc_low = simulated total cloud cover forclouds with CTP greater than 680 hPa
       cfc_mid = simulated total cloud cover forclouds with CTP less than 680 hPa and greater than 440 hPa
      cfc_high = simulated total cloud cover forclouds with CTP less than 440 hPa
    cla_vis006 = Simulated cloud albedo at 0.6 micron
           cot = grid average cloud optical thickness
           cth = simulated cloud top height
           ctp = simulated cloud top pressure. Derived from linear-averaging of the sub-grid CTP
       ctp_log = simulated cloud top pressure. Derived from log-averaging of the sub-grid CTP
hist2d_cot_ctp = CTP--tau hits (lon,lat,n*tau bins,n*pressure bins)  based on   70sub-grids in each grid box and only including clouds that have tau >  0.20
           iwp = simulated grid average ice water path
           lwp = grid average liquid water path


