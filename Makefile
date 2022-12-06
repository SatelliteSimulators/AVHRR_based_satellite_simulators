# ******************************************
#       MAKEFILE
# Salomon.Eliasson@smhi.se
# ******************************************

# some paths are also given in here
include make.macro

##############################
# COSP routines
##############################

p_cosp_kinds = src/tools/from_COSP/cosp_kinds.F90
p_from_COSP2 = src/tools/from_COSP/from_COSP2.F90
p_cosp_constants = src/tools/from_COSP/cosp_constants.F90
p_mo_rng = src/tools/from_COSP/mo_rng.F90
p_cosp_errorHandling = src/tools/from_COSP/cosp_errorHandling.F90
p_scops = src/tools/from_COSP/scops.F90
p_cosp_config = src/tools/from_COSP/cosp_config.F90
p_cosp_stats = src/tools/from_COSP/cosp_stats.F90
p_icarus= src/tools/from_COSP/icarus.F90
p_cosp_isccp_interface= src/tools/from_COSP/cosp_isccp_interface.F90

##############################
# SHARED ROUTINES
##############################

# Types
p_optics = src/types/optics.F90
p_namelist_input = src/types/namelist_input.F90
p_model_input = src/types/model_input.F90
p_internal_simulator = src/types/internal_simulator.F90
p_satellite_specs = src/types/satellite_specs.F90
p_simulator_input_variables = src/types/simulator_input_variables.F90
p_simulator_variables = src/types/simulator_variables.F90

# Tools
p_auxiliary_functions = src/shared/auxiliary_functions.F90
p_data_check = src/shared/data_check.F90
p_handy = src/tools/handy.F90
p_my_maths = src/tools/my_maths.F90
p_my_netcdfTools = src/tools/my_netcdfTools.F90

# Modules
p_simulator_netcdf_module = src/shared/simulator_netcdf_module.F90
p_make_level2B = src/shared/make_level2B.F90
p_subcolumns = src/shared/subcolumns.F90
p_simple_simulator = src/shared/simple_simulator.F90
p_calc_from_model = src/shared/calc_from_model.F90
p_simulate_cloud_microphys = src/shared/simulate_cloud_microphys.F90

##############################
# CLOUD_CCI
##############################
p_Bcucof = src/cloud_cci/from_CC4CL/Bcucof.F90
p_Bcuint = src/cloud_cci/from_CC4CL/Bcuint.F90
p_from_CC4CL = src/cloud_cci/from_CC4CL/from_CC4CL.F90
p_cloud_cci_functions = src/cloud_cci/cloud_cci_functions.F90
p_cloud_cci_netcdf = src/cloud_cci/cloud_cci_netcdf.F90
p_cloud_cci_type = src/cloud_cci/cloud_cci_type.F90
p_cloud_cci = src/cloud_cci/cloud_cci.F90

CLOUD_CCI=\
	obj/cosp_kinds.o\
	obj/handy.o\
	obj/my_maths.o\
	obj/optics.o\
	obj/namelist_input.o\
	obj/my_netcdfTools.o\
	obj/satellite_specs.o\
	obj/model_input.o\
	obj/simulator_input_variables.o\
	obj/simulator_variables.o\
	obj/internal_simulator.o\
	obj/cloud_cci_type.o\
	obj/from_COSP2.o\
	obj/data_check.o\
	obj/from_CC4CL.o\
	obj/simple_simulator.o\
	obj/make_level2B.o\
	obj/cosp_constants.o\
	obj/calc_from_model.o\
	obj/cosp_config.o\
	obj/cosp_stats.o\
	obj/simulate_cloud_microphys.o\
	obj/auxiliary_functions.o\
	obj/mo_rng.o\
	obj/cosp_errorHandling.o\
	obj/scops.o\
	obj/subcolumns.o\
	obj/Bcucof.o\
	obj/Bcuint.o\
	obj/cloud_cci_functions.o\
	obj/cloud_cci_netcdf.o\
	obj/cloud_cci.o\

##############################
# CLARA
##############################

p_brightness_temperatures = src/clara/brightness_temperatures.F90
p_clara_functions = src/clara/clara_functions.F90
p_clara_netcdf = src/clara/clara_netcdf.F90
p_clara_type = src/clara/clara_type.F90
p_clara = src/clara/clara.F90

CLARA=\
	obj/cosp_kinds.o\
	obj/handy.o\
	obj/my_maths.o\
	obj/optics.o\
	obj/namelist_input.o\
	obj/internal_simulator.o\
	obj/my_netcdfTools.o\
	obj/model_input.o\
	obj/satellite_specs.o\
	obj/simulator_variables.o\
	obj/simulator_input_variables.o\
	obj/data_check.o\
	obj/cosp_constants.o\
	obj/auxiliary_functions.o\
	obj/clara_type.o\
	obj/cosp_config.o\
	obj/cosp_stats.o\
	obj/icarus.o\
	obj/cosp_isccp_interface.o\
	obj/brightness_temperatures.o\
	obj/from_COSP2.o\
	obj/simulate_cloud_microphys.o\
	obj/mo_rng.o\
	obj/cosp_errorHandling.o\
	obj/scops.o\
	obj/subcolumns.o\
	obj/calc_from_model.o\
	obj/make_level2B.o\
	obj/simple_simulator.o\
	obj/clara_functions.o\
	obj/clara_netcdf.o\
	obj/clara.o\


.PHONY: cloud_cci.x clara.x

clara:  $(CLARA)
	$(COMPILE) $(CLARA) -o clara.x $(NETCDF_LIB)
	@ echo " --- Finished compiling clara!"
cloud_cci:$(CLOUD_CCI)
	$(COMPILE) $(CLOUD_CCI) -o cloud_cci.x $(NETCDF_LIB)
	@ echo " --- Finished compiling cloud_cci!"

# COSP
obj/cosp_kinds.o: $(p_cosp_kinds) 
	$(COMPILE) -c -o obj/cosp_kinds.o $(p_cosp_kinds)
obj/from_COSP2.o: $(p_from_COSP2) 
	$(COMPILE) -c -o obj/from_COSP2.o $(p_from_COSP2)
obj/cosp_constants.o: $(p_cosp_constants) 
	$(COMPILE) -c -o obj/cosp_constants.o $(p_cosp_constants)
obj/mo_rng.o: $(p_mo_rng) 
	$(COMPILE) -c -o obj/mo_rng.o $(p_mo_rng)
obj/cosp_errorHandling.o: $(p_cosp_errorHandling) 
	$(COMPILE) -c -o obj/cosp_errorHandling.o $(p_cosp_errorHandling)
obj/scops.o: $(p_scops) 
	$(COMPILE) -c -o obj/scops.o $(p_scops)
obj/cosp_config.o: $(p_cosp_config) 
	$(COMPILE) -c -o obj/cosp_config.o $(p_cosp_config)
obj/cosp_stats.o: $(p_cosp_stats) 
	$(COMPILE) -c -o obj/cosp_stats.o $(p_cosp_stats)
obj/icarus.o: $(p_icarus) 
	$(COMPILE) -c -o obj/icarus.o $(p_icarus)
obj/cosp_isccp_interface.o: $(p_cosp_isccp_interface) 
	$(COMPILE) -c -o obj/cosp_isccp_interface.o $(p_cosp_isccp_interface)

#TYPES
obj/namelist_input.o: $(p_namelist_input)
	$(COMPILE) -c -o obj/namelist_input.o $(p_namelist_input)
obj/model_input.o: $(p_model_input)
	$(COMPILE) $(NETCDF_INC) -c -o obj/model_input.o $(p_model_input)
obj/internal_simulator.o: $(p_internal_simulator) 
	$(COMPILE) -c -o obj/internal_simulator.o $(p_internal_simulator)
obj/satellite_specs.o: $(p_satellite_specs)
	$(COMPILE) $(NETCDF_INC) -c -o obj/satellite_specs.o $(p_satellite_specs)
obj/simulator_input_variables.o: $(p_simulator_input_variables)
	$(COMPILE) -c -o obj/simulator_input_variables.o $(p_simulator_input_variables)
obj/simulator_variables.o: $(p_simulator_variables)
	$(COMPILE) -c -o obj/simulator_variables.o $(p_simulator_variables)
obj/optics.o: $(p_optics) 
	$(COMPILE) -c -o obj/optics.o $(p_optics)

# SHARED TOOLS
obj/handy.o: $(p_handy)
	$(COMPILE) -c -o obj/handy.o $(p_handy)
obj/my_maths.o: $(p_my_maths)
	$(COMPILE) -c -o obj/my_maths.o $(p_my_maths)
obj/my_netcdfTools.o: $(p_my_netcdfTools)
	$(COMPILE) $(NETCDF_INC) -c -o obj/my_netcdfTools.o $(p_my_netcdfTools)

# SHARED
obj/simulator_netcdf_module.o: $(p_simulator_netcdf_module)
	$(COMPILE) $(NETCDF_INC) -c -o obj/simulator_netcdf_module.o $(p_simulator_netcdf_module)
obj/make_level2B.o: $(p_make_level2B)
	$(COMPILE) -c -o obj/make_level2B.o $(p_make_level2B)
obj/data_check.o: $(p_data_check)
	$(COMPILE) -c -o obj/data_check.o $(p_data_check)
obj/auxiliary_functions.o: $(p_auxiliary_functions)
	$(COMPILE) $(NETCDF_INC) -c -o obj/auxiliary_functions.o $(p_auxiliary_functions)
obj/subcolumns.o: $(p_subcolumns)
	$(COMPILE) -c -o obj/subcolumns.o $(p_subcolumns)
obj/calc_from_model.o: $(p_calc_from_model)
	$(COMPILE) -c -o obj/calc_from_model.o $(p_calc_from_model)
obj/simple_simulator.o: $(p_simple_simulator)
	$(COMPILE) -c -o obj/simple_simulator.o $(p_simple_simulator)
obj/simulate_cloud_microphys.o: $(p_simulate_cloud_microphys)
	$(COMPILE) -c -o obj/simulate_cloud_microphys.o $(p_simulate_cloud_microphys)

# CLARA
obj/clara_type.o: $(p_clara_type)
	$(COMPILE) -c -o obj/clara_type.o $(p_clara_type)
obj/brightness_temperatures.o: $(p_brightness_temperatures)
	$(COMPILE) -c -o obj/brightness_temperatures.o $(p_brightness_temperatures)
obj/clara_functions.o: $(p_clara_functions)
	$(COMPILE) -c -o obj/clara_functions.o $(p_clara_functions)
obj/clara_netcdf.o: $(p_clara_netcdf)
	$(COMPILE) $(NETCDF_INC) -c -o obj/clara_netcdf.o $(p_clara_netcdf)
obj/clara.o: $(p_clara)
	$(COMPILE) -c -o obj/clara.o $(p_clara)

# CLOUD_CCI
obj/cloud_cci_type.o: $(p_cloud_cci_type) 
	$(COMPILE) -c -o obj/cloud_cci_type.o $(p_cloud_cci_type)
obj/Bcucof.o: $(p_Bcucof)
	$(COMPILE) -c -o obj/Bcucof.o $(p_Bcucof)
obj/Bcuint.o: $(p_Bcuint)
	$(COMPILE) -c -o obj/Bcuint.o $(p_Bcuint)
obj/from_CC4CL.o: $(p_from_CC4CL)
	$(COMPILE) -c -o obj/from_CC4CL.o $(p_from_CC4CL)
obj/cloud_cci_functions.o: $(p_cloud_cci_functions)
	$(COMPILE) -c -o obj/cloud_cci_functions.o $(p_cloud_cci_functions)
obj/cloud_cci_netcdf.o: $(p_cloud_cci_netcdf)
	$(COMPILE) $(NETCDF_INC) -c -o obj/cloud_cci_netcdf.o $(p_cloud_cci_netcdf)
obj/cloud_cci.o: $(p_cloud_cci)
	$(COMPILE) -c -o obj/cloud_cci.o $(p_cloud_cci)

clean:
	rm -rf mod/* obj/* *.x *~ *\#
#
all: clara cloud_cci

install:
	cp -f *.x ${INSTALL_DIR}/
	mv -f *.x bin/
	mv -f *.mod mod/
