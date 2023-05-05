cmake_minimum_required(VERSION 2.9 FATAL_ERROR)
project(nnbar-calo-sim)

set(PROGRAM_PERMISSIONS_DEFAULT
	OWNER_WRITE OWNER_READ OWNER_EXECUTE
	GROUP_WRITE GROUP_READ GROUP_EXECUTE
	WORLD_WRITE WORLD_READ WORLD_EXECUTE)

option(WITH_GEANT4_UIVIS "Build nnbar-calo-sim with Geant4 UI and Vis drivers" ON)

if(WITH_GEANT4_UIVIS)
	find_package(Geant4 REQUIRED ui_all vis_all)
else()
	find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})
MESSAGE(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#

include_directories(${PROJECT_SOURCE_DIR}/include 
	${CRY_INCLUDE_DIRS}
	${CMAKE_CURRENT_BINARY_DIR}
	)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*)

add_executable(nnbar-calo-sim nnbar-calo-sim.cc ${sources} ${headers})

target_link_libraries(nnbar-calo-sim ${Geant4_LIBRARIES} ${CRY_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.

set(SCRIPTS
	gui.mac
	vis.mac
	init_vis.mac
	template.mac
	)

foreach(_script ${SCRIPTS})
	configure_file(
		${PROJECT_SOURCE_DIR}/${_script}
		${PROJECT_BINARY_DIR}/${_script}
		COPYONLY
		)
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS nnbar-calo-sim DESTINATION bin)
install(
	DIRECTORY
	DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/output"
	DIRECTORY_PERMISSIONS ${PROGRAM_PERMISSIONS_DEFAULT}
	)
