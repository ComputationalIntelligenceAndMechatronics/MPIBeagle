cmake_minimum_required(VERSION 2.5)
cmake_policy(SET CMP0003 NEW)

set( CMAKE_VERBOSE_MAKEFILE ON )
 
project (MPIBEAGLE)
set(MPIBEAGLE_VERSION 0.0.1)


# include( FindMPI )
# if( NOT MPI_FOUND )
# 	if( CMAKE_SYSTEM_NAME STREQUAL SunOS )
# 		set( CMAKE_CXX_COMPILER /opt/SUNWhpc/bin/mpCC )
# 		set( MPI_INCLUDE_PATH /opt/SUNWhpc/include )
# 		set( MPI_LINK_FLAGS /opt/SUNWhpc/lib )
# 		set( MPI_LIB_DIR /opt/SUNWhpc/HPC5.0/lib )
# 		set( MPI_LIBRARIES mpi )
# 	else( CMAKE_SYSTEM_NAME STREQUAL SunOS )
		set( MPI_DIR /usr/mpi/gcc/openmpi-1.3.2 CACHE FILEPATH "MPI directory" )
		set( MPI_COMPILER ${MPI_DIR}/bin/mpic++ CACHE FILEPATH "MPI compiler" )
		set( CC ${MPI_COMPILER} )
		set( MPI_INCLUDE_PATH ${MPI_DIR}/include )
		set( MPI_LINK_FLAGS ${MPI_DIR}/lib64 )
		set( MPI_LIB_DIR ${MPI_DIR}/lib64 )
		set( MPI_LIBRARIES mpi mpi_cxx)
#	endif( CMAKE_SYSTEM_NAME STREQUAL SunOS )
#endif( NOT MPI_FOUND )


# Check for PACC
if(NOT PACC_DIR)
    set(PACC_DIR "$ENV{HOME}/local")
    message("++ PACC_DIR set to default ~/local")
endif(NOT PACC_DIR)
if(NOT PACC_INCLUDE_DIR)
    set(PACC_INCLUDE_DIR "${PACC_DIR}/include/PACC")
    message("++ PACC_INCLUDE_DIR set to default ${PACC_DIR}/include/PACC")
endif(NOT PACC_INCLUDE_DIR)
if(NOT PACC_LIB_DIR)
    set(PACC_LIB_DIR "${PACC_DIR}/lib")
    message("++ PACC_LIB_DIR set to default ${PACC_DIR}/lib")
endif(NOT PACC_LIB_DIR)

include_directories( $ENV{HOME}/local/include /sw/include ${PACC_INCLUDE_DIR} ${PACC_INCLUDE_DIR}/PACC ${MPI_INCLUDE_PATH} )
#find_library( $ENV{HOME}/local/lib /sw/lib  ${MPI_LIB_DIR})
link_directories($ENV{HOME}/local/lib /sw/lib  ${MPI_LIB_DIR} ${PACC_LIB_DIR}) 

if( CMAKE_SYSTEM_NAME STREQUAL SunOS )
	set( CMAKE_CXX_FLAGS "-DSunMPI")
endif( CMAKE_SYSTEM_NAME STREQUAL SunOS )

set( MPIBEAGLE_HEADERS 
	Source/CommunicationMPI.h
	Source/MPI_EvaluationOp.hpp
	Source/MPI_Evolver.hpp
	Source/MPI_GA_EvolverBitString.hpp
	Source/MPI_GA_EvolverFloatVector.hpp
	Source/MPI_GP_EvaluationOp.hpp
	Source/MPI_GP_Evolver.hpp
	Source/MPI_Coev_EvaluationOp.hpp
	Source/MPI_Coev_FitnessEvaluationClient.hpp
	Source/VectorUtil.h
)

set( MPIBEAGLE_SRCS 
	Source/CommunicationMPI.cpp
	Source/MPI_EvaluationOp.cpp
	Source/MPI_Evolver.cpp
	Source/MPI_GA_EvolverBitString.cpp
	Source/MPI_GA_EvolverFloatVector.cpp
	Source/MPI_GP_EvaluationOp.cpp
	Source/MPI_GP_Evolver.cpp
	Source/MPI_Coev_EvaluationOp.cpp
	Source/MPI_Coev_FitnessEvaluationClient.cpp
	Source/VectorUtil.cpp
)

add_library (openbeagle-MPI SHARED ${MPIBEAGLE_SRCS})
target_link_libraries(openbeagle-MPI openbeagle openbeagle-GP openbeagle-GA pacc z ${MPI_LIBRARIES})
set_target_properties(openbeagle-MPI PROPERTIES VERSION ${MPIBEAGLE_VERSION})
set_target_properties(openbeagle-MPI PROPERTIES LINKER_LANGUAGE CXX)

# Install dependencies
install(FILES ${MPIBEAGLE_HEADERS} DESTINATION include/beagle)

# Install targets
install(TARGETS openbeagle-MPI DESTINATION lib)

# Test example
set( MAXFCT_SRCS 
	Source/MaxFctEvalOp.cpp
	Source/MaxFctMain.cpp
)

set( MAXFCTLIBS openbeagle-MPI openbeagle openbeagle-GA pacc z ssl pthread ${MPI_LIBRARIES})
add_executable (MaxFct ${MAXFCT_SRCS})
target_link_libraries(MaxFct ${MAXFCTLIBS} )
