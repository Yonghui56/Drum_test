## pthread ##
SET ( CMAKE_THREAD_PREFER_PTHREAD On )
FIND_PACKAGE( Threads )
IF ( CMAKE_USE_PTHREADS_INIT )
        SET (HAVE_PTHREADS TRUE)
        MESSAGE (STATUS "pthread library found." )
ENDIF (CMAKE_USE_PTHREADS_INIT )


# Source files
FILE(GLOB CUI_Files *.h *.cpp)
SOURCE_GROUP( ogs6CUI FILES ${CUI_Files} )


INCLUDE_DIRECTORIES(
        .
        ../external
        ../external/logog/include
        ../external/Eigen
        ../external/tinyxml2
        ../external/muparser/include
        ../ogsCore
        ../ogsNumerics
        ../ogsFileIO
        ../ogs6THMC
)


## Create the library
ADD_EXECUTABLE( ogs6 ${CUI_Files})
TARGET_LINK_LIBRARIES( ogs6 
    ogsCore
    ogsNumerics
    ogs6THMC
    logog
    ${CMAKE_THREAD_LIBS_INIT}
    )

#IF(USE_EXTERNAL_LIS)
#    INCLUDE_DIRECTORIES( ${LIS_INCLUDE_DIR} )
#    TARGET_LINK_LIBRARIES( ogs6 ${LIS_LIBRARIES} )
#ELSE()
#	IF (USE_LIS)
#		INCLUDE_DIRECTORIES( ../external/lis/include )
#		TARGET_LINK_LIBRARIES( ogs6 lis )
#	ENDIF()
#ENDIF ()

IF (USE_BLAS_LAPACK)
    TARGET_LINK_LIBRARIES( ogs6 ${BLAS_LIBRARIES} )
    TARGET_LINK_LIBRARIES( ogs6 ${LAPACK_LIBRARIES} )
ENDIF (USE_BLAS_LAPACK)

IF(USE_PETSC)
        TARGET_LINK_LIBRARIES( ogs6  ${PETSC_LIBRARIES} )
    IF (MSVC)
        SET_TARGET_PROPERTIES(ogs6 PROPERTIES
                LINK_FLAGS "
            /NODEFAULTLIB:\"libcmt.lib\"
                ")
    ENDIF(MSVC)
ENDIF(USE_PETSC)

IF (USE_MPI)
    INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
        TARGET_LINK_LIBRARIES( ogs6 ${MPI_C_LIBRARIES} )
ENDIF (USE_MPI)

IF (WIN32)
    TARGET_LINK_LIBRARIES( ogs6 Winmm.lib )
ENDIF()

SET_TARGET_PROPERTIES(ogs6 PROPERTIES LINKER_LANGUAGE CXX)

