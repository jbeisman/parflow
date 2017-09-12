add_executable(parflow main.c)
target_link_libraries(parflow pfsimulator pfkinsol amps)

if( ${PARFLOW_HAVE_CLM} )
  target_link_libraries(parflow pfclm)
endif( ${PARFLOW_HAVE_CLM} )

if (${PARFLOW_HAVE_HYPRE})
  target_link_libraries (parflow ${HYPRE_LIBRARIES})
endif (${PARFLOW_HAVE_HYPRE})

if (${PARFLOW_HAVE_MPI})
  target_link_libraries (parflow ${MPI_LIBRARIES})
endif (${PARFLOW_HAVE_MPI})

if (${PARFLOW_HAVE_SILO})
  target_link_libraries (parflow ${SILO_LIBRARIES})
endif (${PARFLOW_HAVE_SILO})

if (${PARFLOW_HAVE_HDF5})
  target_link_libraries (parflow ${HDF5_LIBRARIES})
endif (${PARFLOW_HAVE_HDF5})

if (${PARFLOW_HAVE_ZLIB})
  target_link_libraries (parflow ${ZLIB_LIBRARIES})
endif (${PARFLOW_HAVE_ZLIB})

if (${PARFLOW_HAVE_SZLIB})
  target_link_libraries (parflow ${SZLIB_LIBRARIES})
endif (${PARFLOW_HAVE_SZLIB})

if (${PARFLOW_HAVE_SLURM})
  target_link_libraries (parflow ${SLURM_LIBRARIES})
endif (${PARFLOW_HAVE_SLURM})

if( ${PARFLOW_ENABLE_PROFILING} )
  set_target_properties(parflow PROPERTIES LINK_FLAGS ${PARFLOW_PROFILE_OPTS})
endif( ${PARFLOW_ENABLE_PROFILING} )

if ( DEFINED PARFLOW_LINKER_FLAGS)
   set_target_properties(parflow PROPERTIES LINK_FLAGS ${PARFLOW_LINKER_FLAGS})
endif ( DEFINED PARFLOW_LINKER_FLAGS)


  
install(TARGETS parflow DESTINATION bin)