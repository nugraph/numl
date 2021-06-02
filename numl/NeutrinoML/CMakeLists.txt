include_directories( $ENV{HEP_HPC_INC} )
include_directories( $ENV{HDF5_INC} )
cet_find_library( HEP_HPC_HDF5_LIB NAMES hep_hpc_hdf5 PATHS ENV HEP_HPC_LIB NO_DEFAULT_PATH )
cet_find_library( HDF5_LIB NAMES hdf5 PATHS ENV HDF5_LIB NO_DEFAULT_PATH )

#art_make(
#  LIB_LIBRARIES
simple_plugin(HDF5Maker "module"
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_BASIC}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ROOT_CORE}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_PERSISTENCY_COMMON}
    canvas
    ${ART_UTILITIES}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    cetlib
    cetlib_except
    ${ROOT_BASIC_LIB_LIST}
    nusimdata_SimulationBase
    larsim_MCCheater_BackTrackerService_service           
    larsim_MCCheater_ParticleInventoryService_service
    ${HEP_HPC_HDF5_LIB}
    ${HDF5_LIB}
)

simple_plugin(PandoraNuSliceHitsProducer "module"
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_BASIC}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_ROOT_IO_TFILE_SUPPORT}
    ${ROOT_CORE}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_PERSISTENCY_COMMON}
    canvas
    ${ART_UTILITIES}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    cetlib
    cetlib_except
    ${ROOT_BASIC_LIB_LIST}
    nusimdata_SimulationBase
    lardataobj_RecoBase
)

install_headers()
install_fhicl()
install_source()
