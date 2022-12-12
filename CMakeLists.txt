cmake_minimum_required (VERSION 3.12)

project(numl LANGUAGES CXX)

find_package(cetbuildtools REQUIRED)


list(APPEND CMAKE_MODULE_PATH
     $ENV{CANVAS_ROOT_IO_DIR}/Modules
     $ENV{ART_DIR}/Modules
)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

find_ups_product( lardata )
find_ups_product( lardataalg )
find_ups_product( lardataobj )
find_ups_product( larsim )
find_ups_product( nutools )
find_ups_product( art )
find_ups_product( cetbuildtools )
find_ups_product( hep_hpc )
find_ups_product( hdf5 )

include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(numl)
add_subdirectory(test)
add_subdirectory(ups)

include(UseCPack)
