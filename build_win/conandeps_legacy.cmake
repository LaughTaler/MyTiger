message(STATUS "Conan: Using CMakeDeps conandeps_legacy.cmake aggregator via include()")
message(STATUS "Conan: It is recommended to use explicit find_package() per dependency instead")

find_package(cgns)
find_package(opencascade)
find_package(libtiger)

set(CONANDEPS_LEGACY  cgns::cgns  opencascade::opencascade  libtiger::libtiger )