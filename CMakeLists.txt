project(liegroupfem)

cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
  HINTS /usr $ENV{NETGENDIR}/.. /tmp/inst/usr /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)


add_ngsolve_python_module(liegroupfem src/liegroupfe.hpp src/liegroupfem.cpp src/quaternion.hpp src/liegroups.hpp src/SO3coefficient.hpp src/liegroupfem.cpp)




# check if CMAKE_INSTALL_PREFIX is set by user, otherwise install to netgen directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${NETGEN_DIR} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib(1,0,''))" OUTPUT_VARIABLE install_dir OUTPUT_STRIP_TRAILING_WHITESPACE)


message("The python package will be installed to: ${CMAKE_INSTALL_PREFIX}/${install_dir}")
install(TARGETS liegroupfem DESTINATION ${install_dir})
