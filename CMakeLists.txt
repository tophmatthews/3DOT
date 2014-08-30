cmake_minimum_required(VERSION 2.6)
project(mytrim)

# build the mytrim library
# put install prefix somehow into config file
set( DATA_INSTALL ${CMAKE_SOURCE_DIR}/data )
set( DATA_DIR ${DATA_INSTALL} )
configure_file( ${CMAKE_SOURCE_DIR}/config.h.cmake ${CMAKE_SOURCE_DIR}/config.h )
add_library(mytrim invert.cc element.cc material.cc sample.cc ion.cc trim.cc simconf.cc sample_clusters.cc bubble.cc sample_single.cc)


include_directories( "${PROJECT_SOURCE_DIR}/alglib")
add_subdirectory(alglib)

install(TARGETS mytrim ARCHIVE DESTINATION lib)
install(FILES simconf.h element.h material.h ion.h trim.h invert.h sample.h bubble.h sample_single.cc DESTINATION include/mytrim)
install(FILES algplibinternal.h alglibmisc.h ap.h dataanalysis.h diffequations.h fasttransforms.h integration.h interpolation.h linalg.h optimization.h solvers.h specialfunctions.h statistics.h stdafx.h DESTINATION include/mytrim)
install(FILES sample_clusters.h DESTINATION include/mytrim)
install(FILES data/SCOEF.95A data/SLFCTR.dat data/SNUC03.dat data/ELNAME.dat DESTINATION ${DATA_DIR})

# build the trim programs
set(PROGRAMS mytrim_bub mytrim_range )
foreach(program ${PROGRAMS})
  add_executable(${program} ${program}.cc)
  target_link_libraries(${program}  mytrim m)
  target_link_libraries(${program}  alglib)
endforeach(program)
  
# might use boost sometime
#find_package(Boost REQUIRED program_options)
#include_directories(${Boost_INCLUDE_DIRS})
#target_link_libraries(firstpass ${Boost_LIBRARIES} ${R250C_LIB} m)

