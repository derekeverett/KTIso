set (SOURCES_LIB
    KTiso.cpp
    )

set (SOURCES
    RunWrapper.cpp
    )

add_library(ktiso_lib SHARED ${SOURCES_LIB})
add_executable(RunWrapper ${SOURCES})
target_link_libraries (ktiso_lib ${LIBS})
target_link_libraries (RunWrapper ${LIBS})

install(TARGETS RunWrapper DESTINATION ${CMAKE_HOME_DIRECTORY})
