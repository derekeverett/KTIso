set (SOURCES_LIB
    ITA.cpp    
    )

set (SOURCES
    RunWrapper.cpp
    )

add_library(ita_lib SHARED ${SOURCES_LIB})
add_executable(RunWrapper ${SOURCES})
target_link_libraries (ita_lib ${LIBS})
target_link_libraries (RunWrapper ${LIBS})

install(TARGETS RunWrapper DESTINATION ${CMAKE_HOME_DIRECTORY})
