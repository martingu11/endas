
add_executable(SimpleExample Utils.cpp SimpleExample.cpp)
add_executable(Lorenz95Example Utils.cpp Lorenz95Example.cpp)

#target_link_libraries(endas PUBLIC Eigen3::Eigen)

target_include_directories(SimpleExample PRIVATE ${THIRDPARTY_DIR}/matplotlibcpp )
target_include_directories(Lorenz95Example PRIVATE ${THIRDPARTY_DIR}/matplotlibcpp )


if (PYTHON2_FOUND)
    target_link_libraries(SimpleExample PRIVATE endas Python2::Python)
    target_link_libraries(Lorenz95Example PRIVATE endas Python2::Python)

    target_include_directories(SimpleExample PRIVATE ${Python2_INCLUDE_DIRS} )
    target_include_directories(Lorenz95Example PRIVATE ${Python2_INCLUDE_DIRS} )
endif()

install(TARGETS SimpleExample DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS Lorenz95Example DESTINATION ${CMAKE_INSTALL_BINDIR})


