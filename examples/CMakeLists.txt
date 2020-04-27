
add_executable(SimpleExample Utils.cpp SimpleExample.cpp)
add_executable(Lorenz95Example Utils.cpp Lorenz95Example.cpp)
add_executable(KrigingExample Utils.cpp KrigingExample.cpp)

#target_link_libraries(endas PUBLIC Eigen3::Eigen)

target_include_directories(SimpleExample PRIVATE ${THIRDPARTY_DIR}/matplotlibcpp )
target_include_directories(Lorenz95Example PRIVATE ${THIRDPARTY_DIR}/matplotlibcpp )
target_include_directories(KrigingExample PRIVATE ${THIRDPARTY_DIR}/matplotlibcpp )

target_link_libraries(SimpleExample PRIVATE endas)
target_link_libraries(Lorenz95Example PRIVATE endas)
target_link_libraries(KrigingExample PRIVATE endas)


if (PYTHON2_FOUND)
    target_link_libraries(SimpleExample PRIVATE Python2::Python)
    target_link_libraries(Lorenz95Example PRIVATE Python2::Python)
    target_link_libraries(KrigingExample PRIVATE Python2::Python)
endif()



if (ENDAS_USING_OPENMP)
    target_link_libraries(KrigingExample PRIVATE OpenMP::OpenMP_CXX)
endif()



install(TARGETS SimpleExample DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS Lorenz95Example DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS KrigingExample DESTINATION ${CMAKE_INSTALL_BINDIR})


