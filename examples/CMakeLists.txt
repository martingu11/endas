
add_executable(simple Utils.cpp SimpleExample.cpp)
add_executable(lorenz95 Utils.cpp Lorenz95Example.cpp)
add_executable(qg Utils.cpp QGExample.cpp)
#dd_executable(datafusion Utils.cpp DataFusionExample.cpp)


#target_link_libraries(endas PUBLIC Eigen3::Eigen)

target_include_directories(simple PRIVATE ${THIRDPARTY_DIR}/matplotlibcpp )
target_include_directories(lorenz95 PRIVATE ${THIRDPARTY_DIR}/matplotlibcpp )
target_include_directories(qg PRIVATE ${THIRDPARTY_DIR}/matplotlibcpp )
#target_include_directories(datafusion PRIVATE ${THIRDPARTY_DIR}/matplotlibcpp )

target_link_libraries(simple PRIVATE endas)
target_link_libraries(lorenz95 PRIVATE endas)
target_link_libraries(qg PRIVATE endas)
#target_link_libraries(datafusion PRIVATE endas)


if (PYTHON2_FOUND)
    target_link_libraries(simple PRIVATE Python2::Python)
    target_link_libraries(lorenz95 PRIVATE Python2::Python)
    target_link_libraries(qg PRIVATE Python2::Python)
    #target_link_libraries(datafusion PRIVATE Python2::Python)
endif()



#if (ENDAS_USING_OPENMP)
#    target_link_libraries(datafusion PRIVATE OpenMP::OpenMP_CXX)
#endif()



install(TARGETS simple DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS lorenz95 DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS qg DESTINATION ${CMAKE_INSTALL_BINDIR})
#install(TARGETS datafusion DESTINATION ${CMAKE_INSTALL_BINDIR})


