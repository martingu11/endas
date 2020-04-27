add_library(endas SHARED "")

# EnDAS source files
target_sources(endas PRIVATE
    Endas.cpp
    Core/LinAlg.cpp
    Core/Exception.cpp
    Core/Profiling.cpp

    DA/Ensemble.cpp
    DA/Model.cpp
    DA/Sequential.cpp
    DA/ObservationOperator.cpp
    DA/SimpleObservationManager.cpp
    DA/CovarianceOperator.cpp
    DA/StateSpace.cpp
    DA/StateSpaceGrid.cpp
    DA/StateSpaceGeneric.cpp
    DA/Taper.cpp

    DA/Algorithm/KalmanSmoother.cpp
    DA/Algorithm/EnsembleKalmanSmoother.cpp
    DA/Algorithm/EnKSVariantEnKS.cpp
    DA/Algorithm/EnKSVariantESTKS.cpp
    
    Caching/ArrayCache.cpp
    Caching/MemoryArrayCache.cpp
    
    Random/Random.cpp
    Random/MultivariateRandomNormal.cpp
    Random/GaussianRandomField.cpp

    Spatial/CoordinateSystem.cpp
    Spatial/Variogram.cpp

    Model/Lorenz95.cpp
)

target_include_directories(endas PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_include_directories(endas PRIVATE ${THIRDPARTY_DIR}/simplefft/include)


target_link_libraries(endas PUBLIC Eigen3::Eigen)

if (ENDAS_USING_BLAS)
    target_link_libraries(endas PUBLIC ${BLAS_LIBRARIES})
endif()

if (ENDAS_USING_LAPACK)
    #target_link_libraries(endas PUBLIC ${LAPACK_LIBRARIES})
    target_link_libraries(endas PUBLIC ${LAPACKE_LIBRARIES})    
endif()


if (ENDAS_USING_OPENMP)
    target_link_libraries(endas PRIVATE OpenMP::OpenMP_CXX)
endif()



install(TARGETS endas DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})



