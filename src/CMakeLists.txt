add_library(endas SHARED "")

# EnDAS source files
target_sources(endas PRIVATE
    Endas.cpp
    Core/LinAlg.cpp
    Core/Exception.cpp
    Core/Profiling.cpp
    Core/Ensemble.cpp

    DA/Model.cpp
    DA/Sequential.cpp
    DA/ObservationOperator.cpp
    DA/CovarianceOperator.cpp
    DA/Domain.cpp
    DA/Taper.cpp

    DA/GenericDomain.cpp
    DA/GridDomain.cpp
    DA/GridDomainPartitioning.cpp
    DA/SimpleObservationManager.cpp

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

    IO/ArrayIO_NumPy.cpp

    Utils/EnsembleSampler.cpp


    Model/Lorenz95.cpp
    Model/QG/QG.cpp
    Model/QG/calc.f90
    Model/QG/helmholtz.f90
    Model/QG/parameters.f90
    Model/QG/qgflux.f90
    Model/QG/qgstep.f90
    Model/QG/utils.f90
)

if (CMAKE_Fortran_COMPILER)
    target_sources(endas PRIVATE
        Model/QG/calc.f90
        Model/QG/helmholtz.f90
        Model/QG/parameters.f90
        Model/QG/qgflux.f90
        Model/QG/qgstep.f90
        Model/QG/utils.f90
    )
endif()


# EnDAS source files coming from thirdparty libs
target_sources(endas PRIVATE
    ${THIRDPARTY_DIR}/cnpy/cnpy.cpp
)



target_include_directories(endas PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_include_directories(endas PRIVATE ${THIRDPARTY_DIR}/simplefft/include)
target_include_directories(endas PRIVATE ${THIRDPARTY_DIR}/cnpy)
target_include_directories(endas PRIVATE ${THIRDPARTY_DIR}/libspatial)


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



