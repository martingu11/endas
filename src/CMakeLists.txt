add_library(endas SHARED "")

# EnDAS source files
target_sources(endas PRIVATE
    Endas.cpp
    Core/LinAlg.cpp
    Core/Exception.cpp
    Algorithm/Algorithm.cpp
    Algorithm/KalmanSmoother.cpp
    Algorithm/EnsembleKalmanSmoother.cpp
    Caching/ArrayCache.cpp
    Caching/MemoryArrayCache.cpp
    Error/CovarianceOperator.cpp
    Error/DiagonalCovariance.cpp
    Model/Trivial.cpp
    Model/Lorenz95.cpp
    Random/Random.cpp
    Random/MultivariateRandomNormal.cpp
    
)


# EnDAS header files to be included in IDEs
target_sources(endas PRIVATE
    ${CMAKE_SOURCE_DIR}/include/Endas/Endas.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Core/LinAlg.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Core/Exception.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Algorithm/Algorithm.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Algorithm/KalmanSmoother.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Algorithm/EnsembleKalmanSmoother.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Caching/ArrayCache.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Caching/MemoryArrayCache.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Error/CovarianceOperator.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Error/DiagonalCovariance.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Model/Model.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Model/Lorenz95.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Model/Trivial.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Random/Random.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Random/MultivariateRandomNormal.hpp
)


target_include_directories(endas PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(endas PUBLIC Eigen3::Eigen)

if (ENDAS_USING_BLAS)
    target_link_libraries(endas PUBLIC ${BLAS_LIBRARIES})
endif()

if (ENDAS_USING_LAPACK)
    #target_link_libraries(endas PUBLIC ${LAPACK_LIBRARIES})
    target_link_libraries(endas PUBLIC ${LAPACKE_LIBRARIES})    
endif()



install(TARGETS endas DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})



