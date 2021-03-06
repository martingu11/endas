# EnDAS header files 
# Not needed for building but we do want headers to be included in IDE projects

target_sources(endas PRIVATE
    ${CMAKE_SOURCE_DIR}/include/Endas/Endas.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Core/LinAlg.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Core/AABox.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Core/Exception.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Core/Profiling.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Core/Ensemble.hpp

    ${CMAKE_SOURCE_DIR}/include/Endas/DA/Model.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/DA/Sequential.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/DA/ObservationOperator.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/DA/ObservationManager.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/DA/CovarianceOperator.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/DA/Domain.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/DA/DomainPartitioning.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/DA/Taper.hpp

    ${CMAKE_SOURCE_DIR}/include/Endas/DA/GenericDomain.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/DA/SimpleObservationManager.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/DA/GridDomain.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/DA/GridDomainPartitioning.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/DA/IndexedPartitionPointQuery.hpp

    ${CMAKE_SOURCE_DIR}/include/Endas/DA/Algorithm/KalmanSmoother.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/DA/Algorithm/EnsembleKalmanSmoother.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/DA/Algorithm/EnsembleKalmanSmootherVariant.hpp

    ${CMAKE_SOURCE_DIR}/include/Endas/Caching/ArrayCache.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Caching/MemoryArrayCache.hpp

    ${CMAKE_SOURCE_DIR}/include/Endas/Parallel/Parallel.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Parallel/AsyncJobExecutor.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Parallel/SerialJobExecutor.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Parallel/ThreadJobExecutor.hpp

    ${CMAKE_SOURCE_DIR}/include/Endas/Random/Random.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Random/MultivariateRandomNormal.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Random/GaussianRandomField.hpp

    ${CMAKE_SOURCE_DIR}/include/Endas/Spatial/Spatial.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Spatial/Variogram.hpp
    ${CMAKE_SOURCE_DIR}/include/Endas/Spatial/CoordinateSystem.hpp

    ${CMAKE_SOURCE_DIR}/include/Endas/IO/ArrayIO.hpp

    ${CMAKE_SOURCE_DIR}/include/Endas/Utils/EnsembleSampler.hpp
)


target_sources(endas PRIVATE
    ${CMAKE_SOURCE_DIR}/include/EndasModels/Lorenz95.hpp
)

