# Distributed under the MIT License.
# See LICENSE.txt for details.

set(EXECUTABLE SolveXcts)

add_spectre_executable(
  ${EXECUTABLE}
  EXCLUDE_FROM_ALL
  SolveXcts.cpp
  )

target_link_libraries(
  ${EXECUTABLE}
  PRIVATE
  AmrCriteria
  Charmxx::main
  Convergence
  CoordinateMaps
  DataStructures
  DiscontinuousGalerkin
  Domain
  DomainCreators
  Elliptic
  EllipticDg
  EllipticDgSubdomainOperator
  EllipticSubdomainPreconditioners
  ErrorHandling
  Events
  EventsAndTriggers
  FunctionsOfTime
  Hydro
  Informer
  InitialDataUtilities
  Initialization
  LinearOperators
  MathFunctions
  Observer
  Options
  Parallel
  ParallelAmr
  ParallelLinearSolver
  ParallelMultigrid
  ParallelNonlinearSolver
  ParallelSchwarz
  PhaseControl
  Utilities
  Xcts
  XctsAnalyticData
  XctsBoundaryConditions
  XctsEvents
  XctsSolutions
  )

# Vacuum XCTS solver (z4c-CCM mission): same as SolveXcts but without the
# HydroQuantitiesCompute observation machinery, so vacuum free-data backgrounds
# (e.g. Xcts::AnalyticData::TeukolskyWave for the arXiv:2308.10361 X=2 test)
# need not implement matter variables.
set(EXECUTABLE SolveXctsVacuum)

add_spectre_executable(
  ${EXECUTABLE}
  EXCLUDE_FROM_ALL
  SolveXctsVacuum.cpp
  )

target_link_libraries(
  ${EXECUTABLE}
  PRIVATE
  AmrCriteria
  Charmxx::main
  Convergence
  CoordinateMaps
  DataStructures
  DiscontinuousGalerkin
  Domain
  DomainCreators
  Elliptic
  EllipticDg
  EllipticDgSubdomainOperator
  EllipticSubdomainPreconditioners
  ErrorHandling
  Events
  EventsAndTriggers
  FunctionsOfTime
  Hydro
  Informer
  InitialDataUtilities
  Initialization
  LinearOperators
  MathFunctions
  Observer
  Options
  Parallel
  ParallelAmr
  ParallelLinearSolver
  ParallelMultigrid
  ParallelNonlinearSolver
  ParallelSchwarz
  PhaseControl
  Utilities
  Xcts
  XctsAnalyticData
  XctsBoundaryConditions
  XctsEvents
  XctsSolutions
  )
