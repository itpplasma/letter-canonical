cmake_minimum_required(VERSION 3.10)
project(letter-canonical Fortran C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Check whether intel or GNU
if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS " -diag-disable=10448 -g -traceback -fPIC -qopenmp -xHost \
        -warn all,nounused")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 ${CMAKE_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 ${CMAKE_Fortran_FLAGS}")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "-fPIC -g -fbacktrace -fopenmp -march=native -mtune=native \
        -Wall -Wuninitialized -Wno-maybe-uninitialized -Wno-unused-label \
        -Wno-unused-dummy-argument -fmax-errors=1")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 ${CMAKE_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 ${CMAKE_Fortran_FLAGS}")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(CMAKE_Fortran_FLAGS "-fPIC -g -traceback -mp -Minform=inform")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 ${CMAKE_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 ${CMAKE_Fortran_FLAGS}")
else()
    message(FATAL_ERROR "Unsupported compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif()


set(SOURCES
    test_util.f90
    magfie.f90
    magfie_test.f90
    magfie_tok.f90
    magfie_factory.f90
    psi_transform.f90
    biotsavart.f90
    canonical.f90
    field_can_base.f90
    field_can_test.f90
    field_can_cyl.f90
    field_can_new.f90
    field_can.f90
    integrator/integrator_base.f90
    integrator/integrator_rk45.f90
    integrator/integrator_euler0.f90
    integrator/integrator_euler1.f90
    integrator/velo_sub.f90
    integrator/rk45_cyl_integrator.f90
    integrator/rk45_can_integrator.f90
    integrator/expl_impl_euler_integrator.f90
    integrator/integrator.f90
    callback/callback_base.f90
    callback/cut_callback.f90
    callback/callback.f90
    letter_canonical.f90
)

set(LIBNEO_SOURCES
    ../libneo/src/libneo_kinds.f90
    ../libneo/src/math_constants.f90
    ../libneo/src/spl_three_to_five.f90
    ../libneo/src/odeint_rkf45.f90
    ../libneo/src/contrib/rkf45.f90
    ../libneo/src/interpolate.f90
    ../libneo/src/binsrc.f90
    ../libneo/src/plag_coeff.f90
)

set(MAGFIE_SOURCES
    ../libneo/src/magfie/spline5_RZ.f90
    ../libneo/src/magfie/theta_rz_mod.f90
    ../libneo/src/magfie/extract_fluxcoord_mod.f90
    ../libneo/src/magfie/amn_mod.f90
    ../libneo/src/magfie/field_mod.f90
    ../libneo/src/magfie/field_eq_mod.f90
    ../libneo/src/magfie/field_c_mod.f90
    ../libneo/src/magfie/input_files.f90
    ../libneo/src/magfie/inthecore_mod.f90
    ../libneo/src/magfie/field_divB0.f90
    ../libneo/src/magfie/bdivfree_mod.f90
    ../libneo/src/magfie/bdivfree.f90
)

add_library(field SHARED ${LIBNEO_SOURCES} ${MAGFIE_SOURCES})
add_library(canonical SHARED ${SOURCES})
target_link_libraries(canonical field)

add_executable(letter_canonical.x main.f90)
target_link_libraries(letter_canonical.x canonical)

add_executable(test_orbit_sympl.x test_orbit_sympl.f90)
target_link_libraries(test_orbit_sympl.x canonical)

add_executable(test_orbit.x test_orbit.f90)
target_link_libraries(test_orbit.x canonical)

add_executable(test_callback.x test_callback.f90)
target_link_libraries(test_callback.x canonical)

add_executable(test_integrator.x test_integrator.f90)
target_link_libraries(test_integrator.x canonical)

add_executable(test_field_can.x test_field_can.f90)
target_link_libraries(test_field_can.x canonical)

add_executable(test_magfie.x test_magfie.f90)
target_link_libraries(test_magfie.x canonical)

add_executable(test.x test.f90)
target_link_libraries(test.x canonical)

add_executable(test_large.x test_large.f90)
target_link_libraries(test_large.x canonical)

add_executable(test_biotsavart.x test_biotsavart.f90)
target_link_libraries(test_biotsavart.x canonical)
