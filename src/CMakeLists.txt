# CMake file to build mf2005 executable

set(SOURCE
  de47.f
  gmg7_c_interface.f90
  gmg7.f
  gwf2bas7.f
  gwf2bcf7.f
  gwf2chd7.f
  gwf2drn7.f
  gwf2drt7.f
  gwf2ets7.f
  gwf2evt7.f
  gwf2fhb7.f
  gwf2gag7.f
  gwf2ghb7.f
  gwf2hfb7.f
  gwf2huf7.f
  gwf2hydmod7.f
  gwf2ibs7.f
  gwf2lak7.f
  gwf2lpf7.f
  gwf2mnw17.f
  gwf2mnw27.f
  gwf2mnw2i7.f
  gwf2rch7.f
  gwf2res7.f
  gwf2riv7.f
  gwf2sfr7.f
  gwf2str7.f
  gwf2sub7.f
  gwf2swi27.f
  gwf2swr7.f
  gwf2swr7util.f
  gwf2swt7.f
  gwf2uzf1.f
  gwf2wel7.f
  gwflakmodule.f
  gwfsfrmodule.f
  gwfuzfmodule.f
  hufutl7.f
  lmt8.f
  mf2005.f
  mhc7.f
  obs2bas7.f
  obs2chd7.f
  obs2drn7.f
  obs2ghb7.f
  obs2riv7.f
  obs2str7.f
  parutl7.f
  pcg7.f
  pcgn2.f90
  pcgn_solve2.f90
  sip7.f
  utl7.f
  ccfd.c
  mf2kgmg.c
  r_vector.c
  solvers.c
)

add_executable(mf2005 ${SOURCE})
target_compile_options(mf2005
  PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${Fortran_WARN_FLAGS}>
  PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${Fortran_BUILD_FLAGS}>
)

write_openspec_inc(openspec.inc mf2005)

install(TARGETS mf2005
  RUNTIME DESTINATION bin
)
