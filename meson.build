# Build file for zreion, building fortran to use from python.
# Most of this is copied from the auto-generated build file from f2py.
project(
  'zreion',
  ['c', 'fortran'],
)

fc = meson.get_compiler('fortran')

py = import('python').find_installation(pure: false)

# Dependencies
py_dep = py.dependency()

incdir_numpy = run_command(
  py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true,
).stdout().strip()

incdir_f2py = run_command(
  py,
  ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
  check : true,
).stdout().strip()

inc_np = include_directories(incdir_numpy)
incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

# MKL include
# TODO: this is hard-coded for building using ifort on Bridges-2
incdir_mkl = '/opt/packages/oneapi/v2023.2.0/mkl/2023.2.0/include/'

fortran_include = include_directories(incdir_numpy, incdir_f2py, incdir_mkl)

# Libraries
# TODO: this is hard-coded for building using ifort on Bridges-2
mkl_intel_ilp64 = declare_dependency(link_args : ['-lmkl_intel_ilp64'])
mkl_intel_thread = declare_dependency(link_args : ['-lmkl_intel_thread'])
mkl_core = declare_dependency(link_args : ['-lmkl_core'])
iomp5 = declare_dependency(link_args : ['-liomp5'])
lib_dir_0 = declare_dependency(link_args : ['-L/opt/packages/oneapi/v2023.2.0/compiler/2023.2.1/linux/compiler/lib/intel64_lin'])
lib_dir_1 = declare_dependency(link_args : ['-L/opt/packages/oneapi/v2023.2.0/mkl/2023.2.0/lib/intel64'])

# Make f2py do some up-front work
input_files = [  # order is important!
    'fortran_source/fft_tools.f90',
    'fortran_source/general_tools.f90',
    'fortran_source/cosmo_tools.f90',
    'fortran_source/ic_tools.f90',
    'fortran_source/domain_tools.f90',
    'fortran_source/particle_tools.f90',
    'fortran_source/zreion_tools.f90',
    'fortran_source/field_tools.f90',
]
module_name = 'fortran_tools'
f2py_output = custom_target(
  input : input_files,
  output : [module_name + '-f2pywrappers2.f90', module_name + 'module.c'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', module_name, '--lower'],
)

# Build the python extension
py.extension_module(module_name,
                    [input_files, f2py_output, fortranobject_c],
                    include_directories : [
                      fortran_include,
                    ],
                    dependencies : [
                      py_dep,
                      mkl_intel_ilp64,
                      mkl_intel_thread,
                      mkl_core,
                      iomp5,
                      lib_dir_0,
                      lib_dir_1,
                    ],
                    install : true,
                    subdir: 'zreion')

# The actual python install is done via a helper build
subdir('zreion')
