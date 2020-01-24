#      # #    # #    # #    #
#      # ##   # #    #  #  #
#      # # #  # #    #   ##
#      # #  # # #    #   ##
#      # #   ## #    #  #  #
###### # #    #  ####  #    #

WARN    = -Wall
CC      = gcc $(WARN)
CXX     = g++ $(WARN) -std=c++11 -pthread
F90     = gfortran
LIBS3RD = -Llib3rd/lib -Llib3rd/dll -Wl,-rpath,lib3rd/dll -llapack_wrapper_linux_static -lsuperlu_linux_static
LIBS    = -Llib/lib -Llib/dll -lAlglin_linux $(LIBS3RD) -Wl,-rpath,lib/dll -ldl
INC     = -I./src -Ilib3rd/include
CLIBS   = -lc++
DEFS    =

CXXFLAGS = -O2 -funroll-loops -floop-interchange -floop-block -fPIC
#
# activate C++11 for g++ >= 4.9
#
VERSION  = $(shell $(CC) -dumpversion)
AR       = ar rcs
LIBSGCC  = -lstdc++ -lm -pthread

ifneq (,$(findstring ALGLIN_USE_LAPACK,$(USED_LIB)))
  override LIBS += -llapack -lblas
endif

ifneq (,$(findstring ALGLIN_USE_OPENBLAS,$(USED_LIB)))
  FPATH=$(dir $(shell gfortran -print-libgcc-file-name))
  override LIBS += -lopenblas -L$(FPATH) -lgfortran
endif

ifneq (,$(findstring ALGLIN_USE_ATLAS,$(USED_LIB)))
  FPATH      = $(dir $(shell gfortran -print-libgcc-file-name))
  ATLAS_PATH = /usr/lib/atlas-base
  ATLAS_LIBS = -llapack -lf77blas -lcblas -latlas -lgfortran
  override LIBS += -L$(ATLAS_PATH) $(ATLAS_LIBS) -Wl,-rpath,$(ATLAS_PATH) -L$(FPATH) -lgfortran
endif

ifneq (,$(findstring ALGLIN_USE_MKL,$(USED_LIB)))
  # for MKL
  MKL_PATH = /opt/intel/mkl
  MKL_ARCH = intel64
  #MKL_LIB = -lmkl_tbb_thread -lmkl_rt -lmkl_core
  MKL_LIB = -lmkl_sequential -lmkl_rt -lmkl_core
  override LIBS += -L$(MKL_PATH)/lib/$(MKL_ARCH) -Wl,-rpath,$(MKL_PATH)/lib/$(MKL_ARCH) $(MKL_LIB)
  override INC  += -I$(MKL_PATH)/include
endif

ifneq (,$(findstring ALGLIN_USE_ACCELERATE,$(USED_LIB)))
  $(error error is "Accelerate is supported only on Darwin!")
endif

LIBNAME = Alglin_linux

lib: lib/lib/lib$(LIBNAME)_static.a lib/dll/lib$(LIBNAME).so

lib/lib/lib$(LIBNAME)_static.a: $(OBJS)
	$(MKDIR) ./lib
	$(MKDIR) ./lib/lib
	$(AR) lib/lib/lib$(LIBNAME)_static.a $(OBJS)

lib/dll/lib$(LIBNAME).so: $(OBJS)
	$(MKDIR) ./lib
	$(MKDIR) ./lib/dll
	$(CXX) -shared -o ./lib/dll/lib$(LIBNAME).so $(OBJS)

install_local: lib/lib/lib$(LIBNAME).a lib/dll/lib$(LIBNAME).so
	$(MKDIR) ./lib/include
	cp -f -P src/*.h*       ./lib/include
	cp -rf lib3rd/include/* ./lib/include

install: lib/lib/lib$(LIBNAME).a lib/dll/lib$(LIBNAME).so
	$(MKDIR) $(PREFIX)/include
	cp -f -P src/*.h*             $(PREFIX)/include
	cp -f -P lib/lib/libAlglin.a  $(PREFIX)/lib
	cp -f -P lib/dll/libAlglin.so $(PREFIX)/dll
