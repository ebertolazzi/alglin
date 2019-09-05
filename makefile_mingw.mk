#      # #    # #    # #    #
#      # ##   # #    #  #  #
#      # # #  # #    #   ##
#      # #  # # #    #   ##
#      # #   ## #    #  #  #
###### # #    #  ####  #    #

WARN  = -Wall
CC    = gcc $(WARN)
CXX   = g++ $(WARN) -std=c++11 -pthread
F90   = gfortran
LIBS  = -L./lib -lAlglin -Llib3rd/lib
INC   = -I./src -Ilib3rd/include
CLIBS = -lc++
DEFS  =

CXXFLAGS = -m64 -O3 -funroll-loops -fPIC
#
# activate C++11 for g++ >= 4.9
#
VERSION  = $(shell $(CC) -dumpversion)
AR       = ar rcs
LIBSGCC  = -lsuperlu_mingw -lstdc++ -lm -pthread

ifneq (,$(findstring ALGLIN_USE_LAPACK,$(USED_LIB)))
  override LIBS += -llapack -lblas
endif

ifneq (,$(findstring ALGLIN_USE_OPENBLAS,$(USED_LIB)))
  FPATH=$(dir $(shell gfortran -print-libgcc-file-name))
  override LIBS += -Llib3rd/lib -Wl,-rpath,lib3rd/dll -lopenblas -L$(FPATH)/../../.. -Wl,-rpath,$(LIB3RD)/../../..  -lgfortran
endif

ifneq (,$(findstring ALGLIN_USE_ATLAS,$(USED_LIB)))
  ATLAS_PATH = /usr/lib/atlas-base
  ATLAS_LIBS = -llapack -lf77blas -lcblas -latlas -lgfortran
  override LIBS += -L$(ATLAS_PATH) $(ATLAS_LIBS) -Wl,-rpath,$(ATLAS_PATH)
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

LIBNAME = Alglin_mingw

lib: config lib/lib/lib$(LIBNAME)_static.a lib/lib/lib$(LIBNAME).so

lib/lib/lib$(LIBNAME)_static.a: $(OBJS)
	$(MKDIR) ./lib
	$(MKDIR) ./lib/lib
	$(AR) lib/lib/lib$(LIBNAME)_static.a $(OBJS)

lib/lib/lib$(LIBNAME).so: $(OBJS)
	$(MKDIR) ./lib
	$(MKDIR) ./lib/lib
	$(CXX) -shared -o ./lib/lib/lib$(LIBNAME).so $(OBJS)

install_local: lib/lib/lib$(LIBNAME).a lib/lib/lib$(LIBNAME).so
	$(MKDIR) ./lib/include
	cp -f -P src/*.h*       ./lib/include
	cp -rf lib3rd/include/* ./lib/include

install: lib/lib/lib$(LIBNAME).a lib/lib/lib$(LIBNAME).so
	$(MKDIR) $(PREFIX)/include
	cp -f -P src/*.h*             $(PREFIX)/include
	cp -f -P lib/lib/libAlglin.a  $(PREFIX)/lib
	cp -f -P lib/lib/libAlglin.so $(PREFIX)/lib

install_as_framework: lib/$(LIB_ALGLIN)
	$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	$(MKDIR) $(PREFIX)/lib
	cp -f -P src/*.h*          $(PREFIX)/include/$(FRAMEWORK)
	cp -rf lib3rd/include/*    $(PREFIX)/include/$(FRAMEWORK)
	cp -f -P lib/$(LIB_ALGLIN) $(PREFIX)/lib/$(LIB_ALGLIN)
