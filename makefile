# get the type of OS currently running
OS=$(shell uname)
PWD=$(shell pwd)

LIB_ALGLIN = libAlglin.a

CC    = gcc
CXX   = g++
F90   = gfortran
INC   =
LIBS  =
CLIBS = -lc++
DEFS  =

CXXFLAGS = -msse4.2 -msse4.1 -mssse3 -msse3 -msse2 -msse -mmmx -m64 -O3 -funroll-loops -fPIC
override LIBS += -L./lib -lAlglin -L$(PWD)/lib3rd/lib
override INC  += -I./src -I$(PWD)/lib3rd/include

#
# select which version of BLAS/LAPACK use
#
USED_LIB=""
ifeq ($(ATLAS),1)
  USED_LIB = ALGLIN_USE_ATLAS
endif
#
ifeq ($(MKL),1)
  USED_LIB = ALGLIN_USE_MKL
endif
#
ifeq ($(OPENBLAS),1)
  USED_LIB = ALGLIN_USE_OPENBLAS
endif
#
ifeq ($(LAPACK),1)
  USED_LIB = ALGLIN_USE_LAPACK
endif
#
ifeq ($(ACCELERATE),1)
  USED_LIB = ALGLIN_USE_ACCELERATE
endif
#
# if missig setup default
#
ifeq ($(USED_LIB), "")
ifeq (,$(wildcard .alglin_config))
ifneq (,$(findstring Darwin, $(OS)))
  USED_LIB = ALGLIN_USE_ACCELERATE
else
  USED_LIB = ALGLIN_USE_LAPACK
endif
else
 USED_LIB = $(shell cat .alglin_config)
endif
endif

$(shell echo "$(USED_LIB)" > .alglin_config )
$(info $(USED_LIB))

#      # #    # #    # #    #
#      # ##   # #    #  #  #
#      # # #  # #    #   ##
#      # #  # # #    #   ##
#      # #   ## #    #  #  #
###### # #    #  ####  #    #

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  WARN = -Wall
  CC   = gcc
  CXX  = g++
  #
  # activate C++11 for g++ >= 4.9
  #
  VERSION  = $(shell $(CC) -dumpversion)
ifneq (,$(findstring 4.9, $(VERSION)))
  CXX += -std=c++11 -pthread
  THREAD = ALGLIN_USE_THREAD
else
ifneq (,$(findstring 5., $(VERSION)))
  CXX += -std=c++11 -pthread
  THREAD = ALGLIN_USE_THREAD
else
ifneq (,$(findstring 6., $(VERSION)))
  CXX += -std=c++11 -pthread
  THREAD = ALGLIN_USE_THREAD
else
  THREAD = ALGLIN_DO_NOT_USE_CXX11
endif
endif
endif
  #
  #
  #
  CC     += $(WARN)
  CXX    += $(WARN)
  AR      = ar rcs
  LIBSGCC = -lstdc++ -lm -pthread

ifneq (,$(findstring ALGLIN_USE_LAPACK,$(USED_LIB)))
  override LIBS += -llapack -lblas
endif

ifneq (,$(findstring ALGLIN_USE_OPENBLAS,$(USED_LIB)))
  OPENBLAS_PATH = $(PWD)/lib3rd/lib/openblas
  override LIBS += -L$(OPENBLAS_PATH) -Wl,-rpath,$(OPENBLAS_PATH) -lopenblas
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

  #
  override INC  += -I/usr/include/eigen3 -I/usr/include/atlas/
endif

#######  #####  #     #
#     # #     #  #   #
#     # #         # #
#     #  #####     #
#     #       #   # #
#     # #     #  #   #
#######  #####  #     #

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  WARN     = -Weverything -Wno-reserved-id-macro -Wno-padded
  CC       = clang
  CXX      = clang++
  VERSION  = $(shell $(CC) --version 2>&1 | grep -o "Apple LLVM version [0-9]\.[0-9]\.[0-9]" | grep -o " [0-9]\.")
  #---------
  CXX     += -std=c++11 -stdlib=libc++
  THREAD   = ALGLIN_USE_THREAD
  #---------
  CC     += $(WARN)
  CXX    += $(WARN)
  AR      = libtool -static -o
  LIBSGCC = -lstdc++ -lm

ifneq (,$(findstring ALGLIN_USE_LAPACK,$(USED_LIB)))
  override LIBS += -llapack -lblas
endif

ifneq (,$(findstring ALGLIN_USE_OPENBLAS,$(USED_LIB)))
  override LIBS += -L$(PWD)/lib3rd/lib/openblas -Xlinker -rpath -Xlinker $(PWD)/lib3rd/lib/openblas -lopenblas
endif

ifneq (,$(findstring ALGLIN_USE_ATLAS,$(USED_LIB)))
  $(error error is "Atlas is supported only on Linux!")
endif

ifneq (,$(findstring ALGLIN_USE_MKL,$(USED_LIB)))
  # for MKL
  MKL_PATH = /opt/intel/mkl
  #MKL_LIB = -lmkl_tbb_thread -lmkl_intel -lmkl_core
  MKL_LIB = -lmkl_sequential -lmkl_intel -lmkl_core
  override LIBS += -L$(MKL_PATH)/lib -Xlinker -rpath -Xlinker $(MKL_PATH)/lib $(MKL_LIB)
  override INC  += -I$(MKL_PATH)/include
endif

ifneq (,$(findstring ALGLIN_USE_ACCELERATE,$(USED_LIB)))
  override LIBS += -L./lib -lAlglin -framework Accelerate
endif

  override INC += -I/usr/local/include/eigen3
endif

SRCS = \
src/ABD_Arceco.cc \
src/ABD_Block.cc \
src/ABD_Diaz.cc \
src/Alglin++.cc \
src/Alglin.cc \
src/BlockBidiagonal.cc \
src/BABD.cc \
src/BABD_Block.cc \
src/BABD_C_interface.cc \
src/BABD_BorderedCR.cc \
src/KKT_like.cc \
src/BABD_SuperLU.cc \
src/Simplex.cc
OBJS = $(SRCS:.cc=.o)

SRCS_TESTS = \
src_tests/test0-FD.cc \
src_tests/test1-small-factorization.cc \
src_tests/test2-Threads.cc \
src_tests/test3-Timing.cc \
src_tests/test4-KKT.cc \
src_tests/test5-ABD-Diaz.cc \
src_tests/test6-ABD-Block.cc \
src_tests/test7-BorderedCR.cc \
src_tests/test12-BandedMatrix.cc

OBJS_TESTS = $(SRCS_TESTS:.cc=.o)

#src/AlglinConfig.hh
DEPS = \
src/ABD_Arceco.hh \
src/ABD_Diaz.hh \
src/ABD_Block.hh \
src/Alglin++.hh \
src/Alglin.hh \
src/AlglinEigen.hh \
src/AlglinFD.hh \
src/AlglinSuperLU.hh \
src/Alglin_aux.hh \
src/Alglin_threads.hh \
src/Alglin_tmpl.hh \
src/BABD.hh \
src/BABD_Block.hh \
src/BABD_SuperLU.hh \
src/BlockBidiagonal.hh \
src/BABD_BorderedCR.hh \
src/KKT_like.hh \
src/Simplex.hh \
src/TicToc.hh

MKDIR = mkdir -p

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX    = /usr/local
FRAMEWORK = Alglin

all: config lib $(OBJS_TESTS)
	mkdir -p bin
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test0-FD                  src_tests/test0-FD.o $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test1-small-factorization src_tests/test1-small-factorization.o $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test2-Threads             src_tests/test2-Threads.o $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test3-Timing              src_tests/test3-Timing.o $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test4-KKT                 src_tests/test4-KKT.o $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test5-ABD-Diaz            src_tests/test5-ABD-Diaz.o $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test6-ABD-Block           src_tests/test6-ABD-Block.o $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test7-BorderedCR          src_tests/test7-BorderedCR.o $(LIBS)
	$(CC)  $(INC) $(DEFS) $(CXXFLAGS) -o bin/test8-Cinterface          src_tests/test8-Cinterface.c $(LIBS) $(LIBSGCC)
	$(CC)  $(INC) $(DEFS) $(CXXFLAGS) -o bin/test9-Cinterface          src_tests/test9-Cinterface.c $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test12-BandedMatrix       src_tests/test12-BandedMatrix.o $(LIBS) $(LIBSGCC)

all1: config lib
	mkdir -p bin
	$(F90) $(INC) -o bin/test10-FORTRAN src_tests/test10-FORTRAN.f90 $(LIBS) $(LIBSGCC) $(CLIBS)
	$(F90) $(INC) -o bin/test11-FORTRAN src_tests/test11-FORTRAN.f90 $(LIBS) $(LIBSGCC) $(CLIBS)

all_simplex: config lib
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/SimplexTest1 src_tests/SimplexTest1.cc $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/SimplexTest2 src_tests/SimplexTest2.cc $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/SimplexTest3 src_tests/SimplexTest3.cc $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/SimplexTest4 src_tests/SimplexTest4.cc $(LIBS)

lib: config lib/$(LIB_ALGLIN)

src/%.o: src/%.cc $(DEPS)
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@

src/%.o: src/%.c $(DEPS)
	$(CC) $(INC) $(CFLAGS) $(DEFS) -c -o $@ $<

src_tests/%.o: src_tests/%.cc $(DEPS)
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@

src_tests/%.o: src_tests/%.c $(DEPS)
	$(CC) $(INC) $(CFLAGS) $(DEFS) -c -o $@ $<

lib/libAlglin.a: $(OBJS)
	$(MKDIR) lib
	$(AR) lib/libAlglin.a $(OBJS)

lib/libAlglin.dylib: $(OBJS)
	$(MKDIR) lib
	$(CXX) -shared -o lib/libAlglin.dylib $(OBJS)

lib/libAlglin.so: $(OBJS)
	$(MKDIR) lib
	$(CXX) -shared -o lib/libAlglin.so $(OBJS)

install_local: lib/$(LIB_ALGLIN)
	$(MKDIR) ./lib/include
	cp -f -P src/*.h*          ./lib/include
	cp -f -P lib3rd/include/*  ./lib/include

install: lib/$(LIB_ALGLIN)
	$(MKDIR) $(PREFIX)/include
	cp -f -P src/*.h*          $(PREFIX)/include
	#cp -f -P lib3rd/include/*  $(PREFIX)/include
	cp -f -P lib/$(LIB_ALGLIN) $(PREFIX)/lib

install_as_framework: lib/$(LIB_ALGLIN)
	$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	cp -f -P src/*.h*          $(PREFIX)/include/$(FRAMEWORK)
	cp -f -P lib3rd/include/*  $(PREFIX)/include/$(FRAMEWORK)
	cp -f -P lib/$(LIB_ALGLIN) $(PREFIX)/lib

config:
	rm -f src/AlglinConfig.hh
	sed 's/@@ALGLIN_USE@@/#define $(USED_LIB) 1/' <src/AlglinConfig.hh.tmpl | \
	sed 's/@@ALGLIN_THREAD@@/#define $(THREAD) 1/' >src/AlglinConfig.hh

run:
	./bin/test0-FD
	./bin/test1-small-factorization
	./bin/test2-Threads
	./bin/test3-Timing
	./bin/test4-KKT
	./bin/test5-ABD-Diaz
	./bin/test6-ABD-Block
	./bin/test7-BorderedCR
	./bin/test8-Cinterface
	./bin/test12-BandedMatrix

run_simplex:
	./bin/SimplexTest1
	./bin/SimplexTest2
	./bin/SimplexTest3
	./bin/SimplexTest4

doc:
	doxygen

clean:
	rm -rf lib/libAlglin.* src/*.o  src_tests/*.o
	rm -rf bin
