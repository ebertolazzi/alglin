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

CXXFLAGS = -pthread -msse4.2 -msse4.1 -mssse3 -msse3 -msse2 -msse -mmmx -m64 -O3 -funroll-loops -fPIC
override LIBS += -L./lib -lAlglin
override INC  += -I./src

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  WARN = -Wall
  CC   = gcc
  CXX  = g++
  # activate C++11 for g++ >= 4.9
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
  CC     += $(WARN)
  CXX    += $(WARN)
  AR      = ar rcs
  LIBSGCC = -lstdc++ -lm
ifeq ($(ATLAS),1)
  # for ATLAS (default)
  override LIBS += -L/usr/lib/atlas-base -Wl,-rpath,/usr/lib/atlas-base -llapack -lf77blas -lcblas -latlas -lopenblas
  USED_LIB = ALGLIN_USE_ATLAS
else
#
ifeq ($(MKL),1)
  # for MKL
  MKL_PATH = /opt/intel/mkl
  #MKL_LIB = -lmkl_tbb_thread -lmkl_rt -lmkl_core
  MKL_LIB = -lmkl_sequential -lmkl_rt -lmkl_core
  override LIBS += -L$(MKL_PATH)/lib/intel64 -Wl,-rpath,$(MKL_PATH)/lib/intel64 $(MKL_LIB)
  override INC  += -I$(MKL_PATH)/include
  USED_LIB = ALGLIN_USE_MKL
else
  # for OPENBLAS
  override LIBS += -L/usr/lib/openblas-base -Wl,-rpath,/usr/lib/openblas-base -lopenblas
  USED_LIB = ALGLIN_USE_OPENBLAS
endif
#
endif
  override INC  += -I/usr/include/eigen3 -I/usr/include/atlas/
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  WARN     = -Weverything -Wno-reserved-id-macro -Wno-padded
  CC       = clang
  CXX      = clang++
  VERSION  = $(shell $(CC) --version 2>&1 | grep -o "Apple LLVM version [0-9]\.[0-9]\.[0-9]" | grep -o " [0-9]\.")
ifneq (,$(findstring 8., $(VERSION)))
  CXX += -std=c++11 -stdlib=libc++
  THREAD = ALGLIN_USE_THREAD
else
ifneq (,$(findstring 7., $(VERSION)))
  CXX += -std=c++11 -stdlib=libc++
  THREAD = ALGLIN_USE_THREAD
else
  CXX += -std=c++11
  THREAD = ALGLIN_USE_THREAD
endif
endif
  CC     += $(WARN)
  CXX    += $(WARN)
  AR      = libtool -static -o
  LIBSGCC = -lstdc++ -lm
ifeq ($(OPENBLAS),1)
  # for OPENBLAS
  override LIBS += -L/usr/local/opt/openblas/lib -lopenblas
  override INC  += -I/usr/local/opt/openblas/include
  USED_LIB = ALGLIN_USE_OPENBLAS
else
#
ifeq ($(MKL),1)
  # for MKL
  MKL_PATH = /opt/intel/mkl
  #MKL_LIB = -lmkl_tbb_thread -lmkl_intel -lmkl_core
  MKL_LIB = -lmkl_sequential -lmkl_intel -lmkl_core
  override LIBS += -L$(MKL_PATH)/lib -Xlinker -rpath -Xlinker $(MKL_PATH)/lib $(MKL_LIB)
  override INC  += -I$(MKL_PATH)/include
  USED_LIB = ALGLIN_USE_MKL
else
  override LIBS += -L./lib -lAlglin -framework Accelerate
  USED_LIB = ALGLIN_USE_ACCELERATE
endif
#
endif
  override INC += -I/usr/local/include/eigen3
endif

CC  += -O3 -g0
CXX += -O3 -g0
#CC  += -O1 -g3
#CXX += -O1 -g3

SRCS_ALGLIN = \
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
src/KKT_like.cc

SRCS_SUPERLU=\
src/BABD_SuperLU.cc \
src/Simplex.cc

#SRCS = $(SRCS_ALGLIN) $(SRCS_SUPERLU)
SRCS = $(SRCS_ALGLIN)

OBJS  = $(SRCS:.cc=.o)
DEPS  = \
src/ABD_Arceco.hh \
src/ABD_Diaz.hh \
src/ABD_Block.hh \
src/Alglin++.hh \
src/Alglin.hh \
src/AlglinConfig.hh \
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

all: config lib
	mkdir -p bin
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test0-FD                  src_tests/test0-FD.cc $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test1-small-factorization src_tests/test1-small-factorization.cc $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test2-Threads             src_tests/test2-Threads.cc $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test3-Timing              src_tests/test3-Timing.cc $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test4-KKT                 src_tests/test4-KKT.cc $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test5-ABD-Diaz            src_tests/test5-ABD-Diaz.cc $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test6-ABD-Block           src_tests/test6-ABD-Block.cc $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test7-BorderedCR          src_tests/test7-BorderedCR.cc $(LIBS)
	$(CC)  $(INC) $(DEFS) $(CXXFLAGS) -o bin/test8-Cinterface          src_tests/test8-Cinterface.c $(LIBS) $(LIBSGCC)
	$(CC)  $(INC) $(DEFS) $(CXXFLAGS) -o bin/test9-Cinterface          src_tests/test9-Cinterface.c $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test12-BandedMatrix       src_tests/test12-BandedMatrix.cc $(LIBS) $(LIBSGCC)

all1: config lib
	mkdir -p bin
	$(F90) $(INC) -o bin/test10-FORTRAN src_tests/test10-FORTRAN.f90 $(LIBS) $(LIBSGCC) $(CLIBS)
	$(F90) $(INC) -o bin/test11-FORTRAN src_tests/test11-FORTRAN.f90 $(LIBS) $(LIBSGCC) $(CLIBS)

all_simplex: libAlglin
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/SimplexTest1              src_tests/SimplexTest1.cc $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/SimplexTest2              src_tests/SimplexTest2.cc $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/SimplexTest3              src_tests/SimplexTest3.cc $(LIBS)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/SimplexTest4              src_tests/SimplexTest4.cc $(LIBS)

lib: lib/$(LIB_ALGLIN)

src/%.o: src/%.cc $(DEPS)
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@

src/%.o: src/%.c $(DEPS)
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
	cp src/*.hh ./lib/include
	cp src/*.h  ./lib/include

install: lib/$(LIB_ALGLIN)
	$(MKDIR) $(PREFIX)/include
	cp src/*.hh          $(PREFIX)/include
	cp src/*.h           $(PREFIX)/include
	cp lib/$(LIB_ALGLIN) $(PREFIX)/lib

install_as_framework: lib/$(LIB_ALGLIN)
	$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	cp src/*.hh          $(PREFIX)/include/$(FRAMEWORK)
	cp lib/$(LIB_ALGLIN) $(PREFIX)/lib

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
	rm -rf lib/libAlglin.* lib/include src/*.o
	rm -rf bin
