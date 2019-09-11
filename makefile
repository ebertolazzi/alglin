# get the type of OS currently running
OS=$(shell uname)
PWD=$(shell pwd)

include Makefile_config.mk

.SUFFIXES:
.SUFFIXES: .o .c .cc

SRCS = \
src/Alglin_aux.cc \
src/ABD_Arceco.cc \
src/ABD_Block.cc \
src/ABD_Diaz.cc \
src/BABD.cc \
src/BABD_Block.cc \
src/BABD_BorderedCR.cc \
src/BABD_C_interface.cc \
src/BABD_SuperLU.cc \
src/BlockBidiagonal.cc \
src/KKT_like.cc \
src/Simplex.cc

OBJS = $(SRCS:.cc=.o)

#src/AlglinConfig.hh
DEPS = \
src/ABD_Arceco.hh \
src/ABD_Block.hh \
src/ABD_Diaz.hh \
src/Alglin.hh \
src/Alglin_Config.hh \
src/Alglin_Eigen.hh \
src/Alglin_FD.hh \
src/Alglin_SuperLU.hh \
src/Alglin_aux.hh \
src/Alglin_tmpl.hh \
src/BABD.hh \
src/BABD_Block.hh \
src/BABD_BorderedCR.hh \
src/BABD_C_interface.h \
src/BABD_SuperLU.hh \
src/BlockBidiagonal.hh \
src/KKT_like.hh \
src/Simplex.hh

MKDIR     = mkdir -p
PREFIX    = /usr/local
FRAMEWORK = Alglin

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
include Makefile_linux.mk
endif

# check if the OS string contains 'MINGW'
ifneq (,$(findstring MINGW, $(OS)))
include Makefile_mingw.mk
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
include Makefile_osx.mk
endif

all: config lib
	mkdir -p bin
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test0-FD                  src_tests/test0-FD.cc                  $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test1-small-factorization src_tests/test1-small-factorization.cc $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test2-Threads             src_tests/test2-Threads.cc             $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test3-Timing              src_tests/test3-Timing.cc              $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test4-KKT                 src_tests/test4-KKT.cc                 $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test5-ABD-Diaz            src_tests/test5-ABD-Diaz.cc            $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test6-ABD-Block           src_tests/test6-ABD-Block.cc           $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test7-BorderedCR          src_tests/test7-BorderedCR.cc          $(LIBS) $(LIBSGCC)
	$(CC)  $(INC) $(DEFS) $(CXXFLAGS) -o bin/test8-Cinterface          src_tests/test8-Cinterface.c           $(LIBS) $(LIBSGCC)
	$(CC)  $(INC) $(DEFS) $(CXXFLAGS) -o bin/test9-Cinterface          src_tests/test9-Cinterface.c           $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test12-BandedMatrix       src_tests/test12-BandedMatrix.cc       $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test13-BFGS               src_tests/test13-BFGS.cc               $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test14-BLOCKTRID          src_tests/test14-BLOCKTRID.cc          $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/test15-EIGS               src_tests/test15-EIGS.cc               $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/Simplex-Test1             src_tests/Simplex-Test1.cc             $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/Simplex-Test2             src_tests/Simplex-Test2.cc             $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/Simplex-Test3             src_tests/Simplex-Test3.cc             $(LIBS) $(LIBSGCC)
	$(CXX) $(INC) $(DEFS) $(CXXFLAGS) -o bin/Simplex-Test4             src_tests/Simplex-Test4.cc             $(LIBS) $(LIBSGCC)

all1: config lib
	mkdir -p bin
	$(F90) $(INC) -o bin/test10-FORTRAN src_tests/test10-FORTRAN.f90 $(LIBS) $(LIBSGCC) $(CLIBS)
	$(F90) $(INC) -o bin/test11-FORTRAN src_tests/test11-FORTRAN.f90 $(LIBS) $(LIBSGCC) $(CLIBS)

.cc.o: $(DEPS)
	$(CXX) $(INC) $(CXXFLAGS) $(DEFS) -c $< -o $@

.c.o: $(DEPS)
	$(CC) $(INC) $(CFLAGS) $(DEFS) -c -o $@ $<

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
	./bin/test13-BFGS
	./bin/test14-BLOCKTRID
	./bin/test15-EIGS

run_simplex:
	./bin/Simplex-Test1
	./bin/Simplex-Test2
	./bin/Simplex-Test3
	./bin/Simplex-Test4

doc:
	doxygen

clean:
	rm -rf lib/libAlglin.* src/*.o
	rm -rf bin
