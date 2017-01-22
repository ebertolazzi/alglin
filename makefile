# get the type of OS currently running
OS=$(shell uname)
PWD=$(shell pwd)

LIB_ALGLIN = libAlglin.a

CC   = gcc
CXX  = g++
INC  = -I./src -I./include
LIBS = -L./lib -lAlglin
DEFS =

# check if the OS string contains 'Linux'
ifneq (,$(findstring Linux, $(OS)))
  INC     += -I/usr/include/atlas -I/usr/include/eigen3
  CXXFLAGS = -Wall -O3 -fPIC -Wno-sign-compare -std=c++11
  AR       = ar rcs
  LIBS     = -L./lib -lAlglin -llapack -lblas
endif

# check if the OS string contains 'Darwin'
ifneq (,$(findstring Darwin, $(OS)))
  WARN     = -Weverything -Wno-reserved-id-macro -Wno-padded -Wno-documentation-unknown-command -Wno-float-equal -Wimplicit-fallthrough
  CC       = clang
  CXX      = clang++ -std=c++11
  INC     += -I/usr/local/include -I/usr/local/include/eigen3
  CXXFLAGS = -Wall -O3 -fPIC -Wno-sign-compare
  AR       = libtool -static -o
  LIBS     = -L./lib -lAlglin -framework Accelerate
endif

SRCS = \
src/ABD_Arceco.cc \
src/ABD_Diaz.cc \
src/Alglin++.cc \
src/Alglin.cc \
src/BABD.cc \
src/BABD_Amodio.cc \
src/BABD_AmodioN.cc \
src/BABD_Block.cc \
src/BABD_QR.cc \
src/BABD_QR_N.cc \
src/BABD_SuperLU.cc \
src/BlockBidiagonal.cc \
src/BorderedCR.cc \
src/CyclicReduction.cc \
src/CyclicReductionN.cc \
src/CyclicReductionQR.cc \
src/KKT_like.cc

OBJS  = $(SRCS:.cc=.o)
DEPS  = \
src/ABD_Arceco.hh \
src/ABD_Diaz.hh \
src/Alglin++.hh \
src/Alglin.hh \
src/AlglinConfig.hh \
src/AlglinEigen.hh \
src/AlglinFD.hh \
src/AlglinSuperLU.hh \
src/Alglin_aux.hh \
src/Alglin_tmpl.hh \
src/BABD.hh \
src/BABD_Amodio.hh \
src/BABD_AmodioN.hh \
src/BABD_Block.hh \
src/BABD_QR.hh \
src/BABD_QR_N.hh \
src/BABD_SuperLU.hh \
src/BlockBidiagonal.hh \
src/BorderedCR.hh \
src/CyclicReduction.hh \
src/CyclicReductionQR.hh \
src/KKT_like.hh \
src/TicToc.hh

MKDIR = mkdir -p

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX    = /usr/local
FRAMEWORK = Alglin

all: lib
	mkdir -p bin
	$(CXX) $(INC) $(CXXFLAGS) -o bin/AmodioN_test           src_tests/AmodioN_test.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/Amodio_test            src_tests/Amodio_test.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/BABD_test              src_tests/BABD_test.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/BABD_test2             src_tests/BABD_test2.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/Block_test             src_tests/Block_test.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/BorderedCR_test        src_tests/BorderedCR_test.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/CyclicReductionQR_test src_tests/CyclicReductionQR_test.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/CyclicReduction_test   src_tests/CyclicReduction_test.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/Diaz_test              src_tests/Diaz_test.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/KKT-test               src_tests/KKT-test.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/QR_N_test              src_tests/QR_N_test.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/QR_test                src_tests/QR_test.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/SMALL-test             src_tests/SMALL-test.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/TestFD                 src_tests/TestFD.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/Timing                 src_tests/Timing.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test1_Diaz             src_tests/test1_Diaz.cc $(LIBS)

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

install: lib/$(LIB_ALGLIN)
	cp src/*.hh          $(PREFIX)/include
	cp lib/$(LIB_ALGLIN) $(PREFIX)/lib

install_as_framework: lib/$(LIB_ALGLIN)
	$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	cp src/*.hh          $(PREFIX)/include/$(FRAMEWORK)
	cp lib/$(LIB_ALGLIN) $(PREFIX)/lib

run:
	./bin/AmodioN_test
	./bin/Amodio_test
	./bin/Block_test
	./bin/Diaz_test
	./bin/QR_N_test
	./bin/Timing
	./bin/CyclicReductionQR_test
	./bin/CyclicReduction_test
	./bin/QR_test
	./bin/TestFD
	./bin/KKT-test
	./bin/SMALL-test
	./bin/KKT-test
	./bin/BABD_test
	./bin/BABD_test2

doc:
	doxygen
	
clean:
	rm -f lib/libAlglin.* src/*.o
	rm -rf bin
	
