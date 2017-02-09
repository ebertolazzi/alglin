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
  CXXFLAGS = -Wall -O0 -g -fPIC -Wno-sign-compare -std=c++11 -pthread
  AR       = ar rcs
  LIBS     = -L./lib -lAlglin -llapack -lblas
  DEFS     = -DSUPERLU_OLD
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
src/BABD_Block.cc \
src/BABD_SuperLU.cc \
src/BlockBidiagonal.cc \
src/BorderedCR.cc \
src/KKT_like.cc \
src/Simplex.cc

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
src/Alglin_threads.hh \
src/Alglin_tmpl.hh \
src/BABD.hh \
src/BABD_Block.hh \
src/BABD_SuperLU.hh \
src/BlockBidiagonal.hh \
src/BorderedCR.hh \
src/KKT_like.hh \
src/Simplex.hh \
src/TicToc.hh

MKDIR = mkdir -p

# prefix for installation, use make PREFIX=/new/prefix install
# to override
PREFIX    = /usr/local
FRAMEWORK = Alglin

all: lib
	mkdir -p bin
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test0-FD                  src_tests/test0-FD.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test1-small-factorization src_tests/test1-small-factorization.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test2-Threads             src_tests/test2-Threads.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test3-Timing              src_tests/test3-Timing.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test4-KKT                 src_tests/test4-KKT.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test5-Diaz                src_tests/test5-Diaz.cc $(LIBS)
	#$(CXX) $(INC) $(CXXFLAGS) -o bin/test6-Block               src_tests/test6-Block.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/test7-BorderedCR          src_tests/test7-BorderedCR.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/SimplexTest1              src_tests/SimplexTest1.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/SimplexTest2              src_tests/SimplexTest2.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/SimplexTest3              src_tests/SimplexTest3.cc $(LIBS)
	$(CXX) $(INC) $(CXXFLAGS) -o bin/SimplexTest4              src_tests/SimplexTest4.cc $(LIBS)

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
	./bin/test0-FD
	./bin/test1-small-factorization
	./bin/test2-Threads
	./bin/test3-Timing
	./bin/test4-KKT
	./bin/test5-Diaz
	./bin/test7-BorderedCR
	./bin/SimplexTest1
	./bin/SimplexTest2
	./bin/SimplexTest3
	./bin/SimplexTest4

doc:
	doxygen

clean:
	rm -f lib/libAlglin.* src/*.o
	rm -rf bin
