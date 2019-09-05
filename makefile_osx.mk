#######  #####  #     #
#     # #     #  #   #
#     # #         # #
#     #  #####     #
#     #       #   # #
#     # #     #  #   #
#######  #####  #     #

WARN     = -Weverything -Wno-reserved-id-macro -Wno-padded
CC       = clang -fPIC $(WARN)
CXX      = clang++ -fPIC -std=c++11 -stdlib=libc++ $(WARN)
CXXFLAGS = -O3
INC      = -Isrc -Ilib3rd/include
AR       = libtool -static -o
LIBS3RD  = -Llib3rd/lib -llapack_wrapper_osx_static -lsuperlu_osx_static
LIBS     = -L./lib/lib -lAlglin $(LIBS3RD)
LIBSGCC  = -lsuperlu_osx -lstdc++ -lm


ifneq (,$(findstring ALGLIN_USE_LAPACK,$(USED_LIB)))
  override LIBS3RD += -llapack -lblas
endif

ifneq (,$(findstring ALGLIN_USE_OPENBLAS,$(USED_LIB)))
  FPATH=$(dir $(shell gfortran -print-libgcc-file-name))
  override LIBS3RD += -Llib3rd/lib -Wl,-rpath,lib3rd/dll -lopenblas -L$(FPATH)/../../.. -Wl,-rpath,$(FPATH)/../../..  -lgfortran
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
  override LIBS3RD += -framework Accelerate
endif

LIBNAME = Alglin_osx

lib: config lib/lib/lib$(LIBNAME)_static.a lib/lib/lib$(LIBNAME).dylib

lib/lib/lib$(LIBNAME)_static.a: $(OBJS)
	$(MKDIR) ./lib
	$(MKDIR) ./lib/lib
	$(AR) lib/lib/lib$(LIBNAME)_static.a $(OBJS)

lib/lib/lib$(LIBNAME).dylib: $(OBJS)
	$(MKDIR) ./lib
	$(MKDIR) ./lib/lib
	$(CXX) -shared -o ./lib/lib/lib$(LIBNAME).dylib $(OBJS) $(LIBS3RD)

install_local: lib/lib/lib$(LIBNAME).a lib/lib/lib$(LIBNAME).dylib
	$(MKDIR) ./lib/include
	cp -f -P src/*.h*       ./lib/include
	cp -rf lib3rd/include/* ./lib/include

install: lib/lib/lib$(LIBNAME).a lib/lib/lib$(LIBNAME).dylib
	$(MKDIR) $(PREFIX)/include
	cp -f -P src/*.h*                    $(PREFIX)/include
	cp -f -P lib/lib/lib$(LIBNAME).a     $(PREFIX)/lib
	cp -f -P lib/lib/lib$(LIBNAME).dylib $(PREFIX)/lib

install_as_framework: lib/lib/lib$(LIBNAME).a lib/lib/lib$(LIBNAME).dylib
	$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	$(MKDIR) $(PREFIX)/lib
	cp -f -P src/*.h*                    $(PREFIX)/include/$(FRAMEWORK)
	cp -f -P lib/lib/lib$(LIBNAME).a     $(PREFIX)/lib
	cp -f -P lib/lib/lib$(LIBNAME).dylib $(PREFIX)/lib
