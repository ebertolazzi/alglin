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
CXXFLAGS = -O2
INC      = -Isrc -Ilib3rd/include
AR       = libtool -static -o
LIBS3RD  = -Llib3rd/lib -Llib3rd/dll  -Wl,-rpath,lib3rd/dll -llapack_wrapper_osx_static -lsuperlu_osx_static
LIBS     = -Llib/lib -Llib/dll -lfmt_osx_static -lfort_osx_static -lAlglin_osx $(LIBS3RD) -Wl,-rpath,lib/dll
LIBSGCC  = -lstdc++ -lm


ifneq (,$(findstring ALGLIN_USE_LAPACK,$(USED_LIB)))
  override LIBS3RD += -llapack -lblas
endif

ifneq (,$(findstring ALGLIN_USE_OPENBLAS,$(USED_LIB)))
  FPATH=$(dir $(shell gfortran -print-libgcc-file-name))
  override LIBS3RD += -Llib3rd/lib -lopenblas -L$(FPATH)/../../.. -Wl,-rpath,$(FPATH)/../../..  -lgfortran
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

lib: lib/lib/lib$(LIBNAME)_static.a lib/dll/lib$(LIBNAME).dylib

lib/lib/lib$(LIBNAME)_static.a: $(OBJS)
	$(MKDIR) ./lib
	$(MKDIR) ./lib/lib
	$(AR) lib/lib/lib$(LIBNAME)_static.a $(OBJS)

lib/dll/lib$(LIBNAME).dylib: $(OBJS)
	$(MKDIR) ./lib
	$(MKDIR) ./lib/dll
	$(CXX) -shared -o ./lib/dll/lib$(LIBNAME).dylib $(OBJS) $(LIBS3RD)

install_local: lib/lib/lib$(LIBNAME).a lib/dll/lib$(LIBNAME).dylib
	$(MKDIR) ./lib/include
	cp -f -P src/*.h*       ./lib/include
	cp -rf lib3rd/include/* ./lib/include

install: lib/lib/lib$(LIBNAME).a lib/dll/lib$(LIBNAME).dylib
	$(MKDIR) $(PREFIX)/include
	cp -f -P src/*.h*                    $(PREFIX)/include
	cp -f -P lib/lib/lib$(LIBNAME).a     $(PREFIX)/lib
	cp -f -P lib/dll/lib$(LIBNAME).dylib $(PREFIX)/dll

install_as_framework: lib/lib/lib$(LIBNAME).a lib/dll/lib$(LIBNAME).dylib
	$(MKDIR) $(PREFIX)/include/$(FRAMEWORK)
	$(MKDIR) $(PREFIX)/lib
	cp -f -P src/*.h*                    $(PREFIX)/include/$(FRAMEWORK)
	cp -f -P lib/lib/lib$(LIBNAME).a     $(PREFIX)/lib
	cp -f -P lib/dll/lib$(LIBNAME).dylib $(PREFIX)/dll
