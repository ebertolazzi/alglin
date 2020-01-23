#
# select which version of BLAS/LAPACK use
#
USED_LIB=""
SYSTEMOPENBLAS=\/\/
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
  USED_LIB = ALGLIN_USE_OPENBLAS
endif
else
  USED_LIB = $(shell cat .alglin_config)
endif
endif

$(shell echo "$(USED_LIB)" > .alglin_config )
$(info $(USED_LIB))
