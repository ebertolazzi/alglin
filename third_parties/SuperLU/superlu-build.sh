VER=5.2.1
DIR=superlu_$VER
DIR1=SuperLU_$VER
FILE="superlu_$VER.tar.gz"
URL="http://crd-legacy.lbl.gov/~xiaoye/SuperLU/$FILE"
if [ -f $FILE ];
then
  echo "$FILE already downloaded"
else
  curl $URL > $FILE
fi

if [ -f lib/lib/libsuperlu.a ];
then
  echo "SuperLU already compiled"
else
  rm -rf superlu
  gunzip $FILE -c | tar -xvf -

  # rinomina
  if [ -f $DIR ];
  then
    mv -f $DIR superlu
  else
    mv -f $DIR1 superlu
  fi

  sed 's/enable_language (Fortran)/\#enable_language (Fortran)/' superlu/CMakeLists.txt > tmp.txt
  mv -f tmp.txt superlu/CMakeLists.txt
  sed 's/set(NOFORTRAN FALSE)/set(NOFORTRAN TRUE)/' superlu/CMakeLists.txt > tmp.txt
  mv -f tmp.txt superlu/CMakeLists.txt
  sed 's/option(enable_tests     "Build tests" ON)/option(enable_tests     "Build tests" OFF)/' superlu/CMakeLists.txt > tmp.txt
  mv -f tmp.txt superlu/CMakeLists.txt
  sed 's/CMAKE_C_FLAGS \"/CMAKE_C_FLAGS "-fPIC /' superlu/CMakeLists.txt > tmp.txt
  mv -f tmp.txt superlu/CMakeLists.txt
  cd superlu ; cmake -DCMAKE_INSTALL_PREFIX:PATH=../lib . ; make ; make install ; cd ..
fi

PREFIX="../../lib3rd"

if [ ! -d $PREFIX/include/superlu ];
then
  mkdir -p $PREFIX/include/superlu
fi
if [ ! -d $PREFIX/lib/superlu  ];
then
  mkdir -p $PREFIX/lib/superlu
fi

cp -f -R -P lib/include/* $PREFIX/include/superlu
cp -f -R -P lib/lib/*     $PREFIX/lib/superlu
