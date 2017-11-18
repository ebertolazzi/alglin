DIR="xianyi-OpenBLAS-6d2da63"
FILE="OpenBLAS.tar.gz"
URL="https://sourceforge.net/projects/openblas/files/v0.2.20/OpenBLAS%200.2.20%20version.tar.gz/download"
echo "$FILE"
if [ -f $FILE ];
then
  echo "$FILE already downloaded"
else
  curl -L $URL > $FILE
fi

if [ -f libs/lib/libopenblas.dylib ] || [ -f libs/lib/libopenblas.so ];
then
  echo "OpenBlas already compiled"
else
  rm -rf OpenBlas
  tar -zxvf $FILE
  mv $DIR OpenBlas
  cd OpenBlas ; make ; make install PREFIX=../libs ; cd ..
fi

if [ -f libs/lib/libopenblas.dylib ]; then
  LIB=libopenblas.dylib
fi

if [ -f libs/lib/libopenblas.so ]; then
  LIB=libopenblas.so
fi

if [ -x "$(command -v install_name_tool)" ]; then
  cd libs/lib ; install_name_tool -id @rpath/$LIB $LIB ; cd ../..
elif [ -x "$(command -v patchelf)" ]; then
  cd libs/lib ; patchelf --set-rpath @rpath/$LIB $LIB ; cd ../..
elif [ -x "$(command -v chrpath)" ]; then
  cd libs/lib ; chrpath --replace @rpath/$LIB $LIB ; cd ../..
fi

PREFIX="../../lib3rd"
if [ ! -d $PREFIX/include/openblas ];
then
  mkdir -p $PREFIX/include/openblas
fi
if [ ! -d $PREFIX/lib/openblas ];
then
  mkdir -p $PREFIX/lib/openblas
fi

cp -f -R -P libs/include/* $PREFIX/include/openblas
cp -f -R -P libs/lib/*     $PREFIX/lib/openblas
