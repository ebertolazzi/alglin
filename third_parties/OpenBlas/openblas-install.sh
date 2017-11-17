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

rm -rf OpenBlas
tar -zxvf $FILE
mv $DIR OpenBlas

cd OpenBlas ; make ; make install PREFIX=../libs ; cd ..

PREFIX="../../lib3rd"
mkdir -p $PREFIX/include/openblas
mkdir -p $PREFIX/lib

LIB=libopenblas.dylib
#LIB1=libopenblas_haswellp-r0.2.20.dylib

cd libs/lib ; install_name_tool -id @rpath/$LIB $LIB ; cd ../..

cp -f -R -P libs/include/* $PREFIX/include/openblas
cp -f -R -P libs/lib/*     $PREFIX/lib
