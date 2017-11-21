DIR="xianyi-OpenBLAS-6d2da63"
FILE="OpenBLAS.tar.gz"


URL=http://cs2.swfu.edu.cn/~zyl/lapack/LAPACK_Debug_x86.zip
URL=http://cs2.swfu.edu.cn/~zyl/lapack/LAPACK_Release_x86.zip
URL=http://cs2.swfu.edu.cn/~zyl/lapack/LAPACK_Debug_x64.zip
URL=http://cs2.swfu.edu.cn/~zyl/lapack/LAPACK_Release_x64.zip

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

cd libs/lib ; install_name_tool -id @rpath/$LIB $LIB ; cd ../..

cp -f -R -P libs/include/* $PREFIX/include/openblas
cp -f -R -P libs/lib/*     $PREFIX/lib
