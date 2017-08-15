mkdir lib
mkdir lib/win32
mkdir lib/win32/dll
mkdir lib/win32/lib
mkdir lib/win32/include
mkdir lib/win32_int64
mkdir lib/win32_int64/dll
mkdir lib/win32_int64/lib
mkdir lib/win32_int64/include
mkdir lib/win64
mkdir lib/win64/dll
mkdir lib/win64/lib
mkdir lib/win64/include

#
#---------------------------------------------------------
#
DIR="openblas-win64-int64"
FILE="$DIR.zip"
DIROUT="OpenBLAS-v0.2.14-Win64-int64"
URL="https://sourceforge.net/projects/openblas/files/v0.2.14/OpenBLAS-v0.2.14-Win64-int64.zip/download"
if [ -f $FILE ];
then
  echo "$FILE already downloaded"
else
  curl -L $URL > $FILE
fi
unzip $FILE
cp -rf $DIROUT/bin/libopenblas.dll   lib/win64/dll/libopenblas_win64.dll
cp -rf $DIROUT/include/*             lib/win64/include
cp -rf $DIROUT/lib/libopenblas.a     lib/win64/lib/libopenblas_win64_static.lib
cp -rf $DIROUT/lib/libopenblas.dll.a lib/win64/lib/libopenblas_win64.lib
rm -rf $DIROUT

#
#---------------------------------------------------------
#
DIR="openblas-win64-int32"
FILE="$DIR.zip"
DIROUT="OpenBLAS-v0.2.14-Win64-int32"
URL="https://sourceforge.net/projects/openblas/files/v0.2.14/OpenBLAS-v0.2.14-Win64-int32.zip/download"
if [ -f $FILE ];
then
  echo "$FILE already downloaded"
else
  curl -L $URL > $FILE
fi
unzip $FILE
cp -rf $DIROUT/bin/libopenblas.dll   lib/win32_int64/dll/libopenblas_win32_int64.dll
cp -rf $DIROUT/include/*             lib/win32_int64/include
cp -rf $DIROUT/lib/libopenblas.a     lib/win32_int64/lib/libopenblas_win32_int64_static.lib
cp -rf $DIROUT/lib/libopenblas.dll.a lib/win32_int64/lib/libopenblas_win32_int64.lib
rm -rf $DIROUT

#
#---------------------------------------------------------
#
DIR="openblas-win32-int32"
FILE="$DIR.zip"
DIROUT="OpenBLAS-v0.2.9.rc2-x86_64-Win"
URL="https://sourceforge.net/projects/openblas/files/v0.2.9.rc2/OpenBLAS-v0.2.9.rc2-x86_64-Win.zip/download"
if [ -f $FILE ];
then
  echo "$FILE already downloaded"
else
  curl -L $URL > $FILE
fi
unzip $FILE
cp -rf $DIROUT/bin/libopenblas.dll   lib/win32/dll/libopenblas_win32.dll
cp -rf $DIROUT/include/*             lib/win32/include
cp -rf $DIROUT/lib/libopenblas.a     lib/win32/lib/libopenblas_win32_static.lib
cp -rf $DIROUT/lib/libopenblas.dll.a lib/win32/lib/libopenblas_win32.lib
rm -rf $DIROUT

