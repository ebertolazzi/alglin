VER=3.10.3
DIR=ATLAS
FILE="atlas$VER.tar.bz2"
URL="https://sourceforge.net/projects/math-atlas/files/Stable/$VER/atlas$VER.tar.bz2/download"
if [ -f $FILE ];
then
  echo "$FILE already downloaded"
else
  curl $URL > $FILE
fi

if [ -d ATLAS/build ];
then
  echo "ATLAS already compiled"
else
  rm -rf ATLAS
  bunzip2 $FILE -c | tar -xvf -
  mkdir ATLAS/build ; cd ATLAS/build ; ../configure ; cmake ; make ; cd ../..
fi

if [ -f libs/lib/libatlas.a ];
then
  echo "ATLAS already installed"
else
  cd ATLAS/build ; make DESTDIR=../../libs install ; cd ../..
fi

PREFIX="../../lib3rd"

if [ ! -d $PREFIX/include/atlas ];
then
  mkdir -p $PREFIX/include/atlas
fi
if [ ! -d $PREFIX/lib/atlas ];
then
  mkdir -p $PREFIX/lib/atlas
fi

cp -f -R -P libs/include/* $PREFIX/include/atlas
cp -f -R -P libs/lib/*     $PREFIX/lib/atlas

