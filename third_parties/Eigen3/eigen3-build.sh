VER=3.3.4
DIR=$VER
FILE="$DIR.tar.gz"
URL="http://bitbucket.org/eigen/eigen/get/$DIR.tar.gz"
if [ -f $FILE ];
then
  echo "$FILE already downloaded"
else
  curl -L $URL > $FILE
fi

if [ ! -d eigen3 ]; then
  tar -zxvf $FILE
  rm -rf eigen3
  mv -f eigen-eigen* eigen3
fi

PREFIX="../../lib3rd"
if [ ! -d $PREFIX/include ];
then
  mkdir -p $PREFIX/include
fi

cp -r eigen3/Eigen $PREFIX/include
