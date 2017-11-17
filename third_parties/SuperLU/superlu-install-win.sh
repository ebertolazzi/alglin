#    ____                      _ _
#   / ___|___  _ __ ___  _ __ (_) | ___
#  | |   / _ \| '_ ` _ \| '_ \| | |/ _ \
#  | |__| (_) | | | | | | |_) | | |  __/
#   \____\___/|_| |_| |_| .__/|_|_|\___|
#                       |_|
#
PWD=`pwd`
ORIGIN="../superlu"
PREFIX="$PWD/../../LibSources/os_win"
YEAR=$1
BITS=$2

case $YEAR in
  "2010") str="Visual Studio 10 2010" ;;
  "2013") str="Visual Studio 12 2013" ;;
  "2015") str="Visual Studio 14 2015" ;;
   *)     str="UNKNOWN" ;;
esac

case $BITS in
  "x64") str="$str Win64" ;;
   *)      ;;
esac

DIR="vs${YEAR}_${BITS}"
rm -rf $DIR
mkdir $DIR
cd $DIR ; \
cmake -G "$str" $ORIGIN ; \
cmake --build . --config Release ; \
cmake --build . --config Debug ; \
cd ..

#   ___           _        _ _
#  |_ _|_ __  ___| |_ __ _| | |
#   | || '_ \/ __| __/ _` | | |
#   | || | | \__ \ || (_| | | |
#  |___|_| |_|___/\__\__,_|_|_|
#

LIBDIR="$PREFIX/lib"
INCDIR="$PREFIX/include/MechatronixCore"

# cmake di SUPERLU mette librerie sotto SRC!
cp $DIR/SRC/Debug/superlu.lib   $LIBDIR/superlu_vs${YEAR}_${BITS}_debug.lib
cp $DIR/SRC/Release/superlu.lib $LIBDIR/superlu_vs${YEAR}_${BITS}.lib

DIR="superlu/SRC"
mkdir -p $INCDIR
mkdir -p $INCDIR/superlu
cp -f $DIR/slu_*.h               $INCDIR/superlu
cp -f $DIR/colamd.h              $INCDIR/superlu
cp -f $DIR/superlu_enum_consts.h $INCDIR/superlu
cp -f $DIR/supermatrix.h         $INCDIR/superlu
