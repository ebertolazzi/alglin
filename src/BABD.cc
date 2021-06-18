/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Alglin.hh"

namespace alglin {

  #if 0

  static
  void
  dumpOneMatrix(
    ostream_type &  stream,
    char const      name[],
    real_type const M[],
    integer         numRow,
    integer         numCol
  ) {
    stream << "# " << name << " Size: "
           << numRow << " x " << numCol << '\n';
    stream << name << " := <";
    for ( integer nc = 0; nc < numCol; ++nc ) {
      stream << '<';
      for ( integer nr = 0; nr < numRow; ++nr ) {
        stream << M[ nr + nc * numRow ];
        if ( nr+1 < numRow ) stream << ',';
        else                 stream << '>';
      }
      if ( nc+1 < numCol ) stream << "|\n";
      else                 stream << ">;\n";
    }
  }

  void
  BVNLFD_System::dump_to_Maple( ostream_type & stream ) const {

    stream << "interface( rtablesize = 40 );\n";

    integer N       = numNodes-1;
    integer dim_z_z = 4*dim_x*dim_x;
    for ( integer row = 0; row < N; ++row ) {
      real_type const * Ad = AdAu_blk + 2*row*dim_z_z;
      real_type const * Au = Ad + dim_z_z;
      dumpOneMatrix( stream, "Ad", Ad, 2*dim_x, 2*dim_x );
      dumpOneMatrix( stream, "Au", Au, 2*dim_x, 2*dim_x );
    }

    dumpOneMatrix( stream, "H0Np", H0Np, dim_bc, 4*dim_x + dim_omega );

    stream << "with(LinearAlgebra):\n";
    stream << "Determinant(Ad);\n";
    stream << "Determinant(Au);\n";
    stream << "Rank(H0Np);\n";
    stream << "Rank(<H0N|Hp>);\n";

  }
#endif
}
