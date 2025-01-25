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
    string_view     name,
    real_type const M[],
    integer         numRow,
    integer         numCol
  ) {
    fmt::print( stream,
      "# {} Size: {} x {}\n{} := <";
      name, numRow, numCol, name
    );

    for ( integer nc{0}; nc < numCol; ++nc ) {
      fmt::print( stream, "<" );
      for ( integer nr{0}; nr < numRow; ++nr ) {
        fmt::print( stream, "{}", M[ nr + nc * numRow ] );
        if ( nr+1 < numRow ) fmt::print( stream, "," );
        else                 fmt::print( stream, ">" );
      }
      if ( nc+1 < numCol ) fmt::print( stream, "|\n"  );
      else                 fmt::print( stream, ">;\n" );
    }
  }

  void
  BVNLFD_System::dump_to_Maple( ostream_type & stream ) const {

    fmt::print( stream, "interface( rtablesize = 40 );\n" );

    integer N{numNodes-1};
    integer dim_z_z{4*dim_x*dim_x};
    for ( integer row{0}; row < N; ++row ) {
      real_type const * Ad{AdAu_blk + 2*row*dim_z_z};
      real_type const * Au{Ad + dim_z_z};
      dumpOneMatrix( stream, "Ad", Ad, 2*dim_x, 2*dim_x );
      dumpOneMatrix( stream, "Au", Au, 2*dim_x, 2*dim_x );
    }

    dumpOneMatrix( stream, "H0Np", H0Np, dim_bc, 4*dim_x + dim_omega );

    fmt::print( stream,
      "with(LinearAlgebra):\n"
      "Determinant(Ad);\n"
      "Determinant(Au);\n"
      "Rank(H0Np);\n"
      "Rank(<H0N|Hp>);\n"
    );
  }
#endif
}
