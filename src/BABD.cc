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

#if defined(__GCC__) || defined(__GNUC__) 
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wweak-template-vtables"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

#include "BABD.hh"

namespace alglin {

  string
  BABD_Choice_to_string( BABD_Choice c ) {
    string res = "none";
    switch ( c ) {
    case BABD_DIAZ:
      res = "Diaz";
      break;
    case BABD_CYCLIC_REDUCTION_LU:
      res = "CyclicReduction+LU";
      break;
    case BABD_CYCLIC_REDUCTION_QR:
      res = "CyclicReduction+QR";
      break;
    case BABD_CYCLIC_REDUCTION_QRP:
      res = "CyclicReduction+QRP";
      break;
    }
    return res;
  }

  #if 0

  static
  void
  dumpOneMatrix(
    ostream_type &  stream,
    char const      name[],
    valueType const M[],
    indexType       numRow,
    indexType       numCol
  ) {
    stream << "# " << name << " Size: "
           << numRow << " x " << numCol << '\n';
    stream << name << " := <";
    for ( indexType nc = 0; nc < numCol; ++nc ) {
      stream << '<';
      for ( indexType nr = 0; nr < numRow; ++nr ) {
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

    indexType N       = numNodes-1;
    indexType dim_z_z = 4*dim_x*dim_x;
    for ( indexType row = 0; row < N; ++row ) {
      valueType const * Ad = AdAu_blk + 2*row*dim_z_z;
      valueType const * Au = Ad + dim_z_z;
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
