/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2026                                                      |
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
 |      Universita' degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "Alglin.hh"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#endif

#include <algorithm>
#include <cmath>
#include <vector>

static Utils::Console msg( &std::cout, 4 );

using real_type = double;
using integer   = alglin::integer;

namespace {

  real_type constexpr TOL{ 1e-10 };

  struct TopBottomLayout {
    integer nblock{2};
    integer n{2};
    integer row0{2};
    integer col0{3};
    integer rowN{2};
    integer colN{3};
    integer neq{ nblock * n + row0 + rowN };
  };

  struct TopBottomData {
    std::vector<real_type> block0;
    std::vector<real_type> blockN;
    std::vector<real_type> blocks;
  };

  void
  check_close(
    real_type value,
    real_type expected,
    char const where[],
    integer i,
    real_type tol = TOL
  ) {
    UTILS_ASSERT(
      std::abs(value-expected) < tol,
      "{} mismatch at {}: found {}, expected {}\n",
      where, i, value, expected
    );
  }

  template <typename F>
  void
  expect_throw( F && fun, char const where[] ) {
    bool thrown{ false };
    try {
      fun();
    } catch ( ... ) {
      thrown = true;
    }
    UTILS_ASSERT( thrown, "{} expected an exception\n", where );
  }

  TopBottomData
  make_top_bottom_data() {
    TopBottomData data;

    data.block0 = {
      4.2, -0.3,
      0.1,  4.6,
      0.2, -0.1
    };

    data.blockN = {
      3.9,  0.15,
      0.05, 4.4,
      0.2, -0.25
    };

    TopBottomLayout const layout;
    data.blocks.resize( static_cast<size_t>(2 * layout.n * layout.n * layout.nblock) );
    for ( integer k{0}; k < layout.nblock; ++k ) {
      real_type * blk{ data.blocks.data() + 2 * k * layout.n * layout.n };
      blk[0] = 5.0 + 0.3 * k;
      blk[1] = -0.1;
      blk[2] = 0.2;
      blk[3] = 4.8 + 0.2 * k;

      blk += layout.n * layout.n;
      blk[0] = 0.3;
      blk[1] = 0.15;
      blk[2] = -0.1;
      blk[3] = 0.25;
    }

    return data;
  }

  std::vector<real_type>
  make_rhs_from_arceco_structure(
    TopBottomData const & data,
    std::vector<real_type> const & xref
  ) {
    TopBottomLayout       layout;
    std::vector<real_type> A( static_cast<size_t>(layout.neq * layout.neq), 0 );
    std::vector<real_type> rhs( static_cast<size_t>(layout.neq), 0 );

    auto put_block = [&]( integer row, integer col, integer rows, integer cols, real_type const src[], integer ldSrc ) {
      for ( integer j{0}; j < cols; ++j )
        for ( integer i{0}; i < rows; ++i )
          A[static_cast<size_t>( row + i + layout.neq * (col + j) )] = src[i + ldSrc * j];
    };

    integer row{0};
    integer col{0};
    put_block( row, col, layout.row0, layout.col0, data.block0.data(), layout.row0 );
    row += layout.row0;
    col += layout.col0 - layout.n;

    for ( integer k{0}; k < layout.nblock; ++k ) {
      real_type const * blk{ data.blocks.data() + 2 * k * layout.n * layout.n };
      put_block( row, col, layout.n, 2 * layout.n, blk, layout.n );
      row += layout.n;
      col += layout.n;
    }

    put_block( row, col, layout.rowN, layout.colN, data.blockN.data(), layout.rowN );

    for ( integer j{0}; j < layout.neq; ++j )
      for ( integer i{0}; i < layout.neq; ++i )
        rhs[static_cast<size_t>(i)] += A[static_cast<size_t>( i + layout.neq * j )] * xref[static_cast<size_t>(j)];

    return rhs;
  }

  void
  fill_structure_and_array(
    TopBottomLayout const &      layout,
    TopBottomData const &        data,
    std::vector<integer> &       matrix_structure,
    std::vector<real_type> &     array
  ) {
    matrix_structure.resize( static_cast<size_t>(3 * (layout.nblock + 2)) );
    array.resize(
      static_cast<size_t>(
        layout.row0 * layout.col0 +
        2 * layout.n * layout.n * layout.nblock +
        layout.rowN * layout.colN
      )
    );

    integer * mtr{ matrix_structure.data() };
    *mtr++ = layout.row0;
    *mtr++ = layout.col0;
    *mtr++ = layout.n;
    for ( integer i{0}; i < layout.nblock; ++i ) {
      *mtr++ = layout.n;
      *mtr++ = 2 * layout.n;
      *mtr++ = layout.n;
    }
    *mtr++ = layout.rowN;
    *mtr++ = layout.colN;
    *mtr++ = 0;

    real_type * ptr{ array.data() };
    ptr = std::copy_n( data.block0.data(), data.block0.size(), ptr );
    ptr = std::copy_n( data.blocks.data(),  data.blocks.size(),  ptr );
    std::copy_n( data.blockN.data(), data.blockN.size(), ptr );
  }

  void
  test_arceco_factorize_overload() {
    TopBottomLayout              layout;
    TopBottomData const          data{ make_top_bottom_data() };
    alglin::ArcecoLU<real_type>  arceco;
    std::vector<real_type>       xref( static_cast<size_t>(layout.neq) );
    std::vector<real_type>       rhs;
    std::vector<real_type>       x;

    for ( integer i{0}; i < layout.neq; ++i )
      xref[static_cast<size_t>(i)] = 1.0 + 0.4 * i;

    rhs = make_rhs_from_arceco_structure( data, xref );
    x   = rhs;

    arceco.factorize(
      layout.row0,
      layout.col0,
      data.block0.data(),
      layout.nblock,
      layout.n,
      data.blocks.data(),
      layout.rowN,
      layout.colN,
      data.blockN.data()
    );
    arceco.solve( x.data() );

    for ( integer i{0}; i < layout.neq; ++i )
      check_close( x[static_cast<size_t>(i)], xref[static_cast<size_t>(i)], "ArcecoLU::factorize overload", i, 5e-9 );
  }

  void
  test_arceco_load_by_ref_and_checkStructure() {
    TopBottomLayout              layout;
    TopBottomData const          data{ make_top_bottom_data() };
    alglin::ArcecoLU<real_type>  arceco;
    std::vector<integer>         matrix_structure;
    std::vector<real_type>       array;
    std::vector<integer>         pivot( static_cast<size_t>(layout.neq) );
    std::vector<real_type>       xref( static_cast<size_t>(layout.neq) );
    std::vector<real_type>       rhs;
    std::vector<real_type>       x;

    fill_structure_and_array( layout, data, matrix_structure, array );
    arceco.load_by_ref( layout.nblock + 2, matrix_structure.data(), array.data(), pivot.data() );
    arceco.checkStructure( layout.neq );

    for ( integer i{0}; i < layout.neq; ++i )
      xref[static_cast<size_t>(i)] = 0.5 + 0.35 * i;

    rhs = make_rhs_from_arceco_structure( data, xref );
    x   = rhs;

    arceco.factorize();
    arceco.solve( x.data() );

    for ( integer i{0}; i < layout.neq; ++i )
      check_close( x[static_cast<size_t>(i)], xref[static_cast<size_t>(i)], "ArcecoLU::load_by_ref", i, 5e-9 );
  }

  void
  test_arceco_invalid_structures() {
    TopBottomLayout              layout;
    TopBottomData const          data{ make_top_bottom_data() };
    alglin::ArcecoLU<real_type>  arceco;
    std::vector<integer>         matrix_structure;
    std::vector<real_type>       array;
    std::vector<integer>         pivot( static_cast<size_t>(layout.neq) );

    fill_structure_and_array( layout, data, matrix_structure, array );
    arceco.load_by_ref( layout.nblock + 2, matrix_structure.data(), array.data(), pivot.data() );
    expect_throw( [&]() { arceco.checkStructure( layout.neq + 1 ); }, "ArcecoLU::checkStructure neq mismatch" );

    fill_structure_and_array( layout, data, matrix_structure, array );
    matrix_structure.back() = 1;
    arceco.load_by_ref( layout.nblock + 2, matrix_structure.data(), array.data(), pivot.data() );
    expect_throw( [&]() { arceco.checkStructure( layout.neq ); }, "ArcecoLU::checkStructure last overlap" );
  }

} // namespace

int
main() {
  try {
    test_arceco_factorize_overload();
    test_arceco_load_by_ref_and_checkStructure();
    test_arceco_invalid_structures();
    msg.green( "All ArcecoLU coverage checks passed!\n" );
  } catch ( std::exception const & err ) {
    msg.red( fmt::format( "Error: {}\n", err.what() ) );
    return 1;
  } catch ( ... ) {
    msg.red( "Unknown error\n" );
    return 1;
  }

  return 0;
}
