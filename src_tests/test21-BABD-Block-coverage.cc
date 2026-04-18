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

#include <cmath>
#include <vector>

static Utils::Console msg( &std::cout, 4 );

using real_type = double;
using integer   = alglin::integer;

namespace {

  real_type constexpr TOL{ 1e-10 };

  struct CyclicLayout {
    integer nblock{2};
    integer n{2};
    integer q{1};
    integer neq{ (nblock + 1) * n + q };
  };

  struct CyclicData {
    std::vector<real_type> blocks;
    std::vector<real_type> H0;
    std::vector<real_type> HN;
    std::vector<real_type> Hq;
  };

  void
  check_close(
    real_type value,
    real_type expected,
    char const where[],
    integer i,
    integer j,
    real_type tol = TOL
  ) {
    UTILS_ASSERT(
      std::abs(value - expected) < tol,
      "{} mismatch at ({},{}): found {}, expected {}\n",
      where, i, j, value, expected
    );
  }

  void
  check_close(
    real_type value,
    real_type expected,
    char const where[],
    integer i,
    real_type tol = TOL
  ) {
    UTILS_ASSERT(
      std::abs(value - expected) < tol,
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

  CyclicData
  make_cyclic_data() {
    CyclicLayout layout;
    CyclicData   data;

    data.blocks.resize( static_cast<size_t>( 2 * layout.n * layout.n * layout.nblock ) );
    for ( integer k{0}; k < layout.nblock; ++k ) {
      real_type * blk{ data.blocks.data() + 2 * k * layout.n * layout.n };

      blk[0] = 4.5 + 0.4 * k;
      blk[1] = -0.2;
      blk[2] = 0.15;
      blk[3] = 4.0 + 0.3 * k;

      blk += layout.n * layout.n;
      blk[0] = 0.35;
      blk[1] = 0.08;
      blk[2] = -0.12;
      blk[3] = 0.25;
    }

    data.H0 = {
      1.0,  0.0,  0.4,
      0.0,  1.0, -0.2
    };
    data.HN = {
      0.3, -0.2, 1.0,
      0.1,  0.4, 0.7
    };
    data.Hq = { 1.2, -0.3, 2.0 };

    return data;
  }

  void
  load_cyclic_system(
    alglin::BBlockLU<real_type> & solver,
    CyclicLayout const &          layout,
    CyclicData const &            data
  ) {
    solver.allocate( layout.nblock, layout.n, 0, 0, 0, layout.n + layout.q, 0, 0, layout.q );
    solver.load_blocks( data.blocks.data(), layout.n );
    solver.load_bottom(
      data.H0.data(), layout.n + layout.q,
      data.HN.data(), layout.n + layout.q,
      data.Hq.data(), layout.n + layout.q
    );
  }

  std::vector<real_type>
  dense_from_sparse(
    alglin::BBlockLU<real_type> const & solver,
    integer                             neq
  ) {
    integer const nnz{ solver.sparse_nnz() };
    std::vector<integer>   I( static_cast<size_t>(nnz) );
    std::vector<integer>   J( static_cast<size_t>(nnz) );
    std::vector<real_type> V( static_cast<size_t>(nnz) );
    std::vector<real_type> A( static_cast<size_t>(neq * neq), 0 );

    solver.sparse_pattern( I.data(), J.data() );
    solver.sparse_values( V.data() );

    for ( integer k{0}; k < nnz; ++k ) {
      size_t const kk{ static_cast<size_t>(k) };
      A[static_cast<size_t>( I[kk] + neq * J[kk] )] += V[kk];
    }

    return A;
  }

  std::vector<real_type>
  dense_mv(
    std::vector<real_type> const & A,
    integer                        neq,
    real_type const                x[]
  ) {
    std::vector<real_type> y( static_cast<size_t>(neq), 0 );
    for ( integer j{0}; j < neq; ++j )
      for ( integer i{0}; i < neq; ++i )
        y[static_cast<size_t>(i)] += A[static_cast<size_t>( i + neq * j )] * x[j];
    return y;
  }

  void
  test_bblocklu_rejects_unsupported_layouts() {
    alglin::BBlockLU<real_type> solver;

    expect_throw(
      [&]() { solver.allocate_top_bottom( 2, 2, 2, 3, 2, 3, 0 ); },
      "BBlockLU::allocate_top_bottom"
    );

    expect_throw(
      [&]() { solver.allocate( 2, 2, 1, 0, 0, 3, 0, 0, 1 ); },
      "BBlockLU::allocate border"
    );

    expect_throw(
      [&]() { solver.allocate( 2, 2, 0, 1, 0, 2, 0, 0, 1 ); },
      "BBlockLU::allocate non-cyclic BC"
    );
  }

  void
  test_bblocklu_solve() {
    CyclicLayout                  layout;
    CyclicData const              data{ make_cyclic_data() };
    alglin::BBlockLU<real_type>   solver;
    std::vector<real_type>        xref( static_cast<size_t>(layout.neq) );
    std::vector<real_type>        rhs( static_cast<size_t>(layout.neq) );
    std::vector<real_type>        x;
    integer constexpr             nrhs{3};
    std::vector<real_type>        Xref( static_cast<size_t>(layout.neq * nrhs) );
    std::vector<real_type>        RHS( static_cast<size_t>(layout.neq * nrhs) );
    std::vector<real_type>        X;

    load_cyclic_system( solver, layout, data );

    for ( integer i{0}; i < layout.neq; ++i )
      xref[static_cast<size_t>(i)] = 0.8 + 0.35 * i;

    std::vector<real_type> const A{ dense_from_sparse( solver, layout.neq ) };

    solver.Mv( xref.data(), rhs.data() );
    std::vector<real_type> const rhs_expected{ dense_mv( A, layout.neq, xref.data() ) };
    for ( integer i{0}; i < layout.neq; ++i )
      check_close( rhs[static_cast<size_t>(i)], rhs_expected[static_cast<size_t>(i)], "BBlockLU::Mv", i );

    for ( integer c{0}; c < nrhs; ++c ) {
      for ( integer i{0}; i < layout.neq; ++i )
        Xref[static_cast<size_t>( i + c * layout.neq )] = (1 + c) * xref[static_cast<size_t>(i)] + 0.1 * c;
      solver.Mv( Xref.data() + c * layout.neq, RHS.data() + c * layout.neq );
    }

    solver.factorize();

    x = rhs;
    solver.solve( x.data() );
    for ( integer i{0}; i < layout.neq; ++i )
      check_close( x[static_cast<size_t>(i)], xref[static_cast<size_t>(i)], "BBlockLU::solve", i, 5e-9 );

    X = RHS;
    solver.solve( nrhs, X.data(), layout.neq );
    for ( integer c{0}; c < nrhs; ++c )
      for ( integer i{0}; i < layout.neq; ++i )
        check_close(
          X[static_cast<size_t>( i + c * layout.neq )],
          Xref[static_cast<size_t>( i + c * layout.neq )],
          "BBlockLU::solve(nrhs)",
          i, c, 5e-9
        );
  }

} // namespace

int
main() {
  try {
    test_bblocklu_rejects_unsupported_layouts();
    test_bblocklu_solve();
    msg.green( "All BBlockLU coverage checks passed!\n" );
  } catch ( std::exception const & err ) {
    msg.red( fmt::format( "Error: {}\n", err.what() ) );
    return 1;
  } catch ( ... ) {
    msg.red( "Unknown error\n" );
    return 1;
  }

  return 0;
}
