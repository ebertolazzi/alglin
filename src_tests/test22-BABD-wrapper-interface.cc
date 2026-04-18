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
#include <string>
#include <vector>

static Utils::Console msg( &std::cout, 4 );

using real_type = double;
using integer   = alglin::integer;
using BABD      = alglin::BABD<real_type>;
using BCR       = alglin::BorderedCR<real_type>;

namespace {

  real_type constexpr TOL{ 1e-10 };

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

  struct TopBottomLayout {
    integer nblock{2};
    integer n{2};
    integer row0{2};
    integer col0{3};
    integer rowN{2};
    integer colN{3};
    integer nb{1};
    integer neq{ nblock * n + row0 + rowN };
    integer total{ neq + nb };
  };

  struct TopBottomData {
    std::vector<real_type> block0;
    std::vector<real_type> blockN;
    std::vector<real_type> blocks;
    std::vector<real_type> B;
    std::vector<real_type> C;
    std::vector<real_type> D;
  };

  TopBottomData
  make_top_bottom_data( TopBottomLayout const & layout ) {
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

    data.B.resize( static_cast<size_t>(layout.neq * layout.nb) );
    data.C.resize( static_cast<size_t>(layout.nb * layout.neq) );
    data.D.resize( static_cast<size_t>(layout.nb * layout.nb) );

    for ( integer j{0}; j < layout.nb; ++j ) {
      for ( integer i{0}; i < layout.neq; ++i ) {
        data.B[static_cast<size_t>( i + layout.neq * j )] = 0.03 * real_type(1 + i + 2 * j);
        data.C[static_cast<size_t>( j + layout.nb * i )]  = 0.02 * real_type(1 + 2 * i + j);
      }
    }

    data.D[0] = 2.4;
    return data;
  }

  void
  load_diaz_system( BABD & wrap, TopBottomLayout const & layout ) {
    TopBottomData const data{ make_top_bottom_data( layout ) };
    auto & solver = wrap.diaz();

    solver.allocate_top_bottom(
      layout.nblock,
      layout.n,
      layout.row0,
      layout.col0,
      layout.rowN,
      layout.colN,
      layout.nb
    );

    for ( integer k{0}; k < layout.nblock; ++k )
      solver.load_block(
        k,
        data.blocks.data() + 2 * k * layout.n * layout.n,
        layout.n
      );

    solver.load_right_blocks( data.B.data(), layout.neq );
    solver.load_bottom_blocks( data.C.data(), layout.nb );
    solver.load_RB_block( data.D.data(), layout.nb );
    solver.load_top_bottom( data.block0.data(), layout.row0, data.blockN.data(), layout.rowN );
    solver.select_last_block_solver_LU();
    solver.select_last_border_block_solver_LU();
  }

  integer constexpr N_BLOCKS { 6 };
  integer constexpr N_SIZE   { 2 };
  integer constexpr QR_SIZE  { 1 };
  integer constexpr QX_SIZE  { 2 };
  integer constexpr NR_SIZE  { 3 };
  integer constexpr NX_SIZE  { 2 };

  void
  fill_stable_cr_system( BCR & bcr ) {
    integer const nblock{ bcr.number_of_blocks() };
    integer const n     { bcr.block_size()       };
    integer const qr    { bcr.dim_qr()           };
    integer const qx    { bcr.dim_qx()           };
    integer const nr    { bcr.dim_nr()           };
    integer const nx    { bcr.dim_nx()           };
    integer const m     { n + qr                 };

    bcr.fill_zero();

    for ( integer k{0}; k < nblock; ++k ) {
      bcr.D(k,0,0) = 5.0 + 0.2*k;
      bcr.D(k,1,0) = 0.15;
      bcr.D(k,0,1) = -0.08;
      bcr.D(k,1,1) = 4.7 + 0.2*k;

      bcr.E(k,0,0) = -0.35;
      bcr.E(k,1,0) = 0.04;
      bcr.E(k,0,1) = 0.06;
      bcr.E(k,1,1) = -0.28;

      for ( integer j{0}; j < nx; ++j )
        for ( integer i{0}; i < n; ++i )
          bcr.B(k,i,j) = 0.04 * (1 + k + i + 2*j);

      for ( integer j{0}; j < n; ++j )
        for ( integer i{0}; i < nr; ++i )
          bcr.C(k,i,j) = 0.03 * (1 + 2*k + i + j);
    }

    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < nr; ++i )
        bcr.C(nblock,i,j) = 0.05 * (1 + i + j);

    for ( integer j{0}; j < qx; ++j )
      for ( integer i{0}; i < nr; ++i )
        bcr.Cq(i,j) = 0.07 * (1 + i + 2*j);

    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < nr; ++i )
        bcr.F(i,j) = 0.02 * (1 + i + j);

    bcr.F(0,0) += 3.2;
    bcr.F(1,1) += 2.8;

    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < m; ++i )
        bcr.H(i,j) = 0.02 * (1 + i + j);

    for ( integer j{0}; j < n; ++j )
      for ( integer i{0}; i < m; ++i )
        bcr.H(i,n+j) = 0.03 * (1 + i + j);

    bcr.H(0,n+0) += 5.1;
    bcr.H(1,n+1) += 4.9;
    bcr.H(2,n+0) += 0.2;
    bcr.H(2,n+1) -= 0.25;

    for ( integer j{0}; j < qx; ++j )
      for ( integer i{0}; i < m; ++i )
        bcr.H(i,2*n+j) = 0.04 * (1 + i + j);

    for ( integer j{0}; j < nx; ++j )
      for ( integer i{0}; i < m; ++i )
        bcr.H(i,2*n+qx+j) = 0.05 * (1 + i + j);
  }

  void
  load_cr_system( BABD & wrap ) {
    auto & solver = wrap.cr();
    solver.allocate( N_BLOCKS, N_SIZE, QR_SIZE, QX_SIZE, NR_SIZE, NX_SIZE );
    fill_stable_cr_system( solver );
    solver.select_last_LU();
    solver.select_last2_NONE();
  }

  void
  test_diaz_wrapper() {
    TopBottomLayout layout;
    BABD            wrap( BABD::BABD_Choice::DIAZ );

    UTILS_ASSERT0( wrap.using_diaz(), "BABD should start on Diaz backend\n" );
    UTILS_ASSERT0( !wrap.using_cyclic_reduction(), "BABD Diaz backend mismatch\n" );
    expect_throw( [&]() { (void) wrap.cr(); }, "BABD::cr on Diaz backend" );

    load_diaz_system( wrap, layout );

    UTILS_ASSERT0( wrap.nrows() == layout.total, "BABD Diaz nrows mismatch\n" );
    UTILS_ASSERT0( wrap.ncols() == layout.total, "BABD Diaz ncols mismatch\n" );
    UTILS_ASSERT0( wrap.info() == "Diaz", "BABD Diaz info mismatch\n" );
    UTILS_ASSERT0( wrap.sparse_nnz() > 0, "BABD Diaz sparse_nnz should be positive\n" );

    std::vector<real_type> xref( static_cast<size_t>(layout.total) );
    std::vector<real_type> rhs( static_cast<size_t>(layout.total) );
    std::vector<real_type> x;
    integer constexpr      nrhs{2};
    std::vector<real_type> Xref( static_cast<size_t>(layout.total * nrhs) );
    std::vector<real_type> RHS( static_cast<size_t>(layout.total * nrhs) );
    std::vector<real_type> X;

    for ( integer i{0}; i < layout.total; ++i )
      xref[static_cast<size_t>(i)] = 0.8 + 0.3 * i;

    wrap.diaz().Mv( xref.data(), rhs.data() );
    UTILS_ASSERT0( wrap.factorize_system(), "BABD Diaz factorize_system failed\n" );

    x = rhs;
    UTILS_ASSERT0( wrap.solve_system( x.data() ), "BABD Diaz solve_system failed\n" );
    for ( integer i{0}; i < layout.total; ++i )
      check_close( x[static_cast<size_t>(i)], xref[static_cast<size_t>(i)], "BABD Diaz solve_system", i, 5e-9 );

    load_diaz_system( wrap, layout );
    for ( integer c{0}; c < nrhs; ++c ) {
      for ( integer i{0}; i < layout.total; ++i )
        Xref[static_cast<size_t>( i + c * layout.total )] = (1 + c) * xref[static_cast<size_t>(i)] + 0.15 * c;
      wrap.diaz().Mv( Xref.data() + c * layout.total, RHS.data() + c * layout.total );
    }

    UTILS_ASSERT0( wrap.factorize(), "BABD Diaz factorize alias failed\n" );
    X = RHS;
    UTILS_ASSERT0( wrap.solve( nrhs, X.data(), layout.total ), "BABD Diaz solve alias failed\n" );
    for ( integer c{0}; c < nrhs; ++c )
      for ( integer i{0}; i < layout.total; ++i )
        check_close(
          X[static_cast<size_t>( i + c * layout.total )],
          Xref[static_cast<size_t>( i + c * layout.total )],
          "BABD Diaz solve(nrhs)",
          i, c, 5e-9
        );
  }

  void
  test_cr_wrapper() {
    BABD wrap;

    UTILS_ASSERT0( wrap.using_cyclic_reduction(), "BABD default backend should be CR\n" );
    UTILS_ASSERT0(
      wrap.selected_solver() == BABD::BABD_Choice::CYCLIC_REDUCTION_LU,
      "BABD default choice mismatch\n"
    );
    expect_throw( [&]() { (void) wrap.diaz(); }, "BABD::diaz on CR backend" );

    load_cr_system( wrap );
    UTILS_ASSERT0( wrap.nrows() == wrap.ncols(), "BABD CR should be square in this test\n" );
    UTILS_ASSERT0( wrap.sparse_nnz() > 0, "BABD CR sparse_nnz should be positive\n" );

    integer const N{ wrap.ncols() };
    std::vector<integer> I( static_cast<size_t>(wrap.sparse_nnz()) );
    std::vector<integer> J( static_cast<size_t>(wrap.sparse_nnz()) );
    std::vector<real_type> V( static_cast<size_t>(wrap.sparse_nnz()) );
    wrap.sparse_pattern( I.data(), J.data() );
    wrap.sparse_values( V.data() );
    UTILS_ASSERT0( I.front() >= 0 && J.front() >= 0, "BABD CR sparse pattern not filled\n" );

    std::vector<real_type> xref( static_cast<size_t>(N) );
    std::vector<real_type> rhs( static_cast<size_t>(N) );
    std::vector<real_type> x;

    for ( integer i{0}; i < N; ++i )
      xref[static_cast<size_t>(i)] = 1.0 + 0.2 * i;

    wrap.cr().Mv( xref.data(), rhs.data() );

    wrap.selectSolver( BABD::BABD_Choice::CYCLIC_REDUCTION_QRP );
    UTILS_ASSERT0(
      wrap.selected_solver() == BABD::BABD_Choice::CYCLIC_REDUCTION_QRP,
      "BABD CR solver switch mismatch\n"
    );
    UTILS_ASSERT0( wrap.using_cyclic_reduction(), "BABD should remain on CR backend\n" );
    UTILS_ASSERT(
      wrap.info().find("CyclicReduction+QRP") != std::string::npos,
      "BABD CR info should mention QRP\n"
    );

    UTILS_ASSERT0( wrap.factorize_system(), "BABD CR factorize_system failed\n" );
    x = rhs;
    UTILS_ASSERT0( wrap.solve_system( x.data() ), "BABD CR solve_system failed\n" );
    for ( integer i{0}; i < N; ++i )
      check_close( x[static_cast<size_t>(i)], xref[static_cast<size_t>(i)], "BABD CR solve_system", i, 5e-9 );
  }

} // namespace

int
main() {
  try {
    test_diaz_wrapper();
    test_cr_wrapper();
    msg.green( "All BABD wrapper interface checks passed!\n" );
  } catch ( std::exception const & err ) {
    msg.red( fmt::format( "Error: {}\n", err.what() ) );
    return 1;
  } catch ( ... ) {
    msg.red( "Unknown error\n" );
    return 1;
  }

  return 0;
}
