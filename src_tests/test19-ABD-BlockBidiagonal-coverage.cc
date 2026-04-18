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

#include <array>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

static Utils::Console msg( &std::cout, 4 );

using real_type = double;
using integer   = alglin::integer;
using Choice    = alglin::BlockBidiagonal<real_type>::BB_LASTBLOCK_Choice;

namespace {

  real_type constexpr TOL{ 1e-10 };

  struct DummyBlockBidiagonal : public alglin::BlockBidiagonal<real_type> {
    using Base = alglin::BlockBidiagonal<real_type>;

    void
    allocate(
      integer nblock,
      integer n,
      integer nb,
      integer num_initial_BC,
      integer num_final_BC,
      integer num_cyclic_BC,
      integer num_initial_OMEGA,
      integer num_final_OMEGA,
      integer num_cyclic_OMEGA
    ) override {
      Base::allocate(
        nblock, n, nb,
        num_initial_BC, num_final_BC, num_cyclic_BC,
        num_initial_OMEGA, num_final_OMEGA, num_cyclic_OMEGA,
        0, 0
      );
    }

    void
    allocate_top_bottom(
      integer nblock,
      integer n,
      integer row0,
      integer col0,
      integer rowN,
      integer colN,
      integer nb
    ) override {
      Base::allocate_top_bottom(
        nblock, n, row0, col0, rowN, colN, nb, 0, 0
      );
    }

    void factorize() override {}
    void solve( real_type[] ) const override {}
    void solve( integer, real_type[], integer ) const override {}
  };

  struct GenericLayout {
    integer nblock{2};
    integer n{2};
    integer q{1};
    integer nb{2};
    integer neq{ (nblock+1) * n + q };
    integer total{ neq + nb };
  };

  struct TopBottomLayout {
    integer nblock{2};
    integer n{2};
    integer row0{2};
    integer col0{3};
    integer rowN{2};
    integer colN{3};
    integer nb{0};
    integer col00{ col0 - n };
    integer colNN{ colN - n };
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
      std::abs(value-expected) < tol,
      "{} mismatch at ({},{}): found {}, expected {}\n",
      where, i, j, value, expected
    );
  }

  void
  check_close(
    real_type value,
    real_type expected,
    char const where[],
    real_type tol = TOL
  ) {
    UTILS_ASSERT(
      std::abs(value-expected) < tol,
      "{} mismatch: found {}, expected {}\n",
      where, value, expected
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

  std::vector<real_type>
  make_matrix( integer rows, integer cols, integer tag ) {
    std::vector<real_type> M( static_cast<size_t>(rows * cols) );
    for ( integer j{0}; j < cols; ++j )
      for ( integer i{0}; i < rows; ++i )
        M[static_cast<size_t>( i + rows * j )] = real_type(tag) + 0.1 * real_type(10 * j + i + 1);
    return M;
  }

  void
  zero_block(
    std::vector<real_type> & A,
    integer                  ldA,
    integer                  row,
    integer                  col,
    integer                  rows,
    integer                  cols
  ) {
    for ( integer j{0}; j < cols; ++j )
      for ( integer i{0}; i < rows; ++i )
        A[static_cast<size_t>( row + i + ldA * (col + j) )] = 0;
  }

  void
  put_block(
    std::vector<real_type> & A,
    integer                  ldA,
    integer                  row,
    integer                  col,
    integer                  rows,
    integer                  cols,
    real_type const          src[],
    integer                  ldSrc
  ) {
    for ( integer j{0}; j < cols; ++j )
      for ( integer i{0}; i < rows; ++i )
        A[static_cast<size_t>( row + i + ldA * (col + j) )] = src[i + ldSrc * j];
  }

  void
  add_block(
    std::vector<real_type> & A,
    integer                  ldA,
    integer                  row,
    integer                  col,
    integer                  rows,
    integer                  cols,
    real_type const          src[],
    integer                  ldSrc
  ) {
    for ( integer j{0}; j < cols; ++j )
      for ( integer i{0}; i < rows; ++i )
        A[static_cast<size_t>( row + i + ldA * (col + j) )] += src[i + ldSrc * j];
  }

  std::vector<real_type>
  dense_from_sparse(
    alglin::BlockBidiagonal<real_type> const & solver,
    integer                                    total
  ) {
    integer const nnz{ solver.sparse_nnz() };
    std::vector<integer>   I( static_cast<size_t>(nnz) );
    std::vector<integer>   J( static_cast<size_t>(nnz) );
    std::vector<real_type> V( static_cast<size_t>(nnz) );
    std::vector<real_type> A( static_cast<size_t>(total * total), 0 );

    solver.sparse_pattern( I.data(), J.data() );
    solver.sparse_values( V.data() );

    for ( integer k{0}; k < nnz; ++k ) {
      size_t const kk{ static_cast<size_t>(k) };
      A[static_cast<size_t>( I[kk] + total * J[kk] )] += V[kk];
    }

    return A;
  }

  std::vector<real_type>
  dense_mv(
    std::vector<real_type> const & A,
    integer                        total,
    real_type const                x[]
  ) {
    std::vector<real_type> y( static_cast<size_t>(total), 0 );
    for ( integer j{0}; j < total; ++j )
      for ( integer i{0}; i < total; ++i )
        y[static_cast<size_t>(i)] += A[static_cast<size_t>( i + total * j )] * x[j];
    return y;
  }

  void
  check_dense_equal(
    std::vector<real_type> const & A,
    std::vector<real_type> const & B,
    integer                        total,
    char const                     where[]
  ) {
    UTILS_ASSERT( A.size() == B.size(), "{} size mismatch\n", where );
    for ( integer j{0}; j < total; ++j )
      for ( integer i{0}; i < total; ++i )
        check_close(
          A[static_cast<size_t>( i + total * j )],
          B[static_cast<size_t>( i + total * j )],
          where, i, j
        );
  }

  std::vector<real_type>
  rotate_right_copy( std::vector<real_type> v, integer amount ) {
    if ( amount > 0 ) {
      auto const mid{ v.end() - amount };
      std::rotate( v.begin(), mid, v.end() );
    }
    return v;
  }

  std::vector<real_type>
  rotate_columns_right_copy(
    std::vector<real_type> const & X,
    integer                        rows,
    integer                        cols,
    integer                        amount
  ) {
    std::vector<real_type> Y( X );
    for ( integer c{0}; c < cols; ++c ) {
      auto begin{ Y.begin() + c * rows };
      if ( amount > 0 ) std::rotate( begin, begin + rows - amount, begin + rows );
    }
    return Y;
  }

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

    if ( layout.nb > 0 ) {
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
      data.D[1] = -0.15;
      data.D[2] = 0.1;
      data.D[3] = 2.8;
    }

    return data;
  }

  template <typename Solver>
  void
  load_top_bottom_system( Solver & solver, TopBottomLayout const & layout ) {
    TopBottomData const data{ make_top_bottom_data( layout ) };

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

    if ( layout.nb > 0 ) {
      solver.load_right_blocks( data.B.data(), layout.neq );
      solver.load_bottom_blocks( data.C.data(), layout.nb );
      solver.load_RB_block( data.D.data(), layout.nb );
    }

    solver.load_top_bottom(
      data.block0.data(), layout.row0,
      data.blockN.data(), layout.rowN
    );
  }

  void
  test_blockbidiagonal_base_methods() {
    DummyBlockBidiagonal solver;
    GenericLayout const  layout;

    solver.allocate( layout.nblock, layout.n, layout.nb, 1, 1, 1, 0, 0, 1 );

    UTILS_ASSERT0(
      alglin::BlockBidiagonal<real_type>::LastBlock_to_string( Choice::PINV ) == "last block PINV",
      "LastBlock_to_string(PINV) failed\n"
    );

    std::vector<real_type> expected( static_cast<size_t>(layout.total * layout.total), 0 );

    std::vector<real_type> AdAu{ make_matrix( layout.n, 2 * layout.n * layout.nblock, 10 ) };
    solver.load_blocks( AdAu.data(), layout.n );
    for ( integer k{0}; k < layout.nblock; ++k ) {
      real_type const * blk{ AdAu.data() + 2 * k * layout.n * layout.n };
      put_block( expected, layout.total, k * layout.n, k * layout.n, layout.n, layout.n, blk, layout.n );
      put_block( expected, layout.total, k * layout.n, (k + 1) * layout.n, layout.n, layout.n, blk + layout.n * layout.n, layout.n );
    }
    check_dense_equal( dense_from_sparse( solver, layout.total ), expected, layout.total, "load_blocks" );

    std::vector<real_type> block1{ make_matrix( layout.n, 2 * layout.n, 20 ) };
    solver.load_block( 1, block1.data(), layout.n );
    put_block( expected, layout.total, layout.n, layout.n, layout.n, layout.n, block1.data(), layout.n );
    put_block( expected, layout.total, layout.n, 2 * layout.n, layout.n, layout.n, block1.data() + layout.n * layout.n, layout.n );

    std::vector<real_type> left0{ make_matrix( layout.n, layout.n, 30 ) };
    std::vector<real_type> right0{ make_matrix( layout.n, layout.n, 40 ) };
    solver.load_block_left( 0, left0.data(), layout.n );
    solver.load_block_right( 0, right0.data(), layout.n );
    put_block( expected, layout.total, 0, 0, layout.n, layout.n, left0.data(), layout.n );
    put_block( expected, layout.total, 0, layout.n, layout.n, layout.n, right0.data(), layout.n );

    std::vector<real_type> H0{ make_matrix( layout.n + layout.q, layout.n, 50 ) };
    std::vector<real_type> HN{ make_matrix( layout.n + layout.q, layout.n, 60 ) };
    std::vector<real_type> Hq{ make_matrix( layout.n + layout.q, layout.q, 70 ) };
    solver.load_bottom( H0.data(), layout.n + layout.q, HN.data(), layout.n + layout.q, Hq.data(), layout.n + layout.q );
    put_block( expected, layout.total, layout.nblock * layout.n, 0, layout.n + layout.q, layout.n, H0.data(), layout.n + layout.q );
    put_block( expected, layout.total, layout.nblock * layout.n, layout.nblock * layout.n, layout.n + layout.q, layout.n, HN.data(), layout.n + layout.q );
    put_block( expected, layout.total, layout.nblock * layout.n, layout.nblock * layout.n + layout.n, layout.n + layout.q, layout.q, Hq.data(), layout.n + layout.q );

    std::vector<real_type> LR( static_cast<size_t>(2 * layout.n * layout.n) );
    std::vector<real_type> L( static_cast<size_t>(layout.n * layout.n) );
    std::vector<real_type> R( static_cast<size_t>(layout.n * layout.n) );
    std::vector<real_type> outH0( static_cast<size_t>((layout.n + layout.q) * layout.n) );
    std::vector<real_type> outHN( static_cast<size_t>((layout.n + layout.q) * layout.n) );
    std::vector<real_type> outHq( static_cast<size_t>((layout.n + layout.q) * layout.q) );
    solver.getBlock_LR( LR.data(), layout.n );
    solver.getBlock_L( L.data(), layout.n );
    solver.getBlock_R( R.data(), layout.n );
    solver.getBlock_H0( outH0.data(), layout.n + layout.q );
    solver.getBlock_HN( outHN.data(), layout.n + layout.q );
    solver.getBlock_Hq( outHq.data(), layout.n + layout.q );

    for ( integer j{0}; j < layout.n; ++j ) {
      for ( integer i{0}; i < layout.n; ++i ) {
        check_close( L[static_cast<size_t>(i + layout.n * j)], left0[static_cast<size_t>(i + layout.n * j)], "getBlock_L", i, j );
        check_close( R[static_cast<size_t>(i + layout.n * j)], right0[static_cast<size_t>(i + layout.n * j)], "getBlock_R", i, j );
        check_close( LR[static_cast<size_t>(i + layout.n * j)], left0[static_cast<size_t>(i + layout.n * j)], "getBlock_LR/L", i, j );
        check_close( LR[static_cast<size_t>(i + layout.n * (layout.n + j))], right0[static_cast<size_t>(i + layout.n * j)], "getBlock_LR/R", i, j );
      }
    }
    for ( integer j{0}; j < layout.n; ++j )
      for ( integer i{0}; i < layout.n + layout.q; ++i ) {
        check_close( outH0[static_cast<size_t>(i + (layout.n + layout.q) * j)], H0[static_cast<size_t>(i + (layout.n + layout.q) * j)], "getBlock_H0", i, j );
        check_close( outHN[static_cast<size_t>(i + (layout.n + layout.q) * j)], HN[static_cast<size_t>(i + (layout.n + layout.q) * j)], "getBlock_HN", i, j );
      }
    for ( integer j{0}; j < layout.q; ++j )
      for ( integer i{0}; i < layout.n + layout.q; ++i )
        check_close( outHq[static_cast<size_t>(i + (layout.n + layout.q) * j)], Hq[static_cast<size_t>(i + (layout.n + layout.q) * j)], "getBlock_Hq", i, j );

    real_type const * ptrLR{ solver.getPointer_LR() };
    check_close( ptrLR[0], left0[0], "getPointer_LR[0]" );
    check_close( ptrLR[layout.n * layout.n], right0[0], "getPointer_LR[offset]" );

    std::vector<real_type> Bfull{ make_matrix( layout.neq, layout.nb, 80 ) };
    std::vector<real_type> Cfull{ make_matrix( layout.nb, layout.neq, 90 ) };
    std::vector<real_type> Dfull{ make_matrix( layout.nb, layout.nb, 100 ) };
    solver.load_right_blocks( Bfull.data(), layout.neq );
    solver.load_bottom_blocks( Cfull.data(), layout.nb );
    solver.load_RB_block( Dfull.data(), layout.nb );
    put_block( expected, layout.total, 0, layout.neq, layout.neq, layout.nb, Bfull.data(), layout.neq );
    put_block( expected, layout.total, layout.neq, 0, layout.nb, layout.neq, Cfull.data(), layout.nb );
    put_block( expected, layout.total, layout.neq, layout.neq, layout.nb, layout.nb, Dfull.data(), layout.nb );
    check_dense_equal( dense_from_sparse( solver, layout.total ), expected, layout.total, "full border loads" );

    solver.set_zero_right_blocks();
    solver.set_zero_bottom_blocks();
    solver.setZeroRBblock();
    zero_block( expected, layout.total, 0, layout.neq, layout.neq, layout.nb );
    zero_block( expected, layout.total, layout.neq, 0, layout.nb, layout.neq );
    zero_block( expected, layout.total, layout.neq, layout.neq, layout.nb, layout.nb );
    check_dense_equal( dense_from_sparse( solver, layout.total ), expected, layout.total, "zero border blocks" );

    std::vector<real_type> B0{ make_matrix( layout.n, layout.nb, 110 ) };
    std::vector<real_type> B1{ make_matrix( layout.n, layout.nb, 120 ) };
    std::vector<real_type> Blast{ make_matrix( layout.n + layout.q, layout.nb, 130 ) };
    solver.loadRightBlock( 0, B0.data(), layout.n );
    solver.loadRightBlock( 1, B1.data(), layout.n );
    solver.loadRightLastBlock( Blast.data(), layout.n + layout.q );
    put_block( expected, layout.total, 0, layout.neq, layout.n, layout.nb, B0.data(), layout.n );
    put_block( expected, layout.total, layout.n, layout.neq, layout.n, layout.nb, B1.data(), layout.n );
    put_block( expected, layout.total, layout.neq - (layout.n + layout.q), layout.neq, layout.n + layout.q, layout.nb, Blast.data(), layout.n + layout.q );

    std::vector<real_type> C0{ make_matrix( layout.nb, layout.n, 140 ) };
    std::vector<real_type> C1{ make_matrix( layout.nb, layout.n, 150 ) };
    std::vector<real_type> C1add{ make_matrix( layout.nb, layout.n, 160 ) };
    std::vector<real_type> C2add{ make_matrix( layout.nb, 2 * layout.n, 170 ) };
    std::vector<real_type> Clast{ make_matrix( layout.nb, layout.q, 180 ) };
    std::vector<real_type> Dlast{ make_matrix( layout.nb, layout.nb, 190 ) };

    solver.load_bottom_block( 0, C0.data(), layout.nb );
    solver.load_bottom_block( 1, C1.data(), layout.nb );
    solver.add_to_bottom_block( 1, C1add.data(), layout.nb );
    solver.add_to_bottom_block2( 0, C2add.data(), layout.nb );
    solver.load_bottom_last_block( Clast.data(), layout.nb );
    solver.load_RB_block( Dlast.data(), layout.nb );

    put_block( expected, layout.total, layout.neq, 0, layout.nb, layout.n, C0.data(), layout.nb );
    put_block( expected, layout.total, layout.neq, layout.n, layout.nb, layout.n, C1.data(), layout.nb );
    add_block( expected, layout.total, layout.neq, layout.n, layout.nb, layout.n, C1add.data(), layout.nb );
    add_block( expected, layout.total, layout.neq, 0, layout.nb, 2 * layout.n, C2add.data(), layout.nb );
    put_block( expected, layout.total, layout.neq, (layout.nblock + 1) * layout.n, layout.nb, layout.q, Clast.data(), layout.nb );
    put_block( expected, layout.total, layout.neq, layout.neq, layout.nb, layout.nb, Dlast.data(), layout.nb );

    check_dense_equal( dense_from_sparse( solver, layout.total ), expected, layout.total, "piecewise border loads" );

    std::vector<real_type> x( static_cast<size_t>(layout.total) );
    std::vector<real_type> y( static_cast<size_t>(layout.total) );
    for ( integer i{0}; i < layout.total; ++i ) x[static_cast<size_t>(i)] = 0.5 + 0.2 * i;
    solver.Mv( x.data(), y.data() );
    std::vector<real_type> const y_expected{ dense_mv( expected, layout.total, x.data() ) };
    for ( integer i{0}; i < layout.total; ++i )
      check_close( y[static_cast<size_t>(i)], y_expected[static_cast<size_t>(i)], "DummyBlockBidiagonal::Mv", i, 0 );

    std::ostringstream maple;
    solver.dump_to_Maple( maple );
    UTILS_ASSERT0( maple.str().find("interface(") != std::string::npos, "dump_to_Maple missing interface()\n" );

    std::ostringstream ccoord;
    solver.dump_ccoord( ccoord );
    std::istringstream is( ccoord.str() );
    integer dumped_nnz{ -1 };
    is >> dumped_nnz;
    UTILS_ASSERT0( dumped_nnz == solver.sparse_nnz(), "dump_ccoord wrong nnz header\n" );
  }

  template <typename Solver>
  void
  test_solver_no_border( char const where[] ) {
    TopBottomLayout layout;
    layout.nb    = 0;
    layout.total = layout.neq;

    constexpr Choice choices[]{
      Choice::LU, Choice::LUPQ, Choice::QR, Choice::QRP,
      Choice::SVD, Choice::LSS, Choice::LSY, Choice::PINV
    };

    for ( Choice choice : choices ) {
      Solver solver;
      std::string const choice_label{
        fmt::format( "{} [{}]", where, alglin::BlockBidiagonal<real_type>::LastBlock_to_string(choice) )
      };
      load_top_bottom_system( solver, layout );
      solver.select_last_block_solver( choice );

      std::vector<real_type> xref( static_cast<size_t>(layout.neq) );
      std::vector<real_type> rhs( static_cast<size_t>(layout.neq) );
      std::vector<real_type> x( static_cast<size_t>(layout.neq) );
      std::vector<real_type> x_abd_expected;
      std::vector<real_type> x_abd;

      for ( integer i{0}; i < layout.neq; ++i )
        xref[static_cast<size_t>(i)] = 1.0 + 0.25 * i;

      integer constexpr nrhs{3};
      std::vector<real_type> Xref( static_cast<size_t>(layout.neq * nrhs) );
      std::vector<real_type> RHS( static_cast<size_t>(layout.neq * nrhs) );
      std::vector<real_type> X( static_cast<size_t>(layout.neq * nrhs) );
      std::vector<real_type> XabdExpected;
      std::vector<real_type> Xabd;

      for ( integer c{0}; c < nrhs; ++c ) {
        for ( integer i{0}; i < layout.neq; ++i )
          Xref[static_cast<size_t>( i + c * layout.neq )] = (1 + c) * xref[static_cast<size_t>(i)] - 0.2 * c;
        solver.Mv( Xref.data() + c * layout.neq, RHS.data() + c * layout.neq );
      }

      solver.Mv( xref.data(), rhs.data() );
      solver.factorize();

      x = rhs;
      solver.solve( x.data() );
      for ( integer i{0}; i < layout.neq; ++i )
        check_close( x[static_cast<size_t>(i)], xref[static_cast<size_t>(i)], choice_label.c_str(), i, 0, 5e-9 );

      x_abd_expected = rotate_right_copy( xref, layout.col00 );
      x_abd          = rotate_right_copy( rhs, layout.row0 );
      solver.solve_ABD( x_abd.data() );
      for ( integer i{0}; i < layout.neq; ++i )
        check_close( x_abd[static_cast<size_t>(i)], x_abd_expected[static_cast<size_t>(i)], fmt::format("{} solve_ABD", choice_label).c_str(), i, 0, 5e-9 );

      X = RHS;
      solver.solve( nrhs, X.data(), layout.neq );
      for ( integer c{0}; c < nrhs; ++c )
        for ( integer i{0}; i < layout.neq; ++i )
          check_close(
            X[static_cast<size_t>( i + c * layout.neq )],
            Xref[static_cast<size_t>( i + c * layout.neq )],
            fmt::format("{} solve(nrhs)", choice_label).c_str(), i, c, 5e-9
          );

      XabdExpected = rotate_columns_right_copy( Xref, layout.neq, nrhs, layout.col00 );
      Xabd         = rotate_columns_right_copy( RHS,  layout.neq, nrhs, layout.row0 );
      solver.solve_ABD( nrhs, Xabd.data(), layout.neq );
      for ( integer c{0}; c < nrhs; ++c )
        for ( integer i{0}; i < layout.neq; ++i )
          check_close(
            Xabd[static_cast<size_t>( i + c * layout.neq )],
            XabdExpected[static_cast<size_t>( i + c * layout.neq )],
            fmt::format("{} solve_ABD(nrhs)", choice_label).c_str(), i, c, 5e-9
          );
    }
  }

  template <typename Solver>
  void
  test_solver_bordered( char const where[] ) {
    TopBottomLayout layout;
    layout.nb    = 2;
    layout.total = layout.neq + layout.nb;

    constexpr Choice choices[]{
      Choice::LU, Choice::LUPQ, Choice::QR, Choice::QRP,
      Choice::SVD, Choice::LSS, Choice::LSY, Choice::PINV
    };

    integer constexpr nrhs{3};

    for ( Choice last_choice : choices ) {
      for ( Choice border_choice : choices ) {
        Solver solver;
        std::string const choice_label{
          fmt::format(
            "{} [{} | border {}]",
            where,
            alglin::BlockBidiagonal<real_type>::LastBlock_to_string(last_choice),
            alglin::BlockBidiagonal<real_type>::LastBlock_to_string(border_choice)
          )
        };
        load_top_bottom_system( solver, layout );
        solver.select_last_block_solver( last_choice );
        solver.select_last_border_block_solver( border_choice );

        std::vector<real_type> xref( static_cast<size_t>(layout.total) );
        std::vector<real_type> rhs( static_cast<size_t>(layout.total) );
        std::vector<real_type> x( static_cast<size_t>(layout.total) );

        for ( integer i{0}; i < layout.total; ++i )
          xref[static_cast<size_t>(i)] = 0.75 + 0.3 * i;

        std::vector<real_type> Xref( static_cast<size_t>(layout.total * nrhs) );
        std::vector<real_type> RHS( static_cast<size_t>(layout.total * nrhs) );
        std::vector<real_type> X( static_cast<size_t>(layout.total * nrhs) );

        for ( integer c{0}; c < nrhs; ++c ) {
          for ( integer i{0}; i < layout.total; ++i )
            Xref[static_cast<size_t>( i + c * layout.total )] = (1 + c) * xref[static_cast<size_t>(i)] + 0.1 * c;
          solver.Mv( Xref.data() + c * layout.total, RHS.data() + c * layout.total );
        }

        solver.Mv( xref.data(), rhs.data() );
        solver.factorize_bordered();

        x = rhs;
        solver.solve_bordered( x.data() );
        for ( integer i{0}; i < layout.total; ++i )
          check_close( x[static_cast<size_t>(i)], xref[static_cast<size_t>(i)], choice_label.c_str(), i, 0, 5e-9 );

        X = RHS;
        solver.solve_bordered( nrhs, X.data(), layout.total );
        for ( integer c{0}; c < nrhs; ++c )
          for ( integer i{0}; i < layout.total; ++i )
            check_close(
              X[static_cast<size_t>( i + c * layout.total )],
              Xref[static_cast<size_t>( i + c * layout.total )],
              fmt::format("{} solve_bordered(nrhs)", choice_label).c_str(), i, c, 5e-9
            );
      }
    }
  }

} // namespace

int
main() {
  try {
    test_blockbidiagonal_base_methods();
    test_solver_no_border<alglin::BlockLU<real_type>>( "BlockLU" );
    test_solver_no_border<alglin::DiazLU<real_type>>( "DiazLU" );
    test_solver_bordered<alglin::BlockLU<real_type>>( "BlockLU bordered" );
    test_solver_bordered<alglin::DiazLU<real_type>>( "DiazLU bordered" );
    msg.green( "All ABD/BlockBidiagonal coverage checks passed!\n" );
  } catch ( std::exception const & err ) {
    msg.red( fmt::format( "Error: {}\n", err.what() ) );
    return 1;
  } catch ( ... ) {
    msg.red( "Unknown error\n" );
    return 1;
  }

  return 0;
}
