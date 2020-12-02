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

  using namespace std;

  BABD_SuperLU::BABD_SuperLU()
  : m_baseValue("SuperLU_values")
  , m_baseInteger("SuperLU_integes")
  { }

  BABD_SuperLU::~BABD_SuperLU() {
    m_baseValue.free();
    m_baseInteger.free();
  }

  /*\
   |    __            _             _
   |   / _| __ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/
  void
  BABD_SuperLU::factorize(
    integer         _nblock,
    integer         _n,
    integer         q,
    valueType const AdAu[],
    valueType const H0[],
    valueType const HN[],
    valueType const Hq[]
  ) {

    m_nblock = _nblock;
    n        = _n;
    m        = _n+q;
    m_nnz    = 2*m_nblock*n*n + (n+m)*m;
    m_neq    = m_nblock*n + m;

    m_baseInteger.allocate( size_t(2*m_nnz+4*m_neq+1) );
    m_baseValue.allocate( size_t(m_nnz+2*m_neq ));

    set_default_options(&m_slu_options);

    // Initialize the statistics variables.
    StatInit(&m_slu_stats);

    m_perm_r = m_baseInteger( size_t(m_neq) ); /* row permutations from partial pivoting */
    m_perm_c = m_baseInteger( size_t(m_neq) ); /* column permutation vector */
    m_etree  = m_baseInteger( size_t(m_neq) );

    double * values = m_baseValue( size_t(m_nnz) );
    double * sumR   = m_baseValue( size_t(m_neq) );
    double * sumC   = m_baseValue( size_t(m_neq) );
    int    * rowind = m_baseInteger( size_t(m_nnz) );
    int    * colptr = m_baseInteger( size_t(m_neq+1) );

    integer kk = 0;
    integer jj = 0;
    integer ee = m_nblock*n;
    colptr[0] = 0;
    for ( integer k = 0; k <= m_nblock; ++k ) {
      valueType const * Ad = AdAu + 2*n*n*k;
      valueType const * Au = Ad - n*n;
      integer ii = k*n;
      for ( integer j = 0; j < n; ++j ) {
        if ( k > 0         ) for ( integer i = 0; i < n; ++i ) { rowind[kk] = i+ii-n; values[kk] = Au[i+j*n]; ++kk; }
        if ( k < m_nblock  ) for ( integer i = 0; i < n; ++i ) { rowind[kk] = i+ii; values[kk] = Ad[i+j*n]; ++kk; }
        if ( k == 0        ) for ( integer i = 0; i < m; ++i ) { rowind[kk] = i+ee; values[kk] = H0[i+j*m]; ++kk; }
        if ( k == m_nblock ) for ( integer i = 0; i < m; ++i ) { rowind[kk] = i+ee; values[kk] = HN[i+j*m]; ++kk; }
        colptr[++jj] = kk;
      }
    }
    for ( integer j = 0; j < q; ++j ) {
      for ( integer i = 0; i < m; ++i ) { rowind[kk] = i+ee; values[kk++] = Hq[i+j*m]; }
      colptr[++jj] = kk;
    }

    UTILS_ASSERT0( jj == m_neq, "SuperLU::factorize -- bad matrix format\n" );

    // Create matrix A in the format expected by SuperLU.
    //cout << "Create matrix A in the format expected by SuperLU.\n";
    dCreate_CompCol_Matrix(
      &m_A, m_neq, m_neq, m_nnz,
      values, rowind, colptr,
      SLU_NC, SLU_D, SLU_GE
    );

    std::fill( sumR, sumR+m_neq, 0 );
    std::fill( sumC, sumC+m_neq, 0 );
    for ( integer j = 0; j < m_neq; ++j ) {
      for ( integer ii = colptr[j]; ii < colptr[j+1]; ++ii ) {
        integer i = rowind[ii];
        double absAij = std::abs(values[ii]);
        sumR[i] += absAij;
        sumC[j] += absAij;
      }
    }
    m_one_norm_A = *std::max_element( sumC, sumC+m_neq );
    m_inf_norm_A = *std::max_element( sumR, sumR+m_neq );

    #if 0
    ofstream file("Jac.txt");
    file << "A := Matrix(" << neq << ", " << neq << ", fill = 0);\n";
    for ( integer j = 0; j < neq; ++j ) {
      for ( integer ii = colptr[j]; ii < colptr[j+1]; ++ii ) {
        integer i = rowind[ii];
        file << "A[" << i+1 << "," << j+1 << "] := " << values[ii] << ";\n";
      }
    }
    file.close();
    #endif

    factorize();
  }

  /*\
   |    __            _             _
   |   / _| __ _  ___| |_ ___  _ __(_)_______
   |  | |_ / _` |/ __| __/ _ \| '__| |_  / _ \
   |  |  _| (_| | (__| || (_) | |  | |/ /  __/
   |  |_|  \__,_|\___|\__\___/|_|  |_/___\___|
  \*/
  void
  BABD_SuperLU::factorize() {

    SuperMatrix AC;

    #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
    GlobalLU_t  glu;
    #endif

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   ColPerm = 0: natural ordering
     *   ColPerm = 1: minimum degree on structure of A'*A
     *   ColPerm = 2: minimum degree on structure of A'+A
     *   ColPerm = 3: approximate minimum degree for unsymmetric matrices
     */
    //cout << "get_perm_c.\n";
    get_perm_c( m_slu_options.ColPerm, &m_A, m_perm_c );
    //cout << "sp_preorder.\n";
    sp_preorder( &m_slu_options, &m_A, m_perm_c, m_etree, &AC);

    int panel_size = sp_ienv(1);
    int relax      = sp_ienv(2);
    int info = 0;
    //cout << "dgstrf.\n";
    dgstrf(
      &m_slu_options, &AC, relax, panel_size,
      m_etree, nullptr, 0, m_perm_c, m_perm_r, &m_L, &m_U,
    #if defined(SUPERLU_MAJOR_VERSION) && SUPERLU_MAJOR_VERSION >= 5
      &glu,
    #endif
      &m_slu_stats, &info
    );

    // Free un-wanted storage
    Destroy_SuperMatrix_Store(&m_A);
    Destroy_CompCol_Permuted(&AC);
    StatFree(&m_slu_stats);

    UTILS_ASSERT(
      info == 0,
      "BABD_SuperLU::factorize -- dgstrf() error returns INFO = {}\n", info
    );
    //cout << "done\n";
  }

  /*\
   |             _
   |   ___  ___ | |_   _____
   |  / __|/ _ \| \ \ / / _ \
   |  \__ \ (_) | |\ V /  __/
   |  |___/\___/|_| \_/ \___|
  \*/

  void
  BABD_SuperLU::solve( valueType y[] ) const {

    int         info;
    SuperMatrix B;
    trans_t     trans = NOTRANS;

    // Initialize the statistics variables.
    StatInit( &m_slu_stats );
    dCreate_Dense_Matrix(
      &B, m_L.nrow /* neq */,
      1 /* nrhs */,
      y,
      m_L.nrow /* ldy */,
      SLU_DN, SLU_D, SLU_GE
    );

    // Solve the system A*X=B, overwriting B with X.
    dgstrs( trans, &m_L, &m_U, m_perm_c, m_perm_r, &B, &m_slu_stats, &info);

    Destroy_SuperMatrix_Store(&B);
    StatFree(&m_slu_stats);

    UTILS_ASSERT(
      info == 0,
      "BABD_SuperLU::solve -- dgstrs() error returns INFO = {}\n", info
    );
  }

  /*\
   |                       _
   |    ___ ___  _ __   __| |
   |   / __/ _ \| '_ \ / _` |
   |  | (_| (_) | | | | (_| |
   |   \___\___/|_| |_|\__,_|
  \*/
  void
  BABD_SuperLU::cond( valueType & rcond_1, valueType & rcond_inf ) const {
    int info1 = 0, info2 = 0;
    StatInit( &m_slu_stats );
    dgscon(
      const_cast<char*>("1"),
      &m_L, &m_U, m_one_norm_A,
      &rcond_1, &m_slu_stats, &info1
    );
    if ( info1 == 0 )
      dgscon(
        const_cast<char*>("I"),
        &m_L, &m_U, m_inf_norm_A,
        &rcond_inf, &m_slu_stats, &info2
      );
    StatFree( &m_slu_stats );
    UTILS_ASSERT(
      info1 == 0,
      "BABD_SuperLU::cond -- dgscon() error returns INFO = {}\n", info1
    );
    UTILS_ASSERT(
      info2 == 0,
      "BABD_SuperLU::cond -- dgscon() error returns INFO = {}\n", info2
    );
  }

}
