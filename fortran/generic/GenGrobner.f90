module GenGrobner_Module
    use IntModP_Module
    use VecExponent_Module
    use BitPackedExponent_Module
    implicit none
    private

    integer, public :: TERM_ORDER = 0 ! 0 = Lex, 1 = GrLex, 2 = RevLex

    ! ====================================================
    ! Types Instantiation
    ! ====================================================
    public :: term_vec, polynomial_vec, naive_grobner_basis_vec
    public :: term_bit, polynomial_bit, naive_grobner_basis_bit

! Vec Types
#define C_TYPE IntModP
#define E_TYPE VecExponent
#define T_TERM term_vec
#define T_POLY polynomial_vec
#define T_TERM_COMPARE term_compare_vec
#define T_TERM_CAN_REDUCE term_can_reduce_vec
#define T_TERM_LCM term_lcm_vec
#define T_TERM_IS_EQUAL term_is_equal_vec
#define T_POLY_ADD poly_add_vec
#define T_POLY_SUBTRACT poly_subtract_vec
#define T_POLY_MAKE_MONIC poly_make_monic_vec
#define T_POLY_MUL_TERM poly_mul_term_vec
#define T_POLY_REDUCE poly_reduce_vec
#define T_POLY_IS_EQUAL poly_is_equal_vec
#define T_INIT_POLY init_poly_vec
#define T_INIT_POLY_EMPTY init_poly_empty_vec
#include "GenGrobner_Types.inc"
#undef C_TYPE
#undef E_TYPE
#undef T_TERM
#undef T_POLY
#undef T_TERM_COMPARE
#undef T_TERM_CAN_REDUCE
#undef T_TERM_LCM
#undef T_TERM_IS_EQUAL
#undef T_POLY_ADD
#undef T_POLY_SUBTRACT
#undef T_POLY_MAKE_MONIC
#undef T_POLY_MUL_TERM
#undef T_POLY_REDUCE
#undef T_POLY_IS_EQUAL
#undef T_INIT_POLY
#undef T_INIT_POLY_EMPTY

! Bit Types
#define C_TYPE IntModP
#define E_TYPE BitPackedExponent
#define T_TERM term_bit
#define T_POLY polynomial_bit
#define T_TERM_COMPARE term_compare_bit
#define T_TERM_CAN_REDUCE term_can_reduce_bit
#define T_TERM_LCM term_lcm_bit
#define T_TERM_IS_EQUAL term_is_equal_bit
#define T_POLY_ADD poly_add_bit
#define T_POLY_SUBTRACT poly_subtract_bit
#define T_POLY_MAKE_MONIC poly_make_monic_bit
#define T_POLY_MUL_TERM poly_mul_term_bit
#define T_POLY_REDUCE poly_reduce_bit
#define T_POLY_IS_EQUAL poly_is_equal_bit
#define T_INIT_POLY init_poly_bit
#define T_INIT_POLY_EMPTY init_poly_empty_bit
#include "GenGrobner_Types.inc"
#undef C_TYPE
#undef E_TYPE
#undef T_TERM
#undef T_POLY
#undef T_TERM_COMPARE
#undef T_TERM_CAN_REDUCE
#undef T_TERM_LCM
#undef T_TERM_IS_EQUAL
#undef T_POLY_ADD
#undef T_POLY_SUBTRACT
#undef T_POLY_MAKE_MONIC
#undef T_POLY_MUL_TERM
#undef T_POLY_REDUCE
#undef T_POLY_IS_EQUAL
#undef T_INIT_POLY
#undef T_INIT_POLY_EMPTY


contains

! ====================================================
! Procedures Instantiation
! ====================================================

! Vec Procedures
#define C_TYPE IntModP
#define E_TYPE VecExponent
#define T_TERM term_vec
#define T_POLY polynomial_vec
#define T_TERM_COMPARE term_compare_vec
#define T_TERM_CAN_REDUCE term_can_reduce_vec
#define T_TERM_LCM term_lcm_vec
#define T_TERM_IS_EQUAL term_is_equal_vec
#define T_POLY_ADD poly_add_vec
#define T_POLY_SUBTRACT poly_subtract_vec
#define T_POLY_MAKE_MONIC poly_make_monic_vec
#define T_POLY_MUL_TERM poly_mul_term_vec
#define T_POLY_REDUCE poly_reduce_vec
#define T_POLY_IS_EQUAL poly_is_equal_vec
#define T_INIT_POLY init_poly_vec
#define T_INIT_POLY_EMPTY init_poly_empty_vec
#define T_SORT_TERMS sort_terms_vec
#define T_S_POLY s_poly_vec
#define T_GROBNER naive_grobner_basis_vec
#include "GenGrobner_Procedures.inc"
#undef C_TYPE
#undef E_TYPE
#undef T_TERM
#undef T_POLY
#undef T_TERM_COMPARE
#undef T_TERM_CAN_REDUCE
#undef T_TERM_LCM
#undef T_TERM_IS_EQUAL
#undef T_POLY_ADD
#undef T_POLY_SUBTRACT
#undef T_POLY_MAKE_MONIC
#undef T_POLY_MUL_TERM
#undef T_POLY_REDUCE
#undef T_POLY_IS_EQUAL
#undef T_INIT_POLY
#undef T_INIT_POLY_EMPTY
#undef T_SORT_TERMS
#undef T_S_POLY
#undef T_GROBNER

! Bit Procedures
#define C_TYPE IntModP
#define E_TYPE BitPackedExponent
#define T_TERM term_bit
#define T_POLY polynomial_bit
#define T_TERM_COMPARE term_compare_bit
#define T_TERM_CAN_REDUCE term_can_reduce_bit
#define T_TERM_LCM term_lcm_bit
#define T_TERM_IS_EQUAL term_is_equal_bit
#define T_POLY_ADD poly_add_bit
#define T_POLY_SUBTRACT poly_subtract_bit
#define T_POLY_MAKE_MONIC poly_make_monic_bit
#define T_POLY_MUL_TERM poly_mul_term_bit
#define T_POLY_REDUCE poly_reduce_bit
#define T_POLY_IS_EQUAL poly_is_equal_bit
#define T_INIT_POLY init_poly_bit
#define T_INIT_POLY_EMPTY init_poly_empty_bit
#define T_SORT_TERMS sort_terms_bit
#define T_S_POLY s_poly_bit
#define T_GROBNER naive_grobner_basis_bit
#include "GenGrobner_Procedures.inc"
#undef C_TYPE
#undef E_TYPE
#undef T_TERM
#undef T_POLY
#undef T_TERM_COMPARE
#undef T_TERM_CAN_REDUCE
#undef T_TERM_LCM
#undef T_TERM_IS_EQUAL
#undef T_POLY_ADD
#undef T_POLY_SUBTRACT
#undef T_POLY_MAKE_MONIC
#undef T_POLY_MUL_TERM
#undef T_POLY_REDUCE
#undef T_POLY_IS_EQUAL
#undef T_INIT_POLY
#undef T_INIT_POLY_EMPTY
#undef T_SORT_TERMS
#undef T_S_POLY
#undef T_GROBNER

end module GenGrobner_Module

