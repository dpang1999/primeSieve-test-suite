type :: T_TERM
    type(C_TYPE) :: coefficient
    type(E_TYPE) :: exponents
contains
    procedure :: compare => T_TERM_COMPARE
    procedure :: can_reduce => T_TERM_CAN_REDUCE
    procedure :: lcm => T_TERM_LCM
    procedure, private :: is_equal => T_TERM_IS_EQUAL
    generic :: operator(==) => is_equal
end type T_TERM

type :: T_POLY
    type(T_TERM), allocatable :: terms(:)
contains
    procedure :: addPoly => T_POLY_ADD
    procedure :: subtract => T_POLY_SUBTRACT
    procedure :: multiply_by_term => T_POLY_MUL_TERM
    procedure :: reduce_poly => T_POLY_REDUCE
    procedure :: make_monic => T_POLY_MAKE_MONIC
    procedure, private :: is_equal => T_POLY_IS_EQUAL
    generic :: operator(==) => is_equal
end type T_POLY

interface T_POLY
    module procedure T_INIT_POLY
    module procedure T_INIT_POLY_EMPTY
end interface
