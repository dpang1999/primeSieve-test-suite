module IntModP_Module
    implicit none
    private

    public :: IntModP, mod_inverse

    ! Note: Modulus is currently a global parameter, but we could use a PDT
    ! like type IntModP(modulus) but keeping it simple for now as modulus=7
    ! is generally what is used, or a module variable.
    integer(8), save, public :: IntModP_modulus = 7

    type :: IntModP
        integer(8) :: val = 0
    contains
        procedure :: a => int_mod_p_add
        procedure :: s => int_mod_p_sub
        procedure :: m => int_mod_p_mul
        procedure :: d => int_mod_p_div
        procedure :: zero => int_mod_p_zero
        procedure :: coerce_to_f64
        procedure, private :: is_equal
        generic :: operator(==) => is_equal
        generic :: operator(+) => a
        generic :: operator(-) => s
        generic :: operator(*) => m
        generic :: operator(/) => d
    end type IntModP

    interface IntModP
        module procedure init_int_mod_p
    end interface

contains

    function init_int_mod_p(v) result(res)
        integer(8), intent(in) :: v
        type(IntModP) :: res
        res%val = mod(v, IntModP_modulus)
        if (res%val < 0) res%val = res%val + IntModP_modulus
    end function init_int_mod_p

    pure function int_mod_p_add(this, other) result(res)
        class(IntModP), intent(in) :: this
        class(IntModP), intent(in) :: other
        type(IntModP) :: res
        res%val = mod(this%val + other%val, IntModP_modulus)
    end function int_mod_p_add

    pure function int_mod_p_sub(this, other) result(res)
        class(IntModP), intent(in) :: this
        class(IntModP), intent(in) :: other
        type(IntModP) :: res
        res%val = mod(this%val + IntModP_modulus - mod(other%val, IntModP_modulus), IntModP_modulus)
    end function int_mod_p_sub

    pure function int_mod_p_mul(this, other) result(res)
        class(IntModP), intent(in) :: this
        class(IntModP), intent(in) :: other
        type(IntModP) :: res
        res%val = mod(this%val * other%val, IntModP_modulus)
    end function int_mod_p_mul

    pure function int_mod_p_div(this, other) result(res)
        class(IntModP), intent(in) :: this
        class(IntModP), intent(in) :: other
        type(IntModP) :: res
        res%val = mod(this%val * mod_inverse(other%val, IntModP_modulus), IntModP_modulus)
    end function int_mod_p_div

    pure function int_mod_p_zero(this) result(res)
        class(IntModP), intent(in) :: this
        type(IntModP) :: res
        res%val = 0
    end function int_mod_p_zero

    pure function coerce_to_f64(this) result(res)
        class(IntModP), intent(in) :: this
        real(8) :: res
        res = real(this%val, 8)
    end function coerce_to_f64

    pure function is_equal(this, other) result(res)
        class(IntModP), intent(in) :: this
        class(IntModP), intent(in) :: other
        logical :: res
        res = (this%val == other%val)
    end function is_equal

    pure function mod_inverse(a, m) result(res)
        integer(8), intent(in) :: a
        integer(8), intent(in) :: m
        integer(8) :: res
        integer(8) :: a_i64, m_i64, m0, y, x, q, t

        if (m == 1) then
            res = 0
            return
        end if

        a_i64 = a
        m_i64 = m
        m0 = m_i64
        y = 0
        x = 1

        do while (a_i64 > 1)
            q = a_i64 / m_i64
            t = m_i64
            m_i64 = mod(a_i64, m_i64)
            a_i64 = t
            t = y
            y = x - q * y
            x = t
        end do

        if (x < 0) x = x + m0
        res = mod(x, m0)
    end function mod_inverse

end module IntModP_Module