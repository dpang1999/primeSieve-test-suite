module BitPackedExponent_Module
    implicit none
    private

    public :: BitPackedExponent

    type :: BitPackedExponent
        integer(8) :: packed = 0
    contains
        procedure :: lex_compare
        procedure :: degree
        procedure :: can_reduce
        procedure :: lcm => lcm_packed
        procedure :: add => add_packed
        procedure :: sub => sub_packed
        procedure, private :: is_equal
        generic :: operator(==) => is_equal
    end type BitPackedExponent

    interface BitPackedExponent
        module procedure init_bit_packed
        module procedure init_bit_packed_scalar
    end interface

contains

    function init_bit_packed_scalar(p) result(res)
        integer(8), intent(in) :: p
        type(BitPackedExponent) :: res
        res%packed = p
    end function init_bit_packed_scalar

    function init_bit_packed(v) result(res)
        integer, intent(in) :: v(:)
        type(BitPackedExponent) :: res
        integer(8) :: shift, deg, val
        integer :: i, exps(6)

        res%packed = 0
        exps = 0
        do i = 1, min(size(v), 6)
            exps(i) = v(i)
        end do

        do i = 0, 5
            shift = 40 - 8 * i
            val = iand(int(exps(i+1), 8), 255_8)
            res%packed = ior(res%packed, ishft(val, shift))
        end do

        deg = 0
        do i = 1, 6
            deg = deg + exps(i)
        end do

        res%packed = ior(res%packed, ishft(iand(deg, 65535_8), 48))
    end function init_bit_packed

    pure function lex_compare(this, other) result(res)
        class(BitPackedExponent), intent(in) :: this
        class(BitPackedExponent), intent(in) :: other
        integer :: res
        integer(8) :: e1, e2
        
        e1 = iand(this%packed, z'0000FFFFFFFFFFFF')
        e2 = iand(other%packed, z'0000FFFFFFFFFFFF')

        if (e1 < e2) then
            res = -1
        else if (e1 > e2) then
            res = 1
        else
            res = 0
        end if
    end function lex_compare

    pure function degree(this) result(res)
        class(BitPackedExponent), intent(in) :: this
        integer(8) :: res
        res = iand(ishft(this%packed, -48), 65535_8)
    end function degree

    pure function can_reduce(this, other) result(res)
        class(BitPackedExponent), intent(in) :: this
        class(BitPackedExponent), intent(in) :: other
        logical :: res
        integer :: i
        integer(8) :: self_exp, divisor_exp

        do i = 0, 40, 8
            self_exp = iand(ishft(this%packed, -i), 255_8)
            divisor_exp = iand(ishft(other%packed, -i), 255_8)
            if (self_exp < divisor_exp) then
                res = .false.
                return
            end if
        end do
        res = .true.
    end function can_reduce

    pure function lcm_packed(this, other) result(res)
        class(BitPackedExponent), intent(in) :: this
        class(BitPackedExponent), intent(in) :: other
        type(BitPackedExponent) :: res
        integer :: i
        integer(8) :: self_exponents, other_exponents, lcm_exponents, deg
        integer(8) :: self_exp, other_exp, lcm_exp

        self_exponents = iand(this%packed, z'0000FFFFFFFFFFFF')
        other_exponents = iand(other%packed, z'0000FFFFFFFFFFFF')
        lcm_exponents = 0
        deg = 0

        do i = 0, 40, 8
            self_exp = iand(ishft(self_exponents, -i), 255_8)
            other_exp = iand(ishft(other_exponents, -i), 255_8)
            lcm_exp = max(self_exp, other_exp)
            lcm_exponents = ior(lcm_exponents, ishft(lcm_exp, i))
            deg = deg + lcm_exp
        end do

        lcm_exponents = ior(lcm_exponents, ishft(iand(deg, 65535_8), 48))
        res%packed = lcm_exponents
    end function lcm_packed

    pure function add_packed(this, other) result(res)
        class(BitPackedExponent), intent(in) :: this
        class(BitPackedExponent), intent(in) :: other
        type(BitPackedExponent) :: res
        res%packed = this%packed + other%packed
    end function add_packed

    pure function sub_packed(this, other) result(res)
        class(BitPackedExponent), intent(in) :: this
        class(BitPackedExponent), intent(in) :: other
        type(BitPackedExponent) :: res
        res%packed = this%packed - other%packed
    end function sub_packed

    pure function is_equal(this, other) result(res)
        class(BitPackedExponent), intent(in) :: this
        class(BitPackedExponent), intent(in) :: other
        logical :: res
        res = (this%packed == other%packed)
    end function is_equal

end module BitPackedExponent_Module