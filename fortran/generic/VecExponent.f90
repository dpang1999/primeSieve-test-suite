module VecExponent_Module
    implicit none
    private

    public :: VecExponent

    type :: VecExponent
        integer, allocatable :: exps(:)
    contains
        procedure :: lex_compare
        procedure :: degree
        procedure :: can_reduce
        procedure :: lcm => lcm_exps
        procedure :: add => add_exps
        procedure :: sub => sub_exps
        procedure, private :: is_equal
        generic :: operator(==) => is_equal
    end type VecExponent

    interface VecExponent
        module procedure init_vec_exponent
    end interface

contains

    function init_vec_exponent(v) result(res)
        integer, intent(in) :: v(:)
        type(VecExponent) :: res
        allocate(res%exps(size(v)))
        res%exps = v
    end function init_vec_exponent

    pure function lex_compare(this, other) result(res)
        class(VecExponent), intent(in) :: this
        class(VecExponent), intent(in) :: other
        integer :: res
        integer :: i, min_len

        min_len = min(size(this%exps), size(other%exps))
        do i = 1, min_len
            if (this%exps(i) < other%exps(i)) then
                res = -1
                return
            end if
            if (this%exps(i) > other%exps(i)) then
                res = 1
                return
            end if
        end do

        if (size(this%exps) < size(other%exps)) then
            res = -1
            return
        end if
        if (size(this%exps) > size(other%exps)) then
            res = 1
            return
        end if

        res = 0
    end function lex_compare

    pure function degree(this) result(res)
        class(VecExponent), intent(in) :: this
        integer(8) :: res
        integer :: i
        res = 0
        if (allocated(this%exps)) then
            do i = 1, size(this%exps)
                res = res + this%exps(i)
            end do
        end if
    end function degree

    pure function can_reduce(this, other) result(res)
        class(VecExponent), intent(in) :: this
        class(VecExponent), intent(in) :: other
        logical :: res
        integer :: i, other_val

        if (.not. allocated(this%exps)) then
            if (.not. allocated(other%exps) .or. size(other%exps) == 0) then
                res = .true.
                return
            else
                res = .false.
                return
            end if
        end if

        do i = 1, size(this%exps)
            if (allocated(other%exps) .and. i <= size(other%exps)) then
                other_val = other%exps(i)
            else
                other_val = 0
            end if
            if (this%exps(i) < other_val) then
                res = .false.
                return
            end if
        end do
        res = .true.
    end function can_reduce

    pure function lcm_exps(this, other) result(res)
        class(VecExponent), intent(in) :: this
        class(VecExponent), intent(in) :: other
        type(VecExponent) :: res
        integer :: len, i, e1, e2

        if (.not. allocated(this%exps) .and. .not. allocated(other%exps)) then
            allocate(res%exps(0))
            return
        end if

        if (.not. allocated(this%exps)) then
            len = size(other%exps)
        else if (.not. allocated(other%exps)) then
            len = size(this%exps)
        else
            len = max(size(this%exps), size(other%exps))
        end if
        
        allocate(res%exps(len))
        do i = 1, len
            e1 = 0
            e2 = 0
            if (allocated(this%exps)) then
                if (i <= size(this%exps)) e1 = this%exps(i)
            end if
            if (allocated(other%exps)) then
                if (i <= size(other%exps)) e2 = other%exps(i)
            end if
            res%exps(i) = max(e1, e2)
        end do
    end function lcm_exps

    pure function add_exps(this, other) result(res)
        class(VecExponent), intent(in) :: this
        class(VecExponent), intent(in) :: other
        type(VecExponent) :: res
        integer :: len, i, e1, e2

        if (.not. allocated(this%exps) .and. .not. allocated(other%exps)) then
            allocate(res%exps(0))
            return
        end if

        if (.not. allocated(this%exps)) then
            len = size(other%exps)
        else if (.not. allocated(other%exps)) then
            len = size(this%exps)
        else
            len = max(size(this%exps), size(other%exps))
        end if
        
        allocate(res%exps(len))
        do i = 1, len
            e1 = 0
            e2 = 0
            if (allocated(this%exps)) then
                if (i <= size(this%exps)) e1 = this%exps(i)
            end if
            if (allocated(other%exps)) then
                if (i <= size(other%exps)) e2 = other%exps(i)
            end if
            res%exps(i) = e1 + e2
        end do
    end function add_exps

    pure function sub_exps(this, other) result(res)
        class(VecExponent), intent(in) :: this
        class(VecExponent), intent(in) :: other
        type(VecExponent) :: res
        integer :: len, i, e1, e2

        if (.not. allocated(this%exps) .and. .not. allocated(other%exps)) then
            allocate(res%exps(0))
            return
        end if

        if (.not. allocated(this%exps)) then
            len = size(other%exps)
        else if (.not. allocated(other%exps)) then
            len = size(this%exps)
        else
            len = max(size(this%exps), size(other%exps))
        end if
        
        allocate(res%exps(len))
        do i = 1, len
            e1 = 0
            e2 = 0
            if (allocated(this%exps)) then
                if (i <= size(this%exps)) e1 = this%exps(i)
            end if
            if (allocated(other%exps)) then
                if (i <= size(other%exps)) e2 = other%exps(i)
            end if
            res%exps(i) = e1 - e2
        end do
    end function sub_exps

    pure function is_equal(this, other) result(res)
        class(VecExponent), intent(in) :: this
        class(VecExponent), intent(in) :: other
        logical :: res
        integer :: i

        if (.not. allocated(this%exps) .and. .not. allocated(other%exps)) then
            res = .true.
            return
        end if

        if (.not. allocated(this%exps) .or. .not. allocated(other%exps)) then
            res = .false.
            return
        end if

        if (size(this%exps) /= size(other%exps)) then
            res = .false.
            return
        end if

        do i = 1, size(this%exps)
            if (this%exps(i) /= other%exps(i)) then
                res = .false.
                return
            end if
        end do
        res = .true.
    end function is_equal

end module VecExponent_Module