module PolyMod
    implicit none

    integer(8), parameter :: MODULUS = 7

    type :: Term
        integer(8) :: coeff = 0
        integer, allocatable :: exps(:)
    end type Term

    type :: Polynomial
        type(Term), allocatable :: terms(:)
    end type Polynomial

    interface operator(==)
        module procedure poly_equal
    end interface

contains

    function mod_inverse(a, m) result(res)
        integer(8), intent(in) :: a, m
        integer(8) :: res, a_i64, m_i64, m0, q, temp, x, y

        if (m == 1) then
            res = 0
            return
        end if

        a_i64 = a
        m_i64 = m
        m0 = m
        y = 0
        x = 1

        do while (a_i64 > 1)
            q = a_i64 / m_i64
            temp = m_i64
            m_i64 = mod(a_i64, m_i64)
            a_i64 = temp
            temp = y
            y = x - q * y
            x = temp
        end do

        if (x < 0) x = x + m0
        res = mod(x, m0)
    end function mod_inverse

    integer function compare_terms(t1, t2)
        type(Term), intent(in) :: t1, t2
        integer :: i

        compare_terms = 0
        do i = 1, min(size(t1%exps), size(t2%exps))
            if (t1%exps(i) < t2%exps(i)) then
                compare_terms = -1
                return
            else if (t1%exps(i) > t2%exps(i)) then
                compare_terms = 1
                return
            end if
        end do

        if (size(t1%exps) < size(t2%exps)) then
            compare_terms = -1
        else if (size(t1%exps) > size(t2%exps)) then
            compare_terms = 1
        end if
    end function compare_terms

    logical function poly_equal(p1, p2)
        type(Polynomial), intent(in) :: p1, p2
        integer :: i, j

        poly_equal = .false.
        if (.not. allocated(p1%terms) .and. .not. allocated(p2%terms)) then
            poly_equal = .true.
            return
        end if
        if (.not. allocated(p1%terms) .or. .not. allocated(p2%terms)) return
        if (size(p1%terms) /= size(p2%terms)) return

        do i = 1, size(p1%terms)
            if (p1%terms(i)%coeff /= p2%terms(i)%coeff) return
            if (size(p1%terms(i)%exps) /= size(p2%terms(i)%exps)) return
            do j = 1, size(p1%terms(i)%exps)
                if (p1%terms(i)%exps(j) /= p2%terms(i)%exps(j)) return
            end do
        end do
        poly_equal = .true.
    end function poly_equal

    logical function exps_equal(e1, e2)
        integer, intent(in) :: e1(:), e2(:)
        integer :: i
        exps_equal = .false.
        if (size(e1) /= size(e2)) return
        do i = 1, size(e1)
            if (e1(i) /= e2(i)) return
        end do
        exps_equal = .true.
    end function exps_equal

    subroutine sort_terms(terms)
        type(Term), allocatable, intent(inout) :: terms(:)
        type(Term) :: temp
        integer :: i, j, n, max_idx

        if (.not. allocated(terms)) return
        n = size(terms)
        ! Selection sort (descending order)
        do i = 1, n - 1
            max_idx = i
            do j = i + 1, n
                if (compare_terms(terms(j), terms(max_idx)) > 0) then
                    max_idx = j
                end if
            end do
            if (max_idx /= i) then
                temp = terms(i)
                terms(i) = terms(max_idx)
                terms(max_idx) = temp
            end if
        end do
    end subroutine sort_terms

    function copy_poly(p) result(res)
        type(Polynomial), intent(in) :: p
        type(Polynomial) :: res
        integer :: i
        if (allocated(p%terms)) then
            allocate(res%terms(size(p%terms)))
            do i = 1, size(p%terms)
                res%terms(i)%coeff = p%terms(i)%coeff
                allocate(res%terms(i)%exps(size(p%terms(i)%exps)))
                res%terms(i)%exps = p%terms(i)%exps
            end do
        end if
    end function copy_poly

    function add_poly(p1, p2) result(res)
        type(Polynomial), intent(in) :: p1, p2
        type(Polynomial) :: res
        type(Term), allocatable :: temp_terms(:), result_terms(:)
        integer :: i, j, count
        logical :: found

        if (.not. allocated(p1%terms)) then
            res = copy_poly(p2)
            return
        end if
        if (.not. allocated(p2%terms)) then
            res = copy_poly(p1)
            return
        end if

        allocate(temp_terms(size(p1%terms) + size(p2%terms)))
        count = size(p1%terms)
        do i = 1, size(p1%terms)
            temp_terms(i)%coeff = p1%terms(i)%coeff
            allocate(temp_terms(i)%exps(size(p1%terms(i)%exps)))
            temp_terms(i)%exps = p1%terms(i)%exps
        end do

        do i = 1, size(p2%terms)
            found = .false.
            do j = 1, count
                if (exps_equal(temp_terms(j)%exps, p2%terms(i)%exps)) then
                    temp_terms(j)%coeff = mod(temp_terms(j)%coeff + p2%terms(i)%coeff, MODULUS)
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                count = count + 1
                temp_terms(count)%coeff = mod(p2%terms(i)%coeff, MODULUS)
                allocate(temp_terms(count)%exps(size(p2%terms(i)%exps)))
                temp_terms(count)%exps = p2%terms(i)%exps
            end if
        end do

        ! Filter out coeff = 0
        j = 0
        do i = 1, count
            if (temp_terms(i)%coeff /= 0) j = j + 1
        end do
        allocate(res%terms(j))
        j = 0
        do i = 1, count
            if (temp_terms(i)%coeff /= 0) then
                j = j + 1
                res%terms(j)%coeff = temp_terms(i)%coeff
                allocate(res%terms(j)%exps(size(temp_terms(i)%exps)))
                res%terms(j)%exps = temp_terms(i)%exps
            end if
        end do
        
        call sort_terms(res%terms)
    end function add_poly

    function sub_poly(p1, p2) result(res)
        type(Polynomial), intent(in) :: p1, p2
        type(Polynomial) :: res
        type(Term), allocatable :: temp_terms(:), result_terms(:)
        integer :: i, j, count
        logical :: found

        if (.not. allocated(p1%terms)) then
            allocate(res%terms(size(p2%terms)))
            do i = 1, size(p2%terms)
                res%terms(i)%coeff = mod(MODULUS - mod(p2%terms(i)%coeff, MODULUS), MODULUS)
                allocate(res%terms(i)%exps(size(p2%terms(i)%exps)))
                res%terms(i)%exps = p2%terms(i)%exps
            end do
            return
        end if
        if (.not. allocated(p2%terms)) then
            res = copy_poly(p1)
            return
        end if

        allocate(temp_terms(size(p1%terms) + size(p2%terms)))
        count = size(p1%terms)
        do i = 1, size(p1%terms)
            temp_terms(i)%coeff = p1%terms(i)%coeff
            allocate(temp_terms(i)%exps(size(p1%terms(i)%exps)))
            temp_terms(i)%exps = p1%terms(i)%exps
        end do

        do i = 1, size(p2%terms)
            found = .false.
            do j = 1, count
                if (exps_equal(temp_terms(j)%exps, p2%terms(i)%exps)) then
                    temp_terms(j)%coeff = mod(MODULUS + temp_terms(j)%coeff - mod(p2%terms(i)%coeff, MODULUS), MODULUS)
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                count = count + 1
                temp_terms(count)%coeff = mod(MODULUS - mod(p2%terms(i)%coeff, MODULUS), MODULUS)
                allocate(temp_terms(count)%exps(size(p2%terms(i)%exps)))
                temp_terms(count)%exps = p2%terms(i)%exps
            end if
        end do

        ! Filter out coeff = 0
        j = 0
        do i = 1, count
            if (temp_terms(i)%coeff /= 0) j = j + 1
        end do
        allocate(res%terms(j))
        j = 0
        do i = 1, count
            if (temp_terms(i)%coeff /= 0) then
                j = j + 1
                res%terms(j)%coeff = temp_terms(i)%coeff
                allocate(res%terms(j)%exps(size(temp_terms(i)%exps)))
                res%terms(j)%exps = temp_terms(i)%exps
            end if
        end do
        
        call sort_terms(res%terms)
    end function sub_poly

    function multiply_by_term(p, t) result(res)
        type(Polynomial), intent(in) :: p
        type(Term), intent(in) :: t
        type(Polynomial) :: res
        integer :: i

        if (.not. allocated(p%terms)) return
        allocate(res%terms(size(p%terms)))
        do i = 1, size(p%terms)
            res%terms(i)%coeff = mod(p%terms(i)%coeff * t%coeff, MODULUS)
            allocate(res%terms(i)%exps(size(p%terms(i)%exps)))
            res%terms(i)%exps = p%terms(i)%exps + t%exps
        end do
    end function multiply_by_term

    function leading_term(p) result(res)
        type(Polynomial), intent(in) :: p
        type(Term) :: res
        if (allocated(p%terms)) then
            if (size(p%terms) > 0) then
                res%coeff = p%terms(1)%coeff
                allocate(res%exps(size(p%terms(1)%exps)))
                res%exps = p%terms(1)%exps
            else
                res%coeff = 0
                allocate(res%exps(0))
            end if
        else
            res%coeff = 0
            allocate(res%exps(0))
        end if
    end function leading_term

    function s_polynomial(p1, p2) result(res)
        type(Polynomial), intent(in) :: p1, p2
        type(Polynomial) :: res
        type(Term) :: lt1, lt2, s1_term, s2_term
        integer :: i

        if (.not. allocated(p1%terms) .or. .not. allocated(p2%terms)) return
        if (size(p1%terms) == 0 .or. size(p2%terms) == 0) return

        lt1 = leading_term(p1)
        lt2 = leading_term(p2)

        s1_term%coeff = lt2%coeff
        s2_term%coeff = lt1%coeff
        allocate(s1_term%exps(size(lt1%exps)))
        allocate(s2_term%exps(size(lt2%exps)))

        do i = 1, size(lt1%exps)
            s1_term%exps(i) = max(lt1%exps(i), lt2%exps(i)) - lt1%exps(i)
            s2_term%exps(i) = max(lt1%exps(i), lt2%exps(i)) - lt2%exps(i)
        end do

        res = sub_poly(multiply_by_term(p1, s1_term), multiply_by_term(p2, s2_term))
    end function s_polynomial

    function make_monic(p) result(res)
        type(Polynomial), intent(in) :: p
        type(Polynomial) :: res
        integer(8) :: inv_coeff
        integer :: i

        if (.not. allocated(p%terms)) return
        if (size(p%terms) == 0) return

        inv_coeff = mod_inverse(p%terms(1)%coeff, MODULUS)
        res = copy_poly(p)
        do i = 1, size(res%terms)
            res%terms(i)%coeff = mod(res%terms(i)%coeff * inv_coeff, MODULUS)
        end do
    end function make_monic

    function reduce_polynomial(p, basis) result(res)
        type(Polynomial), intent(in) :: p
        type(Polynomial), intent(in) :: basis(:)
        type(Polynomial) :: res, poly, diff, cur_basis
        type(Term) :: lt_poly, lt_b, scale_term
        type(Term), allocatable :: temp_terms(:), temp_rem(:)
        integer :: i, j, k, remainder_cnt
        logical :: reduced, can_reduce
        type(Term), allocatable :: remainder(:)

        poly = copy_poly(p)
        allocate(remainder(0))

        do while (allocated(poly%terms) .and. size(poly%terms) > 0)
            reduced = .false.
            lt_poly = leading_term(poly)

            do i = 1, size(basis)
                if (.not. allocated(basis(i)%terms)) cycle
                if (size(basis(i)%terms) == 0) cycle
                
                lt_b = leading_term(basis(i))
                can_reduce = .true.
                do j = 1, size(lt_poly%exps)
                    if (lt_poly%exps(j) < lt_b%exps(j)) then
                        can_reduce = .false.
                        exit
                    end if
                end do

                if (can_reduce) then
                    scale_term%coeff = mod(lt_poly%coeff * mod_inverse(lt_b%coeff, MODULUS), MODULUS)
                    allocate(scale_term%exps(size(lt_poly%exps)))
                    do j = 1, size(lt_poly%exps)
                        scale_term%exps(j) = lt_poly%exps(j) - lt_b%exps(j)
                    end do
                    
                    poly = sub_poly(poly, multiply_by_term(basis(i), scale_term))
                    deallocate(scale_term%exps)
                    reduced = .true.
                    exit
                end if
            end do

            if (.not. reduced) then
                allocate(temp_rem(size(remainder) + 1))
                do j = 1, size(remainder)
                    temp_rem(j)%coeff = remainder(j)%coeff
                    allocate(temp_rem(j)%exps(size(remainder(j)%exps)))
                    temp_rem(j)%exps = remainder(j)%exps
                end do
                temp_rem(size(remainder) + 1)%coeff = poly%terms(1)%coeff
                allocate(temp_rem(size(remainder) + 1)%exps(size(poly%terms(1)%exps)))
                temp_rem(size(remainder) + 1)%exps = poly%terms(1)%exps
                call move_alloc(temp_rem, remainder)
                
                if (size(poly%terms) > 1) then
                    allocate(temp_terms(size(poly%terms) - 1))
                    do j = 2, size(poly%terms)
                        temp_terms(j-1)%coeff = poly%terms(j)%coeff
                        allocate(temp_terms(j-1)%exps(size(poly%terms(j)%exps)))
                        temp_terms(j-1)%exps = poly%terms(j)%exps
                    end do
                    call move_alloc(temp_terms, poly%terms)
                else
                    deallocate(poly%terms)
                    allocate(poly%terms(0))
                end if
            end if
        end do

        allocate(res%terms(size(remainder)))
        do i = 1, size(remainder)
            res%terms(i)%coeff = remainder(i)%coeff
            allocate(res%terms(i)%exps(size(remainder(i)%exps)))
            res%terms(i)%exps = remainder(i)%exps
        end do
        call sort_terms(res%terms)
    end function reduce_polynomial

    function in_basis(p, basis) result(res)
        type(Polynomial), intent(in) :: p
        type(Polynomial), intent(in) :: basis(:)
        logical :: res
        integer :: i
        res = .false.
        do i = 1, size(basis)
            if (p == basis(i)) then
                res = .true.
                return
            end if
        end do
    end function in_basis

    function grobner_basis(basis_in) result(reduced_basis)
        type(Polynomial), intent(in) :: basis_in(:)
        type(Polynomial), allocatable :: reduced_basis(:)
        type(Polynomial), allocatable :: cur_basis(:), temp_basis(:), basis_excluding_self(:)
        type(Polynomial) :: s_poly, h, temp_poly
        integer, allocatable :: pair_set(:, :), temp_pairs(:, :)
        integer :: i, j, k, n, num_pairs, p1, p2

        allocate(cur_basis(size(basis_in)))
        do i = 1, size(basis_in)
            cur_basis(i) = copy_poly(basis_in(i))
        end do

        num_pairs = (size(basis_in) * (size(basis_in) - 1)) / 2
        allocate(pair_set(2, num_pairs))
        k = 1
        do i = 1, size(basis_in)
            do j = i + 1, size(basis_in)
                pair_set(1, k) = i
                pair_set(2, k) = j
                k = k + 1
            end do
        end do

        do while (size(pair_set, 2) > 0)
            p1 = pair_set(1, 1)
            p2 = pair_set(2, 1)

            ! Remove first pair
            allocate(temp_pairs(2, size(pair_set, 2) - 1))
            if (size(pair_set, 2) > 1) then
                temp_pairs = pair_set(:, 2:)
            end if
            call move_alloc(temp_pairs, pair_set)

            s_poly = s_polynomial(cur_basis(p1), cur_basis(p2))
            h = reduce_polynomial(s_poly, cur_basis)

            if (allocated(h%terms) .and. size(h%terms) > 0) then
                if (.not. in_basis(h, cur_basis)) then
                    n = size(cur_basis)
                    allocate(temp_basis(n + 1))
                    do i = 1, n
                        temp_basis(i) = copy_poly(cur_basis(i))
                    end do
                    temp_basis(n + 1) = copy_poly(h)
                    call move_alloc(temp_basis, cur_basis)

                    allocate(temp_pairs(2, size(pair_set, 2) + n))
                    if (size(pair_set, 2) > 0) then
                        temp_pairs(:, 1:size(pair_set, 2)) = pair_set
                    end if
                    do i = 1, n
                        temp_pairs(1, size(pair_set, 2) + i) = i
                        temp_pairs(2, size(pair_set, 2) + i) = n + 1
                    end do
                    call move_alloc(temp_pairs, pair_set)
                end if
            end if
        end do

        allocate(reduced_basis(0))
        do i = 1, size(cur_basis)
            allocate(basis_excluding_self(size(cur_basis) - 1))
            k = 1
            do j = 1, size(cur_basis)
                if (i /= j) then
                    basis_excluding_self(k) = copy_poly(cur_basis(j))
                    k = k + 1
                end if
            end do

            temp_poly = reduce_polynomial(cur_basis(i), basis_excluding_self)
            deallocate(basis_excluding_self)

            if (allocated(temp_poly%terms) .and. size(temp_poly%terms) > 0) then
                temp_poly = make_monic(temp_poly)
                if (.not. in_basis(temp_poly, reduced_basis)) then
                    allocate(temp_basis(size(reduced_basis) + 1))
                    do j = 1, size(reduced_basis)
                        temp_basis(j) = copy_poly(reduced_basis(j))
                    end do
                    temp_basis(size(reduced_basis) + 1) = copy_poly(temp_poly)
                    call move_alloc(temp_basis, reduced_basis)
                end if
            end if
        end do
    end function grobner_basis

end module PolyMod

program main
    use PolyMod
    implicit none

    integer :: n_vars = 4, i, j, k, v
    character(len=32) :: arg
    type(Polynomial), allocatable :: polynomials(:), reduced_basis(:)
    type(Polynomial) :: curr_poly
    type(Term), allocatable :: t_arr(:)
    type(Term) :: t
    
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        read(arg, *) n_vars
    end if

    if (n_vars == 4) then
        print *, "Fortran specialized vec exponent cyclic 4"
        allocate(polynomials(4))

        allocate(t_arr(4))
        ! P1
        t_arr(1)%coeff = 1; allocate(t_arr(1)%exps(4)); t_arr(1)%exps = [1, 0, 0, 0]
        t_arr(2)%coeff = 1; allocate(t_arr(2)%exps(4)); t_arr(2)%exps = [0, 1, 0, 0]
        t_arr(3)%coeff = 1; allocate(t_arr(3)%exps(4)); t_arr(3)%exps = [0, 0, 1, 0]
        t_arr(4)%coeff = 1; allocate(t_arr(4)%exps(4)); t_arr(4)%exps = [0, 0, 0, 1]
        allocate(polynomials(1)%terms(4))
        do i = 1, 4; polynomials(1)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(1)%terms)

        ! P2
        t_arr(1)%coeff = 1; deallocate(t_arr(1)%exps); allocate(t_arr(1)%exps(4)); t_arr(1)%exps = [1, 1, 0, 0]
        t_arr(2)%coeff = 1; deallocate(t_arr(2)%exps); allocate(t_arr(2)%exps(4)); t_arr(2)%exps = [0, 1, 1, 0]
        t_arr(3)%coeff = 1; deallocate(t_arr(3)%exps); allocate(t_arr(3)%exps(4)); t_arr(3)%exps = [0, 0, 1, 1]
        t_arr(4)%coeff = 1; deallocate(t_arr(4)%exps); allocate(t_arr(4)%exps(4)); t_arr(4)%exps = [1, 0, 0, 1]
        allocate(polynomials(2)%terms(4))
        do i = 1, 4; polynomials(2)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(2)%terms)

        ! P3
        t_arr(1)%coeff = 1; deallocate(t_arr(1)%exps); allocate(t_arr(1)%exps(4)); t_arr(1)%exps = [1, 1, 1, 0]
        t_arr(2)%coeff = 1; deallocate(t_arr(2)%exps); allocate(t_arr(2)%exps(4)); t_arr(2)%exps = [0, 1, 1, 1]
        t_arr(3)%coeff = 1; deallocate(t_arr(3)%exps); allocate(t_arr(3)%exps(4)); t_arr(3)%exps = [1, 0, 1, 1]
        t_arr(4)%coeff = 1; deallocate(t_arr(4)%exps); allocate(t_arr(4)%exps(4)); t_arr(4)%exps = [1, 1, 0, 1]
        allocate(polynomials(3)%terms(4))
        do i = 1, 4; polynomials(3)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(3)%terms)

        ! P4
        deallocate(t_arr)
        allocate(t_arr(2))
        t_arr(1)%coeff = 1; allocate(t_arr(1)%exps(4)); t_arr(1)%exps = [1, 1, 1, 1]
        t_arr(2)%coeff = 6; allocate(t_arr(2)%exps(4)); t_arr(2)%exps = [0, 0, 0, 0]
        allocate(polynomials(4)%terms(2))
        do i = 1, 2; polynomials(4)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(4)%terms)

    else if (n_vars == 5) then
        print *, "Fortran specialized vec exponent cyclic 5"
        allocate(polynomials(5))
        
        allocate(t_arr(5))
        ! P1
        t_arr(1)%coeff = 1; allocate(t_arr(1)%exps(5)); t_arr(1)%exps = [1, 0, 0, 0, 0]
        t_arr(2)%coeff = 1; allocate(t_arr(2)%exps(5)); t_arr(2)%exps = [0, 1, 0, 0, 0]
        t_arr(3)%coeff = 1; allocate(t_arr(3)%exps(5)); t_arr(3)%exps = [0, 0, 1, 0, 0]
        t_arr(4)%coeff = 1; allocate(t_arr(4)%exps(5)); t_arr(4)%exps = [0, 0, 0, 1, 0]
        t_arr(5)%coeff = 1; allocate(t_arr(5)%exps(5)); t_arr(5)%exps = [0, 0, 0, 0, 1]
        allocate(polynomials(1)%terms(5))
        do i = 1, 5; polynomials(1)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(1)%terms)

        ! P2
        t_arr(1)%coeff = 1; deallocate(t_arr(1)%exps); allocate(t_arr(1)%exps(5)); t_arr(1)%exps = [1, 1, 0, 0, 0]
        t_arr(2)%coeff = 1; deallocate(t_arr(2)%exps); allocate(t_arr(2)%exps(5)); t_arr(2)%exps = [0, 1, 1, 0, 0]
        t_arr(3)%coeff = 1; deallocate(t_arr(3)%exps); allocate(t_arr(3)%exps(5)); t_arr(3)%exps = [0, 0, 1, 1, 0]
        t_arr(4)%coeff = 1; deallocate(t_arr(4)%exps); allocate(t_arr(4)%exps(5)); t_arr(4)%exps = [0, 0, 0, 1, 1]
        t_arr(5)%coeff = 1; deallocate(t_arr(5)%exps); allocate(t_arr(5)%exps(5)); t_arr(5)%exps = [1, 0, 0, 0, 1]
        allocate(polynomials(2)%terms(5))
        do i = 1, 5; polynomials(2)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(2)%terms)

        ! P3
        t_arr(1)%coeff = 1; deallocate(t_arr(1)%exps); allocate(t_arr(1)%exps(5)); t_arr(1)%exps = [1, 1, 1, 0, 0]
        t_arr(2)%coeff = 1; deallocate(t_arr(2)%exps); allocate(t_arr(2)%exps(5)); t_arr(2)%exps = [0, 1, 1, 1, 0]
        t_arr(3)%coeff = 1; deallocate(t_arr(3)%exps); allocate(t_arr(3)%exps(5)); t_arr(3)%exps = [0, 0, 1, 1, 1]
        t_arr(4)%coeff = 1; deallocate(t_arr(4)%exps); allocate(t_arr(4)%exps(5)); t_arr(4)%exps = [1, 0, 0, 1, 1]
        t_arr(5)%coeff = 1; deallocate(t_arr(5)%exps); allocate(t_arr(5)%exps(5)); t_arr(5)%exps = [1, 1, 0, 0, 1]
        allocate(polynomials(3)%terms(5))
        do i = 1, 5; polynomials(3)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(3)%terms)

        ! P4
        t_arr(1)%coeff = 1; deallocate(t_arr(1)%exps); allocate(t_arr(1)%exps(5)); t_arr(1)%exps = [1, 1, 1, 1, 0]
        t_arr(2)%coeff = 1; deallocate(t_arr(2)%exps); allocate(t_arr(2)%exps(5)); t_arr(2)%exps = [0, 1, 1, 1, 1]
        t_arr(3)%coeff = 1; deallocate(t_arr(3)%exps); allocate(t_arr(3)%exps(5)); t_arr(3)%exps = [1, 0, 1, 1, 1]
        t_arr(4)%coeff = 1; deallocate(t_arr(4)%exps); allocate(t_arr(4)%exps(5)); t_arr(4)%exps = [1, 1, 0, 1, 1]
        t_arr(5)%coeff = 1; deallocate(t_arr(5)%exps); allocate(t_arr(5)%exps(5)); t_arr(5)%exps = [1, 1, 1, 0, 1]
        allocate(polynomials(4)%terms(5))
        do i = 1, 5; polynomials(4)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(4)%terms)

        ! P5
        deallocate(t_arr)
        allocate(t_arr(2))
        t_arr(1)%coeff = 1; allocate(t_arr(1)%exps(5)); t_arr(1)%exps = [1, 1, 1, 1, 1]
        t_arr(2)%coeff = 6; allocate(t_arr(2)%exps(5)); t_arr(2)%exps = [0, 0, 0, 0, 0]
        allocate(polynomials(5)%terms(2))
        do i = 1, 2; polynomials(5)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(5)%terms)

    else if (n_vars == 6) then
        print *, "Fortran specialized vec exponent cyclic 6"
        allocate(polynomials(6))
        
        allocate(t_arr(6))
        ! P1
        t_arr(1)%coeff = 1; allocate(t_arr(1)%exps(6)); t_arr(1)%exps = [1, 0, 0, 0, 0, 0]
        t_arr(2)%coeff = 1; allocate(t_arr(2)%exps(6)); t_arr(2)%exps = [0, 1, 0, 0, 0, 0]
        t_arr(3)%coeff = 1; allocate(t_arr(3)%exps(6)); t_arr(3)%exps = [0, 0, 1, 0, 0, 0]
        t_arr(4)%coeff = 1; allocate(t_arr(4)%exps(6)); t_arr(4)%exps = [0, 0, 0, 1, 0, 0]
        t_arr(5)%coeff = 1; allocate(t_arr(5)%exps(6)); t_arr(5)%exps = [0, 0, 0, 0, 1, 0]
        t_arr(6)%coeff = 1; allocate(t_arr(6)%exps(6)); t_arr(6)%exps = [0, 0, 0, 0, 0, 1]
        allocate(polynomials(1)%terms(6))
        do i = 1, 6; polynomials(1)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(1)%terms)

        ! P2
        t_arr(1)%coeff = 1; deallocate(t_arr(1)%exps); allocate(t_arr(1)%exps(6)); t_arr(1)%exps = [1, 1, 0, 0, 0, 0]
        t_arr(2)%coeff = 1; deallocate(t_arr(2)%exps); allocate(t_arr(2)%exps(6)); t_arr(2)%exps = [0, 1, 1, 0, 0, 0]
        t_arr(3)%coeff = 1; deallocate(t_arr(3)%exps); allocate(t_arr(3)%exps(6)); t_arr(3)%exps = [0, 0, 1, 1, 0, 0]
        t_arr(4)%coeff = 1; deallocate(t_arr(4)%exps); allocate(t_arr(4)%exps(6)); t_arr(4)%exps = [0, 0, 0, 1, 1, 0]
        t_arr(5)%coeff = 1; deallocate(t_arr(5)%exps); allocate(t_arr(5)%exps(6)); t_arr(5)%exps = [0, 0, 0, 0, 1, 1]
        t_arr(6)%coeff = 1; deallocate(t_arr(6)%exps); allocate(t_arr(6)%exps(6)); t_arr(6)%exps = [1, 0, 0, 0, 0, 1]
        allocate(polynomials(2)%terms(6))
        do i = 1, 6; polynomials(2)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(2)%terms)

        ! P3
        t_arr(1)%coeff = 1; deallocate(t_arr(1)%exps); allocate(t_arr(1)%exps(6)); t_arr(1)%exps = [1, 1, 1, 0, 0, 0]
        t_arr(2)%coeff = 1; deallocate(t_arr(2)%exps); allocate(t_arr(2)%exps(6)); t_arr(2)%exps = [0, 1, 1, 1, 0, 0]
        t_arr(3)%coeff = 1; deallocate(t_arr(3)%exps); allocate(t_arr(3)%exps(6)); t_arr(3)%exps = [0, 0, 1, 1, 1, 0]
        t_arr(4)%coeff = 1; deallocate(t_arr(4)%exps); allocate(t_arr(4)%exps(6)); t_arr(4)%exps = [0, 0, 0, 1, 1, 1]
        t_arr(5)%coeff = 1; deallocate(t_arr(5)%exps); allocate(t_arr(5)%exps(6)); t_arr(5)%exps = [1, 0, 0, 0, 1, 1]
        t_arr(6)%coeff = 1; deallocate(t_arr(6)%exps); allocate(t_arr(6)%exps(6)); t_arr(6)%exps = [1, 1, 0, 0, 0, 1]
        allocate(polynomials(3)%terms(6))
        do i = 1, 6; polynomials(3)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(3)%terms)

        ! P4
        t_arr(1)%coeff = 1; deallocate(t_arr(1)%exps); allocate(t_arr(1)%exps(6)); t_arr(1)%exps = [1, 1, 1, 1, 0, 0]
        t_arr(2)%coeff = 1; deallocate(t_arr(2)%exps); allocate(t_arr(2)%exps(6)); t_arr(2)%exps = [0, 1, 1, 1, 1, 0]
        t_arr(3)%coeff = 1; deallocate(t_arr(3)%exps); allocate(t_arr(3)%exps(6)); t_arr(3)%exps = [0, 0, 1, 1, 1, 1]
        t_arr(4)%coeff = 1; deallocate(t_arr(4)%exps); allocate(t_arr(4)%exps(6)); t_arr(4)%exps = [1, 0, 0, 1, 1, 1]
        t_arr(5)%coeff = 1; deallocate(t_arr(5)%exps); allocate(t_arr(5)%exps(6)); t_arr(5)%exps = [1, 1, 0, 0, 1, 1]
        t_arr(6)%coeff = 1; deallocate(t_arr(6)%exps); allocate(t_arr(6)%exps(6)); t_arr(6)%exps = [1, 1, 1, 0, 0, 1]
        allocate(polynomials(4)%terms(6))
        do i = 1, 6; polynomials(4)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(4)%terms)

        ! P5
        t_arr(1)%coeff = 1; deallocate(t_arr(1)%exps); allocate(t_arr(1)%exps(6)); t_arr(1)%exps = [1, 1, 1, 1, 1, 0]
        t_arr(2)%coeff = 1; deallocate(t_arr(2)%exps); allocate(t_arr(2)%exps(6)); t_arr(2)%exps = [0, 1, 1, 1, 1, 1]
        t_arr(3)%coeff = 1; deallocate(t_arr(3)%exps); allocate(t_arr(3)%exps(6)); t_arr(3)%exps = [1, 0, 1, 1, 1, 1]
        t_arr(4)%coeff = 1; deallocate(t_arr(4)%exps); allocate(t_arr(4)%exps(6)); t_arr(4)%exps = [1, 1, 0, 1, 1, 1]
        t_arr(5)%coeff = 1; deallocate(t_arr(5)%exps); allocate(t_arr(5)%exps(6)); t_arr(5)%exps = [1, 1, 1, 0, 1, 1]
        t_arr(6)%coeff = 1; deallocate(t_arr(6)%exps); allocate(t_arr(6)%exps(6)); t_arr(6)%exps = [1, 1, 1, 1, 0, 1]
        allocate(polynomials(5)%terms(6))
        do i = 1, 6; polynomials(5)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(5)%terms)

        ! P6
        deallocate(t_arr)
        allocate(t_arr(2))
        t_arr(1)%coeff = 1; allocate(t_arr(1)%exps(6)); t_arr(1)%exps = [1, 1, 1, 1, 1, 1]
        t_arr(2)%coeff = 6; allocate(t_arr(2)%exps(6)); t_arr(2)%exps = [0, 0, 0, 0, 0, 0]
        allocate(polynomials(6)%terms(2))
        do i = 1, 2; polynomials(6)%terms(i) = t_arr(i); end do; call sort_terms(polynomials(6)%terms)

    end if

    do k = 1, 10
        reduced_basis = grobner_basis(polynomials)
        print *, "Iteration ", k-1, " complete"
    end do

    print *, "Basis size: ", size(reduced_basis)
    do i = 1, size(reduced_basis)
        print *, "Polynomial ", i, ":"
        do j = 1, size(reduced_basis(i)%terms)
            print *, "  Coeff: ", reduced_basis(i)%terms(j)%coeff, " Exps: ", reduced_basis(i)%terms(j)%exps
        end do
    end do

end program main
