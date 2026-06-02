module PolyMod
    implicit none
    private
    public :: MODULUS, Term, Polynomial, &
              compare_terms, add_poly, sub_poly, multiply_by_term, &
              leading_term, get_lcm, can_reduce, s_polynomial, &
              make_monic, reduce_polynomial, polynomials_equal, &
              grobner_basis, pack_exponents, init_polynomial

    integer(8), save :: MODULUS = 7

    type :: Term
        integer(8) :: coefficient
        integer(8) :: exponents
    end type Term

    type :: Polynomial
        type(Term), allocatable :: terms(:)
    end type Polynomial

contains

    function mod_inverse(a, m) result(res)
        integer(8), intent(in) :: a, m
        integer(8) :: res
        integer(8) :: a_i64, m_i64, m0, q, t, x, y
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

    function compare_terms(t1, t2) result(res)
        type(Term), intent(in) :: t1, t2
        integer :: res
        integer(8) :: e1, e2
        e1 = iand(t1%exponents, Z'0000FFFFFFFFFFFF')
        e2 = iand(t2%exponents, Z'0000FFFFFFFFFFFF')
        if (e1 < e2) then
            res = -1
        else if (e1 > e2) then
            res = 1
        else
            res = 0
        end if
    end function compare_terms

    subroutine sort_terms(terms)
        type(Term), intent(inout) :: terms(:)
        integer :: n, i, j
        type(Term) :: temp
        n = size(terms)
        do i = 1, n - 1
            do j = i + 1, n
                if (compare_terms(terms(i), terms(j)) < 0) then
                    temp = terms(i)
                    terms(i) = terms(j)
                    terms(j) = temp
                end if
            end do
        end do
    end subroutine sort_terms

    function init_polynomial(terms_in) result(p)
        type(Term), intent(in) :: terms_in(:)
        type(Polynomial) :: p
        integer :: i, count, n
        n = size(terms_in)
        count = 0
        do i = 1, n
            if (terms_in(i)%coefficient /= 0) then
                count = count + 1
            end if
        end do
        allocate(p%terms(count))
        count = 0
        do i = 1, n
            if (terms_in(i)%coefficient /= 0) then
                count = count + 1
                p%terms(count) = terms_in(i)
            end if
        end do
        call sort_terms(p%terms)
    end function init_polynomial

    function add_poly(p1, p2) result(p)
        type(Polynomial), intent(in) :: p1, p2
        type(Polynomial) :: p
        type(Term), allocatable :: temp_terms(:)
        integer :: i, j, n1, n2, count
        logical :: found

        if (.not. allocated(p1%terms)) then
            n1 = 0
        else
            n1 = size(p1%terms)
        end if
        if (.not. allocated(p2%terms)) then
            n2 = 0
        else
            n2 = size(p2%terms)
        end if

        allocate(temp_terms(n1 + n2))
        
        count = 0
        do i = 1, n1
            count = count + 1
            temp_terms(count) = p1%terms(i)
        end do

        do i = 1, n2
            found = .false.
            do j = 1, count
                if (temp_terms(j)%exponents == p2%terms(i)%exponents) then
                    temp_terms(j)%coefficient = mod(temp_terms(j)%coefficient + p2%terms(i)%coefficient, MODULUS)
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                count = count + 1
                temp_terms(count) = p2%terms(i)
            end if
        end do

        p = init_polynomial(temp_terms(1:count))
        deallocate(temp_terms)
    end function add_poly

    function sub_poly(p1, p2) result(p)
        type(Polynomial), intent(in) :: p1, p2
        type(Polynomial) :: p
        type(Term), allocatable :: temp_terms(:)
        integer :: i, j, n1, n2, count
        logical :: found

        if (.not. allocated(p1%terms)) then
            n1 = 0
        else
            n1 = size(p1%terms)
        end if
        if (.not. allocated(p2%terms)) then
            n2 = 0
        else
            n2 = size(p2%terms)
        end if

        allocate(temp_terms(n1 + n2))
        
        count = 0
        do i = 1, n1
            count = count + 1
            temp_terms(count) = p1%terms(i)
        end do

        do i = 1, n2
            found = .false.
            do j = 1, count
                if (temp_terms(j)%exponents == p2%terms(i)%exponents) then
                    temp_terms(j)%coefficient = mod(MODULUS + temp_terms(j)%coefficient - p2%terms(i)%coefficient, MODULUS)
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                count = count + 1
                temp_terms(count)%coefficient = mod(MODULUS - mod(p2%terms(i)%coefficient, MODULUS), MODULUS)
                temp_terms(count)%exponents = p2%terms(i)%exponents
            end if
        end do

        p = init_polynomial(temp_terms(1:count))
        deallocate(temp_terms)
    end function sub_poly

    function multiply_by_term(p1, t) result(p)
        type(Polynomial), intent(in) :: p1
        type(Term), intent(in) :: t
        type(Polynomial) :: p
        type(Term), allocatable :: temp_terms(:)
        integer :: i, n

        if (.not. allocated(p1%terms)) then
            n = 0
        else
            n = size(p1%terms)
        end if

        allocate(temp_terms(n))
        do i = 1, n
            temp_terms(i)%coefficient = mod(p1%terms(i)%coefficient * t%coefficient, MODULUS)
            temp_terms(i)%exponents = p1%terms(i)%exponents + t%exponents
        end do

        p = init_polynomial(temp_terms)
        deallocate(temp_terms)
    end function multiply_by_term

    function leading_term(p) result(t)
        type(Polynomial), intent(in) :: p
        type(Term) :: t
        if (allocated(p%terms)) then
            if (size(p%terms) > 0) then
                t = p%terms(1)
                return
            end if
        end if
        t%coefficient = 0
        t%exponents = 0
    end function leading_term

    function get_lcm(a, b) result(res)
        integer(8), intent(in) :: a, b
        integer(8) :: res
        integer(8) :: self_exponents, other_exponents, lcm_exponents, deg
        integer(8) :: self_exp, other_exp, lcm_exp
        integer :: i

        self_exponents = iand(a, Z'0000FFFFFFFFFFFF')
        other_exponents = iand(b, Z'0000FFFFFFFFFFFF')
        lcm_exponents = 0
        deg = 0

        do i = 0, 40, 8
            self_exp = iand(ishft(self_exponents, -i), Z'FF')
            other_exp = iand(ishft(other_exponents, -i), Z'FF')
            lcm_exp = max(self_exp, other_exp)
            lcm_exponents = ior(lcm_exponents, ishft(lcm_exp, i))
            deg = deg + lcm_exp
        end do
        lcm_exponents = ior(lcm_exponents, ishft(iand(deg, Z'FFFF'), 48))
        res = lcm_exponents
    end function get_lcm

    function can_reduce(a, b) result(res)
        integer(8), intent(in) :: a, b
        logical :: res
        integer(8) :: self_exp, divisor_exp
        integer :: i

        res = .true.
        do i = 0, 40, 8
            self_exp = iand(ishft(a, -i), Z'FF')
            divisor_exp = iand(ishft(b, -i), Z'FF')
            if (self_exp < divisor_exp) then
                res = .false.
                return
            end if
        end do
    end function can_reduce

    function s_polynomial(p1, p2) result(p)
        type(Polynomial), intent(in) :: p1, p2
        type(Polynomial) :: p
        type(Term) :: lt1, lt2, s1_term, s2_term
        integer(8) :: lcm_exps, s1_exp, s2_exp

        if (.not. allocated(p1%terms) .or. .not. allocated(p2%terms)) then
            allocate(p%terms(0))
            return
        end if
        if (size(p1%terms) == 0 .or. size(p2%terms) == 0) then
            allocate(p%terms(0))
            return
        end if

        lt1 = leading_term(p1)
        lt2 = leading_term(p2)

        lcm_exps = get_lcm(lt1%exponents, lt2%exponents)
        s1_exp = lcm_exps - lt1%exponents
        s2_exp = lcm_exps - lt2%exponents

        s1_term%coefficient = lt2%coefficient
        s1_term%exponents = s1_exp

        s2_term%coefficient = lt1%coefficient
        s2_term%exponents = s2_exp

        p = sub_poly(multiply_by_term(p1, s1_term), multiply_by_term(p2, s2_term))
    end function s_polynomial

    function make_monic(p1) result(p)
        type(Polynomial), intent(in) :: p1
        type(Polynomial) :: p
        integer(8) :: lt_coeff, inv_coeff
        integer :: i, n
        type(Term) :: temp_t

        if (.not. allocated(p1%terms)) then
            allocate(p%terms(0))
            return
        end if
        n = size(p1%terms)
        if (n == 0) then
            allocate(p%terms(0))
            return
        end if

        temp_t = leading_term(p1)
        lt_coeff = temp_t%coefficient
        inv_coeff = mod_inverse(lt_coeff, MODULUS)

        allocate(p%terms(n))
        do i = 1, n
            p%terms(i)%coefficient = mod(p1%terms(i)%coefficient * inv_coeff, MODULUS)
            p%terms(i)%exponents = p1%terms(i)%exponents
        end do
    end function make_monic

    function reduce_polynomial(p1, basis) result(p)
        type(Polynomial), intent(in) :: p1
        type(Polynomial), intent(in) :: basis(:)
        type(Polynomial) :: p, poly, scaled_b
        type(Term), allocatable :: remainder(:), temp_remainder(:)
        type(Term) :: lt_poly, lt_b, scale_term
        logical :: reduced
        integer :: i, n_basis, r_count, p_size, k
        integer(8) :: inv_b, scale_coeff, scale_exps
        type(Term), allocatable :: new_poly_terms(:)

        poly = p1
        n_basis = size(basis)
        r_count = 0
        allocate(remainder(0))

        do while (allocated(poly%terms) .and. size(poly%terms) > 0)
            reduced = .false.
            lt_poly = leading_term(poly)

            do i = 1, n_basis
                if (.not. allocated(basis(i)%terms)) cycle
                if (size(basis(i)%terms) == 0) cycle
                
                lt_b = leading_term(basis(i))
                if (can_reduce(lt_poly%exponents, lt_b%exponents)) then
                    inv_b = mod_inverse(lt_b%coefficient, MODULUS)
                    scale_coeff = mod(lt_poly%coefficient * inv_b, MODULUS)
                    scale_exps = lt_poly%exponents - lt_b%exponents
                    
                    scale_term%coefficient = scale_coeff
                    scale_term%exponents = scale_exps

                    scaled_b = multiply_by_term(basis(i), scale_term)
                    poly = sub_poly(poly, scaled_b)
                    reduced = .true.
                    exit
                end if
            end do

            if (.not. reduced) then
                if (allocated(remainder)) then
                    allocate(temp_remainder(r_count + 1))
                    if (r_count > 0) temp_remainder(1:r_count) = remainder
                    temp_remainder(r_count + 1) = leading_term(poly)
                    call move_alloc(temp_remainder, remainder)
                    r_count = r_count + 1
                else
                    allocate(remainder(1))
                    remainder(1) = leading_term(poly)
                    r_count = 1
                end if
                
                p_size = size(poly%terms)
                if (p_size > 1) then
                    allocate(new_poly_terms(p_size - 1))
                    do k = 2, p_size
                        new_poly_terms(k - 1) = poly%terms(k)
                    end do
                    call move_alloc(new_poly_terms, poly%terms)
                else
                    if (allocated(poly%terms)) deallocate(poly%terms)
                    allocate(poly%terms(0))
                end if
            end if
        end do

        if (r_count > 0) then
            p = init_polynomial(remainder)
        else
            allocate(p%terms(0))
        end if
    end function reduce_polynomial

    function polynomials_equal(p1, p2) result(res)
        type(Polynomial), intent(in) :: p1, p2
        logical :: res
        integer :: i, n1, n2

        if (.not. allocated(p1%terms)) then
            n1 = 0
        else
            n1 = size(p1%terms)
        end if
        if (.not. allocated(p2%terms)) then
            n2 = 0
        else
            n2 = size(p2%terms)
        end if

        if (n1 /= n2) then
            res = .false.
            return
        end if

        res = .true.
        do i = 1, n1
            if (p1%terms(i)%coefficient /= p2%terms(i)%coefficient .or. p1%terms(i)%exponents /= p2%terms(i)%exponents) then
                res = .false.
                return
            end if
        end do
    end function polynomials_equal

    function contains_poly(arr, p) result(res)
        type(Polynomial), intent(in) :: arr(:)
        type(Polynomial), intent(in) :: p
        logical :: res
        integer :: i, n
        n = size(arr)
        res = .false.
        do i = 1, n
            if (polynomials_equal(arr(i), p)) then
                res = .true.
                return
            end if
        end do
    end function contains_poly

    subroutine append_poly(arr, p)
        type(Polynomial), allocatable, intent(inout) :: arr(:)
        type(Polynomial), intent(in) :: p
        type(Polynomial), allocatable :: temp(:)
        integer :: n
        if (allocated(arr)) then
            n = size(arr)
            allocate(temp(n + 1))
            temp(1:n) = arr
            temp(n + 1) = p
            call move_alloc(temp, arr)
        else
            allocate(arr(1))
            arr(1) = p
        end if
    end subroutine append_poly

    function grobner_basis(basis_in) result(reduced_basis)
        type(Polynomial), intent(in) :: basis_in(:)
        type(Polynomial), allocatable :: reduced_basis(:)
        type(Polynomial), allocatable :: basis(:), basis_set(:)
        integer, allocatable :: pair_set(:,:)
        integer, allocatable :: temp_pair_set(:,:)
        integer :: n_basis, n_pairs, i, j, k
        type(Polynomial) :: s_poly, h, reduced
        type(Polynomial), allocatable :: basis_excluding_self(:)
        
        n_basis = size(basis_in)
        allocate(basis(n_basis))
        basis = basis_in

        allocate(basis_set(0))
        do i = 1, n_basis
            if (.not. contains_poly(basis_set, basis(i))) then
                call append_poly(basis_set, basis(i))
            end if
        end do

        n_pairs = (n_basis * (n_basis - 1)) / 2
        allocate(pair_set(n_pairs, 2))
        k = 1
        do i = 1, n_basis
            do j = i + 1, n_basis
                pair_set(k, 1) = i
                pair_set(k, 2) = j
                k = k + 1
            end do
        end do

        do while (size(pair_set, 1) > 0)
            i = pair_set(1, 1)
            j = pair_set(1, 2)
            
            if (size(pair_set, 1) > 1) then
                allocate(temp_pair_set(size(pair_set, 1) - 1, 2))
                temp_pair_set = pair_set(2:, :)
                call move_alloc(temp_pair_set, pair_set)
            else
                deallocate(pair_set)
                allocate(pair_set(0, 2))
            end if
            
            s_poly = s_polynomial(basis(i), basis(j))
            h = reduce_polynomial(s_poly, basis)
            
            if (allocated(h%terms)) then
                if (size(h%terms) > 0 .and. .not. contains_poly(basis_set, h)) then
                    call append_poly(basis_set, h)
                    n_basis = size(basis)
                    call append_poly(basis, h)
                    
                    if (size(pair_set, 1) == 0) then
                        deallocate(pair_set)
                        allocate(pair_set(n_basis, 2))
                        do k = 1, n_basis
                            pair_set(k, 1) = k
                            pair_set(k, 2) = n_basis + 1
                        end do
                    else
                        allocate(temp_pair_set(size(pair_set, 1) + n_basis, 2))
                        temp_pair_set(1:size(pair_set, 1), :) = pair_set
                        do k = 1, n_basis
                            temp_pair_set(size(pair_set, 1) + k, 1) = k
                            temp_pair_set(size(pair_set, 1) + k, 2) = n_basis + 1
                        end do
                        call move_alloc(temp_pair_set, pair_set)
                    end if
                end if
            end if
        end do

        allocate(reduced_basis(0))
        n_basis = size(basis)
        do i = 1, n_basis
            allocate(basis_excluding_self(n_basis - 1))
            k = 1
            do j = 1, n_basis
                if (i /= j) then
                    basis_excluding_self(k) = basis(j)
                    k = k + 1
                end if
            end do
            
            reduced = reduce_polynomial(basis(i), basis_excluding_self)
            deallocate(basis_excluding_self)
            
            if (allocated(reduced%terms)) then
                if (size(reduced%terms) > 0) then
                    call append_poly(reduced_basis, make_monic(reduced))
                end if
            end if
        end do

    end function grobner_basis

    function pack_exponents(exps) result(res)
        integer, intent(in) :: exps(:)
        integer(8) :: res
        integer(8) :: packed, shift_val, e_val, deg
        integer :: i, n
        integer :: full_exps(6)

        n = size(exps)
        full_exps = 0
        if (n > 6) n = 6
        full_exps(1:n) = exps(1:n)

        packed = 0
        deg = 0
        do i = 1, 6
            shift_val = 40 - 8 * (i - 1)
            e_val = iand(int(full_exps(i), 8), Z'FF')
            packed = ior(packed, ishft(e_val, int(shift_val)))
            deg = deg + full_exps(i)
        end do
        packed = ior(packed, ishft(iand(deg, Z'FFFF'), 48))
        res = packed
    end function pack_exponents

end module PolyMod

program main
    use PolyMod
    implicit none

    integer :: n
    character(len=32) :: arg
    integer :: i, j
    type(Polynomial), allocatable :: polynomials(:), reduced_basis(:)
    type(Term), allocatable :: terms(:)

    n = 4
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        read(arg, *) n
    end if

    MODULUS = 7

    if (n == 4) then
        write(*,*) "Fortran specialized bit-packed cyclic 4"
        allocate(polynomials(4))
        
        allocate(terms(4))
        terms(1) = Term(1_8, pack_exponents([1, 0, 0, 0]))
        terms(2) = Term(1_8, pack_exponents([0, 1, 0, 0]))
        terms(3) = Term(1_8, pack_exponents([0, 0, 1, 0]))
        terms(4) = Term(1_8, pack_exponents([0, 0, 0, 1]))
        polynomials(1) = init_polynomial(terms); deallocate(terms)
        
        allocate(terms(4))
        terms(1) = Term(1_8, pack_exponents([1, 1, 0, 0]))
        terms(2) = Term(1_8, pack_exponents([0, 1, 1, 0]))
        terms(3) = Term(1_8, pack_exponents([0, 0, 1, 1]))
        terms(4) = Term(1_8, pack_exponents([1, 0, 0, 1]))
        polynomials(2) = init_polynomial(terms); deallocate(terms)

        allocate(terms(4))
        terms(1) = Term(1_8, pack_exponents([1, 1, 1, 0]))
        terms(2) = Term(1_8, pack_exponents([0, 1, 1, 1]))
        terms(3) = Term(1_8, pack_exponents([1, 0, 1, 1]))
        terms(4) = Term(1_8, pack_exponents([1, 1, 0, 1]))
        polynomials(3) = init_polynomial(terms); deallocate(terms)

        allocate(terms(2))
        terms(1) = Term(1_8, pack_exponents([1, 1, 1, 1]))
        terms(2) = Term(6_8, pack_exponents([0, 0, 0, 0]))
        polynomials(4) = init_polynomial(terms); deallocate(terms)
    else if (n == 5) then
        write(*,*) "Fortran specialized bit-packed cyclic 5"
        allocate(polynomials(5))

        allocate(terms(5))
        terms(1) = Term(1_8, pack_exponents([1, 0, 0, 0, 0]))
        terms(2) = Term(1_8, pack_exponents([0, 1, 0, 0, 0]))
        terms(3) = Term(1_8, pack_exponents([0, 0, 1, 0, 0]))
        terms(4) = Term(1_8, pack_exponents([0, 0, 0, 1, 0]))
        terms(5) = Term(1_8, pack_exponents([0, 0, 0, 0, 1]))
        polynomials(1) = init_polynomial(terms); deallocate(terms)

        allocate(terms(5))
        terms(1) = Term(1_8, pack_exponents([1, 1, 0, 0, 0]))
        terms(2) = Term(1_8, pack_exponents([0, 1, 1, 0, 0]))
        terms(3) = Term(1_8, pack_exponents([0, 0, 1, 1, 0]))
        terms(4) = Term(1_8, pack_exponents([0, 0, 0, 1, 1]))
        terms(5) = Term(1_8, pack_exponents([1, 0, 0, 0, 1]))
        polynomials(2) = init_polynomial(terms); deallocate(terms)

        allocate(terms(5))
        terms(1) = Term(1_8, pack_exponents([1, 1, 1, 0, 0]))
        terms(2) = Term(1_8, pack_exponents([0, 1, 1, 1, 0]))
        terms(3) = Term(1_8, pack_exponents([0, 0, 1, 1, 1]))
        terms(4) = Term(1_8, pack_exponents([1, 0, 0, 1, 1]))
        terms(5) = Term(1_8, pack_exponents([1, 1, 0, 0, 1]))
        polynomials(3) = init_polynomial(terms); deallocate(terms)

        allocate(terms(5))
        terms(1) = Term(1_8, pack_exponents([1, 1, 1, 1, 0]))
        terms(2) = Term(1_8, pack_exponents([0, 1, 1, 1, 1]))
        terms(3) = Term(1_8, pack_exponents([1, 0, 1, 1, 1]))
        terms(4) = Term(1_8, pack_exponents([1, 1, 0, 1, 1]))
        terms(5) = Term(1_8, pack_exponents([1, 1, 1, 0, 1]))
        polynomials(4) = init_polynomial(terms); deallocate(terms)

        allocate(terms(2))
        terms(1) = Term(1_8, pack_exponents([1, 1, 1, 1, 1]))
        terms(2) = Term(6_8, pack_exponents([0, 0, 0, 0, 0]))
        polynomials(5) = init_polynomial(terms); deallocate(terms)
    else if (n == 6) then
        write(*,*) "Fortran specialized bit-packed cyclic 6"
        allocate(polynomials(6))

        allocate(terms(6))
        terms(1) = Term(1_8, pack_exponents([1, 0, 0, 0, 0, 0]))
        terms(2) = Term(1_8, pack_exponents([0, 1, 0, 0, 0, 0]))
        terms(3) = Term(1_8, pack_exponents([0, 0, 1, 0, 0, 0]))
        terms(4) = Term(1_8, pack_exponents([0, 0, 0, 1, 0, 0]))
        terms(5) = Term(1_8, pack_exponents([0, 0, 0, 0, 1, 0]))
        terms(6) = Term(1_8, pack_exponents([0, 0, 0, 0, 0, 1]))
        polynomials(1) = init_polynomial(terms); deallocate(terms)

        allocate(terms(6))
        terms(1) = Term(1_8, pack_exponents([1, 1, 0, 0, 0, 0]))
        terms(2) = Term(1_8, pack_exponents([0, 1, 1, 0, 0, 0]))
        terms(3) = Term(1_8, pack_exponents([0, 0, 1, 1, 0, 0]))
        terms(4) = Term(1_8, pack_exponents([0, 0, 0, 1, 1, 0]))
        terms(5) = Term(1_8, pack_exponents([0, 0, 0, 0, 1, 1]))
        terms(6) = Term(1_8, pack_exponents([1, 0, 0, 0, 0, 1]))
        polynomials(2) = init_polynomial(terms); deallocate(terms)

        allocate(terms(6))
        terms(1) = Term(1_8, pack_exponents([1, 1, 1, 0, 0, 0]))
        terms(2) = Term(1_8, pack_exponents([0, 1, 1, 1, 0, 0]))
        terms(3) = Term(1_8, pack_exponents([0, 0, 1, 1, 1, 0]))
        terms(4) = Term(1_8, pack_exponents([0, 0, 0, 1, 1, 1]))
        terms(5) = Term(1_8, pack_exponents([1, 0, 0, 0, 1, 1]))
        terms(6) = Term(1_8, pack_exponents([1, 1, 0, 0, 0, 1]))
        polynomials(3) = init_polynomial(terms); deallocate(terms)

        allocate(terms(6))
        terms(1) = Term(1_8, pack_exponents([1, 1, 1, 1, 0, 0]))
        terms(2) = Term(1_8, pack_exponents([0, 1, 1, 1, 1, 0]))
        terms(3) = Term(1_8, pack_exponents([0, 0, 1, 1, 1, 1]))
        terms(4) = Term(1_8, pack_exponents([1, 0, 0, 1, 1, 1]))
        terms(5) = Term(1_8, pack_exponents([1, 1, 0, 0, 1, 1]))
        terms(6) = Term(1_8, pack_exponents([1, 1, 1, 0, 0, 1]))
        polynomials(4) = init_polynomial(terms); deallocate(terms)

        allocate(terms(6))
        terms(1) = Term(1_8, pack_exponents([1, 1, 1, 1, 1, 0]))
        terms(2) = Term(1_8, pack_exponents([0, 1, 1, 1, 1, 1]))
        terms(3) = Term(1_8, pack_exponents([1, 0, 1, 1, 1, 1]))
        terms(4) = Term(1_8, pack_exponents([1, 1, 0, 1, 1, 1]))
        terms(5) = Term(1_8, pack_exponents([1, 1, 1, 0, 1, 1]))
        terms(6) = Term(1_8, pack_exponents([1, 1, 1, 1, 0, 1]))
        polynomials(5) = init_polynomial(terms); deallocate(terms)

        allocate(terms(2))
        terms(1) = Term(1_8, pack_exponents([1, 1, 1, 1, 1, 1]))
        terms(2) = Term(6_8, pack_exponents([0, 0, 0, 0, 0, 0]))
        polynomials(6) = init_polynomial(terms); deallocate(terms)
    end if

    do i = 1, 10
        if (allocated(reduced_basis)) deallocate(reduced_basis)
        reduced_basis = grobner_basis(polynomials)
        write(*, "('Iteration ', I0, ' complete')") i - 1
    end do

    print *, "Basis size: ", size(reduced_basis)
    do i = 1, size(reduced_basis)
        print *, "Polynomial ", i, ":"
        do j = 1, size(reduced_basis(i)%terms)
            print *, "  Coeff: ", reduced_basis(i)%terms(j)%coefficient, " Exps: ", reduced_basis(i)%terms(j)%exponents
        end do
    end do

end program main