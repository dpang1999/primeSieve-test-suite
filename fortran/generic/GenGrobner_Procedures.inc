

    function T_TERM_COMPARE(this, other) result(res)
        class(T_TERM), intent(in) :: this
        class(T_TERM), intent(in) :: other
        integer :: res
        integer(8) :: d1, d2
        if (TERM_ORDER == 0) then 
            res = this%exponents%lex_compare(other%exponents)
        else if (TERM_ORDER == 1) then 
            d1 = this%exponents%degree()
            d2 = other%exponents%degree()
            if (d1 /= d2) then
                if (d1 < d2) then
                    res = -1
                else
                    res = 1
                end if
                return
            end if
            res = this%exponents%lex_compare(other%exponents)
        else 
            d1 = this%exponents%degree()
            d2 = other%exponents%degree()
            if (d1 /= d2) then
                if (d1 < d2) then
                    res = -1
                else
                    res = 1
                end if
                return
            end if
            res = -this%exponents%lex_compare(other%exponents)
        end if
    end function T_TERM_COMPARE

    function T_TERM_CAN_REDUCE(this, divisor_leading) result(res)
        class(T_TERM), intent(in) :: this
        class(T_TERM), intent(in) :: divisor_leading
        logical :: res
        res = this%exponents%can_reduce(divisor_leading%exponents)
    end function T_TERM_CAN_REDUCE

    function T_TERM_LCM(this, other) result(res)
        class(T_TERM), intent(in) :: this
        class(T_TERM), intent(in) :: other
        type(E_TYPE) :: res
        res = this%exponents%lcm(other%exponents)
    end function T_TERM_LCM

    function T_TERM_IS_EQUAL(this, other) result(res)
        class(T_TERM), intent(in) :: this
        class(T_TERM), intent(in) :: other
        logical :: res
        res = (this%coefficient == other%coefficient) .and. (this%exponents == other%exponents)
    end function T_TERM_IS_EQUAL

    function T_INIT_POLY(terms_in) result(res)
        type(T_TERM), intent(in) :: terms_in(:)
        type(T_POLY) :: res
        integer :: i, count
        type(T_TERM), allocatable :: temp_terms(:)
        
        allocate(temp_terms(size(terms_in)))
        count = 0
        do i = 1, size(terms_in)
            if (abs(terms_in(i)%coefficient%coerce_to_f64()) > 0.0d0) then
                count = count + 1
                temp_terms(count) = terms_in(i)
            end if
        end do
        
        allocate(res%terms(count))
        if (count > 0) res%terms = temp_terms(1:count)
        
        call T_SORT_TERMS(res%terms)
    end function T_INIT_POLY

    function T_INIT_POLY_EMPTY() result(res)
        type(T_POLY) :: res
        allocate(res%terms(0))
    end function T_INIT_POLY_EMPTY

    subroutine T_SORT_TERMS(arr)
        type(T_TERM), intent(inout) :: arr(:)
        integer :: i, j
        type(T_TERM) :: key
        
        do i = 2, size(arr)
            key = arr(i)
            j = i - 1
            do while (j >= 1)
                if (arr(j)%compare(key) <= 0) exit
                arr(j + 1) = arr(j)
                j = j - 1
            end do
            arr(j + 1) = key
        end do
    end subroutine T_SORT_TERMS

    function T_POLY_ADD(this, other) result(res)
        class(T_POLY), intent(in) :: this
        class(T_POLY), intent(in) :: other
        type(T_POLY) :: res
        type(T_TERM), allocatable :: temp(:)
        logical :: found
        integer :: i, j, count

        allocate(temp(size(this%terms) + size(other%terms)))
        count = size(this%terms)
        if (count > 0) temp(1:count) = this%terms
        
        do i = 1, size(other%terms)
            found = .false.
            do j = 1, count
                if (temp(j)%exponents == other%terms(i)%exponents) then
                    temp(j)%coefficient = temp(j)%coefficient + other%terms(i)%coefficient
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                count = count + 1
                temp(count) = other%terms(i)
            end if
        end do
        
        res = T_POLY(temp(1:count))
    end function T_POLY_ADD

    function T_POLY_SUBTRACT(this, other) result(res)
        class(T_POLY), intent(in) :: this
        class(T_POLY), intent(in) :: other
        type(T_POLY) :: res
        type(T_TERM), allocatable :: temp(:)
        logical :: found
        integer :: i, j, count
        type(C_TYPE) :: zero_coeff

        allocate(temp(size(this%terms) + size(other%terms)))
        count = size(this%terms)
        if (count > 0) temp(1:count) = this%terms
        
        do i = 1, size(other%terms)
            found = .false.
            do j = 1, count
                if (temp(j)%exponents == other%terms(i)%exponents) then
                    temp(j)%coefficient = temp(j)%coefficient - other%terms(i)%coefficient
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                count = count + 1
                zero_coeff = other%terms(i)%coefficient%zero()
                temp(count) = T_TERM(zero_coeff - other%terms(i)%coefficient, other%terms(i)%exponents)
            end if
        end do
        
        res = T_POLY(temp(1:count))
    end function T_POLY_SUBTRACT

    function T_POLY_MAKE_MONIC(this) result(res)
        class(T_POLY), intent(in) :: this
        type(T_POLY) :: res
        type(T_TERM), allocatable :: temp(:)
        type(C_TYPE) :: lead_coeff
        integer :: i

        if (size(this%terms) == 0) then
            res = T_POLY()
            return
        end if

        lead_coeff = this%terms(1)%coefficient
        allocate(temp(size(this%terms)))
        
        do i = 1, size(this%terms)
            temp(i) = T_TERM(this%terms(i)%coefficient / lead_coeff, this%terms(i)%exponents)
        end do
        
        res = T_POLY(temp)
    end function T_POLY_MAKE_MONIC

    function T_POLY_MUL_TERM(this, term) result(res)
        class(T_POLY), intent(in) :: this
        class(T_TERM), intent(in) :: term
        type(T_POLY) :: res
        type(T_TERM), allocatable :: temp(:)
        integer :: i

        allocate(temp(size(this%terms)))
        do i = 1, size(this%terms)
            temp(i) = T_TERM(this%terms(i)%coefficient * term%coefficient, &
                             this%terms(i)%exponents%add(term%exponents))
        end do
        
        res = T_POLY(temp)
    end function T_POLY_MUL_TERM

    function T_POLY_REDUCE(this, divisors) result(res)
        class(T_POLY), intent(in) :: this
        type(T_POLY), intent(in) :: divisors(:)
        type(T_POLY) :: res
        
        type(T_POLY) :: current_poly, scaled_divisor
        type(T_TERM), allocatable :: remainder_terms(:)
        type(T_TERM), allocatable :: temp_terms(:)
        type(T_TERM) :: reduction_term
        integer :: d, count
        logical :: reduced
        type(C_TYPE) :: coeff
        type(E_TYPE) :: exps

        current_poly = this
        allocate(remainder_terms(0))
        count = 0

        do while (size(current_poly%terms) > 0)
            reduced = .false.
            do d = 1, size(divisors)
                if (size(divisors(d)%terms) == 0) cycle
                
                if (current_poly%terms(1)%can_reduce(divisors(d)%terms(1))) then
                    coeff = current_poly%terms(1)%coefficient / divisors(d)%terms(1)%coefficient
                    exps = current_poly%terms(1)%exponents%sub(divisors(d)%terms(1)%exponents)
                    reduction_term = T_TERM(coeff, exps)
                    
                    scaled_divisor = divisors(d)%multiply_by_term(reduction_term)
                    current_poly = current_poly%subtract(scaled_divisor)
                    reduced = .true.
                    exit
                end if
            end do
            
            if (.not. reduced) then
                count = count + 1
                allocate(temp_terms(count))
                if (count > 1) temp_terms(1:count-1) = remainder_terms
                temp_terms(count) = current_poly%terms(1)
                call move_alloc(temp_terms, remainder_terms)
                
                if (size(current_poly%terms) > 1) then
                    allocate(temp_terms(size(current_poly%terms)-1))
                    temp_terms = current_poly%terms(2:)
                    call move_alloc(temp_terms, current_poly%terms)
                else
                    deallocate(current_poly%terms)
                    allocate(current_poly%terms(0))
                end if
            end if
        end do
        
        res = T_POLY(remainder_terms)
    end function T_POLY_REDUCE

    function T_POLY_IS_EQUAL(this, other) result(res)
        class(T_POLY), intent(in) :: this
        class(T_POLY), intent(in) :: other
        logical :: res
        integer :: i

        if (size(this%terms) /= size(other%terms)) then
            res = .false.
            return
        end if

        do i = 1, size(this%terms)
            if (.not. (this%terms(i) == other%terms(i))) then
                res = .false.
                return
            end if
        end do
        res = .true.
    end function T_POLY_IS_EQUAL

    function T_S_POLY(p1, p2) result(res)
        type(T_POLY), intent(in) :: p1
        type(T_POLY), intent(in) :: p2
        type(T_POLY) :: res
        
        type(E_TYPE) :: lcm_exponents, scale_factor_p1, scale_factor_p2
        type(T_POLY) :: scaled_p1, scaled_p2

        lcm_exponents = p1%terms(1)%lcm(p2%terms(1))
        scale_factor_p1 = lcm_exponents%sub(p1%terms(1)%exponents)
        scale_factor_p2 = lcm_exponents%sub(p2%terms(1)%exponents)

        scaled_p1 = p1%multiply_by_term(T_TERM(p2%terms(1)%coefficient, scale_factor_p1))
        scaled_p2 = p2%multiply_by_term(T_TERM(p1%terms(1)%coefficient, scale_factor_p2))

        res = scaled_p1%subtract(scaled_p2)
    end function T_S_POLY

    function T_GROBNER(basis_in) result(res)
        type(T_POLY), intent(in) :: basis_in(:)
        type(T_POLY), allocatable :: res(:)
        
        type(T_POLY), allocatable :: basis(:), basis_excluding_self(:)
        integer, allocatable :: pairs(:,:)
        integer :: i, j, k, p_idx, b_size
        type(T_POLY) :: s_poly, remainder, reduced
        logical :: contains_remainder, contains_reduced
        type(T_POLY), allocatable :: temp_basis(:)
        integer, allocatable :: temp_pairs(:,:)

        allocate(basis(size(basis_in)))
        basis = basis_in

        b_size = size(basis)
        allocate(pairs(2, (b_size * (b_size - 1)) / 2))
        p_idx = 0
        do i = 1, b_size
            do j = i + 1, b_size
                p_idx = p_idx + 1
                pairs(1, p_idx) = i
                pairs(2, p_idx) = j
            end do
        end do

        do while (size(pairs, 2) > 0)
            i = pairs(1, 1)
            j = pairs(2, 1)
            
            if (size(pairs, 2) > 1) then
                allocate(temp_pairs(2, size(pairs, 2) - 1))
                temp_pairs = pairs(:, 2:)
                call move_alloc(temp_pairs, pairs)
            else
                deallocate(pairs)
                allocate(pairs(2, 0))
            end if

            s_poly = T_S_POLY(basis(i), basis(j))
            remainder = s_poly%reduce_poly(basis)
            
            if (size(remainder%terms) > 0) then
                contains_remainder = .false.
                do k = 1, size(basis)
                    if (basis(k) == remainder) then
                        contains_remainder = .true.
                        exit
                    end if
                end do
                
                if (.not. contains_remainder) then
                    b_size = size(basis)
                    allocate(temp_pairs(2, size(pairs, 2) + b_size))
                    if (size(pairs, 2) > 0) temp_pairs(:, 1:size(pairs,2)) = pairs
                    do k = 1, b_size
                        temp_pairs(1, size(pairs, 2) + k) = k
                        temp_pairs(2, size(pairs, 2) + k) = b_size + 1
                    end do
                    call move_alloc(temp_pairs, pairs)
                    
                    allocate(temp_basis(b_size + 1))
                    temp_basis(1:b_size) = basis
                    temp_basis(b_size + 1) = remainder
                    call move_alloc(temp_basis, basis)
                end if
            end if
        end do

        allocate(res(0))
        do i = 1, size(basis)
            allocate(basis_excluding_self(0))
            do j = 1, size(basis)
                if (i /= j) then
                    allocate(temp_basis(size(basis_excluding_self)+1))
                    if (size(basis_excluding_self) > 0) temp_basis(1:size(basis_excluding_self)) = basis_excluding_self
                    temp_basis(size(temp_basis)) = basis(j)
                    call move_alloc(temp_basis, basis_excluding_self)
                end if
            end do
            
            reduced = basis(i)%reduce_poly(basis_excluding_self)
            
            if (size(reduced%terms) > 0) then
                reduced = reduced%make_monic()
                
                contains_reduced = .false.
                do k = 1, size(res)
                    if (res(k) == reduced) then
                        contains_reduced = .true.
                        exit
                    end if
                end do
                
                if (.not. contains_reduced) then
                    allocate(temp_basis(size(res)+1))
                    if (size(res) > 0) temp_basis(1:size(res)) = res
                    temp_basis(size(temp_basis)) = reduced
                    call move_alloc(temp_basis, res)
                end if
            end if
        end do

    end function T_GROBNER
