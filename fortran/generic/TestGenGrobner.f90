program TestGenGrobner
    use IntModP_Module
    use VecExponent_Module
    use BitPackedExponent_Module
    use GenGrobner_Module
    implicit none

    integer :: n, i
    character(len=32) :: arg
    integer :: mode = 0
    
    n = 4

    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) n
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg, *) mode
    end if

    TERM_ORDER = 0 ! Lex
    IntModP_modulus = 7 ! Ensure mod is 7

    if (mode == 0) then
        print *, "Fortran generic vec exponent cyclic ", n
        call run_benchmark_vec(n)
    else
        print *, "Fortran generic bit packed exponent cyclic ", n
        call run_benchmark_bit(n)
    end if

contains

    subroutine run_benchmark_vec(n)
        integer, intent(in) :: n
        type(polynomial_vec), allocatable :: sys(:), basis(:)
        integer :: iter
        type(term_vec), allocatable :: t(:)
        
        allocate(sys(n))

        if (n == 4) then
            allocate(t(4))
            t(1) = term_vec(IntModP(1_8), VecExponent([1, 0, 0, 0]))
            t(2) = term_vec(IntModP(1_8), VecExponent([0, 1, 0, 0]))
            t(3) = term_vec(IntModP(1_8), VecExponent([0, 0, 1, 0]))
            t(4) = term_vec(IntModP(1_8), VecExponent([0, 0, 0, 1]))
            sys(1) = polynomial_vec(t)
            
            t(1) = term_vec(IntModP(1_8), VecExponent([1, 1, 0, 0]))
            t(2) = term_vec(IntModP(1_8), VecExponent([0, 1, 1, 0]))
            t(3) = term_vec(IntModP(1_8), VecExponent([0, 0, 1, 1]))
            t(4) = term_vec(IntModP(1_8), VecExponent([1, 0, 0, 1]))
            sys(2) = polynomial_vec(t)

            t(1) = term_vec(IntModP(1_8), VecExponent([1, 1, 1, 0]))
            t(2) = term_vec(IntModP(1_8), VecExponent([0, 1, 1, 1]))
            t(3) = term_vec(IntModP(1_8), VecExponent([1, 0, 1, 1]))
            t(4) = term_vec(IntModP(1_8), VecExponent([1, 1, 0, 1]))
            sys(3) = polynomial_vec(t)

            deallocate(t)
            allocate(t(2))
            t(1) = term_vec(IntModP(1_8), VecExponent([1, 1, 1, 1]))
            t(2) = term_vec(IntModP(6_8), VecExponent([0, 0, 0, 0]))
            sys(4) = polynomial_vec(t)
        else if (n == 5) then
            allocate(t(5))
            t(1) = term_vec(IntModP(1_8), VecExponent([1, 0, 0, 0, 0]))
            t(2) = term_vec(IntModP(1_8), VecExponent([0, 1, 0, 0, 0]))
            t(3) = term_vec(IntModP(1_8), VecExponent([0, 0, 1, 0, 0]))
            t(4) = term_vec(IntModP(1_8), VecExponent([0, 0, 0, 1, 0]))
            t(5) = term_vec(IntModP(1_8), VecExponent([0, 0, 0, 0, 1]))
            sys(1) = polynomial_vec(t)
            
            t(1) = term_vec(IntModP(1_8), VecExponent([1, 1, 0, 0, 0]))
            t(2) = term_vec(IntModP(1_8), VecExponent([0, 1, 1, 0, 0]))
            t(3) = term_vec(IntModP(1_8), VecExponent([0, 0, 1, 1, 0]))
            t(4) = term_vec(IntModP(1_8), VecExponent([0, 0, 0, 1, 1]))
            t(5) = term_vec(IntModP(1_8), VecExponent([1, 0, 0, 0, 1]))
            sys(2) = polynomial_vec(t)

            t(1) = term_vec(IntModP(1_8), VecExponent([1, 1, 1, 0, 0]))
            t(2) = term_vec(IntModP(1_8), VecExponent([0, 1, 1, 1, 0]))
            t(3) = term_vec(IntModP(1_8), VecExponent([0, 0, 1, 1, 1]))
            t(4) = term_vec(IntModP(1_8), VecExponent([1, 0, 0, 1, 1]))
            t(5) = term_vec(IntModP(1_8), VecExponent([1, 1, 0, 0, 1]))
            sys(3) = polynomial_vec(t)

            t(1) = term_vec(IntModP(1_8), VecExponent([1, 1, 1, 1, 0]))
            t(2) = term_vec(IntModP(1_8), VecExponent([0, 1, 1, 1, 1]))
            t(3) = term_vec(IntModP(1_8), VecExponent([1, 0, 1, 1, 1]))
            t(4) = term_vec(IntModP(1_8), VecExponent([1, 1, 0, 1, 1]))
            t(5) = term_vec(IntModP(1_8), VecExponent([1, 1, 1, 0, 1]))
            sys(4) = polynomial_vec(t)

            deallocate(t)
            allocate(t(2))
            t(1) = term_vec(IntModP(1_8), VecExponent([1, 1, 1, 1, 1]))
            t(2) = term_vec(IntModP(6_8), VecExponent([0, 0, 0, 0, 0]))
            sys(5) = polynomial_vec(t)
        else
            print *, "Unsupported N"
            return
        end if


        do iter = 1, 10
            basis = naive_grobner_basis_vec(sys)
            print *, "Iteration ", iter-1, " complete"
        end do
    end subroutine run_benchmark_vec

    subroutine run_benchmark_bit(n)
        integer, intent(in) :: n
        type(polynomial_bit), allocatable :: sys(:), basis(:)
        integer :: iter
        type(term_bit), allocatable :: terms_1(:), terms_2(:), terms_3(:), terms_4(:), terms_5(:), terms_6(:)
        
        allocate(sys(n))

        if (n == 4) then
            allocate(terms_1(4))
            terms_1(1) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0001010000000000', 8)))
            terms_1(2) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0001000100000000', 8)))
            terms_1(3) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0001000001000000', 8)))
            terms_1(4) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0001000000010000', 8)))
            sys(1) = polynomial_bit(terms_1)
            allocate(terms_2(4))
            terms_2(1) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0002010100000000', 8)))
            terms_2(2) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0002000101000000', 8)))
            terms_2(3) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0002000001010000', 8)))
            terms_2(4) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0002010000010000', 8)))
            sys(2) = polynomial_bit(terms_2)
            allocate(terms_3(4))
            terms_3(1) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0003010101000000', 8)))
            terms_3(2) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0003000101010000', 8)))
            terms_3(3) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0003010001010000', 8)))
            terms_3(4) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0003010100010000', 8)))
            sys(3) = polynomial_bit(terms_3)
            allocate(terms_4(2))
            terms_4(1) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0004010101010000', 8)))
            terms_4(2) = term_bit(IntModP(6_8), BitPackedExponent(0_8))
            sys(4) = polynomial_bit(terms_4)
        else if (n == 5) then
            allocate(terms_1(5))
            terms_1(1) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0001010000000000', 8)))
            terms_1(2) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0001000100000000', 8)))
            terms_1(3) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0001000001000000', 8)))
            terms_1(4) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0001000000010000', 8)))
            terms_1(5) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0001000000000100', 8)))
            sys(1) = polynomial_bit(terms_1)
            allocate(terms_2(5))
            terms_2(1) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0002010100000000', 8)))
            terms_2(2) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0002000101000000', 8)))
            terms_2(3) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0002000001010000', 8)))
            terms_2(4) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0002000000010100', 8)))
            terms_2(5) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0002010000000100', 8)))
            sys(2) = polynomial_bit(terms_2)
            allocate(terms_3(5))
            terms_3(1) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0003010101000000', 8)))
            terms_3(2) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0003000101010000', 8)))
            terms_3(3) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0003000001010100', 8)))
            terms_3(4) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0003010000010100', 8)))
            terms_3(5) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0003010100000100', 8)))
            sys(3) = polynomial_bit(terms_3)
            allocate(terms_4(5))
            terms_4(1) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0004010101010000', 8)))
            terms_4(2) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0004000101010100', 8)))
            terms_4(3) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0004010001010100', 8)))
            terms_4(4) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0004010100010100', 8)))
            terms_4(5) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0004010101000100', 8)))
            sys(4) = polynomial_bit(terms_4)
            allocate(terms_5(2))
            terms_5(1) = term_bit(IntModP(1_8), BitPackedExponent(int(Z'0005010101010100', 8)))
            terms_5(2) = term_bit(IntModP(6_8), BitPackedExponent(0_8))
            sys(5) = polynomial_bit(terms_5)
        else
            print *, "Unsupported N"
            return
        end if

        do iter = 1, 10
            basis = naive_grobner_basis_bit(sys)
            print *, "Iteration ", iter-1, " complete"
        end do
    end subroutine run_benchmark_bit

end program TestGenGrobner
