! gfortran (file name) -o (name of executable)
! -o option specifies name of output executable file, else default is a.out
module sieve 

contains
function primeSieve(num) result(primes) 
    implicit none
    integer(kind=8), intent(in) :: num
    logical, dimension(num) :: primes
    integer(kind=8) :: i, j, current

    primes = .true.
    primes(1) = .false.
    i = 2
    do while (i < num)
        if (primes(i)) then
            j = i
            current = i * j
            do while (current < num)
                primes(current) = .false.
                j = j + 1
                current = i * j
            end do
        end if
        i = i + 1
    end do
end function primeSieve
end module sieve

program hello
    use sieve
    implicit none
    integer(kind=8) :: max, i
    logical, dimension(:), allocatable :: temp

    character(100) :: arg1
    integer(kind=8) :: arg1int

    max = 42
    if(command_argument_count() > 0) then
        call get_command_argument(1, arg1)
        read(arg1, *) arg1int
        max = arg1int
    end if
    
    allocate(temp(max))
    temp = primeSieve(max)
    i = 2
    do while (i < size(temp)) 
        if (temp(i)) then
            print *, i
        end if
        i = i + 1
    end do
    
    deallocate(temp)
end program hello