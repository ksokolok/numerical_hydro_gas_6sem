module convection
implicit none

real*8, parameter :: pi = acos(-1.d0)

contains



subroutine Upwind(x, u, Nx, Nt, CFL, c, m, L, b0, b1)
    real*8, allocatable, intent(inout) :: x(:), u(:)
    real*8, allocatable :: un(:)
    real*8, intent(in) :: CFL, c, L, b0, b1
    real*8 :: k, dt, h  !, Time
    integer, intent(inout) :: Nx, Nt
    integer, intent(in) :: m
    integer :: j, n

    if (allocated(x)) deallocate(x)
    if (allocated(u)) deallocate(u)

    allocate(x(0:Nx+1), u(0:Nx+1), un(0:Nx+1))

    k = pi * dble(m) / L
    h = L / (Nx - 1)
    dt = CFL * h / c
    ! Time = dt * Nt

    x(0) = -h
    do j = 1, Nx + 1
        x(j) = x(j-1) + h
    end do

    ! Начальные условия
    un(1:Nx) = b0 + b1 * sin(k * x(1:Nx))
    ! Граничные условия периодичности
    un(0) = un(Nx-1)
    un(Nx+1) = un(2)

    u(:) = un(:)

    if (c .gt. 0.d0) then
        do n = 1, Nt
        do j = 1, Nx
            u(j) = un(j) - CFL * (un(j) - un(j-1))
        end do
        u(0) = u(Nx-1)
        u(Nx+1) = u(2)
        un(:) = u(:)
        end do
    end if
    if (c .lt. 0.d0) then
        do n = 1, Nt
        do j = 1, Nx
            u(j) = un(j) - CFL * (un(j+1) - un(j))
        end do
        u(0) = u(Nx-1)
        u(Nx+1) = u(2)
        un(:) = u(:)
        end do
    end if
end subroutine Upwind

subroutine Output(IO, x, u, Nx, Nt, CFL, c, m, L, b0, b1)
    integer, intent(in) :: IO, Nx, Nt, m
    real*8, intent(inout) :: x(0:Nx+1), u(0:Nx+1), CFL, c, L, b0, b1
    real*8 :: h, dt, Time
    integer :: i

    h = L / (Nx - 1)
    dt = CFL * h / c
    Time = dt * Nt

    ! Шапка файла
    write(IO, '(A)') 'TITLE = "Convection"'
    write(IO, '(A)') 'VARIABLES = "x", "u"'
    write(IO, '(A, I2, A, I4)') 'ZONE T = "sin(', m, 'pi / L x)", I =', Nx
    write(IO, '(A, f10.5, A)') 'AUXDATA L = "', L, '"'
    write(IO, '(A, I3, A)') 'AUXDATA m = "', m, '"'
    write(IO, '(A, f17.7, A)') 'AUXDATA CFL = "', CFL, '"'
    write(IO, '(A, f17.7, A)') 'AUXDATA b0 = "', b0, '"'
    write(IO, '(A, f17.7, A)') 'AUXDATA b1 = "', b1, '"'
    write(IO, '(A, I4, A)') 'AUXDATA Nx = "', Nx, '"'
    write(IO, '(A, I4, A)') 'AUXDATA Nt = "', Nt, '"'
    write(IO, '(A, f10.5, A)') 'AUXDATA h = "', h, '"'
    write(IO, '(A, f10.5, A)') 'AUXDATA Time = "', Time, '"'

    do i = 1, size(x) - 1
            write(IO, '(e25.15, 5x, e25.15)') x(i), u(i)
    end do
end subroutine

end module