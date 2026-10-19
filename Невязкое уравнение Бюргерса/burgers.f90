module burgers

real*8, parameter :: pi = acos(-1.d0)

contains

subroutine Exact_Solution_Stepwise()

end subroutine

subroutine Upwind(x, u, Nx, Nt, CFL, m, L, b0, b1, flag)
    real*8, allocatable, intent(inout) :: x(:), u(:)
    real*8, allocatable :: un(:)
    real*8, intent(in) :: CFL, L, b0, b1
    real*8 :: c, k, dt, h, Time
    integer, intent(in) :: m, Nx, Nt
    integer :: j, n, flag

    if (allocated(x)) deallocate(x)
    if (allocated(u)) deallocate(u)

    allocate(x(0:Nx+1), u(0:Nx+1), un(0:Nx+1))

    k = pi * dble(m) / L
    h = L / (Nx - 1)

    x(0) = -h
    do j = 1, Nx + 1
        x(j) = x(j-1) + h
    end do

    if (flag .eq. 1) then
        un(1:Nx) = b0 + b1 * sin(k * x(1:Nx))
    else if (flag .eq. 2) then
        do j = 1, Nx
            if (x(j) .lt. 0.2d0) then
                un(j) = 1.d0
            else if (x(j) .lt. 0.4d0) then
                un(j) = 3.d0
            else if (x(j) .lt. 0.6d0) then
                un(j) = 2.d0
            else
                un(j) = 1.d0
            end if
        end do
    else
        stop 'Initial values: flag != 1 or 2'
    end if
    un(0) = un(Nx-1)
    un(Nx+1) = un(2)

    ! Максимальная скорость используется для расчета шага по времени
    ! Предполагается, что скорость не увеличится, хотя что если осцилляции?
    c = maxval(abs(un(:))) 
    !c = maxval(un(:)) - minval(un(:))

    dt = abs(CFL * h / c)
    Time = dt * Nt

    u(:) = un(:)

    do n = 1, Nt
    do j = 1, Nx
        if ((un(j) + un(j-1)) .gt. 0.d0) then
            u(j) = un(j) - dt / (2.d0 * h) * (un(j)**2 - un(j-1)**2)
        else
            u(j) = un(j) - dt / (2.d0 * h) * (un(j+1)**2 - un(j)**2)
        end if
    end do
    u(0) = u(Nx-1)
    u(Nx+1) = u(2)
    un(:) = u(:)
    end do

    if (flag .eq. 1) then
        write(*, '(A)') "Upwind Sin"
    else if (flag .eq. 2) then
        write(*, '(A)') "Upwind Stepwise"
    else
        stop 'flag != 1 or 2'
    end if
    write(*, '(2(A, I5, 5x))') "Nx =", Nx, "Nt =", Nt
    write(*, '(5(A, f8.3, 5x))') "CFL =", CFL, "dt =", dt, 'h =', h, 'c =', c, 't =', Time
    write(*, '(A)') "----------------------------------------------"
end subroutine Upwind

subroutine Lax_Wendroff(x, u, Nx, Nt, CFL, m, L, b0, b1, flag)
    real*8, allocatable, intent(inout) :: x(:), u(:)
    real*8, allocatable :: un(:), uhalf(:)
    real*8, intent(in) :: CFL, L, b0, b1
    real*8 :: c, k, dt, h, Time
    integer, intent(in) :: m, Nx, Nt, flag
    integer :: j, n

    if (allocated(x)) deallocate(x)
    if (allocated(u)) deallocate(u)

    allocate(x(0:Nx+1), u(0:Nx+1), un(0:Nx+1), uhalf(0:Nx+1))

    k = pi * dble(m) / L
    h = L / (Nx - 1)

    x(0) = -h
    do j = 1, Nx + 1
        x(j) = x(j-1) + h
    end do

    if (flag .eq. 1) then
        un(1:Nx) = b0 + b1 * sin(k * x(1:Nx))
    else if (flag .eq. 2) then
        do j = 1, Nx
            if (x(j) .lt. 0.2d0) then
                un(j) = 1.d0
            else if (x(j) .lt. 0.4d0) then
                un(j) = 3.d0
            else if (x(j) .lt. 0.6d0) then
                un(j) = 2.d0
            else
                un(j) = 1.d0
            end if
        end do
    else
        stop 'Initial values: flag != 1 or 2'
    end if
    un(0) = un(Nx-1)
    un(Nx+1) = un(2)

    ! Максимальная скорость используется для расчета шага по времени
    ! Предполагается, что скорость не увеличится, хотя что если осцилляции?
    !c = maxval(un(:)) - minval(un(:))
    c = maxval(abs(un(:)))

    dt = abs(CFL * h / c)
    Time = dt * Nt

    u(:) = un(:)
    uhalf(:) = 0.d0

    do n = 1, Nt
    ! Первый шаг
    do j = 0, Nx
        uhalf(j) = 0.5d0 * (un(j+1) + un(j)) &
                 - dt / (4.d0 * h) * (un(j+1)**2 - un(j)**2)
    end do
    ! Второй шаг
    do j = 1, Nx
        u(j) = un(j) - dt / (2.d0 * h) * (uhalf(j)**2 - uhalf(j-1)**2)
    end do
    u(0) = u(Nx-1)
    u(Nx+1) = u(2)
    un(:) = u(:)
    end do

    if (flag .eq. 1) then
        write(*, '(A)') "Lax-Wendroff Sin"
    else if (flag .eq. 2) then
        write(*, '(A)') "Lax-Wendroff Stepwise"
    else
        stop 'flag != 1 or 2'
    end if
    write(*, '(2(A, I5, 5x))') "Nx =", Nx, "Nt =", Nt
    write(*, '(5(A, f8.3, 5x))') "CFL =", CFL, "dt =", dt, 'h =', h, 'c =', c, 't =', Time
    write(*, '(A)') "----------------------------------------------"
end subroutine Lax_Wendroff



subroutine Output(filename, x, u)
    character(len=*) :: filename
    real*8, intent(inout) :: x(:), u(:)
    integer :: i, IO = 10

    open(IO, file=filename)
    ! Шапка файла
    write(IO, '(A)') 'TITLE = "Burgers"'
    write(IO, '(A)') 'VARIABLES = "x", "u"'
    write(IO, '(3A, I4)') 'ZONE T = "', filename,'", I =', size(x) - 2

    do i = 2, size(x) - 1
            write(IO, '(e25.15, 5x, e25.15)') x(i), u(i)
    end do
    close(IO)
end subroutine

subroutine Exact_Solution_Stepwise10(x, u, t, Nx)
    real*8, allocatable, intent(inout) :: x(:), u(:)
    real*8, intent(in) :: t
    real*8 :: L
    integer, intent(in) :: Nx
    integer :: j

    if (allocated(x)) deallocate(x)
    if (allocated(u)) deallocate(u)

    allocate(x(0:Nx+1), u(0:Nx+1))

    ! Решение, вычисляемое здесь справделиво только для t < 0.2
    if (t .gt. 0.2d0) stop 'Riemann_Stepwise10: t > 0.2'

    L = 1.d0
    h = L / dble(Nx - 1)

    x(0) = -h
    do j = 1, Nx + 1
        x(j) = x(j-1) + h
    end do

    do j = 1, Nx
        if (x(j) .le. t + 0.2d0) then
            u(j) = 1.d0
        else if (x(j) .lt. 3.d0 * t + 0.2d0) then
            u(j) = (x(j) - 0.2d0) / t
        else if (x(j) .lt. 2.5d0 * t + 0.4d0) then
            u(j) = 3.d0
        else if (x(j) .lt. 1.5d0 * t + 0.6d0) then
            u(j) = 2.d0
        else
            u(j) = 1.d0
        end if
    end do
    u(0) = u(Nx-1)
    u(Nx+1) = u(2)

    write(*, '(A)') 'Exact Solution Stepwise'
    write(*, '(A, f8.3)') 't =', t
    write(*, '(A)') '----------------------------------------------'
end subroutine


end module