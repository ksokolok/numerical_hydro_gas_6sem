module convection
implicit none

real*8, parameter :: pi = acos(-1.d0)

contains

subroutine Exact_Solution_Sin(x, u, Nx, Nt, CFL, c, m, L, b0, b1)
    real*8, allocatable, intent(inout) :: x(:), u(:)
    real*8, intent(in) :: CFL, c, L, b0, b1
    real*8 :: k, dt, h, Time
    integer, intent(in) :: Nx, Nt, m
    integer :: j

    if (allocated(x)) deallocate(x)
    if (allocated(u)) deallocate(u)
    allocate(x(0:Nx+1), u(0:Nx+1))

    k = dble(m) * pi / L
    h = L / (Nx - 1)
    dt = abs(CFL * h / c)
    Time = dt * Nt

    x(0) = -h
    do j = 1, Nx + 1
        x(j) = x(j-1) + h
    end do
    u(:) = b0 + b1 * sin(k * (x(:) - c * Time) )

    write(*, '(A)') "Exact Solution Sin"
    write(*, '(A, f8.3)') 't =', Time
    write(*, '(A)') "----------------------------------------------"
end subroutine Exact_Solution_Sin

subroutine Exact_Solution_Stepwise(x, u, Nx, t, c)
    real*8, allocatable, intent(inout) :: x(:), u(:)
    real*8, intent(in) :: t, c
    real*8 :: x1, x2, x3, arg, h, L
    integer, intent(in) :: Nx
    integer :: j

    L = 1.d0

    x1 = 0.2d0
    x2 = 0.4d0
    x3 = 0.6d0
    h = L / (Nx - 1)

    if (allocated(x)) deallocate(x)
    if (allocated(u)) deallocate(u)
    allocate(x(0:Nx+1), u(0:Nx+1))

    x(0) = -h
    do j = 1, Nx  
        x(j) = x(j-1) + h
        arg = modulo(x(j) - c * t, L)
        if (arg .lt. x1) then
            u(j) = 0.d0
        else if (arg .lt. x2) then
            u(j) = 2.d0
        else if (arg .lt. x3) then
            u(j) = 1.d0
        else
            u(j) = 0.d0
        end if
    end do
    u(0) = u(Nx-1)
    u(Nx+1) = u(2)

    write(*, '(A)') "Exact Solution Stepwise"
    write(*, '(A, f8.3)') 't =', t
    write(*, '(A)') "----------------------------------------------"
end subroutine Exact_Solution_Stepwise

subroutine Upwind(x, u, Nx, Nt, CFL, c, m, L, b0, b1, flag)
    real*8, allocatable, intent(inout) :: x(:), u(:)
    real*8, allocatable :: un(:)
    real*8, intent(in) :: CFL, c, L, b0, b1
    real*8 :: k, dt, h, Time
    integer, intent(in) :: m, Nx, Nt
    integer :: j, n, flag

    if (allocated(x)) deallocate(x)
    if (allocated(u)) deallocate(u)

    allocate(x(0:Nx+1), u(0:Nx+1), un(0:Nx+1))

    k = pi * dble(m) / L
    h = L / (Nx - 1)
    dt = abs(CFL * h / c)
    Time = dt * Nt

    x(0) = -h
    do j = 1, Nx + 1
        x(j) = x(j-1) + h
    end do

    ! Начальные условия
    if (flag .eq. 1) then
        un(1:Nx) = b0 + b1 * sin(k * x(1:Nx))
    else if (flag .eq. 2) then
        call Exact_Solution_Stepwise(x, un, Nx, 0.d0, c)
    else
        stop 'Initial values: flag != 1 or 2'
    end if
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
            u(j) = un(j) + CFL * (un(j+1) - un(j))
        end do
        u(0) = u(Nx-1)
        u(Nx+1) = u(2)
        un(:) = u(:)
        end do
    end if

    if (flag .eq. 1) then
        write(*, '(A)') "Upwind Sin"
    else if (flag .eq. 2) then
        write(*, '(A)') "Upwind Stepwise"
    else
        stop 'flag != 1 or 2'
    end if
    write(*, '(2(A, I5, 5x))') "Nx =", Nx, "Nt =", Nt
    write(*, '(4(A, f8.3, 5x))') "CFL =", CFL, "dt =", dt, 'h =', h, 't =', Time
    write(*, '(A)') "----------------------------------------------"
end subroutine Upwind

subroutine Lax_Wendroff(x, u, Nx, Nt, CFL, c, m, L, b0, b1, flag)
    real*8, allocatable, intent(inout) :: x(:), u(:)
    real*8, allocatable :: un(:), uhalf(:)
    real*8, intent(in) :: CFL, c, L, b0, b1
    real*8 :: k, dt, h, Time
    integer, intent(in) :: m, Nx, Nt, flag
    integer :: j, n

    if (allocated(x)) deallocate(x)
    if (allocated(u)) deallocate(u)
    if (allocated(uhalf)) deallocate(uhalf)

    ! uhalf(j) соответствует u_{j+1/2}^{n+1/2}
    allocate(x(0:Nx+1), u(0:Nx+1), un(0:Nx+1), uhalf(0:Nx+1))

    k = pi * dble(m) / L
    h = L / (Nx - 1)
    dt = abs(CFL * h / c)
    Time = dt * Nt

    x(0) = -h
    do j = 1, Nx + 1
        x(j) = x(j-1) + h
    end do

    ! Начальные условия
    if (flag .eq. 1) then
        un(1:Nx) = b0 + b1 * sin(k * x(1:Nx))
    else if (flag .eq. 2) then
        call Exact_Solution_Stepwise(x, un, Nx, 0.d0, c)
    else
        stop 'Initial values: flag != 1 or 2'
    end if
    ! Граничные условия периодичности
    un(0) = un(Nx-1)
    un(Nx+1) = un(2)

    u(:) = un(:)
    uhalf(:) = 0.d0

    if (c .gt. 0.d0) then
        do n = 1, Nt
        ! Первый шаг: u_{j+1/2}^{n+1/2}
        do j = 0, Nx
            uhalf(j) = 0.5d0 * (un(j+1) + un(j)) &
                     - 0.5d0 * CFL * (un(j+1) - un(j))
        end do
        ! Второй шаг
        do j = 1, Nx
            u(j) = un(j) - CFL * (uhalf(j) - uhalf(j-1))
        end do
        u(0) = u(Nx-1)
        u(Nx+1) = u(2)
        un(:) = u(:)
        end do
    end if

    if (c .lt. 0.d0) then
        do n = 1, Nt
        ! Первый шаг: u_{j+1/2}^{n+1/2}
        do j = 0, Nx
            uhalf(j) = 0.5d0 * (un(j+1) + un(j)) &
                     + 0.5d0 * CFL * (un(j+1) - un(j))
        end do
        ! Второй шаг
        do j = 1, Nx
            u(j) = un(j) + CFL * (uhalf(j) - uhalf(j-1))
        end do
        u(0) = u(Nx-1)
        u(Nx+1) = u(2)
        un(:) = u(:)
        end do
    end if
    if (flag .eq. 1) then
        write(*, '(A)') "Lax-Wendroff Sin"
    else if (flag .eq. 2) then
        write(*, '(A)') "Lax-Wendroff Stepwise"
    else
        stop 'flag != 1 or 2'
    end if
    write(*, '(2(A, I5, 5x))') "Nx =", Nx, "Nt =", Nt
    write(*, '(4(A, f8.3, 5x))') "CFL =", CFL, "dt =", dt, 'h =', h, 't =', Time
    write(*, '(A)') "----------------------------------------------"
end subroutine Lax_Wendroff

subroutine Output(filename, x, u)
    character(len=*) :: filename
    real*8, intent(inout) :: x(:), u(:)
    integer :: i, IO = 10

    open(IO, file=filename)
    ! Шапка файла
    write(IO, '(A)') 'TITLE = "Convection"'
    write(IO, '(A)') 'VARIABLES = "x", "u"'
    write(IO, '(3A, I4)') 'ZONE T = "', filename,'", I =', size(x) - 2

    do i = 2, size(x) - 1
            write(IO, '(e25.15, 5x, e25.15)') x(i), u(i)
    end do
    close(IO)
end subroutine

real*8 function abs_g_Upwind(CFL, Theta)
    real*8 :: CFL, Theta, Re, Im

    Re = 1.d0 - CFL * (1.d0 - cos(Theta))
    Im = - CFL * sin(Theta)
    abs_g_Upwind = sqrt(Re**2 + Im**2)
end function

real*8 function abs_g_Lax_Wendroff(CFL, Theta)
    real*8 :: CFL, Theta, Re, Im

    Re = 1.d0 + CFL**2 * (cos(Theta) - 1.d0)
    Im = - CFL * sin(Theta)
    abs_g_Lax_Wendroff = sqrt(Re**2 + Im**2)
end function 

real*8 function arg_g_Upwind(CFL, Theta)
    real*8 :: CFL, Theta, Re, Im

    Re = 1.d0 - CFL * (1.d0 - cos(Theta))
    Im = - CFL * sin(Theta)

    arg_g_Upwind = atan2( Im, Re )
end function

real*8 function arg_g_Lax_Wendroff(CFL, Theta)
    real*8 :: CFL, Theta, Re, Im

    Re = 1.d0 - CFL**2 * (1.d0 - cos(Theta))
    Im = - CFL * sin(Theta)

    arg_g_Lax_Wendroff = atan2( Im, Re )
end function

subroutine Estimation_Upwind_Sin(x, u, Nx, Nt, CFL, m, L, b0, b1)
    real*8, allocatable, intent(inout) :: x(:), u(:)
    real*8, intent(in) :: CFL, L, b0, b1
    real*8 :: k, h, Theta  !, Time
    integer, intent(in) :: m, Nx, Nt
    integer :: j

    if (allocated(x)) deallocate(x)
    if (allocated(u)) deallocate(u)

    allocate(x(0:Nx+1), u(0:Nx+1))

    k = dble(m) * pi / L
    h = L / (Nx - 1)

    x(0) = -h
    do j = 1, Nx + 1
        x(j) = x(j-1) + h
    end do

    k = pi * m / L
    Theta = k * h
    u(:) = b0 + b1 * abs_g_Upwind(CFL, Theta)**Nt &
            * (sin(k * x(:) + arg_g_Upwind(CFL, Theta) * Nt))
end subroutine 

subroutine Estimation_Lax_Wendroff_Sin(x, u, Nx, Nt, CFL, m, L, b0, b1)
    real*8, allocatable, intent(inout) :: x(:), u(:)
    real*8, intent(in) :: CFL, L, b0, b1
    real*8 :: k, h, Theta  !, Time
    integer, intent(in) :: m, Nx, Nt
    integer :: j

    if (allocated(x)) deallocate(x)
    if (allocated(u)) deallocate(u)

    allocate(x(0:Nx+1), u(0:Nx+1))

    k = dble(m) * pi / L
    h = L / (Nx - 1)

    x(0) = -h
    do j = 1, Nx + 1
        x(j) = x(j-1) + h
    end do

    k = pi * m / L
    Theta = k * h
    u(:) = b0 + b1 * abs_g_Lax_Wendroff(CFL, Theta)**Nt &
            * (sin(k * x(:) + arg_g_Lax_Wendroff(CFL, Theta) * Nt))
end subroutine 


end module