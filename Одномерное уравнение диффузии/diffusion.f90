module diffusion

contains

subroutine Analytical_Solution(r, u, rw, Nr, A, nu)
    real*8, intent(inout), allocatable 	:: 	r(:), u(:)
    real*8, intent(in)					::	rw, A, nu
    integer, intent(inout)				::  Nr
    integer								::	j
    real*8 								::  h
    
    if (allocated(r)) deallocate(r)
    if (allocated(u)) deallocate(u)
    allocate(r(0:Nr), u(0:Nr))
    
    h = rw / Nr
    do j = 0, Nr
        r(j) = h * j
        u(j) = A / (4.d0 * nu) * (rw**2.d0 - r(j)**2.d0)
    end do
end subroutine

subroutine FTCS_Cylinder(r, u, rw, Nr, Nt, VNM, Time, A, nu)
    real*8, intent(inout), allocatable 	:: 	r(:), u(:)
    real*8, allocatable					:: 	un(:)
    real*8, intent(in)					::	rw, VNM, nu, A, Time
    integer, intent(inout)				::  Nr, Nt
    integer								::	j, n
    real*8 								::  dt, h
    
    if (allocated(r)) deallocate(r)
    if (allocated(u)) deallocate(u)
    if (allocated(un)) deallocate(un)
    allocate(r(0:Nr), u(0:Nr), un(0:Nr))
    
    h = rw / Nr
    dt = VNM * h**2 / nu
    Nt = Time / dt
    
    ! Сетка
    r(0) = 0.d0
    do j = 1, Nr - 1
        r(j) = r(j - 1) + h
    end do
    r(Nr) = rw
    
    ! Начальные значения
    u(:) = 0.d0
    un(:) = 0.d0
    
    ! Граничные условия
    un(Nr) = 0.d0
    u(Nr) = 0.d0
    
    ! Решение
    do n = 1, Nt
    do j = 1, Nr - 1
        u(j) = un(j) + VNM * ((1.d0 + h / (2.d0 * r(j)))*un(j+1) - 2.d0 * &
                un(j) + (1.d0 - h / (2.d0 * r(j)))*un(j-1)) + A * dt
    end do
    u(0) = un(0) + 4.d0 * VNM * (un(1) - un(0)) + A * dt
    un(:) = u(:)
    end do
end subroutine FTCS_Cylinder

subroutine Output(IO, filename, r, u)
    integer, intent(in) :: IO
    real*8, intent(inout) :: r(0:), u(0:)
    character(len=*) :: filename
    integer :: i
    
    open(IO, file=filename)
    do i = 0, size(r) - 1
        write(IO, '(I7, 5x, e25.15, 5x, e25.15)') i, r(i), u(i)
    end do
    close(IO)
end subroutine Output

subroutine FTCS_Cylinder_Evolution(r, u, rw, Nr, Nt, VNM, Time, A, nu, &
                                    filename, Nsnapshots)
    real*8, intent(inout), allocatable 	:: 	r(:), u(:) 
    real*8, allocatable					:: 	un(:)
    real*8, intent(in)					::	rw, VNM, nu, A, Time
    integer, intent(inout)				::  Nr, Nt
    integer, intent(in)					::	Nsnapshots
    integer								::	j, n, IO, &
                                            dn_snapshots
    real*8 								::  dt, h
    character(len=*), intent(in)		:: filename
    
    if (allocated(r)) deallocate(r)
    if (allocated(u)) deallocate(u)
    if (allocated(un)) deallocate(un)
    allocate(r(0:Nr), u(0:Nr), un(0:Nr))
    
    h = rw / Nr
    dt = VNM * h**2 / nu
    Nt = Time / dt
    
    ! Сетка
    r(0) = 0.d0
    do j = 1, Nr - 1
        r(j) = r(j - 1) + h
    end do
    r(Nr) = rw
    
    ! Начальные значения
    u(:) = 0.d0
    un(:) = 0.d0
    
    ! Граничные условия
    un(Nr) = 0.d0
    u(Nr) = 0.d0
    
    ! Шапка выходного файла
    IO = 12
    open(IO, file=filename)
    write(IO, '(A)') 'TITLE = "Evolution of u(r,t)"'
    write(IO, '(A)') 'VARIABLES = "r", "u"'
    write(IO, '(A,f9.4,A)') 'ZONE T="t=', 0.d0, '"'
    
    ! Интервал (в индексах) между "снимками"
    dn_snapshots = max(Nt / Nsnapshots, 1)
    
    ! Запись начального профиля как первой зоны
    write(IO_snap, '(A,f7.2,A)') 'ZONE T="t =', 0.d0, '"'
    do j = 0, Nr
        write(IO, '(2(E20.10,5x))') r(j), u(j)
    end do   
    
    ! Решение
    do n = 1, Nt
    do j = 1, Nr - 1
        u(j) = un(j) + VNM * ((1.d0 + h / (2.d0 * r(j)))*un(j+1) - 2.d0 * &
                un(j) + (1.d0 - h / (2.d0 * r(j)))*un(j-1)) + A * dt
    end do
    u(0) = un(0) + 4.d0 * VNM * (un(1) - un(0)) + A * dt
    un(:) = u(:)	
    
    ! Запись профиля
    if (mod(n, dn_snapshots) == 0) then
        write(IO, '(A,f7.2,A)') 'ZONE T="t =', n*dt, '"'
        do j = 0, Nr
            write(IO, '(2(E20.10,5x))') r(j), u(j)
        end do
    end if	
    
    end do
    close(IO)
end subroutine FTCS_Cylinder_Evolution							

subroutine FTCS_Cylinder_A(r, u, rw, Nr, Nt, VNM, Time, A1, A2, NA, nu,&
                                    filename)
    real*8, intent(inout), allocatable 	:: 	r(:), u(:) 
    real*8, allocatable					:: 	un(:)
    real*8, intent(in)					::	rw, VNM, nu, A1, A2, Time
    integer, intent(inout)				::  Nr, Nt, NA
    integer								::	j, n, IO, i
    real*8 								::  dt, h, A, dA
    character(len=*), intent(in)		:: filename
    
    if (allocated(r)) deallocate(r)
    if (allocated(u)) deallocate(u)
    if (allocated(un)) deallocate(un)
    allocate(r(0:Nr), u(0:Nr), un(0:Nr))
    
    h = rw / Nr
    dt = VNM * h**2 / nu
    Nt = Time / dt
    
    ! Сетка
    r(0) = 0.d0
    do j = 1, Nr - 1
        r(j) = r(j - 1) + h
    end do
    r(Nr) = rw
    
    ! Начальные значения
    u(:) = 0.d0
    un(:) = 0.d0
    
    ! Граничные условия
    un(Nr) = 0.d0
    u(Nr) = 0.d0
    
    ! Шапка выходного файла
    IO = 12
    open(IO, file=filename)
    write(IO, '(A)') 'TITLE = "Umax(A)"'
    write(IO, '(A)') 'VARIABLES = "u_max", "A"'
    write(IO, '(A)') 'ZONE T="Umax(A)"'
    
    dA = (A2 - A1) / (NA - 1)
    do i = 1, NA
    
    A = A1 + (i - 1) * dA
    ! Решение
    do n = 1, Nt
    do j = 1, Nr - 1
        u(j) = un(j) + VNM * ((1.d0 + h / (2.d0 * r(j)))*un(j+1) - 2.d0 * &
                un(j) + (1.d0 - h / (2.d0 * r(j)))*un(j-1)) + A * dt
    end do
    u(0) = un(0) + 4.d0 * VNM * (un(1) - un(0)) + A * dt
    un(:) = u(:)
    end do
    write(IO, '(2(E20.10,5x))') u(0), A
    end do
end subroutine FTCS_Cylinder_A



subroutine Alen_Chen_Cylinder(r, u, rw, Nr, Nt, VNM, Time, A, nu)
    real*8, intent(inout), allocatable 	:: 	r(:), u(:)
    real*8, allocatable					:: 	un(:), up(:)
    real*8, intent(in)					::	rw, VNM, nu, A, Time
    integer, intent(inout)				::  Nr, Nt
    integer								::	j, n
    real*8 								::  dt, h
    
    if (allocated(r)) deallocate(r)
    if (allocated(u)) deallocate(u)
    if (allocated(un)) deallocate(un)
    if (allocated(up)) deallocate(up)
    allocate(r(0:Nr), u(0:Nr), un(0:Nr), up(0:Nr))
    
    h = rw / Nr
    dt = VNM * h**2 / nu
    Nt = Time / dt
    
    ! Сетка
    r(0) = 0.d0
    do j = 1, Nr - 1
        r(j) = r(j - 1) + h
    end do
    r(Nr) = rw
    
    ! Начальные значения
    u(:) = 0.d0
    un(:) = 0.d0
    up(:) = 0.d0
    
    ! Граничные условия
    un(Nr) = 0.d0
    u(Nr) = 0.d0
    
    ! Решение
    do n = 1, Nt
    
    ! Первый шаг (предиктор)
    do j = 1, Nr - 1
        up(j) = (1.d0 + 2.d0 * VNM)**(-1.d0) * &
                (un(j) + VNM * ((1.d0 + h / (2.d0 * r(j))) * un(j+1) + &
                                (1.d0 - h / (2.d0 * r(j))) * un(j-1)) + &
                         A * dt)
    end do
    up(0) = (1.d0 + 4.d0 * VNM)**(-1.d0) * (un(0) + 4.d0 * VNM * un(1) &
                                            + A * dt)
    up(Nr) = 0.d0
    
    ! Второй шаг (корректор)
    do j = 1, Nr - 1
        u(j) = (1.d0 + 2.d0 * VNM)**(-1.d0) * &
                (un(j) + VNM * ((1.d0 + h / (2.d0 * r(j))) * up(j+1) + &
                                (1.d0 - h / (2.d0 * r(j))) * up(j-1)) + &
                         A * dt)
    end do
    u(0) = (1.d0 + 4.d0 * VNM)**(-1.d0) * (un(0) + 4.d0 * VNM * up(1) &
                                            + A * dt)
    u(Nr) = 0.d0
    
    un(:) = u(:)
    
    end do
end subroutine Alen_Chen_Cylinder


end module diffusion
