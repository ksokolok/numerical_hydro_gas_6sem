program dispersion_surface
use convection
implicit none

real*8, allocatable :: CFL(:), Theta(:), arg_g_ex(:, :)
real*8 :: Theta_min, Theta_max, CFL_min, CFL_max
real*8 :: dTheta, dCFL
integer ::  n = 101, & ! Количество точек по CFL
            m = 101, &  ! Количество точек по фазе Theta
            j, k, IO = 10
character(len=64) :: filename1, filename2


allocate(CFL(1:n), Theta(1:m), arg_g_ex(1:n, 1:m))

Theta_min = 0.d0
Theta_max = pi
dTheta = (Theta_max - Theta_min) / (m - 1)

CFL_min = 0.d0
CFL_max = 2.d0
dCFL = (CFL_max - CFL_min) / (n - 1)

CFL(1) = CFL_min
Theta(1) = Theta_min
do j = 2, n
    CFL(j) = CFL(j-1) + dCFL
end do
do j = 2, m
    Theta(j) = Theta(j-1) + dTheta
end do

do j = 1, n
do k = 1, m
    arg_g_ex(j, k) = - Theta(k) * CFL(j)
end do
end do

! Для явной противопоточной схемы
write(filename1, '(A)') 'Upwind_phi_ratio.plt'
open(IO, file=trim(filename1))
write(IO, '(A)') 'TITLE = "Upwind phi_m/phi_ex(CFL, Theta)"'
write(IO, '(A)') 'VARIABLES = "CFL", "Theta", "phi_m/phi_ex"'
write(IO, '(A, I3, A, I3)') 'ZONE T="Upwind phi_m/phi_ex", i=', m, ', j=', n

do j = 1, n
do k = 1, m
    if (arg_g_ex(j, k) .eq. 0.d0) then
        write(IO, '(3(f20.10, 10x))') CFL(j), Theta(k), 1.d0
    else
    write(IO, '(3(f20.10, 10x))') CFL(j), Theta(k), &
            arg_g_Upwind(CFL(j), Theta(k)) / arg_g_ex(j, k)
    end if
end do
end do
close(IO)

! Для схемы Лакса-Вендроффа
write(filename2, '(A)') 'Lax-Wendroff_phi_ratio.plt'
open(IO, file=trim(filename2))
write(IO, '(A)') 'TITLE = "Lax-Wendroff phi_m/phi_ex(CFL, Theta)"'
write(IO, '(A)') 'VARIABLES = "CFL", "Theta", "phi_m/phi_ex"'
write(IO, '(A, I3, A, I3)') 'ZONE T="Lax-Wendroff phi_m/phi_ex", i=', m, ', j=', n

do j = 1, n
do k = 1, m
    if (arg_g_ex(j, k) .eq. 0.d0) then
        write(IO, '(3(f20.10, 10x))') CFL(j), Theta(k), 1.d0
    else
    write(IO, '(3(f20.10, 10x))') CFL(j), Theta(k), &
            arg_g_Lax_Wendroff(CFL(j), Theta(k)) / arg_g_ex(j, k)
    end if
end do
end do
close(IO)
end program dispersion_surface