program dissipation_surface
use convection
implicit none

real*8, allocatable :: CFL(:), Theta(:)
real*8 :: Theta_min, Theta_max, CFL_min, CFL_max
real*8 :: dTheta, dCFL
integer ::  n = 101, & ! Количество точек по CFL
            m = 101, &  ! Количество точек по фазе Theta
            j, k, IO = 10
character(len=64) :: filename1, filename2


allocate(CFL(1:n), Theta(1:m))

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
do j = 2, n
    Theta(j) = Theta(j-1) + dTheta
end do

! Для явной противопоточной схемы
write(filename1, '(A)') 'Upwind_mod_g.plt'
open(IO, file=trim(filename1))
write(IO, '(A)') 'TITLE = "Upwind |g|(CFL, Theta)"'
write(IO, '(A)') 'VARIABLES = "CFL", "Theta", "|g|"'
write(IO, '(A, I3, A, I3)') 'ZONE T="Upwind |g|", i=', m, ', j=', n

do j = 1, n
do k = 1, m
    write(IO, '(3(f20.10, 10x))') CFL(j), Theta(k), abs_g_Upwind(CFL(j), Theta(k))
end do
end do
close(IO)

! Для схемы Лакса-Вендроффа
write(filename2, '(A)') 'Lax-Wenroff_mod_g.plt'
open(IO, file=trim(filename2))
write(IO, '(A)') 'TITLE = "Lax-Wendroff |g|(CFL, Theta)"'
write(IO, '(A)') 'VARIABLES = "CFL", "Theta", "|g|"'
write(IO, '(A, I3, A, I3)') 'ZONE T="Lax-Wendroff |g|", i=', m, ', j=', n

do j = 1, n
do k = 1, m
    write(IO, '(3(f20.10, 10x))') CFL(j), Theta(k), abs_g_Lax_Wendroff(CFL(j), Theta(k))
end do
end do
close(IO)
end program