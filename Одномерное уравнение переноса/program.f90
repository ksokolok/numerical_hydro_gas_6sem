program main
use convection
implicit none
real*8, allocatable :: x(:), u(:)
real*8 :: CFL, c, L, b0, b1, Time, h, dt
integer :: Nx, Nt, m, IO = 11
character(len=64) :: filename

open(IO, file='Input.txt')
    READ(IO,*) L
    READ(IO,*) m
    READ(IO,*) c ! Скорость переноса
    READ(IO,*) b0, b1
    READ(IO,*) Nx
    READ(IO,*) Nt
    READ(IO,*) CFL ! Число Куранта      
close(IO)

! Расчет точного решения (sin(kx), m=4)
call Exact_Solution_Sin(x, u, Nx, Nt, CFL, c, m, L, b0, b1)
write(filename, '(A, I3, A)') 'Exact_Solution_m=', m, '.plt'
call Output(trim(filename), x, u, Nx, Nt, CFL, c, L)

! Расчет по противопоточной схеме (sin(kx), m=4)
call Upwind(x, u, Nx, Nt, CFL, c, m, L, b0, b1, 1)
write(filename, '(A, I3, A)') 'Upwind_m=', m, '.plt'
call Output(trim(filename), x, u, Nx, Nt, CFL, c, L)

! Расчет по двухшаговой схеме Лакса-Вендроффа (sin(kx), m=4)
call Lax_Wendroff(x, u, Nx, Nt, CFL, c, m, L, b0, b1, 1)
write(filename, '(A, I3, A)') 'Lax-Wendroff_m=', m, '.plt'
call Output(trim(filename), x, u, Nx, Nt, CFL, c, L)


Nt = 20
h = L / (Nx - 1)
dt = abs(CFL * h / c)
Time = dt * Nt
! Расчет точного решения (ступенчатая функция 10)
m = 10
call Exact_Solution_Stepwise(x, u, 2001, Time, c)
write(filename, '(A, I3, A)') 'Exact_Solution_Stepwise_', m, '.plt'
call Output(trim(filename), x, u, Nx, Nt, CFL, c, L)

! Расчет по противопоточной схеме (ступенчатая функция 10)
call Upwind(x, u, Nx, Nt, CFL, c, m, L, b0, b1, 2)
write(filename, '(A, I3, A)') 'Upwind_Stepwise_', m, '.plt'
call Output(trim(filename), x, u, Nx, Nt, CFL, c, L)

! Расчет по двухшаговой схеме Лакса-Вендроффа (ступенчатая функция 10)
call Lax_Wendroff(x, u, Nx, Nt, CFL, c, m, L, b0, b1, 2)
write(filename, '(A, I3, A)') 'Lax-Wendroff_Stepwise', m, '.plt'
call Output(trim(filename), x, u, Nx, Nt, CFL, c, L)
end program main