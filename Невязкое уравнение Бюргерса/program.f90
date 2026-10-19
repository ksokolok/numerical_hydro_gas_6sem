program main
use burgers
implicit none

real*8, allocatable :: x(:), u(:)
real*8 :: CFL, c, L, b0, b1, Time, h, dt
integer :: Nx, Nt, m, IO = 11
character(len=64) :: filename

open(IO, file='Input.txt')
    READ(IO,*) L
    READ(IO,*) m
    READ(IO,*) b0, b1
    READ(IO,*) Nx
    READ(IO,*) Nt
    READ(IO,*) CFL  
close(IO)

! Расчет по противопоточной схеме (sin(kx), m=4)
call Upwind(x, u, Nx, Nt, CFL, m, L, b0, b1, 1)
write(filename, '(A, I3, A)') 'Burgers_Upwind_sin_m=', m, '.plt'
call Output(trim(filename), x, u)

! Расчет по двухшаговой схеме Лакса-Вендроффа (sin(kx), m=4)
call Lax_Wendroff(x, u, Nx, Nt, CFL, m, L, b0, b1, 1)
write(filename, '(A, I3, A)') 'Burgers_Lax-Wendroff_sin_m=', m, '.plt'
call Output(trim(filename), x, u)


! Расчет по противопоточной схеме (ступенчатая функция 10)
call Upwind(x, u, Nx, Nt, CFL, m, L, b0, b1, 2)
write(filename, '(A, I3, A)') 'Burgers_Upwind_Stepwise_', 10, '.plt'
call Output(trim(filename), x, u)

! Расчет по двухшаговой схеме Лакса-Вендроффа (ступенчатая функция 10)
call Lax_Wendroff(x, u, Nx, Nt, CFL, m, L, b0, b1, 2)
write(filename, '(A, I3, A)') 'Burgers_Lax-Wendroff_Stepwise_', 10, '.plt'
call Output(trim(filename), x, u)
end program