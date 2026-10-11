
program main
use diffusion
implicit none

integer 			:: IO = 12
integer 			:: Nr, Nt, j, n, ISCHEME
real*8, allocatable :: r(:), u(:)
real*8 				:: A, rw, VNM, Time, nu
real*8				:: A1, A2
integer				:: NA

open(IO, file="Input.txt")
read(IO, *) A
read(IO, *) rw
read(IO, *) Nr
read(IO, *) Time
read(IO, *) VNM
read(IO, *) nu
read(IO, *) ISCHEME
close(IO)

! Получение решения (ВВЦП)
call FTCS_Cylinder(r, u, rw, Nr, Nt, VNM, Time, A, nu)
call Output(IO, "Output.txt", r, u)

! Получение решения для разных моментов времени (ВВЦП)
call FTCS_Cylinder_Evolution(r, u, rw, Nr, Nt, VNM, Time, A, nu, &
								'snapshots.plt', 10)

! Получение решения для разных градиентов давления	(ВВЦП)						
Time = 1
A1 = 1
A2 = 10
NA = 100
call FTCS_Cylinder_A(r, u, rw, Nr, Nt, VNM, Time, A1, A2, NA, nu,&
									'umax(A).plt')

end program
