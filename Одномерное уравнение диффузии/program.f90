
program main
use diffusion
implicit none

integer 			:: IO = 12
integer 			:: Nr, Nt, j, n, ISCHEME
real*8, allocatable :: r(:), u(:)
real*8 				:: A, rw, Time, VNM, nu

open(IO, file="Input.txt")
read(IO, *) A
read(IO, *) rw
read(IO, *) Nr
read(IO, *) Time
read(IO, *) VNM
read(IO, *) nu
read(IO, *) ISCHEME
close(IO)

call FTCS_Cylinder(r, u, rw, Nr, Nt, VNM, Time, A, nu)

call Output(IO, "Output.txt", r, u)

end program
