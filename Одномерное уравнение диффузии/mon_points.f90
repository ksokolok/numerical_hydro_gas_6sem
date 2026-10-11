program mon_points
use diffusion
implicit none

integer, parameter :: IO = 12, IO2 = 13
integer :: Nr, Nt, i, m
real*8, allocatable :: r(:), u(:)
real*8 :: A, rw, VNM, Time, nu
real*8 :: h, dt, u_prev1, u_prev2, u_prev3
character(len=64) :: filename

rw = 1.d0
A = 1.d0
nu = 1.d0

dt = 1d-8      ! Постоянный шаг по времени

Time = 0.1d0     ! Неустановившееся течение

write(filename, '(A, f4.2, A)') 'MonPoints_t=', Time, '.plt' 
open(IO, file=trim(filename))
write(IO, '(A)') 'TITLE = "MonPoints"'
write(IO, '(A)') 'VARIABLES = "h", "VNM", "d r=0.1rw", "d r=0.5rw", "d r=0.9rw"'
write(IO, '(A, f3.1, A)') 'ZONE T="FTCS t=', Time, '"'

Nr = 10
h = rw / Nr
VNM = nu * dt / h**2

call FTCS_Cylinder(r, u, rw, Nr, Nt, VNM, Time, A, nu)
u_prev1 = u(Nr/10)
u_prev2 = u(Nr/2)
u_prev3 = u(9*Nr/10)
print *, 'FTCS, i=', 1, 'done'

m = 6
do i = 2, m
    Nr = 10 * 2**(i - 1)
    h = rw / Nr
    VNM = nu * dt / h**2
    
    call FTCS_Cylinder(r, u, rw, Nr, Nt, VNM, Time, A, nu)
    write(IO, '(5(es25.16, 5x))') h, VNM, abs(u_prev1 - u(Nr/10)), &
        abs(u_prev2 - u(Nr/2)), abs(u_prev3 - u(9*Nr/10))

    u_prev1 = u(Nr/10)
    u_prev2 = u(Nr/2)
    u_prev3 = u(9*Nr/10)
    print *, 'FTCS, i=', i, 'done'
end do

write(IO, '(A, f3.1, A)') 'ZONE T="Alen-Chen t=', Time, '"'

Nr = 10
h = rw / Nr
VNM = nu * dt / h**2

call Alen_Chen_Cylinder(r, u, rw, Nr, Nt, VNM, Time, A, nu)
u_prev1 = u(Nr/10)
u_prev2 = u(Nr/2)
u_prev3 = u(9*Nr/10)
print *, 'Alen-Chen, i=', 1, 'done'

m = 6
do i = 2, m
    Nr = 10 * 2**(i - 1)
    h = rw / Nr
    VNM = nu * dt / h**2
    
    call Alen_Chen_Cylinder(r, u, rw, Nr, Nt, VNM, Time, A, nu)

    write(IO, '(5(es25.16, 5x))') h, VNM, abs(u_prev1 - u(Nr/10)), &
            abs(u_prev2 - u(Nr/2)), abs(u_prev3 - u(9*Nr/10))
    u_prev1 = u(Nr/10)
    u_prev2 = u(Nr/2)
    u_prev3 = u(9*Nr/10)
    print *, 'Alen-Chen, i=', i, 'done'
end do
close(IO)


end program
