program schemes_research
use diffusion
implicit none

integer 			:: IO = 12, IO2 = 13
integer 			:: Nr, Nt, i, j, m
real*8, allocatable :: r(:), u(:), us(:), ue(:)
real*8 				:: A, rw, VNM, Time, nu
real*8				:: h, dt
real*8 				:: Times(3) = [0.1d0, 0.7d0, 5.d0]
character(len=*), parameter	:: filename11 = 'FTCS_dt_t=0.1.plt', &
							   filename12 = 'FTCS_dt_t=0.7,5.0.plt', &
							   filename21 = 'Alen-Chen_dt_t=0.1.plt', &
							   filename22 = 'Alen-Chen_dt_t=0.7,5.0.plt'

! Исследование влияния шага по времени на решение
rw = 1.d0
Nr = 100
A = 1.d0
nu = 1.d0
! Получение аналитического решения us(:) для устойчивого течения (стационара)
!	Используется для выражения ошибки при t=0.7, 5.0
call Analytical_Solution(r, us, rw, Nr, A, nu)

h = rw / Nr

! ВВЦП
Time = Times(1) ! t = 0.1
open(IO, file=filename11)
write(IO, '(A)') 'TITLE = "FTCS dt"'
write(IO, '(A)') 'VARIABLES = "VNM", "dt", "Err Ref"'
write(IO, '(A, f3.1, A)') 'ZONE T="FTCS t=', Time, '"'

! Эталонное решение (для выражения ошибки)
VNM = dble(1e-5)
dt = VNM * h**2.d0 / nu
call Alen_Chen_Cylinder(r, ue, rw, Nr, Nt, VNM, Time, A, nu)

VNM = 0.4d0
dt = VNM * h**2.d0 / nu
m = 8
do j = 1, m
	call FTCS_Cylinder(r, u, rw, Nr, Nt, VNM, Time, A, nu)
	write(IO, '(3(e20.10, 10x))') VNM, dt, maxval(abs(u(:) - ue(:)))
	VNM = VNM / 2.d0
	dt = VNM * h**2.d0 / nu
end do
print *, 'FTCS, t=', Time, ' done'
close(IO)

! t = 0.7, 5.0
open(IO2, file=filename12)
write(IO2, '(A)') 'TITLE = "FTCS dt"'
write(IO2, '(A)') 'VARIABLES = "VNM", "dt", "Err Exact"'
do i = 2, 3
	Time = Times(i)

	write(IO2, '(A, f3.1, A)') 'ZONE T="FTCS t=', Time, '"'
	VNM = 0.4d0
	dt = VNM * h**2.d0 / nu
	do j = 1, m
		call FTCS_Cylinder(r, u, rw, Nr, Nt, VNM, Time, A, nu)
		write(IO2, '(4(e20.10, 10x))') VNM, dt, maxval(abs(u(:) - us(:)))
		VNM = VNM / 2.d0
		dt = VNM * h**2.d0 / nu
	end do
	
	print *, 'FTCS, t=', Time, ' done'
end do
close(IO2)

! Ален-Чен
Time = Times(1) ! t = 0.1
open(IO, file=filename21)
write(IO, '(A)') 'TITLE = "Alen-Chen dt"'
write(IO, '(A)') 'VARIABLES = "VNM", "dt", "Err Ref"'
write(IO, '(A, f3.1, A)') 'ZONE T="Alen-Chen t=', Time, '"'

VNM = 0.4d0
dt = VNM * h**2.d0 / nu
m = 12
do j = 1, m
	call Alen_Chen_Cylinder(r, u, rw, Nr, Nt, VNM, Time, A, nu)
	write(IO, '(3(e20.10, 10x))') VNM, dt, maxval(abs(u(:) - ue(:)))
	VNM = VNM / 2.d0
	dt = VNM * h**2.d0 / nu
end do
print *, 'Alen-Chen, t=', Time, ' done'
close(IO)

! t = 0.7, 5.0
open(IO2, file=filename22)
write(IO2, '(A)') 'TITLE = "Alen-Chen dt"'
write(IO2, '(A)') 'VARIABLES = "VNM", "dt", "Err Exact"'
do i = 2, 3
	Time = Times(i)

	write(IO2, '(A, f3.1, A)') 'ZONE T="Alen-Chen t=', Time, '"'
	VNM = 0.4d0
	dt = VNM * h**2.d0 / nu			
	do j = 1, m
		call Alen_Chen_Cylinder(r, u, rw, Nr, Nt, VNM, Time, A, nu)
		write(IO2, '(3(e20.10, 10x))') VNM, dt, maxval(abs(u(:) - us(:)))
		VNM = VNM / 2.d0
		dt = VNM * h**2.d0 / nu
	end do
	
	print *, 'Alen-Chen, t=', Time, ' done'
end do
close(IO2)

end program schemes_research
