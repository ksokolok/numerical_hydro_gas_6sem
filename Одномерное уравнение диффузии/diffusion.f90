module diffusion

contains

subroutine FTCS_Cylinder(r, u, rw, Nr, Nt, VNM, Time, A, nu)
	real*8, intent(inout), allocatable 	:: 	r(:), u(:)
	real*8, allocatable					:: 	un(:)
	real*8, intent(in)					::	rw, VNM, nu, Time, A
	integer, intent(inout)				::  Nr, Nt
	integer								::	j, n
	real*8 								::  dt, h
	
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
		u(j) = un(j) + VNM * ((1 + 1.d0 / (2*j))*un(j+1) - 2.d0 * &
				un(j) + (1 - 1.d0 / (2*j))*un(j-1)) + A * dt
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
end module diffusion
