program Alen_Chen_Stability
	integer :: j, i, k, IO = 10
	integer, parameter :: n = 100, &! Количество отрезков по s
						  m = 100 ! Количество отрезков по Theta
	real*8, parameter :: pi = acos(-1.d0)
	real*8 :: Theta_min, Theta_max, s_min, s_max
	real*8, allocatable :: s(:), &  ! числа фон Неймана
						   Theta(:) ! фаза k_m * dr

	
	
	j = 1000	! Узел, в котором считаем
	allocate(s(0:n), Theta(0:m))
	
	Theta_min = 0.d0
	Theta_max = pi
	dTheta = (Theta_max - Theta_min) / m
	
	s_min = 0.1**4
	s_max = 100.d0
	dLogs = (log10(s_max) - log10(s_min)) / n
	
	Theta(0) = Theta_min
	s(0) = s_min
	do k = 1, 100
		Theta(k) = Theta(k-1) + dTheta
		s(k) = 10.d0**(log10(s_min) + k * dLogs)
	end do
	
	open(IO, file="Alen-Chen_Stability.plt")
	! Шапка файла с результатами
	write(IO, '(A)') 'TITLE = "Alen-Chen Stability"'
	write(IO, '(A)') 'VARIABLES = "s", "Theta", "|g|"'
	write(IO, '(A, I3, A, I3, A, I3)') 'ZONE T="|g|, j=', j, '", i=', n, ', j=', m
	
	! Вычисления и запись в файл
	do i = 1, n
		do k = 1, m
			write(10, '(3(f20.10, 10x))') s(i), Theta(k), g_Alen_Chen(j, s(i), Theta(k))
		end do
	end do
	close(IO)
	
contains
	
	real*8 function g_Alen_Chen(i, s, Theta)
		integer :: i
		real*8 :: j, s, Theta
		real*8 :: c1, c2, c3
		
		j = dble(i)
		c1 = 2.d0 - 1.d0 / (2.d0 * (j**2 - 1.d0))
		c2 = 2.d0 + 1.d0 / (2.d0 * (j**2 - 1.d0))
		c3 = (4.d0 * j**2 - 3.d0) / (2.d0 * j * (j**2 - 1.d0))
		
		g_Alen_Chen = sqrt( &
						((s**2 * c1 * cos(2.d0 * Theta) + 2.d0 * s * cos(Theta) &
						+ 1.d0 + 2.d0 * s + s**2 * c2)**2 &
						+ (s**2 * c3 * sin(2.d0 * Theta) + s / j * sin(Theta))**2) &
					  ) / (1.d0 + 2.d0 * s)**2
	end function
end program
