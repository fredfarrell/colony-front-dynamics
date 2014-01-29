!Code to solve the standard Fisher equation in two dimensions, with noise
!Equation is dphi/dt = alpha phi (1-phi) + D grad^2 phi + sqrt(phi)*noise
!Uses Dickmann discretization scheme

include 'shift.f90'

program fisher_noise

	integer, parameter :: SIZE=250
	real(8) :: psi(SIZE,SIZE), dpsi(SIZE,SIZE) !psi is a continuous field keeping track of phi
	integer :: n(SIZE,SIZE) !n is the discretized version of phi
	real(8) :: D,dt,dx,alpha
	real(8) :: phi_min, noise_strength, sqrtdt
	integer :: i,j,t,maxt,iup,idwn
	character(len=1024) :: filename

	!parameters
	D=1
	alpha=1
	dt=0.001
	dx=0.1
	maxt=5000
	phi_min=0.001
	noise_strength=10;
	sqrtdt=sqrt(dt)

	!initialize random number generator
	call srand(123)

	!initialize
	do i = 1,SIZE
		do j = 1,SIZE
			n(i,j)=(0.5*tanh( (SIZE/4-j)/10.0) + 0.5)/phi_min
			psi(i,j)=0
		enddo
	enddo

	!timestep
	do t = 1, maxt

		call shift_down_integer(n,0,1,SIZE,100)

		write(6,*) t

		do i = 1,SIZE
			do j=2, SIZE-1

				iup = merge(1,i+1,i==SIZE)
				idwn = merge(SIZE,i-1,i==1)

				!calculate change in psi over the timestep
				dpsi(i,j) = dt*alpha*n(i,j)*(1-phi_min*n(i,j)) &
				 + (D*dt)/(dx*dx)*(n(iup,j) + n(idwn,j) + n(i,j+1) + n(i,j-1) - 4*n(i,j)) &
				 + noise_strength*sqrt(real(n(i,j)))*sqrtdt*(rand()-0.5)

			enddo
		enddo

		!simulataneous update of psi array
		do i = 1,SIZE
			do j=2,SIZE-1
				psi(i,j) = psi(i,j) + dpsi(i,j)

				! add integer part of psi to n
				n(i,j) = n(i,j) + int(psi(i,j))
				psi(i,j) = psi(i,j) - int(psi(i,j))

			enddo
		enddo

		!zero flux BC in the j-direction
		do i = 1,SIZE
			psi(i,1)=psi(i,2)
			psi(i,SIZE-1)=psi(i,SIZE)
		enddo

		!print data every 1000 time steps
		if(mod(t,1000)==0) then

		write(filename,"(A3,I0,A4)") 'out', t, '.dat'
			open(1,file=filename)

			do i = 1,SIZE
				do j=2,SIZE-1
					write(1,*) i,j,psi(i,j),n(i,j)
				enddo
				write(1,*)
			enddo

		endif

	enddo



end program