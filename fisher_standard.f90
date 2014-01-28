!Code to solve the standard Fisher equation in two dimensions
!Equation is d\phi/dt = alpha \phi (1-\phi) + D grad^2 \phi

program fisher_standard

	integer, parameter :: SIZE=250
	real(8) :: phi(SIZE,SIZE), dphi(SIZE,SIZE), temp(SIZE,SIZE)
	real(8) :: D,dt,dx,alpha
	integer :: i,j,t,maxt,iup,idwn
	real(8)  :: totalshift

	!parameters
	D=1
	alpha=1
	dt=0.001
	dx=0.1
	maxt=5000
	totalshift=0

	!initialize
	do i = 1,SIZE
		do j = 1,SIZE
			phi(i,j)=0.5*tanh( (SIZE/4-j)/10.0) + 0.5
		enddo
	enddo

	!timestep
	do t = 1, maxt

		if(t==1) then
			call print_grid(1)
		endif

		if(t==1) then
			call print_grid(2)
		endif

		write(6,*) t,totalshift

		do i = 1,SIZE
			do j=2, SIZE-1

				iup = merge(1,i+1,i==SIZE)
				idwn = merge(SIZE,i-1,i==1)

				!calculate change in phi over the timestep
				dphi(i,j) = dt*alpha*phi(i,j)*(1-phi(i,j)) + (D*dt)/(dx*dx)*(phi(iup,j) + phi(idwn,j) + phi(i,j+1) + phi(i,j-1) - 4*phi(i,j))

			enddo
		enddo

		!simulataneous update of phi array
		do i = 1,SIZE
			do j=2,SIZE-1
				phi(i,j) = phi(i,j) + dphi(i,j)
			enddo
		enddo

		!zero flux BC in the j-direction
		do i = 1,SIZE
			phi(i,1)=phi(i,2)
			phi(i,SIZE-1)=phi(i,SIZE)
		enddo

		!print data every 1000 time steps
		if(mod(t,1000)==0) then
			call print_grid(t)
		endif

		call shift_everything()

	enddo

contains

	subroutine shift_everything() !shift all the fields so that the contour is roughly centred in the sim box

		integer :: buffer,ymin
		real(8) :: tolerance

		buffer = 100 !distance of contour from bottom of box
		tolerance = 0.001

		ymin=SIZE

		do j = 1,SIZE
			do i=1,SIZE
				if(phi(i,j)<=tolerance .and. j<ymin) then 
					ymin=j
				endif
			enddo
		enddo

		shift = ymin-buffer
		totalshift = totalshift + shift

		if (shift>0) then
			do i = 1,SIZE
				do j=1,SIZE

					x=j+shift
					if(x>SIZE) then 
						temp(i,j)=0
					else
						temp(i,j)=phi(i,x)
					endif 
				enddo
			enddo

		else if(shift<0) then
			
			do i = 1,SIZE
				do j=1,SIZE

					x=j+shift
					if(x<=0) then 
						temp(i,j)=1
					else
						temp(i,j)=phi(i,x)
					endif
				enddo 
			enddo

		endif

		do i = 1,SIZE
			do j=1,SIZE
				phi(i,j)=temp(i,j)
			enddo
		enddo

		end subroutine shift_everything

		subroutine print_grid(label)

			integer, intent(in) :: label
			character(len=1024) :: filename

			write(filename,"(A7,I0,A4)") 'std_out',label, '.dat'
			open(1,file=filename)

			do i = 1,SIZE
				do j=1,SIZE
					write(1,*) i,j,phi(i,j)
				enddo
				write(1,*)
			enddo

		end subroutine print_grid

end program







