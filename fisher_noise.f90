!Code to solve the standard Fisher equation in two dimensions, with noise
!Equation is dpsi/dt = alpha psi (1-psi) + D grad^2 psi + sqrt(psi)*noise
!Uses Dickmann discretization scheme

program fisher_noise

    implicit none
	integer, parameter :: SIZE=250
	real(8) :: psi(SIZE,SIZE), dpsi(SIZE,SIZE) !psi is a continuous field keeping track of psi
	integer :: n(SIZE,SIZE) !n is the discretized version of psi
	real(8) :: D,dt,dx,alpha
	real(8) :: psi_min, noise_strength, sqrtdt
	integer :: i,j,t,maxt,iup,idwn
    integer :: totalshift

	!parameters
	D=1
	alpha=1
	dt=0.001
	dx=0.1
	maxt=50000
	psi_min=0.0001
	noise_strength=10;
	sqrtdt=sqrt(dt)

	!initialize random number generator
	call srand(123)

	!initialize
	do j = 1,SIZE
		do i = 1,SIZE
			n(i,j)=(0.5*tanh( (SIZE/4-j)/10.0) + 0.5)/psi_min
			psi(i,j)=0
		enddo
	enddo

    totalshift=0

	!timestep
	do t = 1, maxt

		write(6,*) t, totalshift

		do j = 2,SIZE-1
			do i=1, SIZE

				iup = merge(1,i+1,i==SIZE)
				idwn = merge(SIZE,i-1,i==1)

				!calculate change in psi over the timestep
				dpsi(i,j) = dt*alpha*n(i,j)*(1-psi_min*n(i,j)) &
				 + (D*dt)/(dx*dx)*(n(iup,j) + n(idwn,j) + n(i,j+1) + n(i,j-1) - 4*n(i,j)) &
				 + noise_strength*sqrt(real(n(i,j)))*sqrtdt*(rand()-0.5)

			enddo
		enddo

		!simulataneous update of psi array
		do j = 2,SIZE-1
			do i=1,SIZE
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
            call print_grid(t)
		endif

        call shift_everything()

	enddo

    contains
    
        subroutine shift_everything() !shift all the fields so that the contour is roughly centred in the sim box

		integer :: buffer,ymin,x,shift
		real(8) :: tolerance
        integer, dimension(SIZE,SIZE) :: temp

		buffer = 200 !distance of contour from bottom of box
		tolerance = 0.001

		ymin=SIZE


		do j = 1,SIZE
			do i=1,SIZE
				if(n(i,j)+psi(i,j)<tolerance/psi_min .and. j<ymin) then 
					ymin=j
				endif
			enddo
		enddo

		shift = ymin-buffer

		totalshift = totalshift + shift

		if (shift>0) then
			do j = 1,SIZE
				do i=1,SIZE

					x=j+shift
					if(x>SIZE) then 
						temp(i,j)=0
                        dpsi(i,j)=0 !dpsi being used as a temporary array to store psi values
					else
						temp(i,j)=n(i,x)
                        dpsi(i,j)=psi(i,x)
					endif 
				enddo
			enddo

		else if(shift<0) then
			
			do j = 1,SIZE
				do i=1,SIZE

					x=j+shift
					if(x<=0) then 
						temp(i,j)=1/psi_min
                        dpsi(i,j)=0
					else
						temp(i,j)=n(i,x)
                        dpsi(i,j)=psi(i,x)
					endif
				enddo 
			enddo

        else
            return

		endif


		do j = 1,SIZE
			do i=1,SIZE
				psi(i,j)=dpsi(i,j)
                n(i,j)=temp(i,j)
			enddo
		enddo

		end subroutine shift_everything

		subroutine print_grid(label)

			integer, intent(in) :: label
			character(len=1024) :: filename

			write(filename,"(A3,I0,A4)") 'out',label, '.dat'
			open(1,file=filename)

			do j = 1,SIZE
				do i=1,SIZE
					write(1,*) i,j,psi(i,j),n(i,j)
				enddo
				write(1,*)
			enddo

		end subroutine print_grid



end program
