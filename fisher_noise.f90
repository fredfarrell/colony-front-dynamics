!Code to solve the standard Fisher equation in two dimensions, with noise
!Equations are is dphi_i/dt = alpha phi_i (1-phi_1-phi_2) + D grad^2 phi_i + sqrt(phi_i)*noise
!Uses Dickmann discretization scheme

program fisher_noise

    implicit none
	integer, parameter :: SIZE=500
	real(8) :: psi_1(SIZE,SIZE), dpsi_1(SIZE,SIZE), psi_2(SIZE,SIZE), dpsi_2(SIZE,SIZE) !psi_i are continuous fields keeping track of phi_i
	integer :: n(SIZE,SIZE), m(SIZE,SIZE) !n,m discretized versions of psi_1, psi_2 resp.
	real(8) :: D,dt,dx,alpha_1,alpha_2
	integer :: i,j,t,maxt,iup,idwn
    integer :: totalshift
    real(8) :: noise, phi_min, g, noise_max, sqrtdt

	!parameters
	D=.1
	alpha_1=1.0
	alpha_2=1.5
	dt=0.0001
	dx=0.1
	maxt=100000
	g=1
	sqrtdt=sqrt(dt)
	phi_min=0.0001 !size of discretizaton of the fields (i.e. dphi)

	!initialize random number generator
	call srand(123)

	!initialize
	do j = 1,SIZE
		do i = 1,SIZE
			m(i,j)=( (0.5*tanh( (SIZE/4-j)/10.0) + 0.5) * exp (-(i-SIZE/2)*(i-SIZE/2)/100.0) )/phi_min
			n(i,j)=(0.5*tanh( (SIZE/4-j)/10.0) + 0.5)/phi_min - m(i,j)
		enddo
	enddo

	call print_grid(0)

    totalshift=0

	!timestep
	do t = 1, maxt

		write(6,*) t, totalshift

		do j = 2,SIZE-1
			do i=1, SIZE

				iup = merge(1,i+1,i==SIZE)
				idwn = merge(SIZE,i-1,i==1)

				!calculate change in psi_i over the timestep
				!noise term with opposite sign for the two fields so that their total is unaffected by it
				noise  = rand()

				dpsi_1(i,j) = dt*alpha_1*n(i,j)*(1-phi_min*n(i,j)-phi_min*m(i,j)) &
				 + (D*dt)/(dx*dx)*(n(iup,j) + n(idwn,j) + n(i,j+1) + n(i,j-1) - 4*n(i,j)) &
				 + g*sqrt( real(n(i,j)*m(i,j)) )*sqrtdt*(noise-0.5)

				dpsi_2(i,j) = dt*alpha_2*m(i,j)*(1-phi_min*n(i,j)-phi_min*m(i,j)) &
				 + (D*dt)/(dx*dx)*(m(iup,j) + m(idwn,j) + m(i,j+1) + m(i,j-1) - 4*m(i,j)) &
				 - g*sqrt( real(n(i,j)*m(i,j)) )*sqrtdt*(noise-0.5)

			enddo
		enddo

		!simulataneous update of psi arrays
		do j = 2,SIZE-1
			do i=1,SIZE
				psi_1(i,j) = psi_1(i,j) + dpsi_1(i,j)
				psi_2(i,j) = psi_2(i,j) + dpsi_2(i,j)

				! add integer part of psi_1 to n
				n(i,j) = n(i,j) + int(psi_1(i,j))
				psi_1(i,j) = psi_1(i,j) - int(psi_1(i,j))
				m(i,j) = m(i,j) + int(psi_2(i,j))
				psi_2(i,j) = psi_2(i,j) - int(psi_2(i,j))

				if(n(i,j)<0 .or. m(i,j)<0) then
					write(6,*) "Error: negative"
					call exit(0)
				endif

			enddo
		enddo

		!zero flux BC in the j-direction
		do i = 1,SIZE
			psi_1(i,1)=psi_1(i,2)
			psi_1(i,SIZE-1)=psi_1(i,SIZE)
			psi_2(i,1)=psi_2(i,2)
			psi_2(i,SIZE-1)=psi_2(i,SIZE)
		enddo

		!print data every 1000 time steps
		if(mod(t,10000)==0) then
            call print_grid(t)
		endif

        !call shift_everything()

	enddo

contains
 
!  	subroutine shift_everything() !shift all the fields so that the contour is roughly centred in the sim box

! 		integer :: buffer,ymin,x,shift
! 		real(8) :: tolerance
!         integer, dimension(SIZE,SIZE) :: temp

! 		buffer = 200 !distance of contour from bottom of box
! 		tolerance = 0.001

! 		ymin=SIZE


! 		do j = 1,SIZE
! 			do i=1,SIZE
! 				if(n(i,j)+psi_1(i,j)<tolerance/phi_min .and. j<ymin) then 
! 					ymin=j
! 				endif
! 			enddo
! 		enddo

! 		shift = ymin-buffer

! 		totalshift = totalshift + shift

! 		if (shift>0) then
! 			do j = 1,SIZE
! 				do i=1,SIZE

! 					x=j+shift
! 					if(x>SIZE) then 
! 						temp(i,j)=0
!                         dpsi_1(i,j)=0 !dpsi_1 being used as a temporary array to store psi_1 values
! 					else
! 						temp(i,j)=n(i,x)
!                         dpsi_1(i,j)=psi_1(i,x)
! 					endif 
! 				enddo
! 			enddo

! 		else if(shift<0) then
			
! 			do j = 1,SIZE
! 				do i=1,SIZE

! 					x=j+shift
! 					if(x<=0) then 
! 						temp(i,j)=1/phi_min
!                         dpsi_1(i,j)=0
! 					else
! 						temp(i,j)=n(i,x)
!                         dpsi_1(i,j)=psi_1(i,x)
! 					endif
! 				enddo 
! 			enddo

!         else
!             return

! 		endif


! 		do j = 1,SIZE
! 			do i=1,SIZE
! 				psi_1(i,j)=dpsi_1(i,j)
!                 n(i,j)=temp(i,j)
! 			enddo
! 		enddo

! 	end subroutine shift_everything

	subroutine print_grid(label)

		integer, intent(in) :: label
		character(len=1024) :: filename

		write(filename,"(A3,I0,A4)") 'out',label, '.dat'
		open(1,file=filename)

		do j = 1,SIZE
			do i=1,SIZE
				write(1,*) i,j,psi_1(i,j),psi_2(i,j),n(i,j),m(i,j)
			enddo
			write(1,*)
		enddo

	end subroutine print_grid



end program
