!Code to solve a 2-component standard Fisher equation in two dimensions
!Equations are dphi_i/dt = alpha_i phi_i (1-sum_i phi_i) + D grad^2 phi_i
!i=1,2

program fisher_standard

    implicit none
	integer, parameter :: SIZE=500
	real(8) :: phi_1(SIZE,SIZE), dphi_1(SIZE,SIZE), phi_2(SIZE,SIZE), dphi_2(SIZE,SIZE)
	real(8) :: D,dt,dx,alpha_1,alpha_2
	integer :: i,j,t,maxt,iup,idwn
	integer  :: totalshift

	!parameters
	D=.1
	alpha_1=1
    alpha_2=1.5
	dt=0.001
	dx=0.1
	maxt=50000
	totalshift=0

	!initialize
	do j = 1,SIZE
		do i = 1,SIZE
			phi_2(i,j)=( (0.5*tanh( (SIZE/4-j)/10.0) + 0.5) * exp (-(i-SIZE/2)*(i-SIZE/2)/100.0) )
			phi_1(i,j)=0.5*tanh( (SIZE/4-j)/10.0) + 0.5 - phi_2(i,j)
		enddo
	enddo

	!timestep
	do t = 1, maxt

		write(6,*) t,totalshift

		do j = 2,SIZE-1
			do i=1, SIZE

				iup = merge(1,i+1,i==SIZE)
				idwn = merge(SIZE,i-1,i==1)

				!calculate change in phi_i over the timestep
				dphi_1(i,j) = dt*alpha_1*phi_1(i,j)*(1-phi_1(i,j)-phi_2(i,j)) + (D*dt)/(dx*dx)*(phi_1(iup,j) + phi_1(idwn,j) &
                 + phi_1(i,j+1) + phi_1(i,j-1) - 4*phi_1(i,j))
				dphi_2(i,j) = dt*alpha_2*phi_2(i,j)*(1-phi_1(i,j)-phi_2(i,j)) + (D*dt)/(dx*dx)*(phi_2(iup,j) + phi_2(idwn,j) + &
                 phi_2(i,j+1) + phi_2(i,j-1) - 4*phi_2(i,j))

			enddo
		enddo

		!simulataneous update of arrays
		do j = 2,SIZE-1
			do i=1,SIZE
				phi_1(i,j) = phi_1(i,j) + dphi_1(i,j)
				phi_2(i,j) = phi_2(i,j) + dphi_2(i,j)
			enddo
		enddo

		!zero flux BC in the j-direction
		do i = 1,SIZE
			phi_1(i,1)=phi_1(i,2)
			phi_1(i,SIZE-1)=phi_1(i,SIZE)
			phi_2(i,1)=phi_2(i,2)
			phi_2(i,SIZE-1)=phi_2(i,SIZE)
		enddo

		!print data every 1000 time steps
		if(mod(t,1000)==0 .or. t==1) then
			call print_grid(t)
		endif

        call shift_everything()

	enddo

contains

	subroutine shift_everything() !shift all the fields so that the contour is roughly centred in the sim box

		integer :: buffer,ymin,x,shift
		real(8) :: tolerance

		buffer = 150 !distance of contour from bottom of box
		tolerance = 0.001

		ymin=SIZE


		do j = 1,SIZE
			do i=1,SIZE
				if(phi_1(i,j)<=tolerance .and. phi_2(i,j)<=tolerance .and. j<ymin) then 
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
						dphi_1(i,j)=0
                        dphi_2(i,j)=0
					else
						dphi_1(i,j)=phi_1(i,x) !dphi's used here as temp arrays to store new values of phi 
						dphi_2(i,j)=phi_2(i,x)
					endif 
				enddo
			enddo

		else if(shift<0) then
			
			do j = 1,SIZE
				do i=1,SIZE

					x=j+shift
					if(x<=0) then 
						dphi_1(i,j)=1
                        dphi_2(i,j)=1
					else
						dphi_1(i,j)=phi_1(i,x)
                        dphi_2(i,j)=phi_2(i,x)
					endif
				enddo 
			enddo

        else
            return

		endif


		do j = 1,SIZE
			do i=1,SIZE
				phi_1(i,j)=dphi_1(i,j)
				phi_2(i,j)=dphi_2(i,j)
			enddo
		enddo

		end subroutine shift_everything

		subroutine print_grid(label)

			integer, intent(in) :: label
			character(len=1024) :: filename

			write(filename,"(A7,I0,A4)") 'std_out',label, '.dat'
			open(1,file=filename)

			do j = 1,SIZE
				do i=1,SIZE
					write(1,*) i,j,phi_1(i,j),phi_2(i,j)
				enddo
				write(1,*)
			enddo

		end subroutine print_grid

end program







