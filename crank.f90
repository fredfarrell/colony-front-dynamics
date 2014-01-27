	program implicit
	implicit none
	integer, parameter :: ISIZE = 5000
	double precision dx,dt,c0,Diff,R
	double precision c_old(ISIZE),c_new(ISIZE)
	double precision a(ISIZE),b(ISIZE),c(ISIZE),d(ISIZE),x(ISIZE)
	integer t,i,j,tfinal

	dx=0.4
	dt=0.05

    Diff=100.0

    tfinal=1000

    R=Diff*dt/(dx*dx)

!	initial conditions
	do i=1,ISIZE
	do j=1,ISIZE
	
           c_new(i) = exp(-(i-ISIZE/2)*(i-ISIZE/2)/100.d0)

	enddo
	enddo

    !write initial condition to a file
	do i=1,ISIZE

        write(10,*) i,c_new(i)

    enddo

	do t=1,tfinal

        c_old=c_new

        !a,b and c in the traidiagonal matrix are constants
        a=-R
        b=1+2*R
        c=-R

        !fill in the d vector
        do i=2, ISIZE-1

            d(i)=R*c_old(i-1)+(1-2*R)*c_old(i)+R*c_old(i+1)

        enddo

        !zero flux BCs
        a(1)=0
        b(1)=1
        c(1)=-1
        d(1)=0

        a(ISIZE)=-1
        b(ISIZE)=1
        c(ISIZE)=0
        d(ISIZE)=0
        

        !invert the matrix
        call solve_tridiag(a,b,c,d,x,ISIZE)

        do i=1,ISIZE

            c_new(i)=x(i)

        enddo


    enddo
            
	do i=1,ISIZE

        write(6,*) i,c_new(i)

    enddo

	return
	end

      subroutine solve_tridiag(a,b,c,v,x,n) !from wikipedia
      implicit none
!      a - sub-diagonal (means it is the diagonal below the main diagonal)
!      b - the main diagonal
!      c - sup-diagonal (means it is the diagonal above the main diagonal)
!      v - right part
!      x - the answer
!      n - number of equations
 
        integer,intent(in) :: n
        real(8),dimension(n),intent(in) :: a,b,c,v
        real(8),dimension(n),intent(out) :: x
        real(8),dimension(n) :: bp,vp
        real(8) :: m
        integer i
 
! Make copies of the b and v variables so that they are unaltered by this subroutine
        bp(1) = b(1)
        vp(1) = v(1)
 
        !The first pass (setting coefficients):
firstpass: do i = 2,n
         m = a(i)/bp(i-1)
         bp(i) = b(i) - m*c(i-1)
         vp(i) = v(i) - m*vp(i-1)
        end do firstpass
 
         x(n) = vp(n)/bp(n)
        !The second pass (back-substition)
backsub:do i = n-1, 1, -1
          x(i) = (vp(i) - c(i)*x(i+1))/bp(i)
        end do backsub
 	

    	end subroutine solve_tridiag
