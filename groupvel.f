C	BOX SCHEME NUMERICAL GROUP VELOCITY
c	In this example we computed the BOX SCHEME numerical solution of
c	U(t)+b*U(x) = 0; subject to the initial condition
c	v(x)=exp(-1.0d0*epslon*(x(i)-1.0d0)**2)*sin(alpha*pi*x(i))
	program groupvel
	implicit none
	integer i,j,n,zx,zt,zx1
	double precision b, dx, dt, c, u, x, t, pi,epslon,alpha,xmax
	parameter (zx=320+1, zt=80+1)
	dimension x(1:zx),t(1:zt),u(1:zx, 1:zt)
	open(unit=13, file="upwind_groupvel.dat")
	open(unit=15, file="groupvel.dat")

 20	format(1x,3(f25.20,5x))
	
       b = 1.0d0   
	epslon = 10.0d0
	alpha = 16.0d0


       dx = 1.25d-2
       dt = 2.5d-2
	xmax = 4.0d0
	zx1 = xmax/dx +1
	c = b * (dt/dx)
	
	pi = 4.0d0 * atan(1.0d0)
	x(1) = 0.0d0
	t(1) = 0.0d0

	do i=2,zx1+1
	   x(i)=x(i-1)+dx
	end do

	do i=2,zt
	   t(i)=t(i-1)+dt
	end do

c	INITIAL CONDITION

	do i =1,zt
	u(1,i)=0.0d0
	end do

	do i = 1,zx1 
	   u(i,1) = exp(-1.0d0*epslon*(x(i)-1.0d0)**2)*sin(alpha*pi*x(i))
	   write(13,20) x(i),t(1),u(i,1)
	end do

	u(1,2) = 0.0d0

c	BOX SCHEME
	do n=1,zt-1
	   do j=1,zx-2
	   u(j+1,n+1) = u(j,n)+(((1-c)/(1+c))*(u(j+1,n)-u(j,n+1)))
	   end do
	end do

       do n = 1, zt
	   do j = 1, zx1
		write(15,20) x(j),t(n),u(j,n)
	   end do
	end do

	do n = 1, zt-1
	   do j = 2, zx
		u(j,n+1) = u(j,n) - c * ( u(j,n) - u(j-1,n) )
	   enddo
	enddo

       do n = 1, zt
	   do j = 1, zx
		write(13,20) x(j),t(n),u(j,n)
	   end do
	end do

	stop
	end
