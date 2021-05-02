!******************************************************************************************
!>
!  dense matrix utilities
!
!  Copyright (C) 1996 Roger Fletcher
!
!  Current version dated 26 May 2011
!
!  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
!  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
!  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

!    module utilities
    
!    contains
!******************************************************************************************

!******************************************************************************************
!>
!  solves Rx=b where R is nxn upper triangular. Solution overwrites b.
!  R is a single suffix array: the first nmax elements contain the first row
!  of R in positions 1:n, the next nmax-1 elements contain the second row of R,
!  and so on. nn indexes the element R(n,n) (where nn=n*(3-n)/2+(n-1)*nmax)

      subroutine rsol(n,nn,nmax,R,b)
      implicit double precision (a-h,o-z)
      dimension R(*),b(*)
      n1=nmax+1
      ii=nn
      b(n)=b(n)/R(nn)
      do i=n-1,1,-1
        ii=ii-n1+i
        b(i)=-scpr(-b(i),R(ii+1),b(i+1),n-i)/R(ii)
      end do
      end subroutine rsol
!******************************************************************************************

!******************************************************************************************
!>
!  solves Rt.x=b with same conventions as above
!  nn is not required on entry but is set on exit

      subroutine rtsol(n,nn,nmax,R,b)
      implicit double precision (a-h,o-z)
      dimension R(*),b(*)
      n2=nmax+2
      nn=1
      b(1)=b(1)/R(1)
      do i=2,n
        i1=i-1
        call mysaxpy(-b(i1),R(nn+1),b(i),n-i1)
        nn=nn+n2-i
        b(i)=b(i)/R(nn)
      end do
      end subroutine rtsol
!******************************************************************************************

!******************************************************************************************
!>
!  forms b=M.x where Q is nxn, stored by columns, with stride nmax

      subroutine Qprod(n,nmax,Q,x,b)
      implicit double precision (a-h,o-z)
      dimension Q(*),x(*),b(*)
      do i=1,n
        b(i)=0.D0
      end do
      i1=1
      do i=1,n
        call mysaxpy(x(i),Q(i1),b,n)
        i1=i1+nmax
      end do
      end subroutine Qprod
!******************************************************************************************

!******************************************************************************************
!>
!  forms b=M'.x where Q is nxn, stored by columns, with stride nmax

      subroutine Qtprod(n,nmax,Q,x,b)
      implicit double precision (a-h,o-z)
      dimension Q(*),x(*),b(*)
      i1=1
      do i=1,n
        b(i)=scpr(0.D0,Q(i1),x,n)
        i1=i1+nmax
      end do
      end subroutine Qtprod
!******************************************************************************************

!******************************************************************************************
!>
!
      subroutine brots(n,nmax,k,kk,R,v)
      implicit double precision (a-h,o-z)
      dimension R(*),v(*)
      ipip=kk
      do i=k-1,1,-1
        ip=i+1
        ipi=ipip-nmax+i
        ii=ipi-1
        call angle(v(i),v(ip),cos,sin)
        call rot(n-i,R(ipi),R(ipip),cos,sin)
        v(ip)=sin*R(ii)
        R(ii)=cos*R(ii)
        ipip=ii
      end do
      end subroutine brots
!******************************************************************************************

!******************************************************************************************
!>
! nr is either nc or nc+1

      subroutine frots(nr,nc,nmax,R,v)
      implicit double precision (a-h,o-z)
      dimension R(*),v(*)
      ii=1
      do i=1,nc
        ip=i+1
        ipi=ii+1
        ipip=ipi+nmax-i
        call angle(R(ii),v(ip),cos,sin)
        call rot(nr-i,R(ipi),R(ipip),cos,sin)
        ii=ipip
      end do
      end subroutine frots
!******************************************************************************************

!******************************************************************************************
!>
!
      subroutine angle(a,b,cos,sin)
      implicit double precision (a-h,o-z)
      z=sqrt(a**2+b**2)
      if (z==0.D0) then
        cos=1.D0
        sin=0.D0
      else
        cos=a/z
        sin=b/z
        a=z
        b=0.D0
      end if
      end subroutine angle
!******************************************************************************************

!******************************************************************************************
!>
!
      subroutine rot(n,a,b,cos,sin)
      implicit double precision (a-h,o-z)
      dimension a(*),b(*)
      if (sin==0.D0) then
        if (cos>0.D0) then
          do i=1,n
            b(i)=-b(i)
          end do
        else
          do i=1,n
            a(i)=-a(i)
          end do
        end if
      else if (cos==0.D0) then
        if (sin>=0.D0) then
          do i=1,n
            z=a(i)
            a(i)=b(i)
            b(i)=z
          end do
        else
          do i=1,n
            z=a(i)
            a(i)=-b(i)
            b(i)=-z
          end do
        end if
      else
        do i=1,n
          z=a(i)
          a(i)=cos*z+sin*b(i)
          b(i)=sin*z-cos*b(i)
        end do
      end if
      end subroutine rot
!******************************************************************************************

!******************************************************************************************
!>
!
      subroutine mysaxpy(a,x,y,n)
      implicit double precision (a-h,o-z)
      dimension x(*),y(*)
      if (a==0.D0) return
      do i=1,n
        y(i)=y(i)+a*x(i)
      end do
      end subroutine mysaxpy
!******************************************************************************************

!******************************************************************************************
!>
!  saxpy with stride

      subroutine saxpys(a,x,is,y,n)
      implicit double precision (a-h,o-z)
      dimension x(*),y(*)
      if (a==0.D0) return
      ix=1
      do i=1,n
        y(i)=y(i)+a*x(ix)
        ix=ix+is
      end do
      end subroutine saxpys
!******************************************************************************************

!******************************************************************************************
!>
!  saxpy with result in x

      subroutine saxpyx(a,x,y,n)
      implicit double precision (a-h,o-z)
      dimension x(*),y(*)
      if (a==0.D0) then
        do i=1,n
          x(i)=y(i)
        end do
      else
        do i=1,n
          x(i)=y(i)+a*x(i)
        end do
      end if
      end subroutine saxpyx
!******************************************************************************************

!******************************************************************************************
!>
!  saxpy with result in z

      subroutine saxpyz(a,x,y,z,n)
      implicit double precision (a-h,o-z)
      dimension x(*),y(*),z(*)
      if (a==0.D0) then
        do i=1,n
          z(i)=y(i)
        end do
      else
        do i=1,n
          z(i)=y(i)+a*x(i)
        end do
      end if
      end subroutine saxpyz
!******************************************************************************************

!******************************************************************************************
!>
!  saxpy with interchange of x and y

      subroutine saxpyi(a,x,y,n)
      implicit double precision (a-h,o-z)
      dimension x(*),y(*)
      if (a==0.D0) then
        do i=1,n
          call rexch(x(i),y(i))
        end do
      else
        do i=1,n
          z=y(i)
          y(i)=x(i)+a*y(i)
          x(i)=z
        end do
      end if
      end subroutine saxpyi
!******************************************************************************************

!******************************************************************************************
!>
!
      function scpr(a,x,y,n)
      implicit double precision (a-h,o-z)
      dimension x(*),y(*)
      scpr=a
      do i=1,n
        scpr=scpr+x(i)*y(i)
      end do
      end function scpr
!******************************************************************************************


!     function xlen(a,x,n)
!     implicit double precision (a-h,o-z)
!     dimension x(*)
!  finds the l_2 length of [a:x] where a is either 0.D0 or 1.D0
!  if overflow occurs the function is calculated in a less efficient way.
!  Users who cannot trap overflow should either use this method of calculation,
!  or use the alternative routine "xlen" below which is not quite so well
!  protected against overflow.
!     external  ieee_handler, abort
!     integer   ieee_flags, ieeer, ieee_handler
!     external  ieee_flags
!     character out*16
!     out = ''
!     ieeer = ieee_flags ( 'clearall','all','',out )
!     ieeer=ieee_handler('clear','overflow',abort)
!  this call of ieee_handler assumes that
!         ieeer=ieee_handler('set','overflow',abort)
!  has been set in the driver. If not this call of ieee_handler and that below
!  should be removed
!     xlen=a
!     do i=1,n
!       xlen=xlen+x(i)**2
!     end do
!     xlen=sqrt(xlen)
!     ieeer=ieee_flags ( 'get','exception','',out )
!     if (out=='overflow') then
!       call linf(n,x,xmx,i)
!       xmx=max(xmx,1.D0) %this is needed if normalization is always used
!       xlen=(a/xmx)**2
!       do i=1,n
!         xlen=xlen+(x(i)/xmx)**2
!       end do
!       xlen=xmx*sqrt(xlen)
!       ieeer=ieee_flags ( 'clear','overflow','',out )
!     end if
!     ieeer=ieee_handler('set','overflow',abort)
!     return
!     end

!******************************************************************************************
!>
!
      function xlen(a,x,n)
      implicit double precision (a-h,o-z)
      dimension x(*)
      xlen=a
      do i=1,n
        xlen=xlen+x(i)**2
      end do
      xlen=sqrt(xlen)
      end function xlen
!******************************************************************************************

!******************************************************************************************
!>
!
      subroutine linf(n,x,z,iz)
      implicit double precision (a-h,o-z)
      dimension x(*)
      z=0.D0
      do i=1,n
        a=abs(x(i))
        if (a>z) then
          z=a
          iz=i
        end if
      end do
      end subroutine linf
!******************************************************************************************

!******************************************************************************************
!>
!
      subroutine r_shift(r,n,k)
      implicit double precision (a-h,o-z)
      dimension r(*)
      if (k>0) then
        do i=1,n
          r(i)=r(i+k)
        end do
      else if (k<0) then
        do i=n,1,-1
          r(i)=r(i+k)
        end do
      end if
      end subroutine r_shift
!******************************************************************************************

!******************************************************************************************
!>
!
      subroutine ishift(l,n,k)
      implicit double precision (a-h,o-z)
      dimension l(*)
      if (k>0) then
        do i=1,n
          l(i)=l(i+k)
        end do
      else if (k<0) then
        do i=n,1,-1
          l(i)=l(i+k)
        end do
      end if
      end subroutine ishift
!******************************************************************************************

!******************************************************************************************
!>
!
      subroutine rexch(a,b)
      double precision a,b,z
      z=a
      a=b
      b=z
      end subroutine rexch
!******************************************************************************************

!******************************************************************************************
!>
!
      subroutine vexch(a,b,n)
      double precision a,b,z
      dimension a(*),b(*)
      do i=1,n
        z=a(i)
        a(i)=b(i)
        b(i)=z
      end do
      end subroutine vexch
!******************************************************************************************

!******************************************************************************************
!>
!
      subroutine iexch(i,j)
      k=i
      i=j
      j=k
      end subroutine iexch
!******************************************************************************************

!******************************************************************************************
!    end module utilities
!******************************************************************************************