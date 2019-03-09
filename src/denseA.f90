!Christen this file    denseA.f

!  Copyright (C) 1996 Roger Fletcher

!  Current version dated 21 May 1998

!  ******************************************
!  Specification of A in dense matrix format
!  ******************************************

!  The matrix A contains gradients of the linear terms in the objective
!  function (column 0) and the general constraints (columns 1:m).
!  No explicit reference to simple bound constraints is required in A.
!  The information is set in the parameters a(*) and la.

!  In this dense case A is set in standard matrix format as a(la,0:m), where la
!  is the stride between columns. la is an integer which must be greater or
!  equal to n.

!  In the straightforward case that la=n, columns of A follow successively
!  in the space occupied by a(.).


      subroutine saipy(s,a,la,i,y,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),y(*)
!  saxpy with column i of A
      call mysaxpy(s,a(1,i),y,n)
      return
      end

      subroutine isaipy(s,a,la,i,y,n,lr,li)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),y(*),lr(*),li(*)
!  indirectly addressed saxpy with column i of A
      call isaxpy(s,a(1,i),lr,y,n)
      return
      end

      subroutine isaipy1(s,a,la,i,y,n,lr,li,m1)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),y(*),lr(*),li(*)
!  indirectly addressed saxpy with column i of A_1
      call isaxpy(s,a(1,i),lr,y,m1)
      return
      end

!     subroutine isaipy2(s,a,la,i,y,n,lr,li,m1)
!     implicit double precision (a-h,r-z), integer (i-q)
!     dimension a(la,0:*),y(*),lr(*),li(*)
!  indirectly addressed saxpy with column i of A_2
!     call isaxpy(s,a(1,i),lr(m1+1),y(m1+1),n-m1)
!     return
!     end

!     subroutine ssaipy(s,a,la,i,y,n)
!     implicit double precision (a-h,r-z), integer (i-q)
!     dimension a(la,0:*),y(*)
!  ssaxpy with column i of A
!     call ssaxpy(s,a(1,i),y,n)
!     return
!     end

!     subroutine ssaxpy(a,x,y,n)
!     implicit double precision (a-h,r-z), integer (i-q)
!     dimension x(*),y(*)
!  saxpy with squares of x
!     do i=1,n
!       y(i)=y(i)+a*x(i)**2
!     end do
!     return
!     end

      function aiscpr(n,a,la,i,x,b)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*)
!  scalar product with column i of A
      aiscpr=scpr(b,a(1,i),x,n)
      return
      end

      function daiscpr(n,a,la,i,x,b)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*)
      DOUBLE PRECISION daiscpr,dscpr
      daiscpr=dscpr(b,a(1,i),x,n)
      return
      end

      function aiscpri(n,a,la,i,x,b,lr,li)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*),lr(*),li(*)
!  indirectly addressed scalar product with column i of A
      aiscpri=scpri(b,a(1,i),lr,x,n)
      return
      end

      function daiscpri(n,a,la,i,x,b,lr,li)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*),lr(*),li(*)
      DOUBLE PRECISION daiscpri,dscpri
      daiscpri=dscpri(b,a(1,i),lr,x,n)
      return
      end

      function aiscpri1(n,a,la,i,x,b,lr,li,m1)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*),lr(*),li(*)
!  indirectly addressed scalar product with column i of A_1
      aiscpri1=scpri(b,a(1,i),lr,x,m1)
      return
      end

      function aiscpri2(n,a,la,i,x,b,lr,li,m1)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*),lr(*),li(*)
!  indirectly addressed scalar product with column i of A_2
      aiscpri2=scpri(b,a(1,i),lr(m1+1),x(m1+1),n-m1)
      return
      end

      function ailen(n,a,la,i)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*)
!  L2 length of column i of A
      ailen=scpr(0.D0,a(1,i),a(1,i),n)
      ailen=sqrt(ailen)
      return
      end

      subroutine iscatter(a,la,i,li,an,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),li(*),an(*)
!  indirect scatter into vector an
      do j=1,n
        an(li(j))=a(j,i)
      end do
      return
      end

      subroutine iunscatter(a,la,i,li,an,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),li(*),an(*)
!  included for compatibility with sparseA.f
      return
      end

      function aij(i,j,a,la)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*)
!  get element A(i,j)
      aij=a(i,j)
      return
      end

      subroutine setaij(aij,i,j,a,la)
      implicit double precision (a-h,o-z)
      dimension a(la,0:*)
!  set element A(i,j)
      aij=a(i,j)
      end

      subroutine isaxpy(a,x,lr,y,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension x(*),lr(*),y(*)
!  saxpy with x indirectly addressed
      if (a==0.D0) return
      do i=1,n
        y(i)=y(i)+a*x(lr(i))
      end do
      return
      end

      function dscpr(a,x,y,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension x(*),y(*)
      DOUBLE PRECISION dscpr
      dscpr=dble(a)
      do i=1,n
        dscpr=dscpr+dble(x(i))*dble(y(i))
      end do
      return
      end

      function scpri(a,x,lr,y,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension x(*),lr(*),y(*)
!  scpr with x indirectly addressed
      scpri=a
      do i=1,n
        scpri=scpri+x(lr(i))*y(i)
      end do
      return
      end

      function dscpri(a,x,lr,y,n)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension x(*),lr(*),y(*)
      DOUBLE PRECISION dscpri
      dscpri=dble(a)
      do i=1,n
        dscpri=dscpri+dble(x(lr(i)))*dble(y(i))
      end do
      return
      end

      subroutine cscale(n,m,a,la,x,bl,bu,s,menu,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(la,0:*),x(*),bl(*),bu(*),s(*)

!     Constraint scaling procedure for use prior to calling bqpd when using
!     denseA.f

!     The user must set the parameter menu to control how the
!     x-variables are scaled (or equivalently how constraints i = 1:n
!     are scaled), as follows

!     menu = 1 indicates that a unit scaling applies to the x-variables

!     menu = 2 the user provides estimates s(i)>0 of the magnitude of
!              x(i) for i = 1:n. In this case the elements  x(i), bl(i), bu(i)
!              are divided by s(i) for i = 1:n.

!     In all cases, cscale goes on to scale the general constraints, in
!     such a way that the normal vector of each nontrivial constraint in
!     the scaled problem has an l_2 norm of unity. This scaling is also
!     applied to the right hand sides  bl(i), bu(i) for i = n+1:n+m.
!     The scaled data overwrites the original data.

!     cscale also scales the constant vector of the quadratic function,
!     which is found in a(1:n). However if a non-unit x-variable scaling
!     is used, it is necessary for the user to scale the Hessian matrix
!     G appropriately. This can be done by passing the x-variable scale
!     factors s(i) i = 1:n into the subroutine gdotx using the
!     parameter ws, and multiplying G(i,j) by s(i)*s(j) (possibly
!     implicitly).

!     cscale sets ifail = 1 to indicate that some s(i)< = 0,
!             and ifail = 2 to indicate an incorrect setting of menu.
!       Otherwise ifail = 0.

      ifail=2
      if (menu<1 .or. menu>2) return
!     z=1.D0/log(2.D0)
      if (menu==1) then
        do j=1,n
          s(j)=1.D0
        end do
      else
        ifail=1
        do j=1,n
          if (s(j)<=0.D0) return
        end do
!       if (menu==2) then
!         do j=1,n
!           s(j)=2.D0**nint(log(s(j))*z)
!         end do
!       end if
        do j=1,n
          if (s(j)/=1.D0) then
            x(j)=x(j)/s(j)
            bl(j)=bl(j)/s(j)
            bu(j)=bu(j)/s(j)
            a(j,0)=a(j,0)*s(j)
          end if
        end do
      end if
      do i=1,m
        t=0.D0
        do j=1,n
          a(j,i)=a(j,i)*s(j)
          t=t+a(j,i)**2
        end do
        t=sqrt(t)
        if (t==0.D0) then
          s(n+i)=1.D0
        else
!         t=2.D0**nint(log(t)*z)
          s(n+i)=t
          do j=1,n
            a(j,i)=a(j,i)/t
          end do
          bl(n+i)=bl(n+i)/t
          bu(n+i)=bu(n+i)/t
        end if
      end do
      ifail=0
      return
      end
