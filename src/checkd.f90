!Christen this file checkd.f

!  Copyright (C) 2010 Roger Fletcher

!  Current version 20 January 2011

!  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
!  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
!  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

      subroutine checkd(n,m,x,al,ws,lws,maxa,maxla,maxu,maxiu, &
        mxws,mxlws,tol)
      implicit double precision (a-h,o-z)

!  first derivative checking subroutine for use with filterSD code

!  Parameters
!  **********
!  n    set number of variables
!  m    set number of constraints
!  x    set any vector of variables in x(i), i=1,n
!  al   set a vector of differencing increments in al(i), i=1,n
!        (al(i) is the difference interval for x(i) and must be nonzero,
!        say 0.01 times the typical magnitude of x(i))
!  ws   double precision workspace (the amount required for filterSD.f should
!       be plenty)
!  lws  integer workspace (use the amount required for filterSD.f)
!  maxa set as maxa parameter for filterSD
!  maxla set as maxla parameter for filterSD
!  maxu  set as maxu parameter for filterSD
!  maxiu  set as maxiu parameter for filterSD
!  mxws  set the length of ws as provided in the driver
!  mxlws  set the length of lws as provided in the driver
!  tol  tolerace (>0) for reporting inconsistencies in derivatives (eg 1.D-12)

!  Usage
!  *****
!  The user must write subroutines 'functions' and 'gradients' as for filterSD
!  Write a driver program for your problem,  but replace the call of filterSD
!  by a call of checkd (having set differencing increments in al).
!    The program will report any inconsistencies in the derivatives.
!  If the difference quotient estimate lies between the derivatives
!  at x and x+h (h is the perturbation stored in in al) then the
!  derivative is assumed to be correct. Small errors in this
!  comparison may be ignored. If no errors are reported then the
!  call of filter.. may be restored.

      dimension x(*),al(*),ws(*),lws(*)
      print *,'entering checkd'
      m1=m+1
!  set real storage map for ws
!  first maxu locations are user storage for functions and gradients
!  vectors required by checkd: two slots of length maxa for a(*)
      last1=maxu+1
      next1=last1+maxa
!  slot of length m+1 for f,c at x
      ncx0=next1+maxa
      ncx1=ncx0+1
!  slot of length m+1 for f,c at x + h.e_i
      ncxd0=ncx0+m1
      ncxd1=ncxd0+1
!  total length of ws used is
      kk=ncxd0+m
      if(kk.gt.mxws)then
        print 1,'ws not large enough: kk, mxws =',kk,mxws
        stop
      endif

!  set integer storage map for lws
!  first maxiu locations are user storage for functions and gradients
!  storage of length maxla for la(0:*)
      nla1=maxiu+1
!  total storage needed is
      ll=nla1+maxla-1
      if(ll.gt.mxlws)then
        print 1,'lws not large enough: ll, mxlws =',ll,mxlws
        stop
      endif


      call functions(n,m,x,ws(ncx0),ws(ncx1),ws,lws)
      call gradients(n,m,x,ws(last1),ws,lws)
!     print 4,'ws_0',(ws(j),j=last1,last1+7)
      do i=1,n
        xi=x(i)
        x(i)=x(i)+al(i)
        call functions(n,m,x,ws(ncxd0),ws(ncxd1),ws,lws,iflag)
        call gradients(n,m,x,ws(next1),ws,lws,iflag)
        do 10 j=0,m
          dfi=(ws(ncxd0+j)-ws(ncx0+j))/al(i)
          a_ij=aij(i,j,ws(last1),lws(nla1))
          ah_ij=aij(i,j,ws(next1),lws(nla1))
          if((dfi.ge.a_ij-tol.and.dfi.le.ah_ij+tol).or. &
            (dfi.ge.ah_ij-tol.and.dfi.le.a_ij+tol)) goto 10
          print 1,'derivative inconsistency in constraint/variable',j,i
          print *,'deriv at x, diff quotient, deriv at x+h =', &
              a_ij,dfi,ah_ij
          print *,'c at x+h, c at x',ws(ncxd0+j),ws(ncx0+j)
!         print 1,'la(0) =',la(0)
!         print 1,'c/s j pointer =',la(la(0)+j)
!         print 3,'c/s j indices',(la(k),k=la(la(0)+j),la(la(0)+j+1)-1)
!         print 4,'c/s j entries',(a(k),k=la(la(0)+j),la(la(0)+j+1)-1)
          return
10      continue
        x(i)=xi
      enddo
      print *,'exiting checkd'
1     format(A,15I5)
2     format(A,6E15.7)
3     format(A/(20I4))
4     format(A/(6E15.7))
      return
      end
