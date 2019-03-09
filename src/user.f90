!Christen this file user.f

      subroutine initialize(n,m,x,bl,bu,nmax,mmax,ws,lws,iuser)
      implicit double precision (a-h,o-z)
      dimension x(*),bl(*),bu(*),ws(*),lws(*),iuser(*)
      character*10 pname
      common/pnamec/pname
      common/maxac/maxa
!     print *,'enter initialize'
!  open SIF output file
      open(file='OUTSDIF.d',unit=7)
!  read SIF output file and set up the constraints
      call csetup(7,6,n,m,x,bl,bu,nmax,ws,ws,ws, &
        bl(nmax+1),bu(nmax+1),mmax,.false.,.false.,.false.)
      do i=1,m
        bl(n+i)=bl(nmax+i)
        bu(n+i)=bu(nmax+i)
      end do
!     print 4,'bl =',(bl(i),i=1,n+m)
!     print 4,'bu =',(bu(i),i=1,n+m)
      call unames(n,pname,ws)
!  call CUTE's coordinate format into sparse vector format ...
!  ... note nonzero column indices are in ascending order
!  ... followed by n zero column indices (objective function)
      call csgr(n,m,.false.,1,x,x,mxa,maxa,ws,iuser,lws(1))
!     print 1,'mxa =',mxa
!     print 3,'i =',(iuser(k),k=1,mxa)
!     print 3,'j =',(lws(k),k=1,mxa)
!     print 4,'a(ij) =',(ws(k),k=1,mxa)
      if (3*mxa+m+3>maxa) then
        print *,'not enough space for CUTE Jacobian indices'
        print *,'mxa =',mxa,'   maxa =',maxa
        stop
      end if
      maxa=mxa
      lws(2*maxa+1)=maxa+1
      ip=3*maxa+3
      lws(ip-1)=1
      lws(ip)=n+1
      k=1
      do j=1,m
10      continue
        if (lws(k)==j) then
          k=k+1
           goto 10
        end if
        lws(ip+j)=k+n
      end do
!     print 3,'pointers',(lws(ip+i),i=-1,m)
      do i=maxa-n+1,maxa
        lws(maxa+n+1+i)=iuser(i)
      end do
      do i=1,maxa-n
        lws(2*maxa+n+1+i)=iuser(i)
      end do
!     print 3,'indices',(lws(i),i=2*maxa+2,3*maxa+1)
1     format(A,15I4)
2     format(A,5E15.7)
3     format(A/(20I4))
4     format(A/(5E15.7))
      return
      end

      subroutine functions(n,m,x,f,c,user,iuser)
      implicit double precision (a-h,o-z)
      dimension x(*),c(*),user(*),iuser(*)
!     print *,'enter functions'
!  call CUTE's objective and constraint function
!     print 4,'x =',(x(i),i=1,n)
      call cfn(n,m,x,f,m,c)
!     print *,'f =',f
!     print 4,'c =',(c(i),i=1,m)
!     do j=1,m
!       if (c(j)-c(j)/=0.D0) then
!         print 4,'x =',(x(i),i=1,n)
!         print 4,'c =',(c(i),i=1,m)
!         stop
!       end if
!     end do
4     format(A/(5E15.7))
      return
      end

      subroutine gradients(n,m,x,a,user,iuser)
      implicit double precision(a-h,o-z)
      dimension x(*),a(*),user(*),iuser(*)
      common/maxac/maxa
!     print *,'enter gradients'
!  call CUTE's sparse gradient/Jacobian evaluation function
      call csgr(n,m,.false.,1,x,x,mxa,maxa,a,iuser,iuser(maxa+1))
!     print 4,'a(ij) =',(a(k),k=1,mxa)
      do i=1,n
        user(i)=a(maxa-n+i)
      end do
      do i=maxa-n,1,-1
        a(n+i)=a(i)
      end do
      do i=1,n
        a(i)=user(i)
      end do
!     print 4,'new a(ij) =',(a(k),k=1,mxa)
4     format(A/(5E15.7))
      return
      end
