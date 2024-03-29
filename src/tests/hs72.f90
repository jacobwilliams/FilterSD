
      program hs72_driver

! program to drive the HS72 test problem, modified to give linear constraints,
! using sparse matrix format

! reference: http://apmonitor.com/wiki/uploads/Apps/hs072.apm

      implicit double precision (a-h, o-z)

      parameter (maxa=8,n=4,m=2,nm=n+m,mlp=n,mxws=30000,mxlws=5000)

      dimension a(maxa),la(0:maxa+m+2),x(n),bl(nm),bu(nm),g(n),r(nm), &
       w(nm),e(nm),ls(nm),alp(mlp),lp(mlp),ws(mxws),lws(mxlws),v(n)

      character cws

      common/wsc/kk,ll,kkk,lll,mxws_,mxlws_
      common/refactorc/mc,mxmc
      common/infoc/rgnorm,vstep,iter,npv,nfn,ngr

      data a/4.D0,2.25D0,1.D0,0.25D0,0.16D0,0.36D0,2*0.64D0/

      parameter(ainfty=1.D20,tol=1.D-12)

      write(*,*) ''
      write(*,*) 'hs72'
      write(*,*) ''

      mxws_=mxws
      mxlws_=mxlws
      kk=0
      ll=0

      do i=1,n
        la(i)=i
        la(4+i)=i
        x(i)=1.D0
        bl(i)=1.D0/((5-i)*1.D5)
        bu(i)=1.D3
      end do
      la(0)=9
      la(9)=1
      la(10)=1
      la(11)=5
      la(12)=9
      bl(5)=-ainfty
      bl(6)=-ainfty
      bu(5)=4.01D-2
      bu(6)=1.0085D-2

      kmax=4
      maxg=5
      fmin=-ainfty
      rgtol=1.D-5
      mode=0
      mxmc=25
      mxgr=100
      iprint=1
      nout=0
      v(1)=1.D0
      nv=1

!     x(1)=0.5170432D-02
!     x(2)=0.5569570D-02
!     x(3)=0.5404878D-02
!     x(4)=0.5927444D-02

!     do i=1,n
!       g(i)=1.D-2
!     end do
!     call checkg(n,x,g,r,w,ws,lws,ch,tol)
!     stop

      call glcpd(n,m,k,kmax,maxg,a,la,x,bl,bu,f,fmin,g,r,w,e,ls,alp,lp, &
       mlp,ipeq,ws,lws,ch,v,nv,rgtol,mode,ifail,mxgr,iprint,nout)

      write(nout,1)'total number of function and gradient calls =', &
        nfn,ngr
      write(nout,4)'x =',(x(i),i=1,n)
      write(nout,1)'ifail,ipeq,k =',ifail,ipeq,k
      write(nout,4)'al =',(r(abs(ls(j))),j=1,n)

1     format(A,6I5)
4     format(A/(5E15.7))
      stop
      end

      subroutine funct(n,x,f,ws,lws,cws)
      implicit double precision (a-h,o-z)
      dimension x(*),ws(*),lws(*)
      character cws(*)
!     print *,'enter funct'
!     print 4,'x =',(x(i),i=1,n)
      f=1.D0+1.D0/x(1)+1.D0/x(2)+1.D0/x(3)+1.D0/x(4)
!     print *,'f =',f
      return
      end

      subroutine grad(n,x,g,ws,lws,cws)
      implicit double precision(a-h,o-z)
      dimension x(*),g(*),ws(*),lws(*)
      character cws(*)
      do i=1,n
        g(i)=-1.D0/x(i)**2
      end do
      return
      end
