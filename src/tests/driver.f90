
      program cute_driver

      implicit double precision (a-h, o-z)

      parameter(nmax=12005,mmax=11997,nmmax=nmax+mmax,kmx=nmax,mbar=5)
!  this one for GASOIL
!     parameter(mxws=50000000,mxlws=1000000)
!  this one for PINENE
!     parameter(mxws=12000000,mxlws=1000000)
!  this one for MARINE and METHANOL
      parameter(mxws=9000000,mxlws=1000000)
!  this one otherwise
!     parameter(mxws=5000000,mxlws=1000000)
!  mxws >= 2*mxlws should hold
      dimension x(nmmax),al(nmmax),bl(nmmax),bu(nmmax),v(mbar), &
        ws(mxws),lws(mxlws)
      character cstype(mmax)
      character*10 pname
      character ch
      dimension y(6)

      common/epsc/eps,tol,emin
      common/defaultc/ainfty,ubd,mlp,mxf
      common/wsc/kk,ll,kkk,lll,mxws_,mxlws_
      common/statsc/dnorm,h,hJt,hJ,ipeq,k,itn,nft,ngt
      common/refactorc/mc,mxmc
      common/infoc/rgnorm,vstep,iter,npv,ngr,ninf
      common/pnamec/pname
      common/ngrc/mxgr
      common/maxac/maxa
      common/mxm1c/mxm1

!     data y/-0.2878796E-02, -0.1318175E+00,  0.0000000E+00,
!    *  0.3961995E+00, 0.7294942E+00, 0.2000000E+00, 0.2000000E-01,
!    * -0.5380555E+00, -0.1347677E+01, 0.1035782E+01, 0.5364218E+00,
!    *  0.8640571E+00, 0.6447348E+00, 0.1750000E+01, 0.1000000E+01,
!    *  0.4944993E+00, 0.7559711E+00, 0.5500721E-02/

!  HEART6 with 4 near zero c/s
!     data y/0.7575597E-02, -0.2474078E+00, 0.4440983E+01,
!    *  0.4072251E+01, 0.4500000E+01, 0.1080038E+01/

!  HEART6 with 5 near zero c/s
!     data y/0.8121858E-03, -0.6847485E-03, 0.2194098E+02,
!    *  0.2754532E+01, 0.2181134E+02, 0.8713753E+00/

!     data y/-0.9642272E+00, -0.4109042E+00, 0.2592764E+01,
!    *  -0.1516964E+01, -0.5244445E+00, 0.2720367E+01/

!     data y/-0.8073889E+00, -0.2342891E+00, 0.2944044E+01,
!    *  -0.1405982E+01, -0.3410198E+00, 0.3956905E+01/

!  PFIT1 with 2 near zero c/s
!     data y/1.922655D0, 1.94114D0, 1.476798D0/

      mxws_=mxws
      mxlws_=mxlws
      maxa=mxlws

      mxm1=min(m+1,n)

      call initialize(n,m,x,bl,bu,nmax,mmax,ws,lws,ws(mxlws+1))
      print *,pname

!     do i=1,n
!       x(i)=y(i)
!     end do

      v(1)=1.D0
      nv=1

      maxu=n
      maxiu=2*maxa
      maxla=maxa+m+3
      maxg=min(n,mbar)+1
      kmax=min(kmx,n)

      fmin=-ainfty
      maxit=50
      maxit=999
      iprint=1
      nout=0
      mxmc=25
      mxgr=400
      mxgr=100
      mxf=100
      mlp=100
!     ubd=1.D0
      rho=2.D1
      rho=25.D-2
      rho=1.D-1
      rho=1.D1
      rho=1.D0
      htol=1.D-6
      rgtol=1.D-5
!     tol=1.D-8

!     do i=1,n
!       al(i)=1.D-2*x(i)
!       al(i)=1.D-2
!     end do
!     call checkd(n,m,x,al,ws,lws,maxa,maxla,maxu,maxiu,
!    *  mxws,mxlws,1.D-8)
!     stop

!     do i=1,n
!       x(i)=x(i)+dble(i)*0.01
!     end do

10    continue

!     call filterSD2(n,m,x,al,f,fmin,cstype,bl,bu,ws,lws,v,nv,
!     call filter3(n,m,x,al,f,fmin,cstype,bl,bu,ws,lws,v,nv,
      call filterSD(n,m,x,al,f,fmin,cstype,bl,bu,ws,lws,v,nv, &
        maxa,maxla,maxu,maxiu,kmax,maxg,rho,htol,rgtol,maxit,iprint, &
        nout,ifail)


      if (ifail==4 .and. h>ubd) then
        ubd=11.D-1*h
         goto 10
      end if

!     print *,'cstype =',(cstype(i),i=1,m)
!     print 1,'number of function and gradient calls =',nft,ngt
!     print 4,'x =',(x(i),i=1,n)
!     print 4,'al =',(al(i),i=n+1,n+m)

      open(99,status='old',err=998)
997   continue
      read(99,*,end=999)ch
       goto 997
998   continue
      open(99)
999   continue

      nout1=99
      nout1=0
      if (ifail==0) then
        if (abs(f)<1.D5 .and. abs(f)>=1.D0) then
          write(nout1,1111)pname,n,m,f,h,rgnorm,k,itn,nft,ngt
        else
          write(nout1,2222)pname,n,m,f,h,rgnorm,k,itn,nft,ngt
        end if
      else if (ifail==3) then
        write(nout1,3333)pname,n,m,hJt,h,rgnorm,k,itn,nft,ngt,ifail
      else
        write(nout1,3333)pname,n,m,f,h,rgnorm,k,itn,nft,ngt,ifail
      end if

      stop

1     format(A,15I5)
2     format(A,6E15.7)
3     format(A/(20I4))
4     format(A/(5E15.7))

1111  format(A9,I4,I5,' ',G14.7,' ',E9.3,E10.3,2I4,2I6)
2222  format(A9,I4,I5,' ',E14.6,' ',E9.3,E10.3,2I4,2I6)
3333  format(A9,I4,I5,' ',E14.6,' ',E9.3,E10.3,2I4,2I6,' fail',I2)
      end
