
!Christen this file l1sold.f
!ut here >>>>>>>>>>>>>>>>>>

!  Copyright (C) 2010 Roger Fletcher

      subroutine l1sold(n,m,k,kmax,maxg,a,la,x,bl,bu,f,g,r,w,e,ls, &
       alp,lp,mlp,peq,ws,lws,cws,v,nv,rgtol,ifail,iprint,nout)
      implicit double precision (a-h,r-z), integer (i-q)

!  This routine is a post-processor for qlcpd and glcpd. In the case that the
!  constraint set is infeasible (ifail=3), l1sold finds a best l1 solution of
!  the general constraints, subject to the simple bounds being satisfied.

!  Parameters are a subset of those for glcpd and must be passed through
!  unchanged. For qlcpd a dummy parameter cws must be included.

!  A modified form of Wolfe's method is used to resolve degeneracy.
!  If the solution is degenerate, there may be inactive constraints with
!  zero residual and multiplier 1. Such constraints are marked on exit by
!  setting their residual value (in r(*)) to -eps, (see common/epsc for eps)

      parameter (ainfty=1.D100)
      dimension a(*),la(*),x(*),bl(*),bu(*),g(*),r(*),w(*),e(*),ls(*), &
        alp(*),lp(*),ws(*),lws(*),v(*)
      character cws(*)
      character(len=32) spaces
      common/lcpdc/na,na1,nb,nb1,krg,krg1,kr,kr1, &
        ka,ka1,kb,kb1,kc,kc1,kd,kd1,ke,ke1,lu1,ll1
      common/epsc/eps,tol,emin
      common/infoc/vstep,iter,npv,nfn,ngr
      common/repc/sgnf,nrep,npiv,nres
      common/wsc/kk,ll,kkk,lll,mxws,mxlws
      common/refactorc/nup,nfreq
      common/alphac/alpha,rp,pj,qqj,qqj1
      logical plus

1     format(A,15I5)
2     format(A,6E15.7)
3     format(A/(15I5))
4     format(A/(5E15.7))
5     format((6E15.7))

!     if (iprint==3) print 4,'a =',(a(i),i=1,110)
      spaces='         '
      n1=n+1
      nm=n+m
      lp(1)=nm
      lev=1
      npv=0
!  collect simple bound equations
      peq=0
      do j=peq+1,n-k
        i=abs(ls(j))
        if (i<=n .and. bl(i)==bu(i)) then
          peq=peq+1
          call iexch(ls(j),ls(peq))
        end if
      end do
      gnorm=sqrt(scpr(0.D0,g,g,n))
      gtol=sgnf*gnorm
      rgtol=max(rgt0l,gtol)
      if (iprint>=1) write(nout,'(''pivots ='',I5, ''  level = 1    f ='',E16.8)')npv,f
!     print 4,'gradient =',(g(i),i=1,n)
       goto 20
!  start of major iteration
10    continue
      if (iprint>=1) write(nout,'(''pivots ='',I5, ''  level = 1    f ='',E16.8)')npv,f
!  calculate multipliers
!     print 4,'gradient =',(g(i),i=1,n)
      do i=1,nm
        w(i)=0.D0
      end do
      call fbsub(n,1,n,a,la,0,g,w,ls,ws(lu1),lws(ll1),.true.)
      call signst(n,r,w,ls)

20    continue
      if (iprint>=3) then
        write(nout,1001)'costs vector and indices', &
          (ls(j),r(abs(ls(j))),j=1,n)
!       write(nout,1000)'steepest edge coefficients',
!    *    (e(abs(ls(j))),j=1,n)
        write(nout,1)'# of bound equations and free variables = ',peq,k
      end if
!     if (iprint==3) print 4,'gradient =',(g(i),i=1,n)
!     if (iprint==3) write(nout,1001)'residual vector and indices',
!    *    (ls(j),r(abs(ls(j))),j=n1,lp(1))
!     call check1(n,lp(1),nm,k,kmax,g,a,la,x,bl,bu,r,ls,lp,ws(nb1),f,
!    *  ws,lws,cws,1,p,rp)

21    continue
!  l1 optimality test
      call optest1(peq,k,n,bl,bu,r,e,ls,rp,pj)

22    continue
      if (rp<=gtol) then
!  allow for changes to norm(g)
        gnorm=sqrt(scpr(0.D0,g,g,n))
        gtol=sgnf*gnorm
      end if

      if (rp<=gtol) then
        do j=peq+1,n-k
          i=abs(ls(j))
          if (i<=n) then
            if (ls(j)>=0) then
              x(i)=bl(i)
            else
              x(i)=bu(i)
            end if
          end if
        end do
        do i=1,n
          x(i)=max(min(x(i),bu(i)),bl(i))
        end do
        do j=n1,nm
          i=abs(ls(j))
          if (r(i)==0.D0 .and. i<=n) then
            if (ls(j)>=0) then
              x(i)=bl(i)
            else
              x(i)=bu(i)
            end if
          end if
        end do
!       write(nout,1)'# of simple bound equations = ',peq
        do j=peq+1,n-k
          i=abs(ls(j))
          if (bl(i)==bu(i)) then
            peq=peq+1
            call iexch(ls(j),ls(peq))
          end if
        end do
!       write(nout,1001)'costs vector and indices',
!    *    (ls(j),r(abs(ls(j))),j=1,n)
!       write(nout,1)'# of active equations and free variables = ',peq,k
        if (iprint>=2) then
          write(nout,*)'OPTIMAL l1 solution'
          if (iprint>=3) then
!           write(nout,1000)'x variables',(x(i),i=1,n)
            write(nout,1001)'residual vector and indices', &
              (ls(j),r(abs(ls(j))),j=n1,nm)
          end if
        end if
        ifail=0
        return
      end if

      p=abs(ls(pj))
      if (iprint>=2) write(nout,*)'CHOOSE p,pj =',ls(pj),pj,r(p)
!  compute +/- Steepest Edge search direction s in an(.)
      call tfbsub(n,a,la,p,ws(na1),ws(na1),ws(lu1),lws(ll1), &
        e(p),.true.)

        rp=scpr(0.D0,ws(na1),g,n)
        if (ls(pj)<0)rp=-rp
        if (rp*r(p)<=0.D0) then
          print 2,'3rp,r(p),rp-r(p)',rp,r(p),rp-r(p)
          r(p)=0.D0
           goto 98
        end if

      if (pj>n-k) then
        rp=-abs(r(p))
        plus=ls(pj)>=0.eqv.r(p)<0.D0
      else if (r(p)>0.D0) then
        rp=1.D0-r(p)
        plus=ls(pj)<0
      else
        rp=r(p)
        plus=ls(pj)>=0
      end if
      snorm=e(p)
!     print 4,'s (or -s if .not.plus) =',(ws(i),i=na1,na+n)

!  form At.s and denominators
      call form_Ats(n1,lp(1),n,plus,a,la,ws(na1),w,ls,snorm*tol)

!  return from degeneracy code
30    continue
      if (iprint>=3) then
!       write(nout,1000)'x variables',(x(i),i=1,n)
        write(nout,1001)'residual vector and indices', &
          (ls(j),r(abs(ls(j))),j=n1,lp(1))
        write(nout,1000)'denominators',(w(abs(ls(j))),j=n1,lp(1))
      end if
!     print 2,'slope for r(169) =',aiscpr(n,a,la,169-n,ws(na1),0.D0)
!     print 2,'slope for r(217) =',aiscpr(n,a,la,217-n,ws(na1),0.D0)
!     print *,'plus =',plus

40    continue
!  level 1 ratio tests
      amax=ainfty
      qj=0
      do 41 j=n-k+1,n
        i=ls(j)
        si=ws(na+i)
        t=abs(si)
        if (t<=tol) goto 41
        if (si>0.D0.eqv.plus) then
          z=bu(i)-x(i)
          if (abs(z)<tol) then
            z=0.D0
            x(i)=bu(i)
          else
            z=z/t
          end if
        else
          z=x(i)-bl(i)
          if (abs(z)<tol) then
            z=0.D0
            x(i)=bl(i)
          else
            z=z/t
          end if
        end if
        if (z>amax) goto 41
        amax=z
        qj=j
41    continue
      if (pj<=n-k .and. r(p)<0.D0 .and. bu(p)-bl(p)<amax) then
        amax=bu(p)-bl(p)
        qj=pj
      end if
      alpha=amax
      do 44 j=n1,lp(1)
        i=abs(ls(j))
        wi=w(i)
        if (wi==0.D0) goto 44
        ri=r(i)
        if (ri==-eps) then
          if (bl(i)==bu(i) .and. wi<0.D0) then
            alpha=0.D0
            q=i
            qj=j
             goto 45
          end if
           goto 44
        else if (ri<0.D0) then
          if (wi>0.D0) goto 44
          z=(ri-tol)/wi
        else if (wi>0.D0) then
          z=(ri+tol)/wi
        else
          z=(bl(i)-bu(i)+ri-tol)/wi
        end if
        if (z>=alpha) goto 44
        alpha=z
        qj=j
44    continue
      q=abs(ls(qj))
      if (qj>n) then
        if (r(q)<0.D0.eqv.w(q)<0.D0) then
          alpha=r(q)/w(q)
        else
          alpha=(bl(q)-bu(q)+r(q))/w(q)
        end if
      end if
45    continue
      if (iprint>=2) then
        write(nout,2)'r(q),w(q) =',r(q),w(q)
        write(nout,*)'alpha =',alpha,'   q =',q
      end if

      if (alpha==0.D0 .and. pj<=n-k) then
        if (iprint>=2) &
          write(nout,*)'degeneracy block at level 1'
        if (w(q)<0.D0 .and. r(q)/=-eps) then
          w(q)=-w(q)
          ls(qj)=-ls(qj)
        end if
        alp(1)=f
        plev=n
        do j=n1,lp(1)
          i=abs(ls(j))
          if (r(i)==-eps) then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            r(i)=-1.D0
          else if (r(i)==0.D0) then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            r(i)=1.D0
          end if
        end do
        lp(2)=plev
        lev=2
!       write(nout,1001)'costs vector and indices',
!    *    (ls(j),r(abs(ls(j))),j=1,n)
!       write(nout,1)'# of bound equations and free variables = ',peq,k
        if (iprint>=1) write(nout,'(''pivots ='',I5,''     level = 2'', ''    f ='',E16.8)')npv,f
         goto 86
      end if

      if (alpha>0.D0) then
        ff=f
        f=f+alpha*rp
        if (f>=ff) then
          if (pj>n-k .or. r(p)<0.D0) then
            r(p)=0.D0
          else
            r(p)=1.D0
          end if
           goto 20
        end if
        if (plus) then
          call mysaxpy(alpha,ws(na1),x,n)
        else
          call mysaxpy(-alpha,ws(na1),x,n)
        end if
!  update r for inactive c/s
        do 61 j=n1,lp(1)
          i=abs(ls(j))
          if (w(i)==0.D0) goto 61
          ri=r(i)-alpha*w(i)
          if (abs(ri)<=tol)ri=0.D0
          if (r(i)<0.D0 .and. ri>=0.D0) then
!  remove contribution to gradient
            call saipy(sign(1.D0,dble(ls(j))),a,la,i-n,g,n)
            call newg
          end if
          if (w(i)<0.D0) then
            ro=(bu(i)-bl(i))-ri
            if (abs(ro)<=tol)ro=0.D0
            if (ro<ri) then
              ri=ro
!             w(i)=-w(i)
              ls(j)=-ls(j)
            end if
          end if
          if (ri==0.D0 .and. i<=n) then
            if (ls(j)>=0) then
              x(i)=bl(i)
            else
              x(i)=bu(i)
            end if
          end if
          r(i)=ri
61      continue
      end if

70    continue
      if (qj/=pj) then
!  pivot interchange
        if (iprint>=2) write(nout,*)'replace',p,' by',q
        call pivot(p,q,n,nm,a,la,e,ws(lu1),lws(ll1),ifail,npv)
        if (ifail>=1) then
          if (iprint>=1) write(nout,*)'near singularity in pivot (1)'
           goto 98
        end if
        if (pj>n-k) then
          r(p)=x(p)-bl(p)
          rpu=bu(p)-x(p)
          if (rpu<r(p)) then
            r(p)=rpu
            ls(pj)=-p
          else
            ls(pj)=p
          end if
          if (r(p)<=tol)r(p)=0.D0
        else if (r(p)>0.D0) then
          call saipy(-sign(1.D0,dble(ls(pj))),a,la,p-n,g,n)
          call newg
          r(p)=-alpha
        else
          rpu=bu(p)-bl(p)-alpha
          if (abs(rpu)<=tol)rpu=0.D0
          if (alpha<=rpu) then
            r(p)=alpha
          else
            r(p)=rpu
            ls(pj)=-ls(pj)
          end if
        end if
        if (pj>n-k) then
          k=k-1
          call iexch(ls(pj),ls(n-k))
          pj=n-k
        end if
        call iexch(ls(pj),ls(qj))
        if (q<=n .and. bl(q)==bu(q)) then
          peq=peq+1
          call iexch(ls(pj),ls(peq))
        end if
         goto 10
      end if
!  opposite bound comes active
!     r(p)=-rp
      if (pj<=n-k) then
        ls(pj)=-ls(pj)
         goto 10
      end if
!  free variable reaches its bound
      if (r(p)<0.D0)ls(pj)=-p
      k=k-1
      call iexch(ls(pj),ls(n-k))
       goto 10

!  recursive code for resolving degeneracy (Wolfe's method)
80    continue
!  calculate multipliers
      call fbsub(n,1,n,a,la,0,g,w,ls,ws(lu1),lws(ll1),.true.)
      call signst(n,r,w,ls)
!  reset multiplier loop
82    continue
      if (iprint>=3) then
         write(nout,1001)'costs vector and indices', &
          (ls(j),r(abs(ls(j))),j=1,n)
!       write(nout,1000)'steepest edge coefficients',
!    *    (e(abs(ls(j))),j=1,n)
        write(nout,1)'# of bound equations and free variables = ',peq,k
      end if

84    continue
!     call check1(n,lp(1),nm,k,kmax,g,a,la,x,bl,bu,r,ls,lp,ws(nb1),f,
!    *  ws,lws,cws,lev,p,rp)

85    continue
      call optest1(peq,k,n,bl,bu,r,e,ls,rp,pj)

      if (rp<=gtol .or. pj>n-k) then
        if (iprint>=2) write(nout,*)'return to level 1'
        do j=n1,lp(2)
          i=abs(ls(j))
          if (r(i)<0.D0) then
            r(i)=-eps
          else
            r(i)=0.D0
          end if
        end do
        lev=1
        f=alp(1)
         goto 22
      end if
      p=abs(ls(pj))
      if (iprint>=2) write(nout,*)'CHOOSE p,pj =',p,pj
!  compute +/- Steepest Edge (SE) search direction s in an(.)
      call tfbsub(n,a,la,p,ws(na1),ws(na1),ws(lu1),lws(ll1), &
        e(p),.true.)

        rp=scpr(0.D0,ws(na1),g,n)
        if (ls(pj)<0)rp=-rp
        if (rp*r(p)<=0.D0) then
          print 2,'4rp,r(p),rp-r(p)',rp,r(p),rp-r(p)
          do j=n1,lp(2)
            i=abs(ls(j))
            if (r(i)<0.D0) then
              r(i)=-eps
            else
              r(i)=0.D0
            end if
          end do
          f=alp(1)
           goto 98
        end if

      if (r(p)>0.D0) then
        rp=1.D0-r(p)
        plus=ls(pj)<0
      else
        rp=r(p)
        plus=ls(pj)>=0
      end if
      snorm=e(p)
!     print 4,'s (or -s if .not.plus) =',(ws(i),i=na1,na+n)

!  form At.s and denominators
      call form_Ats(n1,lp(lev),n,plus,a,la,ws(na1),w,ls,snorm*tol)
!     print 2,'slope for r(169) =',aiscpr(n,a,la,169-n,ws(na1),0.D0)
!     print 2,'slope for r(217) =',aiscpr(n,a,la,217-n,ws(na1),0.D0)
!     print *,'plus =',plus
86    continue
      if (iprint>=3) then
        write(nout,1001)'residual vector and indices', &
          (ls(j),r(abs(ls(j))),j=n1,lp(lev))
        write(nout,1000)'denominators',(w(abs(ls(j))),j=n1,lp(lev))
      end if
88    continue
!  ratio test at higher levels
      alpha=ainfty
      qj=0
      do 90 j=n1,lp(lev)
        i=abs(ls(j))
        wi=w(i)
        if (wi==0.D0) goto 90
        ri=r(i)
        if (ri==-eps) then
          if (bl(i)==bu(i) .and. wi<0.D0) then
            alpha=0.D0
            q=i
            qj=j
             goto 91
          end if
           goto 90
        else if (ri<0.D0) then
          if (wi>0.D0) goto 90
          z=(ri-tol)/wi
        else if (wi>0.D0) then
          z=(ri+tol)/wi
        else
           goto 90
!         if (bl(i)<bu(i)) goto 90
!         z=(ri-2.D0-tol)/wi
        end if
        if (z>=alpha) goto 90
        alpha=z
        qj=j
90    continue
      if (qj==0) then
        do j=n1,lp(lev)
          i=abs(ls(j))
          if (r(i)<0.D0) then
            r(i)=-eps
          else
            r(i)=0.D0
          end if
        end do
        call form_Ats(lp(lev)+1,lp(lev-1),n,plus,a,la,ws(na1), &
          w,ls,snorm*tol)
        lev=lev-1
        f=alp(lev)
        if (iprint>=2) write(nout,*)'UNBOUNDED:   p =',p, &
          '   return to level',lev
        if (lev>1) goto 86
        if (iprint>=3) then
          write(nout,1001)'costs vector and indices', &
            (ls(j),r(abs(ls(j))),j=1,n)
        write(nout,1)'# of bound equations and free variables = ',peq,k
        end if
         goto 30
      end if
      q=abs(ls(qj))
      alpha=r(q)/w(q)
!     if (alpha<0.D0) then
!       r(q)=r(q)-2.D0
!       w(q)=-w(q)
!       ls(qj)=-ls(qj)
!     end if
!     print *,'alpha =',alpha
91    continue
      if (iprint>=2) then
        write(nout,*)'alpha =',alpha,'   p =',p,'   q =',q
        write(nout,2)'r(p),r(q),w(q) =',r(p),r(q),w(q)
      end if

      if (alpha==0.D0) then
        if (iprint>=2) write(nout,1) &
          'degeneracy block at level',lev
        if (lev+2>mlp) then
          ifail=5
          return
        end if
        if (w(q)<0.D0 .and. r(q)/=-eps) then
          w(q)=-w(q)
          ls(qj)=-ls(qj)
        end if
!       r(q)=0.D0
!       alp(lev)=f
        plev=n
        do j=n1,lp(lev)
          i=abs(ls(j))
          if (r(i)==-eps) then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            r(i)=-1.D0
          else if (r(i)==0.D0) then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            r(i)=1.D0
          end if
        end do
        lev=lev+1
        lp(lev)=plev
        if (iprint>=2) write(nout,*) 'degeneracy: increase level to ',lev
        if (iprint>=1) write(nout,'(''pivots ='',I5,A,''level ='',I2, ''    f ='',E16.8)')npv,spaces(:3*lev-1),lev,f
           goto 86
      end if
!  update r and f
      if (alpha>0.D0) then
!       ff=f
!       f=f+alpha*rp
!       if (f>=ff) then
!         if (r(p)>0.D0) then
!           r(p)=1.D0
!         else
!           r(p)=0.D0
!         end if
!          goto 85
!       end if
        do 92 j=n1,lp(lev)
          i=abs(ls(j))
          if (w(i)==0.D0) goto 92
          ri=r(i)-alpha*w(i)
          if (abs(ri)<=tol)ri=0.D0
          if (r(i)<0.D0 .and. ri>=0.D0) then
!  remove contribution to gradient
            call saipy(sign(1.D0,dble(ls(j))),a,la,i-n,g,n)
            call newg
          end if
!         if (w(i)<0.D0 .and. bl(i)==bu(i)) then
!           ro=2.D0-ri
!           if (abs(ro)<=tol)ro=0.D0
!           if (ro<ri) then
!             ri=ro
!             w(i)=-w(i)
!             ls(j)=-ls(j)
!           end if
!         end if
          r(i)=ri
92      continue
      end if
      if (iprint>=2) write(nout,*)'replace',p,' by',q
      call pivot(p,q,n,nm,a,la,e,ws(lu1),lws(ll1),ifail,npv)
      if (ifail>=1) then
        if (ifail>=2) return
!       call iexch(ls(pj),ls(qj))
        if (iprint>=1) write(nout,*)'near singularity in pivot (4)'
         goto 98
      end if
      if (r(p)>0.D0) then
!  add contribution to g from c/s p
        call saipy(-sign(1.D0,dble(ls(pj))),a,la,p-n,g,n)
        call newg
        r(p)=-alpha
      else
        rpu=2.D0-alpha
        if (alpha<=rpu .or. bl(p)<bu(p)) then
          r(p)=alpha
        else
          r(p)=rpu
          ls(pj)=-ls(pj)
        end if
      end if
      if (abs(r(p))<=tol)r(p)=0.D0
!  exchange a constraint
      call iexch(ls(pj),ls(qj))
      if (q<=n .and. bl(q)==bu(q)) then
        peq=peq+1
        call iexch(ls(pj),ls(peq))
      end if
      if (iprint>=1) write(nout,'(''pivots ='',I5,A,''level ='',I2, ''    f ='',E16.8)')npv,spaces(:3*lev-1),lev,f
       goto 80
!  restart sequence
98    continue
      return
1000  format(a/(e16.5,4e16.5))
1001  format(a/(i4,1x,e12.5,4(i4,1x,e12.5)))
!1000 format(a/(e18.8,3e19.8))
!1001 format(a/(i3,1x,e14.8,3(i4,1x,e14.8)))
      end

      subroutine optest1(peq,k,n,bl,bu,r,e,ls,rp,pj)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension bl(*),bu(*),r(*),e(*),ls(*)
      rp=0.D0
      do 1 j=peq+1,n-k
        i=abs(ls(j))
        ri=-r(i)/e(i)
        if (i<=n) then
          if (ri<=rp) goto 1
        else if (ri<0.D0) then
          ri=(r(i)-1.D0)/e(i)
          if (ri<=rp) goto 1
        else if (ri>rp) then
          if (bl(i)==bu(i)) then
            ri=(-r(i)-1.D0)/e(i)
            if (ri<=rp) goto 1
            r(i)=-r(i)
            ls(j)=-ls(j)
          end if
        else
           goto 1
        end if
        rp=ri
        pj=j
1     continue
!  additional test for free variables
      do j=n-k+1,n
        i=abs(ls(j))
        ri=abs(r(i))/e(i)
        if (ri>=rp) then
          rp=ri
          pj=j
        end if
      end do
      return
      end

      subroutine check1(n,nm,nmi,k,kmax,g,a,la,x,bl,bu,r,ls,lp,an,f, &
        ws,lws,cws,lev,p,alp2)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension g(*),a(*),la(*),x(*),bl(*),bu(*),r(*),ls(*),lp(*), &
        an(*),ws(*),lws(*)
      character cws(*)
      common/noutc/nout
      common/epsc/eps,tol,emin
!     print *,'ENTER check1'
      e=0.D0
      if (lev==1) then
        do j=n+1,nm
          i=abs(ls(j))
          if (i<=n) then
            s=x(i)
          else
            s=aiscpr(n,a,la,i-n,x,0.D0)
          end if
          if (ls(j)>0) then
!           print *,'i,s,r(i),bl(i)',i,s,r(i),bl(i)
            s=r(i)-s+bl(i)
          else
            s=r(i)+s-bu(i)
          end if
          if (abs(s)<=tol*max(1.D0,abs(r(i))))s=0.D0
          if (abs(s)>e) then
            e=abs(s)
            ie=i
          end if
        end do
      else
        do j=n+1,lp(2)
          i=abs(ls(j))
          if (i<=n) then
            s=x(i)
          else
            s=aiscpr(n,a,la,i-n,x,0.D0)
          end if
!         print 2,'s,bl(i),bu(i) =',s,bl(i),bu(i)
          if (ls(j)>0) then
            s=-s+bl(i)
          else
            s=s-bu(i)
          end if
          if (abs(s)<=tol)s=0.D0
          if (abs(s)>e) then
            e=abs(s)
            ie=i
          end if
        end do
      end if
      if (e>tol) write(nout,*)'inactive c/s residual error = ',e,ie
!     if (e>tol) stop
      if (e>1.D-6) print 2,'r(ie)',r(ie)
      if (e>1.D-6) stop
      if (lev==1) then
        ff=0.D0
        do j=n+1,nm
          i=abs(ls(j))
          if (r(i)<0.D0)ff=ff-r(i)
        end do
        e=abs(ff-f)
        if (e>tol*max(1.D0,abs(f))) write(nout,*)'function error = ',e, &
          '   f(x) =',ff
        if (e>tol*max(1.D0,abs(f))) stop
      end if
!       print 4,'g =',(g(i),i=1,n)
!       print 4,'an =',(an(i),i=1,n)
!       err=0.D0
!       do i=1,n
!         err=err+abs(g(i)-an(i))
!       end do
!       print 2,'check err =',err
      do i=1,n
        an(i)=0.D0
      end do
      do j=n+1,nm
        i=abs(ls(j))
        if (r(i)<0.D0) then
          if (i>n) then
            call saipy(-sign(1.D0,dble(ls(j))),a,la,i-n,an,n)
          else
            an(i)=an(i)-sign(1.D0,dble(ls(j)))
          end if
        end if
      end do
      gnm=sqrt(scpr(0.D0,an,an,n))
      e=0.D0
      do j=1,n
!       write(nout,*)'an =',(an(i),i=1,n)
        i=abs(ls(j))
        s=sign(1.D0,dble(ls(j)))
        if (i<=n) then
          an(i)=an(i)-s*r(i)
          if (j>n-k) then
            s=max(0.D0,bl(i)-x(i),x(i)-bu(i))
          else if (ls(j)>0) then
            s=x(i)-bl(i)
          else
            s=bu(i)-x(i)
          end if
        else
          call saipy(-s*r(i),a,la,i-n,an,n)
          if (ls(j)>0) then
            s=aiscpr(n,a,la,i-n,x,-bl(i))
          else
            s=-aiscpr(n,a,la,i-n,x,-bu(i))
          end if
        end if
        if (abs(s)>e) then
          e=abs(s)
          ie=i
        end if
      end do
      if (e>tol) write(nout,*)'active c/s residual error = ',e,ie
!     if (e>tol) stop
!     if (e>1.D-4) print 4,'x =',(x(i),i=1,n)
      if (e>1.D-4) stop
      e=0.D0
      do j=1,n
        if (abs(an(j))>e) then
          e=abs(an(j))
          ie=ls(j)
          je=j
        end if
      end do
      if (e>gnm*tol) write(nout,*)'KT condition error = ',e,je,ie,gnm
!     if (e>gnm*tol) write(nout,4)'KT cond_n errors = ',(an(i),i=1,n)
!     if (e>gnm*tol) stop
      if (e>1.D-4) stop
1     format(A,10I5)
2     format(A,5E15.7)
4     format(A/(5E15.6))
      return
      end
