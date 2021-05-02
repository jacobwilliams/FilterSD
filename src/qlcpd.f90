!Christen this file qlcpd.f
!ut here >>>>>>>>>>>>>>>>>

!  Copyright (C) 2010 Roger Fletcher

!  Current version dated 18 April 2013

!  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
!  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
!  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

      subroutine qlcpd(n,m,k,kmax,maxg,a,la,x,bl,bu,f,fmin,g,r,w,e,ls, &
        alp,lp,mlp,peq,ws,lws,v,nv,linear,rgtol,m0de,ifail,mxgr,iprint, &
        nout)
      implicit double precision (a-h,r-z), integer (i-q)
      logical linear

!  This routine finds a KT point for the Quadratic LCP (QP) problem

!       minimize    f(x) = ct.x + xt.G.x/2

!       subject to  l <= [I : A]t.x <= u                  (t = transpose)

!  where x and c are n-vectors, G is a symmetric n*n matrix, and A is an
!  n*m matrix. If G is also positive semi-definite then the KT point is a
!  global solution, else usually a local solution. The method may also be
!  used efficiently to solve an LP problem (G=0). A recursive form of an
!  active set method is used, using Wolfe's method to resolve degeneracy.
!  A limited memory sweep method is used for minimization in the null space.
!  Matrix information is made available and processed by calls to external
!  subroutines. Details of these are given in an auxiliary file named either
!  'denseL.f' or 'schurQR.f'. (schurQR.f is a more recent replacement for
!  the file sparseL.f)

!  parameter list  (variables in a line starting with C must be set on entry)
!  **************

!  n     number of variables
!  m     number of general constraints (columns of A)
!  k     dimension of the null space obtained by eliminating the active
!        constraints (only to be set if mode>=2). The number of constraints in
!        the active set is n-k
!  kmax  maximum value of k (kmax <= n)
!  maxg  max number of reduced gradient vectors stored in sweep method:
!        (1 < maxg <= kmax+1), typically maxg = min(6,kmax+1)
!  a(*)  storage of reals associated with c and A. This storage may be provided
!        in either dense or sparse format. Refer to either denseA.f or sparseA.f
!        for information on how to set a(*) and la(*).
!  la(*) storage of integers associated with c and A
!  x(n)  contains the vector of variables. Initially an estimate of the solution
!        must be set, replaced by the solution (if it exists) on exit.
!  bl(n+m)  vector of lower bounds for variables and general constraints
!  bu(n+m)  vector of upper bounds (use numbers less than about 1.e30, and
!        where possible supply realistic bounds on the x variables)
!  f     returns the value of f(x) when x is a feasible solution
!        Otherwise f stores the sum of constraint infeasibilities
!  fmin  set a strict lower bound on f(x) (used to identify an unbounded LCP)
!  g(n)  returns the gradient vector of f(x) when x is feasible
!  r(n+m) workspace: stores constraint residuals (or multipliers if the
!        constraint is active). The sign convention is such that these are
!        nonnegative at a solution (except multipliers of equality constraints)
!  w(n+m) workspace: stores denominators for ratio tests
!  e(n+m) stores steepest-edge normalization coefficients: if mode>2 then
!        information in this vector from a previous call should not be changed.
!        (In mode 3 these values provide approximate coefficients)
!  ls(n+m) stores indices of the active constraints in locations 1:n and of
!        the inactive constraints in locations n+1:n+m. The simple bounds
!        on the variables are indexed by 1:n and the general constraints by
!        n+1:n+m. The sign of ls(j) indicates whether the lower bound (+) or
!        the upper bound (-) of constraint ls(j) is currently significant.
!        Within the set of active constraints, locations 1:peq store the indices
!        of any equality constraints, locations peq+1:n-k store the indices of
!        any inequality constraints, and locations n-k+1 store the indices of
!        any free variables (variables not on a bound, which are used to
!        parametrise the null space: ls(j) is always positive in this range)
!          If mode>=2, the first n-k elements of ls must be set on entry
!  alp(mlp) workspace associated with recursion
!  lp(mlp)  list of pointers to recursion information in ls
!  mlp   maximum number of levels of recursion allowed (mlp>2: typically
!        mlp=20 would usually be adequate but mlp=m is an upper bound)
!  peq   pointer to the end of equality constraint indices in ls
!  ws(*) real workspace for gdotx (see below), qlcpd and denseL.f (or schurQR.f)
!          Set the total number in mxws (see "Common" below).
!  lws(*) integer workspace for gdotx, qlcpd and denseL.f (or schurQR.f).
!          Set the total number in mxlws (see "Common" below).
!        The storage maps for ws and lws are set by the routine stmapq below
!  v(maxg) set nv estimates of the eigenvalues of the reduced Hessian of f(x)
!          (for example from a previous run of qlcpd). Set nv=1 and v(1)=1.D0
!          in absence of other information. New values of v are left on exit
!  nv    Number of estimates in v
!  linear  (logical) set linear = .true. iff the problem is an LP
!  rgtol required accuracy in the reduced gradient l2 norm: it is advisable not
!        to seek too high accuracy - rgtol may be increased by the code if it
!        is deemed to be too small, see the definition of sgnf below
!  m0de  mode of operation (larger numbers imply extra information):
!          0 = cold start (no other information available, takes simple
!                bounds for the initial active set)
!          1 = as 0 but includes all equality constraints in initial active set
!          2 = user sets n-k active constraint indices in ls(j), j=1,..,n-k.
!                For a general constraint the sign of ls(j) indicates which
!                bound to use. For a simple bound the current value of x is used
!          3 = takes active set and other information from a previous call.
!                Steepest edge weights are approximated using previous values.
!          4 = as 3 but it is also assumed that A is unchanged so
!                that factors of the basis matrix stored in ws and lws are valid
!                (changes in the vectors c, l, u and the matrix G are allowed)
!        A local copy (mode) of m0de is made and may be changed by qlcpd
!  ifail   outcome of the process
!              0 = solution obtained
!              1 = unbounded problem detected (f(x)<=fmin would occur)
!              2 = bl(i) > bu(i) for some i
!              3 = infeasible problem detected in Phase 1
!              4 = line search cannot improve f (possibly increase rgtol)
!              5 = mxgr gradient calls exceeded (this test is only carried
!                    out at the start of each iteration)
!              6 = incorrect setting of m, n, kmax, maxg, mlp, mode or tol
!              7 = not enough space in ws or lws
!              8 = not enough space in lp (increase mlp)
!              9 = dimension of reduced space too large (increase kmax)
!             10 = maximum number of unsuccessful restarts taken
!          >  10 = crash in pivot call
!  mxgr  maximum number of gradient evaluations
!  iprint  switch for diagnostic printing (0 = off, 1 = summary,
!                 2 = scalar information, 3 = verbose)
!  nout  channel number for output

!  Common
!  ******
!  User information about the lengths of ws and lws is supplied to qlcpd in
!    common/wsc/kk,ll,kkk,lll,mxws,mxlws
!  kk and ll refer to the length of ws and lws needed by gdotx.
!  kkk and lll are the numbers of locations used by qlcpd and are set by qlcpd.
!  the rest of ws and lws is used by the files denseL.f or schurQR.f
!  mxws and mxlws must be set to the total length of ws and lws available: a
!  message will be given if more storage is needed.

!  User subroutine
!  ***************
!  The user must provide a subroutine to calculate the vector v := G.x from a
!  given vector x. The header of the routine is
!            subroutine gdotx(n,x,ws,lws,v)
!            implicit double precision(a-h,o-z)
!            dimension x(*),ws(*),lws(*),v(*)
!  In the case that linear is .true. the subroutine is not called by qlcpd and
!  a dummy routine may be supplied. Otherwise the user may use the parameters
!  ws and lws (see above) for passing real or integer arrays relating to G.
!  Locations ws(1),...,ws(kk) are available to the user for storing any real
!  information to be used by gdotx. Likewise locations lws(1),...,lws(ll) are
!  available for storing any integer information. Default values are kk=ll=0.
!  Any other setting is made by changing  common/wsc/kk,ll,kkk,lll,mxws,mxlws

!  Tolerances and accuracy
!  ***********************
!  qlcpd uses tolerance and accuracy information stored in
!     common/epsc/eps,tol,emin
!     common/repc/sgnf,nrep,npiv,nres
!     common/refactorc/mc,mxmc
!     common/infoc/rgnorm,vstep,iter,npv,ngr
!  eps must be set to the machine precision (unit round-off) and tol is a
!  tolerance such that numbers whose absolute value is less than tol are
!  truncated to zero. This tolerance strategy in the code assumes that the
!  problem is well-scaled and a pre-scaling routine cscale is supplied in
!  denseA.f or sparseA.f. The parameter sgnf is used to measure the maximum
!  allowable relative error in gradient values. If at any stage the accuracy
!  requirement rgtol < sgnf*rgnorm then rgtol is increased to sgnf*rgnorm
!    The code allows one or more refinement steps after the
!  calculation has terminated, to improve the accuracy of the solution,
!  and a fixed number nrep of such repeats is allowed. However the code
!  terminates without further repeats if no more than npiv pivots are taken.
!    In case of any breakdown, the code is restarted, usually in mode 2.
!  The maximum number of unsuccessful restarts allowed is set in nres.
!    The basis matrix may be refactorised on occasions, for example to prevent
!  build-up of round-off in the factors or (when using schurQR.f) to limit
!  the growth in the Schur complement. The maximum interval between
!  refactorizations (or size of Schur complement) is set in mxmc.
!    Default values are set in block data but can be reset by the user.
!    infoc returns information about the progress of the method: rgnorm is the
!  norm of the reduced gradient on exit, and vstep is the length of the vertical
!  step in the warm start process. iter is the total number of iterations taken,
!  npv is the number of pivots, and ngr is the number of calls of gdotx.

      parameter (ainfty=1.D100)
      dimension a(*),x(*),bl(*),bu(*),g(*),r(*),w(*),e(*),alp(*),ws(*), &
        la(*),ls(*),lp(*),lws(*),v(*)
      character(len=32) spaces
      common/lcpdc/na,na1,nb,nb1,krg,krg1,kr,kr1, &
        ka,ka1,kb,kb1,kc,kc1,kd,kd1,ke,ke1,lu1,ll1
      common/epsc/eps,t0l,emin
      common/infoc/rgnorm,vstep,iter,npv,ngr
      common/repc/sgnf,nrep,npiv,nres
      common/wsc/kk,ll,kkk,lll,mxws,mxlws
      common/refactorc/mc,mxmc
      common/alphac/alpha,rp,pj,qqj,qqj1
      logical plus

1     format(A,15I5)
2     format(A,6E15.7)
3     format(A/(15I5))
4     format(A/(5E15.7))
5     format((6E15.7))
6     format(A,I5,2E15.7)

      spaces='         '
      mode=m0de
      tol=t0l
      iter=0
      npv=0
      if (m<0 .or. n<=0 .or. mlp<2 .or. mode<0 .or. mode>4 .or.  &
        kmax<0 .or. (kmax>0 .and. maxg<=1) .or. tol<=D0) then
        ifail=6
        return
      end if
      rgt0l=rgtol
      n1=n+1
      nm=n+m
      nmi=nm
      ngr=0
      nv0=nv
      if (iprint>=3) then
        write(nout,1000)'lower bounds',(bl(i),i=1,nm)
        write(nout,1000)'upper bounds',(bu(i),i=1,nm)
      end if
      irep=0
      ires=0
      mres=0
      bestf=ainfty
      do i=1,nm
        t=bu(i)-bl(i)
        if (t<-tol) then
          ifail=2
          return
        end if
        if (t<=tol .and. t>0.D0) then
          bl(i)=5.D-1*(bl(i)+bu(i))
          bu(i)=bl(i)
        end if
      end do
      vmax=0.D0
      do i=1,n
        x(i)=min(bu(i),max(bl(i),x(i)))
        vmax=max(vmax,bu(i)-bl(i))
      end do
      if (mode<=2) then
        call stmapq(n,nm,kmax,maxg)
        if (mode==0) then
          nk=0
        else if (mode==1) then
!  collect equality c/s
          nk=0
          do i=1,nm
            if (bu(i)==bl(i)) then
              nk=nk+1
              ls(nk)=i
            end if
          end do
!         write(nout,*)'number of eqty c/s =',nk
        else
          nk=n-k
        end if
      end if
!  restarts loop
7     continue
      lp(1)=nm
      lev=1
      if (mode<=3) then
!  set up factors of basis matrix and permutation vectors
        ifail=mode
        call start_up(n,nm,nmi,a,la,nk,e,ls,ws(lu1),lws(ll1),mode,ifail)
        if (ifail>0) return
      end if
8     continue
      peq=0
      ig=0
!  refinement step loop
      mpiv=iter+npiv
      ninf=0
      do i=1,n
        g(i)=0.D0
      end do
      if (mode>0) then
        call warm_start(n,nm,a,la,x,bl,bu,r,ls,ws(lu1), &
          lws(ll1),ws(na1),vstep)
!       print *,'vstep,vmax',vstep,vmax
        if (vstep>2.D0*vmax) then
          mpiv=0
          mode=0
          nk=0
          do i=1,n
            x(i)=min(bu(i),max(bl(i),x(i)))
          end do
           goto 7
        end if
        if (vstep>tol)mpiv=0
      end if
      k=0
!  collect free variables
      do j=n,1,-1
        i=abs(ls(j))
        if (i<=n .and. x(i)>bl(i) .and. x(i)<bu(i)) then
          call iexch(ls(j),ls(n-k))
          k=k+1
        end if
      end do
      if (mode==0) then
        do j=1,n-k
          i=ls(j)
          if (x(i)==bu(i))ls(j)=-i
        end do
        lp(1)=n
         goto 9
      end if
      phase=0
!  move inactive general c/s to the end
      do j=nm,n1,-1
        i=abs(ls(j))
        if (i>n) then
          call iexch(ls(j),ls(lp(1)))
          lp(1)=lp(1)-1
        end if
      end do
      call residuals(n,n1,lp(1),a,la,x,bl,bu,r,ls,f,g,ninf)
      if (ninf>0) then
        gnorm=sqrt(dble(ninf))
        gtol=sgnf*gnorm
        rgtol=max(rgt0l,gtol)
         goto 15
      end if
9     continue
!  enter phase 1
      phase=1
!  collect active equality c/s
      do j=1,n-k
        i=abs(ls(j))
        if (bu(i)==bl(i)) then
          peq=peq+1
          call iexch(ls(j),ls(peq))
        end if
      end do
      call residuals(n,lp(1)+1,nm,a,la,x,bl,bu,r,ls,f,g,ninf)
      lp(1)=nm
      if (ninf>0) then
        gnorm=sqrt(scpr(0.D0,g,g,n))
        gtol=sgnf*gnorm
        rgtol=max(rgt0l,gtol)
         goto 15
      end if
10    continue
      phase=2
      if (iprint>=1) write(nout,*)'FEASIBILITY OBTAINED at level 1'
      n_inf=0
      call setfg2(n,linear,a,la,x,f,g,ws,lws)
      fbase=f
      f=0.D0
      ngr=ngr+1
!     write(nout,4)'g =',(g(i),i=1,n)
      call newg
      gnorm=sqrt(scpr(0.D0,g,g,n))
      gtol=sgnf*gnorm
      rgtol=max(rgt0l,gtol)
      ig=0
      if (iprint>=1) write(nout,'(''pivots ='',I5, ''  level = 1    f ='',E16.8)')npv,fbase
       goto 16
!  start of major iteration
15    continue
      if (iprint>=1) then
        if (ninf==0) then
          if (k>0) then
!           write(nout,'(''pivots ='',I5,
!    *        ''  level = 1    df ='',E16.8,''   k ='',I4)')npv,f,k
            write(nout,'(''pivots ='',I5, ''  level = 1    df ='',E16.8,''   rg ='',E12.4, ''  k ='',I4)')npv,f,rgnorm,k
          else
            write(nout,'(''pivots ='',I5, ''  level = 1    f ='',E16.8)')npv,fbase+f
          end if
        else if (phase==0) then
          write(nout,'(''pivots ='',I5,''  level = 1    f ='', E16.8,''   ninfb ='',I4)')npv,f,ninf
        else
          write(nout,'(''pivots ='',I5,''  level = 1    f ='', E16.8,''   ninf ='',I4)')npv,f,ninf
        end if
      end if
16    continue
!  calculate multipliers
      do i=1,nm
        w(i)=0.D0
      end do
!     write(nout,4)'g =',(g(i),i=1,n)
      call fbsub(n,1,n,a,la,0,g,w,ls,ws(lu1),lws(ll1),.true.)
      call signst(n,r,w,ls)
!  opposite bound or reset multiplier loop
20    continue
      if (iprint>=3) then
        write(nout,1001)'costs vector and indices', &
          (ls(j),r(abs(ls(j))),j=1,n)
!       write(nout,1000)'steepest edge coefficients',
!    *    (e(abs(ls(j))),j=1,n)
        if (peq>0 .or. k>0) write(nout,1) &
          '# active equality c/s and free variables = ',peq,k
      end if
!     if (iphase<=1)fbase=0.D0
!     call checkq(n,lp(1),nmi,kmax,g,a,la,x,bl,bu,r,ls,ws(nb1),fbase+f,
!    *  ws,lws,ninf,peq,k,1,p,rp,linear)

21    continue
      call optest(peq+1,n-k,r,e,ls,rp,pj)
      if (phase==0) then
!  possibly choose an active general c/s to relax (marked by rp>0)
        t=-1.D1*rp
        do 13 j=1,n
          i=abs(ls(j))
          if (i<=n) goto 13
          if (bu(i)==bl(i) .and. r(i)<0.D0) then
            r(i)=-r(i)
            ls(j)=-ls(j)
          end if
          if (r(i)/e(i)<=t) goto 13
          rp=r(i)
          t=rp/e(i)
          pj=j
13      continue
      end if

      if (ig==0) then
        gg=0.D0
        do j=n-k+1,n
          i=ls(j)
          gg=gg+r(i)**2
        end do
        rgnorm=sqrt(gg)
      end if
!     print 2,'rgtol,rgnorm,rp',rgtol,rgnorm,rp

25    continue
      if (rgnorm<=rgtol .and. abs(rp)<=gtol) then
!  allow for changes to norm(g)
        gnorm=sqrt(scpr(0.D0,g,g,n))
        gtol=sgnf*gnorm
        rgtol=max(rgt0l,gtol)
      end if

      if ((rgnorm<=rgtol .and. abs(rp)<=gtol) .or. ngr>mxgr) then
!  optimal at current level: first tidy up x
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
        if (ngr>mxgr) then
          f=fbase+f
          ifail=5
          return
        end if
        if (iprint>=2) then
          write(nout,*)'OPTIMAL at level 1'
          if (iprint>=3) then
!           write(nout,1000)'x variables',(x(i),i=1,n)
            write(nout,1001)'residual vector and indices', &
              (ls(j),r(abs(ls(j))),j=n1,nm)
          end if
        end if
        irep=irep+1
        if (irep<=nrep .and. iter>mpiv) then
          if (iprint>=1) write(nout,1)'refinement step #',irep
          mode=4
           goto 8
        end if
        if (iprint>=2 .and. nrep>0) &
          write(nout,*)'total number of restarts =',ires
        if (ninf>0) then
          ifail=3
          return
        end if
        nv=nv0
        ifail=0
        f=fbase+f
        return
      end if

      if (rgnorm>=abs(rp)) then
!  ignore the multiplier of c/s p and set up or continue SD steps
        p=0
      else
        p=abs(ls(pj))
        if (iprint>=2) print 1,'CHOOSE p =',p
        rp=r(p)
        call iexch(ls(pj),ls(n-k))
        pj=n-k
        ig=0
      end if

      if (p>0) then

!  compute +/- Steepest Edge (SE) search direction s in an(.)
        call tfbsub(n,a,la,p,ws(na1),ws(na1),ws(lu1),lws(ll1), &
          e(p),.true.)
        rp=scpr(0.D0,ws(na1),g,n)
        if (ls(pj)<0)rp=-rp
        if (rp*r(p)<=0.D0) then
          r(p)=0.D0
           goto 21
        end if
        if (abs(rp-r(p))>5.D-1*max(abs(rp),abs(r(p)))) then
!       if (abs(rp-r(p))>1.D-1*gnorm) then
          print 2,'1rp,r(p),rp-r(p)',rp,r(p),rp-r(p)
           goto 98
        end if
        snorm=e(p)
        plus=ls(pj)>=0.eqv.rp<0.D0
        f0=f
        ig=0
      else
        if (ig==0) then
!  start up the limited memory sweep method
!         if (p>0) then
!  transfer c/s p into Z
!           if (ls(pj)<0) then
!             r(p)=-r(p)
!             ls(pj)=-ls(pj)
!           end if
!           k=k+1
!           gg=gg+r(p)**2
!         end if
          ig=1
          ngv=1
          f0=f
          ws(kb1)=gg
          rgnorm=sqrt(gg)
!         print 2,'initial rg =',(r(ls(j)),j=n-k+1,n)
          if (k*ngv>kmax*maxg) then
            f=fbase+f
            ifail=9
            return
          end if
          call store_rg(k,ig,ws(krg1),r,ls(n-k+1))
        end if
!  compute Steepest Descent (SD) search direction s = -Z.rg in an(.)
        call zprod(k,n,a,la,ws(na1),r,w,ls,ws(lu1),lws(ll1))
        rp=scpr(0.D0,ws(na1),g,n)
        if (abs(gg+rp)>5.D-1*max(gg,abs(rp))) then
!       if (abs(gg+rp)>1.D-2*max(gg,abs(rp))) then
          print 2,'gg,rp,gg+rp',gg,rp,gg+rp
           goto 98
        end if
        snorm=sqrt(scpr(0.D0,ws(na1),ws(na1),n))
        plus=.true.
      end if
!     print 4,'s (or -s if .not.plus) =',(ws(i),i=na1,na+n)
!     print *,'plus =',plus

!  form At.s and denominators
      call form_Ats(n1,lp(1),n,plus,a,la,ws(na1),w,ls,snorm*tol)

!  return from degeneracy code
30    continue
      if (iprint>=3) then
        write(nout,1000)'x variables',(x(i),i=1,n)
        write(nout,1001)'residual vector and indices', &
          (ls(j),r(abs(ls(j))),j=n1,lp(1))
        write(nout,1000)'denominators',(w(abs(ls(j))),j=n1,lp(1))
      end if
!     read *,i

40    continue
!  level 1 ratio tests
      amax=ainfty
      qj=0
      qj1=0
      do 41 j=n-k+1,n
        i=ls(j)
        if (i<=0) print *,'i<=0'
        if (i<=0) stop
        si=ws(na+i)
        if (si==0.D0) goto 41
        t=abs(si)
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
      if (ig==0 .and. rp<0.D0 .and. bu(p)-bl(p)<amax) then
        amax=bu(p)-bl(p)
        qj=pj
      end if
      if (ninf>0) then
        alpha1=ainfty
        do 42 j=n1,lp(1)
          i=abs(ls(j))
          wi=w(i)
          if (wi==0.D0) goto 42
          ri=r(i)
          if (wi>0.D0) then
            if (ri<0.D0) goto 42
            z=(ri+tol)/wi
          else
            if (ri<0.D0) then
              z=ri/wi
              if (z<alpha1) then
                alpha1=z
                qj1=j
              end if
            end if
            z=((bl(i)-bu(i))+ri-tol)/wi
          end if
          if (z>=amax) goto 42
          amax=z
          qj=j
42      continue
        if (qj1>0 .and. alpha1<=amax) then
!  find feasible step that zeros most infeasible c/s
          do 43 j=n1,lp(1)
            i=abs(ls(j))
            wi=w(i)
            if (wi>=0.D0) goto 43
            ri=r(i)
            if (ri<0.D0) then
              z=ri/wi
              if (z>alpha1 .and. z<=amax) then
                alpha1=z
                qj1=j
              end if
            end if
43        continue
          amax=alpha1
          qj=qj1
        else
          qj1=0
        end if
      else
        do 44 j=n1,lp(1)
          i=abs(ls(j))
          wi=w(i)
          if (wi==0.D0) goto 44
          ri=r(i)
          if (wi>0.D0) then
            z=(ri+tol)/wi
          else
            z=(bl(i)-bu(i)+ri-tol)/wi
          end if
          if (z>=amax) goto 44
          amax=z
          qj=j
44      continue
      end if
      q=abs(ls(qj))
      if (iprint>=2 .and. q/=p .and. qj>n) &
        write(nout,*)'q,r(q),w(q) =',q,r(q),w(q)
      if (qj>n .and. qj1==0) then
        if (w(q)>0.D0) then
          amax=r(q)/w(q)
        else
          amax=(bl(q)-bu(q)+r(q))/w(q)
        end if
      end if

      if (amax==0.D0 .and. rp<=0.D0) then
        alpha=0.D0
!  potential degeneracy block at level 1
        if (p==0) goto 65
        if (bu(q)==bl(q)) goto 70
        plev=n
        do j=n1,lp(1)
          i=abs(ls(j))
          if (r(i)==0.D0) then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            if (bu(i)>bl(i))r(i)=1.D0
          end if
        end do
        if (plev>n1) then
          lp(2)=plev
          lev=2
          alp(1)=f
          f=0.D0
          qj=pj
          q=p
          if (iprint>=1) write(nout,'(''pivots ='',I5,''     level = 2'', ''    f ='',E16.8)')npv,f
           goto 86
        end if
        qj=n1
        r(q)=0.D0
!       print *,'only one degenerate c/s'
         goto 70
      end if

      if (ninf>0) then
        alpha=amax
      else
        if (linear) then
          alpha=amax
          ff=f+alpha*rp
          if (ff<fmin) goto 75
          f=ff
           goto 60
        end if
        call gdotx(n,ws(na1),ws,lws,ws(nb1))
        ngr=ngr+1
        sgs=scpr(0.D0,ws(na1),ws(nb1),n)
!       print 2,'rp,sgs',rp,sgs
        ggo=gg
        if (p==0) then
          t=v(nv)
          if (t<=0.D0) goto 52
          alpha=1.D0/t
          if (alpha>=amax) goto 52
          nv=nv-1
          ff=f+alpha*(rp+5.D-1*alpha*sgs)
          if (ff>=f0) goto 52
!           print 2,'alphar =',alpha
!  need to set f0 somewhere
          if (iprint>=2) write(nout,*)'Ritz value step:  alpha =', &
            alpha,'   p =',p
           goto 54
        end if
52      continue
        if (sgs>0.D0) then
          alpha=-rp/sgs
          if (alpha<amax) then
!     accept Cauchy step
            if (iprint>=2) write(nout,*)'Cauchy step:  alpha =', &
              alpha,'   p =',p
            ff=f+alpha*(rp+5.D-1*alpha*sgs)
            nv=0
!           print 2,'alphac =',alpha
             goto 54
          end if
        end if
!  Cauchy step infeasible
        alpha=amax
        ff=f+alpha*(rp+5.D-1*alpha*sgs)
        if (ff<fmin) goto 75
        if (ff>=f) then
          if (ires<nres) goto 98
          f=fbase+f
          if (iprint>=1) write(nout,'(''pivots ='',I5, ''  level = 1    f ='',E16.8)')npv,f
          ifail=4
          return
        end if
        f=ff
        if (plus) then
          call mysaxpy(alpha,ws(nb1),g,n)
        else
          call mysaxpy(-alpha,ws(nb1),g,n)
        end if
!       print 4,'new g =',(g(i),i=1,n)
        call newg
         goto 60
54    continue
        if (ff<fmin) goto 75
        if (ff>=f) then
          if (ires<nres) goto 98
          f=fbase+f
          if (iprint>=1) write(nout,'(''pivots ='',I5, ''  level = 1    f ='',E16.8)')npv,f
          ifail=4
          return
        end if
        f=ff
        if (plus) then
          call mysaxpy(alpha,ws(nb1),g,n)
        else
          call mysaxpy(-alpha,ws(nb1),g,n)
        end if
!       print 4,'new g =',(g(i),i=1,n)
        call newg
        if (ig==0) goto 60
        ig1=ig+1
        if (ig1>maxg)ig1=1
        call fbsub(n,1,n,a,la,0,g,w,ls,ws(lu1),lws(ll1),.true.)
!       print 4,'new rg =',(w(ls(j)),j=n-k+1,n)
        if (ngv<maxg)ngv=ngv+1
        if (k*ngv>kmax*maxg) then
          f=fbase+f
          ifail=9
          return
        end if
        call store_rg(k,ig1,ws(krg1),w,ls(n-k+1))
        gpg=0.D0
        gg=0.D0
        do j=n-k+1,n
          i=ls(j)
          gpg=gpg+r(i)*w(i)
          gg=gg+w(i)**2
        end do
        rgnorm=sqrt(gg)
!       print 2,'gpg,gg',gpg,gg
!       print 2,'f =',f
        call signst(n,r,w,ls)
        ws(ka+ig)=1.D0/alpha
        ws(kb+ig1)=gg
        ws(kc+ig)=gpg
        if (nv==0 .or. gg>ggo) then
!  compute new Ritz values
          if (ngv==2) then
            nv=1
            v(1)=1.D0/alpha
          else
            nv=min(ngv-1,k)
            if (nv<=0) print 1,'ngv,k,ig,nv =',ngv,k,ig,nv
            if (nv<=0) stop
!           print 1,'ngv,k,ig,nv =',ngv,k,ig,nv
!           print 4,'G =',(ws(krg+i),i=1,k*ngv)
!           print 4,'a =',(ws(ka+i),i=1,ngv)
!           print 4,'b =',(ws(kb+i),i=1,ngv+1)
!           print 4,'c =',(ws(kc+i),i=1,ngv)
            call formR(nv,k,ig,maxg,ws(ka1),ws(kb1),ws(kc1),ws(kd1), &
              ws(ke1),ws(krg1),ws(kr1))
!           call checkT(nv,maxg,ws(kr1),ws(ke1),ws(kd1))
            call formT(nv,maxg,ws(kr1),v,ws(ke1))
!           print 4,'T matrix',(v(i),i=1,nv)
!             if (nv>1) print 5,(ws(ke+i),i=1,nv-1)
            call trid(v(1),ws(ke1),nv)
!           print 4,'eigenvalues of T',(v(i),i=1,nv)
            call insort(nv,v)
!           print 4,'sorted eigenvalues of T',(v(i),i=1,nv)
          end if
          nv0=nv
          f0=f
        end if
        ig=ig1
      end if

60    continue
      if (alpha>0.D0) then
!  update x
        if (plus) then
          call mysaxpy(alpha,ws(na1),x,n)
        else
          call mysaxpy(-alpha,ws(na1),x,n)
        end if
!  update r for inactive c/s
        iter=iter+1
        if (ninf>0) then
          n_inf=0
          ff=f
          f=0.D0
          do 61 j=n1,lp(1)
            i=abs(ls(j))
            if (w(i)==0.D0) then
              if (r(i)>=0.D0) goto 61
              n_inf=n_inf+1
              f=f-r(i)
               goto 61
            end if
            ri=r(i)-alpha*w(i)
            if (abs(ri)<=tol)ri=0.D0
            if (r(i)<0.D0) then
              if (ri>=0.D0) then
!  remove contribution to gradient
                if (i>n) then
                  call saipy(sign(1.D0,dble(ls(j))),a,la,i-n,g,n)
                else
                  g(i)=0.D0
                end if
              else
                n_inf=n_inf+1
                f=f-ri
              end if
            end if
            if (w(i)<0.D0) then
              ro=(bu(i)-bl(i))-ri
              if (abs(ro)<=tol)ro=0.D0
              if (ro<ri) then
                ri=ro
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
61        continue
          if (n_inf/=ninf) then
            call iexch(ninf,n_inf)
            call newg
!         else if (f>=ff) then
          else if (f>=eps*ff+ff) then
             goto 98
          end if
        else
          n_inf=0
          do 62 j=n1,lp(1)
            i=abs(ls(j))
            if (w(i)==0.D0) goto 62
            ri=r(i)-alpha*w(i)
            if (w(i)<0.D0) then
              ro=(bu(i)-bl(i))-ri
              if (ro<ri) then
                ri=max(ro,0.D0)
                w(i)=-w(i)
                ls(j)=-ls(j)
              end if
            end if
            if (ri<=tol) then
              ri=0.D0
              if (i<=n) then
                if (ls(j)>=0) then
                  x(i)=bl(i)
                else
                  x(i)=bu(i)
                end if
              end if
            end if
            r(i)=ri
62        continue
        end if
      end if

      if (alpha<amax) then
        if (ig>0) then
!  continue limited memory SD iterations
          if (iprint>=1) write(nout,'(''pivots ='',I5, ''  level = 1    df ='',E16.8,''   rg ='',E12.4, ''  k ='',I4)')npv,f,rgnorm,k
          if (alpha>0.D0) goto 20
          print *,'alpha<=0'
           goto 98
        end if
!  Cauchy step with SE iteration
        k=k+1
        if (p<=n) then
          ls(pj)=p
           goto 15
        end if
!  case p>n: find best inactive simple bound to replace p in ls(pj)
        t=0.D0
        do j=n1,lp(1)
          i=abs(ls(j))
          if (i<=n) then
            ti=abs(ws(na+i))
            if (ti>t) then
              t=ti
              qj=j
            end if
          end if
        end do
        if (t<=snorm*tol) then
          print *,'no suitable simple bound available'
           goto 98
        end if
        q=abs(ls(qj))
        ls(qj)=q
        if (iprint>=2) write(nout,1)'New free variable',q
         goto 70
      end if

65    continue
      if (iprint>=2) &
        write(nout,*)'New active c/s:  alpha =',alpha,'   q =',q
      if (ig>0) then
!  case alpha=amax and SD step: find best free variable to relax
        k=k-1
        if (qj<=n) then
!  case: q is a free variable
          if (ws(na+q)>0.D0)ls(qj)=-q
          call iexch(ls(qj),ls(n-k))
          ig=0
          if (n_inf>0 .and. ninf==0) goto 10
           goto 15
        end if
        call fbsub(n,n-k,n,a,la,q,w,w,ls,ws(lu1),lws(ll1),.false.)
!       print 4,'w(n-k:n) =',(w(ls(j)),j=n-k,n)
        t=0.D0
        do j=n-k,n
          i=ls(j)
          ti=abs(w(i))/e(i)
          if (ti>t) then
            t=ti
            pj=j
          end if
        end do
        if (t<=tol) then
          print *,'no suitable free variable to relax'
           goto 98
        end if
        p=ls(pj)
        call iexch(ls(pj),ls(n-k))
        pj=n-k
        if (iprint>=2) write(nout,*)'relax free variable',p
      end if

!  return from degeneracy with an equality c/s
70    continue
      if (qj/=pj) then
!  pivot interchange
        if (iprint>=2) write(nout,*)'replace',p,' by',q
        if (p==0) print *,'p==0'
        if (p==0) goto 98
        call pivot(p,q,n,nmi,a,la,e,ws(lu1),lws(ll1),ifail,npv)
        if (ifail>=1) then
          if (ifail>=2) then
            ifail=11
            return
          end if
          if (iprint>=1) write(nout,*)'failure detected in pivot (1)'
          print *,'r(q),w(q),q',r(q),w(q),q
           goto 98
        end if
        if (rp>0.D0) then
          call iexch(ls(pj),ls(qj))
          call iexch(ls(lp(1)),ls(qj))
          lp(1)=lp(1)-1
          if (ninf>0) goto 15
           goto 9
        end if
        if (ig>0) then
          ri=x(p)-bl(p)
          ro=bu(p)-x(p)
          if (ro<ri) then
            ri=ro
            ls(pj)=-p
          end if
          if (ri<=tol)ri=0.D0
          r(p)=ri
          ig=0
        else
          rpu=max(bu(p)-bl(p)-alpha,0.D0)
          if (alpha<=rpu) then
            rpu=alpha
          else
            ls(pj)=-ls(pj)
          end if
          if (abs(rpu)<=tol)rpu=0.D0
          r(p)=rpu
        end if
!       print 2,'r(p)',r(p)
        call iexch(ls(pj),ls(qj))
        if (phase>0 .and. bu(q)==bl(q)) then
          peq=peq+1
          call iexch(ls(pj),ls(peq))
        end if
        if (ninf==0) then
          if (phase==0) goto 9
          if (phase==1) goto 10
        end if
         goto 15
      end if
!  opposite bound comes active
      if (ninf==0) then
        if (iprint>=1) write(nout,'(''pivots ='',I5, ''  level = 1    f ='',E16.8)')npv,fbase+f
      else if (phase==0) then
        if (iprint>=1) write(nout,'(''pivots ='',I5, ''  level = 1    f ='',E16.8,''   ninfb ='',I4)') &
          npv,f,ninf
      else
        if (iprint>=1) write(nout,'(''pivots ='',I5, ''  level = 1    f ='',E16.8,''   ninf ='',I4)') &
          npv,f,ninf
      end if
      ls(pj)=-ls(pj)
      if (ninf==0 .and. .not.linear) goto 16
      if (ninf>0 .and. ninf/=n_inf) goto 16
      r(p)=-rp
       goto 20

!  unbounded solution case
75    continue
      irep=irep+1
      if (irep<=nrep .and. iter>mpiv) then
        mode=4
        if (iprint>=1) write(nout,1) &
          'unbounded solution identified: refinement step #',irep
         goto 8
      end if
      ifail=1
!  tidy up x
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
      nv=nv0
      f=fbase+f
      return

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
        if (peq>0 .or. k>0) write(nout,1) &
          '# active equality c/s and free variables = ',peq,k
      end if

84    continue
      call optest(peq+1,n-k,r,e,ls,rp,pj)

      if (-rp<=gtol) then
        if (iprint>=2) write(nout,*)'return to level 1'
        lev=1
        f=alp(1)
        do j=n1,lp(2)
          r(abs(ls(j)))=0.D0
        end do
        lev=1
        if (rp==0.D0 .and. phase>0) goto 25
         goto 20
      end if
      call iexch(ls(pj),ls(n-k))
      pj=n-k
      plus=ls(pj)>=0
      p=abs(ls(pj))
      rp=r(p)
!  compute search direction s in an(.)
      call tfbsub(n,a,la,p,ws(na1),ws(na1),ws(lu1),lws(ll1), &
        e(p),.true.)

        rp=scpr(0.D0,ws(na1),g,n)
        if (ls(pj)<0)rp=-rp
        if (rp*r(p)<=0.D0) then
          r(p)=0.D0
           goto 84
        end if
        if (abs(rp-r(p))>5.D-1*max(abs(rp),abs(r(p)))) then
!       if (abs(rp-r(p))>1.D-1*gnorm) then
          print 2,'2rp,r(p),rp-r(p)',rp,r(p),rp-r(p)
           goto 98
        end if

      snorm=e(p)
!  form At.s and denominators
      call form_Ats(n1,lp(lev),n,plus,a,la,ws(na1),w,ls,snorm*tol)
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
        if (wi<=0.D0) goto 90
        if (r(i)<0.D0) goto 90
        z=(r(i)+tol)/wi
        if (z>=alpha) goto 90
        alpha=z
        qj=j
90    continue
      if (qj==0) then
        do j=n1,lp(lev)
          i=abs(ls(j))
          w(i)=min(w(i),0.D0)
          r(i)=0.D0
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
          if (peq>0 .or. k>0) print 1, &
            '# active equality c/s and free variables = ',peq,k
        end if
!       call checkq(n,lp(1),nmi,kmax,g,a,la,x,bl,bu,r,ls,ws(nb1),
!         f,ws,lws,ninf,peq,k,1,p,rp,linear)
         goto 30
      end if
      q=abs(ls(qj))
      alpha=r(q)/w(q)
      ff=f+alpha*rp
      if (iprint>=2) then
        write(nout,*)'alpha =',alpha,'   p =',p,'   q =',q
        write(nout,2)'r(p),r(q),w(q) =',r(p),r(q),w(q)
      end if
!  test for equality c/s
      if (bu(q)==bl(q)) then
        do j=n1,lp(2)
          r(abs(ls(j)))=0.D0
        end do
        lev=1
        f=alp(1)
        alpha=0.D0
        if (iprint>=2) write(nout,*)'EQTY:   p =',p,'   q =',q, &
          '   return to level 1'
         goto 70
      end if
      if (alpha==0.D0) then
!  potential degeneracy block at level lev
        if (lev+2>mlp) then
          ifail=8
          return
        end if
        r(q)=0.D0
        plev=n
        do j=n1,lp(lev)
          i=abs(ls(j))
          if (r(i)==0.D0) then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            if (bu(i)>bl(i))r(i)=1.D0
          end if
        end do
        if (plev>n1) then
          lev=lev+1
          lp(lev)=plev
          alp(lev)=f
          f=0.D0
          if (iprint>=2) write(nout,*) &
            'degeneracy: increase level to ',lev
          if (iprint>=1) write(nout,'(''pivots ='',I5,A,''level ='',I2, ''    f ='',E16.8)')npv,spaces(:3*lev-1),lev,f
           goto 86
        end if
        qj=n1
      end if
      iter=iter+1
      if (iprint>=2) write(nout,*)'replace',p,' by',q
      call pivot(p,q,n,nmi,a,la,e,ws(lu1),lws(ll1),ifail,npv)
      if (ifail>=1) then
        if (ifail>=2) then
          ifail=11
          return
        end if
!       call iexch(ls(pj),ls(qj))
        if (iprint>=1) write(nout,*)'failure detected in pivot (4)'
!       print *,'r(q),w(q),q',r(q),w(q),q
         goto 98
      end if
!  update r and f
      do j=n1,lp(lev)
        i=abs(ls(j))
        ri=r(i)-alpha*w(i)
        if (abs(ri)<=tol)ri=0.D0
        r(i)=ri
      end do
      f=ff
!  exchange a constraint
      r(p)=alpha
      if (r(p)<=tol)r(p)=0.D0
      call iexch(ls(pj),ls(qj))
      if (iprint>=1) write(nout,'(''pivots ='',I5,A,''level ='',I2, ''    f ='',E16.8)')npv,spaces(:3*lev-1),lev,f
       goto 80
!  restart sequence
98    continue
      do i=1,n
        x(i)=min(bu(i),max(bl(i),x(i)))
      end do
      nk=peq
      do j=peq+1,n-k
        i=abs(ls(j))
        if (i>n) then
          nk=nk+1
          ls(nk)=ls(j)
        end if
      end do
      k=n-nk
      mode=2
      ires=ires+1
      if (iprint>=1) write(nout,*)'major restart #',ires
      tol=1.D1*tol
      if (ires<=nres) goto 7
      ifail=10
      return
1000  format(a/(e16.5,4e16.5))
1001  format(a/(i4,1x,e12.5,4(i4,1x,e12.5)))
!1000 format(a/(e18.8,3e19.8))
!1001 format(a/(i3,1x,e14.8,3(i4,1x,e14.8)))
      end

      subroutine stmapq(n,nm,kmax,maxg)
!  set storage map for workspace in qlcpd and auxiliary routines
      implicit double precision (a-h,r-z), integer (i-q)
      common/wsc/kk,ll,kkk,lll,mxws,mxlws
      common/lcpdc/na,na1,nb,nb1,krg,krg1,kr,kr1, &
        ka,ka1,kb,kb1,kc,kc1,kd,kd1,ke,ke1,lu1,ll1
!  double precision storage (ws)
!  locations 1:kk are user workspace for gdotx
!  scratch slots of length n+m and n
      na=kk
      na1=kk+1
      nb=na+nm
      nb1=nb+1
!  workspace of length kmax*maxg for reduced gradient vectors
      krg=nb+n
      krg1=krg+1
!  a slot of length maxg*(maxg+1)/2 and 5 slots of length maxg for sweep method
      kr=krg+kmax*maxg
      kr1=kr+1
      ka=kr+maxg*(maxg+1)/2
      ka1=ka+1
      kb=ka+maxg
      kb1=kb+1
      kc=kb+maxg
      kc1=kc+1
      kd=kc+maxg
      kd1=kd+1
      ke=kd+maxg
      ke1=ke+1
!  remaining space for use by denseL.f or schurQR.f
      lu1=ke1+maxg
!  total number of double precision locations required by qlcpd
      kkk=nm+n+maxg*(maxg+1)/2+maxg*(kmax+5)
!  integer storage (lws)
!  locations 1:ll are user workspace for gdotx
!  number of integer locations required by qlcpd
      lll=0
!  remaining space for use by denseL.f or schurQR.f
      ll1=ll+1
      return
      end

      subroutine setfg2(n,linear,a,la,x,f,g,ws,lws)
      implicit double precision (a-h,o-z)
      logical linear
      dimension a(*),la(*),x(*),g(*),ws(*),lws(*)
      common/noutc/nout
      if (linear) then
        do i=1,n
          g(i)=0.D0
        end do
        call saipy(1.D0,a,la,0,g,n)
        f=scpr(0.D0,x,g,n)
      else
        call gdotx(n,x,ws,lws,g)
        call saipy(1.D0,a,la,0,g,n)
        f=5.D-1*scpr(aiscpr(n,a,la,0,x,0.D0),g,x,n)
      end if
1     format(A,15I4)
2     format(A,5E15.7)
3     format(A/(20I4))
4     format(A/(5E15.7))
      return
      end

      subroutine checkq(n,nm,nmi,kmax,g,a,la,x,bl,bu,r,ls,an,f, &
        ws,lws,ninf,peq,k,lev,p,alp2,linear)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension g(*),a(*),la(*),x(*),bl(*),bu(*),r(*),ls(*), &
        an(*),ws(*),lws(*)
      common/noutc/nout
      common/epsc/eps,tol,emin
      logical linear
!     if (lev==2) then
!       do i=1,n
!         an(i)=g(i)
!       end do
!       e=alp2*sign(1.D0,dble(p))
!       i=abs(p)
!       if (i<=n) then
!         an(i)=an(i)-e
!       else
!         call saipy(-e,a,la,i-n,an,n)
!       end if
!        goto 10
!     end if
      j=nmi*(nmi+1)/2
      do i=1,nmi
        j=j-abs(ls(i))
      end do
      if (j/=0) write(nout,*)'indexing error'
      if (j/=0) stop
      do j=1,peq
        i=abs(ls(j))
        if (bu(i)>bl(i)) then
          write(nout,*)'non-equality constraint i =',i
          write(nout,*)'j,peq =',j,peq
          stop
        end if
      end do
      do j=n-k+1,n
        i=ls(j)
        if (i<=0 .or. i>n) then
          write(nout,*)'faulty free variable: i, j =',i,j
          stop
        end if
      end do
      e=0.D0
      do j=n+1,nm
        i=abs(ls(j))
        if (i<=n) then
          s=x(i)
        else
          s=aiscpr(n,a,la,i-n,x,0.D0)
        end if
        if (ls(j)>0) then
!         print *,'i,s,r(i),bl(i)',i,s,r(i),bl(i)
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
      if (e>tol) write(nout,*)'residual error at level 1 = ',e,ie
!     if (e>tol) stop
      if (ninf==0) then
        call setfg2(n,linear,a,la,x,ff,an,ws,lws)
      else
        do i=1,n
          an(i)=0.D0
        end do
        ff=0.D0
        do j=n+1,nm
          i=abs(ls(j))
          if (r(i)<0.D0) then
            ff=ff-r(i)
            if (i>n) then
              call saipy(-sign(1.D0,dble(ls(j))),a,la,i-n,an,n)
            else
              an(i)=an(i)-sign(1.D0,dble(ls(j)))
            end if
          end if
        end do
      end if
      gnm=sqrt(scpr(0.D0,an,an,n))
      if (lev==1 .and. max(abs(f),abs(ff))<1.D20) then
        e=abs(ff-f)
        if (e>tol*max(1.D0,abs(f))) write(nout,*)'function error = ',e, &
          '   f(x) =',ff
!     if (e>tol) stop
        if (e>tol*max(1.D0,abs(f))) print 4,'x =',(x(j),j=1,n)
        if (e>tol*max(1.D0,abs(f))) stop
      end if
10    continue
      e=0.D0
      do j=1,n
!       write(nout,*)'an =',(an(i),i=1,n)
        i=abs(ls(j))
        s=sign(1.D0,dble(ls(j)))
        if (i<=n) then
!         print *,'i,s,r(i)',i,s,r(i)
          an(i)=an(i)-s*r(i)
          if (j>n-k) then
            s=max(0.D0,bl(i)-x(i),x(i)-bu(i))
          else if (ls(j)>0) then
            s=x(i)-bl(i)
          else
            s=bu(i)-x(i)
          end if
        else
!         print *,'i,s,r(i)',i,s,r(i)
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
      if (e>tol) write(nout,*)'residual error at level 2 = ',e,ie
!     if (e>tol) stop
!     if (e>1.D-6) print 4,'x =',(x(i),i=1,n)
      if (e>1.D-6) stop
      e=0.D0
      do j=1,n
        if (abs(an(j))>e) then
          e=abs(an(j))
          ie=ls(j)
          je=j
        end if
      end do
!     write(nout,*)'KT condition error = ',e,je,ie,gnm
      if (e>gnm*tol) write(nout,*)'KT condition error = ',e,je,ie,gnm
!     if (e>gnm*tol) write(nout,4)'KT cond_n errors = ',(an(i),i=1,n)
      if (e>gnm*tol) stop
!     if (e>1.D-4) stop
1     format(A,15I4)
2     format(A,5E15.7)
3     format(A/(20I4))
4     format(A/(5E15.7))
      return
      end
