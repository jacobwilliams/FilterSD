!Christen this file glcpd.f
!ut here >>>>>>>>>>>>>>>>>

!  Copyright (C) 2010 Roger Fletcher

!  Current version dated 27 March 2013

!  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
!  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
!  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

      subroutine glcpd(n,m,k,kmax,maxg,a,la,x,bl,bu,f,fmin,g,r,w,e,ls, &
       alp,lp,mlp,peq,ws,lws,cws,v,nv,rgtol,m0de,ifail,mxgr,iprint,nout)
      implicit double precision (a-h,r-z), integer (i-q)

!  This routine finds a KT point for the General LCP (Linearly Constrained
!  Problem)

!       minimize    f(x)

!       subject to  l <= [I : A]t.x <= u                  (t = transpose)

!  where f(x) is a given function of n variables x, to be determined.
!  Lower and upper bound constraints on the variables x and the linear
!  functions At.x may be supplied, where A is an n*m matrix.
!  A recursive form of an active set method is used, using Wolfe's method to
!  resolve degeneracy. A limited memory reduced gradient sweep method is used
!  for minimization in the null space, so usually the KT point is a local
!  minimizer. Matrix information is made available and processed by calls to
!  external subroutines. Details of these are given in an auxiliary file
!  named either 'denseL.f' or 'schurQR.f'. (schurQR.f is a more recent
!  replacement for the file sparseL.f)

!  parameter list  (variables in a line starting with C must be set on entry)
!  **************

!  n     number of variables
!  m     number of general constraints (columns of A)
!  k     dimension of the null space obtained by eliminating the active
!        constraints (only to be set if mode>=2). The number of constraints in
!        the active set is n-k
!  kmax  maximum value of k (kmax <= n)
!  maxg  max number of reduced gradient vectors stored in sweep method:
!        (1 < maxg <= kmax+1 when kmax>0), typically maxg = min(6,kmax+1)
!  a(*)  storage of reals associated with A. This storage may be provided
!        in either dense or sparse format. Refer to either denseA.f or sparseA.f
!        for information on how to set a(*) and la(*). The vector c referred to
!        in these files should be set to the zero vector.
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
!        any inequality constraints, and locations n-k+1:n store the indices of
!        any free variables (variables not on a bound, which are used to
!        parametrise the null space: ls(j) is always positive in this range)
!          If mode>=2, the first n-k elements of ls must be set on entry
!  alp(mlp) workspace associated with recursion
!  lp(mlp)  list of pointers to recursion information in ls
!  mlp   maximum number of levels of recursion allowed (mlp>2: typically
!        mlp=50 would usually be adequate but mlp=m is an upper bound)
!  peq   pointer to the end of equality constraint indices in ls
!  ws(*) real workspace for gdotx (see below), qlcpd and denseL.f (or schurQR.f)
!          Set the total number in mxws (see "Common" below).
!  lws(*) integer workspace for gdotx, qlcpd and denseL.f (or schurQR.f).
!          Set the total number in mxlws (see "Common" below).
!        The storage maps for ws and lws are set by the routine stmap below
!  cws(*) character workspace (if any) needed by funct
!  v(maxg) set nv estimates of the eigenvalues of the reduced Hessian of f(x)
!          (for example from a previous run of glcpd). Set nv=1 and v(1)=1.D0
!          in absence of other information. New values of v are left on exit
!  nv    Number of estimates in v
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
!          4 = as 3 but it is also assumed that columns of A are unchanged
!                so that factors of the basis matrix stored in ws and lws are
!                valid (changes in f(x) and the vectors l and u are allowed)
!        A local copy (mode) of m0de is made and may be changed by glcpd
!  ifail   outcome of the process
!              0 = solution obtained
!              1 = unbounded problem (f(x)<fmin has occurred: note grad is not
!                    evaluated in this case)
!              2 = bl(i) > bu(i) for some i
!              3 = infeasible problem detected in Phase 1
!              4 = line search cannot improve f (possibly increase rgtol)
!              5 = mxgr gradient calls exceeded (this test is only carried
!                    out at the start of each iteration)
!              6 = incorrect setting of m, n, kmax, maxg, mlp, m0de or tol
!              7 = not enough space in ws or lws
!              8 = not enough space in lp (increase mlp)
!              9 = dimension of reduced space too large (increase kmax)
!             10 = maximum number of unsuccessful restarts taken
!            >10= possible use by later sparse matrix codes
!  mxgr  maximum number of gradient calls
!  iprint  switch for diagnostic printing (0 = off, 1 = summary,
!                 2 = scalar information, 3 = verbose)
!  nout  channel number for output

!  Storage Allocation
!  ******************
!  User information about the lengths of ws and lws is supplied to glcpd in
!    common/wsc/kk,ll,kkk,lll,mxws,mxlws
!  kk and ll refer to the lengths of ws and lws needed by the user subroutines.
!  kkk and lll are the numbers of locations used by glcpd and are set by glcpd.
!  The rest of ws and lws is used by the files denseL.f or schurQR.f
!  mxws and mxlws must be set to the total lengths of ws and lws available: a
!  message will be given if more storage is needed.

!  User subroutines
!  ****************

!  The user must provide two subroutines as follows

!      subroutine funct(n,x,f,ws,lws,cws)
!      implicit double precision (a-h,o-z)
!      dimension x(*),ws(*),lws(*)
!      character cws(*)
!      ...
!      statements to compute f(x) from x
!      ...
!      return
!      end

!      subroutine grad(n,x,g,ws,lws,cws)
!      implicit double precision (a-h,o-z)
!      dimension x(*),ws(*),lws(*)
!      character cws(*)
!      ...
!      statements to compute grad.f(x) in g from x (the user
!      may assume that a call of grad immediately follows one
!      of funct with the same vector x.)
!      ...
!      return
!      end

!  The parameters ws, lws and cws in the above subroutines enables data to be
!  passed from the user's calling program to these subroutines

!  Tolerances, accuracy and diagnostics
!  ************************************
!  glcpd uses tolerance and accuracy information stored in
!     common/epsc/eps,tol,emin
!     common/repc/sgnf,nrep,npiv,nres
!     common/refactorc/mc,mxmc
!     common/infoc/rgnorm,vstep,iter,npv,nfn,ngr
!  eps must be set to the machine precision (unit round-off) and tol is a
!  tolerance such that numbers whose absolute value is less than tol are
!  truncated to zero. This tolerance strategy in the code assumes that the
!  problem is well-scaled. The parameter sgnf is used to measure the maximum
!  allowable relative error in gradient values. If at any stage the accuracy
!  requirement rgtol < sgnf*rgnorm then rgtol is increased to sgnf*rgnorm
!    The code allows one or more refinement steps after the
!  calculation has terminated, to improve the accuracy of the solution,
!  and a fixed number nrep of such repeats is allowed. However the code
!  terminates without further repeats if no more than npiv pivots are taken.
!    In case of any breakdown, the code is restarted in mode 0.
!  The maximum number of unsuccessful restarts allowed is set in nres.
!    The basis matrix may be refactorised on occasions, for example to prevent
!  build-up of round-off in the factors or (when using schurQR.f) to limit
!  the growth in the Schur complement. The maximum interval between
!  refactorizations (or size of Schur complement) is set in mxmc.
!    Default values are set in block data but can be reset by the user.
!    infoc returns information about the progress of the method: rgnorm is the
!  norm of the reduced gradient on exit, and vstep is the length of the vertical
!  step in the warm start process. iter is the total number of iterations taken,
!  npv is the number of pivots, nfn is the number of function evaluations, and
!  ngr is the number of gradient evaluations.

      parameter (ainfty=1.D100)
      dimension a(*),la(*),x(*),bl(*),bu(*),g(*),r(*),w(*),e(*),ls(*), &
        alp(*),lp(*),ws(*),lws(*),v(*)
      character cws(*)
      character*32 spaces
      common/lcpdc/na,na1,nb,nb1,krg,krg1,kr,kr1, &
        ka,ka1,kb,kb1,kc,kc1,kd,kd1,ke,ke1,lu1,ll1
      common/epsc/eps,t0l,emin
!     common/epsc/eps,tol,emin
      common/infoc/rgnorm,vstep,iter,npv,nfn,ngr
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
      if(m<0 .or. n.le.0 .or. mlp<2 .or. mode<0 .or. mode>4 .or.  &
        kmax<0 .or. (kmax>0 .and. maxg.le.1) .or. tol.le.0.D0) then
        ifail=6
        return
      end if
      rgt0l=rgtol
      n1=n+1
      nm=n+m
      nmi=nm
      nfn=0
      ngr=0
      nv0=nv
      if(iprint>=3) then
        write(nout,1000)'lower bounds',(bl(i),i=1,nm)
        write(nout,1000)'upper bounds',(bu(i),i=1,nm)
      end if
      irep=0
      ires=0
      do i=1,nm
        t=bu(i)-bl(i)
        if(t<-tol) then
          print *,'i,bl(i),bu(i)',i,bl(i),bu(i)
          ifail=2
          return
        else if(t.le.tol) then
          bl(i)=5.D-1*(bl(i)+bu(i))
          bu(i)=bl(i)
        end if
      end do
      vmax=0.D0
      do i=1,n
        x(i)=min(bu(i),max(bl(i),x(i)))
        vmax=max(vmax,bu(i)-bl(i))
      end do
      if(mode.le.2) then
        call stmap(n,nm,kmax,maxg)
        if(mode==0) then
          nk=0
        else if(mode==1) then
!  collect equality c/s
          nk=0
          do i=1,nm
            if(bu(i)==bl(i)) then
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
      if(mode.le.3) then
!  set up factors of basis matrix and permutation vectors
        ifail=mode
        call start_up(n,nm,nmi,a,la,nk,e,ls,ws(lu1),lws(ll1),mode,ifail)
        if(ifail>0) return
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
      if(mode>0) then
        call warm_start(n,nm,a,la,x,bl,bu,r,ls,ws(lu1), &
          lws(ll1),ws(na1),vstep)
!       print *,'vstep,vmax',vstep,vmax
        if(vstep>2.D0*vmax) then
          mpiv=0
          mode=0
          nk=0
          do i=1,n
            x(i)=min(bu(i),max(bl(i),x(i)))
          end do
           goto 7
        end if
        if(vstep>tol)mpiv=0
      end if
      k=0
!  collect free variables
      do j=n,1,-1
        i=abs(ls(j))
        if(i.le.n .and. x(i)>bl(i) .and. x(i)<bu(i)) then
          call iexch(ls(j),ls(n-k))
          k=k+1
        end if
      end do
      if(mode==0) then
        do j=1,n-k
          i=ls(j)
          if(x(i)==bu(i))ls(j)=-i
        end do
        lp(1)=n
         goto 9
      end if
      phase=0
!  move inactive general c/s to the end
      do j=nm,n1,-1
        i=abs(ls(j))
        if(i>n) then
          call iexch(ls(j),ls(lp(1)))
          lp(1)=lp(1)-1
        end if
      end do
      call residuals(n,n1,lp(1),a,la,x,bl,bu,r,ls,f,g,ninf)
      if(ninf>0) then
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
        if(bu(i)==bl(i)) then
          peq=peq+1
          call iexch(ls(j),ls(peq))
        end if
      end do
      call residuals(n,lp(1)+1,nm,a,la,x,bl,bu,r,ls,f,g,ninf)
      lp(1)=nm
      if(ninf>0) then
        gnorm=sqrt(scpr(0.D0,g,g,n))
        gtol=sgnf*gnorm
        rgtol=max(rgt0l,gtol)
         goto 15
      end if
10    continue
      phase=2
      if(iprint>=1) write(nout,*)'FEASIBILITY OBTAINED at level 1'
      n_inf=0
      call funct(n,x,f,ws,lws,cws)
      nfn=nfn+1
      if(f<fmin) goto 75
      call grad(n,x,g,ws,lws,cws)
      ngr=ngr+1
!     write(nout,4)'x =',(x(i),i=1,n)
!     write(nout,4)'g =',(g(i),i=1,n)
      call newg
      gnorm=sqrt(scpr(0.D0,g,g,n))
      gtol=sgnf*gnorm
      rgtol=max(rgt0l,gtol)
      alpha=1.D0
      ig=0
      if(iprint>=1) write(nout,'(''pivots ='',I5, &
        ''  level = 1    f ='',E16.8)')npv,f
       goto 16
!  start of major iteration
15    continue
      if(iprint>=1) then
        if(ninf==0) then
          if(k>0) then
!           write(nout,'(''pivots ='',I5,
!    *        ''  level = 1    f ='',E16.8,''   k ='',I4)')npv,f,k
            write(nout,'(''pivots ='',I5, &
              ''  level = 1    f ='',E16.8,''   rg ='',E12.4, &
              ''  k ='',I4)')npv,f,rgnorm,k
          else
            write(nout,'(''pivots ='',I5, &
              ''  level = 1    f ='',E16.8)')npv,f
          end if
        else if(phase==0) then
          write(nout,'(''pivots ='',I5,''  level = 1    f ='', &
            E16.8,''   ninfb ='',I4)')npv,f,ninf
        else
          write(nout,'(''pivots ='',I5,''  level = 1    f ='', &
            E16.8,''   ninf ='',I4)')npv,f,ninf
        end if
      end if
16    continue
!  calculate multipliers
!     print 4,'gradient =',(g(i),i=1,n)
      do i=1,nm
        w(i)=0.D0
      end do
      call fbsub(n,1,n,a,la,0,g,w,ls,ws(lu1),lws(ll1),.true.)
      call signst(n,r,w,ls)
!  opposite bound or reset multiplier loop
20    continue
      if(iprint>=3) then
        write(nout,1001)'costs vector and indices', &
          (ls(j),r(abs(ls(j))),j=1,n)
!       write(nout,1000)'steepest edge coefficients',
!    *    (e(abs(ls(j))),j=1,n)
        if(peq>0 .or. k>0) write(nout,1) &
          '# active equality c/s and free variables = ',peq,k
      end if
!     call check(n,lp(1),nmi,kmax,g,a,la,x,bl,bu,r,ls,ws(nb1),f,
!    *  ws,lws,cws,ninf,peq,k,1,p,rp)

21    continue
      call optest(peq+1,n-k,r,e,ls,rp,pj)
      if(phase==0) then
!  possibly choose an active general c/s to relax (marked by rp>0)
        t=-1.D1*rp
        do 13 j=1,n
          i=abs(ls(j))
          if(i.le.n) goto 13
          if(bu(i)==bl(i) .and. r(i)<0.D0) then
            r(i)=-r(i)
            ls(j)=-ls(j)
          end if
          if(r(i)/e(i).le.t) goto 13
          rp=r(i)
          t=rp/e(i)
          pj=j
13      continue
      end if

      if(ig==0) then
        gg=0.D0
        do j=n-k+1,n
          i=ls(j)
          gg=gg+r(i)**2
        end do
        rgnorm=sqrt(gg)
      end if
!     print 2,'rgtol,rgnorm,rp',rgtol,rgnorm,rp

25    continue
      if(rgnorm.le.rgtol .and. abs(rp).le.gtol) then
!  allow for changes to norm(g)
        gnorm=sqrt(scpr(0.D0,g,g,n))
        gtol=sgnf*gnorm
        rgtol=max(rgt0l,gtol)
      end if
      if(iprint==3) print 2,'gtol,rgtol,rgnorm,rp',gtol,rgtol,rgnorm,rp

      if((rgnorm.le.rgtol .and. abs(rp).le.gtol) .or. ngr>mxgr) then
!  optimal at current level: first tidy up x
        do j=peq+1,n-k
          i=abs(ls(j))
          if(i.le.n) then
            if(ls(j)>=0) then
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
          if(r(i).le.tol .and. i.le.n) then
            r(i)=0.D0
            if(ls(j)>=0) then
              x(i)=bl(i)
            else
              x(i)=bu(i)
            end if
          end if
        end do
        if(ngr>mxgr) then
          ifail=5
          return
        end if
        if(iprint>=2) then
          write(nout,*)'OPTIMAL at level 1'
          if(iprint>=3) then
!           write(nout,1000)'x variables',(x(i),i=1,n)
            write(nout,1001)'residual vector and indices', &
              (ls(j),r(abs(ls(j))),j=n1,nm)
          end if
        end if
        irep=irep+1
        if(irep.le.nrep .and. iter>mpiv) then
          if(iprint>=1) write(nout,*)'refinement step #',irep
          mode=4
           goto 8
        end if
        if(iprint>=2 .and. nrep>0) &
           write(nout,*)'total number of restarts =',ires
        if(ninf>0) then
          ifail=3
          return
        end if
        nv=nv0
        ifail=0
        return
      end if

      if(rgnorm>=abs(rp)) then
!  ignore the multiplier of c/s p and set up or continue SD steps
        p=0
      else
        p=abs(ls(pj))
        if(iprint>=2) print 1,'CHOOSE p =',p
        rp=r(p)
        call iexch(ls(pj),ls(n-k))
        pj=n-k
        ig=0
      end if

!     if(k==0 .or. p>n) then
      if(p>0) then
!  compute +/- Steepest Edge (SE) search direction s in an(.)
        call tfbsub(n,a,la,p,ws(na1),ws(na1),ws(lu1),lws(ll1), &
          e(p),.true.)
        rp=scpr(0.D0,ws(na1),g,n)
        if(ls(pj)<0)rp=-rp
        if(rp*r(p).le.0.D0) then
          r(p)=0.D0
           goto 21
        end if
        if(abs(rp-r(p))>5.D-1*max(abs(rp),abs(r(p)))) then
!       if(abs(rp-r(p))>1.D-1*gnorm) then
          print 2,'1rp,r(p),rp-r(p)',rp,r(p),rp-r(p)
           goto 98
        end if
        snorm=e(p)
        plus=ls(pj)>=0.eqv.rp<0.D0
        f0=f
        ig=0
      else
        if(ig==0) then
!  start up the limited memory sweep method
!         if(p>0) then
!  transfer c/s p into Z
!           if(ls(pj)<0) then
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
          if(k*ngv>kmax*maxg) then
            ifail=9
            return
          end if
          call store_rg(k,ig,ws(krg1),r,ls(n-k+1))
        end if
!  compute Steepest Descent (SD) search direction s = -Z.rg in an(.)
        call zprod(k,n,a,la,ws(na1),r,w,ls,ws(lu1),lws(ll1))
        rp=scpr(0.D0,ws(na1),g,n)
        if(abs(gg+rp)>5.D-1*max(gg,abs(rp))) then
!       if(abs(gg+rp)>1.D-2*max(gg,abs(rp))) then
          print 2,'gg,rp,gg+rp',gg,rp,gg+rp
           goto 98
        end if
        snorm=sqrt(scpr(0.D0,ws(na1),ws(na1),n))
        plus=.true.
      end if
!     print 4,'s (or -s if .not.plus) =',(ws(i),i=na1,na+n)

!  form At.s and denominators
      call form_Ats(n1,lp(1),n,plus,a,la,ws(na1),w,ls,snorm*tol)

!  return from degeneracy code
30    continue

      if(iprint>=3) then
        write(nout,1000)'x variables',(x(i),i=1,n)
        write(nout,1001)'residual vector and indices', &
          (ls(j),r(abs(ls(j))),j=n1,lp(1))
        write(nout,1000)'denominators',(w(abs(ls(j))),j=n1,lp(1))
      end if

40    continue
!  level 1 ratio tests
      amax=ainfty
      qj=0
      qj1=0
      do 41 j=n-k+1,n
        i=ls(j)
        if(i.le.0) print *,'i.le.0'
        if(i.le.0) goto 98
        si=ws(na+i)
        if(si==0.D0) goto 41
        t=abs(si)
!       if(t.le.tol) goto 41
        if(si>0.D0.eqv.plus) then
          z=bu(i)-x(i)
          if(abs(z)<tol) then
            z=0.D0
            x(i)=bu(i)
          else
            z=z/t
          end if
        else
          z=x(i)-bl(i)
          if(abs(z)<tol) then
            z=0.D0
            x(i)=bl(i)
          else
            z=z/t
          end if
        end if
        if(z>amax) goto 41
        amax=z
        qj=j
41    continue
      if(ig==0 .and. rp<0.D0 .and. bu(p)-bl(p)<amax) then
        amax=bu(p)-bl(p)
        qj=pj
      end if
      if(ninf>0) then
        alpha1=ainfty
        do 42 j=n1,lp(1)
          i=abs(ls(j))
          wi=w(i)
          if(wi==0.D0) goto 42
          ri=r(i)
          if(wi>0.D0) then
            if(ri<0.D0) goto 42
            z=(ri+tol)/wi
          else
            if(ri<0.D0) then
              z=ri/wi
              if(z<alpha1) then
                alpha1=z
                qj1=j
              end if
            end if
            z=(bl(i)-bu(i)+ri-tol)/wi
          end if
          if(z>=amax) goto 42
          amax=z
          qj=j
42      continue
        if(qj1>0 .and. alpha1.le.amax) then
!  find feasible step that zeros most infeasible c/s
          do 43 j=n1,lp(1)
            i=abs(ls(j))
            wi=w(i)
            if(wi>=0.D0) goto 43
            ri=r(i)
            if(ri<0.D0) then
              z=ri/wi
              if(z>alpha1 .and. z.le.amax) then
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
          if(wi==0.D0) goto 44
          ri=r(i)
          if(wi>0.D0) then
            z=(ri+tol)/wi
          else
            z=(bl(i)-bu(i)+ri-tol)/wi
          end if
          if(z>=amax) goto 44
          amax=z
          qj=j
44      continue
      end if
      q=abs(ls(qj))
      if(iprint>=2 .and. q/=p .and. qj>n) &
        write(nout,*)'q,r(q),w(q) =',q,r(q),w(q)
      if(qj>n .and. qj1==0) then
        if(w(q)>0.D0) then
          amax=r(q)/w(q)
        else
          amax=(bl(q)-bu(q)+r(q))/w(q)
        end if
      end if

      if(amax==0.D0 .and. rp.le.0.D0) then
        alpha=0.D0
!  potential degeneracy block at level 1
        if(p==0) goto 65
        if(bu(q)==bl(q)) goto 70
        plev=n
        do j=n1,lp(1)
          i=abs(ls(j))
          if(r(i)==0.D0) then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            if(bu(i)>bl(i))r(i)=1.D0
          end if
        end do
        if(plev>n1) then
          lp(2)=plev
          lev=2
          alp(1)=f
          f=0.D0
          qj=pj
          q=p
          if(iprint>=1) write(nout,'(''pivots ='',I5,''     level = 2'', &
            ''    f ='',E16.8)')npv,f
           goto 86
        end if
        qj=n1
        r(q)=0.D0
!       print *,'only one degenerate c/s'
         goto 70
      end if

      if(ninf>0) then
        alpha=amax
        if(plus) then
          call mysaxpy(alpha,ws(na1),x,n)
        else
          call mysaxpy(-alpha,ws(na1),x,n)
        end if
      else
!  take a Ritz value off the stack
!       print 4,'Ritz values =',(v(i),i=1,nv)
        if(nv>0 .and. v(nv)>0.D0) then
          alpha=min(1.D0/v(nv),amax)
          nv=nv-1
        else
          alpha=amax
          nv=0
        end if
!  line search
        alphar=amax
        alphal=0.D0
        dalpha=alpha
        fi=f
        fr=ainfty
        ggo=gg
        gs=rp
        gsi=gs
!       print 2,'f0,fi,gsi,amax =',f0,fi,gsi,amax
51      continue
!  calculate new x
        if(plus) then
          call saxpyz(alpha,ws(na1),x,ws(nb1),n)
        else
          call saxpyz(-alpha,ws(na1),x,ws(nb1),n)
        end if
        call funct(n,ws(nb1),fp,ws,lws,cws)
        if(fp<fmin) goto 75
        nfn=nfn+1
        df=f-fp
!  check for lack of improvement
        if(fp>=f0) then
!         print 2,'alphal,alpha,fp =',alphal,alpha,fp
          if(dalpha<1.D-10 .and. df<-dalpha*gs) then
!           print *,'alpha too small'
            if(alphal>0.D0) goto 52
            ifail=4
            return
          end if
          fr=fp
          alphar=alpha
          z=5.D-1/(1.D0+df/(gs*dalpha))
!         print 2,'df,z =',df,z
          dalpha=dalpha*max(1.D-1,z)
          alpha=alphal+dalpha
          nv=0
           goto 51
        end if
        f=fp
        call grad(n,ws(nb1),g,ws,lws,cws)
        ngr=ngr+1
!       print 4,'new g =',(g(i),i=1,n)
        call newg
        gps=scpr(0.D0,g,ws(na1),n)
        if(.not.plus)gps=-gps
!       print 2,'fp,gps',fp,gps
!       print 2,'alphal,alpha,alphar',alphal,alpha,alphar
!  check for non-positive curvature
        if(alpha<amax .and. (gps.le.gsi .or. (gps<25.D-2*gsi .and.  &
          (alphal>0.D0 .or. fr<ainfty)))) then
!       if(alpha<amax .and. gps.le.gsi) then
          alphal=alpha
          if(fr==ainfty) then
            alpha=min(alpha*5.D0,amax)
            dalpha=alpha-alphal
          else
            dalpha=alphar-alpha
            z=max(2.D-1,5.D-1/(1.D0+(f-fr)/(gps*dalpha)))
            dalpha=dalpha*z
            alpha=min(alpha+dalpha,amax)
          end if
          gs=gps
          nv=0
           goto 51
        end if
!  end of line search
52      continue
        do i=1,n
          x(i)=ws(nb+i)
        end do
        if(ig==0) goto 60
        ig1=ig+1
        if(ig1>maxg)ig1=1
        call fbsub(n,1,n,a,la,0,g,w,ls,ws(lu1),lws(ll1),.true.)
!       print 4,'new rg =',(w(ls(j)),j=n-k+1,n)
        if(ngv<maxg)ngv=ngv+1
        if(k*ngv>kmax*maxg) then
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
        if(nv==0 .or. gg>ggo) then
!  compute new Ritz values
          if(ngv==0) then
            nv=1
            v(1)=1.D0/alpha
          else
            nv=min(ngv-1,k)
            if(nv.le.0) print 1,'ngv,k,ig,nv =',ngv,k,ig,nv
            if(nv.le.0) goto 98
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
!             if(nv>1) print 5,(ws(ke+i),i=1,nv-1)
            call trid(v,ws(ke1),nv)
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
      if(alpha>0.D0) then
!  update r for inactive c/s
        iter=iter+1
        if(ninf>0) then
          n_inf=0
          ff=f
          f=0.D0
          do 61 j=n1,lp(1)
            i=abs(ls(j))
            if(w(i)==0.D0) then
              if(r(i)>=0.D0) goto 61
              n_inf=n_inf+1
              f=f-r(i)
               goto 61
            end if
            ri=r(i)-alpha*w(i)
            if(abs(ri).le.tol)ri=0.D0
            if(r(i)<0.D0) then
              if(ri>=0.D0) then
!  remove contribution to gradient
                if(i>n) then
                  call saipy(sign(1.D0,dble(ls(j))),a,la,i-n,g,n)
                else
                  g(i)=0.D0
                end if
              else
                n_inf=n_inf+1
                f=f-ri
              end if
            end if
            if(w(i)<0.D0) then
              ro=(bu(i)-bl(i))-ri
              if(abs(ro).le.tol)ro=0.D0
              if(ro<ri) then
                ri=ro
                ls(j)=-ls(j)
              end if
            end if
            if(ri==0.D0 .and. i.le.n) then
              if(ls(j)>=0) then
                x(i)=bl(i)
              else
                x(i)=bu(i)
              end if
            end if
            r(i)=ri
61        continue
          if(n_inf/=ninf) then
            call iexch(ninf,n_inf)
            call newg
!         else if(f>=ff) then
          else if(f>=eps*ff+ff) then
             goto 98
          end if
        else
          n_inf=0
          do 62 j=n1,lp(1)
            i=abs(ls(j))
            if(w(i)==0.D0) goto 62
            ri=r(i)-alpha*w(i)
            if(w(i)<0.D0) then
              ro=(bu(i)-bl(i))-ri
              if(ro<ri) then
                ri=ro
                w(i)=-w(i)
                ls(j)=-ls(j)
              end if
            end if
            if(ri.le.tol) then
              ri=0.D0
              if(i.le.n) then
                if(ls(j)>=0) then
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

      if(alpha<amax) then
        if(ig>0) then
!  continue limited memory SD iterations
          if(iprint>=1) write(nout,'(''pivots ='',I5, &
            ''  level = 1    f ='',E16.8,''   rg ='',E12.4, &
            ''  k ='',I4)')npv,f,rgnorm,k
          if(alpha>0.D0) goto 20
          print *,'alpha.le.0'
           goto 98
        end if
!  Cauchy step with SE iteration
        k=k+1
        if(p.le.n) then
          ls(pj)=p
           goto 15
        end if
!  case p>n: find best inactive simple bound to replace p in ls(pj)
        t=0.D0
        do j=n1,lp(1)
          i=abs(ls(j))
          if(i.le.n) then
            ti=abs(ws(na+i))
            if(ti>t) then
              t=ti
              qj=j
            end if
          end if
        end do
        if(t.le.snorm*tol) then
          print *,'no suitable simple bound available'
           goto 98
        end if
        q=abs(ls(qj))
        ls(qj)=q
        if(iprint>=2) write(nout,1)'New free variable',q
         goto 70
      end if

65    continue
      if(iprint>=2) &
        write(nout,*)'New active c/s:  alpha =',alpha,'   q =',q
      if(ig>0) then
!  case alpha=amax and SD step: find best free variable to relax
        k=k-1
        if(qj.le.n) then
!  case: q is a free variable
          if(ws(na+q)>0.D0)ls(qj)=-q
          call iexch(ls(qj),ls(n-k))
          ig=0
          if(n_inf>0 .and. ninf==0) goto 10
           goto 15
        end if
        call fbsub(n,n-k,n,a,la,q,w,w,ls,ws(lu1),lws(ll1),.false.)
!       print 4,'w(n-k:n) =',(w(ls(j)),j=n-k,n)
        t=0.D0
        do j=n-k,n
          i=ls(j)
          ti=abs(w(i))/e(i)
          if(ti>t) then
            t=ti
            pj=j
          end if
        end do
        if(t.le.tol) then
          print *,'no suitable free variable to relax'
           goto 98
        end if
        p=ls(pj)
        call iexch(ls(pj),ls(n-k))
        pj=n-k
        if(iprint>=2) write(nout,*)'relax free variable',p
      end if


!  return from degeneracy with an equality c/s
70    continue
      if(qj/=pj) then
!  pivot interchange
        if(iprint>=2) write(nout,*)'replace',p,' by',q
        if(p==0) print *,'p==0'
        if(p==0) goto 98
        call pivot(p,q,n,nmi,a,la,e,ws(lu1),lws(ll1),ifail,npv)
        if(ifail>=1) then
!         if(ifail>=2) return
          if(ifail==7) return
          if(iprint>=1) write(nout,*)'failure detected in pivot (1)'
!         print *,'r(q),w(q),q',r(q),w(q),q
           goto 98
        end if
        if(rp>0.D0) then
          if(phase>0) print *,'phase =',phase
          call iexch(ls(pj),ls(qj))
          call iexch(ls(lp(1)),ls(qj))
          lp(1)=lp(1)-1
          if(ninf>0) goto 15
           goto 9
        end if
        if(ig>0) then
          ri=x(p)-bl(p)
          ro=bu(p)-x(p)
          if(ro<ri) then
            ri=ro
            ls(pj)=-p
          end if
          if(ri.le.tol)ri=0.D0
          r(p)=ri
          ig=0
        else
          rpu=max(bu(p)-bl(p)-alpha,0.D0)
          if(alpha.le.rpu) then
            rpu=alpha
          else
            ls(pj)=-ls(pj)
          end if
          if(abs(rpu).le.tol)rpu=0.D0
          r(p)=rpu
        end if
!       print 2,'r(p)',r(p)
        call iexch(ls(pj),ls(qj))
        if(phase>0 .and. bu(q)==bl(q)) then
          peq=peq+1
          call iexch(ls(pj),ls(peq))
        end if
        if(ninf==0) then
          if(phase==0) goto 9
          if(phase==1) goto 10
        end if
         goto 15
      end if
!  opposite bound comes active
      if(ninf==0) then
        if(iprint>=1) write(nout,'(''pivots ='',I5, &
          ''  level = 1    f ='',E16.8)')npv,f
      else if(phase==0) then
        if(iprint>=1) write(nout,'(''pivots ='',I5, &
          ''  level = 1    f ='',E16.8,''   ninfb ='',I4)') &
          npv,f,ninf
      else
        if(iprint>=1) write(nout,'(''pivots ='',I5, &
          ''  level = 1    f ='',E16.8,''   ninf ='',I4)') &
          npv,f,ninf
      end if
      ls(pj)=-ls(pj)
      if(ninf==0 .or. ninf/=n_inf) goto 16
      r(p)=-rp
       goto 20

!  unbounded solution case
75    continue
      irep=irep+1
      if(irep.le.nrep .and. iter>mpiv) then
        mode=4
        if(iprint>=1) write(nout,*) &
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
        if(r(i)==0.D0 .and. i.le.n) then
          if(ls(j)>=0) then
            x(i)=bl(i)
          else
            x(i)=bu(i)
          end if
        end if
      end do
      nv=nv0
      return

!  recursive code for resolving degeneracy (Wolfe's method)
80    continue
!  calculate multipliers
      call fbsub(n,1,n,a,la,0,g,w,ls,ws(lu1),lws(ll1),.true.)
      call signst(n,r,w,ls)
!  reset multiplier loop
82    continue
      if(iprint>=3) then
        write(nout,1001)'costs vector and indices', &
          (ls(j),r(abs(ls(j))),j=1,n)
!       write(nout,1000)'steepest edge coefficients',
!    *    (e(abs(ls(j))),j=1,n)
        if(peq>0 .or. k>0) write(nout,1) &
          '# active equality c/s and free variables = ',peq,k
      end if

84    continue
      call optest(peq+1,n-k,r,e,ls,rp,pj)

      if(-rp.le.gtol) then
        if(iprint>=2) write(nout,*)'return to level 1'
        lev=1
        f=alp(1)
        do j=n1,lp(2)
          r(abs(ls(j)))=0.D0
        end do
        lev=1
        if(rp==0.D0 .and. phase>0) goto 25
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
        if(ls(pj)<0)rp=-rp
        if(rp*r(p).le.0.D0) then
          r(p)=0.D0
           goto 84
        end if
        if(abs(rp-r(p))>5.D-1*max(abs(rp),abs(r(p)))) then
!       if(abs(rp-r(p))>1.D-1*gnorm) then
          print 2,'2rp,r(p),rp-r(p)',rp,r(p),rp-r(p)
           goto 98
        end if

      snorm=e(p)
!  form At.s and denominators
      call form_Ats(n1,lp(lev),n,plus,a,la,ws(na1),w,ls,snorm*tol)
86    continue
      if(iprint>=3) then
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
        if(wi.le.0.D0) goto 90
        if(r(i)<0.D0) goto 90
        z=(r(i)+tol)/wi
        if(z>=alpha) goto 90
        alpha=z
        qj=j
90    continue
      if(qj==0) then
        do j=n1,lp(lev)
          i=abs(ls(j))
          w(i)=min(w(i),0.D0)
          r(i)=0.D0
        end do
        call form_Ats(lp(lev)+1,lp(lev-1),n,plus,a,la,ws(na1), &
          w,ls,snorm*tol)
        lev=lev-1
        f=alp(lev)
        if(iprint>=2) write(nout,*)'UNBOUNDED:   p =',p, &
          '   return to level',lev
        if(lev>1) goto 86
        if(iprint>=3) then
          write(nout,1001)'costs vector and indices', &
            (ls(j),r(abs(ls(j))),j=1,n)
          if(peq>0 .or. k>0) print 1, &
            '# active equality c/s and free variables = ',peq,k
        end if
!       call check(n,lp(1),nmi,kmax,g,a,la,x,bl,bu,r,ls,ws(nb1),f,
!    *    ws,lws,cws,ninf,peq,k,1,p,rp)
         goto 30
      end if
      q=abs(ls(qj))
      alpha=r(q)/w(q)
      ff=f+alpha*rp
      if(iprint>=2) then
        write(nout,*)'alpha =',alpha,'   p =',p,'   q =',q
        write(nout,2)'r(p),r(q),w(q) =',r(p),r(q),w(q)
      end if
!  test for equality c/s
      if(bu(q)==bl(q)) then
        do j=n1,lp(2)
          r(abs(ls(j)))=0.D0
        end do
        lev=1
        f=alp(1)
        alpha=0.D0
        if(iprint>=2) write(nout,*)'EQTY:   p =',p,'   q =',q, &
          '   return to level 1'
         goto 70
      end if
      if(alpha==0.D0) then
!  potential degeneracy block at level lev
        if(lev+2>mlp) then
          ifail=8
          return
        end if
        r(q)=0.D0
        plev=n
        do j=n1,lp(lev)
          i=abs(ls(j))
          if(r(i)==0.D0) then
            plev=plev+1
            call iexch(ls(j),ls(plev))
            if(bu(i)>bl(i))r(i)=1.D0
          end if
        end do
        if(plev>n1) then
          lev=lev+1
          lp(lev)=plev
          alp(lev)=f
          f=0.D0
          if(iprint>=2) write(nout,*) &
            'degeneracy: increase level to ',lev
          if(iprint>=1) write(nout,'(''pivots ='',I5,A,''level ='',I2, &
            ''    f ='',E16.8)')npv,spaces(:3*lev-1),lev,f
           goto 86
        end if
        qj=n1
      end if
      iter=iter+1
      if(iprint>=2) write(nout,*)'replace',p,' by',q
      call pivot(p,q,n,nmi,a,la,e,ws(lu1),lws(ll1),ifail,npv)
      if(ifail>=1) then
!       if(ifail>=2) return
        if(ifail==7) return
!       call iexch(ls(pj),ls(qj))
        if(iprint>=1) write(nout,*)'failure detected in pivot (2)'
!       print *,'r(q),w(q),q',r(q),w(q),q
         goto 98
      end if
!  update r and f
      do j=n1,lp(lev)
        i=abs(ls(j))
        ri=r(i)-alpha*w(i)
        if(abs(ri).le.tol)ri=0.D0
        r(i)=ri
      end do
      f=ff
!  exchange a constraint
      r(p)=alpha
      if(r(p).le.tol)r(p)=0.D0
      call iexch(ls(pj),ls(qj))
      if(iprint>=1) write(nout,'(''pivots ='',I5,A,''level ='',I2, &
        ''    f ='',E16.8)')npv,spaces(:3*lev-1),lev,f
       goto 80
!  restart sequence
98    continue
      do i=1,n
        x(i)=min(bu(i),max(bl(i),x(i)))
      end do
      nk=peq
      do j=peq+1,n-k
        i=abs(ls(j))
        if(i>n) then
          nk=nk+1
          ls(nk)=ls(j)
        end if
      end do
      k=n-nk
      mode=2
      ires=ires+1
      if(iprint>=1) write(nout,*)'major restart #',ires
      tol=1.D1*tol
      if(ires.le.nres) goto 7
      ifail=10
      return
1000  format(a/(e16.5,4e16.5))
1001  format(a/(i4,1x,e12.5,4(i4,1x,e12.5)))
!1000 format(a/(e18.8,3e19.8))
!1001 format(a/(i3,1x,e14.8,3(i4,1x,e14.8)))
      end

!     block data defaults
!     implicit double precision (a-h,o-z)
!     common/epsc/eps,tol,emin
!     common/repc/sgnf,nrep,npiv,nres
!     common/refactorc/mc,mxmc
!     common/wsc/kk,ll,kkk,lll,mxws,mxlws
!     data  eps,    tol,   emin, sgnf, nrep, npiv, nres, mxmc, kk, ll
!    * /1111.D-19, 1.D-12, 0.D0, 1.D-8,  2,    3,   2,   500,   0,  0/
!     end

      subroutine stmap(n,nm,kmax,maxg)
!  set storage map for workspace in glcpd and auxiliary routines
      implicit double precision (a-h,r-z), integer (i-q)
      common/wsc/kk,ll,kkk,lll,mxws,mxlws
      common/lcpdc/na,na1,nb,nb1,krg,krg1,kr,kr1, &
        ka,ka1,kb,kb1,kc,kc1,kd,kd1,ke,ke1,lu1,ll1
!  double precision storage (ws)
!  locations 1:kk are user workspace for funct and grad
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
!  total number of double precision locations required by glcpd
      kkk=nm+n+maxg*(maxg+1)/2+maxg*(kmax+5)
!  integer storage (lws)
!  locations 1:ll are user workspace for funct and grad
!  number of integer locations required by glcpd
      lll=0
!  remaining space for use by denseL.f or schurQR.f
      ll1=ll+1
      return
      end

      subroutine check(n,nm,nmi,kmax,g,a,la,x,bl,bu,r,ls,an,f, &
        ws,lws,cws,ninf,peq,k,lev,p,alp2)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension g(*),a(*),la(*),x(*),bl(*),bu(*),r(*),ls(*), &
        an(*),ws(*),lws(*)
      character cws(*)
      common/noutc/nout
      common/epsc/eps,tol,emin
!     if(lev==2) then
!       do i=1,n
!         an(i)=g(i)
!       end do
!       e=alp2*sign(1.D0,dble(p))
!       i=abs(p)
!       if(i.le.n) then
!         an(i)=an(i)-e
!       else
!         call saipy(-e,a,la,i-n,an,n)
!       end if
!        goto 1
!     end if
      j=nmi*(nmi+1)/2
      do i=1,nmi
        j=j-abs(ls(i))
      end do
      if(j/=0) write(nout,*)'indexing error'
      if(j/=0) stop
      do j=1,peq
        i=abs(ls(j))
        if(bu(i)>bl(i)) then
          write(nout,*)'non-equality constraint i =',i
          stop
        end if
      end do
      do j=n-k+1,n
        i=ls(j)
        if(i.le.0 .or. i>n) then
          write(nout,*)'faulty free variable: i, j =',i,j
          stop
        end if
      end do
      e=0.D0
      do j=n+1,nm
        i=abs(ls(j))
        if(i.le.n) then
          s=x(i)
        else
          s=aiscpr(n,a,la,i-n,x,0.D0)
        end if
        if(ls(j)>0) then
!         print *,'i,s,r(i),bl(i)',i,s,r(i),bl(i)
          s=r(i)-s+bl(i)
        else
          s=r(i)+s-bu(i)
        end if
        if(abs(s).le.tol*max(1.D0,abs(r(i))))s=0.D0
        if(abs(s)>e) then
          e=abs(s)
          ie=i
        end if
      end do
      if(e>tol) write(nout,*)'residual error at level 1 = ',e,ie
!     if(e>tol) stop
      if(e>1.D-6) stop
      if(ninf==0) then
        call funct(n,x,ff,ws,lws,cws)
        call grad(n,x,an,ws,lws,cws)
      else
        do i=1,n
          an(i)=0.D0
        end do
        ff=0.D0
        do j=n+1,nm
          i=abs(ls(j))
          if(r(i)<0.D0) then
            ff=ff-r(i)
            if(i>n) then
              call saipy(-sign(1.D0,dble(ls(j))),a,la,i-n,an,n)
            else
              an(i)=an(i)-sign(1.D0,dble(ls(j)))
            end if
          end if
        end do
      end if
      gnm=sqrt(scpr(0.D0,an,an,n))
      if(lev==1 .and. max(abs(f),abs(ff))<1.D20) then
        e=abs(ff-f)
        if(e>tol*max(1.D0,abs(f))) write(nout,*)'function error = ',e, &
          '   f(x) =',ff
!     if(e>tol) stop
!       if(e>tol*max(1.D0,abs(f))) print 4,'x =',(x(j),j=1,n)
        if(e>tol*max(1.D0,abs(f))) stop
      end if
1     continue
      e=0.D0
      do j=1,n
!       write(nout,*)'an =',(an(i),i=1,n)
        i=abs(ls(j))
        s=sign(1.D0,dble(ls(j)))
        if(i.le.n) then
!         print *,'i,s,r(i)',i,s,r(i)
          an(i)=an(i)-s*r(i)
          if(j>n-k) then
            s=max(0.D0,bl(i)-x(i),x(i)-bu(i))
          else if(ls(j)>0) then
            s=x(i)-bl(i)
          else
            s=bu(i)-x(i)
          end if
        else
!         print *,'i,s,r(i)',i,s,r(i)
          call saipy(-s*r(i),a,la,i-n,an,n)
          if(ls(j)>0) then
            s=aiscpr(n,a,la,i-n,x,-bl(i))
          else
            s=-aiscpr(n,a,la,i-n,x,-bu(i))
          end if
        end if
        if(abs(s)>e) then
          e=abs(s)
          ie=i
        end if
      end do
      if(e>tol) write(nout,*)'residual error at level 2 = ',e,ie
!     if(e>tol) stop
!     if(e>1.D-6) print 4,'x =',(x(i),i=1,n)
      if(e>1.D-6) stop
      e=0.D0
      do j=1,n
        if(abs(an(j))>e) then
          e=abs(an(j))
          ie=ls(j)
          je=j
        end if
      end do
      if(e>gnm*tol) write(nout,*)'KT condition error = ',e,je,ie,gnm
!     if(e>gnm*tol) write(nout,4)'KT cond_n errors = ',(an(i),i=1,n)
!     if(e>gnm*tol) stop
      if(e>1.D-4) stop
2     format(A,5E15.7)
4     format(A/(5E15.6))
      return
      end
