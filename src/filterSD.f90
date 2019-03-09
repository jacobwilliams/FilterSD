
      subroutine filterSD(n,m,x,al,f,fmin,cstype,bl,bu,ws,lws,v,nv, &
        maxa,maxla,maxu,maxiu,kmax,maxg,rho,htol,rgtol,maxit,iprint, &
        nout,ifail)
      implicit double precision (a-h,o-z)
      dimension x(*),al(*),bl(*),bu(*),ws(*),lws(*),v(*)
      character cstype(*)

!  Copyright (C) 2010 Roger Fletcher

!  Current version dated 5 October 2011

!  THE ACCOMPANYING PROGRAM IS PROVIDED UNDER THE TERMS OF THE ECLIPSE PUBLIC
!  LICENSE ("AGREEMENT"). ANY USE, REPRODUCTION OR DISTRIBUTION OF THE PROGRAM
!  CONSTITUTES RECIPIENT'S ACCEPTANCE OF THIS AGREEMENT

!  Solves an NLP problem of the form: find a local solution x to
!
!          minimize    f(x)
!                            [  x   ]
!          subject to  bl <= [      ] <= bu
!                            [ c(x) ]
!
!  where f(x) is a given function of n variables x and c(x) is a vector of m
!  given constraint functions. f(x) is minimized subject to lower and upper
!  bounds bl and bu on x and c(x). f(x) and c(x) are defined by a user supplied
!  subroutine 'functions'. The user is also required to supply a subroutine
!  'gradients' which calculates gradients of f(x) and c(x) with respect to x.
!
!  If linearized constraints and the trust region are incompatible, the code
!  enters 'phase 1' in which an 'l1 feasibility problem' is solved. If this
!  is unsuccessful in resolving the situation then the code exits with
!  ifail=3 and returns a 'locally infeasible point' in x.
!
!  Parameter List  (variables in a line starting with C must be set on entry)
!  ==============
!  n     number of variables
!  m     number of general constraints
!  x(n+m)  x(1:n) stores the vector of variables. Initially an estimate of the
!            solution must be set, replaced by the solution (if it exists)
!            on exit. The rest of x is workspace
!  al(n+m) stores Lagrange multipliers at the solution on exit. A positive
!            multiplier indicates that the lower bound is active, and a
!            negative multiplier indicates that the upper bound is active.
!            Inactive constraints have a zero multiplier.
!  f     returns the value of f(x) when x is a feasible solution
!  fmin  set a strict lower bound on f(x) for feasible x (used to identify
!          an unbounded NLP)
!  cstype(m) character workspace: if ifail=3, cstype indicates constraints
!              that are infeasible in the L1 solution. cstype(i)='A' if the
!              lower bound on c/s i is infeasible, 'Z' if the upper bound is
!              infeasible, else 'N' if feasible.
!  bl(n+m) lower bounds on x and c(x) (use numbers no less than -ainfty
!            (see below) and where possible supply realistic bounds on x)
!  bu(n+m) upper bounds on x and c(x) (use numbers no greater than ainfty)
!  ws(*)   double precision workspace
!  lws(*)  integer workspace
!  v(maxg) stores nv Ritz values (estimates of eigenvalues of reduced Hessian)
!            supply the setting from a previous run of filterSD, or set
!            nv=1 and v(1)=1.D0 in absence of other information
!  nv      number of values set in v
!  maxa    maximum number of entries in the Jacobian a(*) set by gradients
!  maxla   number of entries required for sparse matrix indices and pointers
!            la(0:*) to be set up in lws(*) (maxla>=maxa+m+3).
!            Set maxla=1 if using dense matrix format
!  maxu    length of workspace user(*) passed through to user subroutines
!            'functions' and 'gradients'
!  maxiu   length of workspace iuser(*) passed through to user subroutines
!  kmax    maximum dimension of null space allowed for (kmax<=n)
!  maxg    maximum number of reduced gradient vectors stored by the
!            limited memory method (typically 6 or 7)
!  rho     initial trust region radius (typically 1.D1)
!  htol    tolerance allowed in sum h of constraint feasibilities (e.g. 1.D-6)
!  rgtol   tolerance allowed in reduced gradient l2 norm (typically 1.D-4)
!  maxit   maximum number of major iterations allowed
!  iprint  verbosity of printing (0=none, 1=one line per iteration,
!            2=additional text information given)
!  nout    output channel for printing
!  ifail   returns failure indication as follows
!              0 = successful run
!              1 = unbounded NLP (f <= fmin at an htol-feasible point)
!              2 = bounds on x are inconsistent
!              3 = local minimum of feasibility problem and h > htol
!                  (nonlinear constraints are locally inconsistent)
!              4 = initial point x has h > ubd (reset ubd or x and re-enter)
!              5 = maxit major iterations have been carried out
!              6 = termination with rho <= htol
!              7 = not enough workspace in ws or lws (see message)
!              8 = insufficient space for filter (increase mxf and re-enter)
!             >9 = unexpected fail in LCP solver (10 has been added to ifail)

!  User Routines
!  =============
!  The user must supply two subroutines to calculate f(x), c(x) and their
!  first derivatives as follows
!
!     subroutine functions(n,m,x,f,c,user,iuser)
!     implicit double precision (a-h,o-z)
!     dimension x(*),c(*),user(*),iuser(*)
!     ...
!     Statements to calculate f(x) and the m-vector c(x). The user is
!     responsible for ensuring that any failures such as IEEE errors
!     (overflow, NaN's etc.) are trapped and not returned to filterSD.
!     The same holds for gradients.
!     ...
!     return
!     end
!
!     subroutine gradients(n,m,x,a,user,iuser)
!     implicit double precision(a-h,o-z)
!     dimension x(*),a(*),user(*),iuser(*)
!     ...
!     Statements to calculate gradients of f(x) and c(x) and set in a(*).
!     The column vector grad(f) must be followed by the column vectors
!     grad(c_i), i=1,2,...,m, in the one dimensional array a(*). Either a
!     dense or sparse data structure may be used. If using the sparse data
!     structure, only stucturally non-zero entries are set. Pointers etc. for
!     the data structure are set once and for all in lws as described below.
!     The user may assume that a call of 'gradients' immediately follows one
!      of 'functions' with the same vector x.)
!     ...
!     return
!     end
!
!  The user must also supply a driver routine which calls filterSD. This must
!  set parameters and common blocks of filterSD as appropriate.
!  Space for x,al,bl,bu,ws,lws,v and cstype must be assigned.
!  If using the sparse data structure for setting gradients, indices and
!  pointers la(0:maxla-1) must be set in the driver in lws, immediately
!  following any user workspace in lws(1:maxiu). No changes in this
!  data structure are allowed during the operation of filterSD.
!  More details of the format of a(*) and la(*) are given in the file sparseA.f
!  For dense format just set maxla=1 and set lws(maxiu+1) to the
!  'stride' (>=n) used in setting the columns of grad(f) and grad(c_i).
!  For efficiency, constant entries in the gradients may be set in the driver.
!  However two copies of the gradients are kept by filterSD . These reside in
!  ws(maxu+1:maxu+maxa) and ws(maxu+maxa+1:maxu+2*maxa). Any constant entries
!  must be set in both copies.

!  Common blocks
!  =============
!      common/wsc/kk,ll,kkk,lll,mxws,mxlws
!  The user must specify the length of the workspace arrays ws(*) and lws(*)
!  in mxws and mxlws respectively. It may not be easy to specify a-priori how
!  large these arrays should be. Set a suitable large estimate, and filterSD
!  will prompt if larger values are required. As a guide, ws(*) contains
!  first user workspace, then workspace for filterSD, then workspace for glcpd,
!  and finally workspace for denseL.f or schurQR.f. lws(*) contains user
!  workspace, then maxla+n+m+mlp locations for filterSD and additional locationsc  for denseL.f or schurQR.f.
!     common/defaultc/ainfty,ubd,mlp,mxf
!  Default values of some control parameters are set here. ainfty is used to
!  represent infinity. ubd provides an upper bound on the allowed constraint
!  violation. mlp is the maximum length of arrays used in degeneracy control.
!  mxf is the maximum length of filter arrays. Default values are 1.D20, 1.D4,
!  50,  50 respectively.
!     common/ngrc/mxgr
!  The user can limit the time spent in each call of the LCP solver by setting
!  an upper limit on the number of gradient calls in mxgr (default=1000000)
!     common/mxm1c/mxm1
!  When using denseL.f, mxm1 must be set to the maximum number of general
!  constraints allowed in the active set. mxm1=min(m+1,n) is always sufficient
!     common/epsc/eps,tol,emin
!     common/repc/sgnf,nrep,npiv,nres
!  These common blocks provide default parameters that control glcpd
!     common/statsc/dnorm,h,hJt,hJ,ipeq,k,itn,nft,ngt
!  This common block returns information about the outcome of filterSD.
!  dnorm=final step length, h=final c/s violation, hJt=ditto for 'N' c/s,
!  hJ=ditto for 'A' and 'Z' c/s, ipeq=number of active equations, k=number of
!  free variables, itn=number of iterations, nft=total number of function
!  calls, ngt=total number of gradient calls

      common/wsc/kk,ll,kkk,lll,mxws,mxlws
      common/defaultc/ainfty,ubd,mlp,mxf
      common/functc/fxd,alc,m_,iph,last1,next1,nx,nx1, &
        nal,nal1,naal,naal1,nxd,nxd1,ncx,ncx1,ncxd,ncxd1,nla1

      m_=m
      nm=n+m
!  set real storage map for ws
!  first maxu locations are user storage for functions and gradients
!  vectors required by funct: two slots of length maxa for a(*)
      last1=maxu+1
      next1=last1+maxa
!  slot of length n for x
      nx1=next1+maxa
      nx=nx1-1
!  slot of length m for lambda
      nal=nx+n
      nal1=nal+1
!  slot of length n for Ak.al, or for storing x
      naal=nal+m
      naal1=naal+1
!  slot of length n for x at x+d
      nxd=naal+n
      nxd1=nxd+1
!  slot of length m for c at x
      ncx=nxd+n
      ncx1=ncx+1
!  slot of length m for c at x+d
      ncxd=ncx+m
      ncxd1=ncxd+1
!  local storage for filter_SD
!  slot of length n for d
      id1=ncxd1+m
!  slot of length n+m for dl
      idl1=id1+n
!  slot of length n+m for du
      idu1=idl1+nm
!  slot of length n for g
      ig1=idu1+nm
!  slot of length n+m for e
      ie1=ig1+n
!  slot of length mlp for alp
      ialp1=ie1+nm
!  slot of length mxf for filh
      ifilh1=ialp1+mlp
!  slot of length mxf for filf
      ifilf1=ifilh1+mxf
!  total length of ws so far
      kk=ifilf1+mxf-1

!  set integer storage map for lws
!  first maxiu locations are user storage for functions and gradients
!  storage of length maxla for la(0:*)
      nla1=maxiu+1
!  local storage for filter_SD
!  slot of length n+m for ls
      ils1=nla1+maxla
      ils=ils1-1
!  slot of length mlp for lp
      ilp1=ils1+nm
!  total length of lws so far
      ll=ilp1+mlp-1

      do i=1,n
        if(bl(i).gt.bu(i))then
          if(iprint.gt.1)write(nout,*)'simple bounds infeasible'
          ifail=2
          return
        endif
        ws(nx+i)=min(max(bl(i),x(i)),bu(i))
      enddo

!  note x and al are just used as workspace: the true values are those in ws

      call filter_SD(n,f,fmin,cstype,bl,bu,ws,lws,v,nv, &
        maxa,kmax,maxg, &
        ws(id1),ws(idl1),ws(idu1),ws(ig1),x,al,ws(ie1),lws(ils1), &
        ws(ialp1),lws(ilp1),ws(ifilh1),ws(ifilf1),rho,htol,rgtol, &
        maxit,iprint,nout,ifail)

      if(ifail.ge.7)return
!  scatter ws(nx.. and ws(nal.. and bound multipliers into x and al
      do i=1,n
        al(i)=0.D0
      enddo
      do i=1,m
        al(n+i)=ws(nal+i)
      enddo
      do j=1,n
        i=abs(lws(ils+j))
        if(i.le.n)then
          if(ws(nx+i).eq.bl(i))then
            al(i)=x(i)
          elseif(ws(nx+i).eq.bu(i))then
            al(i)=-x(i)
          endif
        endif
      enddo
      do i=1,n
        x(i)=ws(nx+i)
      enddo

      return
1     format(A,15I5)
      end

      subroutine filter_SD(n,f,fmin,cstype,bl,bu,ws,lws,v,nv, &
        maxa,kmax,maxg, &
        d,dl,du,g,r,w,e,ls,alp,lp,filh,filf,rho,htol,rgtol, &
        maxit,iprint,nout,ifail)
      implicit double precision (a-h,o-z)
      dimension bl(*),bu(*),ws(*),lws(*),v(*), &
        d(*),dl(*),du(*),g(*),r(*),w(*),e(*),ls(*),alp(*),lp(*), &
        filh(*),filf(*)
      character cstype(*)

      parameter (sigma=1.D-1,infty=100000000)

      common/defaultc/ainfty,ubd,mlp,mxf
      common/epsc/eps,tol,emin
      common/repc/sgnf,nrep,npiv,nres
      common/functc/fxd,alc,m,iph,last1,next1,nx,nx1, &
        nal,nal1,naal,naal1,nxd,nxd1,ncx,ncx1,ncxd,ncxd1,nla1
      common/infoc/rgnorm,vstep,iter,npv,nfn,ngr
      common/ngrc/mxgr
      common/statsc/dnorm,h,hJt,hJ,ipeq,k,itn,nft,ngt

1     format(A,15I5)
2     format(A,6E15.7)
3     format(A/(15I5))
!   4 format(A/(6E13.5))
4     format(A/(5E15.7))
5     format((6E15.7))
6     format(A,2E15.7,I2)
1000  format(I4,1X,E14.6,E16.8,'  < reset J          ',11X,E12.4)
1001  format(I4,1X,E14.6,E16.8,'  < LCP',E13.5,2E12.4)
1002  format(I4,1X,E14.6,E16.8,'  < project')
2000  format(I4,E14.6,E16.8,'  << feasible LP    ',13X,E12.4)
2001  format(I4,E14.6,E16.8,'  << LCP',E13.5,2E12.4)
2002  format(I4,E14.6,E16.8,'  << project')

      n1=n+1
      nm=n+m
      mode=0
      nrep=0
      iph=2
      nfil=0
      itn=0
      nft=0
      ngt=0

!  evaluate f,c and a
      call functions(n,m,ws(nx1),f,ws(ncx1),ws,lws)
      call gradients(n,m,ws(nx1),ws(last1),ws,lws)
!  evaluate h
      h=0.D0
      do i=1,m
        ci=ws(ncx+i)
        h=h+max(0.D0,bl(n+i)-ci,ci-bu(n+i))
      enddo
      if(h.gt.ubd)then
        if(iprint.gt.1)write(nout,2)'h.gt.ubd: h =',h
        ifail=4
        return
      endif
      if(iprint.ge.1) &
        write(nout,*)' itn     h/hJt          f/hJ      ', &
        '           rgnorm       dnorm        rho'
10    continue
      if(rho.lt.htol)then
        if(iprint.gt.1)write(nout,2)'rho less than htol: rho =',rho
        ifail=6
        return
      endif
!     print 4,'x =',(ws(nx+i),i=1,n)
!     print 4,'c =',(ws(ncx+i),i=1,m)
!  set up LP subproblem
      do i=1,n
        dl(i)=max(-rho,bl(i)-ws(nx+i))
        du(i)=min(rho,bu(i)-ws(nx+i))
        d(i)=0.D0
      enddo
      if(abs(iph).eq.1)then
        do j=1,n
          i=abs(ls(j))-n
          if(i.gt.0)then
            if(cstype(i).eq.'A')then
              ls(j)=n+i
            elseif(cstype(i).eq.'Z')then
              ls(j)=-n-i
            endif
          endif
        enddo
      endif
      do i=1,m
        ci=ws(ncx+i)
        dl(n+i)=bl(n+i)-ci
        du(n+i)=bu(n+i)-ci
        ws(nal+i)=0.D0
        cstype(i)='N'
      enddo
      iph_=iph
      iph=0
      k=0
!     print *,'solve LP subproblem',itn
      iii=0
!     if(itn.eq.14)iii=1
      call glcpd(n,m,k,kmax,maxg,ws(last1),lws(nla1),d,dl,du,phi, &
        -ainfty,g,r,w,e,ls,alp,lp,mlp,ipeq,ws,lws,cstype,v,nv,rgtol, &
        mode,ifail,infty,iii,0)
!     print 4,'d =',(d(i),i=1,n)
!     print 4,'r =',(r(i),i=1,nm)
!     print 3,'ls =',(ls(i),i=1,nm)
!     print 1,'ipeq,k,ifail',ipeq,k,ifail
      mode=2
      iph=iph_
      if(ifail.eq.0.or.ifail.eq.4)then
        if(iprint.ge.1)write(nout,2000)itn,h,f,rho
        if(abs(iph).eq.1)nfil=nfil1-1
        iph=2
         goto 50
      elseif(ifail.ne.3)then
!       print 1,'itn =',itn
        if(iprint.gt.1)write(nout,*)'unexpected fail in LP subproblem'
         goto 99
      endif
15    continue
      if(h.le.htol)then
        if(iprint.gt.1)write(nout,*)'htol-feasible but LP is infeasible'
        ifail=0
        return
      endif
!  infeasibility: enter feasibility restoration
!     print *,'LP is infeasible: solve l1 subproblem',itn
      iii=0
!     if(itn.eq.9)iii=2
      call l1sold(n,m,k,kmax,maxg,ws(last1),lws(nla1),d,dl,du,phi, &
        g,r,w,e,ls,alp,lp,mlp,ipeq,ws,lws,cstype,v,nv,rgtol,ifail,iii,0)
!     print 4,'d1 =',(d(i),i=1,n)
!     print 4,'r =',(r(i),i=1,nm)
!     print 3,'ls =',(ls(i),i=1,nm)
      if(ifail.ne.0)then
        if(print.gt.1)print *,'unexpected fail in l1 subproblem'
         goto 99
      endif
      if(abs(iph).eq.2)then
        call addfil(h,f,filh,filf,1,nfil,mxf,ifail)
        if(ifail.gt.0)return
        nfil1=nfil+1
      endif
!  relax infeasible c/s
      hJ=0.D0
      hJt=0.D0
      do j=1,n
        i=abs(ls(j))
        if(i.gt.n)hJt=hJt+max(0.D0,dl(i),-du(i))
      enddo
      do j=n1,nm
        i=abs(ls(j))
        if(i.gt.n)then
          if(r(i).lt.0.D0)then
            hJ=hJ+max(0.D0,dl(i),-du(i))
            if(ls(j).ge.0)then
              du(i)=dl(i)
              dl(i)=-ainfty
              cstype(i-n)='A'
            else
              dl(i)=du(i)
              du(i)=ainfty
              cstype(i-n)='Z'
            endif
          else
            hJt=hJt+max(0.D0,dl(i),-du(i))
          endif
        endif
      enddo
!     print *,'phase 1 filter entries followed by (hJt,hJ)'
!     do i=nfil1,nfil
!       print 5,filh(i),filf(i)
!     enddo
!     print 5,hJt,hJ
!     print *,'cstype = ',(cstype(i),i=1,m)
!     print 2,'hJt,hJ',hJt,hJ
      if(iprint.ge.1)write(nout,1000)itn,hJt,hJ,rho
!     if(hJt.gt.tol)then
!       call addfil(hJt,hJ,filh,filf,nfil1,nfil,mxf,ifail)
!       if(ifail.gt.0)return
        call testfil(hJt,hJ,filh,filf,nfil1,nfil,ifail)
        if(ifail.eq.1)then
          if(iprint.gt.1)write(nout,*)'l1 solution not acceptable'
          dnorm=0.D0
          do i=1,n
            dnorm=max(dnorm,abs(d(i)))
          enddo
          rho=5.D-1*dnorm
           goto 10
        endif
!     endif
!  collect multipliers from l1 subproblem
      do j=1,n
        i=abs(ls(j))
        if(i.gt.n)then
          if(ls(j).gt.0)then
            ws(nal+i-n)=r(i)
          else
            ws(nal+i-n)=-r(i)
          endif
        endif
        ws(naal+j)=0.D0
      enddo
!     print 4,'al =',(ws(i),i=nal1,nal+m)
      do i=1,m
        call saipy(ws(nal+i),ws(last1),lws(nla1),i,ws(naal1),n)
      enddo
20    continue
      if(itn.eq.maxit)then
        if(iprint.gt.1)write(nout,*)'itn.ge.maxit'
        ifail=5
        return
      endif
      iph=1
!     mode=2
      k=0
!  solve LCP subproblem
!     print 4,'x =',(ws(i),i=nx1,nx+n)
!     print 4,'c =',(ws(i),i=ncx1,ncx+m)
!     print 4,'al =',(ws(i),i=nal1,nal+m)
      do i=1,n
        dl(i)=max(-rho,bl(i)-ws(nx+i))
        du(i)=min(rho,bu(i)-ws(nx+i))
        d(i)=0.D0
      enddo
      alc=scpr(0.D0,ws(nal1),ws(ncx1),m)
!     print *,'solve phase 1 LCP subproblem',itn
      iii=0
!     if(itn.eq.11)iii=2
      call glcpd(n,m,k,kmax,maxg,ws(last1),lws(nla1),d,dl,du,phi, &
        fmin,g,r,w,e,ls,alp,lp,mlp,ipeq,ws,lws,cstype,v,nv,rgtol, &
        mode,ifail,mxgr,iii,0)
      nft=nft+nfn
      ngt=ngt+ngr
!     print 1,'nfn,ngr',nfn,ngr
!     print 4,'d =',(d(i),i=1,n)
!     print 4,'r =',(r(i),i=1,nm)
!     print 3,'ls =',(ls(i),i=1,nm)
!     print 1,'ipeq,k,ifail =',ipeq,k,ifail
      itn=itn+1
      dnorm=0.D0
      do i=1,n
        dnorm=max(dnorm,abs(d(i)))
      enddo
!     print 2,'dnorm,rho',dnorm,rho
      if(ifail.eq.3)then
        if(iprint.gt.1) &
          write(nout,*)'phase 1 LCP problem is infeasible'
        if(dnorm.le.htol) goto 10
        rho=5.D-1*dnorm
         goto 10
      elseif(ifail.eq.1)then
        if(iprint.gt.1) &
          write(nout,*)'phase 1 LCP subproblem is unbounded'
         goto 99
      elseif(ifail.gt.5)then
        if(iprint.gt.1) &
          write(nout,*)'malfunction in phase 1 LCP subproblem'
         goto 99
      endif
      hxdJt=0.D0
      hxdJ=0.D0
      do i=1,m
        ci=ws(ncxd+i)
        if(cstype(i).eq.'N')then
          hxdJt=hxdJt+max(0.D0,bl(n+i)-ci,ci-bu(n+i))
        else
          hxdJ=hxdJ+max(0.D0,bl(n+i)-ci,ci-bu(n+i))
        endif
      enddo
      if(iprint.ge.1)write(nout,1001)itn,hxdJt,hxdJ,rgnorm,dnorm,rho
!     print 4,'x+d =',(ws(nxd+i),i=1,n)
!     print 4,'c at x+d =',(ws(ncxd+i),i=1,m)
      if(hxdJt.le.htol.and.dnorm.le.htol) goto 40
!     print *,'phase 1 filter entries followed by (hJt,hJ)'
!     do i=nfil1,nfil
!       print 5,filh(i),filf(i)
!     enddo
!     print 5,hJt,hJ
      hxd=hxdJt+hxdJ
      if(hxd.ge.ubd)then
        if(iprint.gt.1)write(nout,*)'upper bound on h exceeded (1)'
        rho=max(1.D-1,5.D-1*h/hxd)*dnorm
         goto 10
      endif
      dq=hJ-phi
      df=hJ-hxdJ
!     print 2,'dq,df',dq,df
!  filter test for LCP solution
      call testfil(hxdJt,hxdJ,filh,filf,nfil1,nfil,ifail)
      if(ifail.eq.0)call testfil(hxdJt,hxdJ,hJt,hJ,1,1,ifail)
!     print 6,'hxdJt,hxdJ,ifail',hxdJt,hxdJ,ifail
      if(ifail.eq.1.or.(dq.ge.tol.and.df.lt.sigma*dq))then
        if(hxdJt.eq.0.D0.or.dq.lt.tol)then
          if(iprint.gt.1)write(nout,*)'hxdJt.eq.0.D0.or.dq.lt.tol'
          rho=max(1.D-1,min(5.D-1*h/hxd,5.D-1))*dnorm
           goto 10
        endif
!  projection step
        nv=1
        v(1)=1.D0
        iph=-1
        do i=1,n
          ws(naal+i)=ws(nx+i)
        enddo
30      continue
        hxJt=hxdJt
        hxJ=hxdJ
        do i=1,n
          ws(nx+i)=ws(nxd+i)
          dl(i)=max(-rho,bl(i)-ws(nx+i))
          du(i)=min(rho,bu(i)-ws(nx+i))
!         dl(i)=bl(i)-ws(nx+i)
!         du(i)=bu(i)-ws(nx+i)
          d(i)=0.D0
        enddo
        do i=1,m
          ci=ws(ncxd+i)
          if(cstype(i).eq.'A')then
            dl(n+i)=-ainfty
            du(n+i)=bl(n+i)-ci
          elseif(cstype(i).eq.'Z')then
            dl(n+i)=bu(n+i)-ci
            du(n+i)=ainfty
          else
            dl(n+i)=bl(n+i)-ci
            du(n+i)=bu(n+i)-ci
          endif
        enddo
!       mode=2
        k=0
!  solve projection subproblem
!       print *,'solve phase 1 projection subproblem',itn
        iii=0
!       if(itn.eq.68)iii=1
        call glcpd(n,m,k,kmax,maxg,ws(next1),lws(nla1),d,dl,du,phi,0.D0, &
          g,r,w,e,ls,alp,lp,mlp,ipeq,ws,lws,cstype,v,nv,rgtol, &
          mode,ifail,mxgr,iii,0)
!       print 4,'d =',(d(i),i=1,n)
!       print 4,'r =',(r(i),i=1,nm)
!       print 3,'ls =',(ls(i),i=1,nm)
!       print 1,'ipeq,k,ifail =',ipeq,k,ifail
        if(ifail.eq.3)then
          if(iprint.gt.1) &
            write(nout,*)'phase 1 projection problem is infeasible'
          do i=1,n
            ws(nx+i)=ws(naal+i)
          enddo
          rho=max(1.D-1,min(5.D-1*h/hxd,5.D-1))*dnorm
           goto 10
        elseif(ifail.gt.5)then
          if(iprint.gt.1) &
            write(nout,*)'malfunction in phase 1 projection subproblem'
           goto 99
        endif
        do i=1,n
          ws(nxd+i)=ws(nx+i)+d(i)
        enddo
        call functions(n,m,ws(nxd1),fxd,ws(ncxd1),ws,lws)
        call gradients(n,m,ws(nxd1),ws(next1),ws,lws)
        hxdJt=0.D0
        hxdJ=0.D0
        do i=1,m
          if(cstype(i).eq.'N')then
            hxdJt=hxdJt+max(0.D0,bl(n+i)-ws(ncxd+i),ws(ncxd+i)-bu(n+i))
          else
            hxdJ=hxdJ+max(0.D0,bl(n+i)-ws(ncxd+i),ws(ncxd+i)-bu(n+i))
          endif
        enddo
        if(iprint.ge.1)write(nout,1002)itn,hxdJt,hxdJ
!       print 4,'c at x+d =',(ws(ncxd+i),i=1,m)
!  filter test for projection solution
        hxd=hxdJt+hxdJ
        if(hxd.ge.ubd)then
          if(iprint.gt.1)write(nout,*)'upper bound on h exceeded (2)'
          rho=max(1.D-1,5.D-1*h/hxd)*dnorm
           goto 10
        endif
        df=hJ-hxdJ
!       print 2,'dq,df',dq,df
        call testfil(hxdJt,hxdJ,filh,filf,nfil1,nfil,ifail)
        if(ifail.eq.0)call testfil(hxdJt,hxdJ,hJt,hJ,1,1,ifail)
!       if(ifail.eq.1)print 2,'hxdJt/hxJt =',hxdJt/hxJt
!       print 6,'project: hxdJt,hxdJ,ifail',hxdJt,hxdJ,ifail
        if(ifail.eq.1.or.df.lt.sigma*dq)then
          if(hxdJt.le.8.D-1*hxJt)then
            df=hJ-(hxJt*hxdJ-hxdJt*hxJ)/(hxJt-hxdJt)
            if(df.ge.sigma*dq) goto 30
          endif
          do i=1,n
            ws(nx+i)=ws(naal+i)
          enddo
          rho=max(1.D-1,min(5.D-1*h/hxd,5.D-1))*dnorm
          if(iprint.gt.1)write(nout,*)'phase 1 projection step fails'
           goto 10
        endif
!       print *,'accept projection step (1)'
        do i=1,n
          ws(nx+i)=ws(naal+i)
        enddo
      endif
40    continue
!  accept LCP (iph=1) or projection (iph=-1) solution
      if(dq.lt.tol)then
        call addfil(hJt,hJ,filh,filf,nfil1,nfil,mxf,ifail)
        if(ifail.gt.0)return
      endif
      do i=1,n
        ws(nx+i)=ws(nxd+i)
        ws(naal+i)=0.D0
      enddo
      do i=1,m
        ws(ncx+i)=ws(ncxd+i)
      enddo
      call iexch(last1,next1)
      h=hxdJt+hxdJ
      f=fxd
      if(dnorm.eq.rho)rho=2.D0*rho
      if(h.le.htol) goto 10
!  check for situations where the l1 partition needs recalculating ...
!  if there are any active relaxed c/s
      do j=1,n
        i=abs(ls(j))-n
!       if(i.gt.0.and.cstype(i).ne.'N')print 1,'active relaxed c/s',i
        if(i.gt.0.and.cstype(i).ne.'N') goto 10
      enddo
!  or any infeasible relaxed c/s
      do i=1,m
!       if((cstype(i).eq.'A'.and.ws(ncx+i).ge.bl(n+i)).or.
!    *    (cstype(i).eq.'Z'.and.ws(ncx+i).le.bu(n+i)))
!    *    print 1,'infeasible relaxed c/s',i
        if((cstype(i).eq.'A'.and.ws(ncx+i).ge.bl(n+i)).or. &
          (cstype(i).eq.'Z'.and.ws(ncx+i).le.bu(n+i))) goto 10
      enddo
      hJt=hxdJt
      hJ=hxdJ
      if(iph.eq.-1) goto 10
      if(hxdJt.le.htol.and.dnorm.le.htol)then
        if(iprint.gt.1)write(nout,*)'locally infeasible problem'
!       print 2,'hxdJt.le.htol.and.dnorm.le.htol'
        ifail=3
        return
      endif
!  collect LCP multipliers
      do i=1,m
        ws(nal+i)=0.D0
      enddo
      do j=1,n
        i=abs(ls(j))
        if(i.gt.n)then
          if(ls(j).gt.0)then
            ws(nal+i-n)=r(i)
          else
            ws(nal+i-n)=-r(i)
          endif
        endif
      enddo
!     print 4,'al =',(ws(i),i=nal1,nal+m)
      do i=1,m
        ci=ws(ncx+i)
        if(cstype(i).eq.'A')then
          dl(n+i)=-ainfty
          du(n+i)=bl(n+i)-ci
        elseif(cstype(i).eq.'Z')then
          dl(n+i)=bu(n+i)-ci
          du(n+i)=ainfty
        else
          dl(n+i)=bl(n+i)-ci
          du(n+i)=bu(n+i)-ci
        endif
        call saipy(ws(nal+i),ws(last1),lws(nla1),i,ws(naal1),n)
      enddo
       goto 20

50    continue
!  Phase 2 code
!  collect multipliers from LP subproblem
      do j=1,n
        i=abs(ls(j))
        if(i.gt.n)then
          if(ls(j).gt.0)then
            ws(nal+i-n)=r(i)
          else
            ws(nal+i-n)=-r(i)
          endif
        endif
        ws(naal+j)=0.D0
      enddo
      do i=1,m
        call saipy(ws(nal+i),ws(last1),lws(nla1),i,ws(naal1),n)
      enddo
60    continue
      if(itn.eq.maxit)then
        if(iprint.gt.1)write(nout,*)'itn.ge.maxit'
        ifail=5
        return
      endif
!     print 2,'h,f =',h,f
      iph=2
      k=0
!     mode=2
!  solve LCP subproblem
!     print 4,'x =',(ws(i),i=nx1,nx+n)
!     print 4,'al =',(ws(i),i=nal1,nal+m)
!     print 4,'c =',(ws(i),i=ncx1,ncx+m)
      do i=1,n
        dl(i)=max(-rho,bl(i)-ws(nx+i))
        du(i)=min(rho,bu(i)-ws(nx+i))
        d(i)=0.D0
      enddo
      alc=scpr(0.D0,ws(nal1),ws(ncx1),m)
!     print *,'solve phase 2 LCP subproblem',itn
      iii=0
!     if(itn.eq.164)iii=1
      call glcpd(n,m,k,kmax,maxg,ws(last1),lws(nla1),d,dl,du,phi,fmin, &
        g,r,w,e,ls,alp,lp,mlp,ipeq,ws,lws,cstype,v,nv,rgtol, &
        mode,ifail,mxgr,iii,0)
      nft=nft+nfn
      ngt=ngt+ngr
!     print 1,'nfn,ngr',nfn,ngr
!     print 4,'d =',(d(i),i=1,n)
!     print 4,'r =',(r(i),i=1,nm)
!     print 3,'ls =',(ls(i),i=1,nm)
!     print 1,'ipeq,k,ifail =',ipeq,k,ifail
      itn=itn+1
      dnorm=0.D0
      do i=1,n
        dnorm=max(dnorm,abs(d(i)))
      enddo
!     print 2,'dnorm,rho',dnorm,rho
      if(ifail.eq.3)then
        if(iprint.gt.1) &
          write(nout,*)'phase 2 LCP problem is infeasible'
!       mode=2
         goto 15
!       if(dnorm.le.htol) goto 10
!       rho=5.D-1*dnorm
!        goto 10
      elseif(ifail.eq.1)then
        if(iprint.gt.1) &
          write(nout,*)'phase 2 LCP subproblem is unbounded'
         goto 99
      elseif(ifail.gt.5)then
        if(iprint.gt.1) &
          write(nout,*)'malfunction in phase 2 LCP subproblem'
         goto 99
      endif
      hxd=0.D0
      do i=1,m
        ci=ws(ncxd+i)
        hxd=hxd+max(0.D0,bl(n+i)-ci,ci-bu(n+i))
      enddo
      if(iprint.ge.1)write(nout,2001)itn,hxd,fxd,rgnorm,dnorm,rho
!     print 4,'c at x+d =',(ws(ncxd+i),i=1,m)
      if(hxd.le.htol.and.(fxd.le.fmin.or.dnorm.le.htol)) goto 80
!     print *,'phase 2 filter entries followed by (h,f)'
!     do i=1,nfil
!       print 5,filh(i),filf(i)
!     enddo
!     print 5,h,f
!  filter test for LCP solution
      if(hxd.ge.ubd)then
        if(iprint.gt.1)write(nout,*)'upper bound on h exceeded (3)'
        rho=max(1.D-1,5.D-1*h/hxd)*dnorm
         goto 10
      endif
      dq=f-phi
      df=f-fxd
!     print 2,'dq,df',dq,df
      call testfil(hxd,fxd,filh,filf,1,nfil,ifail)
      if(ifail.eq.0)call testfil(hxd,fxd,h,f,1,1,ifail)
!     print 6,'hxd,fxd,ifail',hxd,fxd,ifail
      if(ifail.eq.1.or.(dq.ge.tol.and.df.lt.sigma*dq))then
        if(hxd.eq.0.D0.or.dq.lt.tol)then
          rho=5.D-1*dnorm
          if(iprint.gt.1)write(nout,*)'hxd.eq.0.D0.or.dq.lt.tol'
          rho=max(1.D-1,min(5.D-1*h/hxd,5.D-1))*dnorm
           goto 10
        endif
!  projection step
        nv=1
        v(1)=1.D0
        iph=-2
        do i=1,n
          ws(naal+i)=ws(nx+i)
        enddo
70      continue
        hx=hxd
        fx=fxd
        do i=1,n
          ws(nx+i)=ws(nxd+i)
          dl(i)=max(-rho,bl(i)-ws(nx+i))
          du(i)=min(rho,bu(i)-ws(nx+i))
!         dl(i)=bl(i)-ws(nx+i)
!         du(i)=bu(i)-ws(nx+i)
          d(i)=0.D0
        enddo
        do i=1,m
          ci=ws(ncxd+i)
          dl(n+i)=bl(n+i)-ci
          du(n+i)=bu(n+i)-ci
        enddo
!       mode=2
        k=0
!       print 4,'x =',(ws(nx+i),i=1,n)
!       print 4,'v =',(v(i),i=1,nv)
!  solve projection subproblem
!       print *,'solve phase 2 projection subproblem'
        iii=0
!       if(itn.eq.9)iii=3
        call glcpd(n,m,k,kmax,maxg,ws(next1),lws(nla1),d,dl,du,phi,0.D0, &
          g,r,w,e,ls,alp,lp,mlp,ipeq,ws,lws,cstype,v,nv,rgtol, &
          mode,ifail,mxgr,iii,0)
!       print 4,'d =',(d(i),i=1,n)
!       print 4,'r =',(r(i),i=1,nm)
!       print 3,'ls =',(ls(i),i=1,nm)
!       print 1,'ipeq,k,ifail =',ipeq,k,ifail
        if(ifail.eq.3)then
          if(iprint.gt.1) &
            write(nout,*)'phase 2 projection problem is infeasible'
          do i=1,n
            ws(nx+i)=ws(naal+i)
          enddo
          rho=5.D-1*dnorm
           goto 10
        elseif(ifail.gt.5)then
          if(iprint.gt.1) &
           write(nout,*),'malfunction in phase 2 projection subproblem'
           goto 99
        endif
        do i=1,n
          ws(nxd+i)=ws(nx+i)+d(i)
        enddo
!       print 4,'xd =',(ws(nxd+i),i=1,n)
        call functions(n,m,ws(nxd1),fxd,ws(ncxd1),ws,lws)
        call gradients(n,m,ws(nxd1),ws(next1),ws,lws)
        hxd=0.D0
        do i=1,m
          ci=ws(ncxd+i)
          hxd=hxd+max(0.D0,bl(n+i)-ci,ci-bu(n+i))
        enddo
        if(iprint.ge.1)write(nout,2002)itn,hxd,fxd
!       print 4,'x+d =',(ws(nxd+i),i=1,n)
!       print 4,'c at x+d =',(ws(ncxd+i),i=1,m)
!  filter test for projection solution
        if(hxd.ge.ubd)then
          if(iprint.gt.1)write(nout,*)'upper bound on h exceeded (4)'
          rho=max(1.D-1,5.D-1*h/hxd)*dnorm
           goto 10
        endif
        df=f-fxd
        call testfil(hxd,fxd,filh,filf,1,nfil,ifail)
        if(ifail.eq.0)call testfil(hxd,fxd,h,f,1,1,ifail)
!       if(ifail.eq.1)print 2,'hxd/hx =',hxd/hx
!       print 6,'hxd,fxd,ifail',hxd,fxd,ifail
        if(ifail.eq.1.or.df.lt.sigma*dq)then
          if(hxd.le.8.D-1*hx)then
            df=f-(hx*fxd-hxd*fx)/(hx-hxd)
            if(df.ge.sigma*dq) goto 70
          endif
          do i=1,n
            ws(nx+i)=ws(naal+i)
          enddo
          rho=5.D-1*dnorm
          if(iprint.gt.1)write(nout,*)'phase 2 projection step fails'
           goto 10
        endif
!       print *,'accept phase 2 projection step'
        do i=1,n
          ws(nx+i)=ws(naal+i)
        enddo
      endif
80    continue
!  accept LCP (iph=2) or projection (iph=-2) solution
      if(dq.lt.tol)then
        call addfil(h,f,filh,filf,1,nfil,mxf,ifail)
        if(ifail.gt.0)return
      endif
      h=hxd
      f=fxd
      if(dnorm.eq.rho)rho=2.D0*rho
      do i=1,n
        ws(nx+i)=ws(nxd+i)
        ws(naal+i)=0.D0
      enddo
      do i=1,m
        ws(ncx+i)=ws(ncxd+i)
      enddo
      call iexch(last1,next1)
      if(iph.eq.-2) goto 10
      if(h.le.htol)then
        if(f.le.fmin)then
          if(iprint.gt.1) &
            write(nout,*)'phase 2 LCP unbounded'
          ifail=1
          return
        elseif(dnorm.le.htol)then
          if(iprint.gt.1)write(nout,*)'local NLP solution found'
!         print 2,'h.le.htol.and.dnorm.le.htol'
          ifail=0
          return
        endif
      endif
!  collect LCP multipliers
      do i=1,m
        ws(nal+i)=0.D0
      enddo
      do j=1,n
        i=abs(ls(j))
        if(i.gt.n)then
          if(ls(j).gt.0)then
            ws(nal+i-n)=r(i)
          else
            ws(nal+i-n)=-r(i)
          endif
        endif
      enddo
!     print 4,'al =',(ws(i),i=nal1,nal+m)
      do i=1,m
        ci=ws(ncx+i)
        dl(n+i)=bl(n+i)-ci
        du(n+i)=bu(n+i)-ci
        call saipy(ws(nal+i),ws(last1),lws(nla1),i,ws(naal1),n)
      enddo
!     iph=2
       goto 60
99    continue
      if(ifail.eq.7)return
      ifail=ifail+10
      return
      end

      block data nlp_defaults
      implicit double precision (a-h,o-z)
      common/defaultc/ainfty,ubd,mlp,mxf
      common/ngrc/mxgr
      data  ainfty, ubd, mlp, mxf,  mxgr &
          / 1.D20,  1.D4, 50,  50, 1000000/
      end

      subroutine testfil(h,f,filh,filf,nfil1,nfil,ifail)
      implicit double precision (a-h,o-z)
      dimension filh(*),filf(*)
      parameter (beta=99999.D-5,gamma=1.D-5)
      ifail=0
!     if(h.eq.0.D0)return
      hd=h/beta
      fp=f+gamma*h
      do i=nfil1,nfil
        if(hd.ge.filh(i).and.fp.gt.filf(i))then
!       if(hd.gt.filh(i).and.fp.gt.filf(i))then
          ifail=1
          return
        endif
      enddo
      return
      end

      subroutine addfil(h,f,filh,filf,nfil1,nfil,mxf,ifail)
      implicit double precision (a-h,o-z)
      dimension filh(*),filf(*)
      do i=nfil,nfil1,-1
        if(h.le.filh(i).and.f.le.filf(i))then
          filh(i)=filh(nfil)
          filf(i)=filf(nfil)
          nfil=nfil-1
        endif
      enddo
      if(nfil.ge.mxf)then
        ifail=8
        return
      endif
      nfil=nfil+1
      filh(nfil)=h
      filf(nfil)=f
      ifail=0
      return
      end

      subroutine funct(n,d,phi,ws,lws,cstype)
      implicit double precision (a-h,o-z)
      dimension d(*),ws(*),lws(*)
      character cstype(*)
      common/defaultc/ainfty,ubd,mlp,mxf
      common/functc/fxd,alc,m,iph,last1,next1,nx,nx1, &
        nal,nal1,naal,naal1,nxd,nxd1,ncx,ncx1,ncxd,ncxd1,nla1
2     format(A,6E15.7)
4     format(A/(5E15.7))
      if(iph.lt.0)then
!  projection subproblem
        phi=5.D-1*scpr(0.D0,d,d,n)
        return
      elseif(iph.eq.0)then
!  LP subproblem
        phi=aiscpr(n,ws(last1),lws(nla1),0,d,0.D0)
        return
      endif
      do i=1,n
        ws(nxd+i)=ws(nx+i)+d(i)
      enddo
!     print 4,'$x =',(ws(nx+i),i=1,n)
!     print 4,'$d =',(d(i),i=1,n)
!     print 4,'$x+d =',(ws(nxd+i),i=1,n)
      call functions(n,m,ws(nxd1),fxd,ws(ncxd1),ws,lws)
!     print 2,'$fxd =',fxd
!     print 4,'$al =',(ws(nal+i),i=1,m)
!     print 4,'$A.al =',(ws(naal+i),i=1,n)
!     print 4,'$cxd =',(ws(ncxd+i),i=1,m)
!     print *,'$cstype =',(cstype(i),i=1,m)
      if(iph.eq.2)then
        phi=scpr(-scpr(-fxd-alc,ws(nal1),ws(ncxd1),m),ws(naal1),d,n)
      else
        phi=scpr(-scpr(-alc,ws(nal1),ws(ncxd1),m),ws(naal1),d,n)
        do i=1,m
          if(cstype(i).eq.'A')then
            phi=phi-ws(ncxd+i)
          elseif(cstype(i).eq.'Z')then
            phi=phi+ws(ncxd+i)
          endif
        enddo
      endif
!     print 2,'$fxd,cxd,phi =',fxd,ws(ncxd1),phi
!     print 2,'$phi =',phi
      return
      end

      subroutine grad(n,d,g,ws,lws,cstype)
      implicit double precision (a-h,o-z)
      dimension d(*),g(*),ws(*),lws(*)
      character cstype(*)
      common/defaultc/ainfty,ubd,mlp,mxf
      common/functc/fxd,alc,m,iph,last1,next1,nx,nx1, &
        nal,nal1,naal,naal1,nxd,nxd1,ncx,ncx1,ncxd,ncxd1,nla1
      common/maxac/maxa
      if(iph.lt.0)then
        do i=1,n
          g(i)=d(i)
        enddo
        return
      elseif(iph.eq.0)then
        do i=1,n
          g(i)=0.D0
        enddo
        call saipy(1.D0,ws(last1),lws(nla1),0,g,n)
        return
      endif
!     print 4,'£x+d =',(ws(nxd+i),i=1,n)
      call gradients(n,m,ws(nxd1),ws(next1),ws,lws)
!     print 4,'£a =',(ws(next1+i),i=0,maxa-1)
!     print 4,'£al =',(ws(nal+i),i=1,m)
      do i=1,n
        g(i)=ws(naal+i)
      enddo
!     print 4,'£Ak.al =',(g(i),i=1,n)
      do i=1,m
        call saipy(-ws(nal+i),ws(next1),lws(nla1),i,g,n)
      enddo
      if(iph.eq.2)then
        call saipy(1.D0,ws(next1),lws(nla1),0,g,n)
      else
        do i=1,m
          if(cstype(i).eq.'A')then
            call saipy(-1.D0,ws(next1),lws(nla1),i,g,n)
          elseif(cstype(i).eq.'Z')then
            call saipy(1.D0,ws(next1),lws(nla1),i,g,n)
          endif
        enddo
      endif
!     print 4,'£g =',(g(i),i=1,n)
1     format(A,15I5)
3     format(A/(15I5))
4     format(A/(5E15.7))
      return
      end
