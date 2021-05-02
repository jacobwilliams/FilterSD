!Christen this file sparseL.f
!ut here >>>>>>>>>>>>>>>>>
!***************** sparse matrix routines for manipulating L *******************

!           ***************************************************
!           Basis matrix routines for bqpd with sparse matrices
!           ***************************************************

!  These routines form and update L-Implicit-U factors LPB=U of a matrix B
!  whose columns are the normal vectors of the active constraints. In this
!  method only the unit lower triangular matrix L and the diagonal of U (in
!  addition to the row permutation P) is stored. B is represented in block form

!    | I  A_2 |    where the last m1 columns (A_2 and A_1) come from the
!    | 0  A_1 |    general constraint normals (columns of the matrix A in bqpd)

!  and the remaining unit columns come from simple bounds. The matrix A must be
!  specified in sparse format and the user is referred to the file  sparseA.f.

!  The data structure used for L is that of a profile or skyline scheme, in
!  which the nontrivial rows of L are stored as dense row spikes. The use of
!  a Tarjan+spk1 ordering algorithm to control the length of these spikes has
!  proved quite effective. The factors are updated by a variant of the
!  Fletcher-Matthews method, which has proved very reliable in practice.
!  However the B matrix is re-factored every 30 updates to control growth in
!  the total spike length.

!  Workspace
!  *********
!  The user needs to supply storage for the rows of L, although the amount
!  required is unknown a-priori.
!  sparse.f requires
!     5*n+nprof          locations of real workspace, and
!     9*n+m              locations of integer workspace
!  where nprof is the space required for storing the row spikes of the L matrix.
!  Storage for sparseL.f is situated at the end of the workspace arrays ws
!  and lws in bqpd.
!  Allow as much space for nprof as you can afford: the routine will report if
!  there is not enough. So far 10^6 locations has proved adequate for problems
!  of up to 5000 variables.

!  In addition the current version of bqpd.f requires
!     kmax*(kmax+9)/2+2*n+m   locations of real workspace in ws
!     kmax                    locations of integer workspace in lws
!  The user is also allowed to reserve storage in ws and lws, for use in the
!  user-supplied routine gdotx. This storage is situated at the start of the
!  arrays ws and lws. The user specifies the amount required by
!  setting the parameters kk and ll in the common block
!     common/wsc/kk,ll,kkk,lll,mxws,mxlws
!  The user MUST also set mxws and mxlws to be (respectively) the total amount
!  of real and integer workspace for the arrays ws and lws.

!  Other information
!  *****************

!  The methodology behind the L-Implicit-U factors and the row spike storage
!  scheme for L is described in the references
!    Fletcher R., Dense Factors of Sparse Matrices, in "Approximation Theory
!    and Optimization. Tributes to M.J.D. Powell", (M.D. Buhmann and A. Iserles,
!    eds), Cambridge University Press (1997), pp. 145-166.
!  and
!    Fletcher R., Block Triangular Orderings and Factors for Sparse Matrices
!    in LP, in "Numerical analysis 1997" (D.F. Griffiths, D.J. Higham and
!    G.A. Watson, eds.), Pitman Research Notes in Mathematics 380, (1998),
!    Longman, Harlow, pp. 91-110.

!  The file contains routines for solving systems with B or its transpose
!  which might be of use in association with bqpd. These routines are
!  documented below.

!  Steepest edge coefficients e(i) are also updated in these routines

!  Copyright, University of Dundee (R.Fletcher), January 1998
!  Current version dated 16/04/02

      subroutine start_up(n,nm,nmi,a,la,nk,e,ls,aa,ll,mode,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(0:*),e(*),ls(*),aa(*),ll(*)
      common/noutc/nout
      common/wsc/kk,ll_,kkk,lll,mxws,mxlws
      common/epsc/eps,tol,emin
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof, &
        lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/refactorc/nup,nfreq
      nfreq=min(30,nfreq)
      nup=0
      ns=kk+kkk+5*n
      nt=ll_+lll+8*n+nmi
      nprof=mxws-ns
      if(nprof.le.0 .or. nt>mxlws) then
        write(nout,*)'not enough real (ws) or integer (lws) workspace'
        write(nout,*)'you give values for mxws and mxlws as',mxws,mxlws
        write(nout,*)'minimum values for mxws and mxlws are',ns,nt
        ifail=7
        return
      end if
3     format(A/(20I5))
4     format(A/(5E15.7))
!  set storage map for sparse factors
      ns=n
      ns1=ns+1
      nt=ns+n
      nt1=nt+1
      nu=nt+n
      nu1=nu+1
      nx=nu+n
      nx1=nx+1
      np=nx+n
      np1=np+1
      lc=n
      lc1=lc+1
      li=lc+n
      li1=li+1
      lm=li+nmi
      lm1=lm+1
      lp=lm+n
      lp1=lp+1
      lq=lp+n
      lq1=lq+1
      lr=lq+n
      lr1=lr+1
      ls_=lr+n
      ls1=ls_+1
      lt=ls_+n
      lt1=lt+1
      m=nm-n
      mp=-1
      mq=-1
!     write(nout,*)'ls',(ls(ij),ij=1,nk)
      if(mode>=3) then
        call re_order(n,nm,a,la(1),la(la(0)),ll,ll(lc1),ll(li1), &
          ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1), &
          nprof,ifail)
        if(ifail>=1) then
!         write(nout,*)'failure in re_order (1)'
          if(ifail.eq.7) return
          mode=2
           goto 1
        end if
        call re_factor(n,nm,a,la,ll,ll(lc1),ll(li1), &
          ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1), &
          nprof,aa,ifail)
        if(ifail.eq.7) return
        call check_L(n,aa,ll(lp1),ifail)
        if(ifail.eq.1) then
          mode=2
           goto 1
        end if
        if(nk.eq.n) return
!  reset ls from e
        do j=1,nk
          i=-ls(j)
          if(i>0)e(i)=-e(i)
        end do
        j=0
        nk=nmi
        do i=1,nmi
          if(e(i)/=0.D0) then
            j=j+1
            if(e(i)>0.D0) then
              ls(j)=i
            else
              ls(j)=-i
              e(i)=-e(i)
            end if
          else
            ls(nk)=i
            nk=nk-1
          end if
        end do
        if(j/=n) then
          write(nout,*)'malfunction in reset sequence in start_up'
          stop
        end if
        return
      end if
1     continue
      if(emin.eq.0.D0) then
!  set a lower bound on e(i): setting emin=0.D0 will force emin to be recalculated: do this only if mode<3
        emin=1.D0
        do i=1,nmi-n
          emin=max(emin,ailen(n,a,la,i))
        end do
        emin=1.D0/emin
      end if
      do i=1,n
        ll(i)=i
        ll(li+i)=i
        e(i)=1.D0
      end do
      do i=n+1,nmi
        ll(li+i)=0
        e(i)=0.D0
      end do
      nu_=0
      if(mode/=0) then
!  shift designated bounds to end and order the resulting rows and columns
        do j=1,nk
          i=abs(ls(j))
          if(i.le.n) then
            nn=n-nu_
            nu_=nu_+1
            call iexch(ls(nu_),ls(j))
            ii=ll(li+i)
            ll(ii)=ll(nn)
            ll(li+ll(ii))=ii
            ll(nn)=i
            ll(li+i)=nn
          end if
        end do
        call order(n,nu_,nk,la,ll,ls,ll(li1),ll(lp1),ll(lq1),ll(lr1), &
          aa(np1),nprof,ifail)
        if(ifail>0) return
      end if
      call factor(n,nmi,nu_,nk,a,la,e,ls,aa(ns1),aa(nt1),aa(nu1), &
        aa(nx1),ll,ll(lc1),ll(li1),ll(lm1),ll(lp1),ll(lq1),ll(lr1), &
        ll(ls1),aa(np1),nprof,aa,ifail)
      if(ifail>0) return
!     write(nout,*)'steepest edge coefficients',(e(ij),ij=1,nm)
!     emax=0.D0
!     do i=1,nm
!       if(e(i)>0.D0) then
!         call eptsol(n,a,la,i,a,aa(ns1),aa(nt1),aa,aa(np1),
!    *      ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
!         ei=xlen(0.D0,aa(ns1),n)
!         ei=sqrt(scpr(0.D0,aa(ns1),aa(ns1),n))
!         emax=max(emax,abs(ei-e(i)))
!       end if
!     end do
!     if(emax>=tol)
!    *  write(nout,*)'error in steepest edge coefficients =',emax
      return
      end

      subroutine refactor(n,nm,a,la,aa,ll,ifail)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),aa(*),ll(*)
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof, &
        lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/noutc/nout
!     write(nout,*)'refactor'
      m=nm-n
      call re_order(n,nm,a,la(1),la(la(0)),ll,ll(lc1),ll(li1), &
        ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1), &
        nprof,ifail)
      if(ifail>=1) then
!       write(nout,*)'failure in re_order (2)'
        return
      end if
      call re_factor(n,nm,a,la,ll,ll(lc1),ll(li1),ll(lm1), &
        ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1), &
        nprof,aa,ifail)
      if(ifail.eq.7) return
      call check_L(n,aa,ll(lp1),ifail)
      return
      end

      subroutine pivot(p,q,n,nm,a,la,e,aa,ll,ifail,info)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(0:*),e(*),aa(*),ll(*),info(*)
      common/noutc/nout
      common/iprintc/iprint
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof, &
        lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/mxm1c/mxm1
      common/refactorc/nup,nfreq
      common/epsc/eps,tol,emin
!     write(nout,*)'pivot: p,q =',p,q
      ifail=0
      if(p/=mp) then
        call eptsol(n,a,la,p,a,aa(ns1),aa(nt1),aa,aa(np1), &
          ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
        if(p>n) then
          e(p)=xlen(0.D0,aa(ns1+m2),m1)
        else
          e(p)=xlen(1.D0,aa(ns1+m2),m1)
        end if
        epp=e(p)
        mp=p
      end if
      if(q/=mq) then
        call aqsol(n,a,la,q,a,aa(nt1),aa(nx1),aa,aa(np1), &
          ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
        mq=q
      end if
!  update steepest edge coefficients
      tp=aa(nt+ll(li+p))
      if(tp.eq.0.D0)tp=eps
      ep=e(p)
      eq=2.D0/ep
!     do i=1,m2-1
!       aa(nu+i)=0.D0
!     end do
!     do i=m2,n
      do i=1,n
        aa(nu+i)=eq*aa(ns+i)
      end do
      call aqsol(n,a,la,-1,a,aa(nu1),aa(nx1),aa,aa(np1), &
        ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
!     write(nout,*)'row perm',(ll(ij),ij=1,n)
!     write(nout,*)'column perm',(ll(lc+ij),ij=m2+1,n)
!     write(nout,*)'s =',(aa(ns+ij),ij=1,n)
!     write(nout,*)'t =',(aa(nt+ij),ij=1,n)
!     write(nout,*)'u =',(aa(nu+ij),ij=1,n)
      e(p)=0.D0
      eq=ep/tp
      do i=1,nm
        if(e(i)>0.D0) then
          j=ll(li+i)
          ei=e(i)
          wi=aa(nt+j)*eq
          awi=abs(wi)
          if(ei>=awi) then
            wi=wi/ei
            e(i)=max(emin,ei*sqrt(max(0.D0,1.D0+wi*(wi-aa(nu+j)/ei))))
          else
            wi=ei/wi
            e(i)=max(emin,awi*sqrt(max(0.D0,1.D0+wi*(wi-aa(nu+j)/ei))))
          end if
        end if
      end do
      e(q)=max(emin,abs(eq))
      info(1)=info(1)+1
      if(nup>=nfreq) then
!     if(nup>=30) then
!  refactorize L
        ip=ll(li+p)
        if(p>n) then
          m2=m2+1
          qq=ll(lc+m2)
          ll(lc+ip)=qq
          ll(li+qq)=ip
          ll(li+p)=0
        else
          ll(ip)=ll(m2)
          ll(li+ll(ip))=ip
          ll(m2)=p
          ll(li+p)=m2
        end if
        if(q>n) then
          ll(lc+m2)=q
          ll(li+q)=m2
          m2=m2-1
        else
          iq=ll(li+q)
          ll(iq)=ll(m2)
          ll(li+ll(iq))=iq
          ll(m2)=q
          ll(li+q)=m2
        end if
        m1=n-m2
        call re_order(n,nm,a,la(1),la(la(0)),ll,ll(lc1),ll(li1), &
          ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1), &
          nprof,ifail)
        if(ifail>=1) then
!         write(nout,*)'failure in re_order (3)'
          return
        end if
        call re_factor(n,nm,a,la,ll,ll(lc1),ll(li1), &
          ll(lm1),ll(lp1),ll(lq1),ll(lr1),ll(ls1),ll(lt1),aa(np1), &
          nprof,aa,ifail)
      else
!  update L
        call update_L(p,q,n,nm,a,la,ll,ll(lc1),ll(li1),ll(lm1),ll(lp1), &
          ll(lq1),ll(lr1),ll(ls1),aa(np1),nprof,aa,aa(ns1),ifail)
      end if
      if(ifail.eq.7) return
      mp=-1
      mq=-1
      call check_L(n,aa,ll(lp1),ifail)
!     write(nout,*)'steepest edge coefficients',(e(ij),ij=1,nm)
!     emax=0.D0
!     do i=1,nm
!       if(e(i)>0.D0) then
!         call eptsol(n,a,la,i,a,aa(ns1),aa(nt1),aa,aa(np1),
!    *      ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
!         ei=xlen(0.D0,aa(ns1),n)
!         ei=sqrt(scpr(0.D0,aa(ns1),aa(ns1),n))
!         emax=max(emax,abs(ei-e(i)))
!       end if
!     end do
!     if(emax>=tol)
!    *  write(nout,*)'error in steepest edge coefficients =',emax
      return
      end

      subroutine fbsub(n,jmin,jmax,a,la,q,b,x,ls,aa,ll,save)
      implicit double precision (a-h,r-z), integer (i-q)
      logical save
      dimension a(*),la(*),b(*),x(*),ls(*),aa(*),ll(*)

!  solves a system  B.x=b

!  Parameter list
!  **************
!   n   number of variables (as for bqpd)
!   jmin,jmax  (see description of ls below)
!   a,la   specification of QP problem data (as for bqpd)
!   q   an integer which, if in the range 1:n+m, specifies that the rhs vector
!       b is to be column q of the matrix A of general constraint normals.
!       In this case the parameter b is not referenced by fbsub.
!       If q=0 then b is taken as the vector given in the parameter b.
!   b(n)  must be set to the r.h.s. vector b (but only if q=0)
!   x(n+m)  contains the required part of the solution x, set according to the
!       index number of that component (in the range 1:n for a simple bound and
!       n+1:n+m for a general constraint)
!   ls(*)  an index vector, listing the components of x that are required.
!       Only the absolute value of the elements of ls are used (this allows
!       the possibility of using of the contents of the ls parameter of bqpd).
!       Elements of x in the range abs(ls(j)), j=jmin:jmax are set by fbsub.
!       These contortions allow bqpd to be independent of the basis matrix code.
!   aa(*)  real storage used by the basis matrix code (supply the vector
!       ws(lu1) with ws as in the call of bqpd and lu1 as in common/bqpdc/...)
!   ll(*)  integer storage used by the basis matrix code (supply the vector
!       lws(ll1) with lws as in the call of bqpd and ll1 as in common/bqpdc/...)
!   save   indicates if fbsub is to save its copy of the solution for possible
!       future use. We suggest that the user only sets save = .false.

      common/noutc/nout
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof, &
        lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
!     write(nout,*)'fbsub  q =',q
      if(save) then
        if(q/=mq) then
          call aqsol(n,a,la,q,b,aa(nt1),aa(nx1),aa,aa(np1), &
            ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
          mq=q
        end if
        do j=jmin,jmax
          i=abs(ls(j))
          x(i)=aa(nt+ll(li+i))
        end do
      else
        call aqsol(n,a,la,q,b,aa(nu1),aa(nx1),aa,aa(np1), &
          ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
        do j=jmin,jmax
          i=abs(ls(j))
          x(i)=aa(nu+ll(li+i))
        end do
      end if
      return
      end

      subroutine ztg(n,k,rg,lv,aa,ll)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension rg(*),lv(*),aa(*),ll(*)
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof, &
        lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls_,ls1,lt,lt1
!     print *,'aa =',(aa(nu+i),i=1,18)
      do j=1,k
        rg(j)=aa(nu+ll(li+lv(j)))
      end do
      return
      end

      subroutine tfbsub(n,a,la,p,b,x,aa,ll,ep,save)
      implicit double precision (a-h,r-z), integer (i-q)
      logical save
      dimension a(*),la(*),b(*),x(*),aa(*),ll(*)

!  solves a system  Bt.x=b

!  Parameter list
!  **************
!   n   number of variables (as for bqpd)
!   a,la   specification of QP problem data (as for bqpd)
!   p    an integer which, if in the range 1:n+m, specifies that the rhs vector
!        b is a unit vector appropriate to the position of p in the current
!        ordering. In this case b is not referenced by tfbsub.
!   b(n+m)  If p=0, this must be set to the r.h.s. vector b. Only the components
!        of b need be set, according to the index number of each component (in
!        the range 1:n for a simple bound and n+1:n+m for a general constraint)
!   x(n)  contains the solution x (in natural ordering)
!   aa(*)  real storage used by the basis matrix code (supply the vector
!       ws(lu1) with ws as in the call of bqpd and lu1 as in common/bqpdc/...)
!   ll(*)  integer storage used by the basis matrix code (supply the vector
!       lws(ll1) with lws as in the call of bqpd and ll1 as in common/bqpdc/...)
!   ep  if p/=0 and save is true, ep contains the l_2 length of x on exit
!   save  indicates if tfbsub is to save its copy of the solution for possible
!       future use. We suggest that the user only sets save = .false.

      common/noutc/nout
      common/sparsec/ns,ns1,nt,nt1,nu,nu1,nx,nx1,np,np1,nprof, &
        lc,lc1,li,li1,lm,lm1,lp,lp1,lq,lq1,lr,lr1,ls,ls1,lt,lt1
      common/factorc/m1,m2,mp,mq,lastr,irow
!     write(nout,*)'tfbsub  p =',p
      if(save) then
        if(p/=mp) then
          call eptsol(n,a,la,p,b,aa(ns1),aa(nt1),aa,aa(np1), &
            ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
          mp=p
        end if
        do i=1,n
          x(ll(i))=aa(ns+i)
        end do
        if(p>n) then
          ep=xlen(0.D0,aa(ns1+m2),m1)
        else if(p>0) then
          ep=xlen(1.D0,aa(ns1+m2),m1)
        end if
      else
        call eptsol(n,a,la,p,b,aa(nu1),aa(nt1),aa,aa(np1), &
          ll,ll(lc1),ll(li1),ll(lp1),ll(lq1))
        do i=1,n
          x(ll(i))=aa(nu+i)
        end do
      end if
!     write(nout,*)'x =',(x(i),i=1,n)
      return
      end

      subroutine newg
      common/factorc/m1,m2,mp,mq,lastr,irow
      mq=-1
      return
      end

!******** The following routines are internal to sparseL.f **************

      subroutine check_L(n,d,p,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension d(*),p(*)
      common/noutc/nout
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/epsc/eps,tol,emin
!     write(nout,*)'check_L'
      ifail=1
!     dmin=1.D37
      do k=nu+1,n
!       dmin=min(dmin,abs(d(k)))
        if(abs(d(k)).le.tol) return
      end do
!     write(nout,*)'dmin =',dmin
!     len=0
!     do i=1,n
!       len=len+p(i)
!     end do
!     write(nout,*)m1*(m1+1)/2,len+m1
!     write(nout,*)'m1 =',m1,'   file length =',len,'   total =',len+m1
      ifail=0
      return
      end

      subroutine aqsol(n,a,la,q,b,tn,xn,d,ws,lr,lc,li,pp,qq)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),b(*),tn(*),xn(*),d(*),ws(*), &
        lr(*),lc(*),li(*),pp(*),qq(*)
      common/noutc/nout
      common/factorc/m1,m2,mp,mq,lastr,irow
!     write(nout,*)'aqsol  q =',q
      if(q>0) then
        do i=1,n
          tn(i)=0.D0
        end do
        if(q.le.n) then
          tn(li(q))=1.D0
        else
          call iscatter(a,la,q-n,li,tn,n)
        end if
      else if(q.eq.0) then
        do i=1,n
          tn(li(i))=b(i)
        end do
      end if
!     write(nout,*)'tn =',(tn(i),i=1,n)
      do i=n,m2+1,-1
        ir=lr(i)
        pri=pp(ir)
        if(pri.eq.0) then
          xn(i)=tn(i)/d(i)
        else
          xn(i)=(scpr(tn(i),ws(qq(ir)+1),tn(i-pri),pri))/d(i)
        end if
        call isaipy(-xn(i),a,la,lc(i)-n,tn,n,lr,li)
      end do
      do i=m2+1,n
        tn(i)=xn(i)
      end do
!     write(nout,*)'tn =',(tn(i),i=1,n)
      return
      end

      subroutine eptsol(n,a,la,p,b,sn,tn,d,ws,lr,lc,li,pp,qq)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(*),b(*),sn(*),tn(*),d(*),ws(*), &
        lr(*),lc(*),li(*),pp(*),qq(*)
      common/noutc/nout
      common/iprintc/iprint
      common/epsc/eps,tol,emin
      common/factorc/m1,m2,mp,mq,lastr,irow
!     write(nout,*)'eptsol  p =',p
      if(p.eq.0) then
        do i=1,m2
          sn(i)=b(lr(i))
        end do
        do i=m2+1,n
          sn(i)=0.D0
        end do
        do i=m2+1,n
          j=lc(i)
          sn(i)=-aiscpri(n,a,la,j-n,sn,-b(j),lr,li)/d(i)
          ir=lr(i)
          pri=pp(ir)
          if(pri>0) call mysaxpy(sn(i),ws(qq(ir)+1),sn(i-pri),pri)
        end do
      else
        do i=1,n
          sn(i)=0.D0
        end do
        pr=li(p)
        if(p.le.n) then
          if(pr>m2) goto 1
          sn(pr)=1.D0
          do i=m2+1,n
            sn(i)=-aiscpri(n,a,la,lc(i)-n,sn,0.D0,lr,li)/d(i)
            ir=lr(i)
            pri=pp(ir)
            if(pri>0) call mysaxpy(sn(i),ws(qq(ir)+1),sn(i-pri),pri)
          end do
        else
          if(pr.le.m2) goto 1
          do i=m2+1,n
            bi=0.D0
            if(i.eq.pr)bi=-1.D0
            sn(i)=-aiscpri(n,a,la,lc(i)-n,sn,bi,lr,li)/d(i)
            ir=lr(i)
            pri=pp(ir)
            if(pri>0) call mysaxpy(sn(i),ws(qq(ir)+1),sn(i-pri),pri)
          end do
        end if
      end if
!     write(nout,*)'sn =',(sn(i),i=1,n)
      return
1     continue
      write(nout,*)'malfunction detected in eptsol: p =',p
      stop
      end

      subroutine order(n,nu,nc,la,lr,ls,li,p,q,r,ws,mxws,ifail)
      implicit integer (c-t)
      double precision ws
      dimension la(0:*),lr(*),ls(*),li(*),p(*),q(*),r(*),ws(*)
      common/noutc/nout
!     character star(1000,80)
!     write(nout,*)'order'
!  spk1 ordering on full matrix
      ifail=0
      if(nu.eq.n) return
!  set row and column counts and row-wise data structure
      nn=n-nu
      ii=mxws/nn
      do j=1,nn
        rowj=lr(j)
        p(rowj)=(j-1)*ii
        r(rowj)=0
      end do
      do j=nn+1,n
        r(lr(j))=0
      end do
1     continue
      do i=nu+1,nc
        coli=abs(ls(i))
        li(coli)=0
        jp=la(0)+coli-n
        do j=la(jp),la(jp+1)-1
          rowj=la(j)
          if(li(rowj).le.nn) then
            li(coli)=li(coli)+1
            r(rowj)=r(rowj)+1
            ij=p(rowj)+r(rowj)
            if(ij>mxws) then
              ij=mxws
              ifail=1
            end if
            ws(ij)=dble(coli)
          end if
        end do
      end do
!  check for no overlaps
      qrj=0
      do j=1,nn
        rowj=lr(j)
        if(p(rowj)<qrj)ifail=1
        qrj=p(rowj)+r(rowj)
        q(rowj)=qrj
        p(rowj)=p(rowj)+1
      end do
      if(ifail.eq.1 .or. qrj>mxws) then
        qrj=0
        do j=1,nn
          rowj=lr(j)
          p(rowj)=qrj
          qrj=qrj+r(rowj)
          r(rowj)=0
        end do
        if(qrj>mxws) then
          write(nout,*)'not enough space for ws in order:  mxws =',mxws
          ifail=7
          return
        end if
        ifail=0
         goto 1
      end if
      ifirstc=nu+1
      ifirstr=1
2     continue
!  move zero-column-count columns to lhs and find minimum column count
      mcc=n
      do i=ifirstc,nc
        coli=abs(ls(i))
        if(li(coli).eq.0) then
          call iexch(ls(i),ls(ifirstc))
          li(coli)=ifirstr-1
          ifirstc=ifirstc+1
        else
          mcc=min(mcc,li(coli))
        end if
      end do
!     write(nout,*)'ifirstc,ifirstr,mcc',ifirstc,ifirstr,mcc
!     write(nout,*)'lr =',(lr(j),j=1,n)
!     write(nout,*)'ls =',(ls(i),i=nu+1,nc)
!     write(nout,*)'row counts =',(r(lr(j)),j=1,n)
!     write(nout,*)'column counts =',(li(abs(ls(i))),i=nu+1,nc)
      if(ifirstc>nc) goto 4
!  apply tie-break rule
      tie=0
      do i=ifirstc,nc
        coli=abs(ls(i))
        if(li(coli).eq.mcc) then
          ti=0
          jp=la(0)+coli-n
          do j=la(jp),la(jp+1)-1
            rowj=la(j)
            if(li(rowj)>=ifirstr)ti=ti+r(rowj)
          end do
          if(ti>tie) then
            tie=ti
            mccc=coli
          end if
        end if
      end do
!     write(nout,*)'tie,mccc',tie,mccc
!  permute rows of m-c-c column to top and update column counts
      jp=la(0)+mccc-n
      do j=la(jp),la(jp+1)-1
        rowj=la(j)
        jr=li(rowj)
        if(jr<ifirstr) goto 3
        if(jr>nn) goto 3
        lr(jr)=lr(ifirstr)
        li(lr(jr))=jr
        lr(ifirstr)=rowj
        li(rowj)=ifirstr
        ifirstr=ifirstr+1
        do i=p(rowj),q(rowj)
          coli=int(ws(i))
          li(coli)=li(coli)-1
        end do
3       continue
      end do
       goto 2
4     continue
!  print star diagram
!     if(nc-nu>80 .or. n>1000) stop
!     write(nout,*)'spk1 ordering'
!     ij=li(abs(ls(nc)))
!     do i=1,ij
!       do j=1,nc-nu
!         star(i,j)=' '
!       end do
!     end do
!     do j=1,nc-nu
!       jp=la(0)+abs(ls(nu+j))-n
!       do i=la(jp),la(jp+1)-1
!         star(li(la(i)),j)='*'
!       end do
!     end do
!     do i=1,ij
!       write(nout,*)(star(i,j),j=1,nc-nu)
!     end do
!     write(nout,*)'lr =',(lr(i),i=1,n)
!     write(nout,*)'ls =',(ls(i),i=nu+1,nc)
!     write(nout,*)'lower profile =',(li(abs(ls(i))),i=nu+1,nc)
      return
      end

      subroutine factor(n,nm,nu,nc,a,la,e,ls,sn,tn,un,xn,lr,lc,li, &
        mao,p,q,r,s,ws,mxws,d,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      integer coli,r,s,rowi,rowp,tl,tu
      dimension a(*),la(0:*),e(*),ls(*),sn(*),tn(*),un(*),xn(*), &
        lr(*),lc(*),li(*),mao(*),p(*),q(*),r(*),s(*),ws(*),d(*)
!     character star(1000,80)
      common/factorc/m1,m2,mp,mq,lastr,irow
      common/iprintc/iprint
      common/refactorc/nup,nfreq
      common/epsc/eps,tol,emin
      common/noutc/nout
      parameter (thresh=1.D-1)
!  factorize LPA=U when A is rectangular
!    p(row) stores the number of stored elements of a natural row
!    q(row) stores the base address in ws of a natural row
!    r(row) stores the previous row stored in ws (or 0 if the first row in ws)
!    s(row) stores the next row stored in ws (or 0 if the last row in ws)
!    li(n+*) stores the lower profile of the sparse matrix
!    irow stores the natural row number of the initial row stored in ws
!    lastr stores the natural row number of the previous row put into ws
!     write(nout,*)'factor'
      nup=0
      lastr=0
      irow=0
      do i=1,n
        p(i)=0
      end do
      m1=0
      tl=1
      do ii=nu+1,nc
        coli=abs(ls(ii))
!       write(nout,*)'coli =',coli
        tu=li(coli)
        do i=1,n
          tn(i)=0.D0
        end do
        call iscatter(a,la,coli-n,li,tn,n)
        do i=m1,1,-1
          rowi=lr(i)
          pri=p(rowi)
          if(pri.eq.0) then
            xn(i)=tn(i)/d(i)
          else
            xn(i)=(scpr(tn(i),ws(q(rowi)+1),tn(i-pri),pri))/d(i)
          end if
          call isaipy(-xn(i),a,la,lc(i)-n,tn,n,lr,li)
        end do
        do i=1,m1
          tn(i)=xn(i)
        end do
        m1p=m1+1
!       write(nout,*)'lr =',(lr(i),i=1,n)
!       write(nout,*)'tn =',(tn(i),i=1,tu)
!  threshold pivot selection
        call linf(tu-m1,tn(m1p),z,iz)
        if(z.le.tol) then
          li(coli)=0
           goto 2
        end if
        zz=max(tol,z*thresh)
        do i=tl,tu
          q(lr(i))=m1p
        end do
!       write(nout,*)'q =',(q(lr(i)),i=m1p,tu)
        iz=iz+m1
        if(iz<tl) then
          z=0.D0
          qri=m1p
          do j=m1p,tu
            tnj=abs(tn(j))
            if(tnj>=zz) then
              qrj=q(lr(j))
              if(qrj.eq.qri) then
                if(tnj>z) then
                  z=tnj
                  iz=j
                end if
              else if(qrj>qri) then
                z=tnj
                iz=j
                qri=qrj
              end if
            end if
          end do
        end if
        tl=tu+1
!       write(nout,*)'zz,z,iz,m1,qri',zz,z,iz,m1,qri
        if(iz>m1p) then
          call rexch(tn(m1p),tn(iz))
          call iexch(lr(m1p),lr(iz))
          li(lr(m1p))=m1p
          li(lr(iz))=iz
        end if
        rowp=lr(m1p)
!  reset q values
        qrp=q(rowp)
        do i=m1p+1,tu
          if(abs(tn(i))>tol) then
            rowi=lr(i)
            if(qrp<q(rowi))q(rowi)=qrp
          end if
        end do
        tnp=tn(m1p)
        do i=1,n
          sn(i)=0.D0
        end do
        sn(m1p)=1.D0
        do i=1,m1
          sn(i)=-aiscpri(n,a,la,lc(i)-n,sn,0.D0,lr,li)/d(i)
          rowi=lr(i)
          pri=p(rowi)
          if(pri>0) call mysaxpy(sn(i),ws(q(rowi)+1),sn(i-pri),pri)
        end do
!       write(nout,*)'sn =',(sn(i),i=1,m1)
!  update steepest edge coefficients
        ep=e(rowp)
        e(rowp)=0.D0
        eq=2.D0/ep
        do i=1,n
          un(i)=eq*sn(i)
        end do
        do i=m1,1,-1
          rowi=lr(i)
          pri=p(rowi)
          if(pri.eq.0) then
            xn(i)=un(i)/d(i)
          else
            xn(i)=(scpr(un(i),ws(q(rowi)+1),un(i-pri),pri))/d(i)
          end if
          call isaipy(-xn(i),a,la,lc(i)-n,un,n,lr,li)
        end do
        do i=1,m1
          un(i)=xn(i)
        end do
!       write(nout,*)'un =',(un(i),i=1,n)
        eq=ep/tnp
        do i=1,nm
          if(e(i)>0.D0) then
            j=li(i)
            ei=e(i)
            wi=tn(j)*eq
            awi=abs(wi)
            if(ei>=awi) then
              wi=wi/ei
              e(i)=max(emin,ei*sqrt(max(0.D0,1.D0+wi*(wi-un(j)/ei))))
            else
              wi=ei/wi
              e(i)=max(emin,awi*sqrt(max(0.D0,1.D0+wi*(wi-un(j)/ei))))
            end if
          end if
        end do
        e(coli)=max(emin,abs(eq))
        do j=qrp,m1
          if(abs(sn(j))>tol) goto 1
        end do
        j=m1p
1       continue
        pri=m1p-j
        if(pri>0) then
          call newslot(rowp,pri,lastr,irow,p,q,r,s,ws,mxws,i,ifail)
          if(ifail>0) return
          p(rowp)=pri
          i=q(rowp)
          do j=j,m1
            i=i+1
            ws(i)=sn(j)
          end do
        end if
        m1=m1p
        ls(m1)=ls(ii)
        lc(m1)=coli
        li(coli)=m1
        d(m1)=tnp
2       continue
      end do
!  complete ls and reorder lr, lc and d
      do i=m1+1,n
        ls(i)=lr(i)
      end do
      j=n
      do i=1,nm
        if(e(i).eq.0.D0) then
          j=j+1
          ls(j)=i
        end if
      end do
      m2=n-m1
      do i=n,m2+1,-1
        lc(i)=lc(i-m2)
        li(lc(i))=i
        lr(i)=lr(i-m2)
        li(lr(i))=i
        d(i)=d(i-m2)
      end do
      do i=1,m2
        lr(i)=ls(m1+i)
        li(lr(i))=i
      end do
!  reset mao
      ilast=n
      ii=ilast
      do i=ilast,m2+1,-1
        mao(i)=ilast
        ii=min(ii,i-p(lr(i)))
        if(ii.eq.i)ilast=i-1
      end do
!     write(nout,*)'PAQ factors:  m1 =',m1
!     write(nout,*)'d =',(d(ij),ij=m2+1,n)
!     do j=m2+1,n
!       rowp=lr(j)
!       if(p(rowp)/=0) then
!         write(nout,*)'L(',rowp,')',
!    *      (ws(k),k=q(rowp)+1,q(rowp)+p(rowp))
!       end if
!     end do
!  print star diagram
!     write(nout,*)'factored ordering:  m1 =',m1
!     if(m1>80 .or. n>1000) stop
!     do i=1,n
!       do j=1,m1
!         star(i,j)=' '
!       end do
!     end do
!     do j=1,m1
!       jp=la(0)+lc(m2+j)-n
!       do i=la(jp),la(jp+1)-1
!         star(li(la(i)),j)='*'
!       end do
!     end do
!     do i=m2+1,n
!       write(nout,*)(star(i,j),j=1,m1)
!     end do
!     write(nout,*)'ls =',(ls(j),j=1,n)
!     write(nout,*)'s.e. coeffs =',(e(i),i=1,nm)
!     write(nout,*)'lr =',(lr(j),j=1,n)
!     write(nout,*)'lc =',(lc(j),j=m2+1,n)
!     write(nout,*)'mao =',(mao(j),j=m2+1,n)
!     call checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
      return
      end

      subroutine re_order(n,nm,a,la,point,lr,lc,li,mao,p,q,r,s, &
        t,ws,mxws,ifail)
      implicit double precision (a-h,u-z), integer (i-t)
      dimension a(*),la(*),point(0:*),lr(*),lc(*),li(*),mao(*), &
        p(*),q(*),r(*),s(*),t(*),ws(*)
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/noutc/nout
      logical backtrack
!     character star(1000,80)
!  print star diagram
!     if(n-nu>80 .or. n>1000) stop
!     write(nout,*)'initial ordering'
!     do i=1,n
!       do j=1,n-nu
!         star(i,j)=' '
!       end do
!     end do
!     do j=1,n-nu
!       ilp=lc(nu+j)-n
!       do i=point(ilp),point(ilp+1)-1
!         star(li(la(i)),j)='*'
!       end do
!     end do
!     do i=nu+1,n
!       write(nout,*)(star(i,j),j=1,n-nu)
!     end do
!     write(nout,*)'re_order'
      if(nu.eq.n) then
        ifail=0
        return
      end if
      m=nm-n
!  transversal search
      do iq=nu+1,n
        backtrack=.false.
        istack=nu
        inode=iq
        nodec=lc(inode)
        nodec_n=nodec-n
        lap=point(nodec_n+1)-point(nodec_n)
!       write(nout,*)'column node =',nodec,'  look-ahead rows =',
!    *    (la(j),j=point(nodec_n),point(nodec_n)+lap-1)
!  look-ahead loop
1       continue
          lap=lap-1
          nextr=la(point(nodec_n)+lap)
          inext=li(nextr)
          if(inext>=iq) goto 4
          if(lap>0) goto 1
          li(nodec)=0
2       continue
!  reassignment depth first search
        t(inode)=point(nodec_n+1)-point(nodec_n)
!       write(nout,*)'column node =',nodec,'  unfathomed rows =',
!    *    (la(j),j=point(nodec_n),point(nodec_n)+t(inode)-1)
3       continue
!  examine successor nodes
        if(t(inode).eq.0) then
          if(istack.eq.nu) then
            ifail=1
!           ifail=iq
!           write(nout,*)'exit: ifail =',iq
            return
          end if
          istack=istack-1
          backtrack=.true.
          if(istack.eq.nu) then
            inode=iq
          else
            inode=mao(istack)
          end if
!         write(nout,*)'backtrack to node at address =',inode
          nodec=lc(inode)
          nodec_n=nodec-n
!         write(nout,*)'column node =',nodec,'  unfathomed rows =',
!    *      (la(j),j=point(nodec_n),point(nodec_n)+t(inode)-1)
           goto 3
        end if
        t(inode)=t(inode)-1
        nextr=la(point(nodec_n)+t(inode))
        inext=li(nextr)
        if(inext.le.nu) goto 3
        if(t(inext)>=0) goto 3
!  extend depth first search
!       write(nout,*)'nextr,inext',nextr,inext
        inode=inext
!       write(nout,*)'put node address on stack'
        istack=istack+1
        mao(istack)=inode
!       write(nout,*)'stack =',(mao(j),j=nu+1,istack)
        nodec=lc(inode)
        nodec_n=nodec-n
        lap=li(nodec)
        if(lap.eq.0) goto 2
!       write(nout,*)'column node =',nodec,'  look-ahead rows =',
!    *    (la(j),j=point(nodec_n),point(nodec_n)+lap-1)
         goto 1
4       continue
!       write(nout,*)'new assignment found in row',nextr
!       write(nout,*)'istack,inext,nextr',istack,inext,nextr
!       if(istack>nu) write(nout,*)'stack =',(mao(j),j=nu+1,istack)
        li(nodec)=lap
!  perform row permutation
        lr(inext)=lr(iq)
        li(lr(inext))=inext
        inode=iq
        do i=nu+1,istack
          inext=mao(i)
          lr(inode)=lr(inext)
          li(lr(inode))=inode
          inode=inext
        end do
        lr(inode)=nextr
        li(nextr)=inode
!       write(nout,*)'lr =',(lr(j),j=nu+1,n)
!       write(nout,*)'look-ahead lengths =',(li(lc(j)),j=nu+1,iq)
        t(iq)=-1
        if(backtrack .or. istack>nu+1) then
          do i=nu+1,iq-1
            t(i)=-1
          end do
        end if
        do i=1,n
          if(li(i)>n) then
            write(nout,*)'iq =',iq
            stop
          end if
        end do
      end do
!     write(nout,*)'transversal found'
!     write(nout,*)'lr =',(lr(j),j=1,n)
!     write(nout,*)'lc =',(lc(j),j=nu+1,n)
!  print star diagram
!     if(n-nu>80 .or. n>1000) stop
!     write(nout,*)'transversal ordering'
!     do i=1,n
!       do j=1,n-nu
!         star(i,j)=' '
!       end do
!     end do
!     do j=1,n-nu
!       ilp=lc(nu+j)-n
!       do i=point(ilp),point(ilp+1)-1
!         star(li(la(i)),j)='*'
!       end do
!     end do
!     do i=nu+1,n
!       write(nout,*)(star(i,j),j=1,n-nu)
!     end do

!  tarjan ordering
      do i=1,n
        q(i)=0
        r(i)=0
      end do
!  reset li and pair off columns with rows
      do i=nu+1,n
        nodec=lc(i)
        li(nodec)=i
        t(lr(i))=nodec
        s(i)=0
      end do
      do i=nu+1,n
        noder=lr(i)
        nodec=t(noder)
        lc(noder)=point(nodec-n+1)-point(nodec-n)
        li(nodec)=-1
      end do
      ifath=nu
      istack=n+1
!  tarjan loop
10    continue
        istack=istack-1
        inode=istack
        noder=lr(inode)
        if(lc(noder).eq.0) then
          write(nout,*)'malfunction: zero length'
          stop
        end if
        nodec=t(noder)
11      continue
        li(nodec)=lc(noder)
        mao(inode)=istack
!       write(nout,*)'put new node',noder,' on stack'
!       write(nout,*)'active part of lr =',(lr(j),j=ifath+1,n)
!       write(nout,*)'ifath,istack =',ifath,istack
!       write(nout,*)'column node =',nodec,'  unfathomed rows =',
!    *    (la(j),j=point(nodec-n),point(nodec-n)+li(nodec)-1)
12      continue
          if(li(nodec).eq.0) then
!           write(nout,*)'backtrack to previous nodes'
13          continue
              if(inode.eq.n) goto 14
              inext=inode+1
              nextr=lr(inext)
              if(mao(inode)<mao(inext)) goto 14
              inode=inext
              noder=nextr
              nodec=t(noder)
              if(li(nodec).eq.0) goto 13
!           write(nout,*)'stack =',(lr(j),j=istack,n)
!           write(nout,*)'lengths =',(li(t(lr(j))),j=istack,n)
!           write(nout,*)'column node =',nodec,'  unfathomed rows =',
!    *        (la(j),j=point(nodec-n),point(nodec-n)+li(nodec)-1)
             goto 12
          end if
!  examine successors of current node
          li(nodec)=li(nodec)-1
          nextr=la(point(nodec-n)+li(nodec))
          inext=li(nextr)
          if(inext.le.ifath) goto 12
          q(nextr)=q(nextr)+1
          nextc=t(nextr)
!         write(nout,*)'nextc,nextr,inext',nextc,nextr,inext
          if(li(nextc)>=0) then
            mx=mao(inext)
            if(mao(inode)>=mx) goto 12
            do j=istack,n
              if(mao(j).eq.mx) goto 12
              mao(j)=mx
            end do
            write(nout,*)'malfunction'
            stop
          end if
          nodec=nextc
          noder=nextr
          istack=istack-1
          inode=istack
          lr(inext)=lr(inode)
          li(lr(inext))=inext
          lr(inode)=noder
          li(noder)=inode
           goto 11
14      continue
!       write(nout,*)'strong component identified'
!       write(nout,*)'active part of lr =',(lr(j),j=ifath+1,n)
!       write(nout,*)'ifath,istack,inode =',ifath,istack,inode,n
!  shift forward strong component
        inext=istack-1
        ir=inode-inext
        do j=istack,inode
          mao(j)=lr(j)
        end do
        do j=inext+ir,ifath+1+ir,-1
          lr(j)=lr(j-ir)
          li(lr(j))=j
        end do
        mx=ifath+ir
        iq=inext-ifath
        ifath=ifath+1
        do j=ifath,mx
          lr(j)=mao(j+iq)
          li(lr(j))=j
          mao(j)=mx
        end do
        istack=inode+1
        ifath=mx
!       write(nout,*)'active part of lr =',(lr(j),j=ifath+1,n)
!       write(nout,*)'ifath,istack =',ifath,istack
        if(istack.le.n) then
          inode=istack
          noder=lr(inode)
          nodec=t(noder)
          nodec_n=nodec-n
!         write(nout,*)'column node =',nodec,'  unfathomed rows =',
!    *      (la(j),j=point(nodec-n),point(nodec-n)+li(nodec)-1)
           goto 12
        end if
      if(ifath<n) goto 10
!  end of tarjan process
!  reset lc and li
      do i=nu+1,n
        lc(i)=t(lr(i))
        li(lc(i))=i
      end do
!     write(nout,*)'mao =',(mao(j),j=nu+1,n)
!     write(nout,*)'q =',(q(j),j=1,n)
!     write(nout,*)'lr =',(lr(j),j=1,n)
!     write(nout,*)'lc =',(lc(j),j=nu+1,n)
!     write(nout,*)'li =',(li(j),j=1,n+m)
!  print star diagram
!     if(n-nu>80 .or. n>1000) stop
!     write(nout,*)'tarjan ordering'
!     do i=1,n
!       do j=1,n-nu
!         star(i,j)=' '
!       end do
!     end do
!     do j=1,n-nu
!       ilp=lc(nu+j)-n
!       do i=point(ilp),point(ilp+1)-1
!         star(li(la(i)),j)='*'
!       end do
!     end do
!     do i=nu+1,n
!       write(nout,*)(star(i,j),j=1,n-nu)
!     end do
!  set up pointers for row-wise sparse structure
      p(1)=1
      do i=1,n-1
        p(i+1)=p(i)+q(i)
        q(i)=p(i)-1
      end do
      if(p(n)+q(n)>mxws) then
        ifail=7
        return
      end if
      q(n)=p(n)-1
      i=nu+1
20    continue
      if(i.eq.mao(i)) then
        t(i)=i
      else
!  spk1 ordering on tarjan block
!  set row and column counts
        do inode=i,mao(i)
          nodec=lc(inode)
          do j=point(nodec-n),point(nodec-n+1)-1
            noder=la(j)
            if(li(noder)>=i) then
              q(noder)=q(noder)+1
              ws(q(noder))=dble(nodec)
              s(inode)=s(inode)+1
            end if
          end do
        end do
!       print *,'r-c counts: i =',i,'   mao(i) =',mao(i)
!       print *,'q =',(q(j),j=i,mao(i))
!       print *,'s =',(s(j),j=i,mao(i))
!  find minimum-column-count column
        mcc=n
        do inode=i,mao(i)
          noder=lr(inode)
          r(noder)=q(noder)-p(noder)+1
          mcc=min(mcc,s(inode))
        end do
!     write(nout,*)'i,mao(i),mcc',i,mao(i),mcc
!     write(nout,*)'p =',(p(lr(j)),j=i,mao(i))
!     write(nout,*)'q =',(q(lr(j)),j=i,mao(i))
!     write(nout,*)'r =',(r(lr(j)),j=i,mao(i))
!     write(nout,*)'s =',(s(j),j=i,mao(i))
!  check for fully dense block
        if(mcc>mao(i)-i) then
          do inode=i,mao(i)
            t(inode)=mao(i)
          end do
           goto 22
        end if
!  determine spk1 ordering
        ifirstr=i
        ifirstc=i
21      continue
!  apply tie-break rule
        tie=0
        do inode=ifirstc,mao(i)
          if(s(inode).eq.mcc) then
            nodec=lc(inode)-n
            ti=0
            do j=point(nodec),point(nodec+1)-1
              noder=la(j)
              if(li(noder)>=ifirstr)ti=ti+r(noder)
            end do
            if(ti>tie) then
              tie=ti
              mccc=nodec
            end if
          end if
        end do
!       write(nout,*)'tie,mccc',tie,mccc+n
!  permute rows of m-c-c column to top and update column counts
        do j=point(mccc),point(mccc+1)-1
          noder=la(j)
          ir=li(noder)
          if(ir>=ifirstr) then
            lr(ir)=lr(ifirstr)
            li(lr(ir))=ir
            lr(ifirstr)=noder
            li(noder)=ifirstr
            ifirstr=ifirstr+1
            do ir=p(noder),q(noder)
              inode=li(int(ws(ir)))
              s(inode)=s(inode)-1
            end do
          end if
        end do
!       write(nout,*)'s =',(s(ij),ij=i,mao(i))
!       write(nout,*)'lr =',(lr(ij),ij=i,mao(i))
!  move zero-column-count columns to lhs and find minimum column count
        mcc=n
        do inode=ifirstc,mao(i)
          if(s(inode).eq.0) then
            nodec=lc(inode)
            lc(inode)=lc(ifirstc)
            li(lc(inode))=inode
            lc(ifirstc)=nodec
            li(nodec)=ifirstc
            s(inode)=s(ifirstc)
            t(ifirstc)=ifirstr-1
            ifirstc=ifirstc+1
          else
            mcc=min(mcc,s(inode))
          end if
        end do
!       write(nout,*)'lc =',(lc(ij),ij=i,mao(i))
!       write(nout,*)'ifirstc,mcc',ifirstc,mcc
        if(ifirstc<mao(i)) goto 21
      end if
22    continue
      i=mao(i)+1
      if(i.le.n) goto 20
!  print star diagram
!     if(n-nu>80 .or. n>1000) stop
!     write(nout,*)'tarjan + spk1 ordering'
!     do i=1,n
!       do j=1,n-nu
!         star(i,j)=' '
!       end do
!     end do
!     do j=1,n-nu
!       ilp=lc(nu+j)-n
!       do i=point(ilp),point(ilp+1)-1
!         star(li(la(i)),j)='*'
!       end do
!     end do
!     do i=nu+1,n
!       write(nout,*)(star(i,j),j=1,n-nu)
!     end do
!     write(nout,*)'lr =',(lr(j),j=nu+1,n)
!     write(nout,*)'lc =',(lc(j),j=nu+1,n)
!     write(nout,*)'lower profile =',(t(j),j=nu+1,n)
      ifail=0
      return
      end

      subroutine re_factor(n,nm,a,la,lr,lc,li,mao,p,q,r,s, &
        t,ws,mxws,d,ifail)
      implicit double precision (a-h,u-z), integer (i-t)
      dimension a(*),la(0:*),lr(*),lc(*),li(*),mao(*), &
        p(*),q(*),r(*),s(*),t(*),d(*),ws(*)
!     character star(1000,80)
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/iprintc/iprint
      common/refactorc/nup,nfreq
      common/epsc/eps,tol,emin
      common/noutc/nout
      double precision thresh,tol
      parameter (thresh=1.D-1)
!  factorize LPA=U
!    p(row) stores the number of stored elements of a natural row
!    q(row) stores the base address in ws of a natural row
!    r(row) stores the previous row stored in ws (or 0 if the first row in ws)
!    s(row) stores the next row stored in ws (or 0 if the last row in ws)
!    t(*) stores the lower profile of the sparse matrix
!    irow stores the natural row number of the initial row stored in ws
!    lastr stores the natural row number of the previous row put into ws
!     write(nout,*)'re_factor'
      nup=0
      m=nm-n
      lastr=0
      irow=0
      do i=1,n
        p(i)=0
      end do
      if(m1.eq.0) return
      i=nu+1
1     continue
      if(i.eq.mao(i)) then
        d(i)=aij(lr(i),lc(i)-n,a,la)
        if(d(i).eq.0.D0)d(i)=eps
!       write(nout,*)'row,col,d(i) =',lr(i),lc(i),d(i)
      else
!       write(nout,*)'lc =',(lc(j),j=i,mao(i))
        do inode=i,mao(i)-1
          nodec=lc(inode)-n
          im=inode-1
!  form L.a_q
          z=0.
!         write(nout,*)'inode,t(inode)',inode,t(inode)
          do j=inode,t(inode)
            rowj=lr(j)
            prj=p(rowj)
            if(prj>0) then
              d(j)=aiscpri2(n,a,la,rowj,nodec,ws(q(rowj)+1),1.D0,im, &
                prj,li)
            else
              d(j)=aij(rowj,nodec,a,la)
            end if
            z=max(z,abs(d(j)))
          end do
!         write(nout,*)'d =',(d(ij),ij=inode,t(inode))
!  threshold pivot selection
          zz=z*thresh
          z=0.D0
          pri=n
          do j=inode,t(inode)
            dj=abs(d(j))
            if(dj>=zz) then
              prj=p(lr(j))
              if(prj.eq.pri) then
                if(dj>z) then
                  z=dj
                  iz=j
                end if
              else if(prj<pri) then
                z=dj
                iz=j
                pri=prj
              end if
            end if
          end do
!       write(nout,*)'zz,z,iz,pri',zz,z,iz,pri
          if(iz>inode) then
!  pivot interchange
            call rexch(d(inode),d(iz))
            call iexch(lr(inode),lr(iz))
            li(lr(iz))=iz
            li(lr(inode))=inode
          end if
          if(d(inode).eq.0.D0)d(inode)=eps
!  update L
          qri=q(lr(inode))
          zz=-d(inode)
          do j=inode+1,t(inode)
            z=d(j)/zz
            rowj=lr(j)
            prj=p(rowj)
            qrj=q(rowj)
!  find space available in-situ in ws
            if(prj.eq.0) then
              len=0
            else if(s(rowj).eq.0) then
              len=mxws-qrj
            else
              len=q(s(rowj))-qrj
            end if
            if(abs(z).le.tol) then
!  special case of a zero multiplier
              if(prj.eq.0) goto 2
              len_=prj+1
              if(len_>len) then
                call newslot(rowj,len_,lastr,irow,p,q,r,s,ws,mxws,qrj, &
                  ifail)
                if(ifail>0) return
                qrj_=q(rowj)
                do k=1,prj
                  ws(qrj_+k)=ws(qrj+k)
                end do
                ws(qrj_+len_)=z
              else
                ws(qrj+len_)=z
              end if
              p(rowj)=len_
               goto 2
            end if
            len_=max(pri,prj)+1
            if(len_>len .or. pri>prj) then
!  create a new slot and use saxpyz ...
              call newslot(rowj,len_,lastr,irow,p,q,r,s,ws,mxws,qrj, &
                ifail)
              if(ifail>0) return
              qrj_=q(rowj)
              len=prj-pri
              if(len>=0) then
                do k=1,len
                  ws(qrj_+k)=ws(qrj+k)
                end do
                len=len+1
                call saxpyz(z,ws(qri+1),ws(qrj+len),ws(qrj_+len), &
                  len_-len)
              else
                len=-len
                do k=1,len
                  ws(qrj_+k)=z*ws(qri+k)
                end do
                len=len+1
                call saxpyz(z,ws(qri+len),ws(qrj+1),ws(qrj_+len), &
                  len_-len)
              end if
              ws(qrj_+len_)=z
            else
!  ... else saxpy in-situ
              if(pri>0) &
                call mysaxpy(z,ws(qri+1),ws(qrj+prj-pri+1),pri)
              ws(qrj+len_)=z
            end if
            p(rowj)=len_
!           do rj=1,n
!             if(p(rj)/=0) then
!               write(nout,*)'storage for row',rj,'  p,q,r,s =',
!    *            p(rj),q(rj),r(rj),s(rj)
!             end if
!           end do
2           continue
          end do
!         write(nout,*)'lr =',(lr(j),j=i,mao(i))
!         do j=i,mao(i)
!           rowj=lr(j)
!           if(p(rowj)/=0) then
!             write(nout,*)'L(',rowj,')',
!    *          (ws(k),k=q(rowj)+1,q(rowj)+p(rowj))
!           end if
!         end do
        end do
        inode=mao(i)
        noder=lr(inode)
        pri=p(noder)
        if(pri>0) then
         d(inode)=aiscpri2(n,a,la,noder,lc(inode)-n,ws(q(noder)+1), &
           1.D0,inode-1,pri,li)
        else
          d(inode)=aij(noder,lc(inode)-n,a,la)
        end if
        if(d(inode).eq.0.D0)d(inode)=eps
      end if
      i=mao(i)+1
      if(i.le.n) goto 1
!     write(nout,*)'PAQ factors:  nu =',nu
!     write(nout,*)'column perm =',(lc(j),j=nu+1,n)
!     write(nout,*)'row perm =',(lr(j),j=nu+1,n)
!     write(nout,*)'d =',(d(ij),ij=nu+1,n)
!     do j=nu+1,n
!       rowj=lr(j)
!       if(p(rowj)/=0) then
!         write(nout,*)'L(',rowj,')',
!    *      (ws(k),k=q(rowj)+1,q(rowj)+p(rowj))
!       end if
!     end do
!     call checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
!  print star diagram
!     if(m1>80 .or. n>1000) stop
!     write(nout,*)'factored tarjan + spk1 ordering:  nu =',nu
!     do i=1,n
!       do j=1,m1
!         star(i,j)=' '
!       end do
!     end do
!     do j=1,m1
!       jp=la(0)+lc(nu+j)-n
!       do i=la(jp),la(jp+1)-1
!         star(li(la(i)),j)='*'
!       end do
!     end do
!     do i=nu+1,n
!       write(nout,*)(star(i,j),j=1,m1)
!     end do
!     write(nout,*)'lr =',(lr(j),j=nu+1,n)
!     write(nout,*)'lc =',(lc(j),j=nu+1,n)
      mp=-1
      mq=-1
      ifail=0
      return
      end

      function aiscpri2(n,a,la,rowi,coli,ws,di,im,pri,li)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),ws(*),li(*)
      integer rowi,coli,rowj,pri
      aiscpri2=0.D0
      jp=la(0)+coli
      do j=la(jp),la(jp+1)-1
        rowj=la(j)
        if(rowj.eq.rowi) then
          aiscpri2=aiscpri2+di*a(j)
        else
          ir=li(rowj)-im
          if(ir>0) goto 1
          ir=ir+pri
          if(ir>0)aiscpri2=aiscpri2+ws(ir)*a(j)
        end if
1       continue
      end do
      return
      end

      subroutine update_L(pp,qq,n,nm,a,la,lr,lc,li,mao,p,q,r,s, &
        ws,mxws,d,sn,ifail)
      implicit double precision (a-h,r-z), integer (i-q)
      dimension a(*),la(0:*),lr(*),lc(*),li(*),mao(*), &
        p(*),q(*),r(*),s(*),ws(*),d(*),sn(*)
!     character star(1000,80)
      double precision l11,l21
      integer r,s,rowim,rowi,rowj,rrj
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/refactorc/nup,nfreq
      common/iprintc/iprint
      common/epsc/eps,tol,emin
      common/noutc/nout
      parameter (thresh=1.D-1,growth=1.D1)
!     write(nout,*)'update_L:  p,q =',pp,qq
      nup=nup+1
      if(qq>n) then
        ilast=nu
        jp=la(0)+qq-n
        do j=la(jp),la(jp+1)-1
          ip=li(la(j))
          if(ip>nu)ilast=max(ilast,mao(ip))
        end do
        qqq=qq
      else
!  row flma procedure to remove row qq (includes qq amongst the unit vectors)
        iq=li(qq)
        if(iq.le.nu) goto 99
        ilast=mao(iq)
        l11=1.D0
        u11=d(iq)
        ss=-sn(iq)
        nu=nu+1
        do i=iq,nu+1,-1
          lr(i)=lr(i-1)
          li(lr(i))=i
          sn(i)=sn(i-1)
          d(i)=d(i-1)
        end do
        lr(nu)=qq
        li(qq)=nu
!  update mao
        do j=iq-1,nu,-1
          if(mao(j)<ilast) goto 5
        end do
        j=nu-1
5       continue
        do j=j,nu,-1
          mao(j+1)=mao(j)+1
        end do
        prq=p(qq)
        if(prq>0)qrq=q(qq)
        do i=iq+1,ilast
          im=i-1
          rowi=lr(i)
          pri=p(rowi)
          u22=d(i)
          if(prq>0) then
            u12=aiscpri2(n,a,la,qq,lc(i)-n,ws(qrq+1),l11,im,prq,li)
          else
            u12=l11*aij(qq,lc(i)-n,a,la)
          end if
          if(abs(u12).le.tol)u12=0.D0
          if(pri>0) then
            qri=q(rowi)
            is=im-iq
            ii=pri-is
            if(ii.le.0) then
              l21=0.
            else
              l21=ws(qri+ii)
              if(abs(l21).le.tol)l21=0.D0
              if(ii.eq.1) then
                call trim_(rowi,pri,qri,q,ws)
                if(pri.eq.0) call erase(rowi,lastr,irow,r,s)
                if(s(rowi).eq.0) then
                  qr_=mxws
                else
                  qr_=q(s(rowi))
                end if
                if(qri+pri>=qr_) then
                  call r_shift(ws(qri),pri,1)
                  qri=qri-1
                  q(rowi)=qri
                end if
              else
                pri=pri-1
                call r_shift(ws(qri+ii),is,1)
              end if
            end if
            p(rowi)=pri
          else
            l21=0.D0
          end if
          rr=-l21/l11
          del=rr*u12+u22
          test=abs(rr)*max(abs(u11),abs(u22))
!         write(nout,*)'l11,l21,u11,u12,u22,del,test',
!    *      l11,l21,u11,u12,u22,del,test
          is=pri-prq
          if(is<0)test=test*growth
          if(u12.eq.0.D0 .and. is>0)test=test*thresh
!           write(nout,*)'rowi,pri,qri =',rowi,pri,qri
!           write(nout,*)'rowq,prq,qrq =',qq,prq,qrq
!           write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
!           do j=1,n
!             if(p(j)/=0) write(nout,*)j,p(j),q(j),r(j),s(j)
!           end do
!           write(nout,*)'rowq =',(ws(qrq+ij),ij=1,prq)
!           write(nout,*)'rowi =',(ws(qri+ij),ij=1,pri)
          if(abs(del).le.test) then
!  no-perm operation for row flma
!           write(nout,*)'no-perm operation for row flma'
            if(is>0) then
              pr_=prq
              prq=pri+1
              call newslot(qq,prq,lastr,irow,p,q,r,s,ws,mxws,qr_,ifail)
              if(ifail>0) return
              qrq=q(qq)
              qri=q(rowi)
              call r_shift(ws(qrq+1),pri,qri-qrq)
              call mysaxpy(rr,ws(qr_+1),ws(qri+is+1),pr_)
            else
              if(prq.eq.0) then
                call erase(rowi,lastr,irow,r,s)
                p(rowi)=0
                call newslot(qq,1,lastr,irow,p,q,r,s,ws,mxws,qr_,ifail)
                if(ifail>0) return
                prq=1
                qrq=q(qq)
              else
                is=-is
                do j=1,is
                  ws(qrq+j)=rr*ws(qrq+j)
                end do
                if(pri>0) then
                  call saxpyx(rr,ws(qrq+is+1),ws(qri+1),pri)
                else
                  call newslot(rowi,1,lastr,irow,p,q,r,s,ws,mxws,qr_, &
                    ifail)
                  if(ifail>0) return
                  qri=q(rowi)
                  qrq=q(qq)
                end if
                if(abs(ws(qrq+1)).le.tol) call trim_(qq,prq,qrq,q,ws)
!  rename qq as rowi and vice-versa
                if(qri<qrq) then
                  if(s(rowi).eq.qq) then
                    r(qq)=r(rowi)
                    r(rowi)=qq
                    s(rowi)=s(qq)
                    s(qq)=rowi
                  else
                    call iexch(r(qq),r(rowi))
                    call iexch(s(qq),s(rowi))
                    r(s(qq))=qq
                    s(r(rowi))=rowi
                  end if
                  if(r(qq)>0) then
                    s(r(qq))=qq
                  else
                    irow=qq
                  end if
                  if(s(rowi)>0)r(s(rowi))=rowi
                else
                  if(s(qq).eq.rowi) then
                    r(rowi)=r(qq)
                    r(qq)=rowi
                    s(qq)=s(rowi)
                    s(rowi)=qq
                  else
                    call iexch(r(rowi),r(qq))
                    call iexch(s(rowi),s(qq))
                    r(s(rowi))=rowi
                    s(r(qq))=qq
                  end if
                  if(r(rowi)>0) then
                    s(r(rowi))=rowi
                  else
                    irow=rowi
                  end if
                  if(s(qq)>0)r(s(qq))=qq
                end if
                call iexch(pri,prq)
                call iexch(qri,qrq)
                call iexch(q(rowi),q(qq))
                if(pri.eq.0) call erase(rowi,lastr,irow,r,s)
                prq=prq+1
              end if
            end if
            p(rowi)=pri
            p(qq)=prq
            ws(qrq+prq)=1.D0
            d(i)=rr*u11
            u11=u22
            l11=l21
          else
!  perm operation for row flma
!           write(nout,*)'perm operation for row flma'
            if(rr/=0.D0) then
              if(is>=0) then
                if(prq>0) then
                  call mysaxpy(rr,ws(qrq+1),ws(qri+is+1),prq)
                  if(abs(ws(qri+1)).le.tol) call trim_(rowi,pri,qri,q,ws)
                  if(pri.eq.0) call erase(rowi,lastr,irow,r,s)
                end if
                is=pri-prq
              else
                pr_=pri
                pri=prq
                call newslot(rowi,pri,lastr,irow,p,q,r,s,ws,mxws,qr_, &
                  ifail)
                if(ifail>0) return
                qrq=q(qq)
                qri=q(rowi)
                is=-is
                do j=1,is
                  ws(qri+j)=rr*ws(qrq+j)
                end do
                call saxpyz(rr,ws(qrq+is+1),ws(qr_+1),ws(qri+is+1),pr_)
                is=0
              end if
            end if
            p(rowi)=pri
            if(u12/=0.D0) then
              u12=-u12/del
              if(is>0) then
                pr_=prq
                prq=pri+1
                call newslot(qq,prq,lastr,irow,p,q,r,s,ws,mxws,qr_, &
                  ifail)
                if(ifail>0) return
                qrq=q(qq)
                qri=q(rowi)
                do j=1,is
                  ws(qrq+j)=u12*ws(qri+j)
                end do
                call saxpyz(u12,ws(qri+is+1),ws(qr_+1),ws(qrq+is+1),pr_)
                ws(qrq+prq)=u12
                 goto 7
              else
                if(pri>0) then
                  is=-is
                  call mysaxpy(u12,ws(qri+1),ws(qrq+is+1),pri)
                  if(abs(ws(qrq+1)).le.tol) then
                    call trim_(qq,prq,qrq,q,ws)
                    if(prq.eq.0) call erase(qq,lastr,irow,r,s)
                    p(qq)=prq
                  end if
                end if
              end if
            end if
            if(prq>0 .or. u12/=0.D0) then
              if(prq.eq.0) then
                len=0
              else if(s(qq).eq.0) then
                len=mxws-qrq
              else
                len=q(s(qq))-qrq
              end if
              if(len.eq.prq) then
                call newslot(qq,prq+1,lastr,irow,p,q,r,s,ws,mxws,qr_, &
                  ifail)
                if(ifail>0) return
                qrq=q(qq)
                qri=q(rowi)
                call r_shift(ws(qrq+1),prq,qr_-qrq)
              end if
              prq=prq+1
              ws(qrq+prq)=u12
            end if
7           continue
            p(rowi)=pri
            p(qq)=prq
            d(i)=del
            u11=u11*u22/del
            call iexch(lc(i),lc(im))
          end if
!           write(nout,*)'rowi,pri,qri =',rowi,pri,qri
!           write(nout,*)'rowq,prq,qrq =',qq,prq,qrq
!           write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
!           do j=1,n
!             if(p(j)/=0) write(nout,*)j,p(j),q(j),r(j),s(j)
!           end do
!           write(nout,*)'rowq* =',(ws(qrq+ij),ij=1,prq)
!           write(nout,*)'rowi* =',(ws(qri+ij),ij=1,pri)
        end do
        if(prq>0) then
!         write(nout,*)'ss,l11,ilast,n,prq',ss,l11,ilast,n,prq
!         write(nout,*)'sn =',(sn(ij),ij=nu+1,n)
          call mysaxpy(ss/l11,ws(qrq+1),sn(ilast-prq+1),prq)
          call erase(qq,lastr,irow,r,s)
          p(qq)=0
        end if
        qqq=lc(ilast)
        do i=ilast,nu+1,-1
          lc(i)=lc(i-1)
          li(lc(i))=i
        end do
!       if(pp.le.n) then
!         ip=li(pp)
!         write(nout,*)'check sn'
!         do i=nu+1,ilast
!           nodec=lc(i)
!           u12=aiscpri2(n,a,la,pp,lc(i)-n,sn(nu+1),1.D0,ilast,
!             ilast-nu,li)
!           if(abs(u12)>tol) write(nout,*)'error,nodec =',u12,nodec
!         end do
!       end if
!       write(nout,*)'intermediate PAQ factors:  new q =',qqq
!       write(nout,*)'lr =',(lr(j),j=nu+1,n)
!       write(nout,*)'lc =',(lc(j),j=nu+1,n)
!       write(nout,*)'d =',(d(ij),ij=nu+1,n)
!       do j=nu+1,n
!         rowj=lr(j)
!         if(p(rowj)/=0) then
!           write(nout,*)'L(',rowj,')',
!    *        (ws(k),k=q(rowj)+1,q(rowj)+p(rowj))
!         end if
!       end do
!       call checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
      end if
      ip=li(pp)
      if(pp>n) then
        li(pp)=0
        if(pp.eq.qqq) goto 30
        if(ip.le.nu) goto 99
        iout=ip
        rowim=lr(ip)
        prim=p(rowim)
        if(prim>0)qrim=q(rowim)
      else
        if(ip>nu .or. p(pp)>0) goto 99
        lr(ip)=lr(nu)
        li(lr(ip))=ip
!  check for growth in sn
!       write(nout,*)'sn =',(sn(i),i=nu+1,n)
        iout=ilast
        i=nu+1
        if(i>ilast) goto 13
11      continue
          do j=i,mao(i)
            if(abs(sn(j))>growth) then
              iout=i-1
               goto 13
            end if
          end do
          i=mao(i)+1
          if(i.le.ilast) goto 11
13      continue
        do j=nu+1,iout
          if(abs(sn(j))>tol) goto 14
        end do
        j=iout+1
14      continue
        rowim=pp
        prim=iout-j+1
        if(prim>0) then
          call newslot(pp,prim,lastr,irow,p,q,r,s,ws,mxws,qr_,ifail)
          if(ifail>0) return
          p(pp)=prim
          qrim=q(pp)
          ii=qrim
          do j=j,iout
            ii=ii+1
            ws(ii)=sn(j)
          end do
        end if
        do i=nu,iout-1
          lr(i)=lr(i+1)
          li(lr(i))=i
          lc(i)=lc(i+1)
          li(lc(i))=i
          d(i)=d(i+1)
        end do
        lr(iout)=pp
        li(pp)=iout
!       write(nout,*)'lr =',(lr(ij),ij=nu,iout)
!       write(nout,*)'lc =',(lc(ij),ij=nu,iout-1)
!       if(prim>0) write(nout,*)'L(',pp,') =',(ws(qrim+j),j=1,prim)
        nu=nu-1
      end if
!     write(nout,*)'iout,ilast,rowim,prim =',iout,ilast,rowim,prim
!  column flma operations to restore L to triangular form
      iswap=0
      do i=iout+1,ilast
        im=i-1
        lc(im)=lc(i)
        li(lc(im))=im
        rowi=lr(i)
        pri=p(rowi)
!       if(pri>0) write(nout,*)'L(',rowi,') =',(ws(q(rowi)+j),j=1,pri)
        u22=d(i)
        if(prim>0) then
          u12=aiscpri2(n,a,la,rowim,lc(i)-n,ws(qrim+1),1.D0,im-1,prim, &
            li)
          if(abs(u12).le.tol)u12=0.D0
        else
          u12=aij(rowim,lc(i)-n,a,la)
        end if
        if(pri>0) then
!         write(nout,*)'pri,iswap',pri,iswap
          qri=q(rowi)
          ii=pri-iswap
          if(ii.le.0) then
            l21=0.D0
          else
            l21=ws(qri+ii)
            if(abs(l21).le.tol)l21=0.D0
            if(ii.eq.1) then
              call trim_(rowi,pri,qri,q,ws)
              if(pri.eq.0) call erase(rowi,lastr,irow,r,s)
              if(s(rowi).eq.0) then
                qr_=mxws
              else
                qr_=q(s(rowi))
              end if
              if(qri+pri>=qr_) then
                call r_shift(ws(qri),pri,1)
                qri=qri-1
                q(rowi)=qri
              end if
            else
              pri=pri-1
              call r_shift(ws(qri+ii),iswap,1)
            end if
            p(rowi)=pri
!           write(nout,*)'rowi =',(ws(qri+ij),ij=1,pri)
          end if
        else
          l21=0.D0
        end if
        del=u22-l21*u12
        test=abs(u12)*max(1.D0,abs(l21))
!       write(nout,*)'l21,u12,u22,del,test',l21,u12,u22,del,test
        is=pri-prim
        if(is>0)test=growth*test
        if(l21.eq.0.D0 .and. is<0)test=thresh*test
!         write(nout,*)'rowim,prim,qrim =',rowim,prim,qrim
!         write(nout,*)'rowi,pri,qri =',rowi,pri,qri
!         write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
!         do j=1,n
!           if(p(j)/=0) write(nout,*)j,p(j),q(j),r(j),s(j)
!         end do
!         write(nout,*)'rowim =',(ws(qrim+ij),ij=1,prim)
!         write(nout,*)'rowi =',(ws(qri+ij),ij=1,pri)
        if(abs(del).le.test) then
!  no-perm operation for column flma
!         write(nout,*)'no-perm operation for column flma'
          rr=-u22/u12
          l21=l21+rr
          if(abs(l21).le.tol)l21=0.D0
          if(is>=0) then
            if(prim>0) then
              call mysaxpy(rr,ws(qrim+1),ws(qri+is+1),prim)
              if(abs(ws(qri+1)).le.tol) call trim_(rowi,pri,qri,q,ws)
              if(pri.eq.0) then
                call erase(rowi,lastr,irow,r,s)
                p(rowi)=0
              end if
            end if
            if(pri>0 .or. l21/=0.D0) then
              if(pri.eq.0) then
                len=0
              else if(s(rowi).eq.0) then
                len=mxws-qri
              else
                len=q(s(rowi))-qri
              end if
              if(len.eq.pri) then
                call newslot(rowi,pri+1,lastr,irow,p,q,r,s,ws,mxws,qr_, &
                  ifail)
                if(ifail>0) return
                qrim=q(rowim)
                qri=q(rowi)
                call r_shift(ws(qri+1),pri,qr_-qri)
              end if
              pri=pri+1
              ws(qri+pri)=l21
            end if
          else
            pr_=pri
            pri=prim+1
            call newslot(rowi,pri,lastr,irow,p,q,r,s,ws,mxws,qr_,ifail) !
            if(ifail>0) return
            qrim=q(rowim)
            qri=q(rowi)
            is=-is
            do j=1,is
              ws(qri+j)=rr*ws(qrim+j)
            end do
            call saxpyz(rr,ws(qrim+is+1),ws(qr_+1),ws(qri+is+1),pr_)
            ws(qri+pri)=l21
          end if
!           write(nout,*)'rowim,prim,qrim =',rowim,prim,qrim
!           write(nout,*)'rowi,pri,qri =',rowi,pri,qri
!           write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
!           do j=1,n
!             if(p(j)/=0) write(nout,*)j,p(j),q(j),r(j),s(j)
!           end do
!           write(nout,*)'rowim* =',(ws(qrim+ij),ij=1,prim)
!           write(nout,*)'rowi* =',(ws(q(rowi)+ij),ij=1,p(rowi))
          p(rowi)=pri
          rowim=rowi
          prim=pri
          qrim=qri
          d(im)=u12
!  perform accumulated cyclic permutation in subsequent rows
          if(iswap>0) then
            do j=i+1,ilast
              rowj=lr(j)
              prj=p(rowj)
              is=prj-j+i
              if(is>0) then
                qrj=q(rowj)
                if(is>iswap) then
                  ii=is-iswap
                  l21=ws(qrj+ii)
                  call r_shift(ws(qrj+ii),iswap,1)
                  ws(qrj+is)=l21
                  if(abs(ws(qrj+1)).le.tol) call trim_(rowj,prj,qrj,q,ws)
                  if(prj.eq.0) call erase(rowj,lastr,irow,r,s)
                else
                  prj=prj+1
                  rrj=r(rowj)
                  if(rrj.eq.0) then
                    len=qrj
                  else
                    len=qrj-q(rrj)-p(rrj)
                  end if
                  if(len>0) then
                    call r_shift(ws(qrj),is,1)
                    ws(qrj+is)=0.D0
                    qrj=qrj-1
                    q(rowj)=qrj
                  else
                    call newslot(rowj,prj,lastr,irow,p,q,r,s,ws,mxws, &
                      qr_,ifail)
                    if(ifail>0) return
                    qrj=q(rowj)
                    qrim=q(rowim)
                    call r_shift(ws(qrj+1),is,qr_-qrj)
                    ws(qrj+is+1)=0.D0
                    call r_shift(ws(qrj+is+2),j-i,qr_-qrj-1)
                  end if
                end if
                p(rowj)=prj
!               write(nout,*)'L(',rowj,')* =',(ws(qrj+ij),ij=1,prj)
              end if
            end do
          end if
          iswap=0
        else
!  perm operation for column flma
!         write(nout,*)'perm operation for column flma'
          rr=-l21
          if(rr/=0.D0) then
            if(is>=0) then
              if(prim>0) then
                call mysaxpy(rr,ws(qrim+1),ws(qri+is+1),prim)
                if(abs(ws(qri+1)).le.tol) call trim_(rowi,pri,qri,q,ws)
                if(pri.eq.0) call erase(rowi,lastr,irow,r,s)
              end if
              is=pri-prim
            else
              pr_=pri
              pri=prim
              call newslot(rowi,pri,lastr,irow,p,q,r,s,ws,mxws,qr_, &
                ifail)
              if(ifail>0) return
              qrim=q(rowim)
              qri=q(rowi)
              is=-is
              do j=1,is
                ws(qri+j)=rr*ws(qrim+j)
              end do
              call saxpyz(rr,ws(qrim+is+1),ws(qr_+1),ws(qri+is+1),pr_)
              is=0
            end if
          end if
          p(rowi)=pri
          if(u12/=0.D0) then
            u12=-u12/del
            if(is>0) then
              pr_=prim
              prim=pri+1
              call newslot(rowim,prim,lastr,irow,p,q,r,s,ws,mxws,qr_, &
                ifail)
              if(ifail>0) return
              qrim=q(rowim)
              qri=q(rowi)
              do j=1,is
                ws(qrim+j)=u12*ws(qri+j)
              end do
              call saxpyz(u12,ws(qri+is+1),ws(qr_+1),ws(qrim+is+1),pr_)
              ws(qrim+prim)=u12
               goto 27
            else
              if(pri>0) then
                is=-is
                call mysaxpy(u12,ws(qri+1),ws(qrim+is+1),pri)
                if(abs(ws(qrim+1)).le.tol) then
                  call trim_(rowim,prim,qrim,q,ws)
                  if(prim.eq.0) call erase(rowim,lastr,irow,r,s)
                  p(rowim)=prim
                end if
              end if
            end if
          end if
          if(prim>0 .or. u12/=0.D0) then
            if(prim.eq.0) then
              len=0
            else if(s(rowim).eq.0) then
              len=mxws-qrim
            else
              len=q(s(rowim))-qrim
            end if
            if(len.eq.prim) then
              call newslot(rowim,prim+1,lastr,irow,p,q,r,s,ws,mxws,qr_, &
                ifail)
              if(ifail>0) return
              qrim=q(rowim)
              qri=q(rowi)
              call r_shift(ws(qrim+1),prim,qr_-qrim)
            end if
            prim=prim+1
            ws(qrim+prim)=u12
          end if
27        continue
          p(rowim)=prim
          p(rowi)=pri
!           write(nout,*)'rowim,prim,qrim =',rowim,prim,qrim
!           write(nout,*)'rowi,pri,qri =',rowi,pri,qri
!           write(nout,*)'j,p(j),q(j),r(j),s(j)   irow =',irow
!           do j=1,n
!             if(p(j)/=0) write(nout,*)j,p(j),q(j),r(j),s(j)
!           end do
!           write(nout,*)'rowim* =',(ws(qrim+ij),ij=1,prim)
!           write(nout,*)'rowi* =',(ws(q(rowi)+ij),ij=1,p(rowi))
          d(im)=del
          call iexch(lr(i),lr(i-1))
          call iexch(li(lr(i)),li(lr(i-1)))
          iswap=iswap+1
        end if
      end do
      lc(ilast)=qqq
      li(qqq)=ilast
!     write(nout,*)'rowim* =',(ws(qrim+ij),ij=1,prim)
!     write(nout,*)'ilast,prim,qrim',ilast,prim,qrim
      if(prim>0) then
       d(ilast)=aiscpri2(n,a,la,rowim,qqq-n,ws(qrim+1),1.D0,ilast-1, &
          prim,li)
      else
        d(ilast)=aij(rowim,qqq-n,a,la)
      end if
!  reset mao
      iout=ilast
      do i=ilast,nu+1,-1
        mao(i)=ilast
        iout=min(iout,i-p(lr(i)))
        if(iout.eq.i)ilast=i-1
      end do
30    continue
      m1=n-nu
!     write(nout,*)'PAQ factors:  nu =',nu
!     write(nout,*)'d =',(d(ij),ij=nu+1,n)
!     do j=nu+1,n
!       rowj=lr(j)
!       if(p(rowj)/=0) then
!         write(nout,*)'L(',rowj,')',
!    *      (ws(k),k=q(rowj)+1,q(rowj)+p(rowj))
!       end if
!     end do
!     call checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
!  print star diagram
!     if(m1>80 .or. n>1000) stop
!     write(nout,*)'updated ordering:  nu =',nu
!     do i=1,n
!       do j=1,m1
!         star(i,j)=' '
!       end do
!     end do
!     do j=1,m1
!       jp=la(0)+lc(nu+j)-n
!       do i=la(jp),la(jp+1)-1
!         star(li(la(i)),j)='*'
!       end do
!     end do
!     do i=nu+1,n
!       write(nout,*)(star(i,j),j=1,m1)
!     end do
!     write(nout,*)'lr =',(lr(j),j=nu+1,n)
!     write(nout,*)'lc =',(lc(j),j=nu+1,n)
!     write(nout,*)'mao =',(mao(j),j=nu+1,n)
      return
99    continue
      write(nout,*)'malfunction in update_L:  p,q =',pp,qq
      stop
      end

      subroutine newslot(row,len,lastr,irow,p,q,r,s,ws,mxws,qr_, &
        ifail)
      implicit double precision (a-h,u-z), integer (i-t)
      parameter (igap=10)
      dimension p(*),q(*),r(*),s(*),ws(*)
      common/noutc/nout
!     write(nout,*)'newslot: row =',row,'   len =',len
!     write(nout,*)'irow,lastr,mxws =',irow,lastr,mxws
      ifail=0
      if(lastr.eq.0) then
        if(mxws<len) then
          write(nout,*)'insufficient space available for profile'
          ifail=7
        else
          irow=row
          q(row)=0
          r(row)=0
          s(row)=0
          lastr=row
        end if
        return
      end if
      igp=igap
1     continue
      len_=len+igp
      thisr=lastr
2     continue
      qrow=q(thisr)+p(thisr)
      nextr=s(thisr)
!     write(nout,*)'thisr,nextr,qrow,p(thisr),len_',
!    *  thisr,nextr,qrow,p(thisr),len_
      if(nextr/=0) then
        if(q(nextr)>=qrow+len_) then
!  free slot after this row
           goto 4
        else
          thisr=nextr
          if(thisr/=lastr) goto 2
        end if
      else
        if(mxws-qrow>=len_) then
!  free slot at end of ws
           goto 4
        else if(q(irow)>=len_) then
!  free slot at beginning of ws
          qrow=0
          thisr=0
          nextr=irow
          irow=row
          igp=0
           goto 4
        end if
        thisr=irow
        if(thisr/=lastr) goto 2
      end if
!  no free space: try minimum value of len
      if(igp>0) then
        igp=0
         goto 1
      end if
!  compress ws
      thisr=irow
      qrow=0
3     continue
      call r_shift(ws(qrow+1),p(thisr),q(thisr)-qrow)
      q(thisr)=qrow
      qrow=qrow+p(thisr)
      if(s(thisr)/=0) then
        thisr=s(thisr)
         goto 3
      end if
      if(mxws<qrow+len_) then
        write(nout,*)'insufficient space available for profile'
        write(nout,*)'mxws,qrow,len_',mxws,qrow,len_
        ifail=7
        return
      end if
!  insert at end of compressed file
      nextr=0
4     continue
      qr_=q(row)
      q(row)=qrow+igp
      if(p(row)>0) then
        if(r(row).eq.thisr .or. s(row).eq.nextr) return
!  insert after row thisr and take out old row
        call erase(row,lastr,irow,r,s)
      end if
      lastr=row
      r(row)=thisr
      if(thisr>0)s(thisr)=row
      s(row)=nextr
      if(nextr>0)r(nextr)=row
      i=0
      return
      end

      subroutine erase(row,lastr,irow,r,s)
!  remove slot for row from the data file
      implicit integer (i-s)
      dimension r(*),s(*)
      common/noutc/nout
!     write(nout,*)'erase: row,irow,lastr =',row,irow,lastr
      if(r(row).eq.0) then
        if(s(row).eq.0) then
          irow=0
          lastr=0
          return
        end if
        irow=s(row)
        r(irow)=0
      else if(s(row).eq.0) then
        s(r(row))=0
      else
        s(r(row))=s(row)
        r(s(row))=r(row)
      end if
      if(row.eq.lastr)lastr=irow
      return
      end

      subroutine trim_(rowi,pri,qri,q,ws)
!  trim leading zeros off slot for row i
      implicit double precision (a-h,s-z), integer (i-r)
      dimension q(*),ws(*)
      common/epsc/eps,tol,emin
1     continue
      qri=qri+1
      pri=pri-1
      if(pri.eq.0) return
      if(abs(ws(qri+1)).le.tol) goto 1
      q(rowi)=qri
      return
      end

      subroutine checkout(n,a,la,lr,lc,li,p,q,r,s,ws,mxws,d)
      implicit double precision (a-h,r-z), integer (i-q)
      integer r,s,rowj,thisr
      dimension a(*),la(*),lr(*),lc(*),li(*),p(*),q(*),r(*),s(*),ws(*), &
        d(*)
      common/factorc/m1,nu,mp,mq,lastr,irow
      common/noutc/nout
      common/epsc/eps,tol,emin
!  check indexing
      do j=1,nu
        if(p(lr(j))/=0) then
          write(nout,*)'p(lr(j))/=0'
           goto 11
        end if
      end do
      np=0
      do i=nu+1,n
        if(p(lr(i))>0)np=np+1
      end do
      if(irow>0) then
        if(r(irow)/=0) then
          write(nout,*)'r(irow)/=0'
           goto 11
        end if
        thisr=irow
1       continue
        if(p(thisr).le.0) then
          write(nout,*)'p(thisr).le.0'
           goto 11
        end if
        np=np-1
        nextr=s(thisr)
        if(nextr.eq.0) then
          if(q(thisr)+p(thisr)>mxws) then
            write(nout,*)'q(thisr)+p(thisr)>mxws'
             goto 11
          end if
        else
          if(r(nextr)/=thisr) then
            write(nout,*)'r(nextr)/=thisr'
             goto 11
          end if
          if(nextr/=s(thisr)) then
            write(nout,*)'nextr/=s(thisr)'
             goto 11
          end if
          if(q(thisr)+p(thisr)>q(nextr)) then
            write(nout,*)'q(thisr)+p(thisr)>q(nextr)'
             goto 11
          end if
          thisr=nextr
           goto 1
        end if
      end if
      if(np/=0) then
        write(nout,*)'np/=0'
         goto 11
      end if
      last=0
      emax=0.D0
      length=0
      do inode=nu+1,n
        nodec=lc(inode)
!  form L.a_q
        rowj=lr(inode)
        prj=p(rowj)
        length=length+prj
        if(prj<0) then
          write(nout,*)'prj<0'
           goto 11
        else if(prj.eq.0) then
          e=abs(aij(rowj,nodec-n,a,la)-d(inode))
        else
          e=abs(d(inode)-aiscpri2(n,a,la,rowj,nodec-n,ws(q(rowj)+1), &
            1.D0,inode-1,prj,li))
        end if
!       if(e>tol) write(nout,*)'error =',e,
!    *    '  inode,nodec,rowj =',inode,nodec,rowj
        emax=max(emax,e)
        do j=inode+1,n
          rowj=lr(j)
          prj=p(rowj)
          if(prj>0) then
            e=abs(aiscpri2(n,a,la,rowj,nodec-n,ws(q(rowj)+1),1.D0,j-1, &
               prj,li))
          else
            e=abs(aij(rowj,nodec-n,a,la))
          end if
!         if(e>tol) write(nout,*)'error =',e,
!    *      '  inode,nodec,j,rowj =',inode,nodec,j,rowj
          emax=max(emax,e)
        end do
      end do
      write(nout,*)'checkout:  m1 =',m1,'  file length =',length
      if(emax>tol) write(nout,*)'error =',emax
      return
11    continue
      write(nout,*)'thisr,nextr =',thisr,nextr
      write(nout,*)'i,p(i),q(i),r(i),s(i):  irow =',irow
      do i=1,n
        if(p(i)/=0) write(nout,*)i,p(i),q(i),r(i),s(i)
      end do
      stop
      end
