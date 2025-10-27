PROGRAM SF
    implicit none
    integer,parameter::DP=kind(0.D0),nq=301, n=1000, nphi=250
    integer::schermo, j, i, jd, h, iq,itermu, imin, imax, s
    character(LEN=30)::file_out1, file_out2, info,file_out3, file_out5, ds, mus
    real(DP)::pi, hbar, echarge,hbar2, echarge2, die, die0, me0, me, mh, mr, mratio, conv_n, ab_effm, BE
    real(DP)::gs, gv, Ry, ab_eff, db, y1P0ne, y1P0nh, y1P0a, Vinmax, meffmax, meffmin, meff(n), meffd(n-2), dE(n-1)
    real(DP)::kc, kc0, mu_e, mu_h, mu, Dmax, preco, precon, precj, precjn, precos, precjs, muf, lp1, lp2, d, Kin, Corr
    real(DP):: k, k1, k2, k3, int_D, gap, ang_ints, ang_int, ang_intintra, En0, xi0, E(n), Emin, kmin, kmax
    real(DP)::q(nq), P0ne(nq), P0nh(nq), P0a(nq), Vq(nq), Vqint(nq), y2P0ne(nq), y2P0nh(nq), y2P0a(nq)
    real(DP)::vk2, vk2s, uk2, dens_e, dens_h, int_ne, int_nh, int_cf, ne, nh, n_tot, int_lp, lp, eps_ek, eps_hk, eps_k, En, dxi, dk
    real(DP)::rs_n, kf_n, Ef_n, kf_e, Ef_e, kf_h, Ef_h, CF, muvector(3), density, intdensity, pairlen, pairlennum, pairlenden
    real(DP):: xk(n), wk(n), D_Old(n), Delta(n), xphi(nphi), wphi(nphi), nk_Old(n), dens(n), xi(n)
    real(DP):: Pot_scr, Pot_se, Pot_sh, Pot_d, termquad, Ptot_e, Ptot_h, t1eh, t1ee, t1hh, ffeh, ffee, ffhh, intra(n)
    real(DP):: intP_ne, intP_nh, intP_a, l1, l2, qis, intraterm, intraddens, intradgap, intrabubble(n), vk(n), uk(n)
    real :: start, finish
    real(DP), allocatable :: int_intra(:,:)
    allocate(int_intra(n,n))
    
    DO s=0,0

              d=1.0d0+(s*0.90d0)   !interlayer distance in nm
              write(ds,'(F5.2)') d
   
!!!!CONSTANT!!!
    pi=acos(-1.0)
    hbar=6.6261E-34/(2.d0*pi)
    hbar2=hbar**2
    echarge=1.602E-19
    echarge2=echarge**2
    me0=9.1094E-31
    die0=8.8542E-12

!!!!PARAMETER OF THE SYSTEM!!!
    gs=2.d0                                      !SPIN degenerance
    gv=1.d0                                      !Valley degenerance
    die=2.0d0*die0                              !dielectric constant
    

!if you change it, you have to change it consistently in sub BUBBLE
    mratio=1.d0
    me=0.04d0
    mh=me*mratio
    mr=(me*mh)/(me+mh)                           !reduced mass

    ab_effm=hbar2*4.d0*pi*die/(mr*me0*echarge2)  !Bohr radius in m
    ab_eff=ab_effm*(10**9)                       !Bohr radius in nm

    conv_n=1000.d0/(ab_eff**2)                   !factor of conversion 10^11 cm^-2

    Ry=echarge/(4.d0*pi*die*2.d0*ab_effm)        !Rydberg eV
    Ry=Ry*1000.d0                                !Rydberg meV
    
    db=d/ab_eff                               !interlayer distance
    
    CALL BindingE(db, BE)

    schermo=4
    WRITE(schermo,*)'Ry(meV)=',Ry, 'ab(nm)=',ab_eff, 'factox10^11cm^-2=',conv_n
    WRITE(schermo,*)'(me,mh,mr)=',me,mh,mr, 'D*=',db,'BE(Ry)=',BE
    WRITE(schermo,*)'============================================================================'
   


    WRITE(6,*)'Ry(meV)=',Ry, 'ab(nm)=',ab_eff, 'factox10^11cm^-2=',conv_n 
    WRITE(6,*)'(me,mh,mr)=',me,mh,mr, 'D*=',db,'BE(Ry)=',BE

     call gauleg(0.d0, 2.d0*pi, xphi, wphi, nphi)        !Calculate Points and Weights for quadrature of phi in [0,2pi]

     !!!INITIAL parameters!!!
  
    mu=-0.80d0 !chemical potential in Ry

    write(mus,'(F5.1)') mu

    
    kc=30.0d0 !max value wave vector k in ab^-1
    call gauleg(0.d0, kc, xk, wk, n)
    l2=xk(n)
    

   
    D_Old=0.5d0*exp(-0.5d0*xk**2.d0)                   !Delta function trial in Ry
                      
    dens=0.0d0

    do h=1,n						!trial density of states nk
          xi0=0.5*(xk(h)**2) - mu
          En0=sqrt((xi0**2) + (D_Old(h)**2))
          nk_Old(h)= 0.5d0*(1.d0-(xi0/En0))
    enddo   
   
    intrabubble=0.0d0                          !Initial HartreeFock term
    P0ne= 0.0d0
    P0nh= 0.0d0
    P0a= 0.0d0
    y2P0ne= 0.0d0
    y2P0nh= 0.0d0
    y2P0a= 0.0d0
   
   
    file_out1='Deltakvkukd'//Trim(ds)//'YINSBECCR.txt'
    file_out2='Deltamaxd'//Trim(ds)//'YINSBECCR.txt'
    file_out3='Pold'//Trim(ds)//'YINSBECCR.txt'
    file_out5='mapd'//Trim(ds)//'YINSBECCR.txt'
    info='infod'//Trim(ds)//'YINSBECCR.txt'

    OPEN(UNIT=1, FILE=file_out1, STATUS='REPLACE', ACTION='WRITE')
    WRITE(1,"('#',T10,'k',T20,'Deltak',T28,'vk',T36,'nk',T44,'Vintra_dens',T52,'ek')")

    OPEN(UNIT=2, FILE=file_out2, STATUS='REPLACE', ACTION='WRITE')
    WRITE(2,"('#',T3,'mu',T13,'rs',T24,'Ef_e',T35,'Ef_h',T46,'n_tot',T55,'Deltamax',T72,'Vinmax',T82,'CF', & 
              T90, 'kmin', T98, 'pl', T106, 'mu_s', T114, 'K', T122, 'E_corr')")

    OPEN(UNIT=3, FILE=file_out3, STATUS='REPLACE', ACTION='WRITE')
    WRITE(3,"('#',T10,'k',T20,'Pne',T30,'Pna',T40,'V_scr',T50,'Vint_scr',T60,'V_un',T70,'Vint_un')")
    

    OPEN(UNIT=5, FILE=file_out5, STATUS='REPLACE', ACTION='WRITE')
    WRITE(5,"('#',T1,'rs',T20,'d',T30,'CF')")

    OPEN(UNIT=4, FILE=info, STATUS='REPLACE', ACTION='WRITE')
     flush(4)

z
        
       !$OMP PARALLEL DO num_threads(32) PRIVATE(j, intraterm, intradgap) &
       !$OMP SHARED(intrabubble), FIRSTPRIVATE(xphi, wphi, q, k1, k2, P0ne, P0nh)& 
       !$OMP FIRSTPRIVATE(P0a, y2P0ne, y2P0nh, y2P0a, xk, mu, wk, db, nk_Old)
 	      
            do j=1,n                                   
                k1=xk(j)
                 intradgap=0                                         !k fixed to have same mapping of k'
                do h=1,n
             		k2=xk(h) 
                    intraterm=ang_intintra(db, xphi, wphi, k1, k2, q, P0ne, P0nh, P0a, y2P0ne, y2P0nh, y2P0a)
                    intradgap=intradgap+(wk(h)*xk(h)*intraterm*nk_Old(h))
                enddo
            intrabubble(j)= intradgap
           end do
         
           !$OMP END PARALLEL DO
 
        
          !$OMP BARRIER

         
          
    DO                     !For different values of mu I calculate Delta(k)-n(k) 
        write(6,*) "chemical potential", mu
        write(1,*) "# chemical potential:", mu
        write(3,*) "# chemical potential:", mu
        call flush(4)

        write(6,*) "kc=", kc                     !(if I change it I need to change it also inside subroutine)
     
        
          
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!RECURSIVE METHOD to calculate the Gap equation !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do 
            
            !I calculate the polarizability using Delta    
            !initialize the precision to zero
            preco=0.0d0
            precon=0.0d0
 	    
 	    
            precj=0.0d0 
            precjn=0.0d0
           
         !$OMP PARALLEL DO num_threads(32) PRIVATE(j, precj, int_D, ang_ints) &
         !$OMP& REDUCTION(+:preco), SHARED(Delta), FIRSTPRIVATE(k, k1, q, intrabubble, xphi, wphi, gap)& 
         !$OMP& FIRSTPRIVATE(P0ne, P0nh, P0a, y2P0ne, y2P0nh, y2P0a, xk, wk, mu, db, D_Old, nk_Old)
 	 
         
          
            do j=1,n                                    
                k=xk(j)          !k fixed to have same mapping of k'
                int_D=0.d0
                
                do i=1,n
                    k1=xk(i)     !points k' 
                    
                   call gapEq(i, intrabubble, xk(i), xk, wk, mu, nk_Old, D_Old, gap)
           
                   ang_ints=ang_int(db, xphi, wphi, xk(j), xk(i), q, P0ne, P0nh, P0a, y2P0ne, y2P0nh, y2P0a) !integral in phi
                             
                   int_D=int_D + wk(i)*xk(i)*ang_ints*gap
                end do
                     
                Delta(j)=int_D                                  !result in the corrisponding position Delta(k)
                precj=abs(int_D-D_Old(j)) 
                
                preco=preco+precj                               !calculate precision
          end do
         
        
       !$OMP END PARALLEL DO 
       !$OMP BARRIER 
     
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! DENSITY OF STATES
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
       
         !$OMP PARALLEL DO num_threads(32) PRIVATE(jd, precjn, vk2s), REDUCTION(+:precon) &
         !$OMP& SHARED(dens), FIRSTPRIVATE(k, xk, wk, mu, db, intrabubble)& 
         !$OMP& FIRSTPRIVATE(D_Old, nk_Old)
      
          do jd=1,n             
             k=xk(jd)
             call denseq(jd, intrabubble, k, xk, wk, mu, nk_Old, D_Old, vk2s)
              
             dens(jd)=vk2s
             precjn=abs(vk2s-nk_Old(jd))  
             precon=precon+precjn
         end do    
         
         !$OMP END PARALLEL DO
    
         !$OMP BARRIER 
          

           !$OMP PARALLEL DO num_threads(32) PRIVATE(j, intraterm, intradgap) &
           !$OMP SHARED(intrabubble), FIRSTPRIVATE(h, xphi, wphi, q, k1, k2, mu, db)& 
           !$OMP FIRSTPRIVATE(P0ne, P0nh, P0a, y2P0ne, y2P0nh, y2P0a, xk, wk, dens)
 	      
            do j=1,n                                   
                k1=xk(j)
                 intradgap=0                                         !k fixed to have same mapping of k'
               do h=1,n
             		k2=xk(h) 
                    intraterm=ang_intintra(db, xphi, wphi, k1, k2, q, P0ne, P0nh, P0a, y2P0ne, y2P0nh, y2P0a)
                    intradgap=intradgap+(wk(h)*xk(h)*intraterm*dens(h))
                enddo
            intrabubble(j)= intradgap
           end do
         
           !$OMP END PARALLEL DO
 
        
            !$OMP BARRIER

       
            D_Old=Delta-((Delta-D_Old)/10)          !I rewrite Delta in DeltaOld and go on with recursive method
            nk_Old=dens-((dens-nk_Old)/10)
            intrabubble=intra-((intra-intrabubble)/10)
           
            Dmax=MAXVAL(Delta)
            Vinmax=MAXVAL(intra)
         
         
         if ((preco.le.0.0001).and.(precon.le.0.0001)) then
             
      
       !$OMP PARALLEL DO num_threads(32) PRIVATE(j, intraterm, intradgap) &
       !$OMP SHARED(intrabubble), FIRSTPRIVATE(xphi, wphi, q, k1, k2, P0ne, P0nh)& 
       !$OMP FIRSTPRIVATE(P0a, y2P0ne, y2P0nh, y2P0a, xk, mu, wk, db, nk_Old)
 	      
            do j=1,n                                   
                k1=xk(j)
                 intradgap=0                                         !k fixed to have same mapping of k'
                do h=1,n
             		k2=xk(h) 
                    intraterm=ang_intintra(db, xphi, wphi, k1, k2, q, P0ne, P0nh, P0a, y2P0ne, y2P0nh, y2P0a)
                    intradgap=intradgap+(wk(h)*xk(h)*intraterm*nk_Old(h))
                enddo
            intrabubble(j)= intradgap
           end do
         
           !$OMP END PARALLEL DO
 
        
          !$OMP BARRIER

                do iq=1,nq
                WRITE(3,"(F10.6,X,F12.6,X,F10.6,X,F10.6,X,F15.6,X,F15.6,X,F15.6,X,F15.6)") q(iq), P0ne(iq), P0a(iq), Vq(iq), &
      							                         Vqint(iq),  dexp(-db*q(iq))/q(iq) ,1/q(iq)
                end do
            
                exit
          
            else 
              
                write(6,*) preco, precon, Dmax, Vinmax
        end if
       
        end do



!!!!!!DENSITY!!!!!!
  
         int_ne=0.d0
         int_nh=0.d0
         int_cf=0.d0
            pairlennum=0.0d0 
          
            
       do jd=1,n
            int_ne= int_ne + wk(jd)*xk(jd)*dens(jd)
            int_nh= int_nh + wk(jd)*xk(jd)*dens(jd)
            int_cf= int_cf + wk(jd)*xk(jd)*dens(jd)*(1-dens(jd))       !DENSITY OF THE PAIRS
            call pairlength(xk(jd), xk, Delta,intrabubble, mu, lp1)
            pairlennum=pairlennum+ wk(jd)*xk(jd)*lp1
             WRITE(1,"(F10.6,X,F10.6,X,F10.6,X,F16.12,X,F10.6,X,F10.6,X,F10.6)") xk(jd), Delta(jd), sqrt(dens(jd)), dens(jd), &
                                                                               intra(jd), 0.5*(xk(jd)**2)
           
        end do
      
      
        ne= gs*gv*int_ne/(2.d0*pi)
        nh= gs*gv*int_nh/(2.d0*pi)
        n_tot= ne
        CF=gs*gv*int_cf/(2.d0*pi)    !DENSITY OF THE PAIR
        rs_n=1.d0/(sqrt(pi*n_tot))

        kf_n=sqrt(4.d0*pi*n_tot/(gs*gv))
        kf_e=sqrt(4.d0*pi*ne/(gs*gv))
        kf_h=sqrt(4.d0*pi*nh/(gs*gv))

        Ef_n=0.5*(kf_n**2)
        Ef_e=mr*(kf_e**2)/me
        Ef_h=mr*(kf_h**2)/mh
        
       pairlen=sqrt(pairlennum/int_cf)
!kmin
        do i=1,n
           E(i)=sqrt((((0.5*(xk(i)**2))-mu-intra(i))**2) + (Delta(i)**2))
       end do 
       
        Emin=MINVAL(E)
                        
        do i=1,n
        	if (E(i).eq.Emin) then
          	    kmin=xk(i) 
          	    imin=i  
        	exit 
        	end if  
       end do
        
!effective mass       
        do i=1,n
           meff(i)=  (0.5*(xk(i)**2))/((0.5*(xk(i)**2))-intra(i))
        end do
        
        do i=1,n-1
           dE(i)=kf_n*(E(i+1)-E(i))/(xk(i+1)-xk(i))
        end do
        
        do i=1,n-2
           meffd(i)=2*Ef_n/(kf_n*(dE(i+1)-dE(i))/(xk(i+1)-xk(i)))
        end do
        
        meffmax=MAXVAL(meff)
        meffmin=MINVAL(meff)
        
        

       

!Kinetic and correlation 
Kin=0.0d0
Corr=0.0d0

do i=1,n
   
   Kin=Kin+xk(i)*0.5*(xk(i)**2)*dens(i)
   Corr=Corr+xk(i)*intra(i)*dens(i)

end do       


       
!effective chemical potential       
       do i=1,n-1
           xi(i)=(0.5*(xk(i)**2))-mu-intra(i)
       end do
       do i=1,n-1
           	
           	if (xi(i)*xi(i+1).lt.0.0d0) then
              	   j=i
              	   exit 
           	end if
        enddo
           
           if (mu+intra(j).gt.0.0d0) then
           	muf=((mu+intra(j))+(mu+intra(j+1)))/2
           	else
           	muf= mu+intra(1)
           end if
           
           
       
         
        WRITE(2,"(F8.3,X,F9.3,X,F9.4,X,F9.4,X,F11.7,X,F14.10,F11.7,X,F11.7,X,F11.7,X,F11.7,X,F11.7,X,F11.7,X,F11.7)") &
             mu, rs_n, Ef_e, Ef_h, n_tot, Dmax, Vinmax, CF/n_tot, kmin/kf_n, pairlen, muf, Kin, Corr


        WRITE(5,"(F8.3,X,F9.3,X,F9.4,X,F9.4)") rs_n, d, CF/n_tot
         
         
         
        
        call flush(2)  
        call flush(5)
        WRITE(schermo,*) 'mu=',mu, 'mu_e=',mu_e, 'mu_h=',mu_h 
        write(schermo,*) 'Delta_max=',Dmax
        WRITE(schermo,*) 'n_e=',ne, 'n_h=',nh, 'n_tot=',n_tot
        write(schermo,*) 'rs=',rs_n, 'kF=',kf_n  
        write(schermo,*) 'EF=',Ef_n, 'EF_e=',Ef_e, 'EF_h=',Ef_h
        write(schermo,*) '------------------------------------------------------------------------'
        call flush(4)

        if (CF/n_tot.gt.0.2d0) then  
            exit
        else        
          
          mu=mu+0.25d0      
        
        end if
     
   END DO

    CLOSE(UNIT=1)
    CLOSE(UNIT=2)
    CLOSE(UNIT=3)
    CLOSE(UNIT=4)
    CLOSE(UNIT=5)
END DO
END PROGRAM SF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Calculate Points and Weight of N-order Legendre polinomia in [x1,x2] integral extremes 
!!! X is the vector containing the points for the integral
!!! W is the vector with the waits assosiated to each element in  X
!!! N number of points for the integral
     SUBROUTINE gauleg(X1,X2,X,W,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X(N),W(N)
      PARAMETER (EPS=3.D-14)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
        Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
1       CONTINUE
          P1=1.D0
          P2=0.D0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE spline(x, y, n, yp1, ypn, y2)
! x(n) and y(n) containing a tabulated function,with x(1)<x(2)<...<x(n), given values yp1 and ypn for the first derivative  
! at points 1 and n, returns an array y2(n) which contains the second derivatives of the interpolating function at the tabulated points x(i).
! If yp1 and/or ypn are equal to 1.e30 or larger,the routine is signaled to set the 
! corresponding boundary condition for a natural spline with zero second derivative on that boundary.
   INTEGER, PARAMETER :: DP=KIND(1.0D0)
   INTEGER:: n
   INTEGER, PARAMETER:: nmax=500
   REAL(DP):: yp1, ypn, x(n), y(n), y2(n)
   INTEGER:: i, k
   REAL(DP):: p, qn, sig, un, u(nmax)

     if (yp1.gt.0.99e30) then
        y2(1)=0.d0
        u(1)=0.d0
     else
        y2(1)=-0.5
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
     end if

     do i=2, n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
     enddo

     if (ypn.gt.0.99e30) then
        qn=0.
        un=0.
     else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
     endif

     y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

     do k=n-1, 1, -1
        y2(k)=y2(k)*y2(k+1)+u(k)
     enddo
    RETURN

END 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION splint(xa,ya,y2a, n, x)
!Given the arrays xa and ya, which tabulate a function (with the xa(i) in order), 
!and given the array y2a which is the output from spline above, and given a value of x, 
!this routine returns a cubic-spline interpolated value.
   INTEGER, PARAMETER :: DP = KIND(1.D0)
   INTEGER:: n
   REAL(DP):: x, splint, xa(n), y2a(n), ya(n)
   INTEGER:: k, khi, klo
   REAL(DP):: a, b, h

     klo=1
     khi=n
1   if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (xa(k).gt.x) then
           khi=k
        else
           klo=k
        endif
        goto 1
     endif

     h=xa(khi)-xa(klo)
     if (h.eq.0.) stop 'bad xa input in splint'

     a=(xa(khi)-x)/h
     b=(x-xa(klo))/h
     splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
END FUNCTION splint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Calculate integral of angular part!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION ang_intintra(db, xphi, wphi, k, k1, q, P0ne, P0nh, P0a, y2P0ne, y2P0nh, y2P0a)
    implicit none
    integer,parameter::DP=kind(0.D0),nphi=250,nq=301
    real(DP):: db, xphi(nphi), wphi(nphi), k, k1, qab, q(nq), P0ne(nq), P0nh(nq), P0a(nq), y2P0ne(nq), y2P0nh(nq), y2P0a(nq) 
    real(DP):: pi, l1, l2, qi, Pne, Pnh, Pa, ang_int, ang_intintra, splint, t1eh,t1ee,t1hh, ffeh,ffee,ffhh
    real(DP):: intg, Pot_scr, Pot_se, Pot_sh,Pot_d, Pot_eff, Pot_intraeff, intgintra, termquad,Ptot_e,Ptot_h,termquadP
    integer::i
    external splint

    pi=acos(-1.0)
    
    l1=q(1)
    l2=q(nq)
    intgintra=0.d0
    !integral in phi
   do i = 1, nphi
        qi=sqrt(k**2+k1**2-2.0*k*k1*cos(xphi(i)))     !q=|k-k'|
         
        ! I define the polarizabilities
            Pne=0.d0
            Pnh=0.d0
            Pa=0.d0
    

        Pot_se=1.d0/qi
        Pot_sh=1.d0/qi                        !unscreen intralayer Potential V(phi,k,k')
        Pot_d=exp(-db*qi)/qi                  !unscreen interlayer Potential V(phi,k,k')

 
        Ptot_e=Pne+Pa  
        Ptot_h=Pnh+Pa   
        termquad=(Pot_se*Pot_sh)-(Pot_d**2)
        termquadP=(Pne*Pnh)-Pa**2

 
        
         Pot_scr=1.d0-(2.d0*((Pne*Pot_se)+(Pa*Pot_d)))+(termquad*termquadP)
         Pot_intraeff= (Pot_se)
          
	    ! Pot_intraeff=Pot_se  !unscreened
      
   
        intgintra=intgintra+wphi(i)*Pot_intraeff       
      
    end do
	
    
    ang_intintra=intgintra/pi                           !divided by pi to normilize in Ry 
            
END FUNCTION ang_intintra






FUNCTION ang_int(db, xphi, wphi, k, k1, q, P0ne, P0nh, P0a, y2P0ne, y2P0nh, y2P0a)
    implicit none
    integer,parameter::DP=kind(0.D0),nphi=250,nq=301
    real(DP):: db, xphi(nphi), wphi(nphi), k, k1, q(nq), P0ne(nq), P0nh(nq), P0a(nq), y2P0ne(nq), y2P0nh(nq), y2P0a(nq) 
    real(DP):: pi, l1, l2, qi, Pne, Pnh, Pa, ang_int, ang_intintra, splint, t1eh,t1ee,t1hh, ffeh,ffee,ffhh
    real(DP):: intg, Pot_scr, Pot_se, Pot_sh,Pot_d, Pot_eff, Pot_intraeff, intgintra, termquad,Ptot_e,Ptot_h,termquadP
    integer::i
    external splint

    pi=acos(-1.0)
    
    l1=q(1)
    l2=q(nq)
    intg=0.d0
    
    
    !integral in phi
    do i = 1, nphi
        qi=sqrt(k**2 + k1**2 - 2.0*k*k1*cos(xphi(i)))     !q=|k-k'|
        ! I define the polarizabilities
            Pne=0.d0
            Pnh=0.d0
            Pa=0.d0
   
        Pot_se=1.d0/qi
        Pot_sh=1.d0/qi                        !unscreen intralayer Potential V(phi,k,k')
        Pot_d=exp(-db*qi)/qi                  !unscreen interlayer Potential V(phi,k,k')

        Ptot_e=Pne+Pa  
        Ptot_h=Pnh+Pa   
        termquad=(Pot_se*Pot_sh)-(Pot_d**2)
        termquadP=(Pne*Pnh)-Pa**2
        
         Pot_scr=1.d0-(2.d0*((Pne*Pot_se)+(Pa*Pot_d)))+(termquad*termquadP)
         Pot_eff= (Pot_d)
         ! Pot_eff=Pot_d   !unscreened
	
      
    intg=intg+wphi(i)*Pot_eff

    end do  


    ang_int=intg/pi         !divided by pi to normilize in Ry 
            
END FUNCTION ang_int


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Gap equation Delta(k)
!FUNCTION gapEq(i, ang_intintras, k, x, w, mu, db, xphi, wphi, q, P0ne, P0nh, P0a, y2P0ne, y2P0nh, y2P0a, nk, delta)

SUBROUTINE gapEq(i, intra, k, x, w, mu, nk, delta, gap)
     implicit none
    integer,parameter::DP=kind(0.D0),  nq=301, nphi=250, n=1000
    integer::h, i
    real(DP)::delta(n), gap, mu, muintra, intradgap, db, k3
    real(DP)::k, xi, En, nk(n), x(n), w(n), xphi(nphi), wphi(nphi), intra(n)
    real(DP)::xi0, En0, vk20
    real(DP):: ang_intintras(n, n)
    
   
    muintra=mu+intra(i)

    xi=0.5*(k**2) - muintra          !xi=(xie+xih)/2 
    En=sqrt((xi**2) + (delta(i)**2))

    gap=delta(i)/(2.d0*En)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Bogoliubov terms 
!FUNCTION vk2(i, ang_intintras, k, x, w, mu, db, xphi, wphi, q, P0ne, P0nh, P0a, y2P0ne, y2P0nh, y2P0a, nk, delta)
SUBROUTINE denseq(i, intra, k, x, w, mu, nk, delta, vk2)
    implicit none
    integer,parameter::DP=kind(0.D0),nq=301, nphi=250, n=1000
    integer::h, i
    real(DP)::delta(n), vk2, mu, muintra, intrad, db, k3
    real(DP)::k, xi, En, nk(n), x(n), w(n), xphi(nphi), wphi(nphi), intra(n)
    real(DP)::xi0,En0, vk20
    real(DP):: ang_intintras(n, n)
   
   muintra=mu+intra(i)
   
    xi=0.5*(k**2) - muintra
    En=sqrt((xi**2) + (delta(i)**2))                 !xi=(xie+xih)/2 
    

    vk2= 0.5d0*(1.d0-(xi/En))
    

END SUBROUTINE



SUBROUTINE BindingE(db,BE)
    implicit none
    integer,parameter::DP=kind(0.D0),n=1000,na=10000
    integer:: i, j
    real(DP)::pi, func, db
    real(DP)::p1, w1, int_e, f1, BE
    real(DP):: x(n), w(n), alpha(na),En_a(na)
    external func

    pi=acos(-1.d0)
    !!inizializzo i vettori e li azzero
    x=0.d0
    w=0.d0
    call gauleg(0.d0, 50.d0, x, w, n)                    !Calculate Points and Weight for quadrature

    DO j=1,na
    alpha(j)=j*0.0005
    int_e=0.d0
        do i=1,n
        p1=x(i)
        w1=w(i)
        f1=func(p1, alpha(j),db)
        int_e=int_e + f1*w1
        end do
    En_a(j)=2.d0*pi*int_e
    END DO
   
    BE=MINVAL(En_a)

END subroutine

!!!!!!!!!!!!!!!
FUNCTION func(rho, alpha, d)
    implicit none
    integer,parameter::DP=kind(0.D0)
    real(DP)::rho, d, alpha, func, phi, pi, N,N1, N2, H
    real(DP):: A, A1, A2, A3, B, C, term, term2, alpha2, rho2

    pi=acos(-1.d0)
    rho2=rho**2
    alpha2=alpha**2

    term2=(d**2) + (rho**2)
    term=sqrt(term2)

    phi=exp(-term/alpha)

    A1=rho2/(alpha*term*term2) !!!!
    A2=rho2/(alpha2*term2)
    A3=1.d0/(alpha*term)
    A=(A1+A2-A3)*phi

    B=-phi/(alpha*term)
 
    C=2.d0*phi/term

    H=-A-B-C

    N1=exp(2.d0*d/alpha)
    N2=alpha2*(1.d0+ (2.d0*d/alpha))
    N=2.d0*N1/(pi*N2)

    func=N*rho*phi*H 

END FUNCTION
