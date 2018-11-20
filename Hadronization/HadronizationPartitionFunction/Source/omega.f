c$$$      implicit none 
c$$$      double precision massc,phi,ephi
c$$$      integer maxpar
c$$$      parameter(maxpar=50)
c$$$      double precision m(maxpar)
c$$$      integer npointsmc
c$$$      common/point/npointsmc
c$$$      integer n,i
c$$$
c$$$      print*,'enter number of MC points'
c$$$      read*,npointsmc
c$$$      print*,'enter massc (GeV)'
c$$$      read*,massc
c$$$      print*,'enter number of final state particles (max 50)'
c$$$      read*,n
c$$$      do i=1,n
c$$$         m(i)= 0.d0
c$$$      enddo
c$$$      print*,'enter the n values of the particle masses'
c$$$      do i=1,n
c$$$         read*,m(i)
c$$$      enddo
c$$$
c$$$      call omega(massc,n,m,phi,ephi)
c$$$      print*,'omega= ',phi,' +- ',ephi
c$$$      stop
c$$$      end

c      subroutine omega(mass,n,m,phi,ephi)
      subroutine omega()
      logical warn
      integer maxpar
c      parameter(maxpar=50)     !cambiare len
      parameter(maxpar=150)     !cambiare len
c      integer n
      integer n,npointsmc
      double precision mass
      double precision m
c      double precision m(maxpar),mi(maxpar)
      double precision mi(maxpar)
      double precision phi,ephi
      double precision etot,masssum,ttot
      double precision r(maxpar),z(maxpar),x(maxpar),s(maxpar)
      double precision t(maxpar),e(maxpar),p(maxpar),pinw(maxpar)
      double precision tmp,etmp
      double precision w
      double precision wint,wsig
      double precision psi
      double precision pi
      common/constants/pi
c      common/point/npointsmc
      common/inputdata/mass,m(maxpar),npointsmc,n
      common/outputdata/phi,ephi

      integer j
            
c for MC integration
      integer lux,int,k1,k2,len
      parameter (lux=4,int=701849210,k1=0,k2=0,len=150) !len must be=maxpar!!
      real*4 vec(len)
c      integer npointsmc,icnt,nhit,nbias
      integer icnt,icnt2,icnt3,icnt4,oldicnt,nhit,ninst,ninsttot
      double precision fcount,fval,fav,f2av,efav,nmc,sqrtnm1

      pi= 3.1415926536d0

      etot= mass

      do i=1,n
         mi(i)= m(i)
      enddo
      masssum= 0.d0
      do i=1,n
         masssum= masssum+mi(i)
      enddo
      ttot= etot-masssum

      icnt = 0
      icnt2 = 0
      icnt3 = 0
      icnt4 = 0
      oldicnt = npointsmc
      fval= 0.
      fav= 0.
      f2av= 0.
      ninst = 0
      ninsttot = 0

c      print*,'massa cluster ',mass
c      print*,'numero punto integrazione ',npointsmc
c      print*,'numero particelle ',n
      

      if(mass.gt.20.d0)then
         phi = -1.d0
         ephi = -1.d0
         print*,'Massa clu sup a 20 GeV,problema vettori'
         return
      endif

*
*-----random number generator initialization
*
      call rluxgo(lux,int,k1,k2)
*
*-----Monte Carlo loop
*
 21   do i = 1,npointsmc
*     
         icnt = icnt + 1
*              
*-----RANDOM NUMBER GENERATION
*     
         call ranlux(vec,len)
         do j=1,len
            r(j)= vec(j)*1.d0
         enddo
*     
         do j=1,n-1
            z(j)= r(j)**(1.d0/j)
         enddo
         x(n)= 1.d0
         s(n)= ttot
         do j=n-1,1,-1
            x(j)= x(j+1)*z(j)
            s(j)= ttot * x(j)
            t(j+1)= s(j+1) - s(j)
            e(j+1)= t(j+1) + mi(j+1)
            p(j+1)= sqrt(t(j+1)*(t(j+1)+2.d0*mi(j+1)))
            pinw(j+1)= p(j+1)
         enddo
         t(1)= s(1)
         e(1)= t(1) + mi(1)
         p(1)= sqrt(t(1)*(t(1)+2.d0*mi(1)))
         pinw(1)= p(1)
         
c     pinw(1) = 0.0588691d0
c     pinw(2) = 0.129392d0
c     pinw(3) = 0.0293137d0
c     pinw(4) = 0.109888d0
c     pinw(5) = 0.120457d0
c     pinw(6) = 0.0351444d0
c     pinw(7) = 1.52307d0
c     pinw(8) = 0.0476349d0
c     pinw(9) = 0.0635834d0
c     pinw(10)= 0.038534d0
c     pinw(11)= 0.125991d0
c     pinw(12)= 0.126978d0
c     pinw(13)= 0.147783d0
c     pinw(14)= 0.0497111d0
c     pinw(15)= 0.332001d0
         
         tmp= 0.d0
         etmp= 0.d0
         do j=1,n
            tmp= tmp+t(j)
            etmp= etmp+sqrt(p(j)**2+mi(j)**2)
         enddo

c     print*,'ttot= ',tmp,ttot
c     print*,'etot= ',etmp,etot
         
         if(n.le.14) then
c            print*,'Poche part ',n
            call wp1pn_sigma(pinw,n,wsig)
            w= wsig
c            print*,'wsig= ',wsig
         else
            call wp1pn(pinw,n,wint,warn)
            if(warn) then
               ninst = ninst + 1
               ninsttot = ninsttot + 1
               icnt = icnt - 1               
               if(ninst.le.10)then
                  goto 21
               else
                  phi = -1.d0
                  ephi = -1.d0
                  print*,'11 instabilità successive'
                  return
               endif
            else
c     print*,'w integrata= ',wint
               ninst = 0
               w= wint
            endif
         endif
         

c     if(w.gt.0.d0) then 
c     print*,'w= ',w
c     endif
         
         psi= 1.d0
         do j=1,n
            psi= psi * (4.d0*pi)/(j*1.d0) * ttot *p(j)*e(j)
         enddo
         psi= psi*n/ttot*w
         
         fval= psi
         fav= fav+fval
         f2av= f2av+fval*fval
         nmc= float(icnt)
         
C         if(ninsttot.gt.100)then
C            print*,'Warning, più di 100 inst totali'
C            phi = -1.d0
C            ephi = -1.d0
C            return
C         endif

c         if((icnt.ge.10).and.(icnt2.eq.0))then
         if(icnt.ge.10)then
            efav = ((f2av - fav*fav/nmc)/(fav*fav))*nmc/(nmc-1.e0)
c     print*,'Numero punti ',icnt,' ',icnt2
c     print*,'Valore ',efav
            if(efav.lt.1.e-4)then
               goto 23
c            else
c               icnt2 = 10*icnt
c                icnt4 = icnt/2
c                icnt2 = icnt + icnt4
c                icnt3 = 10*icnt
            endif
         endif
         
c$$$         if((icnt.gt.0).and.(icnt.eq.icnt2))then
c$$$            efav = ((f2av - fav*fav/nmc)/(fav*fav))*nmc/(nmc-1.e0)
c$$$cc            print*,'Numero punti ',icnt,' ',icnt2
c$$$c     c            print*,'Valore ',efav
c$$$            if(efav.lt.1.e-4)then
c$$$               goto 23
c$$$            else
c$$$               icnt2 = 10*icnt
c$$$c              if(icnt.lt.icnt3)then
c$$$c                  icnt2 = icnt + icnt4
c$$$c               else
c$$$c                  icnt4 = icnt/2
c$$$c                  icnt2 = icnt + icnt4
c$$$c                  icnt3 = 10*icnt
c$$$c               endif
c$$$            endif
c$$$         endif

      enddo
      
 23   continue

      nmc= float(icnt)
      sqrtnm1= sqrt(nmc-1.e0)
      fav= fav/nmc
      f2av= f2av/nmc
      efav= sqrt(f2av-fav**2) / sqrtnm1

      phi= fav
      ephi= efav
    
c$$$      print*,' '
c$$$      print*,'results from MC integration'
c$$$      print*,'n. of points= ',nmc,'    integral= ',fav,' +- ',efav,' ',
c$$$     &efav/fav

      return
      end

      subroutine wp1pn(p,n,w,warn)
      implicit none
      logical warn
      integer maxpar,n
c      parameter(maxpar=50)
      parameter(maxpar=150)     
      double precision p(maxpar)
      double precision w
      double precision a,b,s
      double precision eps0,eps
      integer jmax,jmax0     !2^jmax-1
      integer i,ncount2,ncount3
      double precision pi
      common/boundary/a,b
      common/constants/pi
      double precision pc(maxpar)
      common/momenta/pc
      integer nc
      common/nparticles/nc
      double precision func
      external func
*
      nc= n
      do i=1,nc
         pc(i)= p(i)
      enddo
      a= 0.d0
      b= 1000.d0
      eps0= 1.d-2
c      jmax0= 12
      jmax0= 15
*
* integration with trapezoidal rule
*
      s= 0.
      eps= eps0
      jmax= jmax0
      call qtrap(func,a,b,s,eps,jmax,ncount2,warn)
c      print*,'n. of points trap= ',ncount2,'    integral= ',s
*
* integration with Sympson rule
*
c      s= 0.
c      eps= eps0
c      jmax= jmax0
c      call qsymp(func,a,b,s,eps,jmax,ncount3)
c      print*,'n. of points symps= ',ncount3,'    integral= ',s


      w= s

      return
      end
*
      subroutine wp1pn_sigma(p,n,w)
      implicit none
c      integer maxpar,maxpart,nsigma,n
      integer maxpar,maxpar2,maxpart,nsigma,n
c      parameter(maxpar=50,maxpart=2**20)
      parameter(maxpar=150)
      parameter(maxpar2=14,maxpart=2**14)
      double precision p(maxpar)
      double precision w
      double precision aux
      double precision sumsigjpj
      double precision sigjpj
      double precision sigjpjprod
      double precision tmp
      double precision pi
      common/constants/pi
      
c      double precision pc(maxpar)
c      common/momenta/pc
      double precision pcs(maxpar2)
      
      integer nc
      common/nparticles/nc

      integer jnu
c      integer sigma(maxpar,maxpart),prodsigma(maxpart)
      integer sigma(maxpar2,maxpart),prodsigma(maxpart)
      integer succession(maxpart)
      double precision sigjpjott(maxpart),sigjpjnonott(maxpart)
      integer a,b,c
      integer i,j
      integer factor
      integer imax,x
      integer f(20)

ccc Nelle repliche cambiare il nome di initial

      integer initial
      data initial/0/
      save initial
      save succession
      save sigma
      save prodsigma

ccc Aggiunto per tenere conto del cambio di dim in sigma

      integer oldnc
      data oldnc/0/
      save oldnc
*
      nc = n
*
* this is to be done just one time
*

      if((initial.eq.0).or.(nc.gt.oldnc))then
         do i=1,nc
            sigma(i,1)= 1
         enddo
         oldnc = nc
         initial = 0
      endif  

      if(initial.eq.0) then

* determination of array jnu labelling the sigma vectors
         do jnu=1,maxpart           !2^14
            x= jnu
            do i=1,10000000
               if(x.eq.1) then 
                  f(i)= 1
                  imax= i
                  goto 20
               endif
               f(i)= mod(x,2)
               x= (x-f(i))/2
            enddo
 20         do i=1,imax
               if(f(i).eq.1) then !f(i) gives jnu in basis 2
                  succession(jnu)=i
                  goto 30
               endif
            enddo
 30         continue
            if(jnu.gt.1) then
               do i=1,nc
                  sigma(i,jnu)= sigma(i,jnu-1)
                  if(i.eq.succession(jnu-1)) then
                     sigma(i,jnu)= -sigma(i,jnu-1)
                  endif
               enddo
            endif
*
            prodsigma(jnu)= 1
            do i=1,nc
               prodsigma(jnu)= prodsigma(jnu)*sigma(i,jnu)
            enddo
         enddo
         initial= 1
      endif
*
* this has to be calculated at every call
*
      do i=1,nc
         pcs(i)= p(i)
      enddo

      sumsigjpj= 0.d0
      nsigma= 2**nc
      sigjpjott(1)= 0.d0
      do i=1,nc
         sigjpjott(1)= sigjpjott(1) + sigma(i,1)*pcs(i)
      enddo
      sumsigjpj= sumsigjpj + sigjpjott(1)**(nc-3) * prodsigma(1)

c      do j= 1,nsigma
c         sigjpj= 0.d0
c        do i=1,nc
c            sigjpj= sigjpj+sigma(i,j)*pc(i)
c         enddo
c         sigjpjnonott(j)= sigjpj
c      enddo

      do j= 2,nsigma
         sigjpjott(j)= sigjpjott(j-1)
     +            - 2.d0*sigma(succession(j-1),j-1)*pcs(succession(j-1))
         if(sigjpjott(j).lt.0.d0) then
            sigjpjprod= 0.d0
         else
            sigjpjprod= (sigjpjott(j)**(nc-3)) * prodsigma(j)
         endif
         sumsigjpj= sumsigjpj + sigjpjprod
      enddo
c      print*,' '
c      do j=1,nc
c         print*,'sumsigjpj= ',sigjpjott(j),sigjpjnonott(j)
c      enddo

      aux= 1.d0
      do j=1,nc
         aux= aux*pcs(j)
      enddo


      w= -1.d0/2.d0**(nc+1)/pi/factor(nc-3)/aux * sumsigjpj

      return
      end



*
* factorial
*
      function factor(x)
      implicit none
      integer x,y
      integer factor
      if(x.eq.0) then
         factor= 1
         return
      endif
      y= x
      factor= y
 10   if(y.gt.1) then
         y= y-1
         factor= factor*y
         goto 10
      else
         return
      endif
      return
      end
*
* integrand
*
      function func(x)
      implicit none
      integer maxpar
c      parameter(maxpar=50)
      parameter(maxpar=150)
      double precision func
      double precision x
      integer k
      double precision a,b
      double precision prod
      double precision pi
      common/boundary/a,b
      common/constants/pi
      double precision pc(maxpar)
      common/momenta/pc
      integer nc
      common/nparticles/nc

      prod= 1.d0
      do k= 1,nc
         if(abs(x).lt.1.d-8) then 
            prod= prod
         else
            prod= prod*sin(pc(k)*x)/(pc(k)*x)
         endif
      enddo

c      func= (b-a) * 1.d0/2.d0/pi/pi *x*x * prod   !b-a (jacobian)
      func= 1.d0/2.d0/pi/pi *x*x * prod   !b-a (jacobian)
 
      return
      end

c
c uses trapzd
c
c returns as s the integral of the function func from a to b. 
c eps is the required fractional accuracy and jmax is a parameter 
c limiting the maximum number of steps to 2^(jmax-1)
c
      subroutine qtrap(func,a,b,s,eps,jmax,ncalls,warn)
      implicit none
      logical warn
      integer jmax,ncalls,j,ncount2
      double precision a,b,eps,s,func
      external func

      double precision olds,old
      integer jcount

      olds= -1.d30 !any number who is unlikely to be the average of 
c                   the function at its endpoints

c     Aggiunto da me per fare il calcolo con le sigma
c     se non cnverge l'integrale

      warn = .false.

      do j=1,jmax
         call trapzd(func,a,b,s,j,jcount)
         ncalls= jcount
         if(j.gt.5) then   !avoids spurious convergence
            if(abs(s-olds).lt.eps*abs(olds).or.
     +        (s.eq.0.d0.and.olds.eq.0.d0))return
         endif
         olds=s
c         print*,'j= ',j,'   int= ',s
         if(j.eq.4) old= olds
      enddo
      warn = .true.
c      print*,'too many steps in qtrap',olds,old
      return
      end

c
c uses trapzd
c
c returns as s the integral of the function func from a to b 
c according to the Sympson rule. 
c eps is the required fractional accuracy and jmax is a parameter 
c limiting the maximum number of steps to 2^(jmax-1)
c
      subroutine qsymp(func,a,b,s,eps,jmax,ncalls)
      implicit none
      integer jmax,ncalls,j,ncount2
      double precision a,b,eps,s,func
      external func

      double precision os,ost,st
      integer jcount

      ost= -1.d30 !any number who is unlikely to be the average of 
c                   the function at its endpoints
      os= -1.d30
      do j=1,jmax
         call trapzd(func,a,b,st,j,jcount)
         ncalls= jcount
         s= (4.d0*st - ost)/3.d0
         if(j.gt.5) then   !avoids spurious convergence
            if(abs(s-os).lt.eps*abs(os).or.
     +        (s.eq.0.d0.and.os.eq.0.d0)) return
         endif
         os=s
         ost= st
      enddo
c      print*,'too many steps in qsymp'
      return
      end

c
c computes the nth stage of refinement of the trapezoidal rule. func 
c is the function to be integrated between a and b, also input. When 
c called with n=1, the routine returns as s the crudest estimate of 
c the integral. Subsequent calls with n=2,3... (in sequential order) 
c will improve the accuracy by adding 2^(n-2) additional interior 
c points. s should not be modified between sequential calls.
c
      subroutine trapzd(func,a,b,s,n,jcount)
      implicit none
      integer n,jcount
      double precision a,b,s,func
      external func

      integer it,j
      integer icount,icount0
      double precision del,sum,tnm,x
      save icount

      if(n.eq.1) then
         s= 0.5d0*(b-a)*(func(a)+func(b))
         icount= 2
      else
         it= 2**(n-2)
         tnm= it
         del= (b-a)/tnm
         x= a+0.5d0*del
         sum= 0.d0
         do j=1,it
            sum= sum+func(x)
            x= x+del
         enddo
         s= 0.5d0*(s+(b-a)*sum/tnm)
         icount= icount+it
      endif
      jcount=icount
      return
      end

*
*-----SUBROUTINE RANLUX: RANDOM NUMBER GENERATOR
*
      SUBROUTINE RANLUX(RVEC,LENV)
C         SUBTRACT-AND-BORROW RANDOM NUMBER GENERATOR PROPOSED BY
C         MARSAGLIA AND ZAMAN, IMPLEMENTED BY F. JAMES WITH THE NAME
C         RCARRY IN 1991, AND LATER IMPROVED BY MARTIN LUESCHER
C         IN 1993 TO PRODUCE "LUXURY PSEUDORANDOM NUMBERS".
C     FORTRAN 77 CODED BY F. JAMES, 1993
C
C   LUXURY LEVELS.
C   ------ ------      THE AVAILABLE LUXURY LEVELS ARE:
C
C  LEVEL 0  (P=24): EQUIVALENT TO THE ORIGINAL RCARRY OF MARSAGLIA
C           AND ZAMAN, VERY LONG PERIOD, BUT FAILS MANY TESTS.
C  LEVEL 1  (P=48): CONSIDERABLE IMPROVEMENT IN QUALITY OVER LEVEL 0,
C           NOW PASSES THE GAP TEST, BUT STILL FAILS SPECTRAL TEST.
C  LEVEL 2  (P=97): PASSES ALL KNOWN TESTS, BUT THEORETICALLY STILL
C           DEFECTIVE.
C  LEVEL 3  (P=223): DEFAULT VALUE.  ANY THEORETICALLY POSSIBLE
C           CORRELATIONS HAVE VERY SMALL CHANCE OF BEING OBSERVED.
C  LEVEL 4  (P=389): HIGHEST POSSIBLE LUXURY, ALL 24 BITS CHAOTIC.
C
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  CALLING SEQUENCES FOR RANLUX:                                  ++
C!!!      CALL RANLUX (RVEC, LEN)   RETURNS A VECTOR RVEC OF LEN     ++
C!!!                   32-BIT RANDOM FLOATING POINT NUMBERS BETWEEN  ++
C!!!                   ZERO (NOT INCLUDED) AND ONE (ALSO NOT INCL.). ++
C!!!      CALL RLUXGO(LUX,INT,K1,K2) INITIALIZES THE GENERATOR FROM  ++
C!!!               ONE 32-BIT INTEGER INT AND SETS LUXURY LEVEL LUX  ++
C!!!               WHICH IS INTEGER BETWEEN ZERO AND MAXLEV, OR IF   ++
C!!!               LUX .GT. 24, IT SETS P=LUX DIRECTLY.  K1 AND K2   ++
C!!!               SHOULD BE SET TO ZERO UNLESS RESTARTING AT A BREAK++
C!!!               POINT GIVEN BY OUTPUT OF RLUXAT (SEE RLUXAT).     ++
C!!!      CALL RLUXAT(LUX,INT,K1,K2) GETS THE VALUES OF FOUR INTEGERS++
C!!!               WHICH CAN BE USED TO RESTART THE RANLUX GENERATOR ++
C!!!               AT THE CURRENT POINT BY CALLING RLUXGO.  K1 AND K2++
C!!!               SPECIFY HOW MANY NUMBERS WERE GENERATED SINCE THE ++
C!!!               INITIALIZATION WITH LUX AND INT.  THE RESTARTING  ++
C!!!               SKIPS OVER  K1+K2*E9   NUMBERS, SO IT CAN BE LONG.++
C!!!   A MORE EFFICIENT BUT LESS CONVENIENT WAY OF RESTARTING IS BY: ++
C!!!      CALL RLUXIN(ISVEC)    RESTARTS THE GENERATOR FROM VECTOR   ++
C!!!                   ISVEC OF 25 32-BIT INTEGERS (SEE RLUXUT)      ++
C!!!      CALL RLUXUT(ISVEC)    OUTPUTS THE CURRENT VALUES OF THE 25 ++
C!!!                 32-BIT INTEGER SEEDS, TO BE USED FOR RESTARTING ++
C!!!      ISVEC MUST BE DIMENSIONED 25 IN THE CALLING PROGRAM        ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (MAXLEV=4, LXDFLT=3)
      DIMENSION NDSKIP(0:MAXLEV)
      DIMENSION NEXT(24)
      PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
      SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
      INTEGER LUXLEV
      LOGICAL NOTYET
      DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/
      DATA I24,J24,CARRY/24,10,0./
C                               DEFAULT
C  LUXURY LEVEL   0     1     2   *3*    4
      DATA NDSKIP/0,   24,   73,  199,  365 /
CORRESPONDS TO P=24    48    97   223   389
C     TIME FACTOR 1     2     3     6    10   ON SLOW WORKSTATION
C                 1    1.5    2     3     5   ON FAST MAINFRAME
C
C  NOTYET IS .TRUE. IF NO INITIALIZATION HAS BEEN PERFORMED YET.
C              DEFAULT INITIALIZATION BY MULTIPLICATIVE CONGRUENTIAL
      IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = JSDFLT
         INSEED = JSEED
         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
         LUXLEV = LXDFLT
         NSKIP = NDSKIP(LUXLEV)
         LP = NSKIP + 24
         IN24 = 0
         KOUNT = 0
         MKOUNT = 0
         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
     +        LUXLEV,'      P =',LP
            TWOM24 = 1.
         DO 25 I= 1, 24
            TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
   25    CONTINUE
         TWOM12 = TWOM24 * 4096.
         DO 50 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
   50    CONTINUE
         NEXT(1) = 24
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
      ENDIF
C
C          THE GENERATOR PROPER: "SUBTRACT-WITH-BORROW",
C          AS PROPOSED BY MARSAGLIA AND ZAMAN,
C          FLORIDA STATE UNIVERSITY, MARCH, 1989
C
      DO 100 IVEC= 1, LENV
      UNI = SEEDS(J24) - SEEDS(I24) - CARRY
      IF (UNI .LT. 0.)  THEN
         UNI = UNI + 1.0
         CARRY = TWOM24
      ELSE
         CARRY = 0.
      ENDIF
      SEEDS(I24) = UNI
      I24 = NEXT(I24)
      J24 = NEXT(J24)
      RVEC(IVEC) = UNI
C  SMALL NUMBERS (WITH LESS THAN 12 "SIGNIFICANT" BITS) ARE "PADDED".
      IF (UNI .LT. TWOM12)  THEN
         RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
C        AND ZERO IS FORBIDDEN IN CASE SOMEONE TAKES A LOGARITHM
         IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
      ENDIF
C        SKIPPING TO LUXURY.  AS PROPOSED BY MARTIN LUSCHER.
      IN24 = IN24 + 1
      IF (IN24 .EQ. 24)  THEN
         IN24 = 0
         KOUNT = KOUNT + NSKIP
         DO 90 ISK= 1, NSKIP
         UNI = SEEDS(J24) - SEEDS(I24) - CARRY
         IF (UNI .LT. 0.)  THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
         ELSE
            CARRY = 0.
         ENDIF
         SEEDS(I24) = UNI
         I24 = NEXT(I24)
         J24 = NEXT(J24)
   90    CONTINUE
      ENDIF
  100 CONTINUE
      KOUNT = KOUNT + LENV
      IF (KOUNT .GE. IGIGA)  THEN
         MKOUNT = MKOUNT + 1
         KOUNT = KOUNT - IGIGA
      ENDIF
      RETURN
C
C           ENTRY TO INPUT AND FLOAT INTEGER SEEDS FROM PREVIOUS RUN
      ENTRY RLUXIN(ISDEXT)
         TWOM24 = 1.
         DO 195 I= 1, 24
         NEXT(I) = I-1
  195    TWOM24 = TWOM24 * 0.5
         NEXT(1) = 24
         TWOM12 = TWOM24 * 4096.
      WRITE(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
      WRITE(6,'(5X,5I12)') ISDEXT
      DO 200 I= 1, 24
      SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  200 CONTINUE
      CARRY = 0.
      IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
      ISD = IABS(ISDEXT(25))
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = MOD(ISD,100)
      ISD = ISD/100
      IN24 = MOD(ISD,100)
      ISD = ISD/100
      LUXLEV = ISD
        IF (LUXLEV .LE. MAXLEV) THEN
          NSKIP = NDSKIP(LUXLEV)
          WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
     +                         LUXLEV
        ELSE  IF (LUXLEV .GE. 24) THEN
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
        ELSE
          NSKIP = NDSKIP(MAXLEV)
          WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
          LUXLEV = MAXLEV
        ENDIF
      INSEED = -1
      RETURN
C
C                    ENTRY TO OUPUT SEEDS AS INTEGERS
      ENTRY RLUXUT(ISDEXT)
      DO 300 I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  300 CONTINUE
      ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
      IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
      RETURN
C
C                    ENTRY TO OUTPUT THE "CONVENIENT" RESTART POINT
      ENTRY RLUXAT(LOUT,INOUT,K1,K2)
      LOUT = LUXLEV
      INOUT = INSEED
      K1 = KOUNT
      K2 = MKOUNT
      RETURN
C
C                    ENTRY TO INITIALIZE FROM ONE OR THREE INTEGERS
      ENTRY RLUXGO(LUX,INS,K1,K2)
         IF (LUX .LT. 0) THEN
            LUXLEV = LXDFLT
         ELSE IF (LUX .LE. MAXLEV) THEN
            LUXLEV = LUX
         ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
            LUXLEV = MAXLEV
            WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
         ELSE
            LUXLEV = LUX
            DO 310 ILX= 0, MAXLEV
              IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
  310       CONTINUE
         ENDIF
      IF (LUXLEV .LE. MAXLEV)  THEN
         NSKIP = NDSKIP(LUXLEV)
C         WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
C     +        LUXLEV,'     P=', NSKIP+24
      ELSE
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
      ENDIF
      IN24 = 0
      IF (INS .LT. 0)  WRITE (6,'(A)')
     +   ' ILLEGAL INITIALIZATION BY RLUXGO, NEGATIVE INPUT SEED'
      IF (INS .GT. 0)  THEN
       JSEED = INS
C        WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
C     +      JSEED, K1,K2
      ELSE
        JSEED = JSDFLT
C        WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
      ENDIF
      INSEED = JSEED
      NOTYET = .FALSE.
      TWOM24 = 1.
         DO 325 I= 1, 24
           TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
  325    CONTINUE
      TWOM12 = TWOM24 * 4096.
         DO 350 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
  350    CONTINUE
      NEXT(1) = 24
      I24 = 24
      J24 = 10
      CARRY = 0.
      IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
C        IF RESTARTING AT A BREAK POINT, SKIP K1 + IGIGA*K2
C        NOTE THAT THIS IS THE NUMBER OF NUMBERS DELIVERED TO
C        THE USER PLUS THE NUMBER SKIPPED (IF LUXURY .GT. 0).
      KOUNT = K1
      MKOUNT = K2
      IF (K1+K2 .NE. 0)  THEN
        DO 500 IOUTER= 1, K2+1
          INNER = IGIGA
          IF (IOUTER .EQ. K2+1)  INNER = K1
          DO 450 ISK= 1, INNER
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY
            IF (UNI .LT. 0.)  THEN
               UNI = UNI + 1.0
               CARRY = TWOM24
            ELSE
               CARRY = 0.
            ENDIF
            SEEDS(I24) = UNI
            I24 = NEXT(I24)
            J24 = NEXT(J24)
  450     CONTINUE
  500   CONTINUE
C         GET THE RIGHT VALUE OF IN24 BY DIRECT CALCULATION
        IN24 = MOD(KOUNT, NSKIP+24)
        IF (MKOUNT .GT. 0)  THEN
           IZIP = MOD(IGIGA, NSKIP+24)
           IZIP2 = MKOUNT*IZIP + IN24
           IN24 = MOD(IZIP2, NSKIP+24)
        ENDIF
C       NOW IN24 HAD BETTER BE BETWEEN ZERO AND 23 INCLUSIVE
        IF (IN24 .GT. 23) THEN
           WRITE (6,'(A/A,3I11,A,I5)')
     +    '  ERROR IN RESTARTING WITH RLUXGO:','  THE VALUES', INS,
     +     K1, K2, ' CANNOT OCCUR AT LUXURY LEVEL', LUXLEV
           IN24 = 0
        ENDIF
      ENDIF
      RETURN
      END
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C **********************************************************************
C
C **********************************************************************
