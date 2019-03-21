	             program proc1
      
c       The code is aimed at the determination
c       of snow and ice properties ( e.g., albedo)
c       using satellite observations (S-3)
      
c           A. KOKHANOVSKY
      
c             VERSION 1
c             19.02.2019
c             SICE
c     *************
      
c                                VERSION 2
c                                12.03.2019
c        DELETING INTERPOLATION ERROR in BBA calculation( See  'ALEX')


           
      real r(21),boar(21),rs(21),rp(21),rpalex(21),rsalex(21)
      real xs(21),astra(21)
c         ALEX
      real xa(168),ya(168)
             EXTERNAL fun1,fun2
			 COMMON as,bs,cs,am1,am2,csza,cvza,r400,r865,r1020,NSOLO,
     +						thv, xa, ya
             open(898,file='ice_index.dat')
              do 67 j=1,168  
             read(898,*) xa(j),an,ya(j)
 67       continue
c**********************end of modification***************ALEX
          
 102          format (i4,190( 4x,e12.4))
 1022        format(i4,7e12.4,4x,6i3)
 1024        format(i4,5x,f7.2, 5x,e12.4,5x,e12.4,5x,8i4)
                 pi=acos(-1.)

             
c              4 input files:

c*****************************************************************
c               sza,vza,saa,vaa, OLCI 21 reflectances :       
                open(1,file='input1_TOA.dat')
             
c               OLCI BOA 21 reflectances (after atmospheric correction):
                open(2,file='input2_BOA.dat')
                
c              maximal number of lines to be processed:
                 open(3,file='nlines.dat')
                 read(3,*) N
                 

                 
c     imaginary part of ice refractive index at 21 OLCI channels
                 open(901,file='kap_olci.dat')
                 
c     conversion coefficient ( EAL->>>snow grain size)
                 open(1001,file='psi.dat')
                 

c     if key=1 then the SNAP atmospheric correction is used
c     otherwise - not
                 open(1003, file='key.dat')
c******************************************************************     

c                    output files:

c spectral albedo: 21 planar and 21 spherical albedo
                  open(100,file='output_sp_albedo.dat')
c      BBA: planar and spherical
                  open(701,file='output_bba.dat')
c       EAL, SGS, SSA, ndsi,ndbi,etc
                  open(111, file='output_flags.dat')
c       IMPURITY TYPE AND LOAD
                  open(5001, file='output_impurity.dat')   
c*****************************************************************
c                STEP 1
c                reading satellite and ice refractive index data
c******************************************************************
c                conversion coefficient EAL-diameter of grains and thv for OLCI 1st channel
                  read(1001,*) psi
              
                  read(1003,*) key
c     reading coefficients for spectral ice refractive index
                  

                  
                   do 334 nkd=1,21
                     read(901,*) xs(nkd),astra(nkd)
 334              continue
c                 reading satellite data     
       do 1983 j=1,N              
          read(1,*,end=200) js,sza,vza,saa,vaa,(r(i),i=1,21)
          
           if (key.eq.1) go to 898
          do 897 ihk=1,21
               boar(ihk)=r(ihk)
 897        continue
            go to 899
 898   continue
       read(2,*,end= 200) js,(boar(i),i=1,21)
    
 899          continue

c            for grain size and pollution:
               r400=boar(1)
			   write (*,*) "r400: ", r400
               r560=boar(6)
               r865=boar(17)
               r1020=boar(21)
               
c             for algae:               
               r681=boar(10)
               r709=boar(11)
               ratio=r709/r681
c*********************END OF DATA READING***********************


                            
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


c     STEP 2
               
c     FLAGS and SNOW GRAIN SIZE
c     SNOW POLLUTION LOAD
               
               indexs=0
               indexi=0
               indexd=0
               indexc=0
               
c         NDSI
               andsi=(r(17)-r(21))/(r(17)+r(21))
               if (aNDSI.gt.0.03.and.r(1).gt.0.5) indexs=1
               
c         NDBI
               andbi=(r(1)-r(21))/(r(1)+r(21))
               if (aNDBI.gt.0.33) indexi=1
               if (aNDBI.gt.0.66) indexd=1
               
c             relative azimuthal angle calculation:
               raa=abs(180.-(vaa-saa))
c********************************************               
               am1=cos(pi*sza/180.)
			   csza=am1
               am2=cos(pi*vza/180.)
			   cvza=am2
               sam1=sin(pi*sza/180.)
               sam2=sin(pi*vza/180.)
c             scattering angle calculation:               
               co=-am1*am2+sam1*sam2*cos(raa*pi/180.)
               scat =acos(co)*180./pi
               
c             angular functions:
               u1=3.*(1.+2.*am1)/7.
               u2=3.*(1.+2.*am2)/7.
               
         alam3=0.865
         alam4=1.02
         akap3=2.4e-7
         akap4=2.25e-6

         
         alpha3=4.*pi*akap3/alam3
         alpha4=4.*pi*akap4/alam4
         
         eps=sqrt(alpha3/alpha4)
         ax1=1./(1.-eps)
         ax2=1./(1.-1./eps)
         
         r0=r865**ax1*r1020**ax2                  
         xx=u1*u2/r0
         
c      effective absorption length(microns):     
         alka=(alog(r1020/r0))**2./xx/xx/alpha4
         
c       mm
         al=alka/1000.
   
c     grain size
         dens=0.917
c         mm
         diam=al*psi
         if (diam.le.0.01) indexc=1
         delta=r(1)-r(21)
         if (delta.le.0.13) indexc=1
c           m*m/kg
         surf=6./diam/dens
         
c     SNOW POLLUTION
         
c        2 wavelengths:        
         al1=0.4
         al2=0.56
         
         ap1=(alog( r400/r0) ) **2.
         ap2=(alog(r560/r0)  ) **2.
         
           ang1=alog(ap1/ap2)
           ang2=alog(al2/al1)
           ang=ang1/ang2
           
c          inverse mm:
           AF=ap1*(0.4)**ang/al/xx/xx
         
            write(111,1022) js,al,diam,surf,scat,andsi,andbi,delta,
     c          indexs,indexi,indexd,indexc
            
            
            ipol=1
            ntype=4
            
            if (ang.gt. 10.) ipol=0
            if (ang.lt.0.5)    ipol=0
            if ( ang.ge.0.85.and.ang.le.1.15.and.ipol.ne.0) ntype=1
            if( ang.gt.1.15.and.ipol.ne.0) ntype=2
            if (ratio.gt.1.0) ntype=3
            bbsoot=0.9
            sootka=0.47
            shape=1.6
c           inverse mm:
            absoot=1000.*bbsoot*4.*pi*sootka
c          inverse mm:
            abdust=600.
            
c                 RELATIVE IMPURITY LOAD (RIL)
c     volumetric concentration of pollutants devided by
c     the volumetric concentration of ice grains:            
            if (ntype.eq.1) conc=AF*shape/absoot
            if (ntype.eq.2) conc=AF*shape/abdust
            
            if (ipol.eq.0) conc=0.0
c         algae:
            a680=1.e+7
            b680=600.
c           cells/ml
            if (ntype.eq.3) conc=a680*(alog(ratio))**2.-b680
            write(5001,1024) j,ang,af,conc,ipol,ntype

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@               
                

               
c++++++++++++++++++++++++++++++++++++
         
c               STEP 3
c     spectral albedo determination
                   r0al=falex1(am1,am2,co)
		write(*, *) "r0al: ", am1,am2,co, r0al
                           do 12 i=1,21
                  
                               if (r(i).gt.r0al)r0al=r(i)
                               rs(i)=(r(i)/r0al)**(r0al/u1/u2)
                               rp(i)=rs(i)**u1
                               
 12                         continue
              
                          thv=r0al-0.1
c               output of spectral albedo                           
        if(r400.le.thv) write(100,102) js,                   
     c          (rp(ix),ix=1,21),(rs(iz),iz=1,21)

        if(r400.le.thv) go to 700
        
        do 500 nk=1,21

              dega=alka*4.*pi*astra(nk)/xs(nk)
              pow=sqrt(dega)
              if (pow.gt.1.e-6)rsd=exp(-pow)
              if (pow.le.1.e-6) rsd=1.
              rsalex(nk)=rsd
              rpalex(nk)=rsd**u1

 500      continue
          
        if(r400.gt.thv) write(100,102) js,                   
     c (rpalex(i),i=1,21),(rsalex(i),i=1,21)

700    continue
c++++++++++++++++++++++++++++++++++++++++

               
c             STEP4
c     broadband  albedo determination

c         Step 4.1 planar BBA
                
c     planar BBA
                                     NSOLO=0
            x1=0.4425
            x2=0.70875
            x3=1.020
            
            y1=rp(3)
            y2=rp(11)
            y3=rp(21)
			write(*,*) "rp: ", y1, y2, y3
            
              d1=(y3-y1)*(x2-x1)-(y2-y1)*(x3-x1)
              d2=(x3*x3-x1*x1)*(x2-x1)-(x2*x2-x1*x1)*(x3-x1)
              
c           second order polynomial coefficients  for planar albedo:
              as=d1/d2
              bs=(y3-y2-as*(x3*x3-x2*x2))/(x3-x2)
              cs=y3-as*x3*x3-bs*x3

c           limits of integration
                   at=0.3
                   bt=2.4
                   aat=0.7
                   
c                          call qsimp(fun1,at,bt,ss1)
c                          call qsimp(fun2,at,bt,ss2)
   
c                           answer1=ss1/ss2
      
c                          call qsimp(fun1,at,aat,ss1)
c                         call qsimp(fun2,at,aat,ss2)
 
c                         answer2=ss1/ss2
      
                         call qsimp(fun1,aat,bt,ss1)
                      call qsimp(fun2,aat,bt,ss2)
 
                           answer3=ss1/ss2
			write(*,*) "answer3: ", answer3, ss1, ss2

      
c            Step 4.2 Spherical BBA
c     spherical BBA
                                  NSOLO=1
             x1=0.4425
             x2=0.70875
             x3=1.020
            
             y1=rs(3)
             y2=rs(11)
             y3=rs(21)
            
              d1=(y3-y1)*(x2-x1)-(y2-y1)*(x3-x1)
              d2=(x3*x3-x1*x1)*(x2-x1)-(x2*x2-x1*x1)*(x3-x1)
              
c           second order polynomial coefficients  for planar albedo:
              as=d1/d2
              bs=(y3-y2-as*(x3*x3-x2*x2))/(x3-x2)
              cs=y3-as*x3*x3-bs*x3

c           limits of integration
                   at=0.3
                   bt=2.4
                   aat=0.7
                   
c                           call qsimp(fun1,at,bt,ss1)
c                         call qsimp(fun2,at,bt,ss2)
   
c                            ans11=ss1/ss2
      
c             call qsimp(fun1,at,aat,ss1)
c             call qsimp(fun2,at,aat,ss2)

c             ans22=ss1/ss2
      
c                 call qsimp(fun1,aat,bt,ss1)
c                  call qsimp(fun2,aat,bt,ss2)
 
c            ans33=ss1/ss2
c            write(701,102) js,answer1,answer2,answer3,ans11,ans22,ans33
                  
1983          continue
200            continue
c+++++++++++++++++++++++++++++++++++++++++++++=               
               stop
               end

      
c          convolution of solar spectrum and spectral albedo:     
                 real function fun1(x)
c            Alex
                 real a(6),xa(168),ya(168),dif(168)
                 integer kss(2)
				 COMMON as,bs,cs,am1,am2,csza,cvza,r400,r865,r1020,NSOLO,
     +						thv, xa, ya
c     ALEX
               ks=1
c*****************************************
         
         pi=acos(-1.)
         
         alam3=0.865
         alam4=1.02
         akap3=2.4e-7
         akap4=2.25e-6

         
         alpha3=4.*pi*akap3/alam3
         alpha4=4.*pi*akap4/alam4
         
         eps=sqrt(alpha3/alpha4)
         ax1=1./(1.-eps)
         ax2=1./(1.-1./eps)
         
         r0=r865**ax1*r1020**ax2
         
         um1=3.*(1.+2.*csza)/7.
         um2=3.*(1.+2.*cvza)/7.
         
         xx=um1*um2/r0
         
         al=(alog(r1020/r0))**2./xx/xx/alpha4
         
c     ALEX: start of modification ( please, see lines above)
c    **** change of code N1****
          do 68 k=1,168
               
               dif(k)=abs(x-xa(k))              
               
               if ( dif(k).le.1.e-6) ks=k
               if ( dif(k).le.1.e-6) go to 69
68                continue
         
                 kss(1)=minloc(dif,1)
                 l=kss(1)
         
                 delta=x-xa(l)
                 if (delta.le.0) LS1=L-1
                 if (delta.le.0) LS2=L
                       if (delta.ge.0) LS1=L
                       if (delta.ge.0) LS2=L+1
            x0=xa(ls1)
            x1=xa(ls2)
            y0=ya(ls1)
            y1=ya(ls2)
c23456			
c			write(*,*) "XA: ", ls1, ls2, xa(ls1), ya(ls1),
c     +               xa(ls2), ya(ls2)  
        y=y0+(x-x0)*(y1-y0)/(x1-x0)
             go to 18
 69   continue
            y=ya(ks)
 18         continue

                 ASTRA=Y
c     end of modification ALEX

            
         dega=al*4.*pi*astra/x
         pow=sqrt(dega)
         if (pow.gt.1.e-6)rsd=exp(-pow)
         if (pow.le.1.e-6) rsd=1.
         
         if (NSOLO.eq.0) f1=rsd**um1
         if (NSOLO.eq.1) f1=rsd
         
          if (r400.le.thv.and.x.le.1.02) f1=as*x*x + bs*x +cs

        p0=32.38
        p1=-160140.33
        p2=7959.53

        t1=85.34*1.e-3
        t2=401.79*1.e-3
        funcs=p0+p1*exp(-x/t1)+p2*exp(-x/t2)
        if (x.lt.0.4) funcs=p0+p1*exp(-0.4/t1)+p2*exp(-0.4/t2)
     
        fun1=f1*funcs
c        if (x.ge.2.4) write(*,*) r400, x,f1,funcs,rsd,x,astra,nsolo
            return
            END

      
c     SOLAR LIGHT SPECTRUM AT THE GROUND:
      
           real function fun2(x)
      
               p0=32.38
               p1=-160140.33
               p2=7959.53
               t1=85.34*1.e-3
               t2=401.79*1.e-3       
               fun2=p0+p1*exp(-x/t1)+p2*exp(-x/t2)
               if (x.lt.0.4) fun2=p0+p1*exp(-0.4/t1)+p2*exp(-0.4/t2)
               
             return
             END

             SUBROUTINE qsimp(func,a,b,s)
              INTEGER JMAX
               REAL a,b,func,s,EPS
                 EXTERNAL func
                    PARAMETER (EPS=1.e-3, JMAX=20)

                       INTEGER j
                       REAL os,ost,st
		COMMON as,bs,cs,am1,am2,csza,cvza,r400,r865,r1020,NSOLO,
     +						thv, xa, ya				   
                        ost=-1.e30
                        os= -1.e30
      
                       do 11 j=1,JMAX
          
                            call trapzd(func,a,b,st,j)
                                s=(4.*st-ost)/3.
c23456
       write (*, *) "QSIMP: j, a, b, st, s, ost, os: ", 
     +     j, ", ", a, ", ", b, ", ", st, ", ", s, ", ", ost, ", ", os
                                  if (j.gt.5) then
          if (abs(s-os).lt.EPS*abs(os).or.(s.eq.0..and.os.eq.0.)) return
          endif
          os=s
          ost=st
11        continue
          
                    END
      
                      SUBROUTINE trapzd(func,a,b,s,n)
                       INTEGER n
                          REAL a,b,s,func
                          EXTERNAL func
                            INTEGER it,j
                            REAL del,sum,tnm,x
          COMMON as,bs,cs,am1,am2,csza,cvza,r400,r865,r1020,NSOLO,
     +						thv, xa, ya
				   
		write (*, *) "TRAPZ IN: n, a, b, s: ", 
     +                            n, ", ", a, ", ", b, ", ", s
                                if (n.eq.1) then
                            s=0.5*(b-a)*(func(a)+func(b))
		write (*, *) "TRAPZ FUNC: ", func(a), func(b) 
		write (*, *) "TRAPZ PARMS: as, bs, cs, csza, thv, r400: ", 
     +                as, ", ", bs, ", ", cs, ", ", csza, thv, r400
        write (*, *) "TRAPZ OUT: n, a, b, s: ", 
     +                                    n, ", ", a, ", ", b, ", ", s
                          else
                          it=2**(n-2)
                         tnm=it
                          del=(b-a)/tnm
                          x=a+0.5*del
                              sum=0.
                       do 11 j=1,it
                        sum=sum+func(x)
                       x=x+del
11                         continue
                              s=0.5*(s+(b-a)*sum/tnm)
        write (*, *) "TRAPZ O: n, a, b, s, del, x, sum: ", 
     +      n, ", ", a, ", ", b, ", ", s, ", ", del, ", ", x, ", ", sum
                     endif

                      END
   
     
            function falex1(am1,am2,co)
      
                pi=acos(-1.)
                a=1.247
                b=1.186
                c=5.157
                a1=0.087
                a2=0.014
                scat=acos(co)*180./pi
                p=11.1*exp(-a1*scat)+1.1*exp(-a2*scat)
                falex1=(a+b*(am1+am2)+c*am1*am2+p)/4./(am1+am2)
                
              return
              end
