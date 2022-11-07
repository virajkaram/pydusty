        program main

        real lam(120000),flx(120000)
        real flam(10000),filt(10000) 

        character*30 flist(10)

        real lcen(10)

        real lwant,lmin,lmax

        print*,'enter av '
        read*,av
        rv = 3.1

        open(unit=13,file='foo.dat',form='formatted',status='old')
        open(unit=14,file='flx_wavelengths.vac',form='formatted',status='old')
c'
        i = 1
100     read(13,*,end=110)flx(i)
          read(14,*)lam(i)
          x      = 1.e4/lam(i)
	  rval   = rv*rl(rv,x)
          aval   = av*rval/3.1
          flx(i) = flx(i)*10.0**(-0.4*aval) 
          i    = i + 1
          go to 100
110     close(unit=13) 
        close(unit=14)
        ns = i - 1
        print*,'read ns ',ns

        nf = 8
        flist(1) = 'HST_ACS_WFC.F555W_81.dat'
        flist(2) = 'HST_ACS_WFC.F814W_81.dat'
        flist(3) = '2MASS_2MASS.J.dat'
        flist(4) = '2MASS_2MASS.H.dat'
        flist(5) = 'Spitzer_IRAC.I1.dat'
        flist(6) = 'Spitzer_IRAC.I2.dat'
        flist(7) = 'Spitzer_IRAC.I3.dat'
        flist(8) = 'Spitzer_IRAC.I4.dat'
        lcen(1)  = 0.550
        lcen(2)  = 0.806
        lcen(3)  = 1.235
        lcen(4)  = 1.649
        lcen(5)  = 3.510
        lcen(6)  = 4.440
        lcen(7)  = 5.630
        lcen(8)  = 7.590

        do k=1,nf

        open(unit=13,file=flist(k),form='formatted',status='old')
        i = 1
200     read(13,*,end=210)flam(i),filt(i)
          i = i + 1
          go to 200
210     close(unit=13)
        nf = i - 1
        print*,'read ns ',nf
        lwant = lcen(k)*1.e4
        lmin  = 0.9*lwant
        lmax  = 1.1*lwant

        bar0a = 0.0
        bar1a = 0.0
        bar2a = 0.0
        bar0b = 0.0
        bar1b = 0.0
        barsb = 0.0
        do i=1,ns
          call locate(flam,nf,lam(i),j)
          if ((j.ge.1).and.(j.lt.nf)) then
            fval  = filt(j) + (filt(j+1)-filt(j))*(lam(i)-flam(j))/(flam(j+1)-flam(j))
            bar0a = bar0a + lam(i)*fval
            bar1a = bar1a + lam(i)*fval*flx(i)
            bar2a = bar2a + lam(i)*lam(i)*fval*flx(i)
            endif
          if ((lam(i).ge.lmin).and.(lam(i).le.lmax)) then
            bar0b = bar0b + lam(i)
            bar1b = bar1b + lam(i)*flx(i)
            bar2b = bar2b + lam(i)*lam(i)*flx(i)
            endif
          enddo
        bar2a = bar2a/bar1a
        bar2b = bar2b/bar1b

        bar1a = bar1a/bar0a
        bar1b = bar1b/bar0b
       
        print*,'  '
        print*,'real filter ',bar1a,bar2a
        print*,'fake filter ',bar1b,bar2b
        print*,' fake/real  ',bar1b/bar1a,bar2b/bar2a
        
        enddo

        end

      SUBROUTINE locate(xx,n,x,j)
      INTEGER j,n
      REAL x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      END

	function rl(rv,x)


          if (x.lt.1.1) then
            a =  0.574*x**1.61
            b = -0.527*x**1.61
          else
            if (x.lt.3.3) then
              y = x-1.82
              a = 1.0+0.17699*y-0.50447*y*y-0.02427*y**3+0.72085*y**4 +
     1            0.01979*y**5-0.77530*y**6+0.32999*y**7
              b = 1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-
     1            0.62251*y**5+5.30260*y**6-2.09002*y**7
            else
              if (x.lt.5.9) then
                fa = 0.0
                fb = 0.0
              else
                fa = -0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3
                fb =  0.21300*(x-5.9)**2 + 0.120700*(x-5.9)**3
                endif
              a =  1.752 - 0.316*x - 0.104/((x-4.67)**2+0.341) + fa
              b = -3.090 + 1.825*x + 1.206/((x-4.67)**2+0.263) + fb
              endif
            endif

          rl = a+b/rv 

          return
          end
