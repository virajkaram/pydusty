	program main

        character*100 fname,tline

        real lamf(100),tf(100)

        real lam(120000),flam(120000)

	print*,'enter filter file name '
        read(*,'(a)')fname
        open(unit=13,file=fname,form='formatted',status='old')
        read(13,*)nf
        do i=1,nf
          read(13,*)lamf(i),tf(i)
          lamf(i) = lamf(i)/1.e4
          enddo
        close(unit=13)
        print*,'filter ',fname,' runs from ',lamf(1),' to ',lamf(nf)

	print*,'enter spectrum name '
        read(*,'(a)')fname
        print*,' MARCS1 (1) or Kurucz (2) '
        read*,ilam
        open(unit=13,file=fname,form='formatted',status='old')
        if (ilam.eq.1) then
          open(unit=14,file='flx_wavelengths.vac',form='formatted',status='old')
        else 
          read(13,'(a)')tline
          endif
        i = 1
100     continue
        if (ilam.eq.1) then
          read(13,*,end=110)flam(i)
          read(14,*,end=110)lam(i)
          lam(i) = lam(i)/1.e4
        else
          read(13,*,end=110)lam(i),flam(i)
          flam(i) = flam(i)/lam(i)**2
          endif
          i = i + 1
          go to 100
110     close(unit=13)
        if (ilam.eq.1) close(unit=14)
        ns = i - 1

        bar0 = 0.0
        bar1 = 0.0
        bar2 = 0.0
        do i=1,ns
          call locate(lamf,nf,lam(i),j)
          tval = 0.0
          if ((j.ge.1).and.(j.lt.nf)) then
            tval = tf(j) + (tf(j+1)-tf(j))*(lam(i)-lamf(j))/(lamf(j+1)-lamf(j))
            endif
          if (i.ne.1) then
            bar0 = bar0 + 1.0
            bar1 = bar1 + 0.5*(lam(i)-lam(i-1))*(tval*lam(i)*       flam(i)+tvalo*lam(i-1)*         flam(i-1)) 
            bar2 = bar2 + 0.5*(lam(i)-lam(i-1))*(tval*lam(i)*lam(i)*flam(i)+tvalo*lam(i-1)*lam(i-1)*flam(i-1)) 
            endif
          enddo

        wavemean = bar2/bar1
        print*,wavemean

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
         
