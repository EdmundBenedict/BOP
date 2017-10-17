 
      subroutine sumrule()
          use mod_precision


          use mod_all_scalar

          use mod_const

          use ab_io

!
!    This is a subroutine to check the sum rule on the recursion
!     coefficients.
!

      implicit none

!
!    Include constants.
!

!      include "Include/ALL.const"

!
!    Include scalars.
!

!      include "Include/ALL.scalar"

!
!    Include arrays.
!

      include "Include/Atom.array"
!      include "Include/Moment.array"
      include "Include/NebList.array"
!      include "Include/RecurCoef.array"

!
!    Declare the simple variables.
!

      real(dp) :: hopint,suma,sumb
      real(dp) :: rulea,ruleb

      integer nma,wrflag
      integer n
      integer nmax
      integer nla,nlb
      integer ia,ib,la,lb,lla
      integer nstta,nsttb
      integer za,zb
      integer ja,ja0,jla,jla0

!
!    Calculate the sum.
!

      wrflag = 1

      if (momflg == 1) then

         do ia = 1,nd,1

            call rdab(wrflag,ia)

            write(9,'(/''Atom # '',I4)') ia

            za = z(ia)
            call states(za,nla,nstta,llista)
            jla = 1
            do la = 1,nla,1
               do lla = 1,2*llista(la)+1,1
                  mlista(jla) = la
                  jla = jla + 1
               enddo
            enddo

            jla0 = 1
            do la = 1,nla,1
               nma = 2*llista(la)+1
               nmax = lchain(la)

               write(9,'(''L = '',I1)') la

               do n = 0,nmax+1,1
                  ja = aptr(ia)
                  ja0 = ja
                  suma = 0.0_dp
                  sumb = 0.0_dp
                  do while (bptr(ja) /= eol)
                     ib = bptr(ja)
                     zb = z(ib)
                     call states(zb,nlb,nsttb,llistb)
                     if (ia /= ib) then
                        jla = jla0
                        do lla = 1,nma,1
                           do lb = 1,nsttb,1
                              hopint = darec(0,jla,lb,ja-ja0+1)
                              suma = suma + hopint* & 
     &                                      darec(n,jla,lb,ja-ja0+1)
                              sumb = sumb + hopint* & 
     &                                      dbrec(n,jla,lb,ja-ja0+1)
                           enddo
                           jla = jla + 1
                        enddo
                     endif
                     ja = ja + 1
                  enddo
                  suma = suma/real(nma, dp)
                  sumb = sumb/real(nma, dp)
                  if (n == nmax+1) then
                     rulea = 0.0_dp
                  elseif (n == nmax) then
                     if (term == 1) then
                        rulea = lbinf(la)**2-brec(n,la)**2
                     elseif ((term == 2).or.(term == 3)) then
                        rulea = -brec(n,la)**2
                     endif
                  else
                     rulea = brec(n+1,la)**2-brec(n,la)**2
                  endif
                  if (n == 0) then
                     ruleb = 0.0_dp
                  elseif (n <= nmax) then
                     ruleb = 0.5_dp*brec(n,la)*(arec(n,la) & 
     &                                        -arec(n-1,la))
                  else
                     ruleb = 0.5_dp*lbinf(la)*(lainf(la) & 
     &                                        -arec(nmax,la))
                  endif
                  write(9,'(''N = '',I2,'' SUMA  = '',G22.15, & 
     &                      '' SUMB  = '',G22.15)') & 
     &                       n,suma,sumb
                  write(9,'(''       RULEA = '',G22.15, & 
     &                      '' RULEB = '',G22.15)') rulea,ruleb
               enddo

               jla0 = jla0 + nma

            enddo

         enddo

      elseif (momflg == 2) then

         do ia = 1,nd,1

            call rdab(wrflag,ia)

            write(9,'(/''Atom # '',I4)') ia

            za = z(ia)
            call states(za,nla,nstta,llista)

            do la = 1,nstta,1
               nmax = lchain(la)

               write(9,'(''L = '',I1)') la

               do n = 0,nmax+1,1
                  ja = aptr(ia)
                  ja0 = ja
                  suma = 0.0_dp
                  sumb = 0.0_dp
                  do while (bptr(ja) /= eol)
                     ib = bptr(ja)
                     zb = z(ib)
                     call states(zb,nlb,nsttb,llistb)
                     do lb = 1,nsttb,1
                        if ((ia /= ib).or.(la /= lb)) then
                           hopint = darec(0,la,lb,ja-ja0+1)
                           suma = suma + hopint* & 
     &                                   darec(n,la,lb,ja-ja0+1)
                           sumb = sumb + hopint* & 
     &                                   dbrec(n,la,lb,ja-ja0+1)
                        endif
                     enddo
                     ja = ja + 1
                  enddo
                  if (n == nmax+1) then
                     rulea = 0.0_dp
                  elseif (n == nmax) then
                     if (term == 1) then
                        rulea = lbinf(la)**2-brec(n,la)**2
                     elseif ((term == 2).or.(term == 3)) then
                        rulea = -brec(n,la)**2
                     endif
                  else
                     rulea = brec(n+1,la)**2-brec(n,la)**2
                  endif
                  if (n == 0) then
                     ruleb = 0.0_dp
                  elseif (n <= nmax) then
                     ruleb = 0.5_dp*brec(n,la)*(arec(n,la) & 
     &                                        -arec(n-1,la))
                  else
                     ruleb = 0.5_dp*lbinf(la)*(lainf(la) & 
     &                                        -arec(nmax,la))
                  endif
                  write(9,'(''N = '',I2,'' SUMA  = '',G22.15, & 
     &                      '' SUMB  = '',G22.15)') & 
     &                       n,suma,sumb
                  write(9,'(''       RULEA = '',G22.15, & 
     &                      '' RULEB = '',G22.15)') rulea,ruleb
               enddo

            enddo

         enddo

      elseif (momflg == 3) then

         do ia = 1,nd,1

            call rdab(wrflag,ia)

            write(9,'(/''Atom # '',I4)') ia

            za = z(ia)
            call states(za,nla,nstta,llista)

            do la = 1,nstta,1
               nmax = lchain(la)

               write(9,'(''L = '',I1)') la

               do n = 0,nmax+1,1
                  ja = aptr(ia)
                  ja0 = ja
                  suma = 0.0_dp
                  sumb = 0.0_dp
                  do while (bptr(ja) /= eol)
                     ib = bptr(ja)
                     zb = z(ib)
                     call states(zb,nlb,nsttb,llistb)
                     if (ia /= ib) then
                        do lb = 1,nsttb,1
                           hopint = darec(0,la,lb,ja-ja0+1)
                           suma = suma + hopint* & 
     &                                   darec(n,la,lb,ja-ja0+1)
                           sumb = sumb + hopint* & 
     &                                   dbrec(n,la,lb,ja-ja0+1)
                        enddo
                     endif
                     ja = ja + 1
                  enddo
                  if (n == nmax+1) then
                     rulea = 0.0_dp
                  elseif (n == nmax) then
                     if (term == 1) then
                        rulea = lbinf(la)**2-brec(n,la)**2
                     elseif ((term == 2).or.(term == 3)) then
                        rulea = -brec(n,la)**2
                     endif
                  else
                     rulea = brec(n+1,la)**2-brec(n,la)**2
                  endif
                  if (n == 0) then
                     ruleb = 0.0_dp
                  elseif (n <= nmax) then
                     ruleb = 0.5_dp*brec(n,la)*(arec(n,la) & 
     &                                        -arec(n-1,la))
                  else
                     ruleb = 0.5_dp*lbinf(la)*(lainf(la) & 
     &                                        -arec(nmax,la))
                  endif
                  write(9,'(''N = '',I2,'' SUMA  = '',G22.15, & 
     &                      '' SUMB  = '',G22.15)') & 
     &                       n,suma,sumb
                  write(9,'(''       RULEA = '',G22.15, & 
     &                      '' RULEB = '',G22.15)') rulea,ruleb
               enddo

            enddo

         enddo

      endif

      close(50)

      end

