 
      subroutine sclock
          use mod_precision


!   This routine starts the clock.

      implicit none

      real tarray(2),time,dtime

      time = dtime(tarray)

      return
      end


      subroutine pclock(tret)
          use mod_precision


!  This routine returns the CPU time elasped since the clock was last
!  examined.

      implicit none

      real tarray(2),time, dtime
      real(dp) :: tret

      time = dtime(tarray)

      tret = tarray(1)

      return
      end


