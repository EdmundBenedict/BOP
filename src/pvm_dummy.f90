 
!   This file will contain DUMMY routines so that the BOP package
!   can be linked without PVM libraries.

      subroutine initcomms(me,parflg)
          use mod_precision

      
      integer me, parflg

      me = 1
      parflg = 0

      return
      end


      subroutine send_dp(var,len,tid)

      return
      end


      subroutine send_int(var,len,tid)

      return
      end


      subroutine send_str(var,len,tid)

      return
      end


      subroutine cast_dp(var,len,np,tid)
          use mod_precision


      return
      end


      subroutine cast_int(var,len,np,tid)
          use mod_precision


      return
      end


      subroutine cast_str(var,len,np,tid)
          use mod_precision


      return
      end


      subroutine recv_dp(var,len,tid)
          use mod_precision


      return
      end


      subroutine recv_int(var,len,tid)
          use mod_precision


      return
      end


      subroutine recv_str(var,len,tid)
          use mod_precision


      return
      end

      subroutine quitpvm()
          use mod_precision


      return
      end

      subroutine findhere(i)
          use mod_precision

      integer i

      return
      end

      subroutine kdivide()
          use mod_precision


      return
      end

      subroutine reportload(i)
          use mod_precision

      integer i

      return
      end

      subroutine gandsdp(i,z)
          use mod_precision


      integer i
      real(dp) :: z

      return
      end

      subroutine fbound(f)
          use mod_precision


      real(dp) :: f(*)
      
      return
      end

      subroutine pullfq()
          use mod_precision


      return
      end


      subroutine findmin(var)
          use mod_precision


      return
      end

      subroutine findmax(var)
          use mod_precision

      
      return
      end

      subroutine ftomaster(f,n)
          use mod_precision

      
      integer n
      real(dp) :: f(*)

      return
      end

      subroutine rtoworker(n)
          use mod_precision

      integer n

      return
      end

      subroutine kforcacc()
          use mod_precision


      return
      end

      subroutine kdqacc()
          use mod_precision


      return
      end
