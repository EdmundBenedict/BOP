 
   subroutine bldclus(ia)
      use mod_precision
      use mod_all_scalar
      use mod_const
      use mod_ham
!
!    This is a routine to build the recursion cluster for one atom.
!

      implicit none

      include "Include/NebList.array"

      integer ia,ib
      integer ja,jn
      integer ptr
      integer n

! pcluster: level -> position in cluster array/ hamiltonian
! cluster: position in the local hamiltonian for atom ia -> atom (1 is ia) 
! decipher: the reverse of cluster; atom -> position in the hamiltonian


!    Set up zeroth shell.

      pcluster(0) = 1  
      ptr = 2
      pcluster(1) = ptr 
      
      cluster(1) = ia
      decipher(ia) = 1
      


!    Build remaining shells.
      do n = 1,nbase

      
!       For all atoms in previous shell do ...
         do jn = pcluster(n-1), pcluster(n)-1

!            print *, 'jn',jn
!          For all neighbors of these atoms do ...
            ja = aptr(cluster(jn))
!             print *, 'ja,cluster(jn)',ja,cluster(jn)
            ib = bptr(ja)
            do while(ib /= eol)
               

!
!             If it is not found in an earlier shell, add it to the
!              current shell.
!

               if (decipher(ib) == 0) then
                  cluster(ptr) = ib
                  decipher(ib) = ptr

                  ptr = ptr + 1

                  if (ptr > mxcls) then
                     write(6,'(''Cluster array overflow. Increase MXCLS.'')')
                     call panic()
                  endif
               endif

               ja = ja + 1
               ib = bptr(ja)
            enddo

         enddo

!
!       Mark start of next shell.
!

         pcluster(n+1) = ptr

      enddo

      
!       print *, ''
!       print *, ''
!       print *, 'ia  :', ia
!       print *, 'pcls:', pcluster(:3)
!       print *, 'cls :', cluster(pcluster(0):pcluster(2))
!       print *, 'bptr:', bptr(aptr(ia):aptr(ia)+pcluster(2)-pcluster(0))
      
     
! ! 
!       print *,'clussize:',ptr
! 
! 
!       write(500+ia,'("aptr:")')
!       write(500+ia,*) aptr
!       
!       write(500+ia,'("bptr:")')
!       write(500+ia,*) bptr
!       
!       write(500+ia,'("decipher:")')
!       write(500+ia,*) decipher
!       
!       write(500+ia,'("cluster:")')
!       write(500+ia,*) cluster
!       
!       write(500+ia,'("pcluster:")')
!       write(500+ia,*) pcluster
!       
!       
! !
!    Re-initialize array.
!

!       do n=1,totnd
!         decipher(n) = 0
!       enddo
       decipher(1:totnd) = 0
!      WRITE(11,'(15I5)') (PCLUSTER(N),N=0,NBASE+2,1)


      end subroutine bldclus

