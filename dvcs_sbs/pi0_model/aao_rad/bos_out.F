      subroutine bos_out
      implicit none
      
#include "bcs.inc"
#include "mc.inc"
     
      integer nbank,ind,ierr,indx,j,output_counter
      logical do_mctk
      
      SAVE output_counter
 
      do_mctk = .false.
      
      output_counter	= output_counter+1
      
      ind		= nbank('HEAD',0,8,1)
      iw(ind+1)		= 1                 	! version
      iw(ind+2)		= 1                 	! run
      iw(ind+3)		= output_counter     	! event number
      iw(ind+4)		= 100                	! event time
      iw(ind+5) 	= -4                 	! event type
      iw(ind+6) 	= 0                 	! ROC
      iw(ind+7) 	= 15                 	! event class
      iw(ind+8)     = 1                  	! weight

      ind		= nbank('TGBI',0,4,1)
      if (hel.gt.0.5) then
        iw(ind+1) = 32768
      else
        iw(ind+1) = 0 
      endif
         
      ind		= nbank('PART',0,13,npart) 
      do j=1,npart
         indx 		= ind+(j-1)*13
         iw(indx+1)	= id(j)			! particle ID (GEANT)
         rw(indx+2)	= v(j,1)		! x vertex position
         rw(indx+3)	= v(j,2)		! y vertex position
         rw(indx+4)	= v(j,3)		! z vertex position
         rw(indx+5)	= p(j,4)		! energy 
         rw(indx+6)	= p(j,1)		! px
         rw(indx+7)	= p(j,2)		! px
         rw(indx+8)	= p(j,3)		! px
         rw(indx+9)	= q(j)			! charge
         iw(indx+10) 	= j			! Track pointer
         rw(indx+11)	= xsec
         rw(indx+12) 	= asym
         iw(indx+13) 	= 1
      enddo
      
      if (do_mctk) then

      ind		= nbank('MCTK',0,11,npart) 
      do j=1,npart
         indx 		= ind+(j-1)*11
         rw(indx+1)	= cs(j,1)		! x dir cosine at track origin
         rw(indx+2)	= cs(j,2)		! y dir cosine
         rw(indx+3)	= cs(j,3)		! z dir cosine
         rw(indx+4)	= ptot(j)		! momentum
         rw(indx+5)	= rm(j) 		! mass 
         rw(indx+6)	= q(j)  		! charge
         iw(indx+7)	= pdgid(j)		! PDG ID
         iw(indx+8)	= 1			    ! track flag
         iw(indx+9)	= 1  			! beginning vertex number
         iw(indx+10)	= 1			! ending vertex number
         iw(indx+11) 	= 1			! parent track
      enddo
      
      ind		= nbank('MCVX',0,5,npart) 
      do j=1,npart
         indx 		= ind+(j-1)*4
         rw(indx+1)	= v(j,1)		! x of vertex
         rw(indx+2)	= v(j,2)		! y of vertex
         rw(indx+3)	= v(j,3)		! z of vertex
         rw(indx+4)	= 1.			! secs of flight
         iw(indx+5)	= 1.			! vertex flag
      enddo
      
      endif
      
      call fwbos(iw,33,'E',ierr)	! write banks to file
      call bdrop(iw,'E')     		! drop the bank, reclaim the space
      call bgarb(iw)	     		! garbage collection
         
      return
      END
