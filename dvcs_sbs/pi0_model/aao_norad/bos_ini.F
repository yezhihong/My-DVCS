      subroutine bos_ini(recname)
      implicit none
      
      character*8 recname
      character*132 ddl_file

#include "bcs.inc"
#include "names.inc"

      CALL bnames(1000)
      call bos(iw,nbcs)

      call revinm ('CLAS_PARMS','clasbanks.ddl',ddl_file)
      CALL txt2bos(ddl_file)

      close (33)

      call fparm ('OPEN ' //recname// 'UNIT=33 '// 
     &     'FILE="'//bosout//'"'// 
     &     ' RECL=32760 ACTION=WRITE STATUS=NEW FORM=BINARY')
     
      RETURN
      END
