# $Id: Makefile,v 1.9 2004/02/20 18:12:43 holtrop Exp $
#
CREATED_F=aao_rad.F
LIBNAMES  = $(MAIN_DIRECTORY)$(ADD_DEBUG) 
LIBNAMES  +=clasutil$(ADD_DEBUG) fputil$(ADD_DEBUG) recutl$(ADD_DEBUG) \
    bos$(ADD_DEBUG) fpack$(ADD_DEBUG) online_dummy$(ADD_DEBUG) \
    c_bos_io$(ADD_DEBUG)

LEPTO=$(shell ls -1 $(CERN)/$(CERN_LEVEL)/lib/liblepto*.a)
LEPTO_LIB=$(subst lib,,$(basename $(notdir $(LEPTO))))

LIBNAMES  +=$(LEPTO_LIB) jetset74 packlib mathlib kernlib $(LEPTO_LIB) $(RECSIS_MD_LIBS)
SHARED_LIBS=$(RECSIS_SHARED_LIBS)
include $(CLAS_CMS)/Makefile
#
# end of makefile
