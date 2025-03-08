#  Copyright (C) 2002 Regents of the University of Michigan,
#  portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# A new line to test the email notification


SHELL=/bin/sh

include ../Makefile.def
include ../Makefile.conf

install:
	cp ModBoostedFrame.f90 ${SHAREDIR}/.
	cd ${SHAREDIR};\
	perl -i'.BAK' -pe \
	's/(OBJECTS =)/OBJECTS = ModBoostedFrame.o /' Makefile
	cp Makefile.testboost ${DIR}/share/Library/test/.
	cp test_boosted_frame.f90 ${DIR}/share/Library/test/.
	cp test_mhd_x.ref ${DIR}/share/Library/test/.
	cd ${DIR}/share/Library/test/; \
	perl -i'.BAK' -pe \
	's/(test_geopack_exe:)/\ninclude Makefile.testboost\n\ntest_geopack_exe:/' Makefile
	cd ${DIR}; perl -i'.BAK' -pe \
	's/(TDSETUP:)/\ninclude Makefile.boost\n\nTDSETUP:/' Makefile
	cp Makefile.boost ${DIR}/.
	cp ModUserBoosted.f90 ${DIR}/GM/BATSRUS/srcUserExtra/.

uninstall:
	mv -f ${SHAREDIR}/ModBoostedFrame.f90 .
	cd ${SHAREDIR}; mv -f Makefile.BAK Makefile
	mv -f ${DIR}/share/Library/test/test_boosted_frame.f90 .
	mv -f ${DIR}/share/Library/test/Makefile.testboost .
	mv -f ${DIR}/share/Library/test/test_mhd_x.ref .
	cd ${DIR}/share/Library/test/; mv -f Makefile.BAK Makefile
	cd ${DIR}; mv -f Makefile.BAK Makefile
	mv -f ${DIR}/Makefile.boost .
	mv -f ${DIR}/GM/BATSRUS/srcUserExtra/ModUserBoosted.f90 .

clean: cleanfiles

distclean: clean
