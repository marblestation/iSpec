# for gfortran 
F = gfortran
Loader  = gfortran

SOURCE = source

include1 = ../$(SOURCE)/spectrum.inc
include2 = ../$(SOURCE)/tsuji.par

#LINKOPT = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -fconvert=big-endian -O3 -fno-automatic -mcmodel=small -std=legacy -Wall -Wsurprising

#LINKOPT =  -fconvert=big-endian -O3 -fno-automatic -Wall -Wsurprising -ffpe-trap=zero,invalid -fbacktrace -std=legacy

LINKOPT =  -fconvert=big-endian -O3 -fno-automatic -fbacktrace -std=legacy -mcmodel=small
# -fno-sign-zero
#-ffpe-trap=zero,invalid 

babsma_luobj =  \
      ../$(SOURCE)/absko.f \
      ../$(SOURCE)/nhop.f \
      ../$(SOURCE)/babsma.f \
      ../$(SOURCE)/detabs.f \
      ../$(SOURCE)/getinfospecies.f \
      ../$(SOURCE)/inabs.f \
      ../$(SOURCE)/injon.f \
      ../$(SOURCE)/input.f \
      ../$(SOURCE)/jon.f \
      ../$(SOURCE)/lenstr.f \
      ../$(SOURCE)/qas.f \
      ../$(SOURCE)/qtrav.f \
      ../$(SOURCE)/rossop.f \
      ../$(SOURCE)/tabs.f \
      ../$(SOURCE)/tryck.f \
      ../$(SOURCE)/zerof.f \
      ../$(SOURCE)/heopac.f \
      ../$(SOURCE)/h2opac.f \
      ../$(SOURCE)/CIAh2h2.f \
      ../$(SOURCE)/CIAh2h.f \
      ../$(SOURCE)/CIAhhe.f \
      ../$(SOURCE)/CIAh2he.f \
      ../$(SOURCE)/heborysowopac.f \
      ../$(SOURCE)/h2borysowopac.f \
      ../$(SOURCE)/hbop.f \
      ../$(SOURCE)/hlinprof.f \
      ../$(SOURCE)/eqmol_pe_lu.f \
      ../$(SOURCE)/die_pe_lu.f \
      ../$(SOURCE)/funcv.f \
      ../$(SOURCE)/moleq.f \
      ../$(SOURCE)/takemolec.f \
      ../$(SOURCE)/makeabund.f \
      ../$(SOURCE)/getisotopmass.f \
      ../$(SOURCE)/molmat.f \
      ../$(SOURCE)/partf.f \
      ../$(SOURCE)/pemake.f \
      ../$(SOURCE)/pemakero.f \
      ../$(SOURCE)/rinteg.f \
      ../$(SOURCE)/molecpartf.f \
      ../$(SOURCE)/ludcmp.f \
      ../$(SOURCE)/lubksb.f \
      ../$(SOURCE)/lowercase_conversion.f90 \
      ../$(SOURCE)/molfys.f 

eqwidt_luobj = \
      ../$(SOURCE)/anstee.f \
      ../$(SOURCE)/bpl.f \
      ../$(SOURCE)/eqwidtb.f \
      ../$(SOURCE)/eqwidtbplatt.f \
      ../$(SOURCE)/eqwidt.f \
      ../$(SOURCE)/clock-gf.f \
      ../$(SOURCE)/depth.f \
      ../$(SOURCE)/gausi.f \
      ../$(SOURCE)/getinfospecies.f \
      ../$(SOURCE)/getlele.f \
      ../$(SOURCE)/hydropac.f \
      ../$(SOURCE)/hbop.f \
      ../$(SOURCE)/hlinprof.f \
      ../$(SOURCE)/inp3.f \
      ../$(SOURCE)/input.f \
      ../$(SOURCE)/lenstr.f \
      ../$(SOURCE)/eqmol_pe_lu.f \
      ../$(SOURCE)/die_pe_lu.f \
      ../$(SOURCE)/funcv.f \
      ../$(SOURCE)/ludcmp.f \
      ../$(SOURCE)/lubksb.f \
      ../$(SOURCE)/makeabund.f \
      ../$(SOURCE)/getisotopmass.f \
      ../$(SOURCE)/molfys.f \
      ../$(SOURCE)/takemolec.f \
      ../$(SOURCE)/readmo.f \
      ../$(SOURCE)/readfilt.f \
      ../$(SOURCE)/tint.f \
      ../$(SOURCE)/lint.f \
      ../$(SOURCE)/traneq.f \
      ../$(SOURCE)/traneqplatt.f \
      ../$(SOURCE)/tranw.f \
      ../$(SOURCE)/newvoigt.f \
      ../$(SOURCE)/molmat.f \
      ../$(SOURCE)/tranfr.f \
      ../$(SOURCE)/tranfrplatt.f \
      ../$(SOURCE)/trrays.f \
      ../$(SOURCE)/trqua2.f \
      ../$(SOURCE)/transc.f \
      ../$(SOURCE)/transcplatt.f \
      ../$(SOURCE)/partf.f \
      ../$(SOURCE)/molecpartf.f \
      ../$(SOURCE)/partffordepth.f \
      ../$(SOURCE)/h2opartf.f \
      ../$(SOURCE)/h2opartf2001.f \
      ../$(SOURCE)/partfamesh2o.f \
      ../$(SOURCE)/partfBarberH2O.f \
      ../$(SOURCE)/weightlimb.f 


bsyn_luobj = \
      ../$(SOURCE)/anstee.f \
      ../$(SOURCE)/bpl.f \
      ../$(SOURCE)/bsyn.f \
      ../$(SOURCE)/bsynb.f \
      ../$(SOURCE)/bsynbplatt.f \
      ../$(SOURCE)/clock-gf.f \
      ../$(SOURCE)/cubint_module.f90 \
      ../$(SOURCE)/depth.f \
      ../$(SOURCE)/gausi.f \
      ../$(SOURCE)/getinfospecies.f \
      ../$(SOURCE)/getlele.f \
      ../$(SOURCE)/hydropac.f \
      ../$(SOURCE)/hbop.f \
      ../$(SOURCE)/hlinprof.f \
      ../$(SOURCE)/inp3.f \
      ../$(SOURCE)/input.f \
      ../$(SOURCE)/Iplus_calc.f \
      ../$(SOURCE)/lenstr.f \
      ../$(SOURCE)/eqmol_pe_lu.f \
      ../$(SOURCE)/die_pe_lu.f \
      ../$(SOURCE)/funcv.f \
      ../$(SOURCE)/ludcmp.f \
      ../$(SOURCE)/lubksb.f \
      ../$(SOURCE)/makeabund.f \
      ../$(SOURCE)/getisotopmass.f \
      ../$(SOURCE)/molfys.f \
      ../$(SOURCE)/takemolec.f \
      ../$(SOURCE)/readmo.f \
      ../$(SOURCE)/readfilt.f \
      ../$(SOURCE)/read_segments.f \
      ../$(SOURCE)/read_modelatom.f \
      ../$(SOURCE)/read_departure.f \
      ../$(SOURCE)/read_nlteinfofile.f \
      ../$(SOURCE)/strings_module.f90 \
      ../$(SOURCE)/tint.f \
      ../$(SOURCE)/lint.f \
      ../$(SOURCE)/traneq.f \
      ../$(SOURCE)/traneqplatt.f \
      ../$(SOURCE)/tranw.f \
      ../$(SOURCE)/newvoigt.f \
      ../$(SOURCE)/molmat.f \
      ../$(SOURCE)/tranfr.f \
      ../$(SOURCE)/tranfrplatt.f \
      ../$(SOURCE)/trrays.f \
      ../$(SOURCE)/trqua2.f \
      ../$(SOURCE)/transc.f \
      ../$(SOURCE)/transcplatt.f \
      ../$(SOURCE)/partf.f \
      ../$(SOURCE)/molecpartf.f \
      ../$(SOURCE)/partffordepth.f \
      ../$(SOURCE)/h2opartf.f \
      ../$(SOURCE)/h2opartf2001.f \
      ../$(SOURCE)/partfamesh2o.f \
      ../$(SOURCE)/partfBarberH2O.f \
      ../$(SOURCE)/weightlimb.f 


#clean : 
#	rm *.o bsyn_lu babsma_lu eqwidt_lu 

all : cubint_module.mod strings_module.mod babsma_lu bsyn_lu
#  eqwidt_lu 

cubint_module.o : ../${SOURCE}/cubint_module.f90
	$(F)  -c ../${SOURCE}/cubint_module.f90

cubint_module.mod : ../$(SOURCE)/cubint_module.f90 cubint_module.o
	$(F)  -c ../$(SOURCE)/cubint_module.f90

strings_module.o : ../${SOURCE}/strings_module.f90
	 $(F)  -c ../${SOURCE}/strings_module.f90

strings_module.mod : ../$(SOURCE)/strings_module.f90 strings_module.o
	 $(F)  -c ../$(SOURCE)/strings_module.f90

babsma_lu : $(babsma_luobj) $(include1) $(include2)
	$(F) $(LINKOPT) -o babsma_lu $(babsma_luobj)

bsyn_lu: $(bsyn_luobj) $(include1) $(include2)
	$(F) $(LINKOPT) -o bsyn_lu  $(bsyn_luobj)
 
eqwidt_lu: $(eqwidt_luobj) $(include1) $(include2)
	$(F) $(LINKOPT) -o eqwidt_lu  $(eqwidt_luobj)
 
