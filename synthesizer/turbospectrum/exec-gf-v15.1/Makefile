# gfortran Makefile

babsma_luobj =  \
      absko.o \
      babsma.o \
      detabs.o \
      getinfospecies.o \
      inabs.o \
      injon.o \
      input.o \
      jon.o \
      lenstr.o \
      qas.o \
      qtrav.o \
      rossop.o \
      tabs.o \
      tryck.o \
      zerof.o \
      heopac.o \
      h2opac.o \
      CIAh2h2.o \
      CIAh2h.o \
      CIAhhe.o \
      CIAh2he.o \
      heborysowopac.o \
      h2borysowopac.o \
      eqmol_pe_lu.o \
      die_pe_lu.o \
      funcv.o \
      moleq.o \
      takemolec.o \
      makeabund.o \
      getisotopmass.o \
      molmat.o \
      partf.o \
      pemake.o \
      pemakero.o \
      rinteg.o \
      molecpartf.o \
      ludcmp.o \
      lubksb.o \
      molfys.o 

eqwidt_luobj = \
      anstee.o \
      bpl.o \
      eqwidtb.o \
      eqwidtbplatt.o \
      eqwidt.o \
      clock-gf.o \
      depth.o \
      gausi.o \
      getinfospecies.o \
      getlele.o \
      hydropac.o \
      hlinop.o \
      inp3.o \
      input.o \
      lenstr.o \
      eqmol_pe_lu.o \
      die_pe_lu.o \
      funcv.o \
      ludcmp.o \
      lubksb.o \
      makeabund.o \
      getisotopmass.o \
      molfys.o \
      takemolec.o \
      readmo.o \
      readfilt.o \
      tint.o \
      lint.o \
      traneq.o \
      traneqplatt.o \
      tranw.o \
      newvoigt.o \
      molmat.o \
      tranfr.o \
      tranfrplatt.o \
      trrays.o \
      trqua2.o \
      transc.o \
      transcplatt.o \
      partf.o \
      molecpartf.o \
      partffordepth.o \
      h2opartf.o \
      h2opartf2001.o \
      partfamesh2o.o \
      partfBarberH2O.o \
      weightlimb.o 


bsyn_luobj = \
      anstee.o \
      bpl.o \
      bsyn.o \
      bsynb.o \
      bsynbplatt.o \
      clock-gf.o \
      depth.o \
      gausi.o \
      getinfospecies.o \
      getlele.o \
      hydropac.o \
      hlinop.o \
      inp3.o \
      input.o \
      Iplus_calc.o \
      lenstr.o \
      eqmol_pe_lu.o \
      die_pe_lu.o \
      funcv.o \
      ludcmp.o \
      lubksb.o \
      makeabund.o \
      getisotopmass.o \
      molfys.o \
      takemolec.o \
      readmo.o \
      readfilt.o \
      tint.o \
      lint.o \
      traneq.o \
      traneqplatt.o \
      tranw.o \
      newvoigt.o \
      molmat.o \
      tranfr.o \
      tranfrplatt.o \
      trrays.o \
      trqua2.o \
      transc.o \
      transcplatt.o \
      partf.o \
      molecpartf.o \
      partffordepth.o \
      h2opartf.o \
      h2opartf2001.o \
      partfamesh2o.o \
      partfBarberH2O.o \
      weightlimb.o 

# For gfortran

F = gfortran
Loader = gfortran
OPT1 = -L/usr/lib -fconvert=big-endian -O3 -fno-automatic
LINKOPT = -L/usr/lib -fconvert=big-endian -O3 -fno-automatic

#SOURCE = source-v15.1
SOURCE = src


#OPT1     = -B/usr/lib -convert big_endian -fp-model precise -ip -ipo  -O3 -fpe3 -ftz -traceback -fp-speculationoff -opt-report  -save  -static-intel
#LINKOPT        = -B/usr/lib -convert big_endian -fp-model precise -ip -ipo  -O3 -fpe3 -ftz -traceback -fp-speculationoff -opt-report  -save -static-intel

#default below
#OPT1     = -B/usr/lib -convert big_endian -fp-model precise   -O3 -fpe3 -ftz -traceback -fp-speculationoff -opt-report  -save  -static-intel
#LINKOPT        = -B/usr/lib -convert big_endian -fp-model precise   -O3 -fpe3 -ftz -traceback -fp-speculationoff -opt-report  -save -static-intel

#OPT1 = -convert big_endian -O0 -g -C -traceback -save
#LINKOPT = -convert big_endian -O0 -g -C -traceback -save

#clean : 
#	rm *.o bsyn_lu babsma_lu eqwidt_lu 

all : babsma_lu bsyn_lu eqwidt_lu 

babsma_lu : $(babsma_luobj)
	$(F) $(LINKOPT) -o babsma_lu $(babsma_luobj)

bsyn_lu: $(bsyn_luobj)
	$(F) $(LINKOPT) -o bsyn_lu  $(bsyn_luobj)
 
eqwidt_lu: $(eqwidt_luobj)
	$(F) $(LINKOPT) -o eqwidt_lu  $(eqwidt_luobj)
 
 
absko.o : ../${SOURCE}/absko.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/absko.f
 
anstee.o : ../${SOURCE}/anstee.f 
	$(F) -c $(OPT1) ../${SOURCE}/anstee.f
 
babsma.o : ../${SOURCE}/babsma.f ../${SOURCE}/spectrum.inc ../${SOURCE}/tsuji.par
	$(F) -c $(OPT1) ../${SOURCE}/babsma.f
 
babsmaFIR.o : ../${SOURCE}/babsmaFIR.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/babsmaFIR.f
 
babsmaIR.o : ../${SOURCE}/babsmaIR.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/babsmaIR.f
 
bpl.o : ../${SOURCE}/bpl.f
	$(F) -c $(OPT1) ../${SOURCE}/bpl.f
 
bsyn.o : ../${SOURCE}/bsyn.f ../${SOURCE}/spectrum.inc ../${SOURCE}/tsuji.par
	$(F) -c $(OPT1) ../${SOURCE}/bsyn.f
 
eqwidt.o : ../${SOURCE}/eqwidt.f ../${SOURCE}/spectrum.inc ../${SOURCE}/tsuji.par
	$(F) -c $(OPT1) ../${SOURCE}/eqwidt.f
 
eqwidtb.o : ../${SOURCE}/eqwidtb.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/eqwidtb.f
 
eqwidtbplatt.o : ../${SOURCE}/eqwidtbplatt.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/eqwidtbplatt.f
 
bsynb.o : ../${SOURCE}/bsynb.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/bsynb.f
 
bsynbplatt.o : ../${SOURCE}/bsynbplatt.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/bsynbplatt.f
 
clock-gf.o : ../${SOURCE}/clock-gf.f
	$(F) -c $(OPT1) ../${SOURCE}/clock-gf.f

depth.o : ../${SOURCE}/depth.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/depth.f
 
detabs.o : ../${SOURCE}/detabs.f ../${SOURCE}/spectrum.inc ../${SOURCE}/tsuji.par
	$(F) -c $(OPT1) ../${SOURCE}/detabs.f
 
ettmin.o : ../${SOURCE}/ettmin.f
	$(F) -c $(OPT1) ../${SOURCE}/ettmin.f
 
gausi.o : ../${SOURCE}/gausi.f
	$(F) -c $(OPT1) ../${SOURCE}/gausi.f
 
getinfospecies.o : ../${SOURCE}/getinfospecies.f
	$(F) -c $(OPT1) ../${SOURCE}/getinfospecies.f
 
getlele.o : ../${SOURCE}/getlele.f
	$(F) -c $(OPT1) ../${SOURCE}/getlele.f
 
hlinop.o : ../${SOURCE}/hlinop.f 
	$(F) -c $(OPT1) ../${SOURCE}/hlinop.f
 
hydropac.o : ../${SOURCE}/hydropac.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/hydropac.f
 
inabs.o : ../${SOURCE}/inabs.f
	$(F) -c $(OPT1) ../${SOURCE}/inabs.f
 
injon.o : ../${SOURCE}/injon.f
	$(F) -c $(OPT1) ../${SOURCE}/injon.f
 
inp3.o : ../${SOURCE}/inp3.f
	$(F) -c $(OPT1) ../${SOURCE}/inp3.f
 
input.o : ../${SOURCE}/input.f
	$(F) -c $(OPT1) ../${SOURCE}/input.f
 
Iplus_calc.o : ../${SOURCE}/Iplus_calc.f
	$(F) -c $(OPT1) ../${SOURCE}/Iplus_calc.f
 
jon.o : ../${SOURCE}/jon.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/jon.f
 
lenstr.o : ../${SOURCE}/lenstr.f
	$(F) -c $(OPT1) ../${SOURCE}/lenstr.f
 
makeabund.o : ../${SOURCE}/makeabund.f
	$(F) -c $(OPT1) ../${SOURCE}/makeabund.f
 
getisotopmass.o : ../${SOURCE}/getisotopmass.f
	$(F) -c $(OPT1) ../${SOURCE}/getisotopmass.f
 
molmat.o : ../${SOURCE}/molmat.f
	$(F) -c $(OPT1) ../${SOURCE}/molmat.f
 
qas.o : ../${SOURCE}/qas.f
	$(F) -c $(OPT1) ../${SOURCE}/qas.f
 
qtrav.o : ../${SOURCE}/qtrav.f
	$(F) -c $(OPT1) ../${SOURCE}/qtrav.f
 
readmo.o : ../${SOURCE}/readmo.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/readmo.f
 
readfilt.o : ../${SOURCE}/readfilt.f
	$(F) -c $(OPT1) ../${SOURCE}/readfilt.f
 
rossop.o : ../${SOURCE}/rossop.f  ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/rossop.f
 
tabs.o : ../${SOURCE}/tabs.f  ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/tabs.f
 
lint.o : ../${SOURCE}/lint.f
	$(F) -c $(OPT1) ../${SOURCE}/lint.f
 
tint.o : ../${SOURCE}/tint.f
	$(F) -c $(OPT1) ../${SOURCE}/tint.f
 
traneq.o : ../${SOURCE}/traneq.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/traneq.f
 
traneqplatt.o : ../${SOURCE}/traneqplatt.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/traneqplatt.f
 
tranfr.o : ../${SOURCE}/tranfr.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/tranfr.f
 
transc.o : ../${SOURCE}/transc.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/transc.f
 
tranfrplatt.o : ../${SOURCE}/tranfrplatt.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/tranfrplatt.f
 
transcplatt.o : ../${SOURCE}/transcplatt.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/transcplatt.f
 
trrays.o : ../${SOURCE}/trrays.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/trrays.f
 
trqua2.o : ../${SOURCE}/trqua2.f
	$(F) -c $(OPT1) ../${SOURCE}/trqua2.f
 
tranw.o : ../${SOURCE}/tranw.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/tranw.f
 
tryck.o : ../${SOURCE}/tryck.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/tryck.f
 
newvoigt.o : ../${SOURCE}/newvoigt.f
	$(F) -c $(OPT1) ../${SOURCE}/newvoigt.f
 
xeps.o : ../${SOURCE}/xeps.f
	$(F) -c $(OPT1) ../${SOURCE}/xeps.f
 
zerof.o : ../${SOURCE}/zerof.f
	$(F) -c $(OPT1) ../${SOURCE}/zerof.f

heopac.o : ../${SOURCE}/heopac.f
	$(F) -c $(OPT1) ../${SOURCE}/heopac.f
 
h2opac.o : ../${SOURCE}/h2opac.f
	$(F) -c $(OPT1) ../${SOURCE}/h2opac.f
 
CIAh2h2.o : ../${SOURCE}/CIAh2h2.f
	$(F) $(OPT1) -c ../${SOURCE}/CIAh2h2.f

CIAh2h.o : ../${SOURCE}/CIAh2h.f
	$(F) $(OPT1) -c ../${SOURCE}/CIAh2h.f

CIAh2he.o : ../${SOURCE}/CIAh2he.f
	$(F) $(OPT1) -c ../${SOURCE}/CIAh2he.f

CIAhhe.o : ../${SOURCE}/CIAhhe.f
	$(F) $(OPT1) -c ../${SOURCE}/CIAhhe.f

heborysowopac.o : ../${SOURCE}/heborysowopac.f
	$(F) -c $(OPT1) ../${SOURCE}/heborysowopac.f
 
h2borysowopac.o : ../${SOURCE}/h2borysowopac.f
	$(F) -c $(OPT1) ../${SOURCE}/h2borysowopac.f
 
eqmol_pe.o : ../${SOURCE}/eqmol_pe.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/eqmol_pe.f
 
die_pe.o : ../${SOURCE}/die_pe.f
	$(F) -c $(OPT1) ../${SOURCE}/die_pe.f
 
eqmol_pe_lu.o : ../${SOURCE}/eqmol_pe_lu.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/eqmol_pe_lu.f
 
die_pe_lu.o : ../${SOURCE}/die_pe_lu.f
	$(F) -c $(OPT1) ../${SOURCE}/die_pe_lu.f

funcv.o : ../${SOURCE}/funcv.f
	$(F) -c $(OPT1) ../${SOURCE}/funcv.f

ludcmp.o : ../${SOURCE}/ludcmp.f
	$(F) -c $(OPT1) ../${SOURCE}/ludcmp.f

lubksb.o : ../${SOURCE}/lubksb.f
	$(F) -c $(OPT1) ../${SOURCE}/lubksb.f
 
takemolec.o : ../${SOURCE}/takemolec.f ../${SOURCE}/tsuji.par ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/takemolec.f
 
moleq.o : ../${SOURCE}/moleq.f 
	$(F) -c $(OPT1) ../${SOURCE}/moleq.f
 
molfys.o : ../${SOURCE}/molfys.f
	$(F) -c $(OPT1) ../${SOURCE}/molfys.f
 
partf.o : ../${SOURCE}/partf.f
	$(F) -c $(OPT1) ../${SOURCE}/partf.f

molecpartf.o : ../${SOURCE}/molecpartf.f
	$(F) -c $(OPT1) ../${SOURCE}/molecpartf.f

partffordepth.o : ../${SOURCE}/partffordepth.f
	$(F) -c $(OPT1) ../${SOURCE}/partffordepth.f

pemake.o : ../${SOURCE}/pemake.f
	$(F) -c $(OPT1) ../${SOURCE}/pemake.f
 
pemakero.o : ../${SOURCE}/pemakero.f
	$(F) -c $(OPT1) ../${SOURCE}/pemakero.f
 
h2opartf.o : ../${SOURCE}/h2opartf.f
	$(F) -c $(OPT1) ../${SOURCE}/h2opartf.f

h2opartf2001.o : ../${SOURCE}/h2opartf2001.f
	$(F) -c $(OPT1) ../${SOURCE}/h2opartf2001.f

partfamesh2o.o : ../${SOURCE}/partfamesh2o.f
	$(F) -c $(OPT1) ../${SOURCE}/partfamesh2o.f

partfBarberH2O.o : ../${SOURCE}/partfBarberH2O.f
	$(F) -c $(OPT1) ../${SOURCE}/partfBarberH2O.f

weightlimb.o : ../${SOURCE}/weightlimb.f ../${SOURCE}/spectrum.inc
	$(F) -c $(OPT1) ../${SOURCE}/weightlimb.f

rinteg.o : ../${SOURCE}/rinteg.f
	$(F) -c $(OPT1) ../${SOURCE}/rinteg.f


