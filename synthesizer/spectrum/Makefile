CC = gcc
#CFLAGS = -O3
#Replace the CFLAGS line above with the CFLAGS line below to compile with 
#gcc on a 64-bit machine
CFLAGS = -O3 -m64 -mtune=native
LDLIBS = -lm
LOCALBIN = /usr/local/bin

all: spectrum

spectrum: abund2.o al1op3.o autoion3.o balmer8.o brackett.o broad12.o \
c1op_av.o ca1op_av.o capnu.o chop.o coolop5.o density9.o depth.o \
depthfl2.o depthmu.o eqtaukap.o fe1op2.o flux.o fluxflx2.o getisotope.o \
he12.o he13.o he14a.o he15a.o he16a.o he17a.o \
he1op_av.o he313.o he314a.o he315a.o he316a.o he317a.o \
he617a.o helines.o helium6.o heprof4.o hotdensity.o hprofl5.o \
humphreys.o inatom2.o infix.o inisotope.o inline8.o inmodel6.o integ4.o \
intensit.o interva4.o invelgrad.o isorelabun.o linelst12b.o lukeop2.o \
lyman3.o maxcharge.o mg1op_av.o mghop.o ohop.o opacity6.o \
optdepth.o opttrap.o partfn5.o paschen3.o pfinit5.o pfunctio.o \
pfund.o planck.o pop13.o qround.o setreset.o si1op3.o \
spaux.o spectrum.o strong8.o tauflx2.o taukap7.o tauref.o \
tauwave.o trapez.o unified.o veryhotdensity.o voigt.o xi7.o
abund2.o: abund2.c spectrum.h
al1op3.o: al1op3.c spectrum.h
autoion3.o: autoion3.c spectrum.h
balmer8.o: balmer8.c spectrum.h
brackett.o: brackett.c spectrum.h
broad12.o: broad12.c spectrum.h
c1op_av.o: c1op_av.c spectrum.h
ca1op_av.o: ca1op_av.c spectrum.h
capnu.o: capnu.c spectrum.h
chop.o: chop.c spectrum.h
coolop5.o: coolop5.c spectrum.h
density9.o: density9.c spectrum.h
depth.o: depth.c spectrum.h
depthfl2.o: depthfl2.c spectrum.h 
depthmu.o: depthmu.c spectrum.h
eqtaukap.o: eqtaukap.c spectrum.h
fe1op2.o: fe1op2.c spectrum.h
flux.o: flux.c spectrum.h
fluxflx2.o: fluxflx2.c spectrum.h
getisotope.o: getisotope.c spectrum.h
he12.o: he12.c spectrum.h
he13.o: he13.c spectrum.h
he14a.o: he14a.c spectrum.h
he15a.o: he15a.c spectrum.h
he16a.o: he16a.c spectrum.h
he17a.o: he17a.c spectrum.h
he1op_av.o: he1op_av.c spectrum.h 
he313.o: he313.c spectrum.h
he314a.o: he314a.c spectrum.h
he315a.o: he315a.c spectrum.h
he316a.o: he316a.c spectrum.h
he317a.o: he317a.c spectrum.h
he617a.o: he617a.c spectrum.h
helines.o: helines.c spectrum.h
helium6.o: helium6.c spectrum.h
heprof4.o: heprof4.c spectrum.h
hotdensity.o: hotdensity.c spectrum.h
hprofl5.o: hprofl5.c spectrum.h
humphreys.o: humphreys.c spectrum.h
inatom2.o: inatom2.c spectrum.h
infix.o: infix.c spectrum.h
inisotope.o: inisotope.c spectrum.h
inline8.o: inline8.c spectrum.h
inmodel6.o: inmodel6.c spectrum.h
integ4.o: integ4.c spectrum.h
intensit.o: intensit.c spectrum.h
interva4.o: interva4.c spectrum.h
invelgrad.o: invelgrad.c spectrum.h
isorelabund.o: isorelabund.c spectrum.h
linelst12b.o: linelst12b.c spectrum.h
lukeop2.o: lukeop2.c spectrum.h
lyman3.o: lyman3.c spectrum.h
maxcharge.o: maxcharge.c spectrum.h
mg1op_av.o: mg1op_av.c spectrum.h
mghop.o: mghop.c spectrum.h
ohop.o: ohop.c spectrum.h
opacity6.o: opacity6.c spectrum.h
optdepth.o: optdepth.c spectrum.h
opttrap.o: opttrap.c spectrum.h
partfn5.o: partfn5.c spectrum.h
paschen3.o: paschen3.c spectrum.h
pfinit5.o: pfinit5.c spectrum.h
pfunctio.o: pfunctio.c spectrum.h
pfund.o: pfund.c spectrum.h
planck.o: planck.c
pop13.o: pop13.c spectrum.h
qround.o: qround.c spectrum.h
setreset.o: setreset.c spectrum.h
si1op3.o: si1op3.c spectrum.h
spaux.o: spaux.c
spectrum.o: spectrum.c spectrum.h
strong8.o: strong8.c spectrum.h
tauflx2.o: tauflx2.c spectrum.h
taukap7.o: taukap7.c spectrum.h
tauref.o: tauref.c spectrum.h
tauwave.o: tauwave.c spectrum.h
trapez.o: trapez.c spectrum.h
unified.o: unified.c spectrum.h
veryhotdensity.o: veryhotdensity.c spectrum.h
voigt.o: voigt.c spectrum.h
xi7.o: xi7.c spectrum.h 

clean:
	rm -f *.o
	rm -f spectrum
