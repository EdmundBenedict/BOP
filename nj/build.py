#!/usr/bin/env python

import os, sys, re
from nj import writefile, laddprefix, laddsuffix, nj_t


nj = nj_t(fppstyle='cpp')


ofln = list(filter(lambda o: not o.startswith('-'), sys.argv[1:]))
ofln = 'build.ninja' if len(ofln) == 0 else ofln[0]
nj.bld(['$srcpath/build.py'], 'nj', ['$bldpath/'+ofln], si='flags.mk $srcpath/nj/nj.py'.split())
writefile(ofln, nj.rules)

mpiavail = '-DMPI' in (nj.flags['fc'].split() + nj.flags['fflags'].split())

nj.emods.update(set(['xc_f90_lib_m', 'mpi', 'elpa1', 'elpa2']))
if not mpiavail: nj.emods -= set(['mpi'])

src = laddprefix('src/','''
   topologia.f90 mod_precision.f90
   mod_const.f90 mod_all_scalar.f90 mod_io.F90 mod_misc.f90
   mod_srt.f90 mod_conf.f90 ab_io.f90 mod_ham.f90 mod_chi.f90
   mod_gsp.f90  sa_link.f90  g0n.f90 mod_funptr.f90 avennb.f90 bldlist.f90 bndinf.f90
   bseatm.f90 constr.f90 dmatml.f90 dradf.f90 dscale.f90 entden.f90 findj.f90
   fndrec.f90 getmass.f90 grdmat.f90 inichm.f90 initvel.f90 insprs.f90 instem.f90
   inv3x3.f90 finv3x3.f90 lineqsolv.f90 maker.f90 matel.f90 mdnve.f90 mdnvt.f90 msd.f90
   mstmin.f90 mul3x3.f90 onebdy.f90 onsite.f90 orbrot.f90 orderp.f90 pbc.f90 scale.f90 radf.f90
   ran1.f90 rdtokn.f90 rescale.f90 rotbo.f90 rpoly.f90 states.f90 theta.f90 zap.f90 zcore.f90
   ai.f90 bldclus.f90 bldh.f90 screenf.f90 scrcut.f90 dabdl.f90 dchinl.f90
   dchisr.f90 ddelta.f90 delta.f90 evalfn.f90 simp.f90 prodcoeff.f90 dot.f90 getb.f90
   srtint.f90 farcdat.f90 fdpoly.f90 locabinf.f90 locsym.f90 dndmfn.f90 mdnpt.f90 pascal.f90
   parfrac.f90 precab.f90 gethu.f90 recab.f90 trncav.f90 trncno.f90 utran.f90 getroots.f90
   poly.f90 numelnull.f90 numelsrt.f90 getchisrt.f90 epromavg.f90 eprommix.f90 epromnoavg.f90
   febond.f90 getptm.f90 evlptm.f90 dscreenf.f90 utranv.f90 getebsnull.f90 getchinull.f90
   getnch.f90 move.f90 getefsrt.f90 ebsos.f90 entropy.f90 fintrm.f90 getefnull.f90 volrelax.f90
   diffuse.f90

   Library/DiagComp.f90 Library/DiagSym.f90 Library/FindRoots.f90 Library/FindZero.f90
   Library/cdiv.f90 Library/dcabs1.f90 Library/dswap.f90 Library/fft2d.f90 Library/four1.f90
   Library/htribk.f90 Library/htridi.f90 Library/integ.f90 Library/iyamax.f90 Library/pythag.f90
   Library/tql1.f90 Library/tql2.f90 Library/tred2.f90 Library/yaxpy.f90 Library/ygefa.f90
   Library/yscal.f90

   repulsive_data_mods.f90 bop.f90 addg3.f90 addgr.f90 addsw.f90
   bldlcell.f90 dump.f90 erasab.f90 gammas.f90 getdos.f90 geteatom.f90
   getelast.f90 get_quick_elast.f90 getetot.F90 getldos.f90
   getocc.f90 getrho.f90 getvib.f90 moldyn.f90 outfil.f90 panic.f90
   rdab.f90 relax.f90 repeng.f90 repf.f90 elcon_bcc.f90 elcon_l10.f90
   tag.f90 gamma.f90 getensurf.f90 getenvol.f90 get_en.f90
   gb.f90 strain.f90 outblock.f90 makeblocks.f90 maxf_gs.f90
   maxf_gb.f90 relax_sb.f90 move_sb.f90 sort_out.f90
   rose_curve.f90 findminhcp.f90 findminfcc.f90 elcon_hcp.f90 elcon_hcp_short.f90
   elcon_fcc.f90 get_en_2.f90 relax_ds.f90 maxf_ds.f90 move_ds.f90 stress.f90
   logplot.f90 report.f90 safemin.f90 setpress.f90 resetup.f90  sumrule.f90
   usrexit.f90 wrtab.f90 wrthis.f90 latgfbc.f90 diffdisp.f90 lgfsafemin.f90
   random.f90 mpython.f90 gfbc_outp.f90 gfbcaf.f90 envscr.f90

   kspace/mod_kspace.F90 kspace/bsf.F90 kspace/kbldh.f90 kspace/shiftons.f90 kspace/spglib_f08.f90
   kspace/atq.f90 kspace/keprom.f90 kspace/kdiag.f90 kspace/iktran.f90 kspace/intersite.F90

   mod_pft.f90 print_forces.f90
   zbrent.f90 zriddr.f90 eff.f90 spnbop.F90 femag.F90 mod_atom_ar.f90 mod_tail.f90
   strs_tetr.f90 strs_cube.f90 bldhdo.f90 recurse.f90 classic.f90 dosplot.f90 eval_rfun.f90
   emb/emb_sc.F90 writecfg.f90 writexyz.f90 print_bond_scalings.f90
   writecell.f90 writeddp.f90


   exit_orderly.f90 kspace/kentropy.f90 kspace/krhodiag.f90 kspace/spnkspc.f90
   kspace/fndocc.f90 getebsfbs.f90 numel.f90  get_dq2chia.f90 kspace/kemag.f90
   eval_bsens.F90 forcebop.f90 kspace/eval_kens.f90 kspace/kebsfbs.f90 forcecheck.f90
   bldnebt.f90 forcedetails.f90 neb.f90

   mod_par.F90

   kspace/kdosplot.F90
   
   mkstrxd.F90 scg.F90 makcgn.F90  mkstrxidx.F90 soldhj.F90  ropyln.F90  ropcsm.F90                 \
    ropyln1.F90  besslr.F90  tbshfl.F90 scglp1.F90 qmix.F90 broyj.F90 dgedi.F90
    dsifa.F90 dsidi.F90 dsmpy.F90 dqinv.F90 
''')

if not mpiavail: src.append('src/nullmpi.f90')

nj.mks(src)

if '-time' in sys.argv: nj.printtimes()

nj.mk(laddsuffix('.o', src), o='libbop.a')
nj.mks(['src/main_nmp.f90'])

nj.mk('src/main_nmp.f90.o libbop.a', o='bin/bop')


#bin/bop_dos bin/bop_elas bin/bop_kdos: bin/bop

if (nj.srcpath != nj.bldpath): nj.mkincludes()


for i in laddprefix('$bldpath/bin/', 'bop_dos bop_elas bop_kdos'):
    nj.bld([], 'cmd', [i], so=['$bldpath/bin/bop'], vs=dict(command='ln -sf bop ' + i))


nj.bld(laddprefix('$bldpath/bin/','bop bop_dos bop_elas bop_kdos'),'phony',['all'])

nj.bld([],'cmd',['clean'],vs=dict(command='ninja -t clean'))

nj.direct += '\ndefault all\n'

nj.write(ofln)

