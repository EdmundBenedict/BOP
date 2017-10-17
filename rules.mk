
modpath := mods

# Main executable and list of sources.
exec := bin/bop bin/bop_dos bin/bop_elas bin/bop_kdos # the Makefile will sort out what are the "main" srcs from the .x file name -> they have to have the same name
libbop.a: $(addsuffix .o,                                                                                                                    \
   src/topologia.f90 src/mod_precision.f90                                                                                                   \
   src/mod_const.f90 src/mod_all_scalar.f90 src/mod_io.F90 src/mod_misc.f90                                                                  \
   src/mod_srt.f90 src/mod_conf.f90 src/ab_io.f90 src/mod_ham.f90 src/mod_chi.f90                                                            \
   src/mod_gsp.f90  src/sa_link.f90  src/g0n.f90 src/mod_funptr.f90 src/avennb.f90 src/bldlist.f90 src/bndinf.f90                            \
   src/bseatm.f90 src/constr.f90 src/dmatml.f90 src/dradf.f90 src/dscale.f90 src/entden.f90 src/findj.f90                                    \
   src/fndrec.f90 src/getmass.f90 src/grdmat.f90 src/inichm.f90 src/initvel.f90 src/insprs.f90 src/instem.f90                                \
   src/inv3x3.f90 src/finv3x3.f90 src/lineqsolv.f90 src/maker.f90 src/matel.f90 src/mdnve.f90 src/mdnvt.f90 src/msd.f90                      \
   src/mstmin.f90 src/mul3x3.f90 src/onebdy.f90 src/onsite.f90 src/orbrot.f90 src/orderp.f90 src/pbc.f90 src/scale.f90 src/radf.f90          \
   src/ran1.f90 src/rdtokn.f90 src/rescale.f90 src/rotbo.f90 src/rpoly.f90 src/states.f90 src/theta.f90 src/zap.f90 src/zcore.f90            \
   src/ai.f90 src/bldclus.f90 src/bldh.f90 src/screenf.f90 src/scrcut.f90 src/dabdl.f90 src/dchinl.f90                                       \
   src/dchisr.f90 src/ddelta.f90 src/delta.f90 src/evalfn.f90 src/simp.f90 src/prodcoeff.f90 src/dot.f90 src/getb.f90                        \
   src/srtint.f90 src/farcdat.f90 src/fdpoly.f90 src/locabinf.f90 src/locsym.f90 src/dndmfn.f90 src/mdnpt.f90 src/pascal.f90                 \
   src/parfrac.f90 src/precab.f90 src/gethu.f90 src/recab.f90 src/trncav.f90 src/trncno.f90 src/utran.f90 src/getroots.f90                   \
   src/poly.f90 src/numelnull.f90 src/numelsrt.f90 src/getchisrt.f90 src/epromavg.f90 src/eprommix.f90 src/epromnoavg.f90                    \
   src/febond.f90 src/getptm.f90 src/evlptm.f90 src/dscreenf.f90 src/utranv.f90 src/getebsnull.f90 src/getchinull.f90                        \
   src/getnch.f90 src/move.f90 src/getefsrt.f90 src/ebsos.f90 src/entropy.f90 src/fintrm.f90 src/getefnull.f90 src/volrelax.f90              \
   src/diffuse.f90 src/scg.F90 src/makcgn.F90 src/mkstrxidx.F90 src/mkstrxd.F90 src/soldhj.F90 src/ropyln.F90 src/ropcsm.F90                 \
   src/ropyln1.F90 src/besslr.F90 src/tbshfl.F90                                                                                             \
                                                                                                                                             \
   src/Library/DiagComp.f90 src/Library/DiagSym.f90 src/Library/FindRoots.f90 src/Library/FindZero.f90                                       \
   src/Library/cdiv.f90 src/Library/dcabs1.f90 src/Library/dswap.f90 src/Library/fft2d.f90 src/Library/four1.f90                             \
   src/Library/htribk.f90 src/Library/htridi.f90 src/Library/integ.f90 src/Library/iyamax.f90 src/Library/pythag.f90                         \
   src/Library/tql1.f90 src/Library/tql2.f90 src/Library/tred2.f90 src/Library/yaxpy.f90 src/Library/ygefa.f90                               \
   src/Library/yscal.f90                                                                                                                     \
                                                                                                                                             \
   src/repulsive_data_mods.f90 src/bop.f90 src/addg3.f90 src/addgr.f90 src/addsw.f90                                                         \
   src/bldlcell.f90 src/dump.f90 src/erasab.f90 src/gammas.f90 src/getdos.f90 src/geteatom.f90                                               \
   src/getelast.f90 src/get_quick_elast.f90 src/getetot.F90 src/getldos.f90                                                                  \
   src/getocc.f90 src/getrho.f90 src/getvib.f90 src/moldyn.f90 src/outfil.f90 src/panic.f90                                                  \
   src/rdab.f90 src/relax.f90 src/repeng.f90 src/repf.f90 src/elcon_bcc.f90 src/elcon_l10.f90                                                \
   src/tag.f90 src/gamma.f90 src/getensurf.f90 src/getenvol.f90 src/get_en.f90                                                               \
   src/gb.f90 src/strain.f90 src/outblock.f90 src/makeblocks.f90 src/maxf_gs.f90                                                             \
   src/maxf_gb.f90 src/relax_sb.f90 src/move_sb.f90 src/sort_out.f90                                                                         \
   src/rose_curve.f90 src/findminhcp.f90 src/findminfcc.f90 src/elcon_hcp.f90 src/elcon_hcp_short.f90                                        \
   src/elcon_fcc.f90 src/get_en_2.f90 src/relax_ds.f90 src/maxf_ds.f90 src/move_ds.f90 src/stress.f90                                        \
   src/logplot.f90 src/report.f90 src/safemin.f90 src/setpress.f90 src/resetup.f90  src/sumrule.f90                                          \
   src/usrexit.f90 src/wrtab.f90 src/wrthis.f90 src/latgfbc.f90 src/diffdisp.f90 src/lgfsafemin.f90                                          \
   src/random.f90 src/mpython.f90 src/gfbc_outp.f90 src/gfbcaf.f90 src/envscr.f90                                                            \
                                                                                                                                             \
   src/kspace/mod_kspace.F90 src/kspace/bsf.F90 src/kspace/kbldh.f90 src/kspace/shiftons.f90 src/kspace/spglib_f08.f90                       \
   src/kspace/atq.f90 src/kspace/keprom.f90 src/kspace/kdiag.f90 src/kspace/iktran.f90 src/kspace/intersite.F90                              \
                                                                                                                                             \
   src/mod_pft.f90 src/print_forces.f90                                                                                                      \
   src/zbrent.f90 src/zriddr.f90 src/eff.f90 src/spnbop.F90 src/femag.F90 src/mod_atom_ar.f90 src/mod_tail.f90                               \
   src/strs_tetr.f90 src/strs_cube.f90 src/bldhdo.f90 src/recurse.f90 src/classic.f90 src/dosplot.f90 src/eval_rfun.f90                      \
   src/emb/emb_sc.F90 src/writecfg.f90 src/writexyz.f90 src/print_bond_scalings.f90                                                          \
   src/writecell.f90 src/writeddp.f90                                                                                                        \
                                                                                                                                             \
                                                                                                                                             \
   src/exit_orderly.f90 src/kspace/kentropy.f90 src/kspace/krhodiag.f90 src/kspace/spnkspc.f90                                               \
   src/kspace/fndocc.f90 src/getebsfbs.f90 src/numel.f90  src/get_dq2chia.f90 src/kspace/kemag.f90                                           \
   src/eval_bsens.f90 src/forcebop.f90 src/kspace/eval_kens.f90 src/kspace/kebsfbs.f90 src/forcecheck.f90                                    \
   src/bldnebt.f90 src/forcedetails.f90 src/neb.f90                                                                                          \
                                                                                                                                             \
   src/mod_par.f90                                                                                                            \
                                                                                                                                             \
   src/kspace/kdosplot.F90                                                                                                                   \
) # only the common srcs

# src/nullmpi.f90

bin/bop: src/main_nmp.f90.o libbop.a


bin/bop_dos bin/bop_elas bin/bop_kdos: bin/bop
	ln -sf $(notdir $<) $@


# List of .mod files produced together with a given .f90.o file.
ab_io.mod : src/ab_io.f90.o
distrib.mod : src/distrib.f90.o
mod_all_scalar.mod : src/mod_all_scalar.f90.o
mod_chi.mod : src/mod_chi.f90.o
mod_clock.mod mod_precision.mod: src/mod_precision.f90.o
mod_conf.mod : src/mod_conf.f90.o
mod_const.mod : src/mod_const.f90.o
mod_fit_model.mod : src/mod_fit_model.f90.o
mod_g0n.mod : src/g0n.f90
mod_gsp.mod : src/mod_gsp.f90.o
mod_ham.mod : src/mod_ham.f90.o
mod_io.mod : src/mod_io.F90.o
mod_kspace.mod : src/kspace/mod_kspace.F90.o
mod_misc.mod : src/mod_misc.f90.o
mod_srt.mod : src/mod_srt.f90.o
pair_coefficients.mod : src/repulsive_data_mods.f90.o
sa_link.mod : src/sa_link.f90.o
sa.mod : src/mod_sa.f90.o
ssa.mod : src/mod_ssa.f90.o
topologia.mod : src/topologia.f90.o
mod_tail.mod : src/mod_tail.f90.o
mod_pft.mod : src/mod_pft.f90.o
#mpi.mod : src/nullmpi.f90.o
mpi.mod :
mod_atom_ar.mod : src/mod_atom_ar.f90.o
mod_funptr.mod : src/mod_funptr.f90.o
spglib_f08.mod : src/kspace/spglib_f08.f90.o
mod_neb.mod : src/neb.f90.o
tbbop_emb.mod : src/emb/emb_sc.F90.o
mod_forcedetails.mod : src/forcedetails.f90.o
mod_par.mod : src/mod_par.f90.o


# List of .f90.o files depending on the given .mod file.
# Multiple targets and prerequisites allowed.
src/main_nmp.f90.o : mod_const.mod mod_conf.mod mod_io.mod mod_clock.mod topologia.mod mod_misc.mod mpi.mod
src/dump.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_conf.mod mod_io.mod
src/ran1.f90.o : mod_precision.mod
src/kspace/eval_kens.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_kspace.mod topologia.mod mpi.mod
src/Library/htridi.f90.o : mod_precision.mod
src/locsym.f90.o : mod_precision.mod
src/writexyz.f90.o : mod_const.mod mod_all_scalar.mod mod_io.mod
src/kspace/bsf.F90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_kspace.mod mod_atom_ar.mod topologia.mod mpi.mod
src/writecfg.f90.o : mod_const.mod mod_all_scalar.mod mod_atom_ar.mod
src/inichm.f90.o : mod_const.mod
src/ai.f90.o : mod_precision.mod
src/febond.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_chi.mod mod_ham.mod
src/screenf.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_atom_ar.mod
src/dot.f90.o : mod_precision.mod
src/getebsfbs.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod topologia.mod mod_ham.mod ab_io.mod mod_chi.mod mod_clock.mod mod_conf.mod mod_atom_ar.mod mpi.mod mod_par.mod
src/erasab.f90.o : mod_precision.mod ab_io.mod
src/onsite.f90.o : mod_precision.mod
src/gammas.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_srt.mod
src/outfil.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_srt.mod ab_io.mod mod_g0n.mod mod_ham.mod mod_chi.mod mod_kspace.mod mod_funptr.mod mod_atom_ar.mod
src/rdab.f90.o : mod_precision.mod ab_io.mod
src/states.f90.o : mod_precision.mod mod_const.mod
src/avennb.f90.o : mod_precision.mod
src/outblock.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_atom_ar.mod
src/gb.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/orbrot.f90.o : mod_precision.mod
src/Library/yaxpy.f90.o : mod_precision.mod
src/dradf.f90.o : mod_precision.mod
src/zcore.f90.o : mod_precision.mod
src/Library/fft2d.f90.o : mod_precision.mod
src/wrtab.f90.o : mod_precision.mod ab_io.mod
src/getvib.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_atom_ar.mod
src/Library/dswap.f90.o : mod_precision.mod
src/logplot.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/geteatom.f90.o : mod_precision.mod
src/elcon_fcc.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/scale.f90.o : mod_precision.mod mod_gsp.mod
src/findminfcc.f90.o : mod_precision.mod mod_all_scalar.mod
src/scrcut.f90.o : mod_precision.mod mod_const.mod
src/mod_all_scalar.f90.o : mod_const.mod
src/mod_conf.f90.o : mod_precision.mod mod_io.mod mod_tail.mod mod_pft.mod mod_const.mod
src/pascal.f90.o : mod_precision.mod mod_const.mod mod_srt.mod
src/mdnpt.f90.o : mod_precision.mod
src/findj.f90.o : mod_precision.mod mod_const.mod
src/move_ds.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/writecell.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_conf.mod mod_atom_ar.mod mod_io.mod
src/Library/tred2.f90.o : mod_precision.mod
src/spnbop.F90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_chi.mod mod_clock.mod topologia.mod mod_atom_ar.mod
src/trncav.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_ham.mod
src/getrho.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_funptr.mod
src/exit_orderly.f90.o : mod_io.mod ab_io.mod mod_ham.mod mod_chi.mod topologia.mod mpi.mod
src/ebsos.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_srt.mod ab_io.mod mod_g0n.mod
src/dchinl.f90.o : mod_precision.mod
src/entden.f90.o : mod_precision.mod
src/latgfbc.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/radf.f90.o : mod_precision.mod
src/mod_gsp.f90.o : mod_precision.mod
src/bndinf.f90.o : mod_precision.mod
src/epromnoavg.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_funptr.mod ab_io.mod
src/getb.f90.o : mod_precision.mod
src/repeng.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_atom_ar.mod
src/safemin.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/zap.f90.o : mod_precision.mod
src/strs_cube.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_conf.mod
src/dscreenf.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_atom_ar.mod
src/finv3x3.f90.o : mod_precision.mod
src/kspace/kentropy.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_kspace.mod topologia.mod
src/addg3.f90.o : mod_precision.mod
src/stress.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/wrthis.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/mod_pft.f90.o : mod_precision.mod mod_const.mod mod_tail.mod
src/zbrent.f90.o : mod_precision.mod mod_all_scalar.mod
src/get_en_2.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/Library/DiagComp.f90.o : mod_precision.mod
src/mod_chi.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod
src/print_forces.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod
src/relax.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/maker.f90.o : mod_precision.mod
src/locabinf.f90.o : mod_precision.mod
src/rescale.f90.o : mod_precision.mod
src/addsw.f90.o : mod_precision.mod
src/print_bond_scalings.f90.o : mod_precision.mod mod_const.mod mod_conf.mod mod_io.mod mod_pft.mod
src/mod_tail.f90.o : mod_precision.mod mod_const.mod
src/kspace/keprom.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_kspace.mod
src/kspace/spnkspc.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_kspace.mod topologia.mod mod_io.mod mod_atom_ar.mod mpi.mod
src/getroots.f90.o : mod_precision.mod
src/dchisr.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_srt.mod ab_io.mod mod_chi.mod mod_g0n.mod
src/Library/cdiv.f90.o : mod_precision.mod
src/rotbo.f90.o : mod_precision.mod
src/forcecheck.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_conf.mod mod_io.mod mod_atom_ar.mod topologia.mod mpi.mod
src/kspace/intersite.F90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_kspace.mod mod_atom_ar.mod topologia.mod
src/utranv.f90.o : mod_precision.mod
src/relax_sb.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/gamma.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/getnch.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod
src/recab.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_ham.mod
src/numelnull.f90.o : mod_precision.mod
src/getchinull.f90.o : mod_precision.mod
src/lineqsolv.f90.o : mod_precision.mod
src/getdos.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_funptr.mod mod_g0n.mod
src/Library/yscal.f90.o : mod_precision.mod
src/initvel.f90.o : mod_precision.mod
src/mod_ham.f90.o : mod_precision.mod mod_const.mod
src/orderp.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/dmatml.f90.o : mod_precision.mod
src/kspace/iktran.f90.o : mod_precision.mod mod_const.mod
src/grdmat.f90.o : mod_precision.mod mod_const.mod mod_conf.mod
src/dscale.f90.o : mod_precision.mod mod_gsp.mod
src/msd.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/gfbc_outp.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/parfrac.f90.o : mod_precision.mod
src/kspace/krhodiag.f90.o : mod_precision.mod mod_const.mod
src/elcon_hcp_short.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_conf.mod mod_io.mod
src/strs_tetr.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_conf.mod
src/dosplot.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_io.mod mod_g0n.mod mod_conf.mod mod_funptr.mod
src/getchisrt.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_srt.mod ab_io.mod mod_chi.mod mod_g0n.mod
src/precab.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_ham.mod mod_clock.mod
src/bldlcell.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_atom_ar.mod
src/numel.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod topologia.mod ab_io.mod mod_funptr.mod mod_atom_ar.mod mpi.mod
src/theta.f90.o : mod_precision.mod
src/get_en.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/maxf_ds.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod topologia.mod
src/makeblocks.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/numelsrt.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_srt.mod ab_io.mod mod_g0n.mod
src/elcon_hcp.f90.o : mod_precision.mod mod_const.mod
src/pbc.f90.o : mod_precision.mod
src/getocc.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/elcon_bcc.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/kspace/mod_kspace.F90.o : mod_io.mod mod_precision.mod mod_const.mod mod_all_scalar.mod mod_conf.mod topologia.mod mod_atom_ar.mod spglib_f08.mod mod_clock.mod
src/getenvol.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/simp.f90.o : mod_precision.mod
src/mul3x3.f90.o : mod_precision.mod
src/get_quick_elast.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/gethu.f90.o : mod_precision.mod mod_const.mod mod_ham.mod
src/emb/emb_sc.F90.o : mod_precision.mod mod_const.mod mod_conf.mod mod_clock.mod topologia.mod mpi.mod mod_all_scalar.mod mod_atom_ar.mod ab_io.mod mod_ham.mod mod_g0n.mod mod_chi.mod
src/bldhdo.f90.o : mod_precision.mod mod_all_scalar.mod mod_ham.mod mod_const.mod mod_atom_ar.mod
src/mdnve.f90.o : mod_precision.mod
src/inv3x3.f90.o : mod_precision.mod mod_all_scalar.mod
src/bldnebt.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_clock.mod topologia.mod mpi.mod
src/relax_ds.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod topologia.mod
src/Library/iyamax.f90.o : mod_precision.mod
src/volrelax.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/Library/htribk.f90.o : mod_precision.mod
src/Library/ygefa.f90.o : mod_precision.mod
src/Library/four1.f90.o : mod_precision.mod
src/eval_rfun.f90.o : mod_precision.mod mod_const.mod mod_conf.mod mod_atom_ar.mod mod_pft.mod
src/delta.f90.o : mod_precision.mod
src/g0n.f90.o : mod_precision.mod
src/kspace/kdiag.f90.o : mod_precision.mod mod_const.mod mod_io.mod
src/Library/integ.f90.o : mod_precision.mod
src/sumrule.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod
src/poly.f90.o : mod_precision.mod
src/constr.f90.o : mod_precision.mod
src/usrexit.f90.o : mod_precision.mod
src/forcedetails.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_conf.mod ab_io.mod mod_ham.mod mod_chi.mod topologia.mod mod_clock.mod mod_atom_ar.mod mod_pft.mod mod_io.mod mpi.mod
src/evalfn.f90.o : mod_precision.mod mod_all_scalar.mod mod_srt.mod ab_io.mod
src/sort_out.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/farcdat.f90.o : mod_precision.mod
src/dndmfn.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_srt.mod ab_io.mod mod_g0n.mod
src/mod_misc.f90.o : mod_precision.mod mod_io.mod
src/getelast.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_conf.mod mod_io.mod
src/kspace/kbldh.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_kspace.mod topologia.mod mod_atom_ar.mod mod_io.mod
src/panic.f90.o : mod_precision.mod
src/getldos.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_funptr.mod mod_g0n.mod
src/moldyn.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_atom_ar.mod
src/Library/pythag.f90.o : mod_precision.mod
src/resetup.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_gsp.mod mod_conf.mod mod_tail.mod mod_pft.mod topologia.mod mod_atom_ar.mod mod_srt.mod pair_coefficients.mod ab_io.mod mod_ham.mod mod_funptr.mod mod_kspace.mod
src/fndrec.f90.o : mod_precision.mod
src/eval_bsens.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_chi.mod ab_io.mod mod_ham.mod mod_funptr.mod topologia.mod mpi.mod
src/setpress.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/get_dq2chia.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_chi.mod ab_io.mod mod_funptr.mod topologia.mod mod_atom_ar.mod
src/Library/dcabs1.f90.o : mod_precision.mod
src/classic.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_clock.mod mod_gsp.mod mod_conf.mod mod_atom_ar.mod mod_pft.mod
src/getefsrt.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/move_sb.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/instem.f90.o : mod_precision.mod
src/gfbcaf.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod topologia.mod
src/repf.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_atom_ar.mod
src/mod_atom_ar.f90.o : mod_precision.mod mod_const.mod
src/bldclus.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_ham.mod
src/recurse.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_srt.mod mod_ham.mod ab_io.mod topologia.mod mod_clock.mod
src/kspace/kemag.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_atom_ar.mod mpi.mod
src/addgr.f90.o : mod_precision.mod
src/mdnvt.f90.o : mod_precision.mod
src/prodcoeff.f90.o : mod_precision.mod
src/evlptm.f90.o : mod_precision.mod
src/insprs.f90.o : mod_precision.mod
src/elcon_l10.f90.o : mod_precision.mod mod_const.mod
src/kspace/shiftons.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_kspace.mod topologia.mod mod_io.mod
src/getebsnull.f90.o : mod_precision.mod
src/onebdy.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/sa_link.f90.o : mod_precision.mod mod_conf.mod mod_all_scalar.mod
src/random.f90.o : mod_precision.mod
src/femag.F90.o : mod_precision.mod mod_const.mod topologia.mod mod_atom_ar.mod
src/bldlist.f90.o : mod_precision.mod mod_const.mod
src/mod_srt.f90.o : mod_precision.mod mod_const.mod
src/utran.f90.o : mod_precision.mod
src/move.f90.o : mod_precision.mod
src/trncno.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_ham.mod
src/diffdisp.f90.o : mod_precision.mod
src/getetot.F90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_conf.mod tbbop_emb.mod mod_forcedetails.mod
src/maxf_gb.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/getensurf.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/fintrm.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_ham.mod
src/mod_const.f90.o : mod_precision.mod
src/mod_par.f90.o : mod_precision.mod mpi.mod
src/epromavg.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_funptr.mod ab_io.mod
src/mstmin.f90.o : mod_precision.mod
src/repulsive_data_mods.f90.o : mod_precision.mod mod_conf.mod
src/entropy.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_srt.mod ab_io.mod mod_g0n.mod
src/matel.f90.o : mod_precision.mod mod_const.mod mod_conf.mod
src/Library/tql2.f90.o : mod_precision.mod
src/Library/FindZero.f90.o : mod_precision.mod
src/rose_curve.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/eff.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/strain.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod topologia.mod
src/writeddp.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_io.mod
src/kspace/fndocc.f90.o : mod_precision.mod mod_const.mod mod_all_scalar.mod mod_kspace.mod topologia.mod mpi.mod
src/neb.f90.o : mod_precision.mod mod_const.mod mod_io.mod mod_atom_ar.mod mod_all_scalar.mod topologia.mod mpi.mod mod_conf.mod
src/srtint.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_srt.mod ab_io.mod
src/kspace/kebsfbs.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod topologia.mod mod_clock.mod mod_conf.mod mod_kspace.mod mod_atom_ar.mod
src/lgfsafemin.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/Library/FindRoots.f90.o : mod_precision.mod
src/bseatm.f90.o : mod_precision.mod
src/mod_funptr.f90.o : mod_precision.mod
src/rpoly.f90.o : mod_precision.mod
src/dabdl.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_ham.mod
src/diffuse.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/Library/tql1.f90.o : mod_precision.mod
src/tag.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_srt.mod mod_atom_ar.mod
src/report.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_chi.mod mod_funptr.mod
src/findminhcp.f90.o : mod_precision.mod mod_all_scalar.mod
src/bldh.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_ham.mod
src/Library/DiagSym.f90.o : mod_precision.mod
src/forcebop.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_ham.mod mod_chi.mod mod_funptr.mod topologia.mod mpi.mod
src/ab_io.f90.o : mod_precision.mod mod_const.mod mod_srt.mod mod_all_scalar.mod
src/mod_io.F90.o : mod_precision.mod
src/maxf_gs.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/ddelta.f90.o : mod_precision.mod
src/zriddr.f90.o : mod_precision.mod
src/fdpoly.f90.o : mod_precision.mod mod_const.mod mod_srt.mod
src/eprommix.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod ab_io.mod mod_funptr.mod mod_chi.mod
src/getmass.f90.o : mod_precision.mod
src/bop.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_conf.mod mod_ham.mod ab_io.mod mod_clock.mod mod_chi.mod mod_kspace.mod mod_atom_ar.mod topologia.mod mod_neb.mod
src/envscr.f90.o : mod_precision.mod mod_const.mod
src/mpython.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/getptm.f90.o : mod_precision.mod
src/kspace/atq.f90.o : mod_precision.mod mod_const.mod mod_kspace.mod
src/getefnull.f90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod
src/kspace/kdosplot.F90.o : mod_precision.mod mod_all_scalar.mod mod_const.mod mod_io.mod topologia.mod mod_kspace.mod


all: $(exec)


.SECONDEXPANSION:
%.f90.o : %.f90 | $(modpath)/ $$(@D)/
	$(fc) $(fflags) -c $< -o $@

.SECONDEXPANSION:
%.F90.o : %.F90 | $(modpath)/ $$(@D)/
	$(fc) $(fflags) -c $< -o $@

%.c.o : %.c
	$(cc) $(cflags) -c $< -o $@

%.cxx.o : %.cxx
	$(cxx) $(cxxflags) -c $< -o $@



%.a: $^
	ar rc $@ $^

%.x: $^
	$(fc) $(fflags) $^ $(ldflags) -o $@

.SECONDEXPANSION:
bin/%: $$^ | bin/
	$(fc) $(fflags) $^ $(ldflags) -o $@

%/:
	mkdir -p $@

.SUFFIXES:
