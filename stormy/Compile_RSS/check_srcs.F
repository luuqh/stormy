#include "cppdefs.h"

      subroutine check_srcs

!!!!!! WARNING: THIS IS A MACHINE GENERATED CODE, DO NOT EDIT! !!!!!!
!!!!!! This file needs to be updated only if the new files     !!!!!!
!!!!!! were introduced into or deleted from the list of source !!!!!!
!!!!!! codes SRCS in the Makefile.                             !!!!!!

      implicit none
      integer i
#include "strings.h"
      do i=1,max_opt_size
        srcs(i:i)=' '
      enddo
      srcs(   1:   6)='main.F'
      srcs(   7:  13)=' step.F'
      srcs(  14:  24)=' read_inp.F'
      srcs(  25:  33)=' timers.F'
      srcs(  34:  48)=' init_scalars.F'
      srcs(  49:  62)=' init_arrays.F'
      srcs(  63:  76)=' set_weights.F'
      srcs(  77:  89)=' set_scoord.F'
      srcs(  90: 100)=' ana_grid.F'
      srcs( 101: 114)=' setup_grid1.F'
      srcs( 115: 128)=' setup_grid2.F'
      srcs( 129: 142)=' set_nudgcof.F'
      srcs( 143: 156)=' ana_initial.F'
      srcs( 157: 169)=' analytical.F'
      srcs( 170: 182)=' setup_kwds.F'
      srcs( 183: 195)=' check_kwds.F'
      srcs( 196: 208)=' check_srcs.F'
      srcs( 209: 226)=' check_switches1.F'
      srcs( 227: 244)=' check_switches2.F'
      srcs( 245: 253)=' step2d.F'
      srcs( 254: 261)=' u2dbc.F'
      srcs( 262: 269)=' v2dbc.F'
      srcs( 270: 278)=' zetabc.F'
      srcs( 279: 292)=' obc_volcons.F'
      srcs( 293: 305)=' pre_step3d.F'
      srcs( 306: 316)=' step3d_t.F'
      srcs( 317: 329)=' step3d_uv1.F'
      srcs( 330: 342)=' step3d_uv2.F'
      srcs( 343: 351)=' prsgrd.F'
      srcs( 352: 359)=' rhs3d.F'
      srcs( 360: 371)=' set_depth.F'
      srcs( 372: 379)=' omega.F'
      srcs( 380: 389)=' uv3dmix.F'
      srcs( 390: 403)=' uv3dmix_spg.F'
      srcs( 404: 412)=' t3dmix.F'
      srcs( 413: 425)=' t3dmix_spg.F'
      srcs( 426: 437)=' hmix_coef.F'
      srcs( 438: 445)=' u3dbc.F'
      srcs( 446: 453)=' v3dbc.F'
      srcs( 454: 461)=' t3dbc.F'
      srcs( 462: 471)=' rho_eos.F'
      srcs( 472: 482)=' ab_ratio.F'
      srcs( 483: 493)=' alfabeta.F'
      srcs( 494: 504)=' ana_vmix.F'
      srcs( 505: 514)=' bvf_mix.F'
      srcs( 515: 525)=' lmd_vmix.F'
      srcs( 526: 536)=' lmd_skpp.F'
      srcs( 537: 547)=' lmd_bkpp.F'
      srcs( 548: 560)=' lmd_swfrac.F'
      srcs( 561: 573)=' lmd_wscale.F'
      srcs( 574: 580)=' diag.F'
      srcs( 581: 589)=' wvlcty.F'
      srcs( 590: 601)=' checkdims.F'
      srcs( 602: 618)=' grid_stiffness.F'
      srcs( 619: 629)=' bio_diag.F'
      srcs( 630: 640)=' get_date.F'
      srcs( 641: 649)=' lenstr.F'
      srcs( 650: 660)=' closecdf.F'
      srcs( 661: 676)=' ana_initracer.F'
      srcs( 677: 690)=' insert_node.F'
      srcs( 691: 701)=' nf_fread.F'
      srcs( 702: 712)=' get_grid.F'
      srcs( 713: 726)=' get_initial.F'
      srcs( 727: 737)=' def_grid.F'
      srcs( 738: 747)=' def_his.F'
      srcs( 748: 757)=' def_rst.F'
      srcs( 758: 767)=' set_avg.F'
      srcs( 768: 778)=' wrt_grid.F'
      srcs( 779: 788)=' wrt_his.F'
      srcs( 789: 798)=' wrt_rst.F'
      srcs( 799: 808)=' wrt_avg.F'
      srcs( 809: 817)=' output.F'
      srcs( 818: 835)=' put_global_atts.F'
      srcs( 836: 847)=' def_diags.F'
      srcs( 848: 859)=' wrt_diags.F'
      srcs( 860: 875)=' wrt_diags_avg.F'
      srcs( 876: 891)=' set_diags_avg.F'
      srcs( 892: 904)=' def_diagsM.F'
      srcs( 905: 917)=' wrt_diagsM.F'
      srcs( 918: 934)=' wrt_diagsM_avg.F'
      srcs( 935: 951)=' set_diagsM_avg.F'
      srcs( 952: 967)=' def_bio_diags.F'
      srcs( 968: 983)=' wrt_bio_diags.F'
      srcs( 984:1003)=' wrt_bio_diags_avg.F'
      srcs(1004:1023)=' set_bio_diags_avg.F'
      srcs(1024:1033)=' get_vbc.F'
      srcs(1034:1045)=' set_cycle.F'
      srcs(1046:1057)=' get_wwave.F'
      srcs(1058:1067)=' get_slp.F'
      srcs(1068:1080)=' get_tclima.F'
      srcs(1081:1093)=' get_uclima.F'
      srcs(1094:1103)=' get_ssh.F'
      srcs(1104:1113)=' get_sss.F'
      srcs(1114:1126)=' get_smflux.F'
      srcs(1127:1139)=' get_stflux.F'
      srcs(1140:1152)=' get_srflux.F'
      srcs(1153:1162)=' get_sst.F'
      srcs(1163:1174)=' get_tides.F'
      srcs(1175:1186)=' clm_tides.F'
      srcs(1187:1197)=' get_bulk.F'
      srcs(1198:1209)=' bulk_flux.F'
      srcs(1210:1219)=' get_bry.F'
      srcs(1220:1233)=' get_bry_bio.F'
      srcs(1234:1247)=' nf_read_bry.F'
      srcs(1248:1259)=' MPI_Setup.F'
      srcs(1260:1272)=' MessPass2D.F'
      srcs(1273:1285)=' MessPass3D.F'
      srcs(1286:1296)=' exchange.F'
      srcs(1297:1306)=' biology.F'
      srcs(1307:1315)=' o2sato.F'
      srcs(1316:1326)=' sediment.F'
      srcs(1327:1332)=' bbl.F'
      srcs(1333:1346)=' init_floats.F'
      srcs(1347:1359)=' wrt_floats.F'
      srcs(1360:1373)=' step_floats.F'
      srcs(1374:1386)=' rhs_floats.F'
      srcs(1387:1399)=' interp_rho.F'
      srcs(1400:1412)=' def_floats.F'
      srcs(1413:1433)=' init_arrays_floats.F'
      srcs(1434:1447)=' random_walk.F'
      srcs(1448:1468)=' get_initial_floats.F'
      srcs(1469:1479)=' init_sta.F'
      srcs(1480:1489)=' wrt_sta.F'
      srcs(1490:1500)=' step_sta.F'
      srcs(1501:1513)=' interp_sta.F'
      srcs(1514:1523)=' def_sta.F'
      srcs(1524:1541)=' init_arrays_sta.F'
      srcs(1542:1548)=' zoom.F'
      srcs(1549:1559)=' update2D.F'
      srcs(1560:1578)=' set_nudgcof_fine.F'
      srcs(1579:1590)=' zoombc_2D.F'
      srcs(1591:1602)=' zoombc_3D.F'
      srcs(1603:1615)=' uv3dpremix.F'
      srcs(1616:1627)=' t3dpremix.F'
      srcs(1628:1638)=' update3D.F'
      srcs(1639:1649)=' ropigrid.F'
      srcs(1650:1660)=' ropistep.F'
      srcs(1661:1671)=' pirostep.F'
      srcs(1672:1680)=' trcini.F'
      srcs(1681:1689)=' trclsm.F'
      srcs(1690:1702)=' pisces_ini.F'
      srcs(1703:1711)=' p4zday.F'
      srcs(1712:1721)=' p4znano.F'
      srcs(1722:1730)=' p4zlys.F'
      srcs(1731:1740)=' p4zdiat.F'
      srcs(1741:1749)=' p4zopt.F'
      srcs(1750:1759)=' p4zsink.F'
      srcs(1760:1768)=' p4zflx.F'
      srcs(1769:1777)=' p4zprg.F'
      srcs(1778:1787)=' p4zslow.F'
      srcs(1788:1796)=' p4zint.F'
      srcs(1797:1806)=' p4zprod.F'
      srcs(1807:1817)=' p4zmicro.F'
      srcs(1818:1826)=' p4zlim.F'
      srcs(1827:1835)=' p4zrem.F'
      srcs(1836:1844)=' p4zche.F'
      srcs(1845:1853)=' p4zbio.F'
      srcs(1854:1863)=' p4zmeso.F'
      srcs(1864:1874)=' p4zsink2.F'
      srcs(1875:1883)=' p4zsed.F'
      srcs(1884:1896)=' autotiling.F'
      srcs(1897:1904)=' debug.F'
      srcs(1905:1918)=' Agrif2Model.F'
      return
      end
