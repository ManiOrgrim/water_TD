# -*- coding: utf-8 -*-
"""
Created on Sun Jul  6 16:59:23 2025

@author: manim
"""

def modify_input32(flux, Tbot, HPperm, HPporo, LPperm, LPporo, dZ):
    texto ="""# Hydrotherm Data-Input Form
# Version 3.2
#      Notes:
#      Syntax of input line is denoted by # ..
#      A suffix letter indicates an exclusive record choice must be made.
#             i.e. A or B or C
#      (O) - Optional data with conditions for requirement
#      [a|b] - Indicates that either A or B must be entered
#      {{nnn}} - Indicates that the default number, nnn, is used if a zero
#              is entered for that variable
#      [T/F] - Indicates a logical variable
#      [I] - Indicates an integer variable
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#   Start of the data file
# TITLE
# ..  TITLE:   title line 1
# ..  VER3.2:  title line 2
TITLE  Nisi_linke system 
VER3.2 hydrothermal system without input
# DIMENSIONS
# ..  nx[I],ny[I],nz[I],tmstr,iyear[I]
16    1   32   0.0   2
# TIME STEP
# ..A  tstepmx[I],tmincr{{1.3}},%PchgTS{{10.}},%HchgTS{{5.}},Wsatn{{0.03}},relax{{1.}},
# ..       deltmin,n_tmcut{{20}}; confined flow
     1000     1.3     10.0     7.0     0.1     1.0     1.0E-10     40
# N-R ITERATION
# ..A  itermx{{10}},resmas{{1.e-12}},reseng{{1.e-2}},nphasechg{{10}},%PchgNR{{0.1}},
# ..        %HchgNR{{0.1}}; confined flow
     40     1.0E-12     1.0E-4     32     0.1     0.1
# LINEAR EQUATION SOLVER
# .. slmeth[I]
1
#    SSOR SOLVER CONTROL
# ..A  ssormx[I]{{10}},tol{{0.05}},ssorw0{{1.0}}; (O) - slmeth = 1
    1   0.01   1.0
# WEIGHTING AND AVERAGING
# ..  ioptupst[T/F],potdif{{0.0002}},theta{{1.0}}
     T     2.0E-5     1.0
# RELATIVE PERMEABILITY
# ..  kodrp[I],wsatn{{0.3}},ssatn{{0.0}}
     0     0.3     0.05
# ROCK PROPERTIES
# ..  heatcpty{{1.e7}},rxden,rxcmprss{{0.}},grav{{981.}},initphi[T/F]
# ..  initphi_filename{{IC_pressporo.xxx}}; (O) - initphi = T
   1.0e7  2.5,    0,   981.0,   F
# CYLINDRICAL COORDINATES
# ..  irad[T/F],wellrad,radiusmax
     T    0.01    16000
# PRINT/PLOT CONTROL
# ..  wide[I],long[I],v-avg[I]
     2     1     0
# PRINT 1
# ..  pres_pr_intrv,enth_pr_intrv,temp_pr_intrv,satn_pr_intrv,iprsat[I]
     100.0     100.0     100.0     100.0     0
# PRINT 2
# ..  dens_pr_intrv,iprden[I],vis_pr_intrv,iprvis[I],potential_pr_intrv,iprpot[I]
     100.    11     0    0     0    0
# PRINT 3
# .. velocity_pr_intrv,iprvel[I],bcflow_pr_intrv,iprbcflow[I],source_pr_intrv,iprsource[I]
     0.0    0     0    0     0    0
# PRINT 4
# .. pm_properties_pr_intrv,iprmprop[I],poros_pr_intrv,permeability_pr_intrv
     0    0     0    0
# PRINT 5
# .. balance_pr_intrv,dimno_pr_intrv,iprdimno[I],residuals_pr_intrv,dump_pr_intrv
     0     0     0     0     10000
# PRINT 6
# .. plotscalar_pr_intrv,plotvector_pr_intrv,plotfile_type[I],time_series_pr_intrv
     100.0     100.0     5     0
# SLICE number
# ..  start_location: [TOP|BOTTOM]
# ..  index_of_rock_type(i,k),i=1 to nx, k=1 to nz (or k=nz to 1)
TOP
16*2
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
# UNCONFINED
# ..  unconfined[T/F],patm
     F     0
# SATURATION FUNCTION
#     LINEAR
# ..A pwb,pwr; (O) -  kodrp = 0 and unconfined
#     COREY
# ..B pwb,lambda[I]; (O) - kodrp = 1 and unconfined
#     COOLEY
# ..C pwb,bb,cc; (O) - kodrp = 3 and unconfined
# ----------------------------------
# Start of Keyword Data Block section
# PARAMETER
# Sample format for data block
# ..  array_identifier,(units)
# ..  format_specifier
# ..  data record set
# ........
XSPACING (m)
CONSTANT
5.0
ZSPACING (m)
CONSTANT
{deltaZ:.1f}
POROSITY (-)
ROCK
1 {HPporo:.2f}
2 {LPporo:.2f}
XPERMEABILITY (cm^2)
ROCK
1 {HPperm:2.1e}
2 {LPperm:2.1e}
ZPERMEABILITY (cm^2)
CALC
1.0e0
THERMAL_COND (erg/s-cm-K)
CONSTANT
1.5e05
CONDUCTIVE Heat Input along Base (mW/m^2)
FREE_FORMAT
16*0.0
PRESSURE (dyne/cm^2)
CALC
1 1.013e6
TEMPERATURE (C)
CALC
1, 30., {Tbot:.1f}
# End of Keyword Data Block section
# ---------------------------------
# Transient data
# TIME PERIOD #1
# .. tchg,delt,nsrce[I],lchgpr[T/F],nparms[I]
 0.001     0.001     0     T     0
# PRINT 1
# ..  pres_pr_intrv,enth_pr_intrv,temp_pr_intrv,satn_pr_intrv,iprsat[I]
    500     500     500.0     500.0     3
# PRINT 2
# ..  dens_pr_intrv,iprden[I],vis_pr_intrv,iprvis[I],potential_pr_intrv,iprpot[I]
     500    11     0    0     0    0
# PRINT 3
# .. velocity_pr_intrv,iprvel[I],bcflow_pr_intrv,iprbcflow[I],source_pr_intrv,iprsource[I]
     0.0    0     0    0     0    0
# PRINT 4
# .. pm_properties_pr_intrv,iprmprop[I],poros_pr_intrv,permeability_pr_intrv
     0    0     0    0
# PRINT 5
# .. balance_pr_intrv,dimno_pr_intrv,iprdimno[I],residuals_pr_intrv,dump_pr_intrv
     0     0     0     0     10000
# PRINT 6
# .. plotscalar_pr_intrv,plotvector_pr_intrv,plotfile_type[I],time_series_pr_intrv
     100   100     5     0
# TIME PERIOD #2  (instantaneous emplacement of magma)
# .. tchg,delt,nsrce[I],lchgpr[T/F],nparms[I]
    2.0     -1     -1     F     1
CONDUCTIVE Heat Input along Base (mW/m^2)
FREE_FORMAT
16*{heat:.1f}
# TIME PERIOD #3  (instantaneous emplacement of magma)
# .. tchg,delt,nsrce[I],lchgpr[T/F],nparms[I]
    1000.0     -1     -1     F     0
# TIME PERIOD #: End of simulation record
# ..  tchg
-1 /
#  .. End of the data file
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
""".format(deltaZ=dZ, heat=flux, Tbot=Tbot, HPporo=HPporo, LPporo=LPporo, LPperm=LPperm, HPperm=HPperm)
    return texto


def modify_input(flux, Tbot, HPperm, HPporo, LPperm, LPporo, dZ):
    texto ="""# Hydrotherm Data-Input Form
# Version 3.2
#      Notes:
#      Syntax of input line is denoted by # ..
#      A suffix letter indicates an exclusive record choice must be made.
#             i.e. A or B or C
#      (O) - Optional data with conditions for requirement
#      [a|b] - Indicates that either A or B must be entered
#      {{nnn}} - Indicates that the default number, nnn, is used if a zero
#              is entered for that variable
#      [T/F] - Indicates a logical variable
#      [I] - Indicates an integer variable
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#   Start of the data file
# TITLE
# ..  TITLE:   title line 1
# ..  VER3.2:  title line 2
TITLE  Nisi_linke system 
VER3.2 hydrothermal system without input
# DIMENSIONS
# ..  nx[I],ny[I],nz[I],tmstr,iyear[I]
16    1   16   0.0   2
# TIME STEP
# ..A  tstepmx[I],tmincr{{1.3}},%PchgTS{{10.}},%HchgTS{{5.}},Wsatn{{0.03}},relax{{1.}},
# ..       deltmin,n_tmcut{{20}}; confined flow
     1000     1.3     10.0     7.0     0.1     1.0     1.0E-10     40
# N-R ITERATION
# ..A  itermx{{10}},resmas{{1.e-12}},reseng{{1.e-2}},nphasechg{{10}},%PchgNR{{0.1}},
# ..        %HchgNR{{0.1}}; confined flow
     40     1.0E-12     1.0E-4     32     0.1     0.1
# LINEAR EQUATION SOLVER
# .. slmeth[I]
1
#    SSOR SOLVER CONTROL
# ..A  ssormx[I]{{10}},tol{{0.05}},ssorw0{{1.0}}; (O) - slmeth = 1
    1   0.01   1.0
# WEIGHTING AND AVERAGING
# ..  ioptupst[T/F],potdif{{0.0002}},theta{{1.0}}
     T     2.0E-5     1.0
# RELATIVE PERMEABILITY
# ..  kodrp[I],wsatn{{0.3}},ssatn{{0.0}}
     0     0.3     0.05
# ROCK PROPERTIES
# ..  heatcpty{{1.e7}},rxden,rxcmprss{{0.}},grav{{981.}},initphi[T/F]
# ..  initphi_filename{{IC_pressporo.xxx}}; (O) - initphi = T
   1.0e7  2.5,    0,   981.0,   F
# CYLINDRICAL COORDINATES
# ..  irad[T/F],wellrad,radiusmax
     T    0.1    1600
# PRINT/PLOT CONTROL
# ..  wide[I],long[I],v-avg[I]
     2     1     0
# PRINT 1
# ..  pres_pr_intrv,enth_pr_intrv,temp_pr_intrv,satn_pr_intrv,iprsat[I]
     100.0     100.0     100.0     100.0     0
# PRINT 2
# ..  dens_pr_intrv,iprden[I],vis_pr_intrv,iprvis[I],potential_pr_intrv,iprpot[I]
     100.    11     0    0     0    0
# PRINT 3
# .. velocity_pr_intrv,iprvel[I],bcflow_pr_intrv,iprbcflow[I],source_pr_intrv,iprsource[I]
     0.0    0     0    0     0    0
# PRINT 4
# .. pm_properties_pr_intrv,iprmprop[I],poros_pr_intrv,permeability_pr_intrv
     0    0     0    0
# PRINT 5
# .. balance_pr_intrv,dimno_pr_intrv,iprdimno[I],residuals_pr_intrv,dump_pr_intrv
     0     0     0     0     10000
# PRINT 6
# .. plotscalar_pr_intrv,plotvector_pr_intrv,plotfile_type[I],time_series_pr_intrv
     100.0     400.0     5     0
# SLICE number
# ..  start_location: [TOP|BOTTOM]
# ..  index_of_rock_type(i,k),i=1 to nx, k=1 to nz (or k=nz to 1)
TOP
16*2
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
16*1
# UNCONFINED
# ..  unconfined[T/F],patm
     F     0
# SATURATION FUNCTION
#     LINEAR
# ..A pwb,pwr; (O) -  kodrp = 0 and unconfined
#     COREY
# ..B pwb,lambda[I]; (O) - kodrp = 1 and unconfined
#     COOLEY
# ..C pwb,bb,cc; (O) - kodrp = 3 and unconfined
# ----------------------------------
# Start of Keyword Data Block section
# PARAMETER
# Sample format for data block
# ..  array_identifier,(units)
# ..  format_specifier
# ..  data record set
# ........
XSPACING (m)
CONSTANT
10.0
ZSPACING (m)
CONSTANT
{deltaZ:.1f}
POROSITY (-)
ROCK
1 {HPporo:.2f}
2 {LPporo:.2f}
XPERMEABILITY (cm^2)
ROCK
1 {HPperm:2.1e}
2 {LPperm:2.1e}
ZPERMEABILITY (cm^2)
CALC
1.0e0
THERMAL_COND (erg/s-cm-K)
CONSTANT
1.5e05
CONDUCTIVE Heat Input along Base (mW/m^2)
FREE_FORMAT
16*0.0
PRESSURE (dyne/cm^2)
CALC
1 1.013e6
TEMPERATURE (C)
CALC
1, 30., {Tbot:.1f}
# End of Keyword Data Block section
# ---------------------------------
# Transient data
# TIME PERIOD #1
# .. tchg,delt,nsrce[I],lchgpr[T/F],nparms[I]
 1.00     0.001     0     T     0
# PRINT 1
# ..  pres_pr_intrv,enth_pr_intrv,temp_pr_intrv,satn_pr_intrv,iprsat[I]
    1     1     1.0     1.0     3
# PRINT 2
# ..  dens_pr_intrv,iprden[I],vis_pr_intrv,iprvis[I],potential_pr_intrv,iprpot[I]
     1    11     0    0     0    0
# PRINT 3
# .. velocity_pr_intrv,iprvel[I],bcflow_pr_intrv,iprbcflow[I],source_pr_intrv,iprsource[I]
     1.0    0     0    0     0    0
# PRINT 4
# .. pm_properties_pr_intrv,iprmprop[I],poros_pr_intrv,permeability_pr_intrv
     1    1     1000    0
# PRINT 5
# .. balance_pr_intrv,dimno_pr_intrv,iprdimno[I],residuals_pr_intrv,dump_pr_intrv
     1     1     0     0     10000
# PRINT 6
# .. plotscalar_pr_intrv,plotvector_pr_intrv,plotfile_type[I],time_series_pr_intrv
     1   1    5     0
# TIME PERIOD #2  (instantaneous emplacement of magma)
# .. tchg,delt,nsrce[I],lchgpr[T/F],nparms[I]
    2.0     -1     -1     F     1
CONDUCTIVE Heat Input along Base (mW/m^2)
FREE_FORMAT
16*{heat:.1f}
# TIME PERIOD #3  (instantaneous emplacement of magma)
# .. tchg,delt,nsrce[I],lchgpr[T/F],nparms[I]
    100.0     -1     -1     F     0
# TIME PERIOD #: End of simulation record
# ..  tchg
-1 /
#  .. End of the data file
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
""".format(deltaZ=dZ, heat=flux, Tbot=Tbot, HPporo=HPporo, LPporo=LPporo, LPperm=LPperm, HPperm=HPperm)
    return texto



# print(modify_input(1000, 300, 1e-11, 0.15, 1e-18,  0.01, 100))