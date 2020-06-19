# GRACEEarthquakeModeling

The code here provides the calculation of coseismic GRACE gravity change, the modeling of GRACE gravity change caused by an earthqauke, and the localized spectral analysis. 

Notice: The shared code is for academic research. They are not user friendly. 
Some of the subroutines are not shared due to copyright limit. 


References:
Dai, C., C. Shum, J. Guo, K. Shang, B. Tapley, R. Wang, Improved source parameter constraints for five undersea earthquakes from north component of GRACE gravity and gravity gradient change measurements, Earth Planet. Sci. Lett., 443, 118-128, 2016.

Dai, C., C. Shum, R. Wang, L. Wang, J. Guo, K. Shang, and B. Tapley, Improved constraints on seismic source parameters of the 2011 Tohoku earthquake from GRACE gravity and gravity gradient changes, Geophys. Res. Lett., 41, doi:10.1002/2013GL059178, 2014.

Li, J., Chen, J.L. and Wilson, C.R., 2016. Topographic effects on coseismic gravity change for the 2011 Tohoku‐Oki earthquake and comparison with GRACE. Journal of Geophysical Research: Solid Earth, 121(7), pp.5509-5537.

Wahr, J., Molenaar, M. and Bryan, F., 1998. Time variability of the Earth's gravity field: Hydrological and oceanic effects and their possible detection using GRACE. Journal of Geophysical Research: Solid Earth, 103(B12), pp.30205-30229.

Guo, J.Y., Shum, C.K., 2009, Application of the cos-Fourier expansion to data transformation between different latitude-longitude grids, Computers & Geosciences, 35, 1439-1444, doi:10.1016/j.cageo.2008.09.010


SH2psd.scr: a shell script for modeling gravity change and calculating localized spectra.


Part one: calculation of coseismic GRACE gravity change.

1\ Preparation

Ggrid.f90: generate a regular grid, e.g. grid_Gauss_bw900WEN in ex1.zip file.
changenameCSR.sh: modify CSR GRACE L2 product name.
            e.g. "./changenameCSR.sh ifilelist renameodir filelist"
            where ifilelist is a text file that have a list of GRACE filenames in it, here are two example lines (see ifilelist in ex1.zip file):
            GSM-2_2003121-2003141_0021_UTCSR_0060_0005
            GSM-2_2012336-2012366_0029_EIGEN_G---_0005
            filelist is the output list of files (see filelist in ex1.zip file).
            
Subtract_Reference_SHCs_wwoH_NMAXncut_wstdM_SLRdC20.f: subtract a reference field (in terms of spherical harmonic coefficients).
            e.g. 
            "./Subtract_Reference_SHCs_wwoH_NMAXncut_wstdM_SLRdC20 >  out_Subtract_Reference_SHCs <<EOF
&parm
   IDIR="./IDIR/", 
   ODIR="./ODIR/",
   filenames="Subtract_filelist",
   NMAX=60,
   ncut=1,
   C20file="TN-07_C20_SLR_match.txt",
   ostdflag=1,
   C20flag=0,
   /
EOF
"
            where, IDIR, ODIR are input output directories. 
            Subtract_filelist is an input file (see Subtract_filelist in ex1.zip file).
            NMAX is the maximum degree, e.g. 60. 
            All degrees and orders less than ncut are set to zero, so ncut can be zero or 1.
            C20file is the file that contains the accurate C20 time series (see TN-07_C20_SLR_match.txt in ex1.zip file). 
            C20flag, a flag (e.g. 0). If C20flag is 1, the code will use the C20file, if it is zero, it will not use the file.
            ostdflag, a flag (e.g. 1) for choosing to output the uncertainties or not. 1 means to write uncercainties.

read_SHstd_wwoH.f: read spherical harmonic coefficients.
plotSHseries_std_io.m: calculating the a posteriori variance of unit weight based on coefficient time series. A subroutine for plotting shaded error bars is here: https://github.com/raacampbell/shadedErrorBar.

2\ Calculate gravity and gravity gradient from coefficients.
GNR_STATC_fast_JP_var_gra_bp.f: calculate North, East, Down components of gravity from coefficients.
AP_gra.f, subAPgra.f: subroutines.
lgdr2.f90: Not shared; calculate the FULLY NORMALIZED LEGENDRE FUNCTION WITH GIVEN ORDER(M), MAXIMUM DEGREE (NMAX), AND GIVEN EVALUATION POINT, T (COSINES OF COLATITUDE).
GNR_STATC_fast_JP_var_grad_bp.f: calculate the six components of gravity gradient from coefficients.
subAPgrad.f, AP.f: subroutines.

3\ Calculate the jump from time series.
plotGRAseriesFIT_wperiods_io.m: Calculate the jump from time series.


Part two: Modeling of GRACE gravity change caused by an earthqauke.

1\ Get crust model.
getCN2point_dai_ave_weight.f: Get the 1D crust model for a local region. Need the latest global CRUST model and software. 
Slipfmtstopscmp.f90: transform earthquake slip models to the format for PSCMP.
Epicenter2CE_io.m: retrieve geocentric distance of the fault plane.

2\ Calculate gravity due to Solid Earth deformation using Dr. Rongjiang Wang's code PSGRN/PSCMP.
psgrn08_io.inp, pscmp08_4SA_crustWei.inp: examples of input.
crust2fmt_io.m: transform the format of crust model.

3\ Get modeled gravity and gravity gradient on R0 (Earth's semi-major axis).
ETOPO1_Bed_g_int.xyz: ETOPO1 topography (optional).
OF_bw900_fmt.txt: global ocean function (1 for ocean, 0 for land).
plotgr_flat_disp_4inv_bw900_topo.f90: Calculate gravity change at a space-fixed location, surface density change by ocean response, and surface density change by topographic effect.
            Please do not use the topography correction term (FYI, J Li, JL Chen, CR Wilson, 2016).
Mass2SHCs.lnx: inversion of coefficients from surface density change. For code, ask Dr. Junyi Guo at OSU. Guo and Shum (2009).
Grid2SHCsLC.exe: inversion of coefficients from gravity. For code, ask Dr. Junyi Guo at OSU.
Trs_graD2potential_bw450_n_a1_4sigma_sum.f: transform of above two sets of coeffients to potential coefficients. For more information, see Wahr et a., 1998.


Part three: localized spectra based on Slepian basis function.
Grid2PSD_io.m: get localized spectra from spatial signals. See Dr. Frederik Simons's code for subroutines (http://geoweb.princeton.edu/people/simons/software.html).


Please acknowledge this code in publications or academic journals by citing:
Dai, C., C. Shum, J. Guo, K. Shang, B. Tapley, R. Wang, Improved source parameter constraints for five undersea earthquakes from north component of GRACE gravity and gravity gradient change measurements, Earth Planet. Sci. Lett., 443, 118-128, 2016.



