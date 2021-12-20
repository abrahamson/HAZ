c  ***** PEER NGASubduction Models (2019/2020) **********

c ----------------------------------------------------------------------
      subroutine S35_AG2020 ( mag, evType, rRup, vs30, z25, ztor, region, mu, sigma, phi, tau,
     1            rockPGA, specT, period2, iFlag, ACadjfac, epiflag )

C     Model Version: March 16, 2021 - PEER Report Version
C     PGA (T=0.0) set equal to T=0.01 coefficients

      implicit none
      integer MAXPER, NPER
      real sigma, ACadjfac
      parameter (MAXPER=25)

c      character*80 RegName(8), RegName1

      integer i, j, k, region, i9, iParam, count1, count2, iflag
      real mag, rRup, ZTOR, evType
      real rockPGA, period2, period(MAXPER), specT
      real mu, r
      real VLin(MAXPER), b_soil(MAXPER), c1, vs
      real vs30, sum1, angle1, taperTheta
      real n, c, T1, T2
      real vs1
      real Z25ref
c      integer iprint
      real c1_slab(8), c1_inter(MAXPER), depthLimit, c4
      real z25
      real term1, term1a, term1b, term1c, term2, term3, term4,
     1    term5, term6, term7, term6a
      real theta(45), part(45), NL_soil, mu1
      real a1(MAXPER), a2(MAXPER), a3(MAXPER), a4(MAXPER), a5(MAXPER),
     1     a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER), a10(MAXPER),
     1     a11(MAXPER), a12(MAXPER), a13(MAXPER), a14(MAXPER), a15(MAXPER),
     1     a16(MAXPER), a17(MAXPER), a18(MAXPER), a19(MAXPER), a20(MAXPER),
     1     a21(MAXPER), a22(MAXPER), a23(MAXPER), a24(MAXPER), a25(MAXPER),
     1     a26(MAXPER), a27(MAXPER), a28(MAXPER), a29(MAXPER), a30(MAXPER),
     1     a31(MAXPER), a32(MAXPER), a33(MAXPER), a34(MAXPER), a35(MAXPER),
     1     a36(MAXPER), a37(MAXPER), a38(MAXPER), a39(MAXPER), a40(MAXPER),
     1     a41(MAXPER), a42(MAXPER), a43(MAXPER), a44(MAXPER), a45(MAXPER)
      real AKFac(MAXPER), CasFac(MAXPER), AKFacT, CasFacT
      real rhow(MAXPER), rhob(MAXPER), rhowT, rhobT
      real c1_inter_T, vLin_T, b_soil_T, temp1
      real d1(MAXPER), d2(MAXPER), d1_T, d2_T
      real alpha1, alpha2, f2, f3, phi, tau, phi1, phi2, phi3
      real d0, d3, d4, d5, d6
      real a1Global(MAXPER), a1G_T
      real dAmp_dPGA, taulinPGA, taulin, f2PGA, f3PGA, philinPGA, phiamp, philin
      real phiB, phiBPGA, phiNL, tauNL
      real e1T, e2T, e3T, e1(MAXPER), e2(MAXPER), e3(MAXPER), Cepi
      integer epiflag
      real T1_phi2, T2_phi2, T3_phi2, T4_phi2, T1_phi3, T2_phi3, T3_phi3, T4_phi3
      real d3_phi2, d4_phi2, d5_phi2, d3_phi3, d4_phi3, d5_phi3, alpha_phi3
      real tau_lin, tau_lin_PGA, phi1_sq_100, phi1_sq_PGA_100, phi1_sq, phi1_sq_PGA
      real A_phi2_100, Alpha_phi2, A_phi2, A_phi3, f2_PGA, f3_PGA, phi_lin, phi_lin_PGA, phi_S2S, phi_S2S_PGA
      real phiSS_lin, phiSS_lin_PGA, partial_f_PGA, phi_Amp, Phi_B_PGA, Phi_B
      real phiSQ_NL, tauSQ_NL, phiSS_B_PGA, PhiSS_B, phiSS_sq_NL, phiSS, sigmaSS

c     slab c1: Alaska, CAS, CenAm, Japan, NewZealand, SouthAm, Taiwan, global
      data c1_slab / 7.9, 7.1, 7.4, 7.6, 8.0, 7.5, 7.7, 7.5 /

c      data RegName / 'Alaska', 'Cascadia', 'Central_America', 'Japan',
c     1       'New_Zealand', 'South_America', 'Taiwan', 'Global' /

C     Updated Coefficients 3/16/21 (Consistent with PEER Report)
      data period / 0.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1              0.5, 0.6, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.5, 10.0 /
      data vLin / 865.1, 865.1, 865.1, 907.8, 1053.5, 1085.7, 1032.5, 877.6, 748.2, 654.3, 587.1,
     1            503.0, 456.6, 430.3, 410.5, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0,
     2            400.0, 400.0, 400.0  /
      data b_soil / -1.186, -1.186, -1.219, -1.273, -1.346, -1.471, -1.624, -1.931, -2.188,
     1              -2.381, -2.518, -2.657, -2.669, -2.599, -2.401, -1.955, -1.025, -0.299,
     2               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /
      data c1_inter / 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2, 8.2,
     1                8.15, 8.1, 8.05, 8.0, 7.95, 7.9, 7.85, 7.8, 7.8, 7.8, 7.8 /
      data A1Global  / 4.596, 4.596, 4.678, 4.773, 5.029, 5.334, 5.455, 5.376, 4.936,
     1                 4.636, 4.423, 4.124, 3.838, 3.562, 3.152, 2.544, 1.636, 1.076,
     2                 0.658, 0.424, 0.093, -0.145, -0.32, -0.556, -0.86 /

      data a2 / -1.45, -1.45, -1.45, -1.45, -1.45, -1.45, -1.45, -1.425, -1.335, -1.275,
     1          -1.231, -1.165, -1.115, -1.071, -1.02, -0.95, -0.86, -0.82, -0.798, -0.793,
     2          -0.793, -0.793, -0.793, -0.793, -0.793 /
      data a6 / -0.0043, -0.0043, -0.0043, -0.0044, -0.0046, -0.0047, -0.0048, -0.0047,
     1          -0.0045, -0.0043, -0.0042, -0.004, -0.0037, -0.0035, -0.0032, -0.0029,
     2          -0.0026, -0.0024, -0.0022, -0.0021, -0.002, -0.002, -0.002, -0.002, -0.002 /
      data a7 / 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21,
     1          3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.13, 2.985, 2.818, 2.682, 2.515, 2.3 /
      data a8 / 0.044, 0.044, 0.044, 0.044, 0.044, 0.044, 0.044, 0.044, 0.043, 0.042, 0.041,
     1          0.04, 0.039, 0.038, 0.037, 0.035, 0.034, 0.032, 0.031, 0.03, 0.029, 0.028,
     2          0.027, 0.026, 0.025 /
      data a10 / 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21,
     1           3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.21, 3.13, 2.985, 2.818, 2.682, 2.515, 2.3 /
      data a11 / 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.0062, 0.0056,
     1           0.0051, 0.0043, 0.0037, 0.0033, 0.0027, 0.0019, 0.0008, 0.0, 0.0,
     2           0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /
      data a12 / 0.9, 0.9, 1.008, 1.127, 1.333, 1.565, 1.679, 1.853, 2.022, 2.181, 2.281,
     1           2.379, 2.339, 2.217, 1.946, 1.416, 0.394, -0.417, -0.725, -0.695, -0.638,
     2          -0.597, -0.561, -0.53, -0.486 /
      data a13 / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.002, -0.007, -0.011,
     1          -0.015, -0.021, -0.028, -0.041, -0.05, -0.057, -0.065, -0.077, -0.088,
     2          -0.098, -0.11, -0.127 /
      data a14 / -0.46, -0.46, -0.46, -0.46, -0.46, -0.46, -0.46, -0.46, -0.46, -0.46,
     1           -0.46, -0.47, -0.48, -0.49, -0.5, -0.51, -0.52, -0.53, -0.54, -0.54,
     2           -0.54, -0.54, -0.54, -0.54, -0.54 /
C      data a16 / 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.08412, 0.08038,
C     1           0.078, 0.07463, 0.07213, 0.06975, 0.06713, 0.06325, 0.05887, 0.05875,
C     2           0.05975, 0.05912, 0.05037, 0.04287, 0.03788, 0.03162, 0.02413 /
      data a16 / 0.090, 0.090, 0.090, 0.090, 0.090, 0.090, 0.090, 0.090, 0.084, 0.080,
     1           0.078, 0.075, 0.072, 0.070, 0.067, 0.063, 0.059, 0.059, 0.060, 0.059,
     2           0.050, 0.043, 0.038, 0.032, 0.024 /
      data a18 / -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.186, -0.15, -0.14, -0.12,
     1           -0.1, -0.08, -0.06, -0.047, -0.035, -0.018, -0.01, -0.005, 0.0, 0.0,
     2            0.0, 0.0, 0.0, 0.0 /
      data a20 / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.055, -0.105, -0.134, -0.15,
     1          -0.15, -0.15, -0.15, -0.15, -0.15, -0.13, -0.11, -0.095, -0.085,
     2          -0.073, -0.065, -0.06, -0.055, -0.045 /
      data a21 / 0.04, 0.04, 0.04, 0.04, 0.04, 0.06, 0.1, 0.135, 0.17, 0.17, 0.17, 0.17,
     1           0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17 /
      data a22 / 0.04, 0.04, 0.04, 0.04, 0.04, 0.06, 0.1, 0.135, 0.17, 0.17, 0.17, 0.17,
     1           0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17 /
      data a23 / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.069, 0.14, 0.164, 0.19, 0.206, 0.22,
     1           0.225, 0.217, 0.185, 0.083, 0.045, 0.026, 0.035, 0.053, 0.072,
     2           0.086, 0.115, 0.151 /
      data a24 / 0.0015, 0.0015, 0.0015, 0.0015, 0.0011, 0.0011, 0.0012, 0.0013, 0.0013,
     1           0.0013, 0.0014, 0.0015, 0.0015, 0.0015, 0.0014, 0.0013, 0.0014, 0.0015,
     2           0.0014, 0.0014, 0.0014, 0.0014, 0.0014, 0.0014, 0.0014 /
      data a25 / 0.0007, 0.0007, 0.0006, 0.0006, 0.0006, 0.0004, 0.0003, -0.0002, -0.0007,
     1          -0.0009, -0.001, -0.001, -0.0011, -0.0012, -0.0011, -0.0008, -0.0004, 0.0002,
     2           0.0004, 0.0007, 0.001, 0.0013, 0.0015, 0.0017, 0.0017 /
      data a26 / 0.0036, 0.0036, 0.0036, 0.0037, 0.0039, 0.0039, 0.0039, 0.0037, 0.0031,
     1           0.0027, 0.002, 0.0013, 0.0009, 0.0006, 0.0003, 0.0001, -0.0001, 0.0,
     2           0.0, 0.0003, 0.0007, 0.0014, 0.0015, 0.0015, 0.0015 /
      data a27 / -0.0004, -0.0004, -0.0005, -0.0007, -0.0009, -0.0009, -0.0008, -0.0009,
     1           -0.001, -0.0011, -0.0009, -0.0007, -0.0007, -0.0007, -0.0007, -0.0008,
     2           -0.0008, -0.0007, -0.0007, -0.0007, -0.0006, -0.0004, -0.0003, -0.0002, -0.0001 /
      data a28 / 0.0025, 0.0025, 0.0025, 0.0025, 0.0026, 0.0026, 0.0026, 0.0022, 0.0018,
     1           0.00155, 0.0014, 0.0011, 0.0008, 0.0006, 0.0004, 0.0002, 0.0001, 0.0002,
     2           0.0002, 0.0004, 0.0006, 0.0008, 0.0011, 0.0014, 0.0017 /
c      data a29 / 0.0006, 0.0006, 0.0005, 0.0005, 0.0004, 0.0003, 0.0003, 0.0001, -0.0001,
c     1          -0.0003, -0.0002, 0.0, 0.0002, 0.0002, 0.0002, 0.0001, 0.0, 0.0, -0.00015,
c     2          -0.0002, -0.00015, -0.00005, 0.0, 0.0001, 0.0002 /
      data a29 / 0.0006, 0.0006, 0.0005, 0.0005, 0.0004, 0.0003, 0.0003, 0.0001, -0.0001,
     1          -0.0003, -0.0002, 0.0, 0.0002, 0.0002, 0.0002, 0.0001, 0.0, 0.0, -0.00015,
     2          -0.0002, -0.0002, -0.0001, 0.0, 0.0001, 0.0002 /
      data a30 / 0.0033, 0.0033, 0.0033, 0.0034, 0.0036, 0.0037, 0.0038, 0.0037, 0.0035,
     1           0.0033, 0.0032, 0.003, 0.0027, 0.0025, 0.0022, 0.0019, 0.0016, 0.0014,
     2           0.0012, 0.0011, 0.001, 0.001, 0.001, 0.001, 0.001 /
      data a31 / 3.77831, 3.77831, 3.82805, 3.89327, 4.2867, 4.59402, 4.70773, 4.6065,
     1           4.1866, 3.85147, 3.57825, 3.24933, 2.98178, 2.77838, 2.47799, 1.92518,
     2           0.99242, 0.46763, 0.0579, -0.13913, -0.30296, -0.40939, -0.50103,
     3          -0.62091, -0.62208 /
      data a32 / 3.34676, 3.34676, 3.44008, 3.50871, 3.65526, 3.97985, 4.1312, 4.27372,
     1           3.96502, 3.68208, 3.5415, 3.32556, 3.13343, 2.92152, 2.53804, 1.96256,
     2           1.35679, 0.81801, 0.43889, 0.10455, -0.15972, -0.20626, -0.32229,
     3          -0.42231, -0.59085 /
      data a33 / 3.80252, 3.80252, 3.90532, 4.01892, 4.2952, 4.54636, 4.61378, 4.52897,
     1           4.1656, 3.91471, 3.78463, 3.57016, 3.35516, 3.09216, 2.65718, 2.14585,
     2           1.3499, 0.81478, 0.39788, 0.1046, -0.23242, -0.57223, -0.86305, -1.17731, -1.40701 /
      data a34 / 5.03612, 5.03612, 5.13749, 5.26986, 5.61569, 6.02041, 6.16254, 5.96141,
     1           5.39202, 5.01168, 4.70573, 4.28959, 3.93215, 3.61493, 3.17852, 2.57218,
     2           1.64993, 1.06575, 0.63103, 0.38818, 0.01644, -0.28023, -0.4822, -0.75661, -1.08704 /
      data a35 / 4.62715, 4.62715, 4.69578, 4.78092, 5.02105, 5.3474, 5.50651, 5.51804, 5.1668,
     1           4.87436, 4.6544, 4.36597, 4.07791, 3.81455, 3.4391, 2.80555, 1.85463, 1.30203,
     2           0.8017, 0.59579, 0.35222, 0.18743, -0.12426, -0.33163, -0.67834 /
      data a36 / 4.80444, 4.80444, 4.89425, 5.00277, 5.28185, 5.61229, 5.76676, 5.7313, 5.29426,
     1           5.00579, 4.75879, 4.3789, 4.03937, 3.73663, 3.29299, 2.64751, 1.68417, 1.10019,
     2           0.67373, 0.41262, 0.0097, -0.27152, -0.45911, -0.68215, -0.91734 /
      data a37 / 3.56691, 3.56691, 3.64249, 3.70633, 3.91843, 4.22072, 4.35356, 4.36641,
     1           4.01685, 3.75904, 3.59143, 3.37039, 3.15641, 2.95835, 2.65564, 2.06665,
     2           1.33157, 0.76066, 0.36483, 0.16882, -0.03231, -0.15158, -0.22172,
     3          -0.33381, -0.54408 /

      data a39  / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.101, 0.184, 0.315, 0.416,
     1            0.499, 0.6, 0.731, 0.748, 0.761, 0.77, 0.778, 0.79, 0.799, 0.807, 0.817, 0.829 /
      data a41  / -0.029, -0.029, -0.024, -0.034, -0.061, -0.076, -0.049, -0.026, -0.011, -0.009,
     1             0.005, 0.04, 0.097, 0.145, 0.197, 0.269, 0.347, 0.384, 0.397, 0.404, 0.397,
     2             0.378, 0.358, 0.333, 0.281 /

      data d1 / 0.325, 0.325, 0.325, 0.325, 0.325, 0.325, 0.325, 0.325, 0.325, 0.325,
     1          0.325, 0.325, 0.325, 0.325, 0.325, 0.325, 0.312, 0.302, 0.295, 0.289,
     2          0.28, 0.273, 0.267, 0.259, 0.25 /
      data d2 / 0.137, 0.137, 0.137, 0.137, 0.137, 0.137, 0.137, 0.137, 0.137, 0.137,
     1          0.137, 0.137, 0.137, 0.137, 0.137, 0.137, 0.113, 0.096, 0.082, 0.072,
     2          0.055, 0.041, 0.03, 0.017, 0.0 /
c      data AKfac / 0.487, 0.487, 0.519, 0.543, 0.435, 0.410, 0.397, 0.428, 0.442,
c     1             0.494, 0.565, 0.625, 0.634, 0.581, 0.497, 0.469, 0.509, 0.478,
c     2             0.492, 0.47, 0.336, 0.228, 0.151, 0.051, -0.251 /
c      data CasFac / 0.828, 0.828, 0.825, 0.834, 0.895, 0.863, 0.842, 0.737, 0.746,
c     1              0.796, 0.782, 0.768, 0.728, 0.701, 0.685, 0.642, 0.325, 0.257,
c     2              0.211, 0.296, 0.232, 0.034, -0.037, -0.178, -0.313 /
      data AKfac / 0.487, 0.487, 0.519, 0.543,  0.435, 0.410, 0.397, 0.428, 0.442, 0.494,
     2             0.565, 0.625, 0.634, 0.581,  0.498, 0.458, 0.499, 0.478, 0.495, 0.451,
     4             0.322, 0.194, 0.127, 0.051, -0.190 /
      data CasFac / 0.828, 0.828,  0.825,  0.834,  0.895, 0.863, 0.842, 0.737, 0.746, 0.796,
     1              0.782, 0.768,  0.728,  0.701,  0.686, 0.631, 0.315, 0.257, 0.214, 0.278,
     2              0.219, 0.001, -0.062, -0.177, -0.252 /

      data rhow / 1.0, 1.0, 0.99, 0.99, 0.97, 0.95, 0.92, 0.9, 0.87, 0.84, 0.82,
     1            0.74, 0.66, 0.59, 0.5, 0.41, 0.33, 0.3, 0.27, 0.25, 0.22, 0.19,
     2            0.17, 0.14, 0.1 /
      data rhob / 1.0, 1.0, 0.99, 0.99, 0.985, 0.98, 0.97, 0.96, 0.94, 0.93, 0.91,
     1            0.86, 0.8, 0.78, 0.73, 0.69, 0.62, 0.56, 0.52, 0.495, 0.43, 0.4,
     2            0.37, 0.32, 0.28 /
      data e1 / 0.550, 0.550, 0.550, 0.550, 0.560, 0.580, 0.590, 0.590, 0.570, 0.530,
     1          0.490, 0.425, 0.375, 0.345, 0.300, 0.240, 0.230, 0.230, 0.230, 0.240,
     2          0.270, 0.300, 0.320, 0.350, 0.350 /
      data e2 / -0.270, -0.270, -0.270, -0.270, -0.270, -0.270, -0.270, -0.270, -0.270,
     1          -0.224, -0.186, -0.126, -0.079, -0.041,  0.005,  0.065,  0.065,  0.065,
     2           0.065,  0.065,  0.065,  0.065,  0.065,  0.065,  0.065 /
      data e3 / 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.043,
     1          0.037, 0.028, 0.022, 0.016, 0.009, 0.000, 0.000, 0.000, 0.000, 0.000,
     2          0.000, 0.000, 0.000, 0.000, 0.000 /

c     set period-independent terms
      depthLimit = 200.
      c4 = 10.
      n = 1.18
      c = 1.88
      theta(1) = 0.
      theta(3) = 0.1
      theta(4) = 0.73
      theta(5) = 0.
      theta(9) = 0.4
c      theta(15) = 0.
      theta(15) = -0.1
      theta(17) = 0.
      theta(19) = 0.
      theta(38) = 0.
      theta(40) = 0.
      theta(42) = 0.
      theta(43) = 0.
      theta(44) = 0.
      theta(45) = 0.34
      d0 = 0.47
      d3 = 0.109
      d4 = 0.062
      d5 = 0.470
      d6 = 0.242
      alpha2 = 0.42

      nper = 25

C First check for the PGA case
      if (specT .eq. 0.0) then
         vlin_T = vLin(1)
         b_soil_T = b_soil(1)
         c1 = c1_inter(1)
         c1_inter_T = c1_inter(1)
c         if ( region .eq. 8 ) theta(1) = a1Global(1)
         if ( region .eq. 8 ) then
            theta(1) = a1Global(1)
         else
            theta(1) = 0.0
         endif
         theta(2) = a2(1)
         theta(6) = a6(1)
         theta(7) = a7(1)
         theta(8) = a8(1)
         theta(10) = a10(1)
         theta(11) = a11(1)
         theta(12) = a12(1)
         theta(13) = a13(1)
         theta(14) = a14(1)
         theta(16) = a16(1)
         theta(18) = a18(1)
         theta(20) = a20(1)
         theta(21) = a21(1)
         theta(22) = a22(1)
         theta(23) = a23(1)
         theta(24) = a24(1)
         theta(25) = a25(1)
         theta(26) = a26(1)
         theta(27) = a27(1)
         theta(28) = a28(1)
         theta(29) = a29(1)
         theta(30) = a30(1)
         theta(31) = a31(1)
         theta(32) = a32(1)
         theta(33) = a33(1)
         theta(34) = a34(1)
         theta(35) = a35(1)
         theta(36) = a36(1)
         theta(37) = a37(1)
         theta(39) = a39(1)
         theta(41) = a41(1)
         d1_T = d1(1)
         d2_T = d2(1)
         AKfacT = AKFac(1)
         CasfacT = CasFac(1)
         rhowT = rhow(1)
         rhobT = rhob(1)
         e1T = e1(1)
         e2T = e2(1)
         e3T = e3(1)
         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=2,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020
         endif
      enddo

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*)
      write (*,*) 'AG (2020) Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*)
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),b_soil(count1),b_soil(count2),
     +                   specT,b_soil_T,iflag)
      call S24_interp (period(count1),period(count2),vlin(count1),vlin(count2),
     +                   specT,vlin_T,iflag)
      call S24_interp (period(count1),period(count2),c1_inter(count1),c1_inter(count2),
     +                   specT,c1_inter_T,iflag)
      if (region .eq. 8) then
         call S24_interp (period(count1),period(count2),a1Global(count1),a1Global(count2),
     +                   specT,theta(1),iflag)
      endif
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,theta(2),iflag)
      call S24_interp (period(count1),period(count2),a6(count1),a6(count2),
     +                   specT,theta(6),iflag)
      call S24_interp (period(count1),period(count2),a7(count1),a7(count2),
     +                   specT,theta(7),iflag)
      call S24_interp (period(count1),period(count2),a8(count1),a8(count2),
     +                   specT,theta(8),iflag)
      call S24_interp (period(count1),period(count2),a10(count1),a10(count2),
     +                   specT,theta(10),iflag)
      call S24_interp (period(count1),period(count2),a11(count1),a11(count2),
     +                   specT,theta(11),iflag)
      call S24_interp (period(count1),period(count2),a12(count1),a12(count2),
     +                   specT,theta(12),iflag)
      call S24_interp (period(count1),period(count2),a13(count1),a13(count2),
     +                   specT,theta(13),iflag)
      call S24_interp (period(count1),period(count2),a14(count1),a14(count2),
     +                   specT,theta(14),iflag)
      call S24_interp (period(count1),period(count2),a16(count1),a16(count2),
     +                   specT,theta(16),iflag)
      call S24_interp (period(count1),period(count2),a18(count1),a18(count2),
     +                   specT,theta(18),iflag)
      call S24_interp (period(count1),period(count2),a20(count1),a20(count2),
     +                   specT,theta(20),iflag)
      call S24_interp (period(count1),period(count2),a21(count1),a21(count2),
     +                   specT,theta(21),iflag)
      call S24_interp (period(count1),period(count2),a22(count1),a22(count2),
     +                   specT,theta(22),iflag)
      call S24_interp (period(count1),period(count2),a23(count1),a23(count2),
     +                   specT,theta(23),iflag)
      call S24_interp (period(count1),period(count2),a24(count1),a24(count2),
     +                   specT,theta(24),iflag)
      call S24_interp (period(count1),period(count2),a25(count1),a25(count2),
     +                   specT,theta(25),iflag)
      call S24_interp (period(count1),period(count2),a26(count1),a26(count2),
     +                   specT,theta(26),iflag)
      call S24_interp (period(count1),period(count2),a27(count1),a27(count2),
     +                   specT,theta(27),iflag)
      call S24_interp (period(count1),period(count2),a28(count1),a28(count2),
     +                   specT,theta(28),iflag)
      call S24_interp (period(count1),period(count2),a29(count1),a29(count2),
     +                   specT,theta(29),iflag)
      call S24_interp (period(count1),period(count2),a30(count1),a30(count2),
     +                   specT,theta(30),iflag)
      call S24_interp (period(count1),period(count2),a31(count1),a31(count2),
     +                   specT,theta(31),iflag)
      call S24_interp (period(count1),period(count2),a32(count1),a32(count2),
     +                   specT,theta(32),iflag)
      call S24_interp (period(count1),period(count2),a33(count1),a33(count2),
     +                   specT,theta(33),iflag)
      call S24_interp (period(count1),period(count2),a34(count1),a34(count2),
     +                   specT,theta(34),iflag)
      call S24_interp (period(count1),period(count2),a35(count1),a35(count2),
     +                   specT,theta(35),iflag)
      call S24_interp (period(count1),period(count2),a36(count1),a36(count2),
     +                   specT,theta(36),iflag)
      call S24_interp (period(count1),period(count2),a37(count1),a37(count2),
     +                   specT,theta(37),iflag)
      call S24_interp (period(count1),period(count2),a39(count1),a39(count2),
     +                   specT,theta(39),iflag)
      call S24_interp (period(count1),period(count2),a41(count1),a41(count2),
     +                   specT,theta(41),iflag)
      call S24_interp (period(count1),period(count2),d1(count1),d1(count2),
     +                   specT,d1_T,iflag)
      call S24_interp (period(count1),period(count2),d2(count1),d2(count2),
     +                   specT,d2_T,iflag)
      call S24_interp (period(count1),period(count2),AKFac(count1),AKFac(count2),
     +                   specT,AKFacT,iflag)
      call S24_interp (period(count1),period(count2),CasFac(count1),CasFac(count2),
     +                   specT,CasFacT,iflag)
      call S24_interp (period(count1),period(count2),rhow(count1),rhow(count2),
     +                   specT,rhowT,iflag)
      call S24_interp (period(count1),period(count2),rhob(count1),rhob(count2),
     +                   specT,rhobT,iflag)
      call S24_interp (period(count1),period(count2),e1(count1),e1(count2),
     +                   specT,e1T,iflag)
      call S24_interp (period(count1),period(count2),e2(count1),e2(count2),
     +                   specT,e2T,iflag)
      call S24_interp (period(count1),period(count2),e3(count1),e3(count2),
     +                   specT,e3T,iflag)

 1011 period2 = specT


c     Set c1_slab
      if (evType .eq. 0. ) then
          c1 = c1_inter_T
      else
          c1 = c1_slab(region)
      endif

      r = rrup +  c4*exp((mag-6.)*theta(9))

c     Set break in slope of soil term
      VS1 = 1000.
      if ( vs30 .gt. VS1 ) then
        vs = VS1
      else
        vs = vs30
      endif
c      write (*,'( 2f10.3)') vs, vlin_T
c      pause 'vs'

c     set fixed non-linear soil
      if ( vs .lt. vLin_T) then
        nl_soil =  ( - b_soil_T*alog(c+rockPGA)
     1            + b_soil_T*alog(rockPga+c*(vs/vLin_T)**(n)) )
      else
        nl_soil = (b_soil_T*n*alog(vs/vLIN_T) )
      endif

c     Set partial derivatives with respect to coefficients
      part(1) = 1.
      part(2) = alog(r)
      part(3) = (mag - 7.) * alog(r)

      if ( mag .le. c1 ) then
        part(4) = mag - c1 + (c1_slab(region)-7.5) * evType
        part(5) = 0.
      else
        part(4) = (c1_slab(region)-7.5) * evType
        part(5) = mag - c1
      endif
      part(6) = rRup
      part(7) = 0.
      part(8) = 0.
      part(9) = 0.

      if ( region .eq. 3 .or. region .eq. 4 .or. region .eq. 6 ) then
        part(7) = 0.
        part(10) = evType
      else
        part(7) = evType
        part(10) = 0.
      endif

      if ( ztor .lt. depthLimit ) then
        if ( ztor .lt. 50. ) then
          part(8) = (ztor-50.) * evType
          part(11) = 0.
        else
          part(8) = 0.
          part(11) = evType * (ztor - 50.)
        endif
      else
         part(11) = evType * (depthLimit - 50.)
      endif

      part(12) = alog(vs/vLin_T)
      part(13) = (10. - mag )**2
      part(14) = alog(r) * evType

c      if ( region .eq. 3 )  then
c        part(15) =  alog(r)
c      else
c        part(15) = 0.
c      endif
C     Mainshock/Aftershock Flag not currently coded
      part(15) = 0.

      if ( region .eq. 7 )  then
        part(16) =  alog(r)
      else
        part(16) = 0.
      endif

c     Initialize regional terms
      do i9=17,45
        part(i9) = 0.
      enddo

c     Region-dependent VS30 terms
      if ( region .eq. 1 ) then
        part(17) =  alog(vs/vLin_T)
      elseif ( region .eq. 2 ) then
        part(18) =  alog(vs/vLin_T)
      elseif ( region .eq. 3 ) then
        part(19) =  alog(vs/vLin_T)
      elseif ( region .eq. 4 ) then
        part(20) =  alog(vs/vLin_T)
      elseif ( region .eq. 5 ) then
        part(21) =  alog(vs/vLin_T)
      elseif ( region .eq. 6 ) then
        part(22) =  alog(vs/vLin_T)
      elseif ( region .eq. 7 ) then
        part(23) =  alog(vs/vLin_T)
      endif

c     Region-dependent R terms (interface)
      if ( region .eq. 1 ) then
        part(24) =  rRup
      elseif ( region .eq. 2 ) then
        part(25) =  rRup
      elseif ( region .eq. 3 ) then
        part(26) =  rRup
      elseif ( region .eq. 4 ) then
        part(27) =  rRup
      elseif ( region .eq. 5 ) then
        part(28) =  rRup
      elseif ( region .eq. 6 ) then
        part(29) =  rRup
      elseif ( region .eq. 7 ) then
        part(30) =  rRup
      endif

c     Region-dependent Constant terms
      if ( region .eq. 1 ) then
        part(31) =  1.
      elseif ( region .eq. 2 ) then
        part(32) =  1.
      elseif ( region .eq. 3 ) then
        part(33) =  1.
      elseif ( region .eq. 4 ) then
        part(34) =  1.
      elseif ( region .eq. 5 ) then
        part(35) =  1.
      elseif ( region .eq. 6 ) then
        part(36) =  1.
      elseif ( region .eq. 7 ) then
        part(37) =  1.
      endif

c     set z2.5_ref
c     Calc reference z25 for given vs30 for Cascadia(Region=2) or Japan (Region=4)
      if ( Region .eq. 2 ) then
	    call S35_z25_interp ( vs30, 200., 570., 8.52, 7.6, z25ref )
      elseif ( Region .eq. 4 ) then
	    call S35_z25_interp ( vs30, 170., 800., 7.3, 4.1, z25ref )
      else
            z25Ref = -1.
      endif

c     Z2.5 scaling
      if ( z25 .ge. 0. .and. z25ref .ge. 0.) then
        temp1 = alog ( ( z25*1000.0 + 50.) /(z25ref+50.) )
      else
        temp1 = 0.
      endif
      if ( z25 .gt. 0. ) then
        if ( region .eq. 1 ) then
           part(38) =  temp1
        elseif ( region .eq. 2 ) then
        if ( temp1 .gt. 0. ) then
           part(39) =  temp1-0.
        else
           part(39) = 0.
        endif
        elseif ( region .eq. 3 ) then
           part(40) =  temp1
        elseif ( region .eq. 4 ) then
          if ( temp1 .gt. -2. ) then
            part(41) =  temp1
          else
c            part(41) = 0.0
            part(41) = -2.0
          endif
        elseif ( region .eq. 5 ) then
          part(42) =  temp1
        elseif ( region .eq. 6 ) then
          part(43) = temp1
        elseif ( region .eq. 7 ) then
          part(44) =  temp1
        endif
      endif

c     Add mag scaling difference for interface and slab
      part(45) = (c1_slab(region)-7.5) * evType
      if ( evType .eq. 1 .and. mag .le. c1 ) then
c      if ( evType .eq. 1 ) then
         part(45) =  mag - c1 + part(45)
      endif

c     comptue median
      mu = NL_soil
      do iParam=1,45
        mu = mu + theta(iParam) * part(iParam)
      enddo

C     Apply Alaska or Cascadia adjustment is requested.
      if (Region .eq. 1 ) then
         ACadjfac = AKfacT
      elseif (Region .eq. 2 ) then
         ACadjfac = CasfacT
      else
         ACadjfac = 0.0
      endif


      term1 = theta(1)*part(1) + theta(2)*part(2) + theta(3)*part(3)
     1            + theta(6)*part(6)

      term2 = theta(4)*part(4) +  theta(5)*part(5) + theta(13)*part(13)

      term3 = theta(10)*part(10) + theta(7)*part(7) + theta(14)*part(14)
      term4 = theta(45)*part(45)
      term5 = theta(8)*part(8) + theta(11)*part(11)
      term6 = theta(12)*part(12) + NL_soil
      term1c = theta(15)*part(15) + theta(16)*part(16)
      term6a = 0.
      do k=17,23
        term6a = term6a + theta(k)*part(k)
      enddo
      term1a = 0.
      do k=24,30
        term1a = term1a + theta(k)*part(k)
      enddo
      term1b = 0.
      do k=31,37
        term1b = term1b + theta(k)*part(k)
      enddo
      term7 = 0.
      do k=38,44
        term7 = term7 + theta(k)*part(k)
      enddo

      mu1 = term1 + term2 + term3 + term4 + term5 + term6 + term6a + term7
     1      + term1a + term1b + term1c

C     Apply Epistemic model only for Global version
      if (Region .eq. 8) then
         cepi = e1T + e2T*(Rrup/100.0) + e3T*(Rrup/100.0)*(Rrup/100.0)
         if (epiflag .eq. -1) then
            mu = mu - cepi
            mu1 = mu1 - cepi
         elseif (epiflag .eq. 1) then
            mu = mu + cepi
            mu1 = mu1 + cepi
         endif
      endif

c      if (specT .eq. 0.0 .and. vs30 .eq. 760.0) then
c          write (*,'( i5,30f10.4))')
c     1     region, mag, ZTOR, evType, rrup, vs30, z25, z25ref,
c     1      rockpga, specT, mu, exp(mu1),
c     1      term1, term1a , term1b , term1c, term2, term3, term4, term5,
c     1      term6, term6a, term7
c
c      endif


c      if (iprint .eq. 1 ) write (33,'( i5,30f10.4))')
c     1     region, mag, ZTOR, evType, rrup, vs30, z25,
c     1      rockpga, mu,
c     1      term1, term2, term3, term4, term5, term6, term6a, term7,
c     1      term1a , term1b , term1c, z25ref, temp1,
c     1      alog ( ( z25 + 50.) /(z25ref+50.) )

c       regName1 = regName(region)

c     Sigma models

c     set coefficents
      d0 = 0.47
      T1_phi2 = 0.03
      T2_Phi2 = 0.075
      T3_phi2 = 0.20
      T4_phi2 = 1.0
      T1_phi3 = 0.03
      T2_Phi3 = 0.075
      T3_phi3 = 0.10
      T4_phi3 = 0.3
      d3_phi2 = 0.109
      d4_phi2 = 0.062
      d5_phi2 = 0.470
      d3_phi3 = 0.242
      d4_phi3 = 0.0
      d5_phi3 = 0.0
      alpha_phi3 = 0.42

C     Set alpha for phi2 (eq 5.6)
      if (rrup .le. 250. ) then
        alpha_phi2 = 1.
      else if ( rrup .lt. 450. ) then
        alpha_phi2 = 1 - 0.0036*(Rrup-250.)
      else
        alpha_phi2 = 0.28
      endif

c     tau model
      tau_LIN = d0
      tau_LIN_PGA = d0

c     phi1 model
      phi1_sq_100 = d1_T
      phi1_sq_PGA_100 = d1(1)
      if ( rrup .lt. 150. ) then
        phi1_SQ = d1_T
        phi1_SQ_PGA = d1(1)
      elseif ( Rrup .lt. 450. ) then
        phi1_SQ = d1_T + d2_T * (Rrup-150)/300.
        phi1_SQ_PGA = d1(1) + d2(1) * (Rrup-150)/300.
      else
        phi1_SQ = d1_T + d2_T
        phi1_SQ_PGA = d1(1) + d2(1)
      endif

c     phi2 model
c     A_phi term (eq 5.4)
      A_phi2_100 = d3_phi2
      if ( Rrup .lt. 225. ) then
        A_phi2 = d3_phi2
      elseif ( Rrup .lt. 450. ) then
        A_phi2 = d3_phi2 + d4_phi2 * (Rrup-225)/225. + d5_phi2 * ((Rrup-225)/225.)**2
      else
        A_phi2 = d3_phi2 + d4_phi2 + d5_phi2
      endif

c     f2 term for phi2 model (eq 5.3)
      if ( specT .lt. T1_phi2 ) then
        f2 = 1 - alpha_phi2
      elseif ( specT .lt. T2_phi2 ) then
        f2 = 1 - alpha_phi2 *alog(specT/T2_phi2)/alog(T1_phi2/T2_phi2)
      elseif ( specT .lt. T3_phi2 ) then
        f2 = 1.0
      elseif ( specT .lt. T4_phi2 ) then
        f2 = alog(specT/T4_phi2)/alog(T3_phi2/T4_phi2)
      else
        f2 = 0.
      endif
      f2_PGA = 1 - alpha_phi2

c     Compute added variance term, phi2
c      phi2_SQ_add = A_phi2 * f2

c     f2 for phi3 (same for as eq 5.3, but different coeff)
      if ( specT .lt. T1_phi3 ) then
        f3 = 1 - alpha_phi3
      elseif ( specT .lt. T2_phi3 ) then
        f3 = 1 - alpha_phi3 *alog(specT/T2_phi3)/alog(T1_phi3/T2_phi3)
      elseif ( specT .lt. T3_phi3 ) then
        f3 = 1.0
      elseif ( specT .lt. T4_phi3 ) then
        f3 = alog(specT/T4_phi3)/alog(T3_phi3/T4_phi3)
      else
        f3 = 0.
      endif
      A_phi3 = d3_phi3
      f3_PGA = 1 - alpha_phi3

c      phi3_SQ_add = A_phi3 * f3

c     compute region-specific phi
      if ( region .eq. 3) then
        phi_LIN = sqrt( phi1_SQ + A_phi3*f3 )
        phi_LIN_PGA = sqrt( phi1_SQ_PGA + A_phi3*f3_PGA )
        phi_S2S = sqrt( phi1_sq_100 + A_phi3*f3 - 0.165)
        phi_S2S_PGA = sqrt( phi1_SQ_PGA_100 + A_phi3*f3_PGA - 0.165)
      elseif ( region .eq. 4 .or. region .eq. 6 ) then
        phi_LIN = sqrt( phi1_SQ + A_phi2*f2 + A_phi3*f3)
        phi_LIN_PGA = sqrt( phi1_SQ_PGA + A_phi2*f2_PGA + A_phi3*f3_PGA)
        phi_S2S = sqrt( phi1_SQ_100 + A_phi2_100*f2 + A_phi3*f3 - 0.18)
        phi_S2S_PGA = sqrt( phi1_SQ_PGA_100 + A_phi2_100*f2 + A_phi3*f3_PGA - 0.18)
      else
        phi_LIN = sqrt(phi1_SQ)
        phi_LIN_PGA = sqrt(phi1_SQ)
        phi_S2S = sqrt( phi1_SQ_100 - 0.165)
        phi_S2S_PGA = sqrt( phi1_SQ_PGA_100 - 0.165)
      endif

c     compute linear single-station sigma
      PhiSS_LIN = sqrt( phi_LIN**2 - phi_S2S**2 )
      PhiSS_LIN_PGA = sqrt( phi_LIN_PGA**2 - phi_S2S**2 )

c     NL site effects on phi and tau

c     eq 5.9
      If ( vs30 .ge. VLIN_T ) then
        partial_f_PGA = 0.
      else
        partial_f_PGA = b_soil_T * rockPGA *
     1       ( -1./(rockPGA+c) + 1. / (rockPGA+c*(vs30/VLIN_T)**n) )
      endif

c     eq 5.7
      phi_AMP = 0.3
      phi_B_PGA = sqrt( phi_LIN_PGA**2 - phi_AMP**2)
      phi_B = sqrt( phi_LIN**2 - phi_AMP**2)

c     eq 5.8
      phiSQ_NL = phi_LIN**2 + partial_f_PGA**2 * phi_B**2
     1           + 2. * partial_f_PGA * phi_B_PGA * phi_B * rhoWT

c     eq 5.10
      tauSQ_NL = tau_LIN**2 + partial_f_PGA**2 * tau_LIN_PGA**2
     1           + 2. * partial_f_PGA * tau_LIN_PGA * tau_LIN * rhoBT

      phi = sqrt(phiSQ_NL)
      tau = sqrt(tauSQ_NL)
      sigma = sqrt (phiSQ_NL + tauSQ_NL)

c     NL effects on phiSS and sigmaSS
c     eq 5.7
      phi_AMP = 0.3
      phiSS_B_PGA = sqrt( phiSS_LIN_PGA**2 - phi_AMP**2)
      phiSS_B = sqrt( phiSS_LIN**2 - phi_AMP**2)

c     eq 5.8
      phiSS_SQ_NL = phiSS_LIN**2 + partial_f_PGA**2 * phiSS_B**2
     1           + 2. * partial_f_PGA * phiSS_B_PGA * phiSS_B * rhoWT

      phiSS = sqrt(phiSS_SQ_NL)
      sigmaSS = sqrt (phiSS_SQ_NL + tauSQ_NL)


c     Convert units spectral acceleration in gal
      mu = mu + 6.89

      return
      end

c -------------------------------------------------------------------

      subroutine S35_z25_interp ( vs30, x1, x2, y1, y2, y)
      real vs30, x1, x2, y1, y2, y, x, x1Log, x2Log

      x = alog(vs30)
      x1Log = alog(x1)
      x2Log = alog(x2)
      if ( x .lt. x1Log ) then
        y = y1
      elseif ( x .gt. x2Log ) then
        y = y2
      else
        y = (x-x1Log)/(x2Log-x1Log) * (y2-y1) + y1
      endif
      y = exp(y)

      return
      end

c ----------------------------------------------------------------------
      subroutine S35_KBCG2019 ( mag, Ftype, rRup, vs30, z25, lnSa, sigma, phi, tau,
     2                     specT, period2, iflag, depth, disthypo, iRegion, mbInter, mbSlab, ztor, CasBas, Z10 )

      implicit none
      integer MAXPER
      parameter (MAXPER=23)
      integer count1, count2, iflag, nPer, i1, i, iRegion, CasBas

      real Period(MAXPER), sigphi(MAXPER), sigtau(MAXPER), theta1a_Global(MAXPER), theta2(MAXPER),
     1     theta2a(MAXPER), theta3(MAXPER), theta4(MAXPER), theta4a(MAXPER), theta5(MAXPER),
     2     theta9(MAXPER), theta9a(MAXPER), delta_zb_if(MAXPER), delta_zb_slab(MAXPER),
     3     nft_1(MAXPER), nft_2(MAXPER), k1(MAXPER), k2(MAXPER), theta1_reg_AL(MAXPER),
     4     theta1_reg_CAS(MAXPER), theta1_reg_CAM(MAXPER), theta1_reg_Ja(MAXPER), theta1_reg_NZ(MAXPER),
     5     theta1_reg_SA(MAXPER), theta1_reg_TW(MAXPER), theta1_global(MAXPER), theta7_reg_AL(MAXPER),
     6     theta7_reg_Cas(MAXPER), theta7_reg_CAM(MAXPER), theta7_reg_Ja(MAXPER), theta7_reg_NZ(MAXPER),
     7     theta7_reg_SA(MAXPER), theta7_reg_Tw(MAXPER), theta7_global(MAXPER), theta62_reg_AL(MAXPER),
     8     theta62_reg_Cas(MAXPER), theta62_reg_CAM(MAXPER), theta62_reg_Ja(MAXPER), theta62_reg_NZ(MAXPER),
     9     theta62_reg_SA(MAXPER), theta62_reg_Tw(MAXPER), theta6_global(MAXPER)
      real sigphiT, sigtauT, theta2T, theta2aT, theta3T, theta4T, theta4aT, theta5T, theta9T,
     1     theta9aT, delta_zb_ifT, delta_zb_slabT, nft_1T, nft_2T, k1T, k2T, theta1_reg_AlT,
     2     theta1_reg_casT, theta1_reg_CAMT, theta1_reg_JaT, theta1_reg_NZT, theta1_reg_SAT,
     3     theta1_reg_TwT, theta1_globalT, theta7_reg_AlT, theta7_reg_CasT, theta7_reg_CAMT, theta7_reg_JaT,
     4     theta7_reg_NZT, theta7_reg_SAT, theta7_reg_TwT, theta7_globalT, theta62_reg_AlT,
     5     theta62_reg_CasT, theta62_reg_CAMT, theta62_reg_JaT, theta62_reg_NZT, theta62_reg_SAT,
     6     theta62_reg_TwT, theta6_globalT
      real theta1a_reg_AL(MAXPER),theta1a_reg_CAS(MAXPER), theta1a_reg_CAM(MAXPER), theta1a_reg_Ja(MAXPER),
     2     theta1a_reg_NZ(MAXPER), theta1a_reg_SA(MAXPER), theta1a_reg_TW(MAXPER)
      real theta1a_globalT, theta1a_reg_ALT, theta1a_reg_CAST, theta1a_reg_CAMT, theta1a_reg_jaT,
     1     theta1a_reg_NZT, theta1a_reg_SAT, theta1a_reg_TWT
      real CasBasinInter(MAXPER), CasBasinSlope(MAXPER), JapBasinInter(MAXPER), JapBasinSlope(MAXPER)
      real CasBasinInterT, CasBasinSlopeT, JapBasinInterT, JapBasinSlopeT, deltaZ25, predZ25, Z10, deltaZ10, predZ10
      real TaiBasinInter(MAXPER), TaiBasinSlope(MAXPER), TaiBasinInterT, TaiBasinSlopeT
      real NZBasinInter(MAXPER), NZBasinSlope(MAXPER), NZBasinInterT, NZBasinSlopeT
      real dmbInter(MAXPER), mbInterPGA, dmbInterT, mbInter0, CasSeaResid(MAXPER), CasSeaResidT

      real c, n, theta10, deltam, deltaz, zbif, zbslab, mref, zifref, zslabref, delta_mb_if, delta_mb_slab
      real delta_mb_BCH
      real theta1RegT, theta6RegT, theta7RegT
      real h, hpga, ztor
      real sigma, lnSa, pgaRock, vs30, rRup, disthypo, mag, depth, Ftype
      real period2, specT, z25, phi, tau, mbinter, mbslab
      real fconst, fmag, fgeom, fdepth, fatten, fsite
      real fconstPga, fmagPGa, fgeomPga, fdepthPGA, fattenPGA, fsitePGA
      real theta6RegPGA, theta7RegPGA, theta1RegPGA, fbasin
      real fsiteRockPGA, fbasinPGA, LnPGA


      data Period /  0.0, -1.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75,
     1               1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0  /
      data k1  / 865.0, 400.0, 865.0, 865.0, 908.0, 1054.0, 1086.0, 1032.0, 878.0, 748.0, 654.0, 587.0,
     1           503.0, 457.0, 410.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0 /
      data k2  / -1.186, -1.955, -1.186, -1.219, -1.273, -1.346, -1.471, -1.624, -1.931, -2.188, -2.381,
     1           -2.518, -2.657, -2.669, -2.401, -1.955, -1.025, -0.299, 0.0, 0.0, 0.0, 0.0, 0.0 /

C     Coefficients 9/02/2021 - Updated Alaska Coefficients Smoothed Long Period
      data theta1_global / 3.715306652, 6.43354268, 3.575653364, 3.716488614, 4.019203939, 4.544768777,
     1                     4.920943189, 5.073489777, 5.025990527, 4.769660602, 4.441347564,
     2                     4.095960014, 3.433449152, 2.846435008, 1.701639568, 0.890491783,
     3                    -0.175096723, -0.854226686, -1.707596773, -2.257679636, -2.663615911,
     4                    -3.368096986, -3.837778998 /
      data theta1a_global / 4.789080359, 9.43573316, 4.559571984, 4.808206646, 5.167197596, 5.640302621,
     1                      5.848237896, 5.830854259, 5.531066721, 5.132677534, 4.737976746,
     2                      4.374731034, 3.757219946, 3.265243755, 2.398836593, 1.829616288,
     3                      1.090154171, 0.586593296, -0.146183991, -0.711213054, -1.182577041,
     4                     -2.100674191, -2.768999266 /
      data theta2 / -2.460896366, -2.217436, -2.428916739, -2.459748358, -2.513834312, -2.589450354,
     1              -2.619896357, -2.608322791, -2.531737908, -2.434298892, -2.338216333,
     2              -2.250204787, -2.102934572, -1.989991219, -1.811384538, -1.719992958,
     3              -1.652688213, -1.646034912, -1.676629343, -1.709830071, -1.732636218,
     4              -1.750019169, -1.735295193 /
      data theta2a / -2.438099383, -2.6695061, -2.395028055, -2.443224205, -2.496444857, -2.541066633,
     1               -2.525555518, -2.477878863, -2.361283158, -2.252963091, -2.16215032,
     2               -2.088032245, -1.979317168, -1.907871091, -1.818965748, -1.791480543,
     3               -1.79780401, -1.825030176, -1.875875935, -1.908284477, -1.925829475,
     4               -1.933397388, -1.918143528 /
      data theta3 / 0.10393058, 0.1278845, 0.104318059, 0.104711329, 0.104531274, 0.10337605, 0.101545988,
     1              0.099788169, 0.096885665, 0.094745404, 0.093177224, 0.092031873, 0.090583773,
     2              0.089848329, 0.089517549, 0.090095535, 0.091896751, 0.09368432, 0.096501829,
     3              0.098438255, 0.099773762, 0.101593855, 0.102312143 /
      data theta4 / 0.952187401, 1.01448335, 0.945373983, 0.958823587, 0.973096474, 0.991235338, 1.000331099,
     1              1.001912667, 0.997173516, 0.990941872, 0.986668856, 0.985056106, 0.989150137,
     2              1.000679554, 1.047320756, 1.103813757, 1.215946727, 1.313705228, 1.463119514,
     3              1.566045928, 1.637882578, 1.738643843, 1.780438096 /
      data theta4a / 1.10501155, 1.30757813, 1.108606759, 1.116235958, 1.109268454, 1.090513308, 1.078532394,
     1               1.078729359, 1.099427554, 1.13154674, 1.166600358, 1.201257777, 1.26504902,
     2               1.320126168, 1.426015163, 1.500926099, 1.60075645, 1.666319535, 1.752683857,
     3               1.811761179, 1.857447102, 1.941357713, 2.00091285 /
      data theta5 / 0.128150984, 0.123716, 0.129146009, 0.132406953, 0.134842835, 0.13730139, 0.137915277,
     1              0.137330854, 0.134973442, 0.132355343, 0.12994531, 0.12784598, 0.124505287,
     2              0.122092052, 0.118613921, 0.117157008, 0.116773256, 0.117650902, 0.120191283,
     3              0.122611682, 0.124654241, 0.128327776, 0.130642589 /
      data theta6_global / -0.002634759, -0.0013477, -0.002699424, -0.002614211, -0.002548149, -0.002550059,
     1                     -0.002672773, -0.002820768, -0.003084744, -0.003273108, -0.003413419,
     2                     -0.003508116, -0.003620243, -0.003627501, -0.003470123, -0.003218289,
     3                     -0.002747866, -0.002339633, -0.001756502, -0.001376542, -0.001123589,
     4                     -0.000826115, -0.000728431 /
      data theta7_global / 0.88767875, 1.687675, 0.892437982, 0.930261469, 1.01275289, 1.1880131, 1.372273896,
     1                     1.531423527, 1.800940611, 2.015482111, 2.182527609, 2.307342334, 2.444039384,
     2                     2.451673734, 2.096480325, 1.516381724, 0.459844405, -0.176024266, -0.592415429,
     3                    -0.585389016, -0.512670766, -0.431586954, -0.428490629 /
      data theta9 / 0.025480655, 0.0161514, 0.025058198, 0.027178223, 0.029343849, 0.031911432, 0.032822575,
     1              0.032447898, 0.030225714, 0.027520137, 0.024917221, 0.02257519, 0.01872305,
     2              0.015822693, 0.011351668, 0.009168065, 0.00781531, 0.008069539, 0.009820042,
     3              0.011723614, 0.013386835, 0.016398519, 0.01828198 /
      data theta9a / 0.022385169, 0.01552129, 0.021963597, 0.023277106, 0.02453943, 0.02575995, 0.02579158,
     1               0.025080709, 0.022995598, 0.02084772, 0.018906568, 0.017206788, 0.014446286,
     2               0.012344736, 0.008820911, 0.006644533, 0.004020303, 0.002398108, 0.000278753,
     3              -0.001183113, -0.002310934, -0.004286542, -0.005541411 /
c      data theta1_reg_Al / 3.324582852, 6.192513888, 3.111678527, 3.272985702, 3.599643065, 4.153519077,
c     1                     4.558982077, 4.744258602, 4.768031065, 4.579911977, 4.311011865, 4.015117027,
c     2                     3.424768965, 2.881877487, 1.77365556, 0.947723395, -0.188551198, -0.934816148,
c     3                    -1.861303948, -2.41549166, -2.78204016, -3.299650373, -3.545156523 /
      data theta1_reg_Al / 3.324582852, 6.192513888, 3.111678527, 3.272985702, 3.599643065, 4.153519077,
     1                     4.558982077, 4.744258602, 4.768031065, 4.579911977, 4.311011865, 4.015117027,
     2                     3.424768965, 2.881877487, 1.77365556, 0.947723395, -0.188551198, -0.934816148,
     3                    -1.861303948, -2.41549166, -2.78204016, -3.376621108, -3.812917708 /
      data theta1_reg_Cas / 3.600304874, 6.43341768, 3.480334924, 3.601501862, 3.881527524, 4.387812562,
     1                      4.768796999, 4.937742037, 4.927519437, 4.702098799, 4.396801049,
     2                      4.068026237, 3.425683137, 2.848400136, 1.70893745, 0.894849762,
     3                     -0.177633126, -0.859883288, -1.712325988, -2.258375988, -2.660212813,
     4                     -3.357720901, -3.824347713 /
      data theta1_reg_CAM / 3.696420974, 6.41128435, 3.521664099, 3.684086811, 4.009646711, 4.549220574,
     1                      4.916524961, 5.053960299, 4.981228686, 4.710481311, 4.375634036,
     2                      4.028587899, 3.370173599, 2.791079115, 1.666433496, 0.869768161,
     3                     -0.180839526, -0.854612751, -1.706495114, -2.257421039, -2.663634851,
     4                     -3.365218976, -3.829128839 /
      data theta1_reg_Ja / 4.008987257, 6.50527566, 3.786327782, 3.966415432, 4.346256007, 4.964188257,
     1                     5.36193617, 5.488463607, 5.350477145, 5.005995045, 4.606384295,
     2                     4.206322057, 3.472405045, 2.847101424, 1.673874536, 0.868523382,
     3                    -0.172629318, -0.83651918, -1.684797118, -2.243241693, -2.659023743,
     4                    -3.375555293, -3.83850338 /
      data theta1_reg_NZ / 3.943337888, 6.46045166, 3.810325438, 3.912658438, 4.214684, 4.772026638,
     1                     5.180892763, 5.346631513, 5.2890084, 5.002106925, 4.6387195,
     2                     4.259727838, 3.5410181, 2.9129436, 1.712003021, 0.880914761,
     3                    -0.186117187, -0.8541594, -1.68967475, -2.233627612, -2.640906075,
     4                    -3.360428112, -3.845874587 /
      data theta1_reg_SA / 4.051210372, 6.43999505, 3.935535147, 4.057270272, 4.354312135, 4.88884236,
     1                     5.274643422, 5.425392947, 5.350816772, 5.05429396, 4.68392416,
     2                     4.29905291, 3.570140297, 2.93325554, 1.71546997, 0.873661997,
     3                    -0.203168403, -0.87287104, -1.703078565, -2.239111678, -2.638954865,
     4                    -3.344667653, -3.821734003 /
      data theta1_reg_Tw / 3.446871055, 6.43448706, 3.31199428, 3.43677063, 3.698230605, 4.16120153,
     1                     4.51271993, 4.676007205, 4.68984193, 4.504351367, 4.23980503,
     2                     3.947873555, 3.364050705, 2.827443326, 1.739913367, 0.94144859,
     3                    -0.135476695, -0.83215377, -1.70617387, -2.261116483, -2.664491845,
     4                    -3.355595383, -3.816492708 /
c      data theta1a_reg_Al / 4.149829797, 9.8146774, 3.923380172, 4.24921056, 4.60123381, 4.99344756,
c     1                      5.120876872, 5.069682685, 4.792279285, 4.470693522, 4.168180172, 3.896536697,
c     2                      3.439476451, 3.071433322, 2.388354905, 1.892909355, 1.167362272, 0.626671222,
c     3                     -0.164460478, -0.724734053, -1.139277078, -1.791549815, -2.13024199 /
      data theta1a_reg_Al / 4.149829797, 9.8146774, 3.923380172, 4.24921056, 4.60123381, 4.99344756,
     1                      5.120876872, 5.069682685, 4.792279285, 4.470693522, 4.168180172, 3.896536697,
     2                      3.439476451, 3.071433322, 2.388354905, 1.892909355, 1.167362272, 0.626671222,
     3                     -0.164460478, -0.724734053, -1.139277078, -2.098557165, -2.735575178 /
      data theta1a_reg_Cas / 4.486160214, 9.44157083, 4.309173564, 4.492272227, 4.780066177, 5.191630102,
     1                       5.409688439, 5.437210652, 5.242301202, 4.932354902, 4.604524402,
     2                       4.290405089, 3.734440359, 3.274077642, 2.429917247, 1.855767138,
     3                       1.097427152, 0.582327377, -0.154870361, -0.714677586, -1.179599661,
     4                      -2.089046048, -2.759944648 /
      data theta1a_reg_CAM / 5.041801411, 9.40842905, 4.809797423, 5.138580211, 5.547287961, 6.041135236,
     1                       6.214868648, 6.146773936, 5.751216648, 5.279849461, 4.832425536,
     2                       4.431417661, 3.767458256, 3.251794002, 2.368038756, 1.801356947,
     3                       1.074441061, 0.579301661, -0.148164002, -0.713874664, -1.186648539,
     4                      -2.103672102, -2.764279914 /
      data theta1a_reg_Ja / 5.36094074, 9.56171671, 5.048204203, 5.35605789, 5.823778528, 6.424629515,
     1                      6.65407559, 6.587999378, 6.135932665, 5.591939103, 5.07930914,
     2                      4.624838865, 3.885418432, 3.324406957, 2.394374697, 1.819634315,
     3                      1.09961974, 0.610200003, -0.12041156, -0.696791947, -1.179522347,
     4                     -2.106138247, -2.756891322 /
      data theta1a_reg_NZ / 4.822427116, 9.44242678, 4.653637454, 4.860414466, 5.176631041, 5.618999054,
     1                      5.835589729, 5.839812216, 5.579030904, 5.203979691, 4.820203141,
     2                      4.459910304, 3.835826743, 3.330662296, 2.430571835, 1.837887843,
     3                      1.077221991, 0.570044391, -0.154812746, -0.710636734, -1.176416246,
     4                     -2.094423559, -2.773091609 /
      data theta1a_reg_SA / 5.079648822, 9.36631119, 4.840730747, 5.154175935, 5.57760931, 6.114167272,
     1                      6.32265626, 6.26824361, 5.86644186, 5.371816747, 4.898129985,
     2                      4.472628647, 3.769547817, 3.22733893, 2.313242356, 1.74213861,
     3                      1.029526247, 0.55235046, -0.15173344, -0.708760028, -1.181004115,
     4                     -2.10961869, -2.784359065 /
      data theta1a_reg_Tw / 4.391529844, 9.33738787, 4.203531294, 4.388374744, 4.666983819, 5.054938069,
     1                      5.254962794, 5.277204731, 5.092471831, 4.803260306, 4.497810356,
     2                      4.204567281, 3.682760571, 3.246846383, 2.434937705, 1.870986279,
     3                      1.110181869, 0.585690969, -0.166118681, -0.732027631, -1.197657881,
     4                     -2.099842756, -2.761835231 /
c      data theta7_reg_AL / 0.551140468, 1.246646108, 0.572924573, 0.606050262, 0.715872694, 0.986759236,
c     1                     1.244194485, 1.399564893, 1.549306356, 1.62948428, 1.699971678, 1.76786809,
c     2                     1.871947148, 1.903775085, 1.675204801, 1.214579429, 0.299504686, -0.29013837,
c     3                    -0.718257827, -0.739439677, -0.682618715, -0.623527877, -0.66026114 /
      data theta7_reg_AL / 0.551140468, 1.246646108, 0.572924573, 0.606050262, 0.715872694, 0.986759236,
     1                     1.244194485, 1.399564893, 1.549306356, 1.62948428, 1.699971678, 1.76786809,
     2                     1.871947148, 1.903775085, 1.675204801, 1.214579429, 0.299504686, -0.29013837,
     3                    -0.718257827, -0.739439677, -0.682618715, -0.631117395, -0.569526682 /
      data theta7_reg_Cas / 0.662527069, 1.55887653, 0.667082832, 0.707218024, 0.782766296, 0.959353513,
     1                      1.183366699, 1.382117651, 1.682451037, 1.889590743, 2.041945092,
     2                      2.156154854, 2.287970004, 2.302687729, 1.969235669, 1.400361498,
     3                      0.356549952, -0.262818326, -0.638106383, -0.592149796, -0.487937496,
     4                     -0.364708696, -0.362123483 /
      data theta7_reg_CAM / 1.147446734, 1.85177447, 1.109680891, 1.161897312, 1.24741899, 1.368367591,
     1                      1.494059219, 1.643536028, 1.973240894, 2.261112938, 2.476115899,
     2                      2.621750749, 2.748278586, 2.718123161, 2.279570523, 1.664982299,
     3                      0.59274564, -0.063796126, -0.544892764, -0.580628139, -0.519162314,
     4                     -0.410008426, -0.387004114 /
      data theta7_reg_Ja / 0.86972574, 1.53952407, 0.870923483, 0.912943812, 1.011507346, 1.257663811,
     1                     1.49028266, 1.630928725, 1.786976198, 1.896749447, 1.995613121,
     2                     2.081211992, 2.186512067, 2.187920655, 1.83732347, 1.264853747,
     3                     0.224215841, -0.39993712, -0.802764383, -0.780815983, -0.68758182,
     4                    -0.549707683, -0.505965695 /
      data theta7_reg_NZ / 1.005782506, 1.78608459, 1.012431679, 1.056041081, 1.141592746, 1.296723954,
     1                     1.462358148, 1.625128475, 1.924311811, 2.15745282, 2.32477949,
     2                     2.438066215, 2.540737377, 2.519785177, 2.140886776, 1.572665626,
     3                     0.544469013, -0.08281136, -0.49804516, -0.483198185, -0.401492923,
     4                    -0.332415348, -0.387987785 /
      data theta7_reg_SA / 1.131304474, 2.098709, 1.147077322, 1.208208624, 1.276473281, 1.367020497,
     1                     1.469249031, 1.606935822, 1.945485117, 2.274055939, 2.542924494,
     2                     2.742512956, 2.955971819, 2.980301156, 2.590644506, 1.980834172,
     3                     0.913346769, 0.274633447, -0.174358569, -0.205313894, -0.159049369,
     4                    -0.113042244, -0.143450381 /
      data theta7_reg_Tw / 0.933688039, 1.67248556, 0.910859521, 0.953422087, 1.033921126, 1.187023224,
     1                     1.35787035, 1.520788001, 1.816593859, 2.057558467, 2.242567296,
     2                     2.376440146, 2.511422896, 2.502295571, 2.089379765, 1.460664727,
     3                     0.34695872, -0.315465991, -0.749482754, -0.741057329, -0.659579766,
     4                    -0.551928766, -0.531306454 /
c      data theta62_reg_AL / -0.002727249, -0.001656003, -0.002760841, -0.002737411, -0.002736235, -0.002810708,
c     1                      -0.002912096, -0.003009239, -0.00314302, -0.003235036, -0.003288512, -0.003332162,
c     2                      -0.003360836, -0.003320488, -0.003111385, -0.002887742, -0.002461284, -0.002101573,
c     3                      -0.001584346, -0.001266269, -0.001058123, -0.000816203, -0.000766118 /
      data theta62_reg_AL / -0.002727249, -0.001656003, -0.002760841, -0.002737411, -0.002736235, -0.002810708,
     1                      -0.002912096, -0.003009239, -0.00314302, -0.003235036, -0.003288512, -0.003332162,
     2                      -0.003360836, -0.003320488, -0.003111385, -0.002887742, -0.002461284, -0.002101573,
     3                      -0.001584346, -0.001266269, -0.001058123, -0.000868331, -0.00083004 /
      data theta62_reg_Cas / -0.00378872, -0.0013097, -0.003919403, -0.003787575, -0.003683177, -0.003680303,
     1                       -0.003883777, -0.004138913, -0.004589887, -0.004916118, -0.005140047,
     2                       -0.005276242, -0.005410199, -0.005393397, -0.005048778, -0.004564926,
     3                       -0.003666813, -0.002920969, -0.001909979, -0.001295891, -0.000913294,
     4                       -0.000494664, -0.000390228 /
      data theta62_reg_CAM / -0.001378838, -0.0009827, -0.00132529, -0.001329181, -0.00137577, -0.001486615,
     1                       -0.001608964, -0.00170333, -0.001825801, -0.001891814, -0.001937343,
     2                       -0.001959577, -0.001983679, -0.001975796, -0.001915739, -0.001837589,
     3                       -0.001701139, -0.001591565, -0.00140805, -0.00127133, -0.001167597,
     4                       -0.001003264, -0.000899137 /
      data theta62_reg_Ja / -0.00281163, -0.0019607, -0.002915696, -0.002888919, -0.002837462, -0.002794009,
     1                      -0.002831023, -0.00288636, -0.00304475, -0.003173064, -0.003288383,
     2                      -0.003360069, -0.0034541, -0.003481557, -0.003410605, -0.003229835,
     3                      -0.002861705, -0.0025148, -0.001984637, -0.001630897, -0.001413493,
     4                      -0.001218557, -0.00125276 /
      data theta62_reg_NZ / -0.001984352, -0.0005497, -0.00193905, -0.001895321, -0.001910774, -0.002021342,
     1                      -0.002193522, -0.002347337, -0.002575118, -0.002722229, -0.002818411,
     2                      -0.002881902, -0.002936009, -0.002924349, -0.002787277, -0.002607331,
     3                      -0.002268416, -0.001981616, -0.001552311, -0.001265052, -0.001061875,
     4                      -0.000774037, -0.000634906 /
      data theta62_reg_SA / -0.00209537, -0.0011235, -0.002071798, -0.002033545, -0.002045657, -0.002130904,
     1                      -0.002262679, -0.002366816, -0.002506983, -0.002564728, -0.002590677,
     2                      -0.002587561, -0.002556041, -0.002491535, -0.002310983, -0.002133191,
     3                      -0.001875352, -0.001684623, -0.001415576, -0.001219844, -0.001080541,
     4                      -0.000875424, -0.00073729 /
      data theta62_reg_Tw / -0.003046717, -0.0018548, -0.003128106, -0.002916352, -0.002687556, -0.002453776,
     1                      -0.002469507, -0.002641826, -0.003111805, -0.003563224, -0.003945293,
     2                      -0.004249917, -0.004668449, -0.004877669, -0.004934195, -0.00467427,
     3                      -0.003939296, -0.003234659, -0.00218645, -0.001522874, -0.001110153,
     4                      -0.000666067, -0.000611445 /
      data delta_zb_if / 16.81605756, 12.5660947, 16.8238131, 16.66218908, 17.18166098, 18.32074494, 18.92369876,
     1                   18.75864934, 17.286575, 15.32167741, 13.37493716, 11.61349725, 8.768199935,
     2                    6.72811905, 3.964759788, 2.99839585, 2.958589638, 3.467766463, 3.92726955,
     3                    3.370655838, 2.129626625, -2.2543041, -6.949083475 /
      data delta_zb_slab / -15.84125345, -13.359066, -16.12506073, -16.28403215, -16.31449001, -16.30017631,
     1                     -16.24007111, -16.14727661, -15.87766608, -15.52663778, -15.12641286,
     2                     -14.69940528, -13.81827288, -12.94797861, -10.95322734, -9.255852251,
     3                      -6.60057452, -4.646262708, -1.994229308, -0.305374608, 0.840704092,
     4                       2.472428092, 3.239206467 /
      data nft_1 / 0.894588882, 0.93637922, 0.896423521, 0.897168407, 0.896405171, 0.893480171, 0.889433915,
     1             0.885756638, 0.879903372, 0.875691026, 0.872635815, 0.870398829, 0.867543383,
     2             0.866017459, 0.864982593, 0.86559616, 0.868035809, 0.870533621, 0.874445533,
     3             0.877064259, 0.87880376, 0.881004891, 0.881692145 /
      data nft_2 / 0.198851305, 0.21094258, 0.199091802, 0.199150631, 0.198909735, 0.198209877, 0.197375325,
     1             0.19669797, 0.195765446, 0.195220279, 0.194916495, 0.194765678, 0.194726462,
     2             0.194876022, 0.195519746, 0.196226914, 0.197458581, 0.198405476, 0.199694436,
     3             0.200494411, 0.201025627, 0.201732235, 0.20200423 /
      data sigphi / 0.595755278, 0.51148649, 0.599595609, 0.59241506, 0.604430227, 0.638018891, 0.666119594,
     1              0.678368904, 0.677739961, 0.664654218, 0.649527824, 0.635787238, 0.615314866,
     2              0.603262958, 0.594773331, 0.59897004, 0.611546271, 0.617348448, 0.609290796,
     3              0.587788163, 0.562101633, 0.500922328, 0.453378144 /
      data sigtau / 0.488744744, 0.45098461, 0.494229731, 0.505458398, 0.511667679, 0.515792855, 0.514872686,
     1              0.511974691, 0.505515371, 0.500076345, 0.495846384, 0.492609188, 0.488214611,
     2              0.485606983, 0.482858615, 0.482329313, 0.482849552, 0.483535832, 0.48402177,
     3              0.483600245, 0.482729987, 0.480025032, 0.477483008 /

C     Cascadia Basin and Seattle Basin Coefficients updated 11/1/19
      data  CasBasinInter / -0.033754301, -0.012278, -0.03437367, -0.034662015, -0.032014189, -0.028009288,
     1                      -0.033581978, -0.025227652, -0.032857328, -0.034587858, -0.031240354,
     2                      -0.027222135, -0.045714985, -0.031599507, -0.032778512, -0.029057663,
     3                      -0.023875663, -0.039217868, 0.00042061, 0.018831799, 0.016382687,
     4                       0.012890241, 0.015757278 /
      data  CasBasinSlope / 0.020800512, 0.11844209, 0.020403813, 0.018233466, 0.014046186, -0.016393706,
     1                     -0.060272212, -0.077067179, -0.051386265, -0.016655504, 0.035523765,
     2                      0.051800012, 0.11065903, 0.129281574, 0.184907514, 0.173915665,
     3                      0.201715295, 0.230981824, 0.210505666, 0.244947111, 0.261637485,
     4                      0.220665408, 0.186336738 /
      data  CasSeaResid / -0.126368667, 0.118442089, -0.126722682, -0.128641453, -0.132860752, -0.17904093,
     1                    -0.261658424, -0.276774462, -0.24668039, -0.174197802, -0.08799574,
     2                    -0.046549387, 0.023017763, 0.09922914, 0.20957875, 0.236757351, 0.286400113,
     3                     0.351803208531182, 0.417364971913835, 0.448002389, 0.521325059, 0.485046352, 0.40273411 /
      data  JapBasinInter / -0.000552945, -0.0162883, 0.000252659, -0.000199039, 0.001563091, 0.009631227,
     1                       0.013762593, 0.010235154, -0.009497614, -0.013801794, -0.008879054,
     2                      -0.010823018, -0.01393108, -0.019981112, -0.0241875, -0.002670865,
     3                       0.004608325, -0.014419342, 0.017050267, 0.052368336, 0.010241569,
     4                      -0.005722435, -0.012161281 /
      data  JapBasinSlope / -0.027347239, 0.10123977, -0.028776069, -0.030456671, -0.03626794, -0.057726275,
     1                      -0.076256191, -0.068037275, -0.050343659, -0.035030582, -0.022592393,
     2                      -0.005371836, 0.030155088, 0.074733416, 0.141358586, 0.18875034,
     3                       0.249536336, 0.280580414, 0.305661189, 0.307163658, 0.292740348,
     4                       0.248739909, 0.214443383 /

      data  TaiBasinInter / 0.013176757, 0.0116407, 0.012098583, 0.011978376, 0.012910306, 0.010570091,
     1                      0.013508881, 0.01025206, 0.015651805, 0.014435965, 0.021983556,
     2                      0.019669811, 0.016888863, 0.018284057, 0.006109487, 0.015144663,
     3                      0.022357196, 0.013477586, 0.025343721, 0.030810534, 0.012982342,
     4                      0.013252225, -0.002628502 /
      data  TaiBasinSlope / 0.047719326, 0.078625, 0.046885023, 0.047268369, 0.047962287, 0.032748082,
     1                      0.040801587, 0.045775043, 0.040260644, 0.042212758, 0.050554267,
     2                      0.053314924, 0.051926286, 0.065895191, 0.085522692, 0.110850338,
     3                      0.15125173, 0.168365385, 0.157175193, 0.13208657, 0.135284696,
     4                      0.125746568, 0.112043062 /

      data NZBasinInter / -0.029314763, 0.00788926, -0.02830015, -0.028717664, -0.032177337, -0.044733843,
     1                    -0.043511236, -0.047968738, -0.039382031, -0.025000918, -0.022148713,
     2                    -0.010104446, 0.007677064, 0.016125544, 0.012059183, 0.030448367,
     3                     0.043018473, 0.046788465, 0.058672133, 0.042133993, 0.044754822,
     4                     0.03714183, 0.020092163 /
      data NZBasinSlope / -0.083452599, -0.0111422, -0.08319884, -0.082965508, -0.083556642, -0.081382979,
     1                    -0.097523701, -0.124030048, -0.118887899, -0.089642052, -0.11874965,
     2                    -0.084007842, -0.045937946, -0.019815379, 0.004010888, 0.044832102,
     3                     0.094098855, 0.112241422, 0.12463767, 0.136546104, 0.159273263,
     4                     0.126676539, 0.109094018 /

C     Mag scaling adjustment for Japan and South America Regions only
      data dmbInter / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1               -0.147628099, -0.252371901, -0.4, -0.4, -0.4, -0.4, -0.4 /

C     Constant parameters

      delta_mb_if = 0.0
      delta_mb_slab = 0.0
      delta_mb_BCH = 0.0
      c = 1.88
      n = 1.18
      theta10 = 0.0
      deltam = 0.1
      deltaz = 1.0
      zbif = 30.0
      zbslab = 80.0

      mref = 6.0
      Zifref = 15.0
      zslabref = 50.0
      mbInter0 = mbInter
      mbInterPGA = mbInter

C Find the requested spectral period and corresponding coefficients
      nPer = 23

C First check for the PGA case
      if (specT .eq. 0.0) then
         i1=1
         period2 = period(i1)
         sigphiT = sigphi(i1)
         sigtauT = sigtau(i1)
         theta2T = theta2(i1)
         theta2aT = theta2a(i1)
         theta3T = theta3(i1)
         theta4T = theta4(i1)
         theta4aT = theta4a(i1)
         theta5T = theta5(i1)
         theta9T = theta9(i1)
         theta9aT = theta9a(i1)
         delta_zb_ifT = delta_zb_if(i1)
         delta_zb_slabT = delta_zb_slab(i1)
         nft_1T = nft_1(i1)
         nft_2T = nft_2(i1)
         k1T = k1(i1)
         k2T = k2(i1)
         theta1_reg_AlT = theta1_reg_Al(i1)
         theta1_reg_CasT = theta1_reg_Cas(i1)
         theta1_reg_CAMT = theta1_reg_CAM(i1)
         theta1_reg_JaT = theta1_reg_Ja(i1)
         theta1_reg_NZT = theta1_reg_NZ(i1)
         theta1_reg_SAT = theta1_reg_SA(i1)
         theta1_reg_TwT = theta1_reg_Tw(i1)
         theta1_GlobalT = theta1_global(i1)

         theta1a_globalT = theta1a_global(i1)
         theta1a_reg_ALT = theta1a_reg_AL(i1)
         theta1a_reg_CasT = theta1a_reg_Cas(i1)
         theta1a_reg_CAMT = theta1a_reg_CAM(i1)
         theta1a_reg_JAT = theta1a_reg_Ja(i1)
         theta1a_reg_NZT = theta1a_reg_NZ(i1)
         theta1a_reg_SAT = theta1a_reg_SA(i1)
         theta1a_reg_TwT = theta1a_reg_Tw(i1)

         theta7_reg_AlT = theta7_reg_Al(i1)
         theta7_reg_CasT = theta7_reg_Cas(i1)
         theta7_reg_CAMT = theta7_reg_CAM(i1)
         theta7_reg_JaT = theta7_reg_Ja(i1)
         theta7_reg_NZT = theta7_reg_NZ(i1)
         theta7_reg_SAT = theta7_reg_SA(i1)
         theta7_reg_TwT = theta7_reg_Tw(i1)
         theta7_GlobalT = theta7_global(i1)

         theta62_reg_AlT = theta62_reg_Al(i1)
         theta62_reg_CasT = theta62_reg_Cas(i1)
         theta62_reg_CAMT = theta62_reg_CAM(i1)
         theta62_reg_JaT = theta62_reg_Ja(i1)
         theta62_reg_NZT = theta62_reg_NZ(i1)
         theta62_reg_SAT = theta62_reg_SA(i1)
         theta62_reg_TwT = theta62_reg_Tw(i1)
         theta6_GlobalT = theta6_global(i1)

         CasBasinInterT = CasBasinInter(i1)
         CasBasinSlopeT = CasBasinSlope(i1)
         CasSeaResidT = CasSeaResid(i1)
         JapBasinInterT = JapBasinInter(i1)
         JapBasinSlopeT = JapBasinSlope(i1)
         TaiBasinInterT = TaiBasinInter(i1)
         TaiBasinSlopeT = TaiBasinSlope(i1)
         NZBasinInterT = NZBasinInter(i1)
         NZBasinSlopeT = NZBasinSlope(i1)

         dmbInterT = 0.0

         goto 1011
      elseif (specT .eq. -1.0) then
         i1=2
         period2 = period(i1)
         sigphiT = sigphi(i1)
         sigtauT = sigtau(i1)
         theta2T = theta2(i1)
         theta2aT = theta2a(i1)
         theta3T = theta3(i1)
         theta4T = theta4(i1)
         theta4aT = theta4a(i1)
         theta5T = theta5(i1)
         theta9T = theta9(i1)
         theta9aT = theta9a(i1)
         delta_zb_ifT = delta_zb_if(i1)
         delta_zb_slabT = delta_zb_slab(i1)
         nft_1T = nft_1(i1)
         nft_2T = nft_2(i1)
         k1T = k1(i1)
         k2T = k2(i1)
         theta1_reg_AlT = theta1_reg_Al(i1)
         theta1_reg_CasT = theta1_reg_Cas(i1)
         theta1_reg_CAMT = theta1_reg_CAM(i1)
         theta1_reg_JaT = theta1_reg_Ja(i1)
         theta1_reg_NZT = theta1_reg_NZ(i1)
         theta1_reg_SAT = theta1_reg_SA(i1)
         theta1_reg_TwT = theta1_reg_Tw(i1)
         theta1_GlobalT = theta1_global(i1)

         theta1a_globalT = theta1a_global(i1)
         theta1a_reg_ALT = theta1a_reg_AL(i1)
         theta1a_reg_CasT = theta1a_reg_Cas(i1)
         theta1a_reg_CAMT = theta1a_reg_CAM(i1)
         theta1a_reg_JAT = theta1a_reg_Ja(i1)
         theta1a_reg_NZT = theta1a_reg_NZ(i1)
         theta1a_reg_SAT = theta1a_reg_SA(i1)
         theta1a_reg_TwT = theta1a_reg_Tw(i1)

         theta7_reg_AlT = theta7_reg_Al(i1)
         theta7_reg_CasT = theta7_reg_Cas(i1)
         theta7_reg_CAMT = theta7_reg_CAM(i1)
         theta7_reg_JaT = theta7_reg_Ja(i1)
         theta7_reg_NZT = theta7_reg_NZ(i1)
         theta7_reg_SAT = theta7_reg_SA(i1)
         theta7_reg_TwT = theta7_reg_Tw(i1)
         theta7_GlobalT = theta7_global(i1)

         theta62_reg_AlT = theta62_reg_Al(i1)
         theta62_reg_CasT = theta62_reg_Cas(i1)
         theta62_reg_CAMT = theta62_reg_CAM(i1)
         theta62_reg_JaT = theta62_reg_Ja(i1)
         theta62_reg_NZT = theta62_reg_NZ(i1)
         theta62_reg_SAT = theta62_reg_SA(i1)
         theta62_reg_TwT = theta62_reg_Tw(i1)
         theta6_GlobalT = theta6_global(i1)

         CasBasinInterT = CasBasinInter(i1)
         CasBasinSlopeT = CasBasinSlope(i1)
         CasSeaResidT = CasSeaResid(i1)
         JapBasinInterT = JapBasinInter(i1)
         JapBasinSlopeT = JapBasinSlope(i1)
         TaiBasinInterT = TaiBasinInter(i1)
         TaiBasinSlopeT = TaiBasinSlope(i1)
         NZBasinInterT = NZBasinInter(i1)
         NZBasinSlopeT = NZBasinSlope(i1)

         dmbInterT = 0.0

         goto 1011

      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=3,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020
         endif
      enddo
      write (*,'( i5,2f12.6)') nper, specT, period(nper)

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*)
      write (*,*) 'KBCG (2019) Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*)
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),sigphi(count1),sigphi(count2),
     +                   specT,sigphiT,iflag)
      call S24_interp (period(count1),period(count2),sigtau(count1),sigtau(count2),
     +                   specT,sigtauT,iflag)
      call S24_interp (period(count1),period(count2),theta2(count1),theta2(count2),
     +                   specT,theta2T,iflag)
      call S24_interp (period(count1),period(count2),theta2a(count1),theta2a(count2),
     +                   specT,theta2aT,iflag)
      call S24_interp (period(count1),period(count2),theta3(count1),theta3(count2),
     +                   specT,theta3T,iflag)
      call S24_interp (period(count1),period(count2),theta4(count1),theta4(count2),
     +                   specT,theta4T,iflag)
      call S24_interp (period(count1),period(count2),theta4a(count1),theta4a(count2),
     +                   specT,theta4aT,iflag)
      call S24_interp (period(count1),period(count2),theta5(count1),theta5(count2),
     +                   specT,theta5T,iflag)
      call S24_interp (period(count1),period(count2),theta9(count1),theta9(count2),
     +                   specT,theta9T,iflag)
      call S24_interp (period(count1),period(count2),theta9a(count1),theta9a(count2),
     +                   specT,theta9aT,iflag)
      call S24_interp (period(count1),period(count2),delta_zb_if(count1),delta_zb_if(count2),
     +                   specT,delta_zb_ifT,iflag)
      call S24_interp (period(count1),period(count2),delta_zb_slab(count1),delta_zb_slab(count2),
     +                   specT,delta_zb_slabT,iflag)
      call S24_interp (period(count1),period(count2),nft_1(count1),nft_1(count2),
     +                   specT,nft_1T,iflag)
      call S24_interp (period(count1),period(count2),nft_2(count1),nft_2(count2),
     +                   specT,nft_2T,iflag)
      call S24_interp (period(count1),period(count2),k1(count1),k1(count2),
     +                   specT,k1T,iflag)
      call S24_interp (period(count1),period(count2),k2(count1),k2(count2),
     +                   specT,k2T,iflag)
      call S24_interp (period(count1),period(count2),theta1_reg_AL(count1),theta1_reg_AL(count2),
     +                   specT,theta1_reg_ALT,iflag)
      call S24_interp (period(count1),period(count2),theta1_reg_CAS(count1),theta1_reg_CAS(count2),
     +                   specT,theta1_reg_CAST,iflag)
      call S24_interp (period(count1),period(count2),theta1_reg_CAM(count1),theta1_reg_CAM(count2),
     +                   specT,theta1_reg_CAMT,iflag)
      call S24_interp (period(count1),period(count2),theta1_reg_Ja(count1),theta1_reg_Ja(count2),
     +                   specT,theta1_reg_JaT,iflag)
      call S24_interp (period(count1),period(count2),theta1_reg_NZ(count1),theta1_reg_NZ(count2),
     +                   specT,theta1_reg_NZT,iflag)
      call S24_interp (period(count1),period(count2),theta1_reg_SA(count1),theta1_reg_SA(count2),
     +                   specT,theta1_reg_SAT,iflag)
      call S24_interp (period(count1),period(count2),theta1_reg_Tw(count1),theta1_reg_Tw(count2),
     +                   specT,theta1_reg_TwT,iflag)
      call S24_interp (period(count1),period(count2),theta1_global(count1),theta1_global(count2),
     +                   specT,theta1_globalT,iflag)
      call S24_interp (period(count1),period(count2),theta1a_global(count1),theta1a_global(count2),
     +                   specT,theta1a_GlobalT,iflag)
      call S24_interp (period(count1),period(count2),theta1a_reg_AL(count1),theta1a_reg_AL(count2),
     +                   specT,theta1a_reg_ALT,iflag)
      call S24_interp (period(count1),period(count2),theta1a_reg_Cas(count1),theta1a_reg_Cas(count2),
     +                   specT,theta1a_reg_CasT,iflag)
      call S24_interp (period(count1),period(count2),theta1a_reg_CAM(count1),theta1a_reg_CAM(count2),
     +                   specT,theta1a_reg_CAMT,iflag)
      call S24_interp (period(count1),period(count2),theta1a_reg_Ja(count1),theta1a_reg_Ja(count2),
     +                   specT,theta1a_reg_JaT,iflag)
      call S24_interp (period(count1),period(count2),theta1a_reg_NZ(count1),theta1a_reg_NZ(count2),
     +                   specT,theta1a_reg_NZT,iflag)
      call S24_interp (period(count1),period(count2),theta1a_reg_SA(count1),theta1a_reg_SA(count2),
     +                   specT,theta1a_reg_SAT,iflag)
      call S24_interp (period(count1),period(count2),theta1a_reg_Tw(count1),theta1a_reg_Tw(count2),
     +                   specT,theta1a_reg_TwT,iflag)
      call S24_interp (period(count1),period(count2),theta7_reg_AL(count1),theta7_reg_AL(count2),
     +                   specT,theta7_reg_ALT,iflag)
      call S24_interp (period(count1),period(count2),theta7_reg_CAS(count1),theta7_reg_CAS(count2),
     +                   specT,theta7_reg_CAST,iflag)
      call S24_interp (period(count1),period(count2),theta7_reg_CAM(count1),theta7_reg_CAM(count2),
     +                   specT,theta7_reg_CAMT,iflag)
      call S24_interp (period(count1),period(count2),theta7_reg_Ja(count1),theta7_reg_Ja(count2),
     +                   specT,theta7_reg_JaT,iflag)
      call S24_interp (period(count1),period(count2),theta7_reg_NZ(count1),theta7_reg_NZ(count2),
     +                   specT,theta7_reg_NZT,iflag)
      call S24_interp (period(count1),period(count2),theta7_reg_SA(count1),theta7_reg_SA(count2),
     +                   specT,theta7_reg_SAT,iflag)
      call S24_interp (period(count1),period(count2),theta7_reg_Tw(count1),theta7_reg_Tw(count2),
     +                   specT,theta7_reg_TwT,iflag)
      call S24_interp (period(count1),period(count2),theta7_global(count1),theta7_global(count2),
     +                   specT,theta7_globalT,iflag)
      call S24_interp (period(count1),period(count2),theta62_reg_AL(count1),theta62_reg_AL(count2),
     +                   specT,theta62_reg_ALT,iflag)
      call S24_interp (period(count1),period(count2),theta62_reg_CAS(count1),theta62_reg_CAS(count2),
     +                   specT,theta62_reg_CAST,iflag)
      call S24_interp (period(count1),period(count2),theta62_reg_CAM(count1),theta62_reg_CAM(count2),
     +                   specT,theta62_reg_CAMT,iflag)
      call S24_interp (period(count1),period(count2),theta62_reg_Ja(count1),theta62_reg_Ja(count2),
     +                   specT,theta62_reg_JaT,iflag)
      call S24_interp (period(count1),period(count2),theta62_reg_NZ(count1),theta62_reg_NZ(count2),
     +                   specT,theta62_reg_NZT,iflag)
      call S24_interp (period(count1),period(count2),theta62_reg_SA(count1),theta62_reg_SA(count2),
     +                   specT,theta62_reg_SAT,iflag)
      call S24_interp (period(count1),period(count2),theta62_reg_Tw(count1),theta62_reg_Tw(count2),
     +                   specT,theta62_reg_TwT,iflag)
      call S24_interp (period(count1),period(count2),theta6_global(count1),theta6_global(count2),
     +                   specT,theta6_globalT,iflag)
      call S24_interp (period(count1),period(count2),CasBasinInter(count1),CasBasinInter(count2),
     +                   specT,CasBasinInterT,iflag)
      call S24_interp (period(count1),period(count2),CasBasinSlope(count1),CasBasinSlope(count2),
     +                   specT,CasBasinSlopeT,iflag)
      call S24_interp (period(count1),period(count2),CasSeaResid(count1),CasSeaResid(count2),
     +                   specT,CasSeaResidT,iflag)
      call S24_interp (period(count1),period(count2),JapBasinInter(count1),JapBasinInter(count2),
     +                   specT,JapBasinInterT,iflag)
      call S24_interp (period(count1),period(count2),JapBasinSlope(count1),JapBasinSlope(count2),
     +                   specT,JapBasinSlopeT,iflag)
      call S24_interp (period(count1),period(count2),TaiBasinInter(count1),TaiBasinInter(count2),
     +                   specT,TaiBasinInterT,iflag)
      call S24_interp (period(count1),period(count2),TaiBasinSlope(count1),TaiBasinSlope(count2),
     +                   specT,TaiBasinSlopeT,iflag)
      call S24_interp (period(count1),period(count2),NZBasinInter(count1),NZBasinInter(count2),
     +                   specT,NZBasinInterT,iflag)
      call S24_interp (period(count1),period(count2),NZBasinSlope(count1),NZBasinSlope(count2),
     +                   specT,NZBasinSlopeT,iflag)
      call S24_interp (period(count1),period(count2),dmbInter(count1),dmbInter(count2),
     +                   specT,dmbInterT,iflag)

 1011 period2 = specT

C     Set the regional terms
c     0=Global
c     1=Alaska
c     2=Cascadia
c     3=Central America
c     4=Japan
c     5=New Zealand
c     6=South America
c     7=Taiwan

      if (iRegion .eq. 0) then
         if (ftype .eq. 0.0) then
            theta1RegT = theta1_globalT
            theta1RegPGA = theta1_global(1)
         else
            theta1RegT = theta1a_globalT
            theta1RegPGA = theta1a_global(1)
         endif
         theta6RegT = theta6_globalT
         theta7RegT = theta7_globalT
         theta6RegPGA = theta6_global(1)
         theta7RegPGA = theta7_global(1)
      elseif (iRegion .eq. 1) then
         if (ftype .eq. 0.0) then
            theta1RegT = theta1_reg_ALT
            theta1RegPGA = theta1_reg_AL(1)
         else
            theta1RegT = theta1a_reg_ALT
            theta1RegPGA = theta1a_reg_AL(1)
         endif
         theta6RegT = theta62_reg_ALT
         theta7RegT = theta7_reg_ALT
         theta6RegPGA = theta62_reg_AL(1)
         theta7RegPGA = theta7_reg_AL(1)
      elseif (iRegion .eq. 2) then
         if (ftype .eq. 0.0) then
            theta1RegT = theta1_reg_CasT
            theta1RegPGA = theta1_reg_Cas(1)
         else
            theta1RegT = theta1a_reg_CasT
            theta1RegPGA = theta1a_reg_Cas(1)
         endif
         theta6RegT = theta62_reg_CasT
         theta7RegT = theta7_reg_CasT
         theta6RegPGA = theta62_reg_Cas(1)
         theta7RegPGA = theta7_reg_Cas(1)
      elseif (iRegion .eq. 3) then
         if (ftype .eq. 0.0) then
            theta1RegT = theta1_reg_CAMT
            theta1RegPGA = theta1_reg_CAM(1)
         else
            theta1RegT = theta1a_reg_CAMT
            theta1RegPGA = theta1a_reg_CAM(1)
         endif
         theta6RegT = theta62_reg_CAMT
         theta7RegT = theta7_reg_CAMT
         theta6RegPGA = theta62_reg_CAM(1)
         theta7RegPGA = theta7_reg_CAM(1)
      elseif (iRegion .eq. 4) then
         if (ftype .eq. 0.0) then
            theta1RegT = theta1_reg_JaT
            theta1RegPGA = theta1_reg_Ja(1)
         else
            theta1RegT = theta1a_reg_JaT
            theta1RegPGA = theta1a_reg_Ja(1)
         endif
         theta6RegT = theta62_reg_JaT
         theta7RegT = theta7_reg_JaT
         theta6RegPGA = theta62_reg_Ja(1)
         theta7RegPGA = theta7_reg_Ja(1)
      elseif (iRegion .eq. 5) then
         if (ftype .eq. 0.0) then
            theta1RegT = theta1_reg_NZT
            theta1RegPGA = theta1_reg_NZ(1)
         else
            theta1RegT = theta1a_reg_NZT
            theta1RegPGA = theta1a_reg_NZ(1)
         endif
         theta6RegT = theta62_reg_NZT
         theta7RegT = theta7_reg_NZT
         theta6RegPGA = theta62_reg_NZ(1)
         theta7RegPGA = theta7_reg_NZ(1)
      elseif (iRegion .eq. 6) then
         if (ftype .eq. 0.0) then
            theta1RegT = theta1_reg_SAT
            theta1RegPGA = theta1_reg_SA(1)
         else
            theta1RegT = theta1a_reg_SAT
            theta1RegPGA = theta1a_reg_SA(1)
         endif
         theta6RegT = theta62_reg_SAT
         theta7RegT = theta7_reg_SAT
         theta6RegPGA = theta62_reg_SA(1)
         theta7RegPGA = theta7_reg_SA(1)
      elseif (iRegion .eq. 7) then
         if (ftype .eq. 0.0) then
            theta1RegT = theta1_reg_TwT
            theta1RegPGA = theta1_reg_Tw(1)
         else
            theta1RegT = theta1a_reg_TwT
            theta1RegPGA = theta1a_reg_Tw(1)
         endif
         theta6RegT = theta62_reg_TwT
         theta7RegT = theta7_reg_TwT
         theta6RegPGA = theta62_reg_Tw(1)
         theta7RegPGA = theta7_reg_Tw(1)
      endif

C     Adjust mbinterface term for Japan and South America
      if (iRegion .eq. 4 .or. iRegion .eq. 6) then
         if (ftype .eq. 0.0) then
            mbInter = mbInter0 + dmbInterT
         endif
      else
         if (ftype .eq. 0.0) then
            mbInter = mbInter0
         endif
      endif

C     Compute the base model including PGARock

      h = 10**(nft_1T+nft_2T*(mag-6.0))
      hpga = 10**(nft_1(1)+nft_2(1)*(mag-6.0))

C     Constant term
      fconst = theta1RegT
      fconstpga = theta1RegPGA

C     f(mag) term
      if (ftype .eq. 0.0) then
         call S35_KBCGLH (mag,mbInter,theta4T*(mbInter-Mref),theta4T,theta5T,deltam,fmag)
         call S35_KBCGLH (mag,mbInter,theta4(1)*(mbInter-Mref),theta4(1),theta5(1),deltam,fmagpga)
      else
         call S35_KBCGLH (mag,mbSlab,theta4aT*(mbslab-Mref),theta4aT,theta5T,deltam,fmag)
         call S35_KBCGLH (mag,mbSlab,theta4a(1)*(mbslab-Mref),theta4a(1),theta5(1),deltam,fmagpga)
      endif

C     f(geom) term
      if (ftype .eq. 0.0) then
         fgeom = (theta2T + theta3T*mag)*alog(Rrup + h)
         fgeompga = (theta2(1) + theta3(1)*mag)*alog(Rrup + hpga)
      else
         fgeom = (theta2aT + theta3T*mag)*alog(Rrup + h)
         fgeompga = (theta2a(1) + theta3(1)*mag)*alog(Rrup + hpga)
      endif

C     f(depth) term
      if (ftype .eq. 0.0) then
         call S35_KBCGLH (ztor,(zbif+delta_zb_ifT),theta9T*(zbif+delta_zb_ifT-Zifref),theta9T,theta10,deltaz,fdepth)
         call S35_KBCGLH (ztor,(zbif+delta_zb_ifT),theta9(1)*(zbif+delta_zb_ifT-Zifref),theta9(1),theta10,deltaz,fdepthpga)
      else
         call S35_KBCGLH (ztor,(zbslab+delta_zb_slabT),theta9aT*(zbslab+delta_zb_slabT-Zslabref),theta9aT,theta10,deltaz,fdepth)
         call S35_KBCGLH (ztor,(zbslab+delta_zb_slabT),theta9a(1)*(zbslab+delta_zb_slabT-Zslabref),theta9a(1),
     +                    theta10,deltaz,fdepthpga)
      endif

C     f(atten) term
C     Currently assumes that R2=Rrup (i.e., does not cross arc region)
      fatten = Rrup*theta6RegT
      fattenpga = Rrup*theta6RegPGA

C     f(site) term
C     For PGA Vs=1100m/s
      fsiteRockpga = (theta7RegPGA + k2(1)*n)*alog(1100.0/k1(1))

C     Compute PGARock
      pgarock = exp(fconstpga + fmagpga + fgeompga + fdepthpga + fattenpga + fsiteRockPga)

C     Now compute full site term
      if (vs30 .gt. k1T) then
         fsite = (theta7RegT + k2T*n)*alog(Vs30/k1T)
         fsitePGA = (theta7RegPGA + k2(1)*n)*alog(Vs30/k1(1))
      else
         fsite =theta7RegT*alog(vs30/k1T) + k2T*(alog(PGArock + c*(vs30/k1T)**n) - alog(PGARock + c))
         fsitePGA =theta7RegPGA*alog(vs30/k1(1)) + k2(1)*(alog(PGArock + c*(vs30/k1(1))**n) - alog(PGARock + c))
      endif

C     Compute the basin term - only for Cascadia, Japan, Taiwan, and New Zealand
      if (iRegion .eq. 2) then
         predZ25 = 8.2940 + (2.3026 - 8.2940) * Exp((alog(Vs30) - 6.3969)/0.2708)/(1 + Exp((alog(Vs30) - 6.3969)/0.2708))
         deltaZ25 = alog(Z25*1000.0) - predZ25
C     Add in Seattle Basin Residual if requested
         if (CasBas .eq. 1) then
            fbasin = CasSeaResidT
            fbasinPGA = CasSeaResid(1)
C     Add in Cascadia non-Seattle Basin Factor not to exceed Seattle Basin Factors
         elseif (CasBas .eq. 2) then
            if ((CasBasinInterT + CasBasinSlopeT * deltaz25) .gt. 0.0) then
               if ( (CasBasinInterT + CasBasinSlopeT * deltaz25) .gt. CasSeaResidT) then
                  fbasin = CasSeaResidT
                  fbasinPGA = CasSeaResid(1)
               else
                  fbasin = CasBasinInterT + CasBasinSlopeT * deltaz25
                  fbasinPGA = CasBasinInter(1) + CasBasinSlope(1) * deltaz25
               endif

           elseif ( (CasBasinInterT + CasBasinSlopeT * deltaz25) .le. 0.0) then
               if ( CasSeaResidT .gt. 0.0) then
                  fbasin = CasBasinInterT + CasBasinSlopeT * deltaz25
                  fbasinPGA = CasBasinInter(1) + CasBasinSlope(1) * deltaz25
               elseif ( (CasBasinInterT + CasBasinSlopeT * deltaz25) .lt. CasSeaResidT ) then
                  fbasin = CasSeaResidT
                  fbasinPGA = CasSeaResid(1)
               else
                  fbasin = CasBasinInterT + CasBasinSlopeT * deltaz25
                  fbasinPGA = CasBasinInter(1) + CasBasinSlope(1) * deltaz25
               endif
           endif
         endif

      elseif (iRegion .eq. 4) then
         predZ25 = 7.6894 + (2.3026 - 7.6894) * Exp((alog(Vs30) - 6.3092)/0.7529)/(1 + Exp((alog(Vs30) - 6.3092)/0.7529))
         deltaZ25 = alog(Z25*1000.0) - predZ25
         fbasin = JapBasinInterT + JapBasinSlopeT * deltaz25
         fbasinPGA = JapBasinInter(1) + JapBasinSlope(1) * deltaz25
      elseif (iRegion .eq. 7) then
         predZ10 = 6.3056 + (2.3026 - 6.3056) * Exp((alog(Vs30) - 6.1105)/0.4367)/(1 + Exp((alog(Vs30) - 6.1105)/0.4367))
         deltaZ10 = alog(Z10*1000.0) - predZ10
         fbasin = TaiBasinInterT + TaiBasinSlopeT * deltaz10
         fbasinPGA = TaiBasinInter(1) + TaiBasinSlope(1) * deltaz10
      elseif (iRegion .eq. 5) then
         predZ10 = 6.8598 + (2.3026 - 6.8598) * Exp((alog(Vs30) - 5.7457)/0.9156)/(1 + Exp((alog(Vs30) - 5.7457)/0.9156))
         deltaZ10 = alog(Z10*1000.0) - predZ10
         fbasin = NZBasinInterT + NZBasinSlopeT * deltaz10
         fbasinPGA = NZBasinInter(1) + NZBasinSlope(1) * deltaz10
      else
         fbasin = 0.0
         fbasinPGA = 0.0
      endif

C     Compute the site specific ground motion
      lnSA = fconst + fmag + fgeom + fdepth + fatten + fsite + fbasin
      lnPGA = fconstPGA + fmagPGA + fgeomPGA + fdepthPGA + fattenPGA + fsitePGA + fbasinPGA

C     Check that PSA > PGA for T < 0.1sec
C     In cases where PSA < PGA, then set PSA = PGA
      if (specT .lt. 0.1) then
         if (lnSA .lt. lnPGA ) then
            lnSA = lnPGA
         endif
      endif

c     Convert units spectral acceleration in gal
      lnSa = lnSa + 6.89

C     Set sigma values to return
      sigma = sqrt(sigPhiT*sigPhiT + sigTauT*sigTauT)
      phi = sigPhiT
      tau = sigTauT

      return
      end

c ----------------------------------------------------------------------
      subroutine S35_kbcgLH (x, x0, a, b0, b1, delta, LH)

C     Kuehn et al. (2019) LH Function

      real x, x0, a, b0, b1, delta, LH

C     For numerical stability for deep slab events check that exp((x-x0)/delta
C         is less than 85.0. If greater set equal to 85.0. These large
C         values do not change the ground motions since depth scaling is constant.
      if ((x - x0)/delta .gt. 85.0) then
         LH = a + b0*(85.0) + (b1 - b0)*delta*alog(1.0 + exp(85.0))
      else
         LH = a + b0*(x - x0) + (b1 - b0)*delta*alog(1.0 + exp((x - x0)/delta))
      endif

      return
      end

c ----------------------------------------------------------------------
      subroutine S35_PSHAB2019 ( mag, Ftype, rRup, vs30, z25, lnSa, sigma, phi, tau,
     2                     specT, period2, iflag, depth, disthypo, iRegion, mbInter, mbSlab, pnwbflag )

C     Model Version: December 16, 2020

      implicit none
      integer MAXPER
      parameter (MAXPER=26)
      integer count1, count2, iflag, nPer, i1, i, iRegion, pnwbflag

      real Period(MAXPER), c1(MAXPER), c1Slab(MAXPER), c4(MAXPER), c5(MAXPER), C6(MAXPER),
     1     f4(MAXPER), f5(MAXPER), Vc(MAXPER),
     2     Alaska_a0(MAXPER), CAM_a0(MAXPER),
     3     Japan_a0(MAXPER), SA_a0(MAXPER), Taiwan_a0(MAXPER), Global_a0(MAXPER),
     4     Global_c0(MAXPER), Alaska_c0(MAXPER), CAM_c0(MAXPER),
     5     Taiwan_C0(MAXPER), Cascadia_c0(MAXPER), Japan_Pac_c0(MAXPER), Japan_Phi_c0(MAXPER), SAN_c0(MAXPER), SAS_c0(MAXPER),
     6     Global_c0s(MAXPER), Alaska_c0s(MAXPER),  Aleutian_c0s(MAXPER), CAM_c0s(MAXPER), Japan_c0s(MAXPER),
     5     SAN_c0s(MAXPER), SAS_c0s(MAXPER), Taiwan_C0s(MAXPER), Cascadia_c0s(MAXPER),
     6     Global_c(MAXPER), Alaska_c(MAXPER), Cascadia_c(MAXPER), Japan_c(MAXPER),
     7     SA_c(MAXPER), Taiwan_c(MAXPER), d(MAXPER), m(MAXPER), db(MAXPER),
     8     Japan_e1(MAXPER), Japan_e2(MAXPER), Japan_e3(MAXPER),
     9     Cas_e1(MAXPER), Cas_e2(MAXPER), Cas_e3(MAXPER), SeaBasinD(MAXPER), BasinD(MAXPER)
      real Japan_c1(MAXPER), Taiwan_c1(MAXPER)
      real Aleutian_c0(MAXPER), Aleutian_c0T

      real c4slab(MAXPER), c5slab(MAXPER), c6slab(MAXPER), Alaska_a0Slab(MAXPER), Cascadia_a0Slab(MAXPER),
     1     CAM_a0Slab(MAXPER), Japan_a0Slab(MAXPER), SA_a0Slab(MAXPER), Taiwan_a0Slab(MAXPER), Global_a0Slab(MAXPER)

      real sigtau(MAXPER), sigphi1(MAXPER), sigphi2(MAXPER), sigphiv(MAXPER), tau, phi
      real c1T, c1SlabT, c4T, c5T, c6T, f4T, f5T, VcT, sigPhi1T, sigphi2T, sigPhivT, sigTauT, a0T, c0T, cT
      real Global_a0T, Alaska_a0T, CAM_a0T, Japan_a0T, SA_a0T, Taiwan_a0T
      real Global_c0T, Alaska_c0T, CAM_c0T, Taiwan_c0T, Cascadia_c0T, Japan_Pac_c0T, Japan_Phi_c0T, SAN_c0T, SAS_c0T
      real Global_c0sT, Alaska_c0sT, Aleutian_c0sT, CAM_c0sT, Japan_c0sT, SAN_c0sT, SAS_c0sT, Taiwan_c0sT, Cascadia_c0sT
      real Global_cT, Alaska_cT, Cascadia_cT, Japan_cT, SA_cT, Taiwan_cT, dT, mT, dbT
      real Japan_e1T, Japan_e2T, Japan_e3T, Cas_e1T, Cas_e2T, Cas_e3T
      real SeaBasinDT, BasinDT
      real Japan_c1T, Taiwan_c1T
      real c4slabT, c5SlabT, c6SlabT
      real Alaska_a0SlabT, Cascadia_a0SlabT, CAM_a0SlabT,
     1     Japan_a0SlabT, SA_a0SlabT, Taiwan_a0SlabT, Global_a0SlabT
      real a0Pga, c0Pga, fpPga, fmPga, fdPga, x

      real sigma, lnSa, pgaRock, vs30, rRup, disthypo, mag, depth, Ftype
      real period2, specT, z25
      real b4, f1, f3, Vb, Vref, mbInter, mbslab
      real Fm, Fp, Fd, Fs, Flin, Fnl, f2, db1, h, R, Rref, muz25, deltaz25, fb, cas_e3TBasin, cas_e2TBasin
      real phirup, phivs, phivar, V1

      data  period  / 0.0, -1.0, 0.01, 0.02, 0.025, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15,
     1                0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0  /

C     Model Constants - Updated (12/16/20) revised C6 for Interface

      Data Global_c0 / 4.082, 8.097, 3.714, 3.762, 3.859, 4.014, 4.223, 4.456, 4.742, 4.952, 5.08,
     1                 5.035, 4.859, 4.583, 4.18, 3.752, 3.085, 2.644, 2.046, 1.556, 1.167, 0.92,
     2                 0.595, 0.465, 0.078, 0.046 /
      Data Alaska_c0 / 4.458796298, 9.283796298, 4.094796298, 4.132796298, 4.246796298, 4.386796298,
     1                 4.553796298, 4.745796298, 4.972796298, 5.160796298, 5.285796298, 5.277796298,
     2                 5.154796298, 4.910796298, 4.548796298, 4.168796298, 3.510796298, 3.067796298,
     3                 2.513796298, 2.061796298, 1.709796298, 1.456796298, 1.207796298, 1.131796298,
     4                 0.758796298, 0.708796298 /
      Data Aleutian_c0 / 3.652796298, 8.374796298, 3.288796298, 3.338796298, 3.392796298, 3.535796298,
     1                   3.747796298, 3.959796298, 4.231796298, 4.471796298, 4.665796298, 4.661796298,
     2                   4.503796298, 4.276796298, 3.919796298, 3.486796298, 2.710796298, 2.238796298,
     3                   1.451796298, 0.906796298, 0.392796298, 0.099796298, -0.356203702, -0.601203702,
     4                  -1.137203702, -1.290203702 /
      Data Cascadia_c0 / 3.856, 7.728, 3.488, 3.536, 3.633, 3.788, 3.997, 4.23, 4.516, 4.726, 4.848,
     1                   4.798, 4.618, 4.34, 3.935, 3.505, 2.837, 2.396, 1.799, 1.31, 0.922, 0.675,
     2                   0.352, 0.223, -0.162, -0.193 /
      Data CAM_c0 / 2.875899908, 7.046899908, 2.564899908, 2.636899908, 2.731899908, 2.890899908,
     1              3.075899908, 3.287899908, 3.560899908, 3.788899908, 3.945899908, 3.943899908,
     2              3.800899908, 3.491899908, 3.128899908, 2.640899908, 1.987899908, 1.553899908,
     3              0.990899908, 0.534899908, 0.186899908, -0.087100092, -0.353100092, -0.491100092,
     4             -0.837100092, -0.864100092 /
      Data Japan_Pac_c0 / 5.373125851, 8.772125851, 5.022125851, 5.066125851, 5.140125851, 5.317125851,
     1                    5.564125851, 5.843125851, 6.146125851, 6.346125851, 6.425125851, 6.288125851,
     2                    5.972125851, 5.582125851, 5.091125851, 4.680125851, 3.906125851, 3.481125851,
     3                    2.870125851, 2.507125851, 2.160125851, 1.969125851, 1.675125851, 1.601125851,
     4                    1.270125851, 1.364125851 /
      Data Japan_Phi_c0 / 4.309125851, 7.579125851, 3.901125851, 3.935125851, 4.094125851, 4.278125851,
     1                    4.531125851, 4.816125851, 5.126125851, 5.333125851, 5.420125851, 5.289125851,
     2                    4.979125851, 4.592125851, 4.089125851, 3.571125851, 2.844125851, 2.371125851,
     3                    1.779125851, 1.293125851, 0.895125851, 0.607125851, 0.303125851, 0.183125851,
     4                   -0.143874149, -0.195874149 /
      Data SAN_c0 / 5.064671414, 8.528671414, 4.673671414, 4.694671414, 4.779671414, 4.935671414,
     1              5.182671414, 5.457671414, 5.788671414, 5.998671414, 6.103671414, 6.013671414,
     2              5.849671414, 5.603671414, 5.151671414, 4.719671414, 3.995671414, 3.512671414,
     3              2.875671414, 2.327671414, 1.950671414, 1.766671414, 1.524671414, 1.483671414,
     4              1.175671414, 1.271671414 /
      Data SAS_c0 / 5.198671414, 8.679671414, 4.807671414, 4.827671414, 4.911671414, 5.066671414,
     1              5.312671414, 5.586671414, 5.917671414, 6.126671414, 6.230671414, 6.140671414,
     2              5.974671414, 5.728671414, 5.277671414, 4.848671414, 4.129671414, 3.653671414,
     3              3.023671414, 2.481671414, 2.111671414, 1.932671414, 1.698671414, 1.665671414,
     4              1.366671414, 1.462671414 /
      Data Taiwan_c0 / 3.032846279, 7.559846279, 2.636846279, 2.698846279, 2.800846279, 2.926846279,
     1                 3.069846279, 3.236846279, 3.446846279, 3.643846279, 3.798846279, 3.827846279,
     2                 3.765846279, 3.602846279, 3.343846279, 3.028846279, 2.499846279, 2.140846279,
     3                 1.645846279, 1.217846279, 0.871846279, 0.596846279, 0.268846279, 0.014846279,
     4                -0.446153721, -0.473153721 /

      data Global_c0s / 9.907, 13.194, 9.962, 10.099, 10.181, 10.311, 10.588, 10.824, 11.084,
     1                 11.232, 11.311, 11.055, 10.803, 10.669, 10.116, 9.579, 8.837, 8.067,
     2                  6.829, 5.871, 5.2, 4.83, 4.173, 3.833, 3.132, 2.72 /
      data Cascadia_c0s / 9.6, 12.874, 9.802, 9.933, 10.009, 10.133, 10.404, 10.634, 10.888, 11.03,
     1                   11.103, 10.841, 10.583, 10.443, 9.884, 9.341, 8.593, 7.817, 6.573, 5.609,
     2                    4.932, 4.556, 3.893, 3.547, 2.84, 2.422 /
      data Alaska_c0s / 9.404, 12.79, 9.451, 9.587, 9.667, 9.808, 10.086, 10.379, 10.65,
     1                 10.816, 10.883, 10.633, 10.322, 10.116, 9.561, 8.973, 8.246, 7.507,
     2                  6.213, 5.206, 4.594, 4.206, 3.517, 3.142, 2.391, 2.031 /
      data Aleutian_c0s / 9.912, 13.6, 9.954, 10.086, 10.172, 10.302, 10.602, 10.862, 11.184,
     1                   11.304, 11.402, 11.183, 10.965, 10.87, 10.411, 9.901, 9.335, 8.68,
     2                    7.581, 6.671, 6.047, 5.667, 4.97, 4.592, 3.65, 2.95 /
      data CAM_c0s / 9.58, 12.81, 9.612, 9.771, 9.85, 9.993, 10.317, 10.563, 10.785, 10.841,
     1              10.809, 10.519, 10.268, 10.134, 9.598, 9.097, 8.324, 7.557, 6.35, 5.434,
     2               4.773, 4.441, 3.849, 3.502, 2.821, 2.408 /
      data Japan_c0s / 10.145, 13.248, 10.162, 10.306, 10.387, 10.498, 10.744, 10.981, 11.25,
     1                 11.466, 11.619, 11.351, 11.063, 10.878, 10.296, 9.711, 8.934, 8.164,
     2                  6.896, 5.935, 5.234, 4.849, 4.074, 3.814, 3.152, 2.791 /
      data SAN_c0s / 9.254, 12.754, 9.293, 9.403, 9.481, 9.592, 9.834, 10.027, 10.265, 10.467,
     1              10.566, 10.33, 10.124, 10.077, 9.539, 9.03, 8.258, 7.467, 6.22, 5.261, 4.567,
     2               4.176, 3.495, 3.038, 2.368, 1.939 /
      data SAS_c0s / 9.991, 12.927, 9.994, 10.152, 10.292, 10.459, 10.818, 11.102, 11.424, 11.49,
     1              11.32, 10.927, 10.555, 10.328, 9.639, 9.03, 8.258, 7.417, 6.18, 5.161, 4.517,
     2               4.076, 3.445, 3.038, 2.368, 1.939 /
      data Taiwan_c0s / 10.071, 13.516, 10.174, 10.273, 10.329, 10.451, 10.678, 10.86, 11.093,
     1                  11.283, 11.503, 11.32, 11.147, 11.079, 10.547, 10.049, 9.327, 8.504,
     2                   7.204, 6.227, 5.517, 5.157, 4.55, 4.229, 3.554, 3.166 /

C     Geometrical Spreading
      data c1 / -1.662, -1.661, -1.587, -1.593, -1.607, -1.63, -1.657, -1.687, -1.715,
     1          -1.737, -1.745, -1.732, -1.696, -1.643, -1.58, -1.519, -1.44, -1.419, -1.4,
     2          -1.391, -1.394, -1.416, -1.452, -1.504, -1.569, -1.676 /
      data c1slab / -2.543, -2.422, -2.554, -2.566, -2.578, -2.594, -2.629, -2.649, -2.65,
     1              -2.647, -2.634, -2.583, -2.539, -2.528, -2.452, -2.384, -2.338, -2.267,
     2              -2.166, -2.077, -2.015, -2.012, -1.989, -1.998, -2.019, -2.047 /

C     Anelastic Attenuation
      data Global_a0 / -0.00657, -0.00395, -0.00657, -0.00657, -0.00657, -0.00657, -0.00657,
     1                 -0.00657, -0.00657, -0.00657, -0.00657, -0.00657, -0.00657, -0.00657,
     2                 -0.00657, -0.00657, -0.00635, -0.0058, -0.00505, -0.00429, -0.00369,
     3                 -0.00321, -0.00244, -0.0016, -0.000766, 0.0 /
      data Alaska_a0 / -0.00541, -0.00404, -0.00541, -0.00541, -0.00541, -0.00541, -0.00541,
     1                 -0.00541, -0.00541, -0.00541, -0.00541, -0.00541, -0.00541, -0.00541,
     2                 -0.00541, -0.00541, -0.00478, -0.00415, -0.00342, -0.0029, -0.0025,
     3                 -0.00217, -0.00165, -0.00125, -0.000519, 0.0 /
      data CAM_a0 / -0.00387, -0.00153, -0.00387, -0.00387, -0.00387, -0.00387, -0.00387,
     1              -0.00387, -0.00387, -0.00387, -0.00387, -0.00387, -0.00387, -0.00387,
     2              -0.00387, -0.00387, -0.00342, -0.00297, -0.00245, -0.00208, -0.00179,
     3              -0.00156, -0.00118, -0.000895, -0.000371, 0.0 /
      data Japan_a0 / -0.00862, -0.00239, -0.00862, -0.00862, -0.00862, -0.00862, -0.00862,
     1                -0.00862, -0.00862, -0.00862, -0.00862, -0.00862, -0.00862, -0.00862,
     2                -0.00862, -0.00862, -0.00763, -0.00663, -0.00546, -0.00463, -0.00399,
     3                -0.00347, -0.00264, -0.002, -0.000828, 0.0 /
      data SA_a0 / -0.00397, -0.000311, -0.00397, -0.00397, -0.00397, -0.00397, -0.00397,
     1             -0.00397, -0.00397, -0.00397, -0.00397, -0.00397, -0.00397, -0.00397,
     2             -0.00397, -0.00397, -0.00351, -0.00305, -0.00252, -0.00214, -0.00184,
     3             -0.0016, -0.00122, -0.000919, -0.000382, 0.0 /
      data Taiwan_a0 / -0.00787, -0.00514, -0.00787, -0.00787, -0.00787, -0.00787, -0.00787,
     1                 -0.00787, -0.00787, -0.00787, -0.00787, -0.00787, -0.00787, -0.00787,
     2                 -0.00787, -0.00787, -0.0068, -0.00605, -0.00498, -0.00423, -0.00364,
     3                 -0.00316, -0.00241, -0.00182, -0.000755, 0.0 /

      data Global_a0slab / -0.00255, -0.0019, -0.00255, -0.00255, -0.00255, -0.00255, -0.00255,
     1                     -0.00255, -0.00255, -0.00255, -0.00255, -0.00255, -0.00255, -0.00255,
     2                     -0.00255, -0.00255, -0.00211, -0.00187, -0.00154, -0.00131, -0.00113,
     3                     -0.000979, -0.000745, -0.000564, -0.000234, 0.0 /
      data Alaska_a0slab / -0.00227, -0.00238, -0.00219, -0.00219, -0.00219, -0.00219, -0.00219,
     1                     -0.00219, -0.00219, -0.00219, -0.00219, -0.00219, -0.00219, -0.00219,
     2                     -0.00219, -0.00219, -0.00189, -0.00168, -0.00139, -0.00118, -0.00101,
     3                     -0.00088, -0.00067, -0.000507, -0.00021, 0.0 /
      data Cascadia_a0slab / -0.00354, -0.00109, -0.00401, -0.00401, -0.00401, -0.00401, -0.00401,
     1                       -0.00401, -0.00401, -0.00401, -0.00401, -0.00401, -0.00401, -0.00401,
     2                       -0.00401, -0.00401, -0.00347, -0.00309, -0.00254, -0.00216, -0.00186,
     3                       -0.00161, -0.00123, -0.000929, -0.000385, 0.0 /
      data CAM_a0slab / -0.00238, -0.00192, -0.00217, -0.00217, -0.00217, -0.00217, -0.00217,
     1                  -0.00217, -0.00217, -0.00217, -0.00217, -0.00217, -0.00217, -0.00217,
     2                  -0.00217, -0.00217, -0.00188, -0.00167, -0.00138, -0.00117, -0.00101,
     3                  -0.000873, -0.000664, -0.000503, -0.000209, 0.0 /
      data Japan_a0slab / -0.00335, -0.00215, -0.00311, -0.00311, -0.00311, -0.00311, -0.00311,
     1                    -0.00311, -0.00311, -0.00311, -0.00311, -0.00311, -0.00311, -0.00311,
     2                    -0.00311, -0.00311, -0.00269, -0.00239, -0.00197, -0.00167, -0.00144,
     3                    -0.00125, -0.000952, -0.00072, -0.000299, 0.0 /
      data SA_a0slab / -0.00238, -0.00192, -0.00217, -0.00217, -0.00217, -0.00217, -0.00217,
     1                 -0.00217, -0.00217, -0.00217, -0.00217, -0.00217, -0.00217, -0.00217,
     2                 -0.00217, -0.00217, -0.00188, -0.00167, -0.00138, -0.00117, -0.00101,
     3                 -0.000873, -0.000664, -0.000503, -0.000209, 0.0 /
      data Taiwan_a0slab / -0.00362, -0.00366, -0.00355, -0.00355, -0.00355, -0.00355, -0.00355,
     1                     -0.00355, -0.00355, -0.00355, -0.00355, -0.00355, -0.00355, -0.00355,
     2                     -0.00355, -0.00355, -0.00307, -0.00273, -0.00225, -0.00191, -0.00164,
     3                     -0.00143, -0.00109, -0.000822, -0.000341, 0.0 /

C     Source Scaling
      data c4 / 1.246, 1.336, 1.246, 1.227, 1.221, 1.215, 1.207, 1.201, 1.19, 1.182, 1.171,
     1          1.163, 1.156, 1.151, 1.143, 1.143, 1.217, 1.27, 1.344, 1.396, 1.437, 1.47,
     2          1.523, 1.564, 1.638, 1.69 /
      data c5 / -0.021, -0.039, -0.021, -0.021, -0.021, -0.021, -0.021, -0.021, -0.021,
     1          -0.021, -0.021, -0.021, -0.021, -0.021, -0.022, -0.023, -0.026, -0.028,
     2          -0.031, -0.034, -0.036, -0.038, -0.044, -0.048, -0.059, -0.067 /

      data c6 / 1.128, 1.336, 1.128, 1.128, 1.128, 1.128, 1.128, 1.128, 1.128, 1.128, 1.162, 1.163,
     1          1.156, 1.151, 1.143, 1.143, 1.217, 1.240, 1.237, 1.232, 1.227, 1.223, 1.216, 1.21,
     2          1.2, 1.194 /
c      data c6 /1.128, 1.844, 1.128, 1.128, 1.128, 1.128, 1.128, 1.128, 1.128, 1.128, 1.162, 1.187,
c     1         1.204, 1.215, 1.227, 1.234, 1.24, 1.24, 1.237, 1.232, 1.227, 1.223, 1.216, 1.21,
c     2         1.2, 1.194 /

      data c4slab / 1.84, 1.84, 1.84, 1.884, 1.884, 1.884, 1.884, 1.884, 1.884, 1.884, 1.884,
     1              1.884, 1.884, 1.884, 1.884, 1.884, 1.884, 1.884, 1.884, 1.884, 1.884, 1.949,
     2              2.031, 2.131, 2.185, 2.35 /
      data c5slab / -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.06,
     1              -0.068, -0.075, -0.082, -0.091, -0.1, -0.115, -0.134, -0.154, -0.154, -0.154,
     2              -0.154, -0.154, -0.154, -0.154, -0.154 /
      data c6slab / 0.4, 0.8, 0.4, 0.415, 0.43, 0.445, 0.46, 0.475, 0.49, 0.505, 0.52, 0.535, 0.55,
     1              0.565, 0.58, 0.595, 0.61, 0.625, 0.64, 0.655, 0.67, 0.685, 0.7, 0.715, 0.73, 0.745 /

      data d / 0.3004, 0.2693, 0.2839, 0.2854, 0.2891, 0.2932, 0.3004, 0.3048, 0.2992, 0.2854, 0.2814,
     1         0.291, 0.2758, 0.2719, 0.2539, 0.2482, 0.2227, 0.1969, 0.1452, 0.06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /
      data m / 0.0314, 0.0252, 0.0296, 0.0298, 0.0302, 0.0306, 0.0313, 0.0316, 0.0321, 0.032, 0.0325, 0.0306,
     1         0.0306, 0.0323, 0.0302, 0.0295, 0.0266, 0.0231, 0.0118, 0.007, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /
      data db / 67.0, 67.0, 67.0, 67.0, 67.0, 67.0, 67.0, 67.0, 67.0, 67.0, 67.0, 67.0,
     1          67.0, 67.0, 67.0, 67.0, 67.0, 67.0, 67.0, 67.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /

C     Linear Site Amplification, Fs
      data Vc / 1350.0, 850.0, 1300.0, 1225.0, 1200.0, 1200.0, 1200.0, 1225.0, 1350.0, 1450.0,
     1          1500.0, 1425.0, 1350.0, 1250.0, 1150.0, 1025.0, 900.0, 800.0, 760.0,
     2          760.0, 760.0, 760.0, 760.0, 760.0, 760.0, 760.0  /

      data Japan_c1 / -0.586, -0.738, -0.604, -0.593, -0.569, -0.539, -0.468, -0.403, -0.325,
     1                -0.264, -0.25, -0.288, -0.36, -0.455, -0.617, -0.757, -0.966, -0.986,
     2                -0.966, -0.901, -0.822, -0.751, -0.68, -0.592, -0.494, -0.395 /
      data Taiwan_c1 / -0.44, -0.454, -0.44, -0.458, -0.454, -0.455, -0.453, -0.452, -0.456, -0.468,
     1                 -0.484, -0.498, -0.511, -0.514, -0.51, -0.506, -0.5, -0.49, -0.486, -0.475,
     2                 -0.453, -0.428, -0.396, -0.353, -0.311, -0.261 /
      data Global_c / -0.498, -0.601, -0.498, -0.478, -0.464, -0.446, -0.431, -0.42, -0.442,
     1                -0.485, -0.546, -0.612, -0.688, -0.748, -0.802, -0.845, -0.911, -0.926,
     2                -0.888, -0.808, -0.743, -0.669, -0.585, -0.506, -0.418, -0.321 /
      data Alaska_c / -0.785, -1.031, -0.803, -0.785, -0.745, -0.69, -0.636, -0.594, -0.586,
     1                -0.629, -0.729, -0.867, -1.011, -1.133, -1.238, -1.321, -1.383, -1.414,
     2                -1.43, -1.421, -1.391, -1.343, -1.297, -1.233, -1.147, -1.06 /
      data Cascadia_c / -0.572, -0.671, -0.571, -0.575, -0.573, -0.565, -0.546, -0.519, -0.497,
     1                  -0.486, -0.499, -0.533, -0.592, -0.681, -0.772, -0.838, -0.922, -0.932,
     2                  -0.814, -0.725, -0.632, -0.57, -0.489, -0.421, -0.357, -0.302 /
      data Japan_c / -0.586, -0.738, -0.604, -0.593, -0.579, -0.561, -0.508, -0.461, -0.452,
     1               -0.498, -0.568, -0.667, -0.781, -0.867, -0.947, -1.003, -1.052, -1.028,
     2               -0.971, -0.901, -0.822, -0.751, -0.68, -0.592, -0.52, -0.395 /
      data SA_c / -0.333, -0.681, -0.333, -0.345, -0.362, -0.38, -0.403, -0.427, -0.458, -0.49,
     1            -0.536, -0.584, -0.654, -0.725, -0.801, -0.863, -0.942, -0.96, -0.942, -0.891,
     2            -0.842, -0.787, -0.706, -0.621, -0.52, -0.42 /
      data Taiwan_c / -0.44, -0.59, -0.44, -0.458, -0.459, -0.464, -0.466, -0.468, -0.473, -0.482,
     1                -0.499, -0.522, -0.555, -0.596, -0.643, -0.689, -0.745, -0.777, -0.79,
     2                -0.765, -0.724, -0.675, -0.613, -0.536, -0.444, -0.352 /

C     Nonlinear Site Amplification, Fnl
      data f4 / -0.44169, -0.31763, -0.4859, -0.4859, -0.4859, -0.4908, -0.49569, -0.49823,
     1          -0.49724, -0.49471, -0.48583, -0.47383, -0.47696, -0.4845, -0.48105, -0.46492,
     2          -0.43439, -0.38484, -0.32318, -0.26577, -0.21236, -0.17807, -0.13729, -0.07733,
     3          -0.05443, -0.03313 /
      data f5 / -0.0052, -0.0052, -0.0052, -0.00518, -0.00515, -0.00511, -0.00505, -0.00497,
     1          -0.00489, -0.00478, -0.0046, -0.00434, -0.00402, -0.0037, -0.00342, -0.00322,
     2          -0.00312, -0.0031, -0.0031, -0.0031, -0.0031, -0.0031, -0.0031, -0.0031,
     3          -0.0031, -0.0031 /

C     Basin Terms
      data Japan_e1 / 0.0, -0.137, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.1, 0.164, 0.164,
     1                0.08, 0.0, -0.13, -0.2, -0.401, -0.488, -0.578, -0.645, -0.678, -0.772,
     2               -0.699, -0.642, -0.524, -0.327 /
      data Japan_e2 / 0.0, 0.137, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.043, -0.085, -0.139, -0.139,
     1               -0.08, 0.0, 0.113, 0.176, 0.284, 0.346, 0.48, 0.579, 0.609, 0.635, 0.709,
     2                0.63, 0.306, 0.182 /
      data Japan_e3 / 1.0, 0.091, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -0.025, -0.05, -0.082, -0.082,
     1               -0.053, 1.0, 0.087, 0.118, 0.167, 0.203, 0.24, 0.254, 0.267, 0.265, 0.259,
     2                0.215, 0.175, 0.121 /
      data Cas_e1 / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.3, 0.333, 0.29, 0.177, 0.1, 0.0,
     1              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /
      data Cas_e2 / 0.0, 0.115, 0, 0, 0, 0, 0, -0.1, -0.34, -0.377, -0.29, -0.192, -0.035, 0.0,
     1              0.05, 0.1, 0.2, 0.245, 0.32, 0.37, 0.4, 0.43, 0.44, 0.45, 0.406, 0.345 /
      data Cas_e3 / 1.0, 0.068, 1.0, 1.0, 1.0, 1.0, 1.0, -0.063, -0.2, -0.222, -0.193, -0.148,
     1             -0.054, 1.0, 0.2, 0.2, 0.125, 0.153, 0.2, 0.239, 0.264, 0.287, 0.303,
     2              0.321, 0.312, 0.265 /
      data BasinD / 0.0, -0.115, 0.0, 0.0, 0.0, 0.0, 0.0, -0.05, -0.075, -0.081, -0.091,
     1             -0.092, 0.0, 0.0, 0.0, 0.0, -0.2, -0.245, -0.32, -0.28, -0.313, -0.355,
     2             -0.417, -0.45, -0.35, -0.331 /
      data SeaBasinD / 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.078, 0.075, 0.064, 0.075, 0.0,
     1                 0.0, 0.0, 0.0, 0.012, 0.037, 0.064, 0.14, 0.19, 0.165, 0.163, 0.132,
     2                 0.15, 0.117 /

C     Sigma Model
      data sigtau / 0.48, 0.477, 0.476, 0.482, 0.49, 0.5, 0.515, 0.528, 0.53, 0.524, 0.51,
     1              0.501, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492, 0.492,
     2              0.492, 0.492, 0.492, 0.492, 0.492 /
      data sigphi1 / 0.396, 0.348, 0.397, 0.401, 0.405, 0.413, 0.439, 0.473, 0.529, 0.517, 0.457,
     1               0.432, 0.45, 0.436, 0.433, 0.428, 0.448, 0.43, 0.406, 0.393, 0.381, 0.367,
     2               0.33, 0.298, 0.254, 0.231 /
      data sigphi2 / 0.565, 0.288, 0.56, 0.563, 0.575, 0.589, 0.616, 0.653, 0.722, 0.712,
     1               0.644, 0.64, 0.633, 0.584, 0.556, 0.51, 0.471, 0.43, 0.406, 0.393, 0.381,
     2               0.367, 0.33, 0.298, 0.254, 0.231 /
      data sigphiv / -0.18, -0.179, -0.18, -0.181, -0.183, -0.188, -0.205, -0.23, -0.262,
     1               -0.239, -0.185, -0.138, -0.185, -0.158, -0.19, -0.186, -0.177, -0.166,
     2               -0.111, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /

C     Constant parameters
      b4 = 0.1
      f1 = 0.0
      f3 = 0.05
      Vb = 200.0
      Vref = 760.0
      db1 = 20.0
      V1 = 270.0

C Find the requested spectral period and corresponding coefficients
      nPer = 26

C First check for the PGA case
      if (specT .eq. 0.0) then
         i1=1
         period2 = period(i1)
         c1T = c1(i1)
         c1slabT = c1slab(i1)
         c4T = c4(i1)
         c5T = c5(i1)
         c6T = c6(i1)
         c4slabT = c4slab(i1)
         c5slabT = c5slab(i1)
         c6slabT = c6slab(i1)
         f4T = f4(i1)
         f5T = f5(i1)
         VcT = Vc(i1)
         dT = d(i1)
         mT = m(i1)
         dbT = db(i1)

         global_a0T = Global_a0(i1)
         Alaska_a0T = Alaska_a0(i1)
         CAM_a0T = CAM_a0(i1)
         Japan_a0T = Japan_a0(i1)
         SA_a0T = SA_a0(i1)
         Taiwan_a0T = Taiwan_a0(i1)

         global_a0slabT = Global_a0slab(i1)
         Alaska_a0slabT = Alaska_a0slab(i1)
         Cascadia_a0slabT = Cascadia_a0slab(i1)
         CAM_a0slabT = CAM_a0slab(i1)
         Japan_a0slabT = Japan_a0slab(i1)
         SA_a0slabT = SA_a0slab(i1)
         Taiwan_a0slabT = Taiwan_a0slab(i1)

         global_c0T = Global_c0(i1)
         Alaska_c0T = Alaska_c0(i1)
         Aleutian_c0T = Aleutian_c0(i1)
         Cascadia_c0T = Cascadia_c0(i1)
         CAM_c0T = CAM_c0(i1)
         Japan_Pac_c0T = Japan_Pac_c0(i1)
         Japan_Phi_c0T = Japan_Phi_c0(i1)
         SAN_c0T = SAN_c0(i1)
         SAS_c0T = SAS_c0(i1)
         Taiwan_c0T = Taiwan_c0(i1)

         global_c0sT = Global_c0s(i1)
         Cascadia_c0sT = Cascadia_c0s(i1)
         Alaska_c0sT = Alaska_c0s(i1)
         Aleutian_c0sT = Aleutian_c0s(i1)
         CAM_c0sT = CAM_c0s(i1)
         Japan_c0sT = Japan_c0s(i1)
         SAN_c0sT = SAN_c0s(i1)
         SAS_c0sT = SAS_c0s(i1)
         Taiwan_c0sT = Taiwan_c0s(i1)

         global_cT = Global_c(i1)
         Alaska_cT = Alaska_c(i1)
         Cascadia_cT = Cascadia_c(i1)
         Japan_cT = Japan_c(i1)
         Japan_c1T = Japan_c1(i1)
         SA_cT = SA_c(i1)
         Taiwan_cT = Taiwan_c(i1)
         Taiwan_c1T = Taiwan_c1(i1)

         Japan_e1T = Japan_e1(i1)
         Japan_e2T = Japan_e2(i1)
         Japan_e3T = Japan_e3(i1)
         Cas_e1T = Cas_e1(i1)
         Cas_e2T = Cas_e2(i1)
         Cas_e3T = Cas_e3(i1)
         SeaBasinDT = SeaBasinD(i1)
         BasinDT = BasinD(i1)

         sigPhi1T = sigphi1(i1)
         sigPhi2T = sigphi2(i1)
         sigPhivT = sigphiv(i1)
         sigtauT = sigtau(i1)

         goto 1011

      elseif (specT .eq. -1.0) then
         i1=2
         period2 = period(i1)
         c1T = c1(i1)
         c1SlabT = c1Slab(i1)
         c4T = c4(i1)
         c5T = c5(i1)
         c6T = c6(i1)
         c4slabT = c4slab(i1)
         c5slabT = c5slab(i1)
         c6slabT = c6slab(i1)
         f4T = f4(i1)
         f5T = f5(i1)
         VcT = Vc(i1)
         dT = d(i1)
         mT = m(i1)
         dbT = db(i1)

         global_a0T = Global_a0(i1)
         Alaska_a0T = Alaska_a0(i1)
         CAM_a0T = CAM_a0(i1)
         Japan_a0T = Japan_a0(i1)
         SA_a0T = SA_a0(i1)
         Taiwan_a0T = Taiwan_a0(i1)

         global_a0slabT = Global_a0slab(i1)
         Alaska_a0slabT = Alaska_a0slab(i1)
         Cascadia_a0slabT = Cascadia_a0slab(i1)
         CAM_a0slabT = CAM_a0slab(i1)
         Japan_a0slabT = Japan_a0slab(i1)
         SA_a0slabT = SA_a0slab(i1)
         Taiwan_a0slabT = Taiwan_a0slab(i1)

         global_c0T = Global_c0(i1)
         Alaska_c0T = Alaska_c0(i1)
         Aleutian_c0T = Aleutian_c0(i1)
         Cascadia_c0T = Cascadia_c0(i1)
         CAM_c0T = CAM_c0(i1)
         Japan_Pac_c0T = Japan_Pac_c0(i1)
         Japan_Phi_c0T = Japan_Phi_c0(i1)
         SAN_c0T = SAN_c0(i1)
         SAS_c0T = SAN_c0(i1)
         Taiwan_c0T = Taiwan_c0(i1)

         global_c0sT = Global_c0s(i1)
         Cascadia_c0sT = Cascadia_c0s(i1)
         Alaska_c0sT = Alaska_c0s(i1)
         Aleutian_c0sT = Aleutian_c0s(i1)
         CAM_c0sT = CAM_c0s(i1)
         Japan_c0sT = Japan_c0s(i1)
         SAN_c0sT = SAN_c0s(i1)
         SAS_c0sT = SAS_c0s(i1)
         Taiwan_c0sT = Taiwan_c0s(i1)

         global_cT = Global_c(i1)
         Alaska_cT = Alaska_c(i1)
         Cascadia_cT = Cascadia_c(i1)
         Japan_cT = Japan_c(i1)
         Japan_c1T = Japan_c1(i1)
         SA_cT = SA_c(i1)
         Taiwan_cT = Taiwan_c(i1)
         Taiwan_c1T = Taiwan_c1(i1)

         Japan_e1T = Japan_e1(i1)
         Japan_e2T = Japan_e2(i1)
         Japan_e3T = Japan_e3(i1)
         Cas_e1T = Cas_e1(i1)
         Cas_e2T = Cas_e2(i1)
         Cas_e3T = Cas_e3(i1)
         SeaBasinDT = SeaBasinD(i1)
         BasinDT = BasinD(i1)

         sigPhi1T = sigphi1(i1)
         sigPhi2T = sigphi2(i1)
         sigPhivT = sigphiv(i1)
         sigtauT = sigtau(i1)

         goto 1011
      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=3,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020
         endif
      enddo
      write (*,'( i5,2f12.6)') nper, specT, period(nper)

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*)
      write (*,*) 'PSHAB (2020) Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*)
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),c1(count1),c1(count2),
     +                   specT,c1T,iflag)
      call S24_interp (period(count1),period(count2),c1slab(count1),c1slab(count2),
     +                   specT,c1slabT,iflag)
      call S24_interp (period(count1),period(count2),c4(count1),c4(count2),
     +                   specT,c4T,iflag)
      call S24_interp (period(count1),period(count2),c5(count1),c5(count2),
     +                   specT,c5T,iflag)
      call S24_interp (period(count1),period(count2),c6(count1),c6(count2),
     +                   specT,c6T,iflag)
      call S24_interp (period(count1),period(count2),c4slab(count1),c4slab(count2),
     +                   specT,c4slabT,iflag)
      call S24_interp (period(count1),period(count2),c5slab(count1),c5slab(count2),
     +                   specT,c5slabT,iflag)
      call S24_interp (period(count1),period(count2),c6slab(count1),c6slab(count2),
     +                   specT,c6slabT,iflag)
      call S24_interp (period(count1),period(count2),f4(count1),f4(count2),
     +                   specT,f4T,iflag)
      call S24_interp (period(count1),period(count2),f5(count1),f5(count2),
     +                   specT,f5T,iflag)
      call S24_interp (period(count1),period(count2),Vc(count1),Vc(count2),
     +                   specT,VcT,iflag)
      call S24_interp (period(count1),period(count2),Global_a0(count1),Global_a0(count2),
     +                   specT,Global_a0T,iflag)
      call S24_interp (period(count1),period(count2),Alaska_a0(count1),Alaska_a0(count2),
     +                   specT,Alaska_a0T,iflag)
      call S24_interp (period(count1),period(count2),CAM_a0(count1),CAM_a0(count2),
     +                   specT,CAM_a0T,iflag)
      call S24_interp (period(count1),period(count2),Japan_a0(count1),Japan_a0(count2),
     +                   specT,Japan_a0T,iflag)
      call S24_interp (period(count1),period(count2),SA_a0(count1),SA_a0(count2),
     +                   specT,SA_a0T,iflag)
      call S24_interp (period(count1),period(count2),Taiwan_a0(count1),Taiwan_a0(count2),
     +                   specT,Taiwan_a0T,iflag)
      call S24_interp (period(count1),period(count2),Global_a0slab(count1),Global_a0slab(count2),
     +                   specT,Global_a0slabT,iflag)
      call S24_interp (period(count1),period(count2),Alaska_a0slab(count1),Alaska_a0slab(count2),
     +                   specT,Alaska_a0slabT,iflag)
      call S24_interp (period(count1),period(count2),Cascadia_a0slab(count1),Cascadia_a0slab(count2),
     +                   specT,Cascadia_a0slabT,iflag)
      call S24_interp (period(count1),period(count2),CAM_a0slab(count1),CAM_a0slab(count2),
     +                   specT,CAM_a0slabT,iflag)
      call S24_interp (period(count1),period(count2),Japan_a0slab(count1),Japan_a0slab(count2),
     +                   specT,Japan_a0slabT,iflag)
      call S24_interp (period(count1),period(count2),SA_a0slab(count1),SA_a0slab(count2),
     +                   specT,SA_a0slabT,iflag)
      call S24_interp (period(count1),period(count2),Taiwan_a0slab(count1),Taiwan_a0slab(count2),
     +                   specT,Taiwan_a0slabT,iflag)
      call S24_interp (period(count1),period(count2),Global_c0(count1),Global_c0(count2),
     +                   specT,Global_c0T,iflag)
      call S24_interp (period(count1),period(count2),Alaska_c0(count1),Alaska_c0(count2),
     +                   specT,Alaska_c0T,iflag)
      call S24_interp (period(count1),period(count2),Aleutian_c0(count1),Aleutian_c0(count2),
     +                   specT,Aleutian_c0T,iflag)
      call S24_interp (period(count1),period(count2),Cascadia_c0(count1),Cascadia_c0(count2),
     +                   specT,Cascadia_c0T,iflag)
      call S24_interp (period(count1),period(count2),CAM_c0(count1),CAM_c0(count2),
     +                   specT,CAM_c0T,iflag)
      call S24_interp (period(count1),period(count2),Japan_Pac_c0(count1),Japan_Pac_c0(count2),
     +                   specT,Japan_Pac_c0T,iflag)
      call S24_interp (period(count1),period(count2),Japan_Phi_c0(count1),Japan_Phi_c0(count2),
     +                   specT,Japan_Phi_c0T,iflag)
      call S24_interp (period(count1),period(count2),SAN_c0(count1),SAN_c0(count2),
     +                   specT,SAN_c0T,iflag)
      call S24_interp (period(count1),period(count2),SAS_c0(count1),SAS_c0(count2),
     +                   specT,SAS_c0T,iflag)
      call S24_interp (period(count1),period(count2),Taiwan_c0(count1),Taiwan_c0(count2),
     +                   specT,Taiwan_c0T,iflag)
      call S24_interp (period(count1),period(count2),Global_c0s(count1),Global_c0s(count2),
     +                   specT,Global_c0sT,iflag)
      call S24_interp (period(count1),period(count2),Cascadia_c0s(count1),Cascadia_c0s(count2),
     +                   specT,Cascadia_c0sT,iflag)
      call S24_interp (period(count1),period(count2),Alaska_c0s(count1),Alaska_c0s(count2),
     +                   specT,Alaska_c0sT,iflag)
      call S24_interp (period(count1),period(count2),Aleutian_c0s(count1),Aleutian_c0s(count2),
     +                   specT,Aleutian_c0sT,iflag)
      call S24_interp (period(count1),period(count2),CAM_c0s(count1),CAM_c0s(count2),
     +                   specT,CAM_c0sT,iflag)
      call S24_interp (period(count1),period(count2),Japan_c0s(count1),Japan_c0s(count2),
     +                   specT,Japan_c0sT,iflag)
      call S24_interp (period(count1),period(count2),SAN_c0s(count1),SAN_c0s(count2),
     +                   specT,SAN_c0sT,iflag)
      call S24_interp (period(count1),period(count2),SAS_c0s(count1),SAS_c0s(count2),
     +                   specT,SAS_c0sT,iflag)
      call S24_interp (period(count1),period(count2),Taiwan_c0s(count1),Taiwan_c0s(count2),
     +                   specT,Taiwan_c0sT,iflag)
      call S24_interp (period(count1),period(count2),Global_c(count1),Global_c(count2),
     +                   specT,Global_cT,iflag)
      call S24_interp (period(count1),period(count2),Alaska_c(count1),Alaska_c(count2),
     +                   specT,Alaska_cT,iflag)
      call S24_interp (period(count1),period(count2),Cascadia_c(count1),Cascadia_c(count2),
     +                   specT,Cascadia_cT,iflag)
      call S24_interp (period(count1),period(count2),Japan_c(count1),Japan_c(count2),
     +                   specT,Japan_cT,iflag)
      call S24_interp (period(count1),period(count2),Japan_c1(count1),Japan_c1(count2),
     +                   specT,Japan_c1T,iflag)
      call S24_interp (period(count1),period(count2),SA_c(count1),SA_c(count2),
     +                   specT,SA_cT,iflag)
      call S24_interp (period(count1),period(count2),Taiwan_c(count1),Taiwan_c(count2),
     +                   specT,Taiwan_cT,iflag)
      call S24_interp (period(count1),period(count2),Taiwan_c1(count1),Taiwan_c1(count2),
     +                   specT,Taiwan_c1T,iflag)
      call S24_interp (period(count1),period(count2),d(count1),d(count2),
     +                   specT,dT,iflag)
      call S24_interp (period(count1),period(count2),m(count1),m(count2),
     +                   specT,mT,iflag)
      call S24_interp (period(count1),period(count2),db(count1),db(count2),
     +                   specT,dbT,iflag)
      call S24_interp (period(count1),period(count2),Japan_e1(count1),japan_e1(count2),
     +                   specT,Japan_e1T,iflag)
      call S24_interp (period(count1),period(count2),Japan_e2(count1),japan_e2(count2),
     +                   specT,Japan_e2T,iflag)
      call S24_interp (period(count1),period(count2),Japan_e3(count1),japan_e3(count2),
     +                   specT,Japan_e3T,iflag)
      call S24_interp (period(count1),period(count2),Cas_e1(count1),Cas_e1(count2),
     +                   specT,Cas_e1T,iflag)
      call S24_interp (period(count1),period(count2),Cas_e2(count1),Cas_e2(count2),
     +                   specT,Cas_e2T,iflag)
      call S24_interp (period(count1),period(count2),CAS_e3(count1),Cas_e3(count2),
     +                   specT,Cas_e3T,iflag)
      call S24_interp (period(count1),period(count2),SeaBasinD(count1),SeaBasinD(count2),
     +                   specT,SeaBasinDT,iflag)
      call S24_interp (period(count1),period(count2),BasinD(count1),BasinD(count2),
     +                   specT,BasinDT,iflag)
      call S24_interp (period(count1),period(count2),sigPhi1(count1),sigPhi1(count2),
     +                   specT,sigPhi1T,iflag)
      call S24_interp (period(count1),period(count2),sigPhi2(count1),sigPhi2(count2),
     +                   specT,sigPhi2T,iflag)
      call S24_interp (period(count1),period(count2),sigPhiv(count1),sigPhiv(count2),
     +                   specT,sigPhivT,iflag)
      call S24_interp (period(count1),period(count2),sigtau(count1),sigtau(count2),
     +                   specT,sigtauT,iflag)

 1011 period2 = specT

C     Set the regional terms
c     0=Global
c     1=Alaska
c     2=Aleutian (a0=Alaska)
c     3=Cascadia (a0=Global(Interface))
c     4=Central America / Mexico (Global Linear site term)
c     5=Japan - Pacific Plate
c     6=Japan - Philippine Plate
c     7=Northern South America
c     8=Southern South America
c     9=Taiwan

      if (iRegion .eq. 0) then
         if (ftype .eq. 0.0) then
            c0T = Global_c0T
            c0Pga = Global_c0(1)
            a0T = Global_a0T
            a0Pga = Global_a0(1)
         else
            c0T = Global_c0sT
            c0Pga = Global_c0s(1)
            a0T = Global_a0slabT
            a0Pga = Global_a0slab(1)
         endif
         cT = Global_cT
      elseif (iRegion .eq. 1) then
         if (ftype .eq. 0.0) then
            c0T = Alaska_c0T
            c0Pga = Alaska_c0(1)
            a0T = Alaska_a0T
            a0Pga = Alaska_a0(1)
         else
            c0T = Alaska_c0sT
            c0Pga = Alaska_c0s(1)
            a0T = Alaska_a0slabT
            a0Pga = Alaska_a0slab(1)
         endif
         cT = Alaska_cT
      elseif (iRegion .eq. 2) then
         if (ftype .eq. 0.0) then
            c0T = Aleutian_c0T
            c0Pga = Aleutian_c0(1)
            a0T = Alaska_a0T
            a0Pga = Alaska_a0(1)
         else
            c0T = Aleutian_c0sT
            c0Pga = Aleutian_c0s(1)
            a0T = Alaska_a0slabT
            a0Pga = Alaska_a0slab(1)
         endif
         cT = Alaska_cT
      elseif (iRegion .eq. 3) then
         if (ftype .eq. 0.0) then
            c0T = Cascadia_c0T
            c0Pga = Cascadia_c0(1)
            a0T = Global_a0T
            a0Pga = Global_a0(1)
         else
            c0T = Cascadia_c0sT
            c0Pga = Cascadia_c0s(1)
            a0T = Cascadia_a0SlabT
            a0Pga = Cascadia_a0Slab(1)
         endif
         cT = Cascadia_cT
      elseif (iRegion .eq. 4) then
         if (ftype .eq. 0.0) then
            c0T = CAM_c0T
            c0Pga = CAM_c0(1)
            a0T = CAM_a0T
            a0Pga = CAM_a0(1)
         else
            c0T = CAM_c0sT
            c0Pga = CAM_c0s(1)
            a0T = CAM_a0slabT
            a0Pga = CAM_a0slab(1)
         endif
         cT = Global_cT
      elseif (iRegion .eq. 5) then
         if (ftype .eq. 0.0) then
            c0T = Japan_Pac_c0T
            c0Pga = Japan_Pac_c0(1)
            a0T = Japan_a0T
            a0Pga = Japan_a0(1)
         else
            c0T = Japan_c0sT
            c0Pga = Japan_c0s(1)
            a0T = Japan_a0slabT
            a0Pga = Japan_a0slab(1)
         endif
         cT = Japan_cT
      elseif (iRegion .eq. 6) then
         if (ftype .eq. 0.0) then
            c0T = Japan_Phi_c0T
            c0Pga = Japan_Phi_c0(1)
            a0T = Japan_a0T
            a0Pga = Japan_a0(1)
         else
            c0T = Japan_c0sT
            c0Pga = Japan_c0s(1)
            a0T = Japan_a0slabT
            a0Pga = Japan_a0slab(1)
         endif
         cT = Japan_cT
      elseif (iRegion .eq. 7) then
         if (ftype .eq. 0.0) then
            c0T = SAN_c0T
            c0Pga = SAN_c0(1)
            a0T = SA_a0T
            a0Pga = SA_a0(1)
         else
            c0T = SAN_c0sT
            c0Pga = SAN_c0s(1)
            a0T = SA_a0slabT
            a0Pga = SA_a0slab(1)
         endif
         cT = SA_cT
      elseif (iRegion .eq. 8) then
         if (ftype .eq. 0.0) then
            c0T = SAS_c0T
            c0Pga = SAS_c0(1)
            a0T = SA_a0T
            a0Pga = SA_a0(1)
         else
            c0T = SAS_c0sT
            c0Pga = SAS_c0s(1)
            a0T = SA_a0slabT
            a0Pga = SA_a0slab(1)
         endif
         cT = SA_cT
      elseif (iRegion .eq. 9) then
         if (ftype .eq. 0.0) then
            c0T = Taiwan_c0T
            c0Pga = Taiwan_c0(1)
            a0T = Taiwan_a0T
            a0Pga = Taiwan_a0(1)
         else
            c0T = Taiwan_c0sT
            c0Pga = Taiwan_c0s(1)
            a0T = Taiwan_a0slabT
            a0Pga = Taiwan_a0slab(1)
         endif
         cT = Taiwan_cT
      endif

      if (ftype .eq. 0.0) then
         h = 10.0**(-0.82+0.252*mag)
      else
C     New Slab Model for H
          if (mag .ge. mbslab) then
             h = 35.0
          else
             h = 10.0**((1.05/(mbslab-4.0))*(mag-mbslab) + 1.544)
          endif
      endif
      R = sqrt(Rrup*Rrup + h*h)
      Rref = sqrt(1.0+h*h)

C     Compute Fp term
      if (Ftype .eq. 0.0) then
         Fp = c1T*alog(R) + b4*mag*alog(R/Rref) + a0T*R
         FpPga = c1(1)*alog(R) + b4*mag*alog(R/Rref) + a0Pga*R
      else
         Fp = c1slabT*alog(R) + b4*mag*alog(R/Rref) + a0T*R
         FpPga = c1slab(1)*alog(R) + b4*mag*alog(R/Rref) + a0Pga*R
      endif

C     Compute the Fm term
      if (ftype .eq. 0.0) then
         if (mag .le. mbinter) then
            Fm = c4T*(mag-mbinter) + c5T*(mag-mbinter)**2.0
            FmPga = c4(1)*(mag-mbinter) + c5(1)*(mag-mbinter)**2.0
         else
            Fm = c6T*(mag-mbinter)
            FmPga = c6(1)*(mag-mbinter)
         endif
      else
         if (mag .le. mbslab) then
            Fm = c4slabT*(mag-mbslab) + c5slabT*(mag-mbslab)**2.0
            FmPga = c4slab(1)*(mag-mbslab) + c5slab(1)*(mag-mbslab)**2.0
         else
            Fm = c6slabT*(mag-mbslab)
            FmPga = c6slab(1)*(mag-mbslab)
         endif
      endif

C     Compute the Fd term
      if (Ftype .eq. 0.0) then
         Fd = 0.0
         FdPga = 0.0
      else
         if (depth .le. 20.0) then
            Fd = mT*(20.0 - dbT) + dT
            FdPga = m(1)*(20.0 - db(1)) + d(1)
         elseif (depth .gt. dbT) then
            Fd = dT
            FdPga = d(1)
         else
            Fd = mT*(depth-dbT) + dT
            FdPga = m(1)*(depth-db(1)) + d(1)
         endif
      endif

C     Compute PGARock

      pgarock = exp(c0Pga + FpPga + FmPga + FdPga)

C     Now compute the ground motion for given Vs30 value.

C     Compute the linear Flin term.
      if (Vs30 .le. V1) then
         if (iRegion .eq. 5)then
            Flin = Japan_c1T*alog(Vs30/V1) + cT*alog(V1/Vref)
         elseif (iRegion .eq. 8) then
            Flin = Taiwan_c1T*alog(Vs30/V1) + cT*alog(V1/Vref)
         else
            Flin = cT*alog(Vs30/V1) + cT*alog(V1/Vref)
         endif
      elseif (Vs30 .gt. VcT) then
         Flin = cT*alog(VcT/Vref)
      else
         Flin = cT*alog(Vs30/Vref)
      endif

C     Now Compute the Non-linear term (Fnl)
      f2 = f4T*(exp(f5T*(min(vs30,760.0)-200)) - exp(f5T*(760-200)))
      if (specT .ge. 3.0) then
         Fnl = 0.0
      else
         Fnl = f1 + f2*alog((pgarock+f3)/f3)
      endif

      Fs = Flin + Fnl

C     Compute the basin term
C     Cascadia Basin Terms
c      z25 = z25
      if (iRegion .eq. 3) then
          x = (alog10(vs30) - alog10(500.0) ) / (0.42*sqrt(2.0))
          muz25 = 10**(3.75-0.74*(1+erf(x)))
          deltaz25 = alog(z25*1000.0) - alog(muz25)
C     Adjust the e3 term based on the PNW Basin Flag
C     Note values of 0 and 1 are for previous models not currently recommended.
          if (pnwbflag .eq. 0) then
             cas_e3Tbasin = cas_e3T + basinDT
             cas_e2Tbasin = cas_e2T + basinDT
          elseif (pnwbflag .eq. 1) then
             cas_e3Tbasin = cas_e3T + SeabasinDT
             cas_e2Tbasin = cas_e2T + SeabasinDT
          else
             cas_e3Tbasin = cas_e3T
             cas_e2Tbasin = cas_e2T
          endif

          if (deltaz25 .le. (Cas_e1T/Cas_e3Tbasin) ) then
             fb = Cas_e1T
          elseif (deltaz25 .ge. (cas_e2Tbasin/Cas_e3Tbasin) ) then
             fb = cas_e2Tbasin
          else
             fb = Cas_e3TBasin*deltaz25
          endif

C     Japan Basin Terms
      elseif (iRegion .eq. 5 .or. iRegion .eq. 6) then
          x = (alog10(vs30) - alog10(500.0) ) / (0.33*sqrt(2.0))
          muz25 = 10**(3.05-0.8*(1+erf(x)))
          deltaz25 = alog(z25*1000.0) - alog(muz25)

          if (deltaz25 .le. (Japan_e1T/Japan_e3T) ) then
             fb = Japan_e1T
          elseif (deltaz25 .ge. (Japan_e2T/Japan_e3T) ) then
             fb = Japan_e2T
          else
             fb = Japan_e3T*deltaz25
          endif
      else
         fb = 0.0
      endif
C     Compute the site specific ground motion
      lnSA = c0T + Fp + Fm + Fd + Fs + Fb

c     Convert units spectral acceleration in gal
      if (spect .ne. -1.0) then
         lnSa = lnSa + 6.89
      endif

C     Compute sigma values
      if (rRup .le. 200.0) then
         phirup = sigphi1T
      elseif (rRup .ge. 500.0) then
         phirup = sigphi2T
      else
         phirup = ((sigphi2T - sigphi1T)/0.9163)*alog(rRup/200.0) + sigphi1T
      endif

      if (Vs30 .le. 200.0) then
         phivs = sigphivT*(alog(500.0/max(200.0,min(500.0,rRup)))/alog(500.0/200.0))
      elseif (Vs30 .ge. 500.0) then
         phivs = 0.0
      else
         phivs = sigphivT*(alog(500.0/Vs30)/alog(500.0/200.0))*(alog(500.0/(max(200.0,min(500.0,rRup))))/(alog(500.0/200.0)))
      endif

      phivar = phirup + phiVs
      phi = sqrt(phivar)

      tau = sigtauT

      sigma = sqrt(phi*phi + tau*tau)

      return
      end

c ----------------------------------------------------------------------

      subroutine S35_SMK2020 ( mag, Ftype, rRup, vs30, z25, lnSa, sigma, phi, tau,
     2                     specT, period2, iflag, depth, mohodepth )

      implicit none
      integer MAXPER
      parameter (MAXPER=23)
      integer count1, count2, iflag, nPer, i1, i

      real Period(MAXPER), e(MAXPER), a1(MAXPER), d1(MAXPER),
     1     a2(MAXPER), h(MAXPER), cd(MAXPER), Dd(MAXPER), C(MAXPER),
     2     Vc(MAXPER), f4(MAXPER), f5(MAXPER), sigphi(MAXPER), sigtau(MAXPER)
      real eT, a1T, d1T, a2T, hT, cdT, DdT, CT, VcT, f4T, f5T, sigPhiT, sigTauT
      real Vref, f1, f3
      real sigma, lnSa, pgaRock, vs30, rRup, mag, depth, Ftype, mohodepth
      real period2, specT, z25, phi, tau
      real kterm, ktermpga, cterm, ctermpga, cutDD, gterm, gtermpga, gdterm, gdtermpga
      real bterm, btermpga, flin, fnl, gsterm, f2, median


C     Coefficients 08/25/2021

      data period / 0.00, -1.00, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
     1              0.5, 0.75, 1.00, 1.5, 2.00, 3.00, 4.00, 5.00, 7.00, 10.00 /
      data e  / -2.7267439, -1.933444, -2.7079904, -2.7080162, -2.6903921, -2.6167275, -2.525663,
     1          -2.4914814, -2.5189582, -2.5989579, -2.7092901, -2.8172584, -3.0714905, -3.3136917,
     2          -3.785221385, -4.16659131, -4.901364443, -5.470871899, -6.337799541, -6.921473929,
     3          -7.319266537, -7.945565121, -8.779899021 /
      data a1 / 0.4908892, 0.6445959, 0.50175701, 0.503327193, 0.504897379, 0.508037773, 0.511178213,
     1          0.515888791, 0.523737608, 0.531579061, 0.539407049, 0.547215792, 0.56275672,
     2          0.578170887, 0.615999628, 0.652553714, 0.720559114, 0.77991566, 0.868274098,
     3          0.923481828, 0.96007551, 1.025041319, 1.137022039 /
      data d1 / 0.2118015, 0.1805972, 0.213370189, 0.223607193, 0.233844197, 0.246864583, 0.25625647,
     1          0.263091069, 0.259885833, 0.2515129, 0.243631882, 0.232014968, 0.217884443, 0.204814643,
     2          0.187899494, 0.176347416, 0.171592549, 0.174051474, 0.181359962, 0.187096027,
     3          0.193129321, 0.197171869, 0.203235691 /
      data a2 / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1348718,
     1          0.1302177, 0.2549352, 0.5137211, 0.5672251, 0.6083083, 0.6428719, 0.6942251, 0.5671919 /
      data h / 0.007210868, 0.005419657, 0.007202214, 0.007430087, 0.007657961, 0.007979789,
     1         0.008247439, 0.008447657, 0.008211716, 0.007787228, 0.00738587, 0.006845832,
     2         0.006224042, 0.005800522, 0.005048717, 0.004544581, 0.004168906, 0.004056775,
     3         0.004123161, 0.003972515, 0.003648234, 0.002182624, 0.0 /
      data Cd / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1          0.002032605, 0.025594005, 0.031828072, 0.030465723, -0.001921725, -0.042352667,
     2         -0.077322458, -0.100115037 /
      data Dd / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1          0.085030533, 0.097461559, 0.11418673, 0.133516203, 0.163964771, 0.196531265,
     2          0.199814063, 0.172679121 /
      data c / -0.6, -0.84, -0.6037, -0.5739, -0.5341, -0.458, -0.4441, -0.4872, -0.5796,
     1         -0.6876, -0.7718, -0.8417, -0.9109, -0.9693, -1.0154, -1.05, -1.0454, -1.0392,
     2         -1.0112, -0.9694, -0.9195, -0.8046, -0.6558 /
      data Vc / 1500.0, 1300.0, 1500.2, 1500.36, 1502.95, 1501.42, 1494.0, 1479.12, 1442.85,
     1          1392.61, 1356.21, 1308.47, 1252.66, 1203.91, 1147.59, 1109.95, 1072.39,
     2          1009.49, 922.43, 844.48, 793.13, 772.68, 775.0 /
      data f4 / -0.15, -0.1, -0.1483, -0.1471, -0.1549, -0.1963, -0.2287, -0.2492, -0.2571,
     1          -0.2466, -0.2357, -0.2191, -0.1958, -0.1704, -0.1387, -0.1052, -0.0679,
     2          -0.0361, -0.0136, -0.0032, -0.0003, 0.0, 0.0 /
      data f5 / -0.00701, -0.00844, -0.00701, -0.00728, -0.00735, -0.00647, -0.00573,
     1          -0.0056, -0.00585, -0.00614, -0.00644, -0.0067, -0.00713, -0.00744, -0.00812,
     2          -0.00844, -0.00771, -0.00479, -0.00183, -0.00152, -0.00144, -0.00137, -0.00136 /
      data sigphi / 0.602511929, 0.554426633, 0.587138878, 0.586636948, 0.591884361, 0.616636506,
     1              0.662828733, 0.696266316, 0.712217671, 0.703333328, 0.683114937, 0.664253842,
     2              0.635447131, 0.624947234, 0.63600237, 0.654737997, 0.692422766, 0.714571508,
     3              0.702827819, 0.685656575, 0.653101909, 0.606218533, 0.53921879 /
      data sigtau / 0.678907779, 0.458422498, 0.651303557, 0.654062863, 0.672310989, 0.738544673,
     1              0.788641603, 0.76633797, 0.651840928, 0.610739167, 0.567120468, 0.537991628,
     2              0.506166155, 0.486000441, 0.481555738, 0.484183713, 0.4689282, 0.45925422,
     3              0.47022716, 0.500884413, 0.510577956, 0.516643095, 0.493780604 /


c      data period / 0.000, -1.000, 0.010, 0.020, 0.030, 0.050, 0.075, 0.100, 0.150,
c     1              0.200, 0.250, 0.300, 0.400, 0.500, 0.750, 1.000, 1.500, 2.000,
c     2              3.000, 4.000, 5.000, 7.000, 10.000 /
c      data e / -2.7319602, -1.9377316, -2.6899752, -2.6889864, -2.668915, -2.5900945,
c     1         -2.5040933, -2.4680936, -2.5026276, -2.5946368, -2.7083263, -2.8238059,
c     2         -3.0871773, -3.3423891, -3.844262719, -4.26099559, -5.052034618,
c     3         -5.708858434, -6.731732648, -7.432432738, -7.871883886,
c     4         -8.274318056, -8.49984373 /
c      data a1 / 0.4917519, 0.6449087, 0.498552103, 0.500247032, 0.501941964, 0.505331864,
c     1          0.509569357, 0.513806919, 0.52228081, 0.530749161, 0.539206716, 0.547648573,
c     2          0.564469723, 0.581187948, 0.622416389, 0.662622688, 0.738788056, 0.807385789,
c     3          0.915403426, 0.98511283, 1.025265403, 1.056000115, 1.071204328 /
c      data d1 / 0.2114314, 0.1818077, 0.216457821, 0.223950035, 0.23144225, 0.241116032,
c     1          0.250082418, 0.255436924, 0.255226857, 0.25035302, 0.245524039, 0.235806898,
c     2          0.221382562, 0.208054723, 0.191214688, 0.180722026, 0.177385417, 0.180865794,
c     3          0.179712346, 0.176580661, 0.174348944, 0.167460034, 0.157126668  /
c      data a2 / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
c     1          0.1118103, 0.1339921, 0.2557821, 0.5168177, 0.5906224, 0.6243605, 0.6601545,
c     2          0.7026209, 0.617846 /
c      data h / 0.007203427, 0.005427202, 0.007287751, 0.007438551, 0.007589351, 0.007853247,
c     1         0.008122348, 0.008257189, 0.008109089, 0.007767407, 0.007413394, 0.006935728,
c     2         0.006314757, 0.005913315, 0.00514768, 0.004683866, 0.00441547, 0.004410852,
c     3         0.004658857, 0.004796962, 0.004725937, 0.003719351, 0.002209471 /
c      data Cd / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
c     1         -0.0008447, 0.01982, 0.02578, 0.02093, -0.01318, -0.05295, -0.08429, -0.10673 /
c      data Dd / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
c     1          0.0, 0.0849857, 0.10041, 0.1174, 0.13777, 0.16772, 0.20112, 0.20962, 0.19631 /
c      data c / -0.84, -0.6, -0.604, -0.574, -0.534, -0.458, -0.444, -0.487, -0.58,
c     1         -0.688, -0.772, -0.842, -0.911, -0.969, -1.02, -1.05, -1.05, -1.04,
c     2         -1.01, -0.969, -0.92, -0.805, -0.656 /
c      data Vc / 1300.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1490.0, 1480.0, 1440.0,
c     1          1390.0, 1360.0, 1310.0, 1250.0, 1200.0, 1150.0, 1110.0, 1070.0, 1010.0,
c     2           922.0, 844.0, 793.0, 773.0, 775.0 /
c      data f4 / -0.1, -0.15, -0.148, -0.147, -0.155, -0.196, -0.229, -0.249, -0.257,
c     1          -0.247, -0.236, -0.219, -0.196, -0.17, -0.139, -0.105, -0.0679, -0.0361,
c     2          -0.0136, -0.00321, -0.000255, -0.0000455, 0.0 /
c      data f5 / -0.00844, -0.00701, -0.00701, -0.00728, -0.00735, -0.00647, -0.00573,
c     1          -0.0056, -0.00585, -0.00614, -0.00644, -0.0067, -0.00713, -0.00744,
c     2          -0.00812, -0.00844, -0.00771, -0.00479, -0.00183, -0.00152, -0.00144,
c     3          -0.00137, -0.00136 /
C     Updated 11/16/2020
c      data sigphi / 0.698487036, 0.698567591, 0.707423348, 0.745178656, 0.801901481,
c     1              0.824902609, 0.811419926, 0.79259586, 0.765145342, 0.742104729,
c     2              0.711555585, 0.699247374, 0.707484704, 0.724389736, 0.763758646,
c     3              0.795142376, 0.799800905, 0.796173644, 0.774495359, 0.741211134,
c     4              0.699833142, 0.719869168, 0.638371395 /
c      data sigtau / 0.468873526, 0.469119244, 0.482129367, 0.528220944, 0.56902824,
c     1              0.560503577, 0.487844172, 0.45530486, 0.433521089, 0.400213051,
c     2              0.363497143, 0.337707584, 0.297657575, 0.295706859, 0.307184299,
c     3              0.267014442, 0.302249866, 0.342410131, 0.354900285, 0.356886341,
c     4              0.333790792, 0.484805644, 0.335441147 /


C     Constant parameters
      Vref = 760.0
      f1 = 0.0
      f3 = 0.1

C Find the requested spectral period and corresponding coefficients
      nPer = 23

C First check for the PGA case
      if (specT .eq. 0.0) then
         i1=1
         period2 = period(i1)
         eT = e(i1)
         a1T = a1(i1)
         d1T = d1(i1)
         a2T = a2(i1)
         hT  = h(i1)
         CdT = Cd(i1)
         DdT = Dd(i1)
         CT  = C(i1)
         VcT = Vc(i1)
         f4T = f4(i1)
         f5T = f5(i1)
         sigphiT = sigphi(i1)
         sigtauT = sigtau(i1)

         goto 1011

C First check for the PGV case
      elseif (specT .eq. -1.0) then
         i1=2
         eT = e(i1)
         a1T = a1(i1)
         d1T = d1(i1)
         a2T = a2(i1)
         hT  = h(i1)
         CdT = Cd(i1)
         DdT = Dd(i1)
         CT  = C(i1)
         VcT = Vc(i1)
         f4T = f4(i1)
         f5T = f5(i1)
         sigphiT = sigphi(i1)
         sigtauT = sigtau(i1)

         goto 1011

      endif

C   For other periods, loop over the spectral period range of the attenuation relationship.
      do i=3,nper-1
         if (specT .ge. period(i) .and. specT .le. period(i+1) ) then
            count1 = i
            count2 = i+1
            goto 1020
         endif
      enddo
      write (*,'( i5,2f12.6)') nper, specT, period(nper)

C Selected spectral period is outside range defined by attenuaton model.
      write (*,*)
      write (*,*) 'SMK (2020) Horizontal'
      write (*,*) 'attenuation model is not defined for a '
      write (*,*) ' spectral period of: '
      write (*,'(a10,f10.5)') ' Period = ',specT
      write (*,*) 'This spectral period is outside the defined'
      write (*,*) 'period range in the code or beyond the range'
      write (*,*) 'of spectral periods for interpolation.'
      write (*,*) 'Please check the input file.'
      write (*,*)
      stop 99

C Interpolate the coefficients for the requested spectral period.
 1020 call S24_interp (period(count1),period(count2),sigphi(count1),sigphi(count2),
     +                   specT,sigphiT,iflag)
      call S24_interp (period(count1),period(count2),sigtau(count1),sigtau(count2),
     +                   specT,sigtauT,iflag)
      call S24_interp (period(count1),period(count2),e(count1),e(count2),
     +                   specT,eT,iflag)
      call S24_interp (period(count1),period(count2),a1(count1),a1(count2),
     +                   specT,a1T,iflag)
      call S24_interp (period(count1),period(count2),d1(count1),d1(count2),
     +                   specT,d1T,iflag)
      call S24_interp (period(count1),period(count2),a2(count1),a2(count2),
     +                   specT,a2T,iflag)
      call S24_interp (period(count1),period(count2),h(count1),h(count2),
     +                   specT,hT,iflag)
      call S24_interp (period(count1),period(count2),Cd(count1),Cd(count2),
     +                   specT,CdT,iflag)
      call S24_interp (period(count1),period(count2),Dd(count1),Dd(count2),
     +                   specT,DdT,iflag)
      call S24_interp (period(count1),period(count2),C(count1),C(count2),
     +                   specT,CT,iflag)
      call S24_interp (period(count1),period(count2),Vc(count1),Vc(count2),
     +                   specT,VcT,iflag)
      call S24_interp (period(count1),period(count2),f4(count1),f4(count2),
     +                   specT,f4T,iflag)
      call S24_interp (period(count1),period(count2),f5(count1),f5(count2),
     +                   specT,f5T,iflag)

 1011 period2 = specT

C     Compute the base model including PGARock

C     Path Terms
C     k Function Term
      if (specT .eq. -1.0) then
         kterm = 0.002
      elseif (specT .lt. 0.3) then
         kterm = 0.003
      elseif (specT .le. 0.6) then
         kterm = 0.00126 - 0.00332*log10(specT)
      else
         kterm = 0.002
      endif
      ktermpga = 0.003

c     c function term
      if (specT .eq. -1.0) then
         cterm = 0.0028*10**(0.5*min(mag,8.3))
      elseif (specT .lt. 0.3) then
         cterm = 0.0055*10**(0.5*min(mag,8.3))
      elseif (specT .le. 0.6) then
         cterm = (0.000810-0.00897*log10(specT))*10**(0.5*min(mag,8.3))
      else
         cterm = 0.0028*10**(0.5*min(mag,8.3))
      endif
      ctermpga = 0.0055*10**(0.5*min(mag,8.3))

C     Depth Cut off function
      cutDD = 1.7*depth

C     g Function term
      if (depth .gt. mohodepth .and. Rrup .ge. CutDD) then
         gterm = 0.6*log10(CutDD+Cterm) - 1.6*log10(Rrup+Cterm)
         gtermpga = 0.6*log10(CutDD+Ctermpga) - 1.6*log10(Rrup+Ctermpga)
      else
         gterm = -1.0*log10(Rrup + Cterm)
         gtermpga = -1.0*log10(Rrup + Ctermpga)
      endif

C     Basin terms.
c     gd function term
      gdterm = CdT + DdT*z25
      gdtermpga = Cd(1) + Dd(1)*z25

C     Source Terms
c     b function term
      if (specT .lt. 2.0 ) then
         if (mag .lt. 8.3) then
            bterm = a1T*mag + d1T*ftype + hT*depth + eT
            btermpga = a1(1)*mag + d1(1)*ftype + h(1)*depth + e(1)
         else
            bterm = a1T*mag + (a2T-a1T)*(mag-8.3) + d1T*ftype + hT*depth + eT
            btermpga = a1(1)*mag + (a2(1)-a1(1))*(mag-8.3) + d1(1)*ftype + h(1)*depth + e(1)
         endif
      elseif (specT .ge. 2.0) then
         if (mag .lt. 7.5) then
            bterm = a1T*mag + d1T*ftype + hT*depth + eT
            btermpga = a1(1)*mag + d1(1)*ftype + h(1)*depth + e(1)
         else
            bterm = a1T*mag + (a2T-a1T)*(mag-7.5) + d1T*ftype + hT*depth + eT
            btermpga = a1(1)*mag + (a2(1)-a1(1))*(mag-7.5) + d1(1)*ftype + h(1)*depth + e(1)
         endif
      endif

C     Reference Rock PGA

      pgarock = 10**(btermpga + gtermpga - ktermpga*Rrup + gdtermpga)

C     Site Response Terms
c     gs function term
c     linear term
      if (Vs30 .le. VcT) then
         flin = CT*alog(Vs30/Vref)
      else
         flin = CT*alog(VcT/Vref)
      endif
      f2 = 0.5*f4T*(exp(f5T*(min(Vs30,760.0)-360))-exp(f5T*(760.0-360.0)))
      fnl = f2*alog((pgarock+f3)/f3)
      gsterm = (flin + fnl)/alog(10.0)

C     Median Ground motion values
      median = 10**(bterm + gterm - kterm*Rrup + gsterm + gdterm)

C     Convert to LN units
      lnSa = alog(Median)
c     Convert units spectral acceleration in gal
      lnSa = lnSa + 6.89

C     Set sigma values to return
      sigma = sqrt(sigPhiT*sigPhiT + sigTauT*sigTauT)
      phi = sigPhiT
      tau = sigTauT

      return
      end
