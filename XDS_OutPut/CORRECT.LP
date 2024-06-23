
 ***** CORRECT ***** (VERSION Jun 30, 2023  BUILT=20230630)  23-Jun-2024

 INPUT PARAMETER VALUES
 ----------------------
 SPACE_GROUP_NUMBER=    1 as used in the INTEGRATE step
 UNIT_CELL_CONSTANTS=    58.032    58.011   150.510  90.019  90.017  89.989 as used by INTEGRATE
 FRIEDEL'S_LAW=TRUE 
 PROFILE_FITTING= TRUE 
 OVERLOAD=     65535     MINPK=  75.00000    WFAC1=   1.0
 INCLUDE_RESOLUTION_RANGE=    20.000     0.000
 NAME_TEMPLATE_OF_DATA_FRAMES=/uhome/a01431/XDS_example/thaumatin_9_1_??????.h5 GENERIC
 DATA_RANGE=       1    3600
 ROTATION_AXIS=  0.999983 -0.005774  0.000760
 OSCILLATION_RANGE=  0.150000
 STARTING_ANGLE=     0.000     STARTING_FRAME=       1
 X-RAY_WAVELENGTH=  0.979500
 INCIDENT_BEAM_DIRECTION=  0.000651  0.002481  0.999997
 FRACTION_OF_POLARIZATION=   0.999
 POLARIZATION_PLANE_NORMAL=  0.000000  1.000000  0.000000
 AIR=   0.00032
 DETECTOR=PILATUS
 SILICON=  3.706910 SENSOR_THICKNESS=  0.450000
 NUMBER OF DETECTOR SEGMENTS   1
 NX=  4148  NY=  4362    QX=  0.075000  QY=  0.075000
 ORGX=   2214.85  ORGY=   2301.00
 DETECTOR_DISTANCE=   200.000
 DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
 DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
 BEAM_DIVERGENCE_E.S.D.=   0.019
 REFLECTING_RANGE_E.S.D.=   0.056
 MINIMUM_ZETA=     0.050
 MAXIMUM_ERROR_OF_SPOT_POSITION=   3.0
 MAXIMUM_ERROR_OF_SPINDLE_POSITION=   2.0
 MINIMUM_I/SIGMA=   3.0
 REFLECTIONS/CORRECTION_FACTOR=      50
 SNRC=  0.10 I/sigma for use in scaling factor determination
 STRICT_ABSORPTION_CORRECTION=FALSE 
 CORRECTIONS= DECAY MODULATION ABSORPTION                      
 TEST_RESOLUTION_RANGE=    10.000     5.000     MIN_RFL_Rmeas=    50     MAX_FAC_Rmeas=    2.00
 REFERENCE_DATA_SET=                                                  
 FIT_B-FACTOR_TO_REFERENCE_DATA_SET=FALSE 

 1503618 REFLECTIONS ON FILE "INTEGRATE.HKL"
       0 CORRUPTED REFLECTION RECORDS (IGNORED)
       0 MATCHING REFLECTIONS SPECIFIED IN FILE "FILTER.HKL"
       0 REFLECTIONS INCOMPLETE OR OUTSIDE IMAGE RANGE       1 ...    3600
       0 OVERLOADED REFLECTIONS (IGNORED)
     107 REFLECTIONS OUTSIDE ACCEPTED RESOLUTION RANGES
                OR TOO CLOSE TO ROTATION AXIS (IGNORED)
 1503511 REFLECTIONS ACCEPTED



 ******************************************************************************
                      AUTOMATIC SPACE GROUP ASSIGNMENT
 ******************************************************************************

 XDS adopts the following approach.
 (1) it looks for possible symmetries of the crystal lattice
 (2) it computes a redundancy independent R-factor for all enantiomorphous
     point groups compatible with the observed lattice symmetry.
 (3) it selects the group which explains the intensity data at an acceptable,
     redundancy-independent R-factor (Rmeas, Rrim) using a minimum number of
     unique reflections.

 This approach does not test for the presence of screw axes. Consequently,
 orthorhombic cell axes will be specified in increasing length (following
 conventions), despite the possibility that different assignments for the
 cell axes could become necessary for space groups P222(1) and P2(1)2(1)2
 containing one or two screw axes, respectively.

 The user can always override the automatic decisions by specifying the
 correct space group number and unit cell constants in XDS.INP and repeating
 the CORRECT step of XDS. This provides a simple way to rename orthorhombic
 cell constants if screw axes are present.

 In addition, the user has the option to specify in XDS.INP
 (a) a reference data set or
 (b) a reindexing transformation or
 (c) the three basis vectors (if known from processing a previous data set
 taken at the same crystal orientation in a multi-wavelength experiment).
 These features of XDS are useful for resolving the issue of alternative
 settings of polar or rhombohedral cells (like P4, P6, R3).

 =====> Specifications provided by the user in XDS.INP:
 SPACE_GROUP_NUMBER=    0 (0: unknown symmetry and unit cell)
 UNIT_CELL_CONSTANTS=     0.00     0.00     0.00   0.000   0.000   0.000
 UNIT_CELL_A-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_B-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 UNIT_CELL_C-AXIS=     0.000     0.000     0.000 (0 0 0 : unknown)
 REIDX=  0  0  0  0  0  0  0  0  0  0  0  0 (all 0 : not specified)
 TEST_RESOLUTION_RANGE=   10.00    5.00 (Angstrom) for space group determination
 NUMBER OF ACCEPTED REFLECTIONS FOR CALCULATING Rmeas            26853
 NUMBER OF ACCEPTED UNIQUE REFLECTIONS FROM REFERENCE DATA SET       0


 *********** DETERMINATION OF LATTICE CHARACTER AND BRAVAIS LATTICE ***********

 The CHARACTER OF A LATTICE is defined by the metrical parameters of its
 reduced cell as described in the INTERNATIONAL TABLES FOR CRYSTALLOGRAPHY
 Volume A, p. 746 (KLUWER ACADEMIC PUBLISHERS, DORDRECHT/BOSTON/LONDON, 1989).
 Note that more than one lattice character may have the same BRAVAIS LATTICE.

 A lattice character is marked "*" to indicate a lattice consistent with the
 observed locations of the diffraction spots. These marked lattices must have
 low values for the QUALITY OF FIT and their implicated UNIT CELL CONSTANTS
 should not violate the ideal values by more than
 MAX_CELL_AXIS_ERROR=  0.03
 MAX_CELL_ANGLE_ERROR=           2.0 (Degrees)

 The REINDEXING TRANSFORMATION REIDX() consists of 12 integers that relate
 the original indices H,K,L from file INTEGRATE.HKL to the indices H',K',L'
 with respect to the new cell.
    H'=(REIDX(1)*H+REIDX( 2)*K+REIDX( 3)*L)/IDXV+REIDX( 4)
    K'=(REIDX(5)*H+REIDX( 6)*K+REIDX( 7)*L)/IDXV+REIDX( 8)
    L'=(REIDX(9)*H+REIDX(10)*K+REIDX(11)*L)/IDXV+REIDX(12)
 The value of the integer IDXV depends on the lattice type used for
 specifying reflections on file INTEGRATE.HKL; IDXV is
 1 for a primitive,
 2 for a face or body centred,
 3 for a rhombohedral,
 4 for a lattice centred on all faces.
 IDXV is set by XDS and cannot be input by the user.

  LATTICE-  BRAVAIS-   QUALITY  UNIT CELL CONSTANTS (ANGSTROEM & DEGREES)    REINDEXING TRANSFORMATION
 CHARACTER  LATTICE     OF FIT      a      b      c   alpha  beta gamma

 *  31        aP          0.0      58.0   58.0  150.5  90.0  90.0  90.0    0  1  0  0  1  0  0  0  0  0 -1  0
 *  44        aP          0.1      58.0   58.0  150.5  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
 *  10        mC          0.4      82.1   82.0  150.5  90.0  90.0  90.0    1  1  0  0 -1  1  0  0  0  0  1  0
 *  14        mC          0.5      82.1   82.0  150.5  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  33        mP          0.6      58.0   58.0  150.5  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
 *  35        mP          0.6      58.0   58.0  150.5  90.0  90.0  90.0    1  0  0  0  0  1  0  0  0  0  1  0
 *  34        mP          0.9      58.0  150.5   58.0  90.0  90.0  90.0    0  1  0  0  0  0  1  0  1  0  0  0
 *  32        oP          1.0      58.0   58.0  150.5  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
 *  13        oC          1.3      82.1   82.0  150.5  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
 *  11        tP          1.4      58.0   58.0  150.5  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    37        mC        249.8     306.5   58.0   58.0  90.0  90.0  79.1    0 -1 -2  0  0 -1  0  0 -1  0  0  0
    28        mC        249.9      58.0  306.5   58.0  90.0  90.0  79.1    0 -1  0  0  0 -1 -2  0  1  0  0  0
    41        mC        250.0     306.5   58.0   58.0  90.0  90.0  79.1    1  0  2  0  1  0  0  0  0  1  0  0
    29        mC        250.1      58.0  129.7  150.5  90.0  90.0  63.5    0  1  0  0 -2  1  0  0  0  0  1  0
    36        oC        250.1      58.0  306.5   58.0  90.0  90.0 100.9    0 -1  0  0  0  1  2  0 -1  0  0  0
    30        mC        250.2      58.0  306.5   58.0  90.0  90.0  79.1    1  0  0  0  1  0  2  0  0 -1  0  0
    40        oC        250.4      58.0  306.5   58.0  90.0  90.0 100.9    1  0  0  0 -1  0 -2  0  0  1  0  0
    39        mC        250.6     129.8   58.0  150.5  90.0  90.0  63.4    2  1  0  0  0  1  0  0  0  0  1  0
    38        oC        251.0      58.0  129.8  150.5  90.0  90.0 116.6    0  1  0  0 -2 -1  0  0  0  0  1  0
    12        hP        251.4      58.0   58.0  150.5  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
    27        mC        499.5     129.7   58.0  161.3  90.0 108.8  63.5    2 -1  0  0  0 -1  0  0 -1  0 -1  0
    42        oI        499.6      58.0   58.0  312.0 100.7 100.7  90.0    0  1  0  0  1  0  0  0 -1 -1 -2  0
    15        tI        499.7      58.0   58.0  312.0  79.3  79.3  90.0    0 -1  0  0 -1  0  0  0 -1 -1 -2  0
    17        mC        500.0      82.1   82.0  150.5  90.0  90.0  90.0   -1 -1  0  0 -1  1  0  0  0  0 -1  0
    26        oF        624.1      58.0  129.7  306.5  85.1 100.9 116.5    0  1  0  0  2 -1  0  0  0 -1 -2  0
     9        hR        749.5      58.0   82.0  458.9  90.0  97.2 135.0    0  1  0  0  1 -1  0  0 -1 -1 -3  0
     1        cF        999.0     171.4  171.4  171.4 140.4  57.2 140.4   -1  1 -1  0  1  1  1  0  1 -1 -1  0
     2        hR        999.0      82.0  161.3  171.5 134.2  90.0 104.7   -1  1  0  0  0 -1 -1  0 -1 -1  1  0
     3        cP        999.0      58.0   58.0  150.5  90.0  90.0  90.0    0 -1  0  0 -1  0  0  0  0  0 -1  0
     5        cI        999.0     161.3   82.1  161.3  75.3  29.5  75.3    0 -1 -1  0 -1 -1  0  0 -1  0 -1  0
     4        hR        999.0      82.0  161.3  171.4 134.2  90.0 104.7    1 -1  0  0  0  1 -1  0  1  1  1  0
     6        tI        999.0     161.3  161.3   82.1  75.3  75.3  29.5   -1  0 -1  0  0 -1 -1  0 -1 -1  0  0
     7        tI        999.0     161.3   82.1  161.3  75.3  29.5  75.3    0 -1 -1  0 -1 -1  0  0 -1  0 -1  0
     8        oI        999.0      82.0  150.5  171.4  28.6  90.0  90.0   -1  1  0  0  0  0 -1  0 -1 -1 -1  0
    16        oF        999.0      82.0   82.1  312.0 105.2  90.0  90.0   -1  1  0  0 -1 -1  0  0  1  1  2  0
    18        tI        999.0     161.3  171.4   58.0  70.2  90.0 134.2   -1  0 -1  0 -1  1  1  0  0  1  0  0
    19        oI        999.0      58.0  161.3  171.4  45.8  70.2  90.0    0 -1  0  0 -1  0 -1  0  1 -1 -1  0
    20        mC        999.0     161.3  161.3   58.0  90.0  90.0 137.8    1  0 -1  0  1  0  1  0  0 -1  0  0
    21        tP        999.0      58.0  150.5   58.0  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
    22        hP        999.0      58.0  150.5   58.0  90.0  90.0  90.0   -1  0  0  0  0  0 -1  0  0 -1  0  0
    23        oC        999.0     161.3  161.3   58.0  90.0  90.0  42.2   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    24        hR        999.0     198.7  161.3   58.0  90.0  73.0  60.2   -2 -1 -1  0  1  0 -1  0  0 -1  0  0
    25        mC        999.0     161.3  161.3   58.0  90.0  90.0  42.2   -1  0 -1  0  1  0 -1  0  0 -1  0  0
    43        mI        999.0      82.0  312.0   58.0 100.7 135.0  90.0    1 -1  0  0  1  1  2  0 -1  0  0  0

 ****** LATTICE SYMMETRY IMPLICATED BY SPACE GROUP SYMMETRY ******

BRAVAIS-            POSSIBLE SPACE-GROUPS FOR PROTEIN CRYSTALS
  TYPE                     [SPACE GROUP NUMBER,SYMBOL]
  aP      [1,P1]
  mP      [3,P2] [4,P2(1)]
 mC,mI    [5,C2]
  oP      [16,P222] [17,P222(1)] [18,P2(1)2(1)2] [19,P2(1)2(1)2(1)]
  oC      [21,C222] [20,C222(1)]
  oF      [22,F222]
  oI      [23,I222] [24,I2(1)2(1)2(1)]
  tP      [75,P4] [76,P4(1)] [77,P4(2)] [78,P4(3)] [89,P422] [90,P42(1)2]
          [91,P4(1)22] [92,P4(1)2(1)2] [93,P4(2)22] [94,P4(2)2(1)2]
          [95,P4(3)22] [96,P4(3)2(1)2]
  tI      [79,I4] [80,I4(1)] [97,I422] [98,I4(1)22]
  hP      [143,P3] [144,P3(1)] [145,P3(2)] [149,P312] [150,P321] [151,P3(1)12]
          [152,P3(1)21] [153,P3(2)12] [154,P3(2)21] [168,P6] [169,P6(1)]
          [170,P6(5)] [171,P6(2)] [172,P6(4)] [173,P6(3)] [177,P622]
          [178,P6(1)22] [179,P6(5)22] [180,P6(2)22] [181,P6(4)22] [182,P6(3)22]
  hR      [146,R3] [155,R32]
  cP      [195,P23] [198,P2(1)3] [207,P432] [208,P4(2)32] [212,P4(3)32]
          [213,P4(1)32]
  cF      [196,F23] [209,F432] [210,F4(1)32]
  cI      [197,I23] [199,I2(1)3] [211,I432] [214,I4(1)32]


 ********************* SYMMETRY OF REFLECTION INTENSITIES *********************

 Analysis of symmetry related reflections is based on strong data within a
 TEST_RESOLUTION_RANGE=   10.00    5.00 Angstrom (defined in XDS.INP).

 The symmetry test is carried out for the enantiomorphous space groups above
 (omitting those with screw axes) that conform with the observed lattice
 symmetry. The test is quantified by the redundancy independent R-factor,
 Rmeas, defined by Diederichs & Karplus, Nature Struct. Biol. 4, 269-275 (1997).

 For computation of Rmeas a minimum number of symmetry-related reflections is
 required and set by the user (in XDS.INP) as the parameter
                      MIN_RFL_Rmeas=   50

 Low values for Rmeas should be obtained if the crystal possesses the symmetry
 elements of the tested space group. The decision whether the crystal's data
 satisfy the requirements of the tested space group is based on a constant set
 by the user (in XDS.INP) as the parameter
                      MAX_FAC_Rmeas=      2.
 which provides an upper limit for acceptable Rmeas as a multiple of the
 lowest Rmeas found for all tested space groups. If Rmeas is larger than this
 upper limit, the corresponding space group symmetry is considered violated
 by the observed data.

 XDS will select (*) that space group that requires a minimum number of unique
 reflections to account for the observed intensities at an acceptable Rmeas.

 UNIQUE   = number of unique reflections
 Rmeas    = redundancy independent R-factor (intensities)
 COMPARED = number of reflections used for calculating Rmeas

 SPACE-GROUP         UNIT CELL CONSTANTS            UNIQUE   Rmeas  COMPARED  LATTICE-
   NUMBER      a      b      c   alpha beta gamma                            CHARACTER

       5      82.1   82.0  150.5  90.0  90.0  90.0    3891    15.2    22962    10 mC
      75      58.0   58.0  150.5  90.0  90.0  90.0    1923    16.0    24930    11 tP
      89      58.0   58.0  150.5  90.0  90.0  90.0    1141    16.0    25712    11 tP
      21      82.1   82.0  150.5  90.0  90.0  90.0    2063    15.8    24790    13 oC
       5      82.1   82.0  150.5  90.0  90.0  90.0    3891    15.2    22962    14 mC
       1      58.0   58.0  150.5  90.0  90.0  90.0    7403     6.9    19450    31 aP
      16      58.0   58.0  150.5  90.0  90.0  90.0    2127    14.8    24726    32 oP
       3      58.0   58.0  150.5  90.0  90.0  90.0    3955     9.6    22898    35 mP
   *   3      58.0   58.0  150.5  90.0  90.0  90.0    3953    13.3    22900    33 mP
       3      58.0  150.5   58.0  90.0  90.0  90.0    3831    14.6    23022    34 mP
       1      58.0   58.0  150.5  90.0  90.0  90.0    7403     6.9    19450    44 aP
 

 ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************

 SPACE_GROUP_NUMBER=    3
 UNIT_CELL_CONSTANTS=    58.01    58.03   150.51  90.000  90.019  90.000
 UNIT_CELL_A-AXIS=   -50.048    18.083   -23.098
 UNIT_CELL_B-AXIS=   -15.089   -55.056   -10.435
 UNIT_CELL_C-AXIS=   -65.236    -7.739   135.417
 REFLECTIONS H,K,L FROM "INTEGRATE.HKL" WILL BE REINDEXED BY
 H' = (  0*H + -1*K +  0*L)/ 1 +  0
 K' = ( -1*H +  0*K +  0*L)/ 1 +  0
 L' = (  0*H +  0*K + -1*L)/ 1 +  0

 *****  2 EQUIVALENT POSITIONS IN SPACE GROUP #  3 *****

    If x',y',z' is an equivalent position to x,y,z, then
        x'=x*ML(1)+y*ML( 2)+z*ML( 3)+ML( 4)/12.0
        y'=x*ML(5)+y*ML( 6)+z*ML( 7)+ML( 8)/12.0
        z'=x*ML(9)+y*ML(10)+z*ML(11)+ML(12)/12.0

    #    1  2  3  4    5  6  7  8    9 10 11 12 
    1    1  0  0  0    0  1  0  0    0  0  1  0
    2   -1  0  0  0    0  1  0  0    0  0 -1  0

 ******************************************************************************
      MEAN DISCREPANCIES BETWEEN OBSERVED AND CALCULATED SPOT LOCATIONS
 ******************************************************************************

 The discrepancies in X- and Y-coordinates of the spots are depicted in the
 two images DX-CORRECTIONS.cbf and DY-CORRECTIONS.cbf for inspection with
 the XDS-Viewer.



 ******************************************************************************
  REFINEMENT OF DIFFRACTION PARAMETERS USING ALL IMAGES
 ******************************************************************************


 REFINED VALUES OF DIFFRACTION PARAMETERS DERIVED FROM    210551 INDEXED SPOTS
 REFINED PARAMETERS:   POSITION BEAM AXIS ORIENTATION CELL
 STANDARD DEVIATION OF SPOT    POSITION (PIXELS)     0.44
 STANDARD DEVIATION OF SPINDLE POSITION (DEGREES)    0.04
 SPACE GROUP NUMBER      3
 UNIT CELL PARAMETERS     57.947    57.953   150.296  90.000  90.000  90.000
 E.S.D. OF CELL PARAMETERS  6.5E-03 9.8E-03 1.6E-02 0.0E+00 4.9E-03 0.0E+00
 REC. CELL PARAMETERS   0.017257  0.017255  0.006654  90.000  90.000  90.000
 COORDINATES OF UNIT CELL A-AXIS   -50.164    18.072   -22.691
 COORDINATES OF UNIT CELL B-AXIS   -15.149   -54.982   -10.299
 COORDINATES OF UNIT CELL C-AXIS   -64.165    -7.737   135.690
 CRYSTAL MOSAICITY (DEGREES)     0.056
 LAB COORDINATES OF ROTATION AXIS  0.999959 -0.005812 -0.006879
 DIRECT BEAM COORDINATES (REC. ANGSTROEM)   0.008512  0.009725  1.020847
 DETECTOR COORDINATES (PIXELS) OF DIRECT BEAM    2215.91   2306.79
 DETECTOR ORIGIN (PIXELS) AT                     2193.71   2281.42
 CRYSTAL TO DETECTOR DISTANCE (mm)       199.71
 LAB COORDINATES OF DETECTOR X-AXIS  1.000000  0.000000  0.000000
 LAB COORDINATES OF DETECTOR Y-AXIS  0.000000  1.000000  0.000000


 THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:
 SPACE_GROUP_NUMBER=    3
 UNIT_CELL_CONSTANTS=    57.95    57.95   150.30  90.000  90.000  90.000


 ******************************************************************************
        MEAN INTENSITY AS FUNCTION OF SPINDLE POSITION WITHIN DATA IMAGE
 ******************************************************************************

 This statistics could serve as a diagnostic tool for detecting
 shutter problems (suggested by Kay Diederichs).
 Data are corrected for this effect if PATCH_SHUTTER_PROBLEM=TRUE in XDS.INP.
 ===> Selected : PATCH_SHUTTER_PROBLEM=FALSE

 INTERVAL = Angular interval in units of the oscillation range
            covered by a data image.
 NUMBER   = Number of reflections in each interval
 INTENSITY= Mean intensity of the reflections in each interval
 FACTOR   = Correction factor applied to intensities;
            INTENSITY(overall)/INTENSITY(interval)

   INTERVAL     NUMBER      INTENSITY      FACTOR
 -0.05  0.05    150703          0.558       1.002
  0.05  0.15    149690          0.564       0.992
  0.15  0.25    150249          0.558       1.001
  0.25  0.35    151084          0.558       1.002
  0.35  0.45    150236          0.560       0.998
  0.45  0.55    150679          0.559       1.000
  0.55  0.65    150173          0.558       1.002
  0.65  0.75    149821          0.556       1.005
  0.75  0.85    150625          0.558       1.003
  0.85  0.95    150251          0.562       0.996
 
 -0.05  0.95   1503511          0.559       1.000
 

 ******************************************************************************
           CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & RESOLUTION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED     1440
 DEGREES OF FREEDOM OF CHI^2 FIT            940991.9
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.917
 NUMBER OF CYCLES CARRIED OUT                      5

 CORRECTION FACTORS for visual inspection by XDS-Viewer DECAY.cbf           
 XMIN=     0.2 XMAX=  2399.8 NXBIN=   72
 YMIN= 0.00256 YMAX= 0.75107 NYBIN=   20
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS    1162448


 ******************************************************************************
  CORRECTION FACTORS AS FUNCTION OF X (fast) & Y(slow) IN THE DETECTOR PLANE
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED    23244
 DEGREES OF FREEDOM OF CHI^2 FIT            940359.1
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.909
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer MODPIX.cbf          
 XMIN=     6.7 XMAX=  4141.7 NXBIN=  149
 YMIN=     6.7 YMAX=  4355.6 NYBIN=  156
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS    1162448


 ******************************************************************************
   CORRECTION FACTORS AS FUNCTION OF IMAGE NUMBER & DETECTOR SURFACE POSITION
 ******************************************************************************

 RECIPROCAL CORRECTION FACTORS FOR INPUT DATA SETS MERGED TO
 OUTPUT FILE: XDS_ASCII.HKL                                     

 THE CALCULATIONS ASSUME         FRIEDEL'S_LAW= TRUE
 TOTAL NUMBER OF CORRECTION FACTORS DEFINED      936
 DEGREES OF FREEDOM OF CHI^2 FIT            940996.0
 CHI^2-VALUE OF FIT OF CORRECTION FACTORS      0.903
 NUMBER OF CYCLES CARRIED OUT                      3

 CORRECTION FACTORS for visual inspection by XDS-Viewer ABSORP.cbf          
 XMIN=     0.2 XMAX=  2399.8 NXBIN=   72
 DETECTOR_SURFACE_POSITION=    2074    2181
 DETECTOR_SURFACE_POSITION=    2776    2920
 DETECTOR_SURFACE_POSITION=    1372    2920
 DETECTOR_SURFACE_POSITION=    1372    1443
 DETECTOR_SURFACE_POSITION=    2776    1443
 DETECTOR_SURFACE_POSITION=    3664    2874
 DETECTOR_SURFACE_POSITION=    2733    3853
 DETECTOR_SURFACE_POSITION=    1416    3853
 DETECTOR_SURFACE_POSITION=     485    2874
 DETECTOR_SURFACE_POSITION=     485    1489
 DETECTOR_SURFACE_POSITION=    1416     510
 DETECTOR_SURFACE_POSITION=    2733     510
 DETECTOR_SURFACE_POSITION=    3664    1489
 NUMBER OF REFLECTIONS USED FOR DETERMINING CORRECTION FACTORS    1162448


 ******************************************************************************
    CORRECTION PARAMETERS FOR THE STANDARD ERROR OF REFLECTION INTENSITIES
 ******************************************************************************

 The variance v0(I) of the intensity I obtained from counting statistics is
 replaced by v(I)=a*(v0(I)+b*I^2). The model parameters a, b are chosen to
 minimize the discrepancies between v(I) and the variance estimated from
 sample statistics of symmetry related reflections. This model implicates
 an asymptotic limit ISa=1/SQRT(a*b) for the highest I/Sigma(I) that the
 experimental setup can produce (Diederichs (2010) Acta Cryst D66, 733-740).

     a        b          ISa
 6.815E-01  1.553E-03   30.74


 ******************************************************************************
      STANDARD ERROR OF REFLECTION INTENSITIES AS FUNCTION OF RESOLUTION
      FOR DATA SET  XDS_ASCII.HKL                                     
 ******************************************************************************

  I/Sigma  = mean intensity/Sigma of a reflection in shell
  Chi^2    = goodness of fit between sample variances of
             symmetry-related intensities and their errors
             (Chi^2 = 1 for perfect agreement)
  R-FACTOR
  observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
  expected = expected R-FACTOR derived from Sigma(I)

   NUMBER  = number of reflections in resolution shell
             used for calculation of R-FACTOR
  ACCEPTED = number of accepted reflections
  REJECTED = number of rejected reflections (MISFITS),
             recognized by comparison with symmetry-related
             reflections.

 RESOLUTION RANGE  I/Sigma  Chi^2  R-FACTOR  R-FACTOR  NUMBER ACCEPTED REJECTED
                                   observed  expected

   19.767   9.978    16.66   1.08      3.19      3.31    3069    3073      25
    9.978   7.553    13.97   1.17      3.85      3.83    5241    5241      34
    7.553   6.323    11.77   1.23      4.94      4.61    6379    6379      58
    6.323   5.547    12.88   1.19      4.81      4.49    6577    6578      57
    5.547   5.001    13.68   1.20      4.50      4.22    8229    8229      64
    5.001   4.590    14.75   1.18      4.05      3.96    9029    9029      84
    4.590   4.266    15.15   1.18      4.03      3.89   10033   10033      75
    4.266   4.002    14.51   1.20      4.32      4.11   10821   10821     114
    4.002   3.782    13.95   1.23      4.46      4.15    9997    9997     120
    3.782   3.594    13.57   1.24      4.55      4.31   10950   10950     137
    3.594   3.432    12.75   1.25      4.92      4.50   11995   11995     161
    3.432   3.290    11.82   1.27      5.48      4.90   13077   13077     147
    3.290   3.164    10.96   1.26      6.09      5.43   13559   13559     147
    3.164   3.052    10.20   1.25      6.32      5.69   13853   13853     134
    3.052   2.951     9.46   1.27      6.98      6.18   14843   14843     140
    2.951   2.859     8.31   1.27      8.09      7.07   15184   15184     137
    2.859   2.776     7.45   1.25      8.62      7.54   13976   13977     136
    2.776   2.699     7.39   1.25      9.11      8.03   14902   14906     138
    2.699   2.628     6.52   1.23     10.08      8.85   15246   15246     134
    2.628   2.563     6.41   1.25     10.44      9.13   16693   16693     156
    2.563   2.502     6.06   1.24     11.30      9.87   16871   16872     137
    2.502   2.445     5.72   1.22     11.74     10.29   17138   17139     152
    2.445   2.392     5.77   1.21     11.57     10.19   18129   18129     168
    2.392   2.343     5.40   1.22     12.13     10.65   18003   18003     144
    2.343   2.296     5.15   1.20     12.55     11.14   18553   18553     146
    2.296   2.252     4.95   1.22     13.59     11.95   17134   17135     123
    2.252   2.210     5.01   1.18     13.24     11.92   17747   17754     116
    2.210   2.171     4.73   1.19     14.21     12.61   19069   19069     141
    2.171   2.134     4.33   1.17     14.88     13.21   19285   19286     148
    2.134   2.098     4.28   1.17     15.16     13.65   20048   20049     120
    2.098   2.065     3.98   1.15     16.25     14.60   20051   20058     114
    2.065   2.032     3.56   1.12     17.81     16.20   21288   21288     132
    2.032   2.002     3.55   1.14     18.13     16.50   20835   20835     126
    2.002   1.972     3.56   1.13     17.97     16.46   22042   22042     116
    1.972   1.944     3.14   1.12     20.92     19.21   21758   21758     143
    1.944   1.917     2.90   1.08     22.21     20.75   21450   21452     116
    1.917   1.891     2.84   1.10     22.44     20.84   20749   20754     108
    1.891   1.867     2.45   1.05     24.84     23.45   20754   20755      95
    1.867   1.843     2.28   1.03     27.14     26.02   23094   23094      89
    1.843   1.820     2.20   1.03     27.70     26.17   22010   22010      89
    1.820   1.798     2.06   1.02     30.01     28.94   23352   23354      84
    1.798   1.776     1.99   1.03     31.10     29.90   23497   23500      78
    1.776   1.756     1.84   1.00     33.51     32.26   23570   23574      74
    1.756   1.736     1.72   0.99     35.09     34.14   25146   25147      94
    1.736   1.717     1.68   0.98     35.62     34.49   24476   24477      84
    1.717   1.698     1.52   0.98     39.70     38.57   24707   24708      70
    1.698   1.680     1.42   0.97     42.02     41.12   23719   23721      56
    1.680   1.662     1.40   0.96     41.72     40.98   23282   23283      75
    1.662   1.645     1.34   0.95     43.98     43.78   24235   24235      78
    1.645   1.629     1.21   0.94     47.78     47.06   25571   25572      55
    1.629   1.613     1.21   0.94     47.70     47.74   24970   24971      53
    1.613   1.598     1.11   0.90     49.91     50.07   26382   26383      62
    1.598   1.583     1.05   0.93     54.04     53.89   26137   26139      43
    1.583   1.568     0.97   0.90     58.15     58.31   25722   25729      60
    1.568   1.554     0.95   0.89     57.41     58.33   27126   27132      46
    1.554   1.540     0.91   0.90     60.40     61.46   26561   26573      49
    1.540   1.526     0.82   0.87     66.67     68.58   26242   26248      34
    1.526   1.513     0.81   0.86     66.20     67.78   25522   25543      27
    1.513   1.500     0.75   0.84     68.74     71.74   24283   24310      32
    1.500   1.488     0.72   0.83     72.05     75.11   22941   22968      29
    1.488   1.476     0.68   0.81     74.94     79.54   23150   23181      20
    1.476   1.464     0.61   0.82     81.49     86.52   22395   22434      17
    1.464   1.452     0.61   0.80     81.59     88.40   22289   22325      12
    1.452   1.441     0.60   0.80     80.75     87.79   21325   21367      16
    1.441   1.430     0.50   0.77     93.60    102.77   21727   21761      10
    1.430   1.419     0.50   0.77     92.63    100.75   20472   20519      22
    1.419   1.409     0.46   0.76     97.60    107.46   20545   20581       5
    1.409   1.398     0.44   0.76    101.42    110.76   20143   20198      14
    1.398   1.388     0.38   0.72    109.18    124.94   18215   18290       9
    1.388   1.378     0.42   0.73    101.53    115.19   16467   16585      12
    1.378   1.368     0.38   0.73    109.27    123.52   15418   15578       5
    1.368   1.359     0.36   0.72    110.97    126.35   14004   14194       8
    1.359   1.350     0.35   0.72    114.38    128.69   12924   13083       9
    1.350   1.341     0.29   0.69    126.38    146.73   12341   12572       8
    1.341   1.332     0.29   0.69    127.20    149.84   11309   11540       4
    1.332   1.323     0.27   0.71    131.19    148.34   10692   10888       1
    1.323   1.314     0.24   0.69    139.07    159.04    9992   10187       0
    1.314   1.306     0.22   0.66    143.31    170.65    9012    9179       0
    1.306   1.298     0.21   0.66    148.66    180.81    8473    8635       3
    1.298   1.290     0.22   0.66    144.56    170.36    7979    8183       2
    1.290   1.282     0.19   0.63    152.91    186.91    7075    7303       0
    1.282   1.274     0.20   0.63    146.18    179.04    6481    6793       3
    1.274   1.266     0.16   0.65    165.52    197.96    6047    6265       1
    1.266   1.259     0.18   0.66    152.60    186.32    5248    5492       2
    1.259   1.251     0.14   0.61    163.47    207.43    4878    5133       0
    1.251   1.244     0.15   0.62    158.15    196.79    4207    4434       0
    1.244   1.237     0.15   0.61    153.99    192.55    3653    3903       0
    1.237   1.230     0.12   0.64    167.04    201.75    3183    3537       0
    1.230   1.223     0.11   0.59    174.88    222.72    2608    3013       0
    1.223   1.216     0.11   0.62    177.56    226.18    2307    2800       2
    1.216   1.209     0.09   0.57    171.81    225.37    1636    2238       0
    1.209   1.203     0.10   0.59    160.04    202.45    1270    1973       0
    1.203   1.196     0.09   0.57    169.12    214.67     720    1520       2
    1.196   1.190     0.10   0.62    130.27    158.22     188    1089       2
    1.190   1.184     0.08   0.76     93.70    101.11      56     933       0
    1.184   1.178     0.10   0.35     53.70     87.48      40     739       0
    1.178   1.172     0.04   0.95    110.50    102.74      32     655       0
    1.172   1.166     0.04   0.89    134.31    126.02      12     486       0
    1.166   1.160     0.02   0.81     66.70     77.37       8     345       0
    1.160   1.154     0.08   0.79     90.97     99.35       6     222       0
   --------------------------------------------------------------------------
   19.767   1.154     3.19   0.99     12.59     12.06 1485629 1497418    6093
 


 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR VARIOUS SUBSETS OF INCLUDED DATA IMAGES
 ******************************************************************************

 Data set statistics is reported below several times, each with a different
 upper limit on the number of images included. This provides the user with
 the information for deciding which data images should be excluded from the
 final data set because of radiation damage or other defects. If the user
 decides for a subset of "good" images that differs from the specification
 as given by the input parameter (see DATA_RANGE= in XDS.INP) the CORRECT-
 step of XDS must be repeated with the new parameter values for DATA_RANGE=.


 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1     401
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     3.42       15579   10243     13827       74.1%       3.5%      3.4%    10442   16.05      4.9%    99.6*     0    0.000       1
     2.43       28107   19803     24044       82.4%       6.4%      5.8%    16280    8.72      9.0%    98.8*     0    0.000       0
     1.99       35578   26082     30720       84.9%      12.2%     11.0%    18729    4.78     17.1%    95.5*     0    0.000       0
     1.73       41876   31052     36307       85.5%      21.1%     20.6%    21405    2.50     29.8%    88.4*     0    0.000       0
     1.55       47000   35541     41066       86.5%      38.5%     39.7%    22711    1.28     54.2%    67.4*     0    0.000       0
     1.41       42593   32511     45168       72.0%      63.3%     72.0%    20025    0.64     89.3%    43.3*     0    0.000       0
     1.31       26334   21752     49241       44.2%     101.4%    119.0%     9133    0.29    143.2%    16.2*     0    0.000       0
     1.22       10742   10213     52713       19.4%     124.0%    155.5%     1058    0.13    175.4%    12.7      0    0.000       0
     1.15        2203    2192     56057        3.9%     245.5%    329.0%       22    0.01    347.2%   -20.5      0    0.000       0
    total      250012  189389    349143       54.2%       8.8%      8.5%   119805    3.24     12.3%    99.4*     0    0.000       1


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES    251148
 NUMBER OF REJECTED MISFITS                              1136
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    0
 NUMBER OF ACCEPTED OBSERVATIONS                       250012
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                 189389



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1     801
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     3.42       31282   12802     13827       92.6%       4.0%      3.7%    28376   20.48      4.9%    99.7*    -1    0.879    2964
     2.43       56091   23340     24044       97.1%       7.1%      6.3%    49890   11.66      8.8%    99.2*     1    0.916    4801
     1.99       71177   29644     30720       96.5%      12.5%     11.2%    63465    6.65     15.4%    97.4*    -1    0.900    5280
     1.73       83744   35568     36307       98.0%      21.9%     20.7%    74376    3.44     27.1%    92.9*     0    0.864    5583
     1.55       93821   39618     41066       96.5%      38.6%     38.0%    83682    1.73     47.8%    81.0*    -2    0.812    5634
     1.41       85264   41352     45168       91.6%      59.3%     63.1%    72495    0.83     76.2%    59.2*    -4    0.744    3228
     1.31       52558   32546     49241       66.1%      87.8%     99.1%    35041    0.36    116.6%    33.0*     2    0.705     912
     1.22       21470   17827     52713       33.8%     124.5%    148.4%     6625    0.16    168.2%    13.3*    -3    0.716      75
     1.15        4437    4350     56057        7.8%     198.8%    286.7%      174    0.09    281.1%   -14.2      0    0.000       0
    total      499844  237047    349143       67.9%      10.6%     10.1%   414124    4.10     13.3%    99.5*    -1    0.852   28477


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES    502022
 NUMBER OF REJECTED MISFITS                              2178
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    0
 NUMBER OF ACCEPTED OBSERVATIONS                       499844
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                 237047



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1    1201
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     3.42       46869   13733     13827       99.3%       4.1%      3.8%    46560   25.29      4.8%    99.8*    -3    0.876    6625
     2.43       84136   24012     24044       99.9%       7.3%      6.4%    83736   14.86      8.6%    99.4*    -1    0.911   11283
     1.99      106600   30707     30720      100.0%      13.0%     11.7%   106054    8.35     15.4%    97.9*    -1    0.897   13258
     1.73      125506   36219     36307       99.8%      23.2%     22.1%   124883    4.31     27.4%    94.1*     0    0.855   14940
     1.55      140663   40698     41066       99.1%      41.9%     41.8%   139938    2.12     49.6%    83.9*    -1    0.811   16081
     1.41      128009   42600     45168       94.3%      66.3%     71.3%   126604    0.98     80.8%    61.0*    -2    0.757    8760
     1.31       78762   37775     49241       76.7%      97.8%    112.4%    67449    0.39    125.8%    31.0*    -1    0.720    1872
     1.22       32207   22878     52713       43.4%     132.1%    161.3%    17429    0.15    180.9%    14.2*   -13    0.662     117
     1.15        6651    6112     56057       10.9%     114.1%    144.3%     1078    0.10    161.4%    18.5*     0    0.000       0
    total      749403  254734    349143       73.0%      11.2%     10.7%   713731    4.96     13.4%    99.6*    -1    0.848   72936


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES    752526
 NUMBER OF REJECTED MISFITS                              3123
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    0
 NUMBER OF ACCEPTED OBSERVATIONS                       749403
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                 254734



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1    1601
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     3.42       62445   13734     13827       99.3%       4.4%      4.0%    62206   28.69      4.9%    99.8*    -3    0.877   10412
     2.43      112189   24035     24044      100.0%       7.9%      6.9%   111967   16.70      8.9%    99.5*    -1    0.915   19132
     1.99      142119   30709     30720      100.0%      14.1%     12.7%   141819    9.31     16.0%    98.3*     0    0.898   23396
     1.73      167325   36227     36307       99.8%      25.3%     24.1%   167014    4.78     28.6%    95.1*     0    0.849   27896
     1.55      187526   40776     41066       99.3%      45.6%     45.7%   187136    2.32     51.5%    86.6*    -1    0.795   29655
     1.41      170518   43859     45168       97.1%      73.4%     78.9%   168486    1.04     84.8%    65.4*    -1    0.731   23375
     1.31      105022   39548     49241       80.3%     111.2%    127.9%    95863    0.38    135.0%    31.2*     0    0.688    7702
     1.22       42909   24805     52713       47.1%     151.7%    186.1%    31176    0.13    199.4%    11.3*    -3    0.643     587
     1.15        8871    7584     56057       13.5%     152.6%    189.5%     2499    0.07    212.7%    13.8*     0    0.000       0
    total      998924  261277    349143       74.8%      12.2%     11.7%   968166    5.41     13.9%    99.7*     0    0.828  142155


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES   1003093
 NUMBER OF REJECTED MISFITS                              4169
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    0
 NUMBER OF ACCEPTED OBSERVATIONS                       998924
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                 261277



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1    2000
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     3.42       78139   13735     13827       99.3%       4.4%      4.2%    78061   32.12      4.9%    99.9*    -2    0.881   11413
     2.43      140111   24038     24044      100.0%       8.0%      7.0%   140023   18.70      8.8%    99.6*    -1    0.917   20920
     1.99      177565   30710     30720      100.0%      14.4%     12.8%   177409   10.45     15.8%    98.7*    -1    0.900   25969
     1.73      209047   36241     36307       99.8%      25.7%     24.3%   208937    5.34     28.3%    96.2*     0    0.857   31277
     1.55      234161   40838     41066       99.4%      46.6%     46.2%   233990    2.56     51.2%    88.8*    -1    0.804   33968
     1.41      213021   44571     45168       98.7%      74.1%     78.5%   211773    1.14     83.0%    70.4*    -1    0.741   30625
     1.31      131209   40292     49241       81.8%     109.7%    124.5%   125528    0.40    129.2%    35.1*     0    0.694   13624
     1.22       53628   25332     52713       48.1%     145.1%    177.5%    45854    0.12    185.3%    16.7*     0    0.639    2327
     1.15       11055    8628     56057       15.4%     153.7%    191.1%     4564    0.07    211.7%    12.6*     1    0.499      31
    total     1247936  264385    349143       75.7%      12.5%     12.0%  1226139    5.98     13.9%    99.7*    -1    0.825  170154


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES   1253158
 NUMBER OF REJECTED MISFITS                              5222
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    0
 NUMBER OF ACCEPTED OBSERVATIONS                      1247936
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                 264385



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1    2400
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     3.42       93737   13737     13827       99.3%       4.4%      4.2%    93732   35.80      4.8%    99.9*    -2    0.881   12441
     2.43      168181   24038     24044      100.0%       8.0%      7.0%   168174   20.93      8.6%    99.7*    -2    0.911   22552
     1.99      212993   30711     30720      100.0%      14.4%     12.8%   212976   11.71     15.5%    99.0*    -1    0.899   28978
     1.73      250880   36243     36307       99.8%      25.8%     24.5%   250862    5.99     27.9%    97.0*    -1    0.857   34441
     1.55      280908   40841     41066       99.5%      47.1%     46.8%   280879    2.87     50.9%    90.9*    -1    0.806   38726
     1.41      255658   44649     45168       98.9%      75.6%     80.1%   255321    1.27     83.1%    74.4*     0    0.737   39602
     1.31      157417   40752     49241       82.8%     113.3%    128.5%   155655    0.44    131.1%    35.8*     0    0.675   22279
     1.22       64370   25904     52713       49.1%     153.6%    187.2%    61524    0.13    195.2%    13.8*     0    0.624    4608
     1.15       13274    9841     56057       17.6%     165.8%    211.6%     6506    0.06    228.7%    11.6*    23    0.667      44
    total     1497418  266716    349143       76.4%      12.6%     12.1%  1485629    6.62     13.7%    99.8*    -1    0.812  203671


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES   1503511
 NUMBER OF REJECTED MISFITS                              6093
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    0
 NUMBER OF ACCEPTED OBSERVATIONS                      1497418
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                 266716



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1    2800
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     3.42       93737   13737     13827       99.3%       4.4%      4.2%    93732   35.80      4.8%    99.9*    -2    0.882   12441
     2.43      168181   24038     24044      100.0%       8.0%      7.0%   168174   20.93      8.6%    99.7*    -1    0.915   22552
     1.99      212993   30711     30720      100.0%      14.4%     12.8%   212976   11.71     15.5%    99.0*    -1    0.898   28978
     1.73      250880   36243     36307       99.8%      25.8%     24.5%   250862    5.99     27.9%    97.0*    -1    0.857   34441
     1.55      280908   40841     41066       99.5%      47.1%     46.8%   280879    2.87     50.9%    90.9*    -1    0.809   38726
     1.41      255658   44649     45168       98.9%      75.6%     80.1%   255321    1.27     83.1%    74.7*     0    0.736   39602
     1.31      157417   40752     49241       82.8%     113.3%    128.5%   155655    0.44    131.1%    35.6*     0    0.679   22279
     1.22       64370   25904     52713       49.1%     153.6%    187.2%    61524    0.13    195.2%    13.6*     0    0.629    4608
     1.15       13274    9841     56057       17.6%     165.8%    211.6%     6506    0.06    228.7%    11.9*    23    0.667      44
    total     1497418  266716    349143       76.4%      12.6%     12.1%  1485629    6.62     13.7%    99.8*    -1    0.814  203671


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES   1503511
 NUMBER OF REJECTED MISFITS                              6093
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    0
 NUMBER OF ACCEPTED OBSERVATIONS                      1497418
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                 266716



 ******************************************************************************
  SUMMARY OF DATA SET STATISTICS FOR IMAGE   DATA_RANGE=       1    3200
 ******************************************************************************


          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     3.42       93737   13737     13827       99.3%       4.4%      4.2%    93732   35.80      4.8%    99.9*    -2    0.877   12441
     2.43      168181   24038     24044      100.0%       8.0%      7.0%   168174   20.93      8.6%    99.7*    -2    0.913   22552
     1.99      212993   30711     30720      100.0%      14.4%     12.8%   212976   11.71     15.5%    99.0*    -1    0.898   28978
     1.73      250880   36243     36307       99.8%      25.8%     24.5%   250862    5.99     27.9%    97.0*    -1    0.857   34441
     1.55      280908   40841     41066       99.5%      47.1%     46.8%   280879    2.87     50.9%    91.0*    -1    0.806   38726
     1.41      255658   44649     45168       98.9%      75.6%     80.1%   255321    1.27     83.1%    74.6*    -1    0.737   39602
     1.31      157417   40752     49241       82.8%     113.3%    128.5%   155655    0.44    131.1%    35.8*     0    0.676   22279
     1.22       64370   25904     52713       49.1%     153.6%    187.2%    61524    0.13    195.2%    13.7*     0    0.626    4608
     1.15       13274    9841     56057       17.6%     165.8%    211.6%     6506    0.06    228.7%    11.6*    25    0.667      44
    total     1497418  266716    349143       76.4%      12.6%     12.1%  1485629    6.62     13.7%    99.8*    -1    0.812  203671


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES   1503511
 NUMBER OF REJECTED MISFITS                              6093
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    0
 NUMBER OF ACCEPTED OBSERVATIONS                      1497418
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                 266716



 ******************************************************************************
  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       1    3600)
 FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=TRUE 
 ******************************************************************************

   REFLECTIONS OF TYPE H,0,0  0,K,0  0,0,L OR EXPECTED TO BE ABSENT (*)
   --------------------------------------------------------------------

   H    K    L  RESOLUTION  INTENSITY     SIGMA    INTENSITY/SIGMA  #OBSERVED

    0    0    8    18.787   0.1124E+02  0.2867E+00        39.20           2 
    0    0    9    16.700  -0.4473E-02  0.4634E-02        -0.97           2 
    0    0   10    15.030  -0.4229E-03  0.4041E-02        -0.10           4 
    0    0   11    13.663  -0.5964E-02  0.3776E-02        -1.58           4 
    0    0   12    12.525   0.7467E+01  0.1493E+00        50.01           4 
    0    0   13    11.561   0.3132E-02  0.4870E-02         0.64           4 
    0    0   14    10.735   0.3649E-02  0.4967E-02         0.73           4 
    0    0   15    10.020  -0.4434E-03  0.5482E-02        -0.08           4 
    0    0   16     9.393   0.5384E-01  0.1047E-01         5.14           4 
    0    0   17     8.841   0.6187E-02  0.5704E-02         1.08           4 
    0    0   18     8.350  -0.5711E-02  0.6994E-02        -0.82           3 
    0    0   19     7.910  -0.1464E-02  0.9449E-02        -0.15           2 
    0    0   20     7.515   0.1169E-01  0.7844E-02         1.49           4 
    0    0   21     7.157   0.3395E-02  0.6573E-02         0.52           4 
    0    0   22     6.832   0.1749E-02  0.7364E-02         0.24           4 
    0    0   23     6.535   0.1664E-02  0.8444E-02         0.20           4 
    0    0   24     6.262   0.2857E+00  0.2586E-01        11.04           4 
    0    0   25     6.012  -0.2006E-03  0.8080E-02        -0.02           4 
    0    0   26     5.781  -0.1250E-01  0.1108E-01        -1.13           2 
    0    0   27     5.567   0.1988E-01  0.1472E-01         1.35           2 
    0    0   28     5.368   0.4181E+02  0.8649E+00        48.34           3 
    0    0   29     5.183  -0.6418E-02  0.1017E-01        -0.63           4 
    0    0   30     5.010   0.1498E-01  0.1085E-01         1.38           4 
    0    0   31     4.848  -0.1781E-02  0.1175E-01        -0.15           4 
    0    0   32     4.697   0.3036E+02  0.5712E+00        53.16           4 
    0    0   33     4.554  -0.8372E-02  0.9884E-02        -0.85           4 
    0    0   34     4.420   0.3389E-02  0.1101E-01         0.31           4 
    0    0   35     4.294   0.3354E-01  0.1610E-01         2.08           4 
    0    0   36     4.175   0.4329E+00  0.4006E-01        10.81           4 
    0    0   37     4.062  -0.2134E-01  0.1179E-01        -1.81           4 
    0    0   38     3.955  -0.1857E-01  0.1116E-01        -1.66           4 
    0    0   39     3.854  -0.4478E-02  0.1263E-01        -0.35           4 
    0    0   40     3.757   0.1810E+02  0.3840E+00        47.13           4 
    0    0   41     3.666  -0.8949E-02  0.1903E-01        -0.47           3 
    0    0   42     3.578  -0.1045E-01  0.2551E-01        -0.41           2 
    0    0   43     3.495  -0.1341E-01  0.2389E-01        -0.56           2 
    0    0   44     3.416   0.7860E+01  0.2131E+00        36.88           4 
    0    0   45     3.340  -0.1340E-01  0.1452E-01        -0.92           4 
    0    0   46     3.267   0.4026E-02  0.1632E-01         0.25           4 
    0    0   47     3.198   0.9117E-02  0.1520E-01         0.60           4 
    0    0   48     3.131   0.7860E+00  0.6019E-01        13.06           4 
    0    0   49     3.067   0.1839E-01  0.1825E-01         1.01           4 
    0    0   50     3.006  -0.1203E-01  0.1515E-01        -0.79           4 
    0    0   51     2.947  -0.1837E-01  0.1391E-01        -1.32           4 
    0    0   52     2.890   0.1112E+02  0.2837E+00        39.19           4 
    0    0   53     2.836  -0.1873E-01  0.1310E-01        -1.43           4 
    0    0   54     2.783  -0.1341E-01  0.1388E-01        -0.97           4 
    0    0   55     2.733  -0.3781E-03  0.1817E-01        -0.02           4 
    0    0   56     2.684   0.1884E+02  0.4261E+00        44.22           4 
    0    0   57     2.637   0.5761E-02  0.1690E-01         0.34           4 
    0    0   58     2.591   0.1854E-02  0.1690E-01         0.11           2 
    0    0   59     2.547  -0.2439E-01  0.1583E-01        -1.54           2 
    0    0   60     2.505   0.1197E+02  0.3143E+00        38.09           4 
    0    0   61     2.464  -0.5732E-02  0.1317E-01        -0.44           4 
    0    0   62     2.424  -0.1738E-01  0.1376E-01        -1.26           4 
    0    0   63     2.386  -0.6453E-02  0.1683E-01        -0.38           4 
    0    0   64     2.348   0.4265E+01  0.1688E+00        25.27           4 
    0    0   65     2.312   0.6649E-02  0.1686E-01         0.39           4 
    0    0   66     2.277  -0.1633E-01  0.1577E-01        -1.04           4 
    0    0   67     2.243  -0.8157E-03  0.1811E-01        -0.05           4 
    0    0   68     2.210   0.8015E+00  0.7410E-01        10.82           4 
    0    0   69     2.178  -0.4181E-02  0.1921E-01        -0.22           4 
    0    0   70     2.147  -0.1512E-01  0.1807E-01        -0.84           4 
    0    0   71     2.117  -0.8621E-02  0.2052E-01        -0.42           4 
    0    0   72     2.087   0.4215E+00  0.9233E-01         4.57           2 
    0    0   73     2.059  -0.3764E-01  0.2808E-01        -1.34           2 
    0    0   74     2.031  -0.7583E-02  0.1756E-01        -0.43           4 
    0    0   75     2.004  -0.7843E-02  0.1757E-01        -0.45           4 
    0    0   76     1.978   0.8903E+00  0.8233E-01        10.81           4 
    0    0   77     1.952  -0.2431E-01  0.1835E-01        -1.32           4 
    0    0   78     1.927   0.1784E-02  0.2353E-01         0.08           4 
    0    0   79     1.902   0.1289E-01  0.1961E-01         0.66           4 
    0    0   80     1.879   0.1349E+00  0.3720E-01         3.63           4 
    0    0   81     1.856  -0.2046E-01  0.1517E-01        -1.35           4 
    0    0   82     1.833  -0.3550E-02  0.2007E-01        -0.18           4 
    0    0   83     1.811   0.1664E-02  0.1990E-01         0.08           4 
    0    0   84     1.789   0.4671E+00  0.6395E-01         7.30           4 
    0    0   85     1.768  -0.1952E-01  0.2023E-01        -0.97           3 
    0    0   86     1.748  -0.1517E-01  0.1749E-01        -0.87           2 
    0    0   87     1.728  -0.1116E-01  0.1494E-01        -0.75           4 
    0    0   88     1.708  -0.1019E-01  0.1457E-01        -0.70           4 
    0    0   89     1.689  -0.2643E-02  0.1711E-01        -0.15           4 
    0    0   90     1.670  -0.1021E-01  0.1389E-01        -0.73           4 
    0    0   91     1.652  -0.2161E-02  0.1934E-01        -0.11           4 
    0    0   92     1.634   0.1727E-01  0.2035E-01         0.85           4 
    0    0   93     1.616   0.3240E-02  0.2091E-01         0.15           4 
    0    0   94     1.599   0.1603E-03  0.1654E-01         0.01           4 
    0    0   95     1.582  -0.6294E-03  0.2162E-01        -0.03           3 
    0    0   96     1.566   0.4836E+00  0.7269E-01         6.65           4 
    0    0   97     1.549  -0.1365E-01  0.2992E-01        -0.46           2 
    0    0   98     1.534   0.1095E-01  0.3124E-01         0.35           3 
    0    0   99     1.518  -0.1367E-01  0.1876E-01        -0.73           4 
    0    0  100     1.503   0.1599E+01  0.1248E+00        12.81           4 
    0    0  101     1.488  -0.1424E-01  0.2075E-01        -0.69           4 
    0    0  102     1.473  -0.1562E-01  0.1947E-01        -0.80           4 
    0    0  103     1.459  -0.1135E-01  0.1932E-01        -0.59           4 
    0    0  104     1.445   0.2333E+01  0.1681E+00        13.88           3 
    0    0  105     1.431  -0.2514E-01  0.1842E-01        -1.36           4 
    0    0  106     1.418  -0.1124E-01  0.2039E-01        -0.55           4 
    0    0  107     1.405  -0.1330E-01  0.2062E-01        -0.64           2 
    0    0  108     1.392   0.1595E+00  0.6190E-01         2.58           2 
    0    0  109     1.379  -0.5482E-02  0.2644E-01        -0.21           2 
    0    0  110     1.366  -0.1859E-01  0.3454E-01        -0.54           2 
    0    0  111     1.354  -0.2717E-01  0.2519E-01        -1.08           2 
    0    0  112     1.342   0.1697E+00  0.6232E-01         2.72           2 
    0    0  113     1.330   0.5834E-02  0.3813E-01         0.15           2 
    0    0  114     1.318   0.6311E-01  0.6192E-01         1.02           1 
    0    0  115     1.307   0.3187E-01  0.3602E-01         0.88           2 
    0    0  116     1.296  -0.1482E-01  0.5700E-01        -0.26           1 
    0    3    0    19.318  -0.3813E-02  0.4634E-02        -0.82           2 
    0    4    0    14.488   0.1089E+00  0.9418E-02        11.57           4 
    0    5    0    11.591  -0.2089E-02  0.4239E-02        -0.49           4 
    0    6    0     9.659   0.1144E+01  0.3869E-01        29.56           4 
    0    7    0     8.279   0.1283E-01  0.6824E-02         1.88           4 
    0    8    0     7.244   0.1405E+00  0.1562E-01         9.00           4 
    0    9    0     6.439   0.2099E-02  0.8052E-02         0.26           4 
    0   10    0     5.795   0.1809E+01  0.8922E-01        20.28           2 
    0   11    0     5.268   0.6878E-02  0.1088E-01         0.63           4 
    0   12    0     4.829   0.2482E+01  0.1149E+00        21.60           2 
    0   13    0     4.458  -0.9791E-02  0.1091E-01        -0.90           4 
    0   14    0     4.140   0.1036E+01  0.5457E-01        18.98           4 
    0   15    0     3.864   0.3516E-01  0.2616E-01         1.34           2 
    0   16    0     3.622   0.4312E+01  0.1290E+00        33.43           4 
    0   17    0     3.409  -0.1870E-01  0.1635E-01        -1.14           4 
    0   18    0     3.220   0.7038E+01  0.1853E+00        37.97           4 
    0   19    0     3.050  -0.1539E-01  0.1570E-01        -0.98           4 
    0   20    0     2.898   0.9054E+01  0.2307E+00        39.25           4 
    0   21    0     2.760  -0.3943E-02  0.2542E-01        -0.16           2 
    0   22    0     2.634   0.3370E+01  0.1295E+00        26.02           4 
    0   23    0     2.520   0.3677E-03  0.1978E-01         0.02           4 
    0   24    0     2.415   0.2462E+00  0.4070E-01         6.05           4 
    0   25    0     2.318   0.2616E-01  0.2114E-01         1.24           4 
    0   26    0     2.229   0.1077E+01  0.1039E+00        10.37           2 
    0   27    0     2.146  -0.1319E-01  0.2241E-01        -0.59           4 
    0   28    0     2.070   0.2036E+01  0.1108E+00        18.38           4 
    0   29    0     1.998  -0.6184E-02  0.2230E-01        -0.28           4 
    0   30    0     1.932   0.2554E+01  0.1315E+00        19.43           4 
    0   31    0     1.869  -0.8931E-02  0.3535E-01        -0.25           2 
    0   32    0     1.811   0.1021E+01  0.8729E-01        11.69           4 
    0   33    0     1.756  -0.2626E-02  0.2278E-01        -0.12           4 
    0   34    0     1.705   0.1193E+00  0.3602E-01         3.31           4 
    0   35    0     1.656  -0.2106E-01  0.1932E-01        -1.09           4 
    0   36    0     1.610   0.5268E-01  0.3018E-01         1.75           4 
    0   37    0     1.566  -0.1821E-01  0.2103E-01        -0.87           4 
    0   38    0     1.525   0.2223E-01  0.2692E-01         0.83           4 
    0   39    0     1.486  -0.4922E-02  0.1493E-01        -0.33           2 
    0   40    0     1.449   0.6238E-01  0.5340E-01         1.17           2 
    0   41    0     1.413  -0.3057E-01  0.3526E-01        -0.87           2 
    0   42    0     1.380   0.4490E-01  0.5600E-01         0.80           2 
    3    0    0    19.316  -0.6609E-03  0.2274E-02        -0.29           4 
    4    0    0    14.487   0.9640E-01  0.8371E-02        11.52           4 
    5    0    0    11.589   0.1450E-01  0.6825E-02         2.12           2 
    6    0    0     9.658   0.1204E+01  0.5366E-01        22.43           2 
    7    0    0     8.278  -0.5292E-02  0.4045E-02        -1.31           4 
    8    0    0     7.243   0.1852E+00  0.1626E-01        11.39           4 
    9    0    0     6.439  -0.5583E-02  0.5568E-02        -1.00           4 
   10    0    0     5.795   0.2011E+01  0.6324E-01        31.80           4 
   11    0    0     5.268  -0.9174E-02  0.6527E-02        -1.41           4 
   12    0    0     4.829   0.2627E+01  0.8029E-01        32.72           4 
   13    0    0     4.457   0.6379E-02  0.9041E-02         0.71           4 
   14    0    0     4.139   0.1147E+01  0.5234E-01        21.91           4 
   15    0    0     3.863  -0.1035E-01  0.1080E-01        -0.96           4 
   16    0    0     3.622   0.4944E+01  0.1326E+00        37.27           4 
   17    0    0     3.409  -0.2111E-01  0.1095E-01        -1.93           4 
   18    0    0     3.219   0.7977E+01  0.2736E+00        29.15           2 
   19    0    0     3.050  -0.1604E-01  0.1741E-01        -0.92           2 
   20    0    0     2.897   0.9641E+01  0.2292E+00        42.06           4 
   21    0    0     2.759   0.2637E-02  0.1339E-01         0.20           4 
   22    0    0     2.634   0.3557E+01  0.1223E+00        29.09           4 
   23    0    0     2.519  -0.7124E-02  0.1349E-01        -0.53           4 
   24    0    0     2.414   0.2520E+00  0.3476E-01         7.25           4 
   25    0    0     2.318  -0.9483E-02  0.1436E-01        -0.66           4 
   26    0    0     2.229   0.1468E+01  0.9165E-01        16.02           3 
   27    0    0     2.146  -0.4387E-02  0.1541E-01        -0.28           4 
   28    0    0     2.070   0.2175E+01  0.1034E+00        21.03           4 
   29    0    0     1.998   0.1704E-01  0.2354E-01         0.72           3 
   30    0    0     1.932   0.2463E+01  0.1640E+00        15.02           2 
   31    0    0     1.869   0.1909E-01  0.2078E-01         0.92           3 
   32    0    0     1.811   0.1415E+01  0.8779E-01        16.12           4 
   33    0    0     1.756  -0.1101E-02  0.1556E-01        -0.07           4 
   34    0    0     1.704   0.7819E-01  0.2537E-01         3.08           4 
   35    0    0     1.656  -0.1042E-01  0.1246E-01        -0.84           4 
   36    0    0     1.610   0.3665E-01  0.2241E-01         1.64           4 
   37    0    0     1.566   0.5681E-03  0.1482E-01         0.04           4 
   38    0    0     1.525   0.6134E-01  0.2416E-01         2.54           4 
   39    0    0     1.486   0.5161E-02  0.2618E-01         0.20           2 
   40    0    0     1.449   0.6310E-02  0.2852E-01         0.22           2 
   41    0    0     1.413   0.2123E-01  0.3040E-01         0.70           2 
   42    0    0     1.380   0.8585E-01  0.4233E-01         2.03           2 

          COMPLETENESS AND QUALITY OF DATA SET
          ------------------------------------

 R-FACTOR
 observed = (SUM(ABS(I(h,i)-I(h))))/(SUM(I(h,i)))
 expected = expected R-FACTOR derived from Sigma(I)

 COMPARED = number of reflections used for calculating R-FACTOR
 I/SIGMA  = mean of intensity/Sigma(I) of unique reflections
            (after merging symmetry-related observations)
 Sigma(I) = standard deviation of reflection intensity I
            estimated from sample statistics

 R-meas   = redundancy independent R-factor (intensities)
            Diederichs & Karplus (1997), Nature Struct. Biol. 4, 269-275.

 CC(1/2)  = percentage of correlation between intensities from
            random half-datasets. Correlation significant at
            the 0.1% level is marked by an asterisk.
            Karplus & Diederichs (2012), Science 336, 1030-33
 Anomal   = percentage of correlation between random half-sets
  Corr      of anomalous intensity differences. Correlation
            significant at the 0.1% level is marked.
 SigAno   = mean anomalous difference in units of its estimated
            standard deviation (|F(+)-F(-)|/Sigma). F(+), F(-)
            are structure factor estimates obtained from the
            merged intensity observations in each parity class.
  Nano    = Number of unique reflections used to calculate
            Anomal_Corr & SigAno. At least two observations
            for each (+ and -) parity are required.


 SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION
 RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano
   LIMIT     OBSERVED  UNIQUE  POSSIBLE     OF DATA   observed  expected                                      Corr

     3.42       93737   13737     13827       99.3%       4.4%      4.2%    93732   35.80      4.8%    99.9*    -2    0.877   12441
     2.43      168181   24038     24044      100.0%       8.0%      7.0%   168174   20.93      8.6%    99.7*    -2    0.915   22552
     1.99      212993   30711     30720      100.0%      14.4%     12.8%   212976   11.71     15.5%    98.9*    -2    0.897   28978
     1.73      250880   36243     36307       99.8%      25.8%     24.5%   250862    5.99     27.9%    97.0*    -1    0.856   34441
     1.55      280908   40841     41066       99.5%      47.1%     46.8%   280879    2.87     50.9%    91.0*    -1    0.808   38726
     1.41      255658   44649     45168       98.9%      75.6%     80.1%   255321    1.27     83.1%    74.6*    -1    0.736   39602
     1.31      157417   40752     49241       82.8%     113.3%    128.5%   155655    0.44    131.1%    35.7*     0    0.676   22279
     1.22       64370   25904     52713       49.1%     153.6%    187.2%    61524    0.13    195.2%    13.8*     0    0.629    4608
     1.15       13274    9841     56057       17.6%     165.8%    211.6%     6506    0.06    228.7%    11.5*    24    0.667      44
    total     1497418  266716    349143       76.4%      12.6%     12.1%  1485629    6.62     13.7%    99.8*    -1    0.813  203671


 NUMBER OF REFLECTIONS IN SELECTED SUBSET OF IMAGES   1503511
 NUMBER OF REJECTED MISFITS                              6093
 NUMBER OF SYSTEMATIC ABSENT REFLECTIONS                    0
 NUMBER OF ACCEPTED OBSERVATIONS                      1497418
 NUMBER OF UNIQUE ACCEPTED REFLECTIONS                 266716



 ******************************************************************************
    WILSON STATISTICS OF DATA SET  "XDS_ASCII.HKL"
 ******************************************************************************

 Data is divided into resolution shells and a straight line 
 A - 2*B*SS is fitted to log<I>, where
   RES    = mean resolution (Angstrom) in shell
   SS     = mean of (sin(THETA)/LAMBDA)**2 in shell
   <I>    = mean reflection intensity in shell
   BO     = (A - log<I>)/(2*SS)
    #     = number of reflections in resolution shell

   WILSON LINE (using all data) : A=   1.307 B=  17.090 CORRELATION=  0.99
      #      RES      SS        <I>       log(<I>)       BO
    4739     6.290   0.006  2.3570E+00       0.857      35.6
    8023     3.994   0.016  3.9962E+00       1.385      -2.5
   10129     3.136   0.025  2.3267E+00       0.844       9.1
   11961     2.665   0.035  1.1087E+00       0.103      17.1
   13497     2.356   0.045  8.2202E-01      -0.196      16.7
   14920     2.135   0.055  6.2759E-01      -0.466      16.2
   16158     1.965   0.065  4.0701E-01      -0.899      17.0
   17276     1.831   0.075  2.3124E-01      -1.464      18.6
   18374     1.721   0.084  1.5142E-01      -1.888      18.9
   19350     1.629   0.094  1.0840E-01      -2.222      18.7
   20217     1.550   0.104  7.8031E-02      -2.551      18.5
   21110     1.481   0.114  5.7241E-02      -2.860      18.3
   21846     1.421   0.124  4.2749E-02      -3.152      18.0
   20044     1.369   0.133  3.2418E-02      -3.429      17.7
   16704     1.321   0.143  2.4015E-02      -3.729      17.6
   13538     1.278   0.153  2.0244E-02      -3.900      17.0
    9931     1.239   0.163  1.6939E-02      -4.078      16.5
    6266     1.204   0.173  1.4289E-02      -4.248      16.1
    2632     1.173   0.182  1.2063E-02      -4.418      15.8


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF  CENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         3<I>**2     15<I>**3    105<I>**4 

     531     6.290        1.329        2.068        3.467
     536     3.994        1.280        1.773        2.680
     540     3.136        0.820        0.663        0.533
     541     2.665        1.035        1.328        2.072
     530     2.356        1.772        2.575        3.323
     549     2.135        1.550        1.864        2.141
     530     1.965        1.689        2.191        2.548
     498     1.831        1.923        2.686        3.389
     501     1.721        1.136        1.064        0.923
     481     1.629        2.701        6.922       18.804
     477     1.550        1.802        1.965        1.791
     442     1.481        3.683       10.094       26.607
     447     1.421        4.742       29.889      205.504
     301     1.369        2.164        5.060       12.985
     223     1.321        3.271        6.707       13.829
     150     1.278        4.748       13.263       35.029
     103     1.239        3.405        4.517        5.548
      56     1.204        5.358        9.873       17.402
      27     1.173       68.612      690.243     5475.288
    7463   overall        2.293        7.409       38.358


 HIGHER ORDER MOMENTS OF WILSON DISTRIBUTION OF ACENTRIC DATA
    AS COMPARED WITH THEORETICAL VALUES. (EXPECTED: 1.00)    
      #      RES        <I**2>/      <I**3>/      <I**4>/  
                         2<I>**2      6<I>**3     24<I>**4 

    4208     6.290        1.213        1.644        2.332
    7487     3.994        1.076        1.317        1.935
    9589     3.136        1.094        1.238        1.421
   11420     2.665        1.097        1.276        1.565
   12967     2.356        1.068        1.216        1.446
   14371     2.135        1.109        1.375        1.882
   15628     1.965        1.115        1.512        2.595
   16778     1.831        1.100        1.342        1.784
   17873     1.721        1.135        1.651        3.703
   18869     1.629        1.174        1.769        3.498
   19740     1.550        1.133        1.476        2.188
   20668     1.481        1.202        1.721        3.182
   21399     1.421        1.271        1.781        3.006
   19743     1.369        1.486        2.073        3.637
   16481     1.321        1.883        2.769        5.615
   13388     1.278        2.456        3.756        7.330
    9828     1.239        3.554        6.146       14.175
    6210     1.204        5.680       11.799       36.525
    2605     1.173        8.883       20.471       72.555
  259252   overall        1.563        2.370        4.996

   ======= CUMULATIVE INTENSITY DISTRIBUTION =======
 DEFINITIONS:
   <I>    = mean reflection intensity
 Na(Z)exp = expected number of acentric reflections with I <= Z*<I>
 Na(Z)obs = observed number of acentric reflections with I <= Z*<I>
 Nc(Z)exp = expected number of  centric reflections with I <= Z*<I>
 Nc(Z)obs = observed number of  centric reflections with I <= Z*<I>



 Nc(Z)obs/Nc(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

     531     6.290  1.26  1.14  1.07  1.04  1.03  1.03  1.03  1.01  0.99  0.99
     536     3.994  1.12  1.08  1.03  1.05  1.04  1.01  1.01  1.01  1.00  0.98
     540     3.136  1.10  1.04  1.02  0.97  1.00  1.00  1.00  1.00  1.00  0.99
     541     2.665  1.04  1.03  1.00  0.97  0.98  0.99  1.00  0.98  0.97  0.99
     530     2.356  1.01  1.00  0.99  1.02  1.00  1.00  0.98  0.97  0.98  0.98
     549     2.135  0.96  0.90  0.90  0.94  0.95  0.96  0.93  0.92  0.92  0.93
     530     1.965  0.92  0.93  0.93  0.94  0.93  0.95  0.95  0.92  0.93  0.93
     498     1.831  0.87  0.87  0.88  0.91  0.91  0.90  0.90  0.89  0.90  0.92
     501     1.721  0.77  0.83  0.90  0.92  0.89  0.89  0.92  0.92  0.93  0.93
     481     1.629  0.68  0.75  0.80  0.87  0.88  0.90  0.91  0.91  0.92  0.92
     477     1.550  0.95  0.88  0.89  0.89  0.89  0.90  0.87  0.88  0.87  0.88
     442     1.481  1.10  0.97  0.88  0.90  0.94  0.98  0.98  0.96  0.95  0.97
     447     1.421  1.07  0.89  0.80  0.79  0.79  0.80  0.81  0.79  0.81  0.82
     301     1.369  1.34  1.13  1.01  0.98  0.94  0.95  0.95  0.96  0.94  0.94
     223     1.321  1.41  1.14  1.02  0.98  0.94  0.93  0.95  0.93  0.89  0.89
     150     1.278  1.42  1.08  0.99  0.90  0.85  0.83  0.79  0.78  0.76  0.79
     103     1.239  1.17  0.93  0.84  0.78  0.75  0.73  0.75  0.73  0.71  0.73
      56     1.204  1.15  0.83  0.77  0.68  0.62  0.60  0.60  0.60  0.63  0.60
      27     1.173  1.79  1.29  1.07  0.94  0.85  0.79  0.74  0.71  0.68  0.65
    7463   overall  1.03  0.97  0.94  0.94  0.94  0.94  0.94  0.93  0.93  0.93


 Na(Z)obs/Na(Z)exp versus resolution and Z (0.1-1.0)
      #      RES     0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9   1.0

    4208     6.290  1.52  1.33  1.26  1.18  1.14  1.10  1.08  1.08  1.07  1.05
    7487     3.994  1.18  1.07  1.08  1.07  1.05  1.05  1.03  1.02  1.01  1.01
    9589     3.136  1.22  1.11  1.09  1.09  1.08  1.06  1.05  1.03  1.02  1.02
   11420     2.665  1.24  1.11  1.10  1.09  1.06  1.05  1.05  1.04  1.03  1.02
   12967     2.356  1.20  1.14  1.08  1.08  1.07  1.06  1.05  1.04  1.03  1.02
   14371     2.135  1.21  1.18  1.13  1.10  1.09  1.08  1.07  1.06  1.06  1.04
   15628     1.965  1.06  1.08  1.08  1.07  1.06  1.06  1.06  1.05  1.04  1.04
   16778     1.831  1.22  1.15  1.11  1.09  1.08  1.06  1.05  1.04  1.03  1.03
   17873     1.721  1.14  1.07  1.06  1.05  1.04  1.04  1.03  1.03  1.03  1.03
   18869     1.629  1.19  1.08  1.08  1.08  1.07  1.07  1.06  1.05  1.05  1.04
   19740     1.550  1.35  1.10  1.05  1.05  1.05  1.04  1.04  1.03  1.03  1.03
   20668     1.481  1.67  1.22  1.09  1.05  1.04  1.03  1.03  1.03  1.03  1.02
   21399     1.421  2.08  1.36  1.16  1.07  1.03  1.01  1.01  1.00  1.00  1.00
   19743     1.369  2.69  1.62  1.28  1.13  1.05  1.00  0.97  0.95  0.94  0.94
   16481     1.321  3.22  1.86  1.41  1.21  1.09  1.02  0.97  0.94  0.92  0.90
   13388     1.278  3.70  2.07  1.54  1.28  1.13  1.04  0.99  0.95  0.91  0.89
    9828     1.239  4.16  2.30  1.68  1.38  1.20  1.09  1.01  0.95  0.91  0.88
    6210     1.204  4.69  2.52  1.82  1.46  1.25  1.12  1.03  0.96  0.92  0.88
    2605     1.173  5.18  2.76  1.96  1.56  1.33  1.19  1.08  1.00  0.94  0.90
  259252   overall  1.92  1.37  1.20  1.12  1.08  1.05  1.03  1.01  1.00  0.99


 List of    457 reflections *NOT* obeying Wilson distribution (Z> 10.0)

   h    k    l     RES      Z     Intensity    Sigma

    9    9  -81    1.72   26.86  0.4067E+01  0.1403E+00 "alien"
    6   10  120    1.21   26.08  0.3726E+00  0.9279E-01 "alien"
   10    6  120    1.21   25.69  0.3670E+00  0.7597E-01 "alien"
    0   18  102    1.34   25.52  0.6129E+00  0.1134E+00 "alien"
   10    4  125    1.17   24.40  0.2943E+00  0.1302E+00 "alien"
   15    6  121    1.17   24.00  0.2895E+00  0.1212E+00 "alien"
    0   10  105    1.39   23.86  0.7735E+00  0.7349E-01 "alien"
    0    3  102    1.47   23.58  0.1350E+01  0.8400E-01 "alien"
    9    9   81    1.72   23.41  0.3545E+01  0.1206E+00 "alien"
    0    8   90    1.63   23.38  0.2534E+01  0.1106E+00 "alien"
    0    4  112    1.34   23.18  0.5567E+00  0.7618E-01 "alien"
   10    2  124    1.19   23.05  0.3293E+00  0.9700E-01 "alien"
    0    1  104    1.44   21.14  0.9038E+00  0.7733E-01 "alien"
    0   10  116    1.26   19.75  0.3999E+00  0.1437E+00 "alien"
    0    9  116    1.27   19.60  0.3968E+00  0.1390E+00 "alien"
   27   35  -36    1.25   19.57  0.3315E+00  0.9202E-01 "alien"
   22   19  103    1.18   19.43  0.2343E+00  0.1223E+00 "alien"
    8   29  102    1.17   19.30  0.2328E+00  0.1029E+00 "alien"
    0    1  100    1.50   19.00  0.1088E+01  0.8109E-01 "alien"
    9    2  125    1.18   18.96  0.2287E+00  0.1056E+00 "alien"
    8   14  117    1.21   18.34  0.2621E+00  0.1126E+00 "alien"
    7   12  116    1.24   18.18  0.3079E+00  0.8304E-01 "alien"
   43    4   -8    1.34   18.05  0.4336E+00  0.1200E+00 "alien"
   44    8  -52    1.18   17.95  0.2166E+00  0.1089E+00 "alien"
   16   12  113    1.21   17.85  0.2550E+00  0.1322E+00 "alien"
   19   12  111    1.20   17.78  0.2540E+00  0.1154E+00 "alien"
   13    7  121    1.18   17.65  0.2130E+00  0.1076E+00 "alien"
   12    0  123    1.18   17.59  0.7301E+00  0.1892E+00 "alien"
    7   25  109    1.17   17.42  0.2101E+00  0.9990E-01 "alien"
    6   16  118    1.19   17.36  0.2480E+00  0.1081E+00 "alien"
   11    9  -65    2.01   17.34  0.7057E+01  0.1758E+00 "alien"
   26   24   86    1.20   17.32  0.2475E+00  0.1099E+00 "alien"
    5   12  120    1.21   17.21  0.2459E+00  0.1143E+00 "alien"
    7   24  109    1.19   17.18  0.2073E+00  0.9918E-01 "alien"
    5   10  116    1.26   16.97  0.2875E+00  0.8219E-01 "alien"
    0   19  106    1.29   16.97  0.3435E+00  0.8561E-01 "alien"
   38   26  -50    1.16   16.96  0.2046E+00  0.1141E+00 "alien"
    0    0  104    1.45   16.87  0.2348E+01  0.1685E+00 "alien"
   28   26  -63    1.28   16.86  0.3413E+00  0.1264E+00 "alien"
    9   11  -65    2.01   16.77  0.6826E+01  0.1713E+00 "alien"
    9   12  120    1.19   16.67  0.2382E+00  0.1095E+00 "alien"
   11    9   65    2.01   16.66  0.6782E+01  0.1516E+00 "alien"
   13   11  -83    1.60   16.65  0.1805E+01  0.9641E-01 "alien"
   11   13  -83    1.60   16.63  0.1803E+01  0.9651E-01 "alien"
    6   20  -96    1.36   16.24  0.5265E+00  0.1052E+00 "alien"
   15    8  122    1.16   16.11  0.1943E+00  0.1007E+00 "alien"
    9   11   65    2.01   16.04  0.6530E+01  0.1522E+00 "alien"
   17    6  -95    1.42   16.04  0.6856E+00  0.9482E-01 "alien"
    3   17  113    1.24   15.97  0.2704E+00  0.7906E-01 "alien"
   10   21  102    1.27   15.84  0.3207E+00  0.8564E-01 "alien"
   16    4   79    1.67   15.65  0.1697E+01  0.8317E-01 "alien"
   18   45    5    1.19   15.57  0.2225E+00  0.1077E+00 "alien"
   37   28  -43    1.18   15.56  0.1877E+00  0.1040E+00 "alien"
    2    2 -104    1.44   15.53  0.6639E+00  0.6153E-01 "alien"
    8   13  119    1.20   15.50  0.2215E+00  0.1052E+00 "alien"
    0    8   27    4.41   15.29  0.6112E+02  0.7588E+00 "alien"
   16   15  101    1.30   15.28  0.3092E+00  0.5351E-01 "alien"
    8    6  117    1.25   15.24  0.2581E+00  0.9032E-01 "alien"
    3   21  107    1.25   15.22  0.2578E+00  0.8235E-01 "alien"
   23   24   95    1.17   15.20  0.1834E+00  0.9420E-01 "alien"
   30    5  -62    1.50   15.04  0.8611E+00  0.7342E-01 "alien"
    4   11  118    1.23   14.92  0.2527E+00  0.8263E-01 "alien"
   32   21   70    1.24   14.92  0.2527E+00  0.1228E+00 "alien"
   30   31   54    1.21   14.90  0.2130E+00  0.1132E+00 "alien"
    4    9  120    1.23   14.87  0.2519E+00  0.1235E+00 "alien"
   30   37  -33    1.18   14.84  0.1791E+00  0.1074E+00 "alien"
    8   19  112    1.21   14.84  0.2121E+00  0.9707E-01 "alien"
   26   42  -12    1.17   14.81  0.1787E+00  0.1007E+00 "alien"
   11   13   83    1.60   14.68  0.1592E+01  0.8767E-01 "alien"
   21   18   99    1.23   14.68  0.2487E+00  0.1123E+00 "alien"
   38   24  -51    1.18   14.63  0.1765E+00  0.9684E-01 "alien"
   10   22  108    1.20   14.63  0.2090E+00  0.7069E-01 "alien"
   20   17  108    1.18   14.62  0.1764E+00  0.1065E+00 "alien"
   29   30   63    1.20   14.61  0.2087E+00  0.9867E-01 "alien"
    3    6  114    1.30   14.59  0.3504E+00  0.6221E-01 "alien"
   13   11   83    1.60   14.50  0.1572E+01  0.7994E-01 "alien"
   16    3  100    1.38   14.47  0.4691E+00  0.5761E-01 "alien"
   13    8  120    1.19   14.45  0.2064E+00  0.1135E+00 "alien"
   31   27   71    1.17   14.45  0.1743E+00  0.1054E+00 "alien"
   43   22   -9    1.20   14.35  0.2051E+00  0.1067E+00 "alien"
   31   29   65    1.18   14.31  0.1727E+00  0.1006E+00 "alien"
   26   40  -10    1.21   14.21  0.2030E+00  0.1095E+00 "alien"
   11    0  126    1.16   14.19  0.4576E+00  0.1562E+00 "alien"
    9    6  121    1.21   14.18  0.2026E+00  0.7910E-01 "alien"
   16    4  -79    1.67   14.16  0.1535E+01  0.1348E+00 "alien"
   28   25   81    1.19   14.15  0.2022E+00  0.9499E-01 "alien"
   41   25   20    1.19   14.15  0.2021E+00  0.9701E-01 "alien"
    9    9  121    1.20   14.11  0.2017E+00  0.1002E+00 "alien"
   11   34   79    1.23   14.11  0.2389E+00  0.7544E-01 "alien"
   33   33  -31    1.20   14.09  0.2014E+00  0.1002E+00 "alien"
    2    2 -106    1.41   14.08  0.6020E+00  0.8829E-01 "alien"
    7    4   94    1.56   14.07  0.1098E+01  0.8166E-01 "alien"
   14    4  122    1.18   14.07  0.1697E+00  0.9787E-01 "alien"
   10    5  116    1.26   14.02  0.2375E+00  0.5013E-01 "alien"
   21   21  -72    1.43   13.96  0.5966E+00  0.8109E-01 "alien"
    5    8  119    1.24   13.93  0.2359E+00  0.7763E-01 "alien"
    7    4  -94    1.56   13.92  0.1086E+01  0.7648E-01 "alien"
   13   12  -99    1.38   13.92  0.4512E+00  0.1395E+00 "alien"
   19   13  112    1.18   13.90  0.1677E+00  0.1167E+00 "alien"
    4   16   79    1.67   13.80  0.1496E+01  0.8084E-01 "alien"
   20   42  -17    1.23   13.78  0.2334E+00  0.1186E+00 "alien"
    5   11 -110    1.31   13.75  0.3302E+00  0.1187E+00 "alien"
    4    7  -94    1.56   13.74  0.1072E+01  0.8191E-01 "alien"
   10    7  119    1.22   13.74  0.2327E+00  0.6202E-01 "alien"
   30   26   72    1.20   13.72  0.1960E+00  0.1028E+00 "alien"
   15   11  112    1.23   13.70  0.2321E+00  0.7686E-01 "alien"
   12   31   83    1.26   13.59  0.2302E+00  0.1000E+00 "alien"
   32   27   63    1.20   13.59  0.1941E+00  0.1056E+00 "alien"
   30   29   66    1.19   13.58  0.1940E+00  0.1022E+00 "alien"
   25   24   86    1.21   13.51  0.1931E+00  0.1025E+00 "alien"
   10   12  122    1.17   13.46  0.1623E+00  0.9332E-01 "alien"
    3   16  100    1.38   13.45  0.4361E+00  0.6193E-01 "alien"
    3    4  -68    2.17   13.43  0.8431E+01  0.2480E+00 "alien"
    7    7  116    1.27   13.39  0.2711E+00  0.6379E-01 "alien"
    7    7  -76    1.87   13.39  0.3096E+01  0.1357E+00 "alien"
    4    7   94    1.56   13.37  0.1043E+01  0.7185E-01 "alien"
   21   45    7    1.17   13.37  0.1612E+00  0.9293E-01 "alien"
    6   11  116    1.25   13.32  0.2257E+00  0.7151E-01 "alien"
   21   21   72    1.43   13.32  0.5694E+00  0.5461E-01 "alien"
   13   14  -91    1.45   13.31  0.7618E+00  0.8725E-01 "alien"
   15    9  118    1.19   13.31  0.1902E+00  0.1090E+00 "alien"
   46    5  -46    1.17   13.28  0.1601E+00  0.8640E-01 "alien"
   24   24   89    1.20   13.25  0.1894E+00  0.9854E-01 "alien"
   10   16  107    1.28   13.24  0.2681E+00  0.6219E-01 "alien"
    1   13 -104    1.37   13.24  0.4292E+00  0.6640E-01 "alien"
   27   38  -20    1.23   13.18  0.2233E+00  0.1156E+00 "alien"
   20   19  102    1.21   13.17  0.1882E+00  0.1125E+00 "alien"
   34   32  -38    1.18   13.17  0.1588E+00  0.1119E+00 "alien"
   34   33  -42    1.16   13.15  0.1587E+00  0.1019E+00 "alien"
   33   32  -44    1.18   13.11  0.1582E+00  0.9378E-01 "alien"
   22   44   -3    1.18   13.10  0.1581E+00  0.1009E+00 "alien"
   37   14  -65    1.24   13.10  0.2218E+00  0.1029E+00 "alien"
    4   16  117    1.21   13.10  0.1871E+00  0.1040E+00 "alien"
   19   13  113    1.18   13.09  0.1579E+00  0.1064E+00 "alien"
   11    4  125    1.17   13.07  0.1576E+00  0.1101E+00 "alien"
   12    5  120    1.21   13.04  0.1863E+00  0.7290E-01 "alien"
   22   12  105    1.22   13.01  0.1859E+00  0.1204E+00 "alien"
   26   41   -3    1.19   12.98  0.1855E+00  0.1143E+00 "alien"
    4    3  -68    2.17   12.96  0.8134E+01  0.2170E+00 "alien"
   20   24   94    1.21   12.92  0.1846E+00  0.1023E+00 "alien"
   20    5  -85    1.50   12.90  0.7384E+00  0.8606E-01 "alien"
    0   21   58    1.89   12.89  0.2980E+01  0.1256E+00 "alien"
    7   18  118    1.17   12.84  0.1549E+00  0.9090E-01 "alien"
   31   27   70    1.18   12.82  0.1547E+00  0.9119E-01 "alien"
    3    4   68    2.17   12.81  0.8038E+01  0.1760E+00 "alien"
   22   19  102    1.18   12.80  0.1544E+00  0.9356E-01 "alien"
   10   15  119    1.18   12.78  0.1542E+00  0.8621E-01 "alien"
    9    2  122    1.21   12.78  0.1827E+00  0.5107E-01 "alien"
   28    5  -79    1.39   12.75  0.4134E+00  0.6719E-01 "alien"
   44   15   -9    1.24   12.75  0.2160E+00  0.1041E+00 "alien"
   14   14  -14    2.82   12.71  0.1409E+02  0.2310E+00 "alien"
   11   40   72    1.16   12.70  0.1532E+00  0.8737E-01 "alien"
   24   11  -76    1.47   12.69  0.7264E+00  0.8390E-01 "alien"
   41   14  -54    1.21   12.68  0.1811E+00  0.1078E+00 "alien"
   10    5  125    1.17   12.66  0.1527E+00  0.9969E-01 "alien"
    6    8  117    1.25   12.64  0.2141E+00  0.6496E-01 "alien"
   22    7  -76    1.55   12.62  0.9849E+00  0.7893E-01 "alien"
    4   16  -79    1.67   12.56  0.1362E+01  0.9754E-01 "alien"
   43    6  -59    1.18   12.55  0.1514E+00  0.8874E-01 "alien"
   27   41  -20    1.17   12.52  0.1511E+00  0.1088E+00 "alien"
   32   18   66    1.30   12.52  0.2535E+00  0.7729E-01 "alien"
   23   20   99    1.19   12.49  0.1785E+00  0.1035E+00 "alien"
   45    6  -54    1.16   12.49  0.1506E+00  0.9368E-01 "alien"
    9    9  -75    1.83   12.48  0.2886E+01  0.1215E+00 "alien"
    0    5   71    2.08   12.48  0.7832E+01  0.1809E+00 "alien"
   23   41   -5    1.23   12.47  0.2113E+00  0.1081E+00 "alien"
   15   15    0    2.73   12.47  0.1382E+02  0.2194E+00 "alien"
   11    7  120    1.21   12.45  0.1779E+00  0.5676E-01 "alien"
   32   12  -73    1.31   12.41  0.2980E+00  0.8782E-01 "alien"
   35   24  -57    1.21   12.41  0.1773E+00  0.1003E+00 "alien"
   13   42   59    1.17   12.40  0.1496E+00  0.8470E-01 "alien"
   32   19   70    1.26   12.40  0.2510E+00  0.1095E+00 "alien"
    5    5  -90    1.64   12.38  0.1341E+01  0.9985E-01 "alien"
   17    4  -75    1.72   12.36  0.1872E+01  0.9702E-01 "alien"
   18   15  105    1.24   12.35  0.2092E+00  0.7467E-01 "alien"
    8    6   31    3.72   12.35  0.4935E+02  0.6091E+00 "alien"
   36   30  -38    1.18   12.35  0.1490E+00  0.9608E-01 "alien"
   25   23   77    1.28   12.35  0.2500E+00  0.5616E-01 "alien"
    7    7   76    1.87   12.35  0.2855E+01  0.1079E+00 "alien"
    0   12   55    2.38   12.34  0.1014E+02  0.1996E+00 "alien"
   20   19  107    1.17   12.31  0.1485E+00  0.8852E-01 "alien"
    9    9   75    1.83   12.28  0.2840E+01  0.1162E+00 "alien"
    1   14 -113    1.27   12.28  0.2485E+00  0.1177E+00 "alien"
    2    2  104    1.44   12.28  0.5248E+00  0.5498E-01 "alien"
    6    3  114    1.30   12.24  0.2940E+00  0.5825E-01 "alien"
    6    8   31    3.72   12.24  0.4891E+02  0.7067E+00 "alien"
    0   11   10    4.97   12.24  0.2885E+02  0.3734E+00 "alien"
    3   14  117    1.22   12.24  0.2073E+00  0.9722E-01 "alien"
   42   25   13    1.18   12.23  0.1475E+00  0.9052E-01 "alien"
    8   15  115    1.22   12.19  0.2065E+00  0.9974E-01 "alien"
    6    8  -31    3.72   12.19  0.4872E+02  0.6315E+00 "alien"
   43   13  -53    1.17   12.17  0.1468E+00  0.9359E-01 "alien"
    3   11  117    1.25   12.17  0.2061E+00  0.6843E-01 "alien"
   15    4  118    1.21   12.14  0.1735E+00  0.9723E-01 "alien"
   35   27  -44    1.22   12.14  0.2056E+00  0.1017E+00 "alien"
   16   16   79    1.53   12.12  0.9460E+00  0.6112E-01 "alien"
   14   45   14    1.22   12.09  0.2048E+00  0.1091E+00 "alien"
   17    4   75    1.72   12.08  0.1830E+01  0.7985E-01 "alien"
   16   16  -79    1.53   12.08  0.9425E+00  0.7232E-01 "alien"
   24   44   -6    1.16   12.08  0.1457E+00  0.1140E+00 "alien"
   11    4  118    1.23   12.07  0.2045E+00  0.5396E-01 "alien"
   31   28   52    1.25   12.07  0.2045E+00  0.6915E-01 "alien"
   14   14   14    2.82   12.06  0.1337E+02  0.2205E+00 "alien"
    1    9 -115    1.28   12.04  0.2437E+00  0.1164E+00 "alien"
    8    6  -31    3.72   12.03  0.4808E+02  0.6278E+00 "alien"
    8   26  103    1.20   12.03  0.1718E+00  0.8829E-01 "alien"
   16   16   94    1.36   12.02  0.3898E+00  0.5228E-01 "alien"
   27   34  -46    1.24   12.02  0.2036E+00  0.1064E+00 "alien"
   20    5   54    1.98   12.01  0.4887E+01  0.1122E+00 "alien"
   25   42   -5    1.18   12.01  0.1448E+00  0.9001E-01 "alien"
   21   17  105    1.19   12.00  0.1715E+00  0.9730E-01 "alien"
   27   26   86    1.16   11.99  0.1447E+00  0.9750E-01 "alien"
    0   23   81    1.49   11.98  0.6856E+00  0.7864E-01 "alien"
    2   15  111    1.28   11.97  0.2423E+00  0.7419E-01 "alien"
    2    2  106    1.41   11.96  0.5114E+00  0.5884E-01 "alien"
    0   13  112    1.28   11.95  0.2420E+00  0.7849E-01 "alien"
   29   10  -32    1.75   11.95  0.1809E+01  0.1170E+00 "alien"
    0    0  100    1.50   11.94  0.1564E+01  0.1268E+00 "alien"
   11    6  123    1.18   11.94  0.1441E+00  0.8757E-01 "alien"
    9    7  117    1.25   11.94  0.2022E+00  0.5000E-01 "alien"
    8   13  121    1.18   11.94  0.1440E+00  0.9133E-01 "alien"
    7    1  120    1.24   11.88  0.2013E+00  0.5088E-01 "alien"
    4    9  -90    1.61   11.88  0.1287E+01  0.8386E-01 "alien"
    8   15  122    1.16   11.87  0.1432E+00  0.8682E-01 "alien"
   37   27   27    1.23   11.86  0.2010E+00  0.1029E+00 "alien"
   25   43  -11    1.16   11.86  0.1431E+00  0.9353E-01 "alien"
    6    5 -109    1.36   11.86  0.3844E+00  0.7282E-01 "alien"
    0    5   75    1.97   11.84  0.4818E+01  0.1388E+00 "alien"
    8    0  -90    1.63   11.83  0.2910E+01  0.1728E+00 "alien"
   22   19   91    1.27   11.80  0.2388E+00  0.6622E-01 "alien"
   12    3  121    1.20   11.79  0.1685E+00  0.7167E-01 "alien"
    4    3   68    2.17   11.78  0.7393E+01  0.1656E+00 "alien"
   23   19  101    1.18   11.73  0.1415E+00  0.8388E-01 "alien"
    8   17  118    1.18   11.71  0.1413E+00  0.9737E-01 "alien"
   13   33   78    1.25   11.68  0.1978E+00  0.9793E-01 "alien"
   16    3  -95    1.45   11.66  0.4983E+00  0.7274E-01 "alien"
   41   27   15    1.17   11.64  0.1404E+00  0.9841E-01 "alien"
    5   20  114    1.19   11.64  0.1663E+00  0.9157E-01 "alien"
   20   19  105    1.18   11.63  0.1403E+00  0.9871E-01 "alien"
   45   17  -28    1.18   11.61  0.1401E+00  0.9528E-01 "alien"
    3   16  -95    1.45   11.61  0.4963E+00  0.6010E-01 "alien"
   41   20  -52    1.16   11.58  0.1397E+00  0.8850E-01 "alien"
    4   13 -102    1.39   11.57  0.3751E+00  0.6462E-01 "alien"
    4    9   90    1.61   11.57  0.1254E+01  0.8864E-01 "alien"
    5   20  -54    1.98   11.56  0.4704E+01  0.1343E+00 "alien"
   44   17  -17    1.22   11.56  0.1651E+00  0.8772E-01 "alien"
   22   18  104    1.18   11.55  0.1394E+00  0.8798E-01 "alien"
    7   13  122    1.18   11.54  0.1393E+00  0.8306E-01 "alien"
    5    4  -98    1.51   11.52  0.6593E+00  0.6144E-01 "alien"
   26    5   87    1.36   11.51  0.3733E+00  0.7715E-01 "alien"
   17   45    9    1.20   11.50  0.1644E+00  0.9962E-01 "alien"
   19   15  111    1.18   11.48  0.1385E+00  0.9630E-01 "alien"
   14   13  111    1.24   11.47  0.1942E+00  0.6495E-01 "alien"
   15   18  108    1.21   11.46  0.1637E+00  0.9652E-01 "alien"
    5    5  115    1.29   11.46  0.2319E+00  0.5243E-01 "alien"
   17   43   28    1.22   11.45  0.1940E+00  0.1071E+00 "alien"
   45   12  -38    1.19   11.45  0.1636E+00  0.9721E-01 "alien"
   20    5  -54    1.98   11.44  0.4657E+01  0.1414E+00 "alien"
    0    7   60    2.40   11.43  0.9397E+01  0.1916E+00 "alien"
   22   12  -87    1.38   11.43  0.3705E+00  0.8891E-01 "alien"
   12   42   49    1.22   11.40  0.1630E+00  0.9504E-01 "alien"
    1    1   83    1.81   11.40  0.2636E+01  0.1045E+00 "alien"
   17   22  100    1.22   11.37  0.1624E+00  0.9540E-01 "alien"
   32   29   51    1.22   11.36  0.1925E+00  0.7593E-01 "alien"
   39   22  -54    1.17   11.36  0.1370E+00  0.8785E-01 "alien"
    3   12  116    1.25   11.35  0.1923E+00  0.7365E-01 "alien"
   42   12  -52    1.21   11.34  0.1620E+00  0.1100E+00 "alien"
    4    7  114    1.30   11.33  0.2294E+00  0.5300E-01 "alien"
    7   22  -76    1.55   11.33  0.8843E+00  0.6825E-01 "alien"
   14    2  -82    1.67   11.32  0.1228E+01  0.1229E+00 "alien"
   13    9  119    1.19   11.31  0.1616E+00  0.1019E+00 "alien"
    5   42   49    1.25   11.31  0.1915E+00  0.9999E-01 "alien"
   41   10  -61    1.20   11.30  0.1615E+00  0.1107E+00 "alien"
   10   10   90    1.55   11.30  0.8816E+00  0.7324E-01 "alien"
    0   17   43    2.44   11.30  0.9287E+01  0.2159E+00 "alien"
   15   14  -86    1.49   11.29  0.6461E+00  0.9187E-01 "alien"
   11    5  121    1.20   11.27  0.1611E+00  0.5786E-01 "alien"
   26   25   83    1.20   11.27  0.1610E+00  0.9066E-01 "alien"
   31   33  -42    1.21   11.27  0.1610E+00  0.7536E-01 "alien"
   22    7   76    1.55   11.27  0.8792E+00  0.5188E-01 "alien"
    9    1  124    1.19   11.26  0.1610E+00  0.6602E-01 "alien"
   12    4  123    1.18   11.24  0.1356E+00  0.5594E-01 "alien"
   19   40   37    1.25   11.24  0.1904E+00  0.1033E+00 "alien"
   14   16  108    1.24   11.23  0.1903E+00  0.7997E-01 "alien"
   10    2  122    1.20   11.23  0.1604E+00  0.5811E-01 "alien"
   10    4  -68    2.04   11.22  0.4568E+01  0.1369E+00 "alien"
    0   10   98    1.48   11.22  0.6424E+00  0.5883E-01 "alien"
    5   28  103    1.19   11.22  0.1603E+00  0.9314E-01 "alien"
    3   16 -100    1.38   11.22  0.3638E+00  0.6348E-01 "alien"
   13   44   46    1.18   11.22  0.1353E+00  0.9122E-01 "alien"
    1    1  -83    1.81   11.21  0.2593E+01  0.1053E+00 "alien"
   33   25   62    1.21   11.21  0.1602E+00  0.6439E-01 "alien"
    7    6  124    1.19   11.21  0.1601E+00  0.8611E-01 "alien"
    4    4  -96    1.55   11.20  0.8742E+00  0.6964E-01 "alien"
    4   10  -68    2.04   11.20  0.4557E+01  0.1331E+00 "alien"
   23   40   12    1.25   11.18  0.1893E+00  0.9950E-01 "alien"
   23   44    4    1.17   11.18  0.1348E+00  0.9431E-01 "alien"
   10    0  116    1.26   11.17  0.4916E+00  0.1014E+00 "alien"
   16    7  120    1.17   11.16  0.1347E+00  0.9880E-01 "alien"
   10    4   68    2.04   11.16  0.4542E+01  0.1185E+00 "alien"
    4    4  115    1.30   11.16  0.2259E+00  0.5109E-01 "alien"
   18   14  111    1.19   11.15  0.1594E+00  0.9460E-01 "alien"
   14   13  -91    1.45   11.14  0.6379E+00  0.7905E-01 "alien"
    5   20   54    1.98   11.12  0.4528E+01  0.1278E+00 "alien"
   38   25  -20    1.26   11.10  0.1881E+00  0.9100E-01 "alien"
   12   40   67    1.18   11.10  0.1339E+00  0.9201E-01 "alien"
   44   19  -13    1.20   11.09  0.1585E+00  0.8381E-01 "alien"
   36   16  -65    1.24   11.08  0.1878E+00  0.9565E-01 "alien"
    4   20  105    1.28   11.08  0.2244E+00  0.8177E-01 "alien"
   41   19  -49    1.18   11.08  0.1337E+00  0.8926E-01 "alien"
   43   23  -13    1.18   11.07  0.1336E+00  0.8276E-01 "alien"
   17   19  107    1.19   11.06  0.1580E+00  0.9188E-01 "alien"
    5   23  -97    1.31   11.03  0.2649E+00  0.1167E+00 "alien"
    5   11  110    1.31   11.02  0.2647E+00  0.4762E-01 "alien"
   10   10  -99    1.42   11.02  0.4711E+00  0.6876E-01 "alien"
   45   10  -42    1.19   11.02  0.1574E+00  0.8688E-01 "alien"
   17   43   29    1.22   11.01  0.1574E+00  0.9220E-01 "alien"
   26   24   87    1.19   10.99  0.1571E+00  0.9916E-01 "alien"
   41   26    6    1.19   10.99  0.1570E+00  0.8684E-01 "alien"
    9   13  122    1.17   10.98  0.1325E+00  0.8152E-01 "alien"
    8   31   98    1.17   10.97  0.1323E+00  0.8817E-01 "alien"
   10    0 -105    1.39   10.95  0.7603E+00  0.1232E+00 "alien"
   10   10  -90    1.55   10.94  0.8534E+00  0.6728E-01 "alien"
   39   27   23    1.20   10.93  0.1562E+00  0.9930E-01 "alien"
   38   23  -56    1.17   10.90  0.1315E+00  0.9915E-01 "alien"
    4    5  -98    1.51   10.90  0.6237E+00  0.6161E-01 "alien"
   29   40  -22    1.16   10.89  0.1313E+00  0.9534E-01 "alien"
    3    2  120    1.25   10.89  0.1844E+00  0.9024E-01 "alien"
   12   12  -93    1.46   10.83  0.6198E+00  0.9405E-01 "alien"
   19    3  -86    1.51   10.82  0.6196E+00  0.7858E-01 "alien"
    4   14  118    1.21   10.82  0.1546E+00  0.9408E-01 "alien"
   28   28   78    1.17   10.82  0.1305E+00  0.9231E-01 "alien"
   13    5  125    1.16   10.81  0.1304E+00  0.8341E-01 "alien"
    9    4   86    1.68   10.80  0.1636E+01  0.8116E-01 "alien"
   38   27   24    1.22   10.79  0.1542E+00  0.8734E-01 "alien"
    0   15   97    1.44   10.78  0.4609E+00  0.5787E-01 "alien"
   26   19  -54    1.51   10.78  0.6171E+00  0.6944E-01 "alien"
   20   37  -41    1.29   10.78  0.2182E+00  0.1157E+00 "alien"
    4   10   68    2.04   10.77  0.4383E+01  0.1465E+00 "alien"
   44   16  -15    1.23   10.75  0.1822E+00  0.7250E-01 "alien"
    4    9  -86    1.68   10.73  0.1625E+01  0.8813E-01 "alien"
    8    8    6    5.02   10.73  0.2528E+02  0.3197E+00 "alien"
    8    8   -6    5.02   10.72  0.2526E+02  0.3308E+00 "alien"
   37   27  -46    1.18   10.72  0.1293E+00  0.9271E-01 "alien"
   11    5  125    1.17   10.72  0.1293E+00  0.8624E-01 "alien"
    5    9  109    1.34   10.71  0.2573E+00  0.4503E-01 "alien"
   42   27    8    1.16   10.71  0.1292E+00  0.8882E-01 "alien"
   36   27   48    1.19   10.69  0.1528E+00  0.8728E-01 "alien"
   29   36  -37    1.20   10.69  0.1527E+00  0.1050E+00 "alien"
   36   19  -63    1.22   10.68  0.1809E+00  0.9609E-01 "alien"
    8    0   90    1.63   10.68  0.2438E+01  0.1449E+00 "alien"
    7   22   76    1.55   10.68  0.8330E+00  0.6257E-01 "alien"
   40   21  -54    1.16   10.67  0.1288E+00  0.8371E-01 "alien"
    4    5   98    1.51   10.67  0.6107E+00  0.6396E-01 "alien"
   11    1  127    1.15   10.65  0.1285E+00  0.9504E-01 "alien"
   29   25   73    1.22   10.65  0.1522E+00  0.8937E-01 "alien"
    5    3  114    1.31   10.65  0.2557E+00  0.5417E-01 "alien"
    3   13 -107    1.34   10.63  0.2553E+00  0.7436E-01 "alien"
   31   26   71    1.19   10.63  0.1519E+00  0.9531E-01 "alien"
    3    4  -95    1.57   10.62  0.8284E+00  0.7969E-01 "alien"
   17    5  116    1.20   10.60  0.1515E+00  0.6492E-01 "alien"
   19   43   -9    1.23   10.60  0.1796E+00  0.1323E+00 "alien"
   14    7  115    1.23   10.60  0.1796E+00  0.7108E-01 "alien"
   31   35  -26    1.21   10.60  0.1514E+00  0.9676E-01 "alien"
   35   31  -41    1.17   10.60  0.1278E+00  0.9225E-01 "alien"
    7   11  118    1.22   10.59  0.1794E+00  0.6902E-01 "alien"
    5   10  115    1.27   10.59  0.2144E+00  0.7722E-01 "alien"
    9    5   96    1.51   10.56  0.6047E+00  0.6347E-01 "alien"
   30   25  -62    1.27   10.56  0.2138E+00  0.1064E+00 "alien"
   34   30  -40    1.21   10.55  0.1507E+00  0.1025E+00 "alien"
   37   29   22    1.21   10.54  0.1506E+00  0.6501E-01 "alien"
    7    9  117    1.25   10.54  0.1785E+00  0.5608E-01 "alien"
    6   26  -89    1.33   10.54  0.2530E+00  0.7754E-01 "alien"
   43   14  -38    1.22   10.53  0.1505E+00  0.6053E-01 "alien"
   30    5   62    1.50   10.53  0.6027E+00  0.8009E-01 "alien"
   39   26   21    1.22   10.52  0.1503E+00  0.9348E-01 "alien"
   15   10  104    1.32   10.50  0.2520E+00  0.4859E-01 "alien"
    0    3  101    1.48   10.49  0.6004E+00  0.5822E-01 "alien"
    1    0 -104    1.44   10.49  0.9328E+00  0.9985E-01 "alien"
    3   21  -96    1.36   10.49  0.3400E+00  0.8457E-01 "alien"
   28   38  -33    1.19   10.48  0.1498E+00  0.1124E+00 "alien"
   11   20  108    1.22   10.48  0.1775E+00  0.6458E-01 "alien"
   18   29   55    1.44   10.48  0.4480E+00  0.4723E-01 "alien"
   18   12   96    1.35   10.47  0.3394E+00  0.5317E-01 "alien"
   29   39  -25    1.17   10.46  0.1261E+00  0.9234E-01 "alien"
    1   28   94    1.27   10.44  0.2114E+00  0.7015E-01 "alien"
   14   14  -90    1.45   10.43  0.4460E+00  0.6880E-01 "alien"
   37   10  -67    1.25   10.43  0.1766E+00  0.9251E-01 "alien"
   38   22  -54    1.19   10.42  0.1489E+00  0.8921E-01 "alien"
   26   14   82    1.34   10.41  0.2500E+00  0.6995E-01 "alien"
   13    4  102    1.39   10.41  0.3374E+00  0.4605E-01 "alien"
    8    5  119    1.24   10.41  0.1763E+00  0.4857E-01 "alien"
    7   16  121    1.16   10.40  0.1255E+00  0.9235E-01 "alien"
   17   19  101    1.25   10.40  0.1762E+00  0.7139E-01 "alien"
   38   23  -53    1.19   10.37  0.1251E+00  0.8609E-01 "alien"
   29   18   55    1.44   10.36  0.4430E+00  0.6284E-01 "alien"
    9    3 -101    1.45   10.36  0.4429E+00  0.6648E-01 "alien"
    4    0 -112    1.34   10.36  0.5134E+00  0.1026E+00 "alien"
    4    3  -95    1.57   10.36  0.8083E+00  0.6584E-01 "alien"
    1   11  -97    1.49   10.35  0.5926E+00  0.5629E-01 "alien"
   29   18  -55    1.44   10.35  0.4423E+00  0.6122E-01 "alien"
   15    8  114    1.23   10.35  0.1752E+00  0.9771E-01 "alien"
   14    8  119    1.19   10.33  0.1476E+00  0.9814E-01 "alien"
    6    7 -110    1.34   10.33  0.2481E+00  0.1014E+00 "alien"
    9    5  124    1.18   10.33  0.1246E+00  0.8433E-01 "alien"
   11    1  -97    1.49   10.31  0.5901E+00  0.7576E-01 "alien"
    9    4   90    1.61   10.31  0.1117E+01  0.7132E-01 "alien"
    0   14  113    1.27   10.30  0.2085E+00  0.8075E-01 "alien"
    9    3   87    1.66   10.28  0.1114E+01  0.6909E-01 "alien"
   12    2  124    1.17   10.27  0.1238E+00  0.5818E-01 "alien"
    5   30   62    1.50   10.25  0.5869E+00  0.6287E-01 "alien"
    9   31   98    1.17   10.25  0.1236E+00  0.7908E-01 "alien"
    4    3   95    1.57   10.25  0.7995E+00  0.6392E-01 "alien"
   23   21  -67    1.43   10.24  0.4377E+00  0.7138E-01 "alien"
   24    2  -70    1.60   10.24  0.1110E+01  0.7709E-01 "alien"
    6    2  -94    1.58   10.24  0.7988E+00  0.6793E-01 "alien"
    3    9  -87    1.66   10.23  0.1109E+01  0.7800E-01 "alien"
   15   16  101    1.30   10.23  0.2071E+00  0.5897E-01 "alien"
    1    0  100    1.50   10.23  0.1197E+01  0.1122E+00 "alien"
   11    5  110    1.31   10.22  0.2455E+00  0.4114E-01 "alien"
   14    6 -101    1.39   10.22  0.3312E+00  0.1189E+00 "alien"
   12   43   49    1.20   10.22  0.1460E+00  0.8614E-01 "alien"
    5    4   98    1.51   10.21  0.5845E+00  0.5631E-01 "alien"
   19   17  106    1.20   10.20  0.1458E+00  0.9024E-01 "alien"
   13    1  123    1.18   10.20  0.1231E+00  0.6509E-01 "alien"
   39   29   19    1.18   10.20  0.1231E+00  0.8273E-01 "alien"
   37   27   30    1.23   10.20  0.1728E+00  0.9731E-01 "alien"
    9   20  114    1.18   10.19  0.1230E+00  0.9210E-01 "alien"
   28   33  -47    1.24   10.19  0.1727E+00  0.9416E-01 "alien"
   32   28   63    1.18   10.19  0.1230E+00  0.8284E-01 "alien"
    1   21 -106    1.26   10.18  0.2061E+00  0.9411E-01 "alien"
   19   20  -68    1.52   10.18  0.7940E+00  0.8183E-01 "alien"
    1    0  104    1.44   10.15  0.8842E+00  0.9720E-01 "alien"
    5    5   90    1.64   10.15  0.1100E+01  0.7100E-01 "alien"
   40   22   21    1.25   10.14  0.1718E+00  0.8820E-01 "alien"
    2   28   93    1.27   10.14  0.2052E+00  0.6826E-01 "alien"
   17    2  -97    1.41   10.13  0.4330E+00  0.6738E-01 "alien"
    0   22   59    1.83   10.13  0.2342E+01  0.9520E-01 "alien"
   27   24   80    1.22   10.12  0.1715E+00  0.9860E-01 "alien"
   38   24  -41    1.22   10.11  0.1445E+00  0.6592E-01 "alien"
   13   13  110    1.25   10.11  0.1713E+00  0.5162E-01 "alien"
    3    5 -114    1.31   10.10  0.2426E+00  0.1071E+00 "alien"
   23    4   87    1.42   10.09  0.4314E+00  0.5893E-01 "alien"
   46    3  -48    1.17   10.09  0.1217E+00  0.7389E-01 "alien"
    0    4   78    1.91   10.09  0.4105E+01  0.1287E+00 "alien"
    5    5  -42    3.28   10.08  0.2345E+02  0.3494E+00 "alien"
    4   17 -103    1.34   10.08  0.2420E+00  0.7381E-01 "alien"
    3   14  112    1.27   10.07  0.2039E+00  0.6686E-01 "alien"
   42   14  -51    1.20   10.07  0.1439E+00  0.8741E-01 "alien"
   18   15  111    1.19   10.07  0.1439E+00  0.1042E+00 "alien"
   13    1  104    1.37   10.06  0.3260E+00  0.4406E-01 "alien"
   16   16    0    2.56   10.04  0.1113E+02  0.2183E+00 "alien"
    3    0 -102    1.47   10.03  0.1161E+01  0.1134E+00 "alien"
   23   22   91    1.22   10.03  0.1699E+00  0.9279E-01 "alien"
   12   43   59    1.16   10.03  0.1210E+00  0.8465E-01 "alien"
   29   10   32    1.75   10.01  0.1515E+01  0.9011E-01 "alien"
   10    2 -105    1.39   10.01  0.3244E+00  0.5896E-01 "alien"


 List of    457 reflections *NOT* obeying Wilson distribution (sorted by resolution)
 Ice rings could occur at (Angstrom):
 3.897,3.669,3.441, 2.671,2.249,2.072, 1.948,1.918,1.883,1.721

   h    k    l     RES      Z     Intensity    Sigma

   11    1  127    1.15   10.65  0.1285E+00  0.9504E-01
   24   44   -6    1.16   12.08  0.1457E+00  0.1140E+00
   13    5  125    1.16   10.81  0.1304E+00  0.8341E-01
   29   40  -22    1.16   10.89  0.1313E+00  0.9534E-01
   12   43   59    1.16   10.03  0.1210E+00  0.8465E-01
   34   33  -42    1.16   13.15  0.1587E+00  0.1019E+00
   27   26   86    1.16   11.99  0.1447E+00  0.9750E-01
   42   27    8    1.16   10.71  0.1292E+00  0.8882E-01
   15    8  122    1.16   16.11  0.1943E+00  0.1007E+00
    8   15  122    1.16   11.87  0.1432E+00  0.8682E-01
   45    6  -54    1.16   12.49  0.1506E+00  0.9368E-01
   25   43  -11    1.16   11.86  0.1431E+00  0.9353E-01
   38   26  -50    1.16   16.96  0.2046E+00  0.1141E+00
   11   40   72    1.16   12.70  0.1532E+00  0.8737E-01
   41   20  -52    1.16   11.58  0.1397E+00  0.8850E-01
    7   16  121    1.16   10.40  0.1255E+00  0.9235E-01
   11    0  126    1.16   14.19  0.4576E+00  0.1562E+00
   40   21  -54    1.16   10.67  0.1288E+00  0.8371E-01
   21   45    7    1.17   13.37  0.1612E+00  0.9293E-01
   28   28   78    1.17   10.82  0.1305E+00  0.9231E-01
    9   31   98    1.17   10.25  0.1236E+00  0.7908E-01
   27   41  -20    1.17   12.52  0.1511E+00  0.1088E+00
   11    5  125    1.17   10.72  0.1293E+00  0.8624E-01
   46    3  -48    1.17   10.09  0.1217E+00  0.7389E-01
   23   44    4    1.17   11.18  0.1348E+00  0.9431E-01
   20   19  107    1.17   12.31  0.1485E+00  0.8852E-01
    9   13  122    1.17   10.98  0.1325E+00  0.8152E-01
   26   42  -12    1.17   14.81  0.1787E+00  0.1007E+00
   11    4  125    1.17   13.07  0.1576E+00  0.1101E+00
   10   12  122    1.17   13.46  0.1623E+00  0.9332E-01
   46    5  -46    1.17   13.28  0.1601E+00  0.8640E-01
   29   39  -25    1.17   10.46  0.1261E+00  0.9234E-01
    8   31   98    1.17   10.97  0.1323E+00  0.8817E-01
    8   29  102    1.17   19.30  0.2328E+00  0.1029E+00
   13   42   59    1.17   12.40  0.1496E+00  0.8470E-01
   10    5  125    1.17   12.66  0.1527E+00  0.9969E-01
   23   24   95    1.17   15.20  0.1834E+00  0.9420E-01
   16    7  120    1.17   11.16  0.1347E+00  0.9880E-01
   41   27   15    1.17   11.64  0.1404E+00  0.9841E-01
    7   18  118    1.17   12.84  0.1549E+00  0.9090E-01
    7   25  109    1.17   17.42  0.2101E+00  0.9990E-01
   31   27   71    1.17   14.45  0.1743E+00  0.1054E+00
   38   23  -56    1.17   10.90  0.1315E+00  0.9915E-01
   10    4  125    1.17   24.40  0.2943E+00  0.1302E+00
   39   22  -54    1.17   11.36  0.1370E+00  0.8785E-01
   15    6  121    1.17   24.00  0.2895E+00  0.1212E+00
   35   31  -41    1.17   10.60  0.1278E+00  0.9225E-01
   43   13  -53    1.17   12.17  0.1468E+00  0.9359E-01
   12    2  124    1.17   10.27  0.1238E+00  0.5818E-01
   30   37  -33    1.18   14.84  0.1791E+00  0.1074E+00
    7   13  122    1.18   11.54  0.1393E+00  0.8306E-01
   45   17  -28    1.18   11.61  0.1401E+00  0.9528E-01
   10   15  119    1.18   12.78  0.1542E+00  0.8621E-01
   31   29   65    1.18   14.31  0.1727E+00  0.1006E+00
   19   13  113    1.18   13.09  0.1579E+00  0.1064E+00
   37   28  -43    1.18   15.56  0.1877E+00  0.1040E+00
   14    4  122    1.18   14.07  0.1697E+00  0.9787E-01
   20   17  108    1.18   14.62  0.1764E+00  0.1065E+00
    8   17  118    1.18   11.71  0.1413E+00  0.9737E-01
   22   19  103    1.18   19.43  0.2343E+00  0.1223E+00
   22   44   -3    1.18   13.10  0.1581E+00  0.1009E+00
   13   44   46    1.18   11.22  0.1353E+00  0.9122E-01
   13    1  123    1.18   10.20  0.1231E+00  0.6509E-01
   31   27   70    1.18   12.82  0.1547E+00  0.9119E-01
   19   15  111    1.18   11.48  0.1385E+00  0.9630E-01
   22   18  104    1.18   11.55  0.1394E+00  0.8798E-01
   39   29   19    1.18   10.20  0.1231E+00  0.8273E-01
   42   25   13    1.18   12.23  0.1475E+00  0.9052E-01
    9   20  114    1.18   10.19  0.1230E+00  0.9210E-01
   37   27  -46    1.18   10.72  0.1293E+00  0.9271E-01
   12   40   67    1.18   11.10  0.1339E+00  0.9201E-01
   36   30  -38    1.18   12.35  0.1490E+00  0.9608E-01
    8   13  121    1.18   11.94  0.1440E+00  0.9133E-01
   12    4  123    1.18   11.24  0.1356E+00  0.5594E-01
    9    2  125    1.18   18.96  0.2287E+00  0.1056E+00
   38   24  -51    1.18   14.63  0.1765E+00  0.9684E-01
   23   19  101    1.18   11.73  0.1415E+00  0.8388E-01
   11    6  123    1.18   11.94  0.1441E+00  0.8757E-01
   43   23  -13    1.18   11.07  0.1336E+00  0.8276E-01
   43    6  -59    1.18   12.55  0.1514E+00  0.8874E-01
   44    8  -52    1.18   17.95  0.2166E+00  0.1089E+00
   33   32  -44    1.18   13.11  0.1582E+00  0.9378E-01
   20   19  105    1.18   11.63  0.1403E+00  0.9871E-01
   41   19  -49    1.18   11.08  0.1337E+00  0.8926E-01
   32   28   63    1.18   10.19  0.1230E+00  0.8284E-01
   19   13  112    1.18   13.90  0.1677E+00  0.1167E+00
   34   32  -38    1.18   13.17  0.1588E+00  0.1119E+00
   13    7  121    1.18   17.65  0.2130E+00  0.1076E+00
   12    0  123    1.18   17.59  0.7301E+00  0.1892E+00
   25   42   -5    1.18   12.01  0.1448E+00  0.9001E-01
    9    5  124    1.18   10.33  0.1246E+00  0.8433E-01
   22   19  102    1.18   12.80  0.1544E+00  0.9356E-01
    7   24  109    1.19   17.18  0.2073E+00  0.9918E-01
   38   23  -53    1.19   10.37  0.1251E+00  0.8609E-01
   10    2  124    1.19   23.05  0.3293E+00  0.9700E-01
   28   38  -33    1.19   10.48  0.1498E+00  0.1124E+00
   30   29   66    1.19   13.58  0.1940E+00  0.1022E+00
   45   10  -42    1.19   11.02  0.1574E+00  0.8688E-01
   31   26   71    1.19   10.63  0.1519E+00  0.9531E-01
    5   28  103    1.19   11.22  0.1603E+00  0.9314E-01
   23   20   99    1.19   12.49  0.1785E+00  0.1035E+00
   28   25   81    1.19   14.15  0.2022E+00  0.9499E-01
   45   12  -38    1.19   11.45  0.1636E+00  0.9721E-01
   18   15  111    1.19   10.07  0.1439E+00  0.1042E+00
   26   24   87    1.19   10.99  0.1571E+00  0.9916E-01
   15    9  118    1.19   13.31  0.1902E+00  0.1090E+00
   13    8  120    1.19   14.45  0.2064E+00  0.1135E+00
    7    6  124    1.19   11.21  0.1601E+00  0.8611E-01
   21   17  105    1.19   12.00  0.1715E+00  0.9730E-01
    9    1  124    1.19   11.26  0.1610E+00  0.6602E-01
   36   27   48    1.19   10.69  0.1528E+00  0.8728E-01
    9   12  120    1.19   16.67  0.2382E+00  0.1095E+00
   41   25   20    1.19   14.15  0.2021E+00  0.9701E-01
   14    8  119    1.19   10.33  0.1476E+00  0.9814E-01
   41   26    6    1.19   10.99  0.1570E+00  0.8684E-01
    6   16  118    1.19   17.36  0.2480E+00  0.1081E+00
   38   22  -54    1.19   10.42  0.1489E+00  0.8921E-01
   26   41   -3    1.19   12.98  0.1855E+00  0.1143E+00
    5   20  114    1.19   11.64  0.1663E+00  0.9157E-01
   13    9  119    1.19   11.31  0.1616E+00  0.1019E+00
   18   45    5    1.19   15.57  0.2225E+00  0.1077E+00
   17   19  107    1.19   11.06  0.1580E+00  0.9188E-01
   18   14  111    1.19   11.15  0.1594E+00  0.9460E-01
   26   24   86    1.20   17.32  0.2475E+00  0.1099E+00
   12   43   49    1.20   10.22  0.1460E+00  0.8614E-01
   42   14  -51    1.20   10.07  0.1439E+00  0.8741E-01
   30   26   72    1.20   13.72  0.1960E+00  0.1028E+00
   43   22   -9    1.20   14.35  0.2051E+00  0.1067E+00
   32   27   63    1.20   13.59  0.1941E+00  0.1056E+00
   29   36  -37    1.20   10.69  0.1527E+00  0.1050E+00
    9    9  121    1.20   14.11  0.2017E+00  0.1002E+00
    8   13  119    1.20   15.50  0.2215E+00  0.1052E+00
   19   12  111    1.20   17.78  0.2540E+00  0.1154E+00
   41   10  -61    1.20   11.30  0.1615E+00  0.1107E+00
   29   30   63    1.20   14.61  0.2087E+00  0.9867E-01
   12    3  121    1.20   11.79  0.1685E+00  0.7167E-01
   24   24   89    1.20   13.25  0.1894E+00  0.9854E-01
   39   27   23    1.20   10.93  0.1562E+00  0.9930E-01
   17   45    9    1.20   11.50  0.1644E+00  0.9962E-01
   26   25   83    1.20   11.27  0.1610E+00  0.9066E-01
   11    5  121    1.20   11.27  0.1611E+00  0.5786E-01
   44   19  -13    1.20   11.09  0.1585E+00  0.8381E-01
   33   33  -31    1.20   14.09  0.2014E+00  0.1002E+00
   19   17  106    1.20   10.20  0.1458E+00  0.9024E-01
   10   22  108    1.20   14.63  0.2090E+00  0.7069E-01
    8   26  103    1.20   12.03  0.1718E+00  0.8829E-01
   10    2  122    1.20   11.23  0.1604E+00  0.5811E-01
   17    5  116    1.20   10.60  0.1515E+00  0.6492E-01
   31   33  -42    1.21   11.27  0.1610E+00  0.7536E-01
   15    4  118    1.21   12.14  0.1735E+00  0.9723E-01
   11    7  120    1.21   12.45  0.1779E+00  0.5676E-01
   41   14  -54    1.21   12.68  0.1811E+00  0.1078E+00
   42   12  -52    1.21   11.34  0.1620E+00  0.1100E+00
   12    5  120    1.21   13.04  0.1863E+00  0.7290E-01
    5   12  120    1.21   17.21  0.2459E+00  0.1143E+00
   20   19  102    1.21   13.17  0.1882E+00  0.1125E+00
    4   16  117    1.21   13.10  0.1871E+00  0.1040E+00
   25   24   86    1.21   13.51  0.1931E+00  0.1025E+00
   16   12  113    1.21   17.85  0.2550E+00  0.1322E+00
    9    2  122    1.21   12.78  0.1827E+00  0.5107E-01
    8   14  117    1.21   18.34  0.2621E+00  0.1126E+00
   30   31   54    1.21   14.90  0.2130E+00  0.1132E+00
   34   30  -40    1.21   10.55  0.1507E+00  0.1025E+00
    9    6  121    1.21   14.18  0.2026E+00  0.7910E-01
   26   40  -10    1.21   14.21  0.2030E+00  0.1095E+00
    8   19  112    1.21   14.84  0.2121E+00  0.9707E-01
   20   24   94    1.21   12.92  0.1846E+00  0.1023E+00
   31   35  -26    1.21   10.60  0.1514E+00  0.9676E-01
   33   25   62    1.21   11.21  0.1602E+00  0.6439E-01
   35   24  -57    1.21   12.41  0.1773E+00  0.1003E+00
   15   18  108    1.21   11.46  0.1637E+00  0.9652E-01
    4   14  118    1.21   10.82  0.1546E+00  0.9408E-01
   37   29   22    1.21   10.54  0.1506E+00  0.6501E-01
   10    6  120    1.21   25.69  0.3670E+00  0.7597E-01
    6   10  120    1.21   26.08  0.3726E+00  0.9279E-01
   38   24  -41    1.22   10.11  0.1445E+00  0.6592E-01
   44   17  -17    1.22   11.56  0.1651E+00  0.8772E-01
   22   12  105    1.22   13.01  0.1859E+00  0.1204E+00
   12   42   49    1.22   11.40  0.1630E+00  0.9504E-01
   17   43   29    1.22   11.01  0.1574E+00  0.9220E-01
   39   26   21    1.22   10.52  0.1503E+00  0.9348E-01
   43   14  -38    1.22   10.53  0.1505E+00  0.6053E-01
   17   22  100    1.22   11.37  0.1624E+00  0.9540E-01
   38   27   24    1.22   10.79  0.1542E+00  0.8734E-01
   29   25   73    1.22   10.65  0.1522E+00  0.8937E-01
   27   24   80    1.22   10.12  0.1715E+00  0.9860E-01
    8   15  115    1.22   12.19  0.2065E+00  0.9974E-01
   11   20  108    1.22   10.48  0.1775E+00  0.6458E-01
   17   43   28    1.22   11.45  0.1940E+00  0.1071E+00
   10    7  119    1.22   13.74  0.2327E+00  0.6202E-01
   32   29   51    1.22   11.36  0.1925E+00  0.7593E-01
   14   45   14    1.22   12.09  0.2048E+00  0.1091E+00
   36   19  -63    1.22   10.68  0.1809E+00  0.9609E-01
   23   22   91    1.22   10.03  0.1699E+00  0.9279E-01
   35   27  -44    1.22   12.14  0.2056E+00  0.1017E+00
    3   14  117    1.22   12.24  0.2073E+00  0.9722E-01
    7   11  118    1.22   10.59  0.1794E+00  0.6902E-01
    4    9  120    1.23   14.87  0.2519E+00  0.1235E+00
   27   38  -20    1.23   13.18  0.2233E+00  0.1156E+00
   37   27   30    1.23   10.20  0.1728E+00  0.9731E-01
   44   16  -15    1.23   10.75  0.1822E+00  0.7250E-01
   21   18   99    1.23   14.68  0.2487E+00  0.1123E+00
   19   43   -9    1.23   10.60  0.1796E+00  0.1323E+00
   15    8  114    1.23   10.35  0.1752E+00  0.9771E-01
   23   41   -5    1.23   12.47  0.2113E+00  0.1081E+00
   14    7  115    1.23   10.60  0.1796E+00  0.7108E-01
   15   11  112    1.23   13.70  0.2321E+00  0.7686E-01
   11    4  118    1.23   12.07  0.2045E+00  0.5396E-01
    4   11  118    1.23   14.92  0.2527E+00  0.8263E-01
   20   42  -17    1.23   13.78  0.2334E+00  0.1186E+00
   37   27   27    1.23   11.86  0.2010E+00  0.1029E+00
   11   34   79    1.23   14.11  0.2389E+00  0.7544E-01
   28   33  -47    1.24   10.19  0.1727E+00  0.9416E-01
   27   34  -46    1.24   12.02  0.2036E+00  0.1064E+00
   14   13  111    1.24   11.47  0.1942E+00  0.6495E-01
    3   17  113    1.24   15.97  0.2704E+00  0.7906E-01
    8    5  119    1.24   10.41  0.1763E+00  0.4857E-01
    5    8  119    1.24   13.93  0.2359E+00  0.7763E-01
   32   21   70    1.24   14.92  0.2527E+00  0.1228E+00
    7   12  116    1.24   18.18  0.3079E+00  0.8304E-01
   37   14  -65    1.24   13.10  0.2218E+00  0.1029E+00
    7    1  120    1.24   11.88  0.2013E+00  0.5088E-01
   18   15  105    1.24   12.35  0.2092E+00  0.7467E-01
   14   16  108    1.24   11.23  0.1903E+00  0.7997E-01
   36   16  -65    1.24   11.08  0.1878E+00  0.9565E-01
   44   15   -9    1.24   12.75  0.2160E+00  0.1041E+00
   17   19  101    1.25   10.40  0.1762E+00  0.7139E-01
    9    7  117    1.25   11.94  0.2022E+00  0.5000E-01
    7    9  117    1.25   10.54  0.1785E+00  0.5608E-01
    3   11  117    1.25   12.17  0.2061E+00  0.6843E-01
   19   40   37    1.25   11.24  0.1904E+00  0.1033E+00
   13   33   78    1.25   11.68  0.1978E+00  0.9793E-01
    6   11  116    1.25   13.32  0.2257E+00  0.7151E-01
    3    2  120    1.25   10.89  0.1844E+00  0.9024E-01
    3   12  116    1.25   11.35  0.1923E+00  0.7365E-01
    3   21  107    1.25   15.22  0.2578E+00  0.8235E-01
   23   40   12    1.25   11.18  0.1893E+00  0.9950E-01
   40   22   21    1.25   10.14  0.1718E+00  0.8820E-01
   31   28   52    1.25   12.07  0.2045E+00  0.6915E-01
   27   35  -36    1.25   19.57  0.3315E+00  0.9202E-01
    5   42   49    1.25   11.31  0.1915E+00  0.9999E-01
   13   13  110    1.25   10.11  0.1713E+00  0.5162E-01
   37   10  -67    1.25   10.43  0.1766E+00  0.9251E-01
    8    6  117    1.25   15.24  0.2581E+00  0.9032E-01
    6    8  117    1.25   12.64  0.2141E+00  0.6496E-01
   12   31   83    1.26   13.59  0.2302E+00  0.1000E+00
   38   25  -20    1.26   11.10  0.1881E+00  0.9100E-01
   10    5  116    1.26   14.02  0.2375E+00  0.5013E-01
    5   10  116    1.26   16.97  0.2875E+00  0.8219E-01
   32   19   70    1.26   12.40  0.2510E+00  0.1095E+00
    1   21 -106    1.26   10.18  0.2061E+00  0.9411E-01
   10    0  116    1.26   11.17  0.4916E+00  0.1014E+00
    0   10  116    1.26   19.75  0.3999E+00  0.1437E+00
    1   28   94    1.27   10.44  0.2114E+00  0.7015E-01
    7    7  116    1.27   13.39  0.2711E+00  0.6379E-01
   30   25  -62    1.27   10.56  0.2138E+00  0.1064E+00
    1   14 -113    1.27   12.28  0.2485E+00  0.1177E+00
    0   14  113    1.27   10.30  0.2085E+00  0.8075E-01
    5   10  115    1.27   10.59  0.2144E+00  0.7722E-01
   10   21  102    1.27   15.84  0.3207E+00  0.8564E-01
    0    9  116    1.27   19.60  0.3968E+00  0.1390E+00
   22   19   91    1.27   11.80  0.2388E+00  0.6622E-01
    2   28   93    1.27   10.14  0.2052E+00  0.6826E-01
    3   14  112    1.27   10.07  0.2039E+00  0.6686E-01
    2   15  111    1.28   11.97  0.2423E+00  0.7419E-01
   10   16  107    1.28   13.24  0.2681E+00  0.6219E-01
    4   20  105    1.28   11.08  0.2244E+00  0.8177E-01
   28   26  -63    1.28   16.86  0.3413E+00  0.1264E+00
    1    9 -115    1.28   12.04  0.2437E+00  0.1164E+00
   25   23   77    1.28   12.35  0.2500E+00  0.5616E-01
    0   13  112    1.28   11.95  0.2420E+00  0.7849E-01
    0   19  106    1.29   16.97  0.3435E+00  0.8561E-01
   20   37  -41    1.29   10.78  0.2182E+00  0.1157E+00
    5    5  115    1.29   11.46  0.2319E+00  0.5243E-01
    4    4  115    1.30   11.16  0.2259E+00  0.5109E-01
   16   15  101    1.30   15.28  0.3092E+00  0.5351E-01
   15   16  101    1.30   10.23  0.2071E+00  0.5897E-01
    4    7  114    1.30   11.33  0.2294E+00  0.5300E-01
   32   18   66    1.30   12.52  0.2535E+00  0.7729E-01
    6    3  114    1.30   12.24  0.2940E+00  0.5825E-01
    3    6  114    1.30   14.59  0.3504E+00  0.6221E-01
    5    3  114    1.31   10.65  0.2557E+00  0.5417E-01
    3    5 -114    1.31   10.10  0.2426E+00  0.1071E+00
   32   12  -73    1.31   12.41  0.2980E+00  0.8782E-01
    5   23  -97    1.31   11.03  0.2649E+00  0.1167E+00
   11    5  110    1.31   10.22  0.2455E+00  0.4114E-01
    5   11 -110    1.31   13.75  0.3302E+00  0.1187E+00
    5   11  110    1.31   11.02  0.2647E+00  0.4762E-01
   15   10  104    1.32   10.50  0.2520E+00  0.4859E-01
    6   26  -89    1.33   10.54  0.2530E+00  0.7754E-01
    6    7 -110    1.34   10.33  0.2481E+00  0.1014E+00
    4   17 -103    1.34   10.08  0.2420E+00  0.7381E-01
    4    0 -112    1.34   10.36  0.5134E+00  0.1026E+00
    0    4  112    1.34   23.18  0.5567E+00  0.7618E-01
    3   13 -107    1.34   10.63  0.2553E+00  0.7436E-01
   43    4   -8    1.34   18.05  0.4336E+00  0.1200E+00
    5    9  109    1.34   10.71  0.2573E+00  0.4503E-01
   26   14   82    1.34   10.41  0.2500E+00  0.6995E-01
    0   18  102    1.34   25.52  0.6129E+00  0.1134E+00
   18   12   96    1.35   10.47  0.3394E+00  0.5317E-01
    6    5 -109    1.36   11.86  0.3844E+00  0.7282E-01
   26    5   87    1.36   11.51  0.3733E+00  0.7715E-01
   16   16   94    1.36   12.02  0.3898E+00  0.5228E-01
    3   21  -96    1.36   10.49  0.3400E+00  0.8457E-01
    6   20  -96    1.36   16.24  0.5265E+00  0.1052E+00
   13    1  104    1.37   10.06  0.3260E+00  0.4406E-01
    1   13 -104    1.37   13.24  0.4292E+00  0.6640E-01
   13   12  -99    1.38   13.92  0.4512E+00  0.1395E+00
   22   12  -87    1.38   11.43  0.3705E+00  0.8891E-01
   16    3  100    1.38   14.47  0.4691E+00  0.5761E-01
    3   16  100    1.38   13.45  0.4361E+00  0.6193E-01
    3   16 -100    1.38   11.22  0.3638E+00  0.6348E-01
   14    6 -101    1.39   10.22  0.3312E+00  0.1189E+00
   10    2 -105    1.39   10.01  0.3244E+00  0.5896E-01
   10    0 -105    1.39   10.95  0.7603E+00  0.1232E+00
    0   10  105    1.39   23.86  0.7735E+00  0.7349E-01
   28    5  -79    1.39   12.75  0.4134E+00  0.6719E-01
   13    4  102    1.39   10.41  0.3374E+00  0.4605E-01
    4   13 -102    1.39   11.57  0.3751E+00  0.6462E-01
   17    2  -97    1.41   10.13  0.4330E+00  0.6738E-01
    2    2  106    1.41   11.96  0.5114E+00  0.5884E-01
    2    2 -106    1.41   14.08  0.6020E+00  0.8829E-01
   23    4   87    1.42   10.09  0.4314E+00  0.5893E-01
   17    6  -95    1.42   16.04  0.6856E+00  0.9482E-01
   10   10  -99    1.42   11.02  0.4711E+00  0.6876E-01
   21   21   72    1.43   13.32  0.5694E+00  0.5461E-01
   21   21  -72    1.43   13.96  0.5966E+00  0.8109E-01
   23   21  -67    1.43   10.24  0.4377E+00  0.7138E-01
    0   15   97    1.44   10.78  0.4609E+00  0.5787E-01
    2    2 -104    1.44   15.53  0.6639E+00  0.6153E-01
    2    2  104    1.44   12.28  0.5248E+00  0.5498E-01
   29   18   55    1.44   10.36  0.4430E+00  0.6284E-01
   29   18  -55    1.44   10.35  0.4423E+00  0.6122E-01
   18   29   55    1.44   10.48  0.4480E+00  0.4723E-01
    0    1  104    1.44   21.14  0.9038E+00  0.7733E-01
    1    0 -104    1.44   10.49  0.9328E+00  0.9985E-01
    1    0  104    1.44   10.15  0.8842E+00  0.9720E-01
    0    0  104    1.45   16.87  0.2348E+01  0.1685E+00
   16    3  -95    1.45   11.66  0.4983E+00  0.7274E-01
    3   16  -95    1.45   11.61  0.4963E+00  0.6010E-01
    9    3 -101    1.45   10.36  0.4429E+00  0.6648E-01
   14   14  -90    1.45   10.43  0.4460E+00  0.6880E-01
   14   13  -91    1.45   11.14  0.6379E+00  0.7905E-01
   13   14  -91    1.45   13.31  0.7618E+00  0.8725E-01
   12   12  -93    1.46   10.83  0.6198E+00  0.9405E-01
   24   11  -76    1.47   12.69  0.7264E+00  0.8390E-01
    3    0 -102    1.47   10.03  0.1161E+01  0.1134E+00
    0    3  102    1.47   23.58  0.1350E+01  0.8400E-01
    0   10   98    1.48   11.22  0.6424E+00  0.5883E-01
    0    3  101    1.48   10.49  0.6004E+00  0.5822E-01
   11    1  -97    1.49   10.31  0.5901E+00  0.7576E-01
    1   11  -97    1.49   10.35  0.5926E+00  0.5629E-01
   15   14  -86    1.49   11.29  0.6461E+00  0.9187E-01
    0   23   81    1.49   11.98  0.6856E+00  0.7864E-01
   20    5  -85    1.50   12.90  0.7384E+00  0.8606E-01
   30    5   62    1.50   10.53  0.6027E+00  0.8009E-01
   30    5  -62    1.50   15.04  0.8611E+00  0.7342E-01
    5   30   62    1.50   10.25  0.5869E+00  0.6287E-01
    1    0  100    1.50   10.23  0.1197E+01  0.1122E+00
    0    1  100    1.50   19.00  0.1088E+01  0.8109E-01
    0    0  100    1.50   11.94  0.1564E+01  0.1268E+00
    9    5   96    1.51   10.56  0.6047E+00  0.6347E-01
   26   19  -54    1.51   10.78  0.6171E+00  0.6944E-01
   19    3  -86    1.51   10.82  0.6196E+00  0.7858E-01
    5    4   98    1.51   10.21  0.5845E+00  0.5631E-01
    5    4  -98    1.51   11.52  0.6593E+00  0.6144E-01
    4    5  -98    1.51   10.90  0.6237E+00  0.6161E-01
    4    5   98    1.51   10.67  0.6107E+00  0.6396E-01
   19   20  -68    1.52   10.18  0.7940E+00  0.8183E-01
   16   16   79    1.53   12.12  0.9460E+00  0.6112E-01
   16   16  -79    1.53   12.08  0.9425E+00  0.7232E-01
   10   10   90    1.55   11.30  0.8816E+00  0.7324E-01
   10   10  -90    1.55   10.94  0.8534E+00  0.6728E-01
    4    4  -96    1.55   11.20  0.8742E+00  0.6964E-01
   22    7   76    1.55   11.27  0.8792E+00  0.5188E-01
   22    7  -76    1.55   12.62  0.9849E+00  0.7893E-01
    7   22   76    1.55   10.68  0.8330E+00  0.6257E-01
    7   22  -76    1.55   11.33  0.8843E+00  0.6825E-01
    7    4   94    1.56   14.07  0.1098E+01  0.8166E-01
    7    4  -94    1.56   13.92  0.1086E+01  0.7648E-01
    4    7   94    1.56   13.37  0.1043E+01  0.7185E-01
    4    7  -94    1.56   13.74  0.1072E+01  0.8191E-01
    4    3   95    1.57   10.25  0.7995E+00  0.6392E-01
    4    3  -95    1.57   10.36  0.8083E+00  0.6584E-01
    3    4  -95    1.57   10.62  0.8284E+00  0.7969E-01
    6    2  -94    1.58   10.24  0.7988E+00  0.6793E-01
   13   11   83    1.60   14.50  0.1572E+01  0.7994E-01
   13   11  -83    1.60   16.65  0.1805E+01  0.9641E-01
   11   13   83    1.60   14.68  0.1592E+01  0.8767E-01
   11   13  -83    1.60   16.63  0.1803E+01  0.9651E-01
   24    2  -70    1.60   10.24  0.1110E+01  0.7709E-01
    9    4   90    1.61   10.31  0.1117E+01  0.7132E-01
    4    9   90    1.61   11.57  0.1254E+01  0.8864E-01
    4    9  -90    1.61   11.88  0.1287E+01  0.8386E-01
    8    0  -90    1.63   11.83  0.2910E+01  0.1728E+00
    8    0   90    1.63   10.68  0.2438E+01  0.1449E+00
    0    8   90    1.63   23.38  0.2534E+01  0.1106E+00
    5    5   90    1.64   10.15  0.1100E+01  0.7100E-01
    5    5  -90    1.64   12.38  0.1341E+01  0.9985E-01
    9    3   87    1.66   10.28  0.1114E+01  0.6909E-01
    3    9  -87    1.66   10.23  0.1109E+01  0.7800E-01
   16    4   79    1.67   15.65  0.1697E+01  0.8317E-01
   16    4  -79    1.67   14.16  0.1535E+01  0.1348E+00
    4   16   79    1.67   13.80  0.1496E+01  0.8084E-01
    4   16  -79    1.67   12.56  0.1362E+01  0.9754E-01
   14    2  -82    1.67   11.32  0.1228E+01  0.1229E+00
    9    4   86    1.68   10.80  0.1636E+01  0.8116E-01
    4    9  -86    1.68   10.73  0.1625E+01  0.8813E-01
   17    4   75    1.72   12.08  0.1830E+01  0.7985E-01
   17    4  -75    1.72   12.36  0.1872E+01  0.9702E-01
    9    9  -81    1.72   26.86  0.4067E+01  0.1403E+00
    9    9   81    1.72   23.41  0.3545E+01  0.1206E+00
   29   10  -32    1.75   11.95  0.1809E+01  0.1170E+00
   29   10   32    1.75   10.01  0.1515E+01  0.9011E-01
    1    1   83    1.81   11.40  0.2636E+01  0.1045E+00
    1    1  -83    1.81   11.21  0.2593E+01  0.1053E+00
    0   22   59    1.83   10.13  0.2342E+01  0.9520E-01
    9    9   75    1.83   12.28  0.2840E+01  0.1162E+00
    9    9  -75    1.83   12.48  0.2886E+01  0.1215E+00
    7    7   76    1.87   12.35  0.2855E+01  0.1079E+00
    7    7  -76    1.87   13.39  0.3096E+01  0.1357E+00
    0   21   58    1.89   12.89  0.2980E+01  0.1256E+00
    0    4   78    1.91   10.09  0.4105E+01  0.1287E+00
    0    5   75    1.97   11.84  0.4818E+01  0.1388E+00
   20    5   54    1.98   12.01  0.4887E+01  0.1122E+00
   20    5  -54    1.98   11.44  0.4657E+01  0.1414E+00
    5   20   54    1.98   11.12  0.4528E+01  0.1278E+00
    5   20  -54    1.98   11.56  0.4704E+01  0.1343E+00
   11    9   65    2.01   16.66  0.6782E+01  0.1516E+00
   11    9  -65    2.01   17.34  0.7057E+01  0.1758E+00
    9   11   65    2.01   16.04  0.6530E+01  0.1522E+00
    9   11  -65    2.01   16.77  0.6826E+01  0.1713E+00
   10    4   68    2.04   11.16  0.4542E+01  0.1185E+00
   10    4  -68    2.04   11.22  0.4568E+01  0.1369E+00
    4   10   68    2.04   10.77  0.4383E+01  0.1465E+00
    4   10  -68    2.04   11.20  0.4557E+01  0.1331E+00
    0    5   71    2.08   12.48  0.7832E+01  0.1809E+00
    4    3  -68    2.17   12.96  0.8134E+01  0.2170E+00
    4    3   68    2.17   11.78  0.7393E+01  0.1656E+00
    3    4   68    2.17   12.81  0.8038E+01  0.1760E+00
    3    4  -68    2.17   13.43  0.8431E+01  0.2480E+00
    0   12   55    2.38   12.34  0.1014E+02  0.1996E+00
    0    7   60    2.40   11.43  0.9397E+01  0.1916E+00
    0   17   43    2.44   11.30  0.9287E+01  0.2159E+00
   16   16    0    2.56   10.04  0.1113E+02  0.2183E+00
   15   15    0    2.73   12.47  0.1382E+02  0.2194E+00
   14   14   14    2.82   12.06  0.1337E+02  0.2205E+00
   14   14  -14    2.82   12.71  0.1409E+02  0.2310E+00
    5    5  -42    3.28   10.08  0.2345E+02  0.3494E+00
    8    6   31    3.72   12.35  0.4935E+02  0.6091E+00
    8    6  -31    3.72   12.03  0.4808E+02  0.6278E+00
    6    8   31    3.72   12.24  0.4891E+02  0.7067E+00
    6    8  -31    3.72   12.19  0.4872E+02  0.6315E+00
    0    8   27    4.41   15.29  0.6112E+02  0.7588E+00
    0   11   10    4.97   12.24  0.2885E+02  0.3734E+00
    8    8    6    5.02   10.73  0.2528E+02  0.3197E+00
    8    8   -6    5.02   10.72  0.2526E+02  0.3308E+00

 NUMBER OF UNIQUE ALIEN REFLECTIONS WITH A Z-SCORE ABOVE LIMIT        13
 (ALIENS ABOVE LIMIT (REJECT_ALIEN=      20.0) ARE MARKED INVALID)

 NUMBER OF REFLECTION RECORDS ON OUTPUT FILE "XDS_ASCII.HKL"     1503511
 NUMBER OF ACCEPTED OBSERVATIONS (INCLUDING SYSTEMATIC ABSENCES) 1497358
 NUMBER OF REJECTED MISFITS & ALIENS (marked by -1*SIGMA(IOBS))     6153


 cpu time used                559.5 sec
 elapsed wall-clock time       18.5 sec
