*----
*  TEST CASE TCM02
*  MACROSCOPIC CROSS SECTIONS
*  FIXED SOURCE PROBLEM
*  FOR 1/8 7X7 PWR ASSEMBLY
*
*  REF: Z. Stankovski, Nucl. Sci. Eng. 92, 255 (1986)
*       R. Roy et al. Advances in Mathematics, Computation
*       and Reactor Physics, April 28 - May 2 1991, Pittsburgh
*
*----
*  Define STRUCTURES and MODULES used
*----
LINKED_LIST 
  PWR TRACK MACRO ;
SEQ_ASCII
  Fig.ps ;
MODULE 
  GEO: EXCELT: MAC: END: PSP: ;
*----
* Macroscopic XS
*----
MACRO := MAC: ::
  NGRO 1 NMIX 20
  READ INPUT
  MIX 1 TOTAL 14.000  SCAT 1 1  0.000  FIXE 0.000
  MIX 2 TOTAL  1.250  SCAT 1 1  1.242  FIXE 1.000
  MIX 3 TOTAL  0.625  SCAT 1 1  0.355  FIXE 0.000
  MIX 4 TOTAL  1.250  SCAT 1 1  1.242  FIXE 1.000
  MIX 5 TOTAL  0.625  SCAT 1 1  0.355  FIXE 0.000
  MIX 6 TOTAL  1.250  SCAT 1 1  1.242  FIXE 1.000
  MIX 7 TOTAL  0.625  SCAT 1 1  0.355  FIXE 0.000
  MIX 8 TOTAL  1.250  SCAT 1 1  1.242  FIXE 1.000
  MIX 9 TOTAL  0.625  SCAT 1 1  0.355  FIXE 0.000
  MIX 10 TOTAL  1.250  SCAT 1 1  1.242  FIXE 1.000
  MIX 11 TOTAL  0.625  SCAT 1 1  0.355  FIXE 0.000
  MIX 12 TOTAL  1.250  SCAT 1 1  1.242  FIXE 1.000
  MIX 13 TOTAL  0.625  SCAT 1 1  0.355  FIXE 0.000
  MIX 14 TOTAL  1.250  SCAT 1 1  1.242  FIXE 1.000
  MIX 15 TOTAL  0.625  SCAT 1 1  0.355  FIXE 0.000
  MIX 16 TOTAL  1.250  SCAT 1 1  1.242  FIXE 1.000
  MIX 17 TOTAL  0.625  SCAT 1 1  0.355  FIXE 0.000
  MIX 18 TOTAL  1.250  SCAT 1 1  1.242  FIXE 1.000
  MIX 19 TOTAL  0.625  SCAT 1 1  0.355  FIXE 0.000
  MIX 20 TOTAL  1.250  SCAT 1 1  1.242  FIXE 1.000
  ;
*----
*  Geometry : PWR - Cartesian 4X4
*  Tracking : EXCELT
*----
PWR := GEO: :: CAR2D 4 4
  X- DIAG  X+ REFL Y- SYME  Y+ DIAG
  CELL   P F1 F2 F3
           F4 F5 F6
              F7 F8
                 F9 
  ::: F1 := GEO: CARCEL 1
    MIX 3  4
    RADIUS 0.000  0.450 SPLITR 4
    MESHX -0.625  0.625 SPLITX 8
    MESHY -0.625  0.625 SPLITY 8          
    ;
  ::: P := GEO: F1
    MIX  1  2  SPLITR 5 
    ;
  ::: F2 := GEO: F1
    MIX  5  6
    ;
  ::: F3 := GEO: F1
    MIX  7  8
    ;
  ::: F4 := GEO: F1
    MIX  9  10
    ;
  ::: F5 := GEO: F1
    MIX  11  12
    ;
  ::: F6 := GEO: F1
    MIX  13  14
    ;
  ::: F7 := GEO: F1
    MIX  15  16
    ;
  ::: F8 := GEO: F1
    MIX  17  18
    ;
  ::: F9 := GEO: F1
    MIX  19  20
    ;
  ;
*----
*   Plotting
*----
Fig.ps := PSP: PWR :: TYPE REGION ;
Fig.ps := PSP: Fig.ps PWR :: TYPE MIXTURE ;
END: ;
QUIT "LIST" .
