      subroutine aaosub_1(epq2,epw,T0,epcos,epphi,epirea,sigma0,sigu
     *  ,sigt,sigl,sigi)
*************************************************************
*                                                           *
*                Amplitudes and Observables for the         *
*                    CLAS  N* Monte Carlo Program           *
*                                                           *
*                   Author: V.D. Burkert                    *
*                    Begin: February 1991                   *
*                                                           *
*                 Modified:                                 *
*                                                           *
*************************************************************
*                                                           *
**********************************************************************
*    All the amplitudes are multiplied by the factor 19.73           *
*    to convert the cross section from inverse GeV**2 to microbarns. *
*                                                                    *
**********************************************************************
*                                                                    *
*     e0=Beam Energy (GeV)                                           *
*     epQ2=negative of four-momentum transfer (GEV**2)')             *
*                                                                    *
**********************************************************************

       implicit none

         COMMON /ANREG2/
     * BGA0P0,BGA1P0,BGA2P0,BGA3P0,
     * BGA1M0,BGA2M0,BGA3M0,
     * BGB1P0,BGB2P0,BGB3P0,
     * BGB2M0,BGB3M0,
     * BGC0P0,BGC1P0,BGC2P0,BGC3P0,
     * BGC1M0,BGC2M0,BGC3M0,
     * BGA0PP,BGA1PP,BGA2PP,BGA3PP,
     * BGA1MP,BGA2MP,BGA3MP,
     * BGB1PP,BGB2PP,BGB3PP,
     * BGB2MP,BGB3MP,
     * BGC0PP,BGC1PP,BGC2PP,BGC3PP,
     * BGC1MP,BGC2MP,BGC3MP,
     * RAS11_1535,             RCS11_1535,
     * RAS11_1650,             RCS11_1650,
     * RAP11_1440,             RCP11_1440,
     * RAP11_1710,             RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,
     * RAP13_1910, RBP13_1910, RCP13_1910,
     * RAD13_1520, RBD13_1520, RCD13_1520,
     * RAD13_1700, RBD13_1700, RCD13_1700,
     * RAD15_1675, RBD15_1675, RCD15_1675,
     * RAF15_1680, RBF15_1680, RCF15_1680,
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,             RCS31_1620,
     * RAS31_1900,             RCS31_1900,
     * RAP31_1910,             RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232,
     * RAP33_1920, RBP33_1920, RCP33_1920,
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420,
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710,
     * WP13_1720, WD13_1520, WD13_1700,
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600,
     * WS31_1620, WS31_1900, WP31_1910,
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710,
     * LP13_1720, LD13_1520, LD13_1700,
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600,
     * LS31_1620, LS31_1900, LP31_1910,
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420



      common /debug/ibug

      real sig0, e0, epq2, epw, sigma0
      real sigu,sigt,sigl,sigi,fkt
      real x1_0,xu_0,xt_0,xl_0,xi_0
      real xmpi, elchar, enue
      real e_elec_scatt, theta_elec, epeps, q2evf, q2evf0
      real dip_evf,qme11,qme22,qmm11,qmm12,qmm21,qmm22,qmm23
      real fact, pvec_gam, pvec_pi,pvec_gam0,q0_gam
      real epcos,epphi,angfac,ffac
      REAL PI, XMPROT,XMPIP,XMPI0,THETA_S11
      real T0,s2
      real sigmao
      real si,st
      real asind,sind,cosd,tand
      integer ibug
      integer sigtot,icos
      integer epirea

      REAL*4
     * BGA0P0,BGA1P0,BGA2P0,BGA3P0,
     * BGA1M0,BGA2M0,BGA3M0,
     * BGB1P0,BGB2P0,BGB3P0,
     * BGB2M0,BGB3M0,
     * BGC0P0,BGC1P0,BGC2P0,BGC3P0,
     * BGC1M0,BGC2M0,BGC3M0,
     * BGA0PP,BGA1PP,BGA2PP,BGA3PP,
     * BGA1MP,BGA2MP,BGA3MP,
     * BGB1PP,BGB2PP,BGB3PP,
     * BGB2MP,BGB3MP,
     * BGC0PP,BGC1PP,BGC2PP,BGC3PP,
     * BGC1MP,BGC2MP,BGC3MP

      REAL*4
     * PIBS11_1535, PIBS11_1650, PIBP11_1440, PIBP11_1710,
     * PIBP13_1720, PIBD13_1520, PIBD13_1700,
     * PIBD15_1675, PIBF15_1680,
c     * PIBG17_2190,
c     * PIBG19_2250, PIBH19_2220,
c     * PIBI111_2600,
     * PIBS31_1620, PIBS31_1900, PIBP31_1910,
     * PIBP33_1232, PIBP33_1920, PIBD33_1700,
     * PIBD35_1930, PIBF35_1905,
     * PIBF37_1950,
     * PIBH311_2420

      REAL
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710,
     * WP13_1720, WD13_1520, WD13_1700,
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600,
     * WS31_1620, WS31_1900, WP31_1910,
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420

      REAL
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710,
     * LP13_1720, LD13_1520, LD13_1700,
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600,
     * LS31_1620, LS31_1900, LP31_1910,
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420

      real
     * RAS11_1535,             RCS11_1535,
     * RAS11_1650,             RCS11_1650,
     * RAP11_1440,             RCP11_1440,
     * RAP11_1710,             RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,
     * RAP13_1910, RBP13_1910, RCP13_1910,
     * RAD13_1520, RBD13_1520, RCD13_1520,
     * RAD13_1700, RBD13_1700, RCD13_1700,
     * RAD15_1675, RBD15_1675, RCD15_1675,
     * RAF15_1680, RBF15_1680, RCF15_1680,
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,             RCS31_1620,
     * RAS31_1900,             RCS31_1900,
     * RAP31_1910,             RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232,
     * RAP33_1920, RBP33_1920, RCP33_1920,
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420

       real
     * a12p33_1232,a32p33_1232, s12p33_1232,
     * a12p11_1440,s12p11_1440,
     * a12d13_1520, a32d13_1520, s12d13_1520,
     * A12S11_1535,s12s11_1535,
     * a12s31_1620,s12s31_1620,
     * A12S11_1650,s12s11_1650,
     * a12d15_1675,a32d15_1675,s12d15_1675,a33d15_1675,
     * a12f15_1680,a32f15_1680,s12f15_1680,
     * a12d13_1700, a32d13_1700,s12d13_1700,
     * a12d33_1700,a32d33_1700, s12d33_1700,
     * a12p11_1710,s12p11_1710,
     * a12p13_1720, a32p13_1720,s12p13_1720,
     * a12f35_1905, a32f35_1905, s12f35_1905,
     * a12s31_1900,s12s31_1900,
     * a12p31_1910,s12p31_1910,
     * a12p33_1920,a32p33_1920, s12p33_1920,
     * a12d35_1930, a32d35_1930, s12d35_1930,
     * a12f37_1950, a32f37_1950, s12f37_1950,
     * a12g17_2190,a32g17_2190,s12g17_2190

      DATA A12S11_1535,              S12S11_1535  /0.  ,0.   /
      DATA A12S11_1650,              S12S11_1650  /0.  ,0.   /
      DATA A12P11_1440,              S12P11_1440  /-0.065 ,0. /
      DATA A12P11_1710,              S12P11_1710  /0.  ,0.   /
      DATA A12P13_1720, A32P13_1720, S12P13_1720  /0.  ,0.   ,0./
c      DATA A12P13_1910, A32P13_1910, S12P13_1910  /0.  ,0.   ,0./
      DATA A12D13_1520, A32D13_1520, S12D13_1520  /0.  ,0.   ,0./
      DATA A12D13_1700, A32D13_1700, S12D13_1700  /0.  ,0.   ,0./
      DATA A12D15_1675, A32D15_1675, S12D15_1675  /0.  ,0.   ,0./
      DATA A12F15_1680, A32F15_1680, S12F15_1680  /0.  ,0.   ,0./
      DATA A12G17_2190, A32G17_2190, S12G17_2190  /0.  ,0.   ,0./
c      DATA A12G19_2250, A32G19_2250, S12G19_2250  /0.  ,0.   ,0./
c      DATA A12H19_2220, A32H19_2220, S12H19_2220  /0.  ,0.   ,0./
c      DATA A12I111_2600,A32I111_2600,S12I111_2600 /0.  ,0.   ,0./
      DATA A12S31_1620,              S12S31_1620  /0.  ,0.   /
      DATA A12S31_1900,              S12S31_1900  /0.  ,0.   /
      DATA A12P31_1910,              S12P31_1910  /0.  ,0.   /
      DATA A12P33_1232, A32P33_1232, S12P33_1232  /-0.15,-0.27,0./
      DATA A12P33_1920, A32P33_1920, S12P33_1920  /0.  ,0.   ,0.   /
      DATA A12D33_1700, A32D33_1700, S12D33_1700  /0.  ,0.   ,0.   /
      DATA A12D35_1930, A32D35_1930, S12D35_1930  /0.  ,0.   ,0.   /
      DATA A12F35_1905, A32F35_1905, S12F35_1905  /0.  ,0.   ,0.   /
      DATA A12F37_1950, A32F37_1950, S12F37_1950  /0.  ,0.   ,0.   /
c      DATA A12H311_2420,A32H311_2420,S12H311_2420 /0.  ,0.   ,0.   /

      DATA LS11_1535      /1.532    /
      DATA LS11_1650      /1.650    /
      DATA LP11_1440      /1.440    /
      DATA LP11_1710      /1.710    /
      DATA LP13_1720      /1.740    /
      DATA LD13_1520      /1.520    /
      DATA LD13_1700      /1.700    /
      DATA LD15_1675      /1.675    /
      DATA LF15_1680      /1.688    /
      DATA LG17_2190      /2.175    /
      DATA LG19_2250      /2.250    /
      DATA LH19_2220      /2.220    /
      DATA LI111_2600     /2.640    /
      DATA LS31_1620      /1.620    /
      DATA LS31_1900      /1.925    /
      DATA LP31_1910      /1.905    /
      DATA LP33_1232      /1.232    /
      DATA LP33_1920      /1.920    /
      DATA LD33_1700      /1.670    /
      DATA LD35_1930      /1.930    /
      DATA LF35_1905      /1.905    /
      DATA LF37_1950      /1.950    /
      DATA LH311_2420     /2.420    /

      DATA WS11_1535      /.150     /
      DATA WS11_1650      /.150     /
      DATA WP11_1440      /.300     /
      DATA WP11_1710      /.110     /
      DATA WP13_1720      /.200     /
      DATA WD13_1520      /.125     /
      DATA WD13_1700      /.100     /
      DATA WD15_1675      /.155     /
      DATA WF15_1680      /.125     /
      DATA WG17_2190      /.350     /
      DATA WG19_2250      /.300     /
      DATA WH19_2220      /.400     /
      DATA WI111_2600     /.400     /
      DATA WS31_1620      /.140     /
      DATA WS31_1900      /.150     /
      DATA WP31_1910      /.220     /
      DATA WP33_1232      /.115     /
      DATA WP33_1920      /.250     /
      DATA WD33_1700      /.250     /
      DATA WD35_1930      /.250     /
      DATA WF35_1905      /.300     /
      DATA WF37_1950      /.240     /
      DATA WH311_2420     /.300     /

      DATA THETA_S11 / 38. /
c      DATA THETA_D13 /  0. /
c      DATA THETA_D15 /  0. /

c      DATA SQ2        /1.41421  /
      DATA PI         /3.1415926  /
      DATA XMPROT     /.938     /
      DATA XMPIP      /.1395    /
      DATA XMPI0      /.1349    /
c      DATA XMETA      /.5488    /
C
C      Branching ratios into N-pion channel
C
      DATA PIBS11_1535      /0.40    /
      DATA PIBS11_1650      /0.60    /
      DATA PIBP11_1440      /0.60    /
      DATA PIBP11_1710      /0.15    /
      DATA PIBP13_1720      /0.15    /
      DATA PIBD13_1520      /0.55    /
      DATA PIBD13_1700      /0.10    /
      DATA PIBD15_1675      /0.40    /
      DATA PIBF15_1680      /0.60    /
c      DATA PIBG17_2190      /0.14    /
c      DATA PIBG19_2250      /0.10    /
c      DATA PIBH19_2220      /0.18    /
c      DATA PIBI111_2600     /0.05    /

      DATA PIBS31_1620      /0.30    /
      DATA PIBS31_1900      /0.10    /
      DATA PIBP31_1910      /0.20    /
      DATA PIBP33_1232      /0.994   /
      DATA PIBP33_1920      /0.20    /
      DATA PIBD33_1700      /0.15    /
      DATA PIBD35_1930      /0.10    /
      DATA PIBF35_1905      /0.10    /
      DATA PIBF37_1950      /0.40    /
      DATA PIBH311_2420     /0.10    /


       PI=3.1415926
       X1_0=0.
       sigma0=0.
c
c     set up the kinematics from Q**2, W, and theta_elec
         ENUE = (EPW**2 - XMPROT**2 + EPQ2)/2./XMPROT
         theta_elec=180.*T0/pi
         s2=sin(T0/2.)**2
         E0=(ENUE+sqrt(epq2/s2+enue**2))/2.
         E_ELEC_SCATT = E0 - ENUE
         if ( E0 .lt. 0.1 .or. E_ELEC_SCATT .lt. 0.01)then
            write(6,*)' aaosub:, e0,e_elec_scatt',e0,e_elec_scatt
            sigma0=0.1e-20
            return
         endif
       
C
C       Now calculate the resonant helicity amplitudes
C
         IF(EPIREA.EQ.1)XMPI=XMPI0
         IF(EPIREA.EQ.3)XMPI=XMPIP
          ELCHAR = SQRT(4.*PI/137.)
         if (epq2 .le. 0.)then
            write(6,*)' aaosub: epq2',epq2
            sigma0=0.1e-26
            return
         endif
         EPEPS = 1./(1.+ 2.*(1.+ENUE**2/EPQ2)*(TAND(THETA_ELEC/2.))**2.)
c         GAMFLX = 1./137./4./PI**2/XMPROT/EPQ2/E0/
c     *            (1.-EPEPS)*(EPW**2 - XMPROT**2)*E_ELEC_SCATT

c           THETA_GAM = ATAND(SIND(THETA_ELEC)/
c     *                 (E0/E_ELEC_SCATT - COSD(THETA_ELEC)))
         Q2EVF=((EPW**2-XMPROT**2)**2+EPQ2*(EPW+XMPROT)**2)/
     *   (4.*EPW*XMPROT)
         Q2EVF0=((EPW**2-XMPROT**2)**2)/(4.*EPW*XMPROT)
         if (q2evf .le. 0. .or. q2evf0 .le. 0.)then
            write(6,*)' aaosub: q2evf, q2evf0',q2evf, q2evf0
            sigma0=0.1e-28
            return
         endif

          DIP_EVF = (1.+ Q2EVF/0.71)**(-2.)
            QME11 = 3.* DIP_EVF
            QMM11 =  3.8 * (0.4 - Q2EVF) * DIP_EVF
         IF(Q2EVF.LE.1.0)THEN
            QMM12 = 5.* Q2EVF * DIP_EVF
          ELSE
            QMM12 = ( 6.0 - 1.0 * Q2EVF ) * DIP_EVF
         ENDIF

        A12S11_1535 = ELCHAR/2./SQRT(XMPROT*(LS11_1535**2-XMPROT**2))*
     *       (SQRT(1./3.)*QME11 - SQRT(2./3.)*QMM11)*COSD(THETA_S11)
        A12S11_1650 = ELCHAR/2./SQRT(XMPROT*(LS11_1650**2-XMPROT**2)) *
     *        (SQRT(1./3.)*QME11 - SQRT(2./3.)*QMM11)*SIND(THETA_S11)

C
C        Assume no mixing for D13 and D15
C
        A12D13_1520=ELCHAR/2./SQRT(XMPROT*(LD13_1520**2-XMPROT**2))*
     *   (SQRT(1./6.)*QME11 + SQRT(1./12.)*QMM11 - SQRT(3./4.)*QMM12)
        A32D13_1520=ELCHAR/2./SQRT(XMPROT*(LD13_1520**2-XMPROT**2))*
     *           (SQRT(1./2.)*QME11 + 0.5*QMM11 + 0.5*QMM12)


         A12D13_1700 = 0.
         A32D13_1700 = 0.

         A12D15_1675 = 0.
         A33D15_1675 = 0.

        A12D33_1700 = -ELCHAR/2./SQRT(XMPROT*(LD33_1700**2-XMPROT**2))*
     *   (-SQRT(1./12.) * QMM12 - SQRT(1./6.) * QME11 +
     *    0.5*SQRT(1./27.)*QMM11)
        A32D33_1700 = -ELCHAR/2./SQRT(XMPROT*(LD33_1700**2-XMPROT**2))*
     *        (1./6.*QMM12 - SQRT(1./2.)*QME11 + 1./6.*QMM11)

        A12S31_1620 = ELCHAR/2./SQRT(XMPROT*(LS31_1620**2-XMPROT**2))*
     *        ( SQRT(1./3.)*QME11 + SQRT(2./27.)*QMM11)

C
C         Now do the {56,0+} ---> {56,2+} transition.
C         Here we use the quark multipole moments of Cottingham & Dunbar
C         directly. Note that quark multipole QMM21 is assumed = 0.
C         This assumption has no justification (data do not allow its
C         determination).
C
          QME22 = 0.99*SQRT(Q2EVF)*DIP_EVF
          QMM21 = 0.* SQRT(Q2EVF)*DIP_EVF
          QMM22 = (0.75*SQRT(5.)-1.5*SQRT(5.)*Q2EVF)*
     *            SQRT(Q2EVF)*DIP_EVF
          IF(Q2EVF.LE.1.0)THEN
          QMM23 = 5.*Q2EVF*SQRT(Q2EVF)*DIP_EVF
          ELSE
          QMM23 = (5.7-0.7*Q2EVF)*SQRT(Q2EVF)*DIP_EVF
          ENDIF

          FACT = ELCHAR/2./SQRT(XMPROT*(LF15_1680**2 - XMPROT**2))
          A12F15_1680 = FACT*(-SQRT(2./3.)*QMM23 + SQRT(1./5.)*QME22-
     *                  SQRT(2./15.)*QMM22)
          A32F15_1680 = FACT*(SQRT(1./3.)*QMM23 + SQRT(2./5.)*QME22+
     *                  SQRT(4./15.)*QMM22)

          FACT = ELCHAR/2./SQRT(XMPROT*(LP13_1720**2 - XMPROT**2))
          A12P13_1720 = -FACT*(SQRT(3./4.)*QMM21 + SQRT(1./10.)*QME22-
     *                   SQRT(3./20.)*QMM22)
          A32P13_1720 =  FACT*(0.5*QMM21 - SQRT(3./10.)*QME22 +
     *                   SQRT(9./20.)*QMM22)

          FACT = ELCHAR/2./SQRT(XMPROT*(LF35_1905**2 - XMPROT**2))
          A12F35_1905 =  FACT*(SQRT(16./189.)*QMM23 +
     *                  SQRT(28./135.)*QMM22)
          A32F35_1905 = -FACT*(SQRT(8./189.)*QMM23 -
     *                  SQRT(56./135.)*QMM22)

          FACT = ELCHAR/2./SQRT(XMPROT*(LP33_1920**2 - XMPROT**2))
          A12P33_1920 = -FACT*(SQRT(3./15.)*QMM22 + 1./3.*QMM21)
          A32P33_1920 =  FACT*(SQRT(1./15.)*QMM22 - SQRT(1./3.)*QMM21)

          FACT = ELCHAR/2./SQRT(XMPROT*(LF37_1950**2 - XMPROT**2))
          A12F37_1950 = -FACT*SQRT(2./7.)*QMM23
          A32F37_1950 = -FACT*SQRT(10./21.)*QMM23

          FACT = ELCHAR/2./SQRT(XMPROT*(LP31_1910**2 - XMPROT**2))
          A12P31_1910 = -FACT*3./2.*QMM21
C
         BGA0P0 =            0.0000D+00

         BGA1P0 =            0.0D+00
         BGA2P0 =            0.0D+00
ccc         BGA1P0 =            0.38
ccc         BGA2P0 =            0.64
         BGA1P0 =            0.
         BGA2P0 =            0.
         BGA3P0 =            0.0000D+00

         BGB1P0 =            0.0D+00
         BGB2P0 =            0.0D+00
ccc         BGB1P0 =            0.63
ccc         BGB2P0 =            0.60
         BGB1P0 =            0.0
         BGB2P0 =            0.0
         BGB3P0 =            0.0000D+00

         BGC0P0 =            0.0D+00
ccc         BGC0P0 =            0.76
         BGC0P0 =            0.0

         BGC1P0 =            0.0000D+00
         BGC2P0 =            0.0D+00
ccc         BGC2P0 =            0.15
         BGC2P0 =            0.0
         BGC3P0 =            0.0000D+00

         BGA1M0 =            0.0D+00
         BGA2M0 =            0.0D+00
ccc         BGA1M0 =            0.76
ccc         BGA2M0 =           -0.36
         BGA1M0 =            0.0
         BGA2M0 =           0.0
         BGA3M0 =            0.0000D+00

         BGB2M0 =            0.0
         BGB3M0 =            0.0000D+00

         BGC1M0 =            0.0D+00
         BGC2M0 =            0.0D+00
ccc         BGC1M0 =           -0.59
ccc         BGC2M0 =            0.47
         BGC1M0 =            0.
         BGC2M0 =            0.
         BGC3M0 =            0.0000D+00

         BGA0PP =            0.0000D+00
         BGA1PP =            0.0000D+00
         BGA2PP =            0.0000D+00
         BGA3PP =            0.0000D+00

         BGB1PP =            0.0000D+00
         BGB2PP =            0.0000D+00
         BGB3PP =            0.0000D+00

         BGC0PP =            0.0000D+00
         BGC1PP =            0.0000D+00
         BGC2PP =            0.0000D+00
         BGC3PP =            0.0000D+00

         BGA1MP =            0.0000D+00
         BGA2MP =            0.0000D+00
         BGA3MP =            0.0000D+00

         BGB2MP =            0.0000D+00
         BGB3MP =            0.0000D+00

         BGC1MP =            0.0000D+00
         BGC2MP =            0.0000D+00
         BGC3MP =            0.0000D+00
C
C       Now calculate the resonant partial wave helicity elements
C

         PVEC_GAM = SQRT((EPQ2 + (LS11_1535+XMPROT)**2)*
     *   (EPQ2+(LS11_1535-XMPROT)**2))/(2.*LS11_1535)
         PVEC_PI= SQRT((LS11_1535**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LS11_1535**2) - XMPI**2)
         FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LS11_1535*PIBS11_1535/WS11_1535)
         RAS11_1535 = FACT*A12S11_1535*19.73*SQRT(2./3.)
         RCS11_1535 = 0.

         PVEC_GAM = SQRT((EPQ2 + (LS11_1650+XMPROT)**2)*
     *   (EPQ2+(LS11_1650-XMPROT)**2))/(2.*LS11_1650)
         PVEC_PI= SQRT((LS11_1650**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LS11_1650**2) - XMPI**2)
         FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LS11_1650*PIBS11_1650/WS11_1650)
         RAS11_1650 = FACT*A12S11_1650*19.73*SQRT(2./3.)
         RCS11_1650 = 0.

         PVEC_GAM = SQRT((EPQ2 + (LP11_1440+XMPROT)**2)*
     *   (EPQ2+(LP11_1440-XMPROT)**2))/(2.*LP11_1440)
         PVEC_PI= SQRT((LP11_1440**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LP11_1440**2) - XMPI**2)
         FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LP11_1440*PIBP11_1440/WP11_1440)

         RAP11_1440 = -FACT*A12P11_1440*
     *                 DIP_EVF*19.73*SQRT(2./3.)
         RCP11_1440 =               0.0000D+00

         RAP11_1710 =               0.0000D+00
         RCP11_1710 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LP13_1720+XMPROT)**2)*
     *   (EPQ2+(LP13_1720-XMPROT)**2))/(2.*LP13_1720)
         PVEC_PI= SQRT((LP13_1720**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LP13_1720**2) - XMPI**2)
         FACT = SQRT(1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LP13_1720*PIBP13_1720/WP13_1720)
         RAP13_1720 =  FACT*A12P13_1720*19.73*SQRT(2./3.)
         RBP13_1720 = -FACT*SQRT(16./12.)*A32P13_1720*19.73*
     *                 SQRT(2./3.)
         RCP13_1720 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LD13_1520+XMPROT)**2)*
     *   (EPQ2+(LD13_1520-XMPROT)**2))/(2.*LD13_1520)
         PVEC_PI= SQRT((LD13_1520**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LD13_1520**2) - XMPI**2)
         FACT = SQRT(1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LD13_1520*PIBD13_1520/WD13_1520)
         RAD13_1520 = -FACT*A12D13_1520*19.73*SQRT(2./3.)
         RBD13_1520 =  FACT*SQRT(16./12.)*A32D13_1520*19.73*
     *                 SQRT(2./3.)
         RCD13_1520 =               0.0000D+00

         RAD13_1700 =               0.0000D+00
         RBD13_1700 =               0.0000D+00
         RCD13_1700 =               0.0000D+00

         RAD15_1675 =               0.0000D+00
         RBD15_1675 =               0.0000D+00
         RCD15_1675 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LF15_1680+XMPROT)**2)*
     *   (EPQ2+(LF15_1680-XMPROT)**2))/(2.*LF15_1680)
         PVEC_PI= SQRT((LF15_1680**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LF15_1680**2) - XMPI**2)
         FACT = SQRT(1./6./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LF15_1680*PIBF15_1680/WF15_1680)
         RAF15_1680 = -FACT*A12F15_1680*19.73*SQRT(2./3.)
         RBF15_1680 =  FACT*SQRT(16./32.)*A32F15_1680*19.73*
     *                 SQRT(2./3.)
         RCF15_1680 =              0.0000D+00

         RAG17_2190 =              0.0000D+00
         RBG17_2190 =              0.0000D+00
         RCG17_2190 =              0.0000D+00

         RAG19_2250 =              0.0000D+00
         RBG19_2250 =              0.0000D+00
         RCG19_2250 =              0.0000D+00

         RAH19_2220 =              0.0000D+00
         RBH19_2220 =              0.0000D+00
         RCH19_2220 =              0.0000D+00

         RAI111_2600 =             0.0000D+00
         RBI111_2600 =             0.0000D+00
         RCI111_2600 =             0.0000D+00
C
         PVEC_GAM = SQRT((EPQ2 + (LS31_1620+XMPROT)**2)*
     *   (EPQ2+(LS31_1620-XMPROT)**2))/(2.*LS31_1620)
         PVEC_PI= SQRT((LS31_1620**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LS31_1620**2) - XMPI**2)
         FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LS31_1620*PIBS31_1620/WS31_1620)
         RAS31_1620 =  FACT*A12S31_1620*19.73*SQRT(1./3.)
         RCS31_1620 =              0.0000D+00

         RAS31_1900 =              0.0000D+00
         RCS31_1900 =              0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LP31_1910+XMPROT)**2)*
     *   (EPQ2+(LP31_1910-XMPROT)**2))/(2.*LP31_1910)
         PVEC_PI= SQRT((LP31_1910**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LP31_1910**2) - XMPI**2)
         FACT = SQRT(1./2./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LP31_1910*PIBP31_1910/WP31_1910)*SQRT(1./3.)
         RAP31_1910 = -FACT*A12P31_1910*19.73*SQRT(1./3.)
         RCP31_1910 =              0.0000D+00


         PVEC_GAM = SQRT((EPQ2 + (LP33_1232+XMPROT)**2)*
     *   (EPQ2+(LP33_1232-XMPROT)**2))/(2.*LP33_1232)
         PVEC_GAM0 = SQRT(((LP33_1232+XMPROT)**2)*
     *   ((LP33_1232-XMPROT)**2))/(2.*LP33_1232)
         Q0_GAM=(LP33_1232**2 - XMPROT**2+EPQ2)/(2.*LP33_1232)
         PVEC_PI= SQRT((LP33_1232**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LP33_1232**2) - XMPI**2)
             FACT = SQRT(1./137.*LP33_1232/6./XMPROT/PVEC_PI/
     *              XMPROT**2/WP33_1232)*PVEC_GAM
         FFAC = 1./(1+EPQ2/0.71)**2*EXP(-0.21   *EPQ2)
         RAP33_1232 =  -0.5*3.*FFAC*FACT*SQRT(1./3.)*19.73
         RBP33_1232 =   3.*FFAC*FACT*SQRT(1./3.)*19.73
         RCP33_1232 =0.0D+00
cccc         RCP33_1232 =  -0.07*2.*SQRT(EPQ2)/Q0_GAM*3.*FFAC*
cccc     *                  FACT*SQRT(1./3.)*19.73

         PVEC_GAM = SQRT(EPQ2+ (LP33_1920+XMPROT)**2)*
     *   (EPQ2+(LP33_1920-XMPROT)**2)/(2.*LP33_1920)
         PVEC_PI= SQRT(((LP33_1920**2+XMPI**2-XMPROT**2)/
     *   (4*LP33_1920**2) - XMPI**2))
         FACT = SQRT((1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LP33_1920*PIBP33_1920/WP33_1920))
         RAP33_1920 =  FACT*A12P33_1920*SQRT(1./3.)
         RBP33_1920 = -FACT*SQRT(16./12.)*A32P33_1920*SQRT(1./3.)
         RCP33_1920 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LD33_1700+XMPROT)**2)*
     *   (EPQ2+(LD33_1700-XMPROT)**2))/(2.*LD33_1700)
         PVEC_PI= SQRT((LD33_1700**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LD33_1700**2) - XMPI**2)
         FACT = SQRT(1./4./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LD33_1700*PIBD33_1700/WD33_1700)
         RAD33_1700 = -FACT*A12D33_1700*19.73*SQRT(1./3.)
         RBD33_1700 =  FACT*SQRT(16./12.)*A32D33_1700*19.73*
     *                SQRT(1./3.)
         RCD33_1700 =            0.0000D+00

         RAD35_1930 =               0.0000D+00
         RBD35_1930 =               0.0000D+00
         RCD35_1930 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LF35_1905+XMPROT)**2)*
     *   (EPQ2+(LF35_1905-XMPROT)**2))/(2.*LF35_1905)
         PVEC_PI= SQRT((LF35_1905**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LF35_1905**2) - XMPI**2)
         FACT = SQRT(1./6./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LF35_1905*PIBF35_1905/WF35_1905)
         RAF35_1905 = -FACT*A12F35_1905*19.73*SQRT(1./3.)
         RBF35_1905 =  FACT*SQRT(16./32.)*A32F35_1905*19.73*
     *                 SQRT(1./3.)
         RCF35_1905 =               0.0000D+00

         PVEC_GAM = SQRT((EPQ2 + (LF37_1950+XMPROT)**2)*
     *   (EPQ2+(LF37_1950-XMPROT)**2))/(2.*LF37_1950)
         PVEC_PI= SQRT((LF37_1950**2+XMPI**2-XMPROT**2)**2/
     *   (4.*LF37_1950**2) - XMPI**2)
         FACT = SQRT(1./8./PI*PVEC_GAM/PVEC_PI*XMPROT/
     *   LF37_1950*PIBF37_1950/WF37_1950)
         RAF37_1950 =  FACT*A12F37_1950*19.73*SQRT(1./3.)
         RBF37_1950 = -FACT*SQRT(16./60.)*A32F37_1950*19.73*
     *                 SQRT(1./3.)
         RCF37_1950 =               0.0000D+00

         RAH311_2420 =              0.0000D+00
         RBH311_2420 =              0.0000D+00
         RCH311_2420 =              0.0000D+00


         sigma0 = sigmao(epw, epq2, epeps,epcos,epphi,epirea,
     *           0.,0.,0.,0.,sig0,sigu,sigt,sigl,sigi,fkt)
         xu_0=fkt*sigu
         xl_0=fkt*sigl*epeps
         xt_0=fkt*sigt*epeps*cos(2.*epphi*pi/180.)
         xi_0=fkt*sigi*sqrt(epeps*(1.+epeps)/2)*cos(epphi*pi/180.)
      return
      END



         COMPLEX FUNCTION EPRES(IREA,XMP,XMPI,W,W0,GAM,X,L,J,ISO,
     *   PIMP,GENER,TSTS11)
         common /debug/ibug

         COMPLEX QETA,GAMMAT,Q0ETA
         INTEGER L,J,ISO,IREA
         REAL XMP,XMPI,W,W0,GAM,X,PIMP0
         REAL PENER0,GENER0,XMETA,GAMMA
         LOGICAL TSTS11
         COMPLEX IM
         DATA IM/(0.,1.)/
         DATA XMETA/.5488/

         PENER0=(W0**2-XMP**2+XMPI**2)/(2.*W0)
         PIMP0 =  PENER0**2 - XMPI**2
         if (pimp0 .le. 0.)then
            write(6,*)' epres: pimp0 =',pimp0
            pimp0=abs(pimp0)
         endif
         PIMP0 = SQRT ( PIMP0 )
         GENER0 = (W0**2-XMP**2)/(2.*W0)

         GAMMA=GAM * ((PIMP/PIMP0)**(2.*L+1))*
     *   (((PIMP0**2+X**2)/(PIMP**2+X**2))**L)

         GAMMAG=GAM * ((GENER/GENER0)**(2.*J))*
     *   (((GENER0**2+X**2)/(GENER**2+X**2))**J)

         IF(TSTS11)THEN
           XETA=(W**2-XMP**2+XMETA**2)/(2.*W)
           XETA=XETA**2-XMETA**2
           IF(XETA.GE.0.)THEN
             QETA=CMPLX(SQRT(XETA),0.)
           ELSE
             QETA=CMPLX(0.,SQRT(-XETA))
           END IF
           Q0ETA=(W0**2-XMP**2+XMETA**2)/(2.*W0)
           Q0ETA=SQRT(Q0ETA**2-XMETA**2)
           GAMMAT=GAM*(PIMP/PIMP0)
         ELSE
           GAMMAT=CMPLX(GAMMA,0.)
         END IF
         EPRES=SQRT((GENER0*PIMP0)/(GENER*PIMP))
         EPRES=EPRES*W0*SQRT(GAMMA*GAMMAG)
         EPRES=EPRES/(W0**2-W**2-IM*W0*GAMMAT)
         IF(IREA.EQ.1)THEN
           IF(ISO.EQ.1)THEN
             EPRES=-EPRES/SQRT(2.)
           ELSE
             EPRES=EPRES*SQRT(2.)
           END IF
         END IF
         RETURN
         END


      SUBROUTINE BORNT(Q2,W,SINX,COSX,SINX2,COSX2,BORN,IREA)
      DIMENSION BORN(6),HBORN(6)
      DIMENSION
     *   XM(6,6),A(6),HI(6),IEPSI(-1:+1),XN(-1:1,6),
     *   AI(-1:1,6),IETA(6),F1(-1:1),F2(-1:1),FP(-1:1)
      DATA IEPSI/-1,1,1/,IETA/1,1,-1,1,-1,-1/
      DATA PI/3.1415927/,XMP/.9382796/,
     + XMPIP/0.139563/,XMPI0/0.1349630/,
     + MAG/0/
      DATA ER/.302862/,GR/13.5656/,SQ2/1.41421/

      IF(IREA.EQ.1.OR.IREA.EQ.2)XMPI=XMPI0
      IF(IREA.EQ.3.OR.IREA.EQ.4)XMPI=XMPIP
      if (cosx .ge. 1.0)cosx=.9999
      if (cosx .le. -1.0)cosx=-.9999
      S=W*W
      SIGMA=S-XMP**2
      Q0=(SIGMA-Q2)/(2.*W)
      Q =SQRT(Q0**2 + Q2)
      XK0PI=(SIGMA+XMPI**2)/(2.*W)
      XKPI =XK0PI**2-XMPI**2
      if (xkpi .lt. 0.)then
         DO  I=1,6
           BORN(I)=0.D0
         enddo
c         write(6,*)' bornt: xkpi',xkpi
         return
      endif
      XKPI =SQRT(XKPI)
      BETA=SIGMA+Q2/2.
      XNUE=-(W+XMP)
      OMEGA=W-XMP
      E1=W+XMP- Q0
      E2=W+XMP-XK0PI
      G1=-W+XMP+ Q0
      G2=-W+XMP+XK0PI
      T = 2.*XKPI*Q*COSX-2.*XK0PI*Q0+XMPI**2-Q2
      U = S - 2.*BETA-T+XMPI**2
      X1S = XMP*XMP - S
      X1U = XMP*XMP - U
      P2  = T - XMPI*XMPI
      FP(-1)=ER/(1.+Q2/0.5)
      FP(0)=0.D0
      FP(1)=0.D0
      GEP=1./(1.+3.04*Q2+1.54*Q2*Q2+.068*Q2*Q2*Q2)
      CAP=Q2/(4.*XMP*XMP)
        F1P=GEP*(1.+2.793*CAP)/(1.+CAP)
        F1N=GEP*(  -1.913*CAP)/(1.+CAP)
      IF(MAG.EQ.0)THEN
        F2P=0.D0
        F2N=0.D0
      ELSE
        F2P=GEP*1.793/(1.+CAP)
        F2N=GEP*(-1.913)/(1.+CAP)
        F2P=F2P/(2.*XMP)
        F2N=F2N/(2.*XMP)
      END IF
      F1(+1)=ER*(F1P-F1N)
      F1(-1)=F1(+1)
      F1(0 )=ER*(F1P+F1N)
      F2(+1)=ER*(F2P-F2N)
      F2(-1)=F2(+1)
      F2(0 )=ER*(F2P+F2N)
      P1=XK0PI-Q0
      XM(1,1)=-4.*OMEGA*E1
      XM(1,2)=BETA*(Q2+2.*XNUE*P1-P2)
      XM(1,3)=-2.*E1*OMEGA*OMEGA+(E1-Q0)*P2+Q2*(-XNUE+2.*P1)
      XM(1,4)=-2.*E1*OMEGA*OMEGA-(E1-Q0)*P2-Q2*(-XNUE+2.*P1)
      XM(1,5)=-Q2*(Q2+2.*XNUE*P1-P2)
      XM(1,6)=-4.*Q2*E1
      DO 16 IXM=1,6
16      XM(1,IXM)=XM(1,IXM)/(2.*Q*Q)
C      END DO
      HI(1)=G1*SQRT(E1*E2)/(8.*PI*W)
      P1=-(XK0PI-Q0)
      XM(2,1)=-4.*XNUE*G1
      XM(2,2)=BETA*(Q2+2.*OMEGA*P1-P2)
      XM(2,3)=-2.*G1*XNUE*XNUE+(G1+Q0)*P2+Q2*(-OMEGA+2.*P1)
      XM(2,4)=-2.*G1*XNUE*XNUE-(G1+Q0)*P2-Q2*(-OMEGA+2.*P1)
      XM(2,5)=-Q2*(Q2+2.*OMEGA*P1-P2)
      XM(2,6)=-4.*Q2*G1
      DO 26 IXM=1,6
26      XM(2,IXM)=XM(2,IXM)/(2.*Q*Q)
C     END DO
      HI(2)=-E1*SQRT(G1*G2)/(8.*PI*W)
      XM(3,1)=0.D0
      XM(3,2)=-BETA
      XM(3,3)=XNUE
      XM(3,4)=-XNUE
      XM(3,5)=Q2
      XM(3,6)=0.D0
      HI(3)=E2*SQRT(G1*G2)/(8.*PI*W)
      XM(4,1)=0.D0
      XM(4,2)=-BETA
      XM(4,3)=OMEGA
      XM(4,4)=-OMEGA
      XM(4,5)=Q2
      XM(4,6)=0.D0
      HI(4)=G2*SQRT(E1*E2)/(8.*PI*W)
      P1=2.*XK0PI-Q0
      XM(5,1)=2.*E1
      XM(5,2)=BETA*P1+(2.*W-1.5*Q0)*P2
      XM(5,3)=E1*OMEGA-XNUE*P1+P2
      XM(5,4)=E1*OMEGA+XNUE*P1-P2
      XM(5,5)=-Q2*P1+Q0*P2
      XM(5,6)=-2.*OMEGA*E1
      DO 56 IXM=1,6
56      XM(5,IXM)=XM(5,IXM)/(2.*Q*Q)
C     END DO
      HI(5)=-G1*SQRT(E1*E2)/(8*PI*W)
      P1=-2.*XK0PI+Q0
      XM(6,1)=2.*G1
      XM(6,2)=BETA*P1-(2.*W-1.5*Q0)*P2
      XM(6,3)=G1*XNUE-OMEGA*P1+P2
      XM(6,4)=G1*XNUE+OMEGA*P1-P2
      XM(6,5)=-Q2*P1-Q0*P2
      XM(6,6)=-2.*XNUE*G1
      DO 66 IXM=1,6
66      XM(6,IXM)=XM(6,IXM)/(2.*Q*Q)
C     END DO
      HI(6)=E1*SQRT(G1*G2)/(8.*PI*W)
      I1=0
      DO 71 I2=1,6
      XN(I1,I2)=GR*(1./X1S + IEPSI(I1)*IETA(I2)/X1U )
71    CONTINUE
      AI(I1,1)=-.5*F1(I1)*XN(I1,1)
      AI(I1,2)=F1(I1)*XN(I1,2)/P2
      AI(I1,3)=F2(I1)*XN(I1,3)/2.
      AI(I1,4)=F2(I1)*XN(I1,4)/2.
      IF(Q2.NE.0.D0)AI(I1,5)=F1(I1)*XN(I1,5)/(2.*P2)
     + + (1-IEPSI(I1))*(FP(I1)-F1(I1))*GR/(Q2*P2)
      IF(Q2.EQ.0.D0)AI(I1,5)=0.D0
      AI(I1,6)=0.D0
      IF(IREA.EQ.1.OR.IREA.EQ.2)I1=1
      IF(IREA.EQ.3.OR.IREA.EQ.4)I1=-1
      DO 72 I2=1,6
      XN(I1,I2)=GR*(1./X1S + IEPSI(I1)*IETA(I2)/X1U )
72    CONTINUE
      AI(I1,1)=-.5*F1(I1)*XN(I1,1)
      AI(I1,2)=F1(I1)*XN(I1,2)/P2
      AI(I1,3)=F2(I1)*XN(I1,3)/2.
      AI(I1,4)=F2(I1)*XN(I1,4)/2.
      IF(Q2.NE.0.D0)AI(I1,5)=F1(I1)*XN(I1,5)/(2.*P2)
     + + (1-IEPSI(I1))*(FP(I1)-F1(I1))*GR/(Q2*P2)
      IF(Q2.EQ.0.D0)AI(I1,5)=0.D0
      AI(I1,6)=0.D0
      IF(IREA.EQ.1)THEN
      DO 73 I2=1,6
73    A(I2)=AI(0,I2)+AI(1,I2)
      END IF
      IF(IREA.EQ.2)THEN
      DO 74 I2=1,6
74    A(I2)=-AI(0,I2)+AI(1,I2)
      END IF
      IF(IREA.EQ.3)THEN
      DO 75 I2=1,6
75    A(I2)=-SQ2*(AI(0,I2)+AI(-1,I2))
      END IF
      IF(IREA.EQ.4)THEN
      DO 76 I2=1,6
76    A(I2)= SQ2*(AI(0,I2)-AI(-1,I2))
      END IF

      DO 2 I=1,6
      HBORN(I)=0.D0
      DO 1 J=1,6
1     HBORN(I)=HBORN(I)+A(J)*XM(I,J)
2     HBORN(I)=HBORN(I)*HI(I)*19.732857D0
      BORN(1)=SINX*COSX2*(HBORN(3)+HBORN(4))/SQ2
      BORN(2)=-COSX2*(HBORN(1)+HBORN(2))/SQ2
      BORN(3)=-SINX*SINX2*(HBORN(3)-HBORN(4))/SQ2
      BORN(4)=SINX2   *  (HBORN(1)-HBORN(2))/SQ2
      BORN(5)=- SQRT(Q2) *SQ2*COSX2  *(HBORN(5)+HBORN(6))
      BORN(6)=  SQRT(Q2) *SQ2*SINX2  *(HBORN(5)-HBORN(6))
      RETURN
      END
      


         REAL FUNCTION SIGMAO(EPW,EPQ2,EPEPS,EPCOS,EPPHI,EPIREA,
     *   POL_ELEC,POL_TARG,POL_TARG_THETA,POL_TARG_PHI,sig0
     *   ,sigu,sigt,sigl,sigi,fkt)

         implicit none
         common /debug/ibug
         integer ibug

         integer sigtot,epirea
         real sig0,st,si
         real epw,epq2,epeps,epcos,epphi,pol_elec,pol_targ
         real pol_targ_theta,pol_targ_phi

      REAL*4
     * BGA0P0,BGA1P0,BGA2P0,BGA3P0,
     * BGA1M0,BGA2M0,BGA3M0,
     * BGB1P0,BGB2P0,BGB3P0,
     * BGB2M0,BGB3M0,
     * BGC0P0,BGC1P0,BGC2P0,BGC3P0,
     * BGC1M0,BGC2M0,BGC3M0,
     * BGA0PP,BGA1PP,BGA2PP,BGA3PP,
     * BGA1MP,BGA2MP,BGA3MP,
     * BGB1PP,BGB2PP,BGB3PP,
     * BGB2MP,BGB3MP,
     * BGC0PP,BGC1PP,BGC2PP,BGC3PP,
     * BGC1MP,BGC2MP,BGC3MP
      REAL
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710,
     * LP13_1720, LD13_1520, LD13_1700,
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600,
     * LS31_1620, LS31_1900, LP31_1910,
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420
         COMMON /ANREG2/
     * BGA0P0,BGA1P0,BGA2P0,BGA3P0,
     * BGA1M0,BGA2M0,BGA3M0,
     * BGB1P0,BGB2P0,BGB3P0,
     * BGB2M0,BGB3M0,
     * BGC0P0,BGC1P0,BGC2P0,BGC3P0,
     * BGC1M0,BGC2M0,BGC3M0,
     * BGA0PP,BGA1PP,BGA2PP,BGA3PP,
     * BGA1MP,BGA2MP,BGA3MP,
     * BGB1PP,BGB2PP,BGB3PP,
     * BGB2MP,BGB3MP,
     * BGC0PP,BGC1PP,BGC2PP,BGC3PP,
     * BGC1MP,BGC2MP,BGC3MP,
     * RAS11_1535,             RCS11_1535,
     * RAS11_1650,             RCS11_1650,
     * RAP11_1440,             RCP11_1440,
     * RAP11_1710,             RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,
     * RAP13_1910, RBP13_1910, RCP13_1910,
     * RAD13_1520, RBD13_1520, RCD13_1520,
     * RAD13_1700, RBD13_1700, RCD13_1700,
     * RAD15_1675, RBD15_1675, RCD15_1675,
     * RAF15_1680, RBF15_1680, RCF15_1680,
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,             RCS31_1620,
     * RAS31_1900,             RCS31_1900,
     * RAP31_1910,             RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232,
     * RAP33_1920, RBP33_1920, RCP33_1920,
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420,
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710,
     * WP13_1720, WD13_1520, WD13_1700,
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600,
     * WS31_1620, WS31_1900, WP31_1910,
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420,
     * LS11_1535, LS11_1650, LP11_1440, LP11_1710,
     * LP13_1720, LD13_1520, LD13_1700,
     * LD15_1675, LF15_1680,
     * LG17_2190,
     * LG19_2250, LH19_2220,
     * LI111_2600,
     * LS31_1620, LS31_1900, LP31_1910,
     * LP33_1232, LP33_1920, LD33_1700,
     * LD35_1930, LF35_1905,
     * LF37_1950,
     * LH311_2420


      COMPLEX
     * QAS11_1535, QCS11_1535,
     * QAS11_1650, QCS11_1650,
     * QAP11_1440, QCP11_1440,
     * QAP11_1710, QCP11_1710,
     * QAP13_1720, QBP13_1720, QCP13_1720,
     * QAD13_1520, QBD13_1520, QCD13_1520,
     * QAD13_1700, QBD13_1700, QCD13_1700,
     * QAD15_1675, QBD15_1675, QCD15_1675,
     * QAF15_1680, QBF15_1680, QCF15_1680,
     * QAG17_2190, QBG17_2190, QCG17_2190,
     * QAG19_2250, QBG19_2250, QCG19_2250,
     * QAH19_2220, QBH19_2220, QCH19_2220,
     * QAI111_2600,QBI111_2600,QCI111_2600,
     * QAS31_1620, QCS31_1620,
     * QAS31_1900, QCS31_1900,
     * QAP31_1910, QCP31_1910,
     * QAP33_1232, QBP33_1232, QCP33_1232,
     * QAP33_1920, QBP33_1920, QCP33_1920,
     * QAD33_1700, QBD33_1700, QCD33_1700,
     * QAD35_1930, QBD35_1930, QCD35_1930,
     * QAF35_1905, QBF35_1905, QCF35_1905,
     * QAF37_1950, QBF37_1950, QCF37_1950,
     * QAH311_2420,QBF311_2420,QCF311_2420,
     * QBH311_2420,QCH311_2420,
     * XA0P  ,XA1P  ,XA2P  ,XA3P  , XA4P, XA5P, XA6P,
     * XA1M  ,XA2M  ,XA3M  ,XA4M,  XA5M,   XA6M,
     * XB1P  ,XB2P  ,XB3P  , XB4P,  XB5P,
     * XB2M  ,XB3M  , XB4M,  XB5M,  XB6M,
     * XC0P  ,XC1P  ,XC2P  ,XC3P  , XC4P, XC5P,
     * XC1M  ,XC2M  ,XC3M, XC4M, XC5M,  XC6M
      COMPLEX
     * RS11_1535, RS11_1650, RP11_1440, RP11_1710,
     * RP13_1720, RD13_1520, RD13_1700,
     * RD15_1675, RF15_1680,
     * RG17_2190,
     * RG19_2250, RH19_2220,
     * RI111_2600,
     * RS31_1620, RS31_1900, RP31_1910,
     * RP33_1232, RP33_1920, RD33_1700,
     * RD35_1930, RF35_1905,
     * RF37_1950,
     * RH311_2420
      COMPLEX EPH1,EPH2,EPH3,EPH4,EPH5,EPH6
      COMPLEX ADDA0,ADDA1,ADDA2,ADDA3,ADDA4,ADDA5
      COMPLEX ADDB1,ADDB2,ADDB3,ADDB4,ADDB5
      COMPLEX ADDC0,ADDC1,ADDC2,ADDC3,ADDC4,ADDC5
      COMPLEX SUBA0,SUBA1,SUBA2,SUBA3,SUBA4,SUBA5
      COMPLEX SUBB1,SUBB2,SUBB3,SUBB4,SUBB5
      COMPLEX SUBC0,SUBC1,SUBC2,SUBC3,SUBC4,SUBC5
      COMPLEX HNP,HNM,HFP,HFM,HN0,HF0
      COMPLEX X1,X2,Y1,Y2,Y3,Y4,Z1,Z2
      COMPLEX R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11
      REAL HB(6)

      real
     * RAS11_1535,             RCS11_1535,
     * RAS11_1650,             RCS11_1650,
     * RAP11_1440,             RCP11_1440,
     * RAP11_1710,             RCP11_1710,
     * RAP13_1720, RBP13_1720, RCP13_1720,
     * RAP13_1910, RBP13_1910, RCP13_1910,
     * RAD13_1520, RBD13_1520, RCD13_1520,
     * RAD13_1700, RBD13_1700, RCD13_1700,
     * RAD15_1675, RBD15_1675, RCD15_1675,
     * RAF15_1680, RBF15_1680, RCF15_1680,
     * RAG17_2190, RBG17_2190, RCG17_2190,
     * RAG19_2250, RBG19_2250, RCG19_2250,
     * RAH19_2220, RBH19_2220, RCH19_2220,
     * RAI111_2600,RBI111_2600,RCI111_2600,
     * RAS31_1620,             RCS31_1620,
     * RAS31_1900,             RCS31_1900,
     * RAP31_1910,             RCP31_1910,
     * RAP33_1232, RBP33_1232, RCP33_1232,
     * RAP33_1920, RBP33_1920, RCP33_1920,
     * RAD33_1700, RBD33_1700, RCD33_1700,
     * RAD35_1930, RBD35_1930, RCD35_1930,
     * RAF35_1905, RBF35_1905, RCF35_1905,
     * RAF37_1950, RBF37_1950, RCF37_1950,
     * RAH311_2420,RBH311_2420,RCH311_2420

      REAL
     * WS11_1535, WS11_1650, WP11_1440, WP11_1710,
     * WP13_1720, WD13_1520, WD13_1700,
     * WD15_1675, WF15_1680,
     * WG17_2190,
     * WG19_2250, WH19_2220,
     * WI111_2600,
     * WS31_1620, WS31_1900, WP31_1910,
     * WP33_1232, WP33_1920, WD33_1700,
     * WD35_1930, WF35_1905,
     * WF37_1950,
     * WH311_2420

       COMPLEX ZERO
       INTEGER IR,IREA
       COMPLEX EPRES
       REAL YMPI,W,PIM,GK

       real sq2,pi,xmprot,xmpip,xmpi0,xmeta,xmpi
       real dum,epsin,epsin2,epcos2
       real pener,pimp,gener,gimp,gamk,x,p1s,p2s,p3s,p4s,p5s,p6s,p7s
       real p2s2,p3s2,p4s2,p5s2,p6s2,p7s2,x0,F,fkt
       real bpipy0,bpipy1,bpipy2,bpipy3,q2,sinx,cosx,sinx2
       real cosx2,eph1b,eph2b,eph3b,eph4b,eph5b,eph6b
       real sind,cosd
       real sigu,sigt,sigl,sigi,sige,sigpt,siget,xila
       real theta_gam,enue,theta_elec,e_elec_scatt
       real e_elec_prim,pol_x,pol_y,pol_z

      DATA ZERO/(0.,0.)/
      DATA SQ2        /1.41421  /
      DATA PI         /3.1415926  /
      DATA XMPROT     /.938     /
      DATA XMPIP      /.1395    /
      DATA XMPI0      /.1349    /
      DATA XMETA      /.5488    /

        sig0=0.
         DUM=ACOS(EPCOS)
         EPSIN=SIN(DUM)
         EPSIN2=SIN(DUM/2.)
         EPCOS2=COS(DUM/2.)
         IF(EPIREA.EQ.1)XMPI=XMPI0
         IF(EPIREA.EQ.3)XMPI=XMPIP
           PENER  =(EPW**2-XMPROT**2+XMPI   **2)/(2.*EPW   )
           PIMP   =SQRT(PENER   **2 - XMPI   **2)
           GENER  =(EPW   **2-XMPROT**2-EPQ2   )/(2.*EPW   )
           GIMP   =SQRT(GENER**2 + EPQ2   )
           GAMK   =(EPW**2-XMPROT**2)/(2.*EPW   )
           X=EPCOS
           P1S  =1.
           P2S  =3.*X
           P3S  =(15.*X*X-3.)/2.
           P4S  =(35.*X*X*X-15.*X)/2.
           P5S  =(315.*X*X*X*X-210.*X*X+15.)/8.
           P6S  =(693.*X*X*X*X*X-630.*X*X*X)/8.
           P7S  =(3003.*X*X*X*X*X*X-3465.*X*X*X*X)/16.
           P2S2 =3.
           P3S2 =15.*X
           P4S2 =(105.*X*X-15.)/2.
           P5S2 =(315.*X*X*X-105.*X)/2.
           P6S2 =(3465.*X*X*X*X-1890.*X*X)/8.
           P7S2 =(18018.*X*X*X*X*X-13860.*X*X*X)/16.
           X0=1./SQRT((PIMP   **2)+(.35**2))
           F=(1.+(PIMP   **2)/.71)
           BPIPY0   =1./F
           BPIPY1   =X0*PIMP/F
           BPIPY2   =X0*X0*PIMP*PIMP/F
           BPIPY3   =X0*X0*X0*PIMP*PIMP*PIMP/F
           FKT=PIMP/GAMK
           Q2=EPQ2
           W =EPW
           SINX =EPSIN
           COSX =EPCOS
           SINX2=EPSIN2
           COSX2=EPCOS2
           IREA=EPIREA
           CALL BORNT(Q2,W,SINX,COSX,SINX2,COSX2,HB,IREA)
           EPH1B   =HB(1)
           EPH2B   =HB(2)
           EPH3B   =HB(3)
           EPH4B   =HB(4)
           EPH5B   =HB(5)
           EPH6B   =HB(6)
           IR=EPIREA
           YMPI  =  XMPI
           W     =  EPW
           PIM   =  PIMP
           GK    =  GAMK


       RS11_1535 = EPRES(IR,XMPROT,YMPI,W,LS11_1535,WS11_1535,
     * .35,0,1,1,PIM,GK,.FALSE.)
       RS11_1650  = EPRES(IR,XMPROT,YMPI,W,LS11_1650,WS11_1650,
     * .35,0,1,1,PIM,GK,.FALSE.)
       RP11_1440 = EPRES(IR,XMPROT,YMPI,W,LP11_1440,WP11_1440,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP11_1710 = EPRES(IR,XMPROT,YMPI,W,LP11_1710,WP11_1710,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RP13_1720 = EPRES(IR,XMPROT,YMPI,W,LP13_1720,WP13_1720,
     * .35,1,1,1,PIM,GK,.FALSE.)
       RD13_1520 = EPRES(IR,XMPROT,YMPI,W,LD13_1520,WD13_1520,
     * .35,2,1,1,PIM,GK,.FALSE.)
       RD13_1700 = EPRES(IR,XMPROT,YMPI,W,LD13_1700,WD13_1700,
     * .35,2,1,1,PIM,GK,.FALSE.)
       RD15_1675 = EPRES(IR,XMPROT,YMPI,W,LD15_1675,WD15_1675,
     * .35,2,2,1,PIM,GK,.FALSE.)
       RF15_1680 = EPRES(IR,XMPROT,YMPI,W,LF15_1680,WF15_1680,
     * .35,3,2,1,PIM,GK,.FALSE.)
       RG17_2190 = EPRES(IR,XMPROT,YMPI,W,LG17_2190,WG17_2190,
     * .35,4,3,1,PIM,GK,.FALSE.)
       RG19_2250 = EPRES(IR,XMPROT,YMPI,W,LG19_2250,WG19_2250,
     * .35,4,4,1,PIM,GK,.FALSE.)
       RH19_2220 = EPRES(IR,XMPROT,YMPI,W,LH19_2220,WH19_2220,
     * .35,5,4,1,PIM,GK,.FALSE.)
       RI111_2600 = EPRES(IR,XMPROT,YMPI,W,LI111_2600,WI111_2600,
     * .35,6,5,1,PIM,GK,.FALSE.)


       RS31_1620 = EPRES(IR,XMPROT,YMPI,W,LS31_1620,WS31_1620,
     * .35,0,1,3,PIM,GK,.FALSE.)
       RS31_1900 = EPRES(IR,XMPROT,YMPI,W,LS31_1900,WS31_1900,
     * .35,0,1,3,PIM,GK,.FALSE.)
       RP31_1910 = EPRES(IR,XMPROT,YMPI,W,LP31_1910,WP31_1910,
     * .35,1,1,3,PIM,GK,.FALSE.)
       RP33_1232 = EPRES(IR,XMPROT,YMPI,W,LP33_1232,WP33_1232,
     * .185,1,1,3,PIM,GK,.FALSE.)
       RP33_1920 = EPRES(IR,XMPROT,YMPI,W,LP33_1920,WP33_1920,
     * .35,1,1,3,PIM,GK,.FALSE.)
       RD33_1700 = EPRES(IR,XMPROT,YMPI,W,LD33_1700,WD33_1700,
     * .35,2,1,3,PIM,GK,.FALSE.)
       RD35_1930 = EPRES(IR,XMPROT,YMPI,W,LD35_1930,WD35_1930,
     * .35,2,2,3,PIM,GK,.FALSE.)
       RF35_1905 = EPRES(IR,XMPROT,YMPI,W,LF35_1905,WF35_1905,
     * .35,3,2,3,PIM,GK,.FALSE.)
       RF37_1950 = EPRES(IR,XMPROT,YMPI,W,LF37_1950,WF37_1950,
     * .35,3,3,3,PIM,GK,.FALSE.)
       RH311_2420 = EPRES(IR,XMPROT,YMPI,W,LH311_2420,WH311_2420,
     * .35,5,5,3,PIM,GK,.FALSE.)


           IF(EPIREA   .EQ.1)THEN
             XA0P    = BGA0P0*BPIPY0   +ZERO
             XA1P    = BGA1P0*BPIPY1   +ZERO
             XA2P    = BGA2P0*BPIPY2   +ZERO
             XA3P    = BGA3P0*BPIPY3   +ZERO
             XA4P    =                  ZERO
             XA5P    =                  ZERO

             XB1P    = BGB1P0*BPIPY1   +ZERO
             XB2P    = BGB2P0*BPIPY2   +ZERO
             XB3P    = BGB3P0*BPIPY3   +ZERO
             XB4P    =                  ZERO
             XB5P    =                  ZERO

             XC0P    = BGC0P0*BPIPY0   +ZERO
             XC1P    = BGC1P0*BPIPY1   +ZERO
             XC2P    = BGC2P0*BPIPY2   +ZERO
             XC3P    = BGC3P0*BPIPY3   +ZERO
             XC4P    =                  ZERO
             XC5P    =                  ZERO

             XA1M    = BGA1M0*BPIPY1   +ZERO
             XA2M    = BGA2M0*BPIPY2   +ZERO
             XA3M    = BGA3M0*BPIPY3   +ZERO
             XA4M    =                  ZERO
             XA5M    =                  ZERO
             XA6M    =                  ZERO

             XB2M    = BGB2M0*BPIPY2   +ZERO
             XB3M    = BGB3M0*BPIPY3   +ZERO
             XB4M    =                  ZERO
             XB5M    =                  ZERO
             XB6M    =                  ZERO

             XC1M    = BGC1M0*BPIPY1   +ZERO
             XC2M    = BGC2M0*BPIPY2   +ZERO
             XC3M    = BGC3M0*BPIPY3   +ZERO
             XC4M    =                  ZERO
             XC5M    =                  ZERO
             XC6M    =                  ZERO

           ELSE
             XA0P    = BGA0PP*BPIPY0   +ZERO
             XA1P    = BGA1PP*BPIPY1   +ZERO
             XA2P    = BGA2PP*BPIPY2   +ZERO
             XA3P    = BGA3PP*BPIPY3   +ZERO
             XA4P    =                  ZERO
             XA5P    =                  ZERO
             XA6P    =                  ZERO

             XB1P    = BGB1PP*BPIPY1   +ZERO
             XB2P    = BGB2PP*BPIPY2   +ZERO
             XB3P    = BGB3PP*BPIPY3   +ZERO
             XB4P    =                  ZERO
             XB5P    =                  ZERO

             XC0P    = BGC0PP*BPIPY0   +ZERO
             XC1P    = BGC1PP*BPIPY1   +ZERO
             XC2P    = BGC2PP*BPIPY2   +ZERO
             XC3P    = BGC3PP*BPIPY3   +ZERO
             XC4P    =                  ZERO
             XC5P    =                  ZERO

             XA1M    = BGA1MP*BPIPY1   +ZERO
             XA2M    = BGA2MP*BPIPY2   +ZERO
             XA3M    = BGA3MP*BPIPY3   +ZERO
             XA4M    =                  ZERO
             XA5M    =                  ZERO
             XA6M    =                  ZERO

             XB2M    = BGB2MP*BPIPY2   +ZERO
             XB3M    = BGB3MP*BPIPY3   +ZERO
             XB4M    =                  ZERO
             XB5M    =                  ZERO
             XB6M    =                  ZERO

             XC1M    = BGC1MP*BPIPY1   +ZERO
             XC2M    = BGC2MP*BPIPY2   +ZERO
             XC3M    = BGC3MP*BPIPY3   +ZERO
             XC4M    =                  ZERO
             XC5M    =                  ZERO
             XC6M    =                  ZERO

           END IF


         QAS11_1535   =RS11_1535   * RAS11_1535
         QCS11_1535   =RS11_1535   * RCS11_1535

         QAS11_1650   =RS11_1650   * RAS11_1650
         QCS11_1650   =RS11_1650   * RCS11_1650

         QAP11_1440   =RP11_1440   * RAP11_1440
         QCP11_1440   =RP11_1440   * RCP11_1440

         QAP11_1710   =RP11_1710   * RAP11_1710
         QCP11_1710   =RP11_1710   * RCP11_1710

         QAP13_1720   =RP13_1720   * RAP13_1720
         QBP13_1720   =RP13_1720   * RBP13_1720
         QCP13_1720   =RP13_1720   * RCP13_1720

         QAD13_1520   =RD13_1520   * RAD13_1520
         QBD13_1520   =RD13_1520   * RBD13_1520
         QCD13_1520   =RD13_1520   * RCD13_1520

         QAD13_1700   =RD13_1700   * RAD13_1700
         QBD13_1700   =RD13_1700   * RBD13_1700
         QCD13_1700   =RD13_1700   * RCD13_1700

         QAD15_1675   =RD15_1675   * RAD15_1675
         QBD15_1675   =RD15_1675   * RBD15_1675
         QCD15_1675   =RD15_1675   * RCD15_1675

         QAF15_1680   =RF15_1680   * RAF15_1680
         QBF15_1680   =RF15_1680   * RBF15_1680
         QCF15_1680   =RF15_1680   * RCF15_1680

         QAG17_2190   =RG17_2190   * RAG17_2190
         QBG17_2190   =RG17_2190   * RBG17_2190
         QCG17_2190   =RG17_2190   * RCG17_2190

         QAG19_2250   =RG19_2250   * RAG19_2250
         QBG19_2250   =RG19_2250   * RBG19_2250
         QCG19_2250   =RG19_2250   * RCG19_2250

         QAH19_2220   =RH19_2220   * RAH19_2220
         QBH19_2220   =RH19_2220   * RBH19_2220
         QCH19_2220   =RH19_2220   * RCH19_2220

         QAI111_2600  =RI111_2600  * RAI111_2600
         QBI111_2600  =RI111_2600  * RBI111_2600
         QCI111_2600  =RI111_2600  * RCI111_2600

         QAS31_1620   =RS31_1620   * RAS31_1620
         QCS31_1620   =RS31_1620   * RCS31_1620

         QAS31_1900   =RS31_1900   * RAS31_1900
         QCS31_1900   =RS31_1900   * RCS31_1900

         QAP31_1910   =RP31_1910   * RAP31_1910
         QCP31_1910   =RP31_1910   * RCP31_1910

         QAP33_1232   =RP33_1232   * RAP33_1232
         QBP33_1232   =RP33_1232   * RBP33_1232
         QCP33_1232   =RP33_1232   * RCP33_1232

         QAP33_1920   =RP33_1920   * RAP33_1920
         QBP33_1920   =RP33_1920   * RBP33_1920
         QCP33_1920   =RP33_1920   * RCP33_1920


         QAD33_1700   =RD33_1700   * RAD33_1700
         QBD33_1700   =RD33_1700   * RBD33_1700
         QCD33_1700   =RD33_1700   * RCD33_1700

         QAD35_1930   =RD35_1930   * RAD35_1930
         QBD35_1930   =RD35_1930   * RBD35_1930
         QCD35_1930   =RD35_1930   * RCD35_1930

         QAF35_1905   =RF35_1905   * RAF35_1905
         QBF35_1905   =RF35_1905   * RBF35_1905
         QCF35_1905   =RF35_1905   * RCF35_1905

         QAF37_1950   =RF37_1950   * RAF37_1950
         QBF37_1950   =RF37_1950   * RBF37_1950
         QCF37_1950   =RF37_1950   * RCF37_1950

         QAH311_2420  =RH311_2420  * RAH311_2420
         QBH311_2420  =RH311_2420  * RBH311_2420
         QCH311_2420  =RH311_2420  * RCH311_2420

         XA0P=XA0P + QAS11_1535 + QAS11_1650 + QAS31_1620 + QAS31_1900
         XC0P=XC0P + QCS11_1535 + QCS11_1650 + QCS31_1620 + QCS31_1900

         XA1M=XA1M + QAP11_1440 + QAP11_1710 + QAP31_1910
         XC1M=XC1M + QCP11_1440 + QCP11_1710 + QCP31_1910

         XA1P=XA1P + QAP13_1720 + QAP33_1232 + QAP33_1920
         XB1P=XB1P + QBP13_1720 + QBP33_1232 + QBP33_1920
         XC1P=XC1P + QCP13_1720 + QCP33_1232 + QCP33_1920

         XA2M = XA2M + QAD13_1520 + QAD13_1700 + QAD33_1700
         XB2M = XB2M + QBD13_1520 + QBD13_1700 + QBD33_1700
         XC2M = XC2M + QCD13_1520 + QCD13_1700 + QCD33_1700

         XA2P = XA2P + QAD15_1675 + QAD35_1930
         XB2P = XB2P + QBD15_1675 + QBD35_1930
         XC2P = XC2P + QCD15_1675 + QCD35_1930

         XA3M = XA3M + QAF15_1680 + QAF35_1905
         XB3M = XB3M + QBF15_1680 + QBF35_1905
         XC3M = XC3M + QCF15_1680 + QCF35_1905

         XA3P = XA3P + QAF37_1950
         XB3P = XB3P + QBF37_1950
         XC3P = XC3P + QCF37_1950

         XA4M = XA4M + QAG17_2190
         XB4M = XB4M + QBG17_2190
         XC4M = XC4M + QCG17_2190

         XA4P = XA4P + QAG19_2250
         XB4P = XB4P + QBG19_2250
         XC4P = XC4P + QCG19_2250

         XA5M = XA5M + QAH19_2220
         XB5M = XB5M + QBH19_2220
         XC5M = XC5M + QCH19_2220

         XA5P = XA5P + QAH311_2420
         XB5P = XB5P + QBH311_2420
         XC5P = XC5P + QCH311_2420

         XA6M = XA6M + QAI111_2600
         XB6M = XB6M + QBI111_2600
         XC6M = XC6M + QCI111_2600
c         write(6,*)' xa0p-xa5p',xa0p,xa1p,xa2p,xa3p,xa4p,xa5p
c         write(6,*)' xa1m-xa6m',xa1m,xa2m,xa3m,xa4m,xa5m,xa6m
c         write(6,*)' xb1p-xb5p',xb1p,xb2p,xb3p,xb4p,xb5p
c         write(6,*)' xb2m-xb6m',xb2m,xb3m,xb4m,xb5m,xb6m
c         write(6,*)' xc0p-xc5p',xc0p,xc1p,xc2p,xc3p,xc4p,xc5p
c         write(6,*)' xc1m-xc6m',xc1m,xc2m,xc3m,xc4m,xc5m,xc6m


         SUBA0=XA0P-XA1M
         SUBA1=XA1P-XA2M
         SUBA2=XA2P-XA3M
         SUBA3=XA3P-XA4M
         SUBA4=XA4P-XA5M
         SUBA5=XA5P-XA6M

         SUBC0=XC0P-XC1M
         SUBC1=XC1P-XC2M
         SUBC2=XC2P-XC3M
         SUBC3=XC3P-XC4M
         SUBC4=XC4P-XC5M
         SUBC5=XC5P-XC6M

         SUBB1=XB1P-XB2M
         SUBB2=XB2P-XB3M
         SUBB3=XB3P-XB4M
         SUBB4=XB4P-XB5M
         SUBB5=XB5P-XB6M

         ADDA0=XA0P+XA1M
         ADDA1=XA1P+XA2M
         ADDA2=XA2P+XA3M
         ADDA3=XA3P+XA4M
         ADDA4=XA4P+XA5M
         ADDA5=XA5P+XA6M

         ADDC0=XC0P+XC1M
         ADDC1=XC1P+XC2M
         ADDC2=XC2P+XC3M
         ADDC3=XC3P+XC4M
         ADDC4=XC4P+XC5M
         ADDC5=XC5P+XC6M

         ADDB1=XB1P+XB2M
         ADDB2=XB2P+XB3M
         ADDB3=XB3P+XB4M
         ADDB4=XB4P+XB5M
         ADDB5=XB5P+XB6M
c
C   Now do the partial wave expansion for the Walker amplitudes
C
       EPH1   =SUBB1*(       -P2S2   )
     *       + SUBB2*(P2S2   -P3S2   )
     *       + SUBB3*(P3S2   -P4S2   )
     *       + SUBB4*(P4S2   -P5S2   )
     *       + SUBB5*(P5S2   -P6S2   )

       EPH1   = EPH1   * EPSIN   * EPCOS2   /SQ2

       EPH2  = SUBA0   *(       -P1S   )
     *       + SUBA1   *(P1S    -P2S   )
     *       + SUBA2   *(P2S    -P3S   )
     *       + SUBA3   *(P3S    -P4S   )
     *       + SUBA4   *(P4S    -P5S   )
     *       + SUBA5   *(P5S    -P6S   )
       EPH2  = EPH2    * EPCOS2   * SQ2


       EPH3   = ADDB1    *(       +P2S2   )
     *       +  ADDB2    *(P2S2   +P3S2   )
     *       +  ADDB3    *(P3S2   +P4S2   )
     *       +  ADDB4    *(P4S2   +P5S2   )
     *       +  ADDB5    *(P5S2   +P6S2   )
       EPH3   = EPH3   *  EPSIN   * EPSIN2   /SQ2


       EPH4  = ADDA0     * (       +P1S   )
     *       + ADDA1     * (P1S    +P2S   )
     *       + ADDA2     * (P2S    +P3S   )
     *       + ADDA3     * (P3S    +P4S   )
     *       + ADDA4     * (P4S    +P5S   )
     *       + ADDA5     * (P5S    +P6S   )
       EPH4   = EPH4   * EPSIN2   * SQ2


       EPH5   =SUBC0     *(       -P1S   )
     *       + SUBC1     *(P1S    -P2S   )
     *       + SUBC2     *(P2S    -P3S   )
     *       + SUBC3     *(P3S    -P4S   )
     *       + SUBC4     *(P4S    -P5S   )
     *       + SUBC5     *(P5S    -P6S   )
       EPH5   = EPH5   * EPCOS2   * SQ2


       EPH6  = ADDC0     * (       +P1S   )
     *       + ADDC1     * (P1S    +P2S   )
     *       + ADDC2     * (P2S    +P3S   )
     *       + ADDC3     * (P3S    +P4S   )
     *       + ADDC4     * (P4S    +P5S   )
     *       + ADDC5     * (P5S    +P6S   )
        EPH6   =  EPH6   * EPSIN2   * SQ2

C
C      Now add the Born term contributions. Remember, these are real functions!
C
         EPH1   =EPH1   +CMPLX(EPH1B   ,0.)
         EPH2   =EPH2   +CMPLX(EPH2B   ,0.)
         EPH3   =EPH3   +CMPLX(EPH3B   ,0.)
         EPH4   =EPH4   +CMPLX(EPH4B   ,0.)
         EPH5   =EPH5   +CMPLX(EPH5B   ,0.)
         EPH6   =EPH6   +CMPLX(EPH6B   ,0.)
C
C       Now we are ready to calculate all the observables
C
C          First the unpolarized transverse term
C

           SIGU   =(CABS(EPH1   )**2
     *            +CABS(EPH2   )**2
     *            +CABS(EPH3   )**2
     *            +CABS(EPH4   )**2)/2
C
C            Next the transverse polarized term
C
           SIGT   =REAL(EPH2   *CONJG(EPH3   )
     *            -     EPH1   *CONJG(EPH4   ))
C
C             And the longitudinal term
C
           SIGL   =CABS(EPH5   )**2
     *            +CABS(EPH6   )**2
C
C             Finally, the longitudinal-transverse interference term
C
           SIGI   =SQ2*REAL(EPH5   *CONJG(EPH4   -EPH1   )
     *            +         EPH6   *CONJG(EPH3   +EPH2   ))
C
C              Sum up all the unpolarized cross section terms
C

           st = EPEPS*SIGT
           si = SQRT(EPEPS*(1.+EPEPS)/2.)*SIGI
             sig0  = SIGU       + EPEPS*SIGL
           SIGMAO   = sig0
     *             + st*COS(EPPHI*PI/90.)
     *             + si*COS(EPPHI*PI/180.)
C
         if(POL_ELEC .eq. 0. .and. POL_TARG .eq. 0.)then
            sigmao=sigmao*FKT

            return
         endif

C      Convert Walker amplitudes to Bartl & Majerotto spin-flip and non
C      spin-flip amplitudes
C
           HNP=(EPH4+EPH1)/SQ2
           HNM=(EPH4-EPH1)/SQ2
           HFP=(EPH3-EPH2)/SQ2
           HFM=(EPH3+EPH2)/SQ2
           HN0=EPH5
           HF0=EPH6
C
C     Define the appreviations for the various polarized cross section terms
C     Polarized beam, polarized target, polarized beam - polarized target
C
           X1=HF0*CONJG(HNP) + HN0*CONJG(HFP)
           X2=HFM*CONJG(HNP) + HNM*CONJG(HFP)
           Y1=HNP*CONJG(HFP) + HNM*CONJG(HFM)
           Y2=HNM*CONJG(HFM) - HNP*CONJG(HFP)
           Y3=HN0*CONJG(HF0)
           Y4=HN0*CONJG(HFM) - HF0*CONJG(HNM)
           Z1=HN0*CONJG(HNP) - HF0*CONJG(HFP)
           Z2=HNM*CONJG(HNP) - HFM*CONJG(HFP)
C
C               These are the terms for the recoil polarizations
C
           R1=CONJG(HN0)*HF0
           R2=CONJG(HNP)*HFP - CONJG(HNM)*HFM
           R3=CONJG(HN0)*HFM + CONJG(HF0)*HNP
           R4=CONJG(HF0)*HNM - CONJG(HN0)*HFM
           R5=CONJG(HN0)*HNP + CONJG(HF0)*HFP
           R6=CONJG(HNP)*HFP + CONJG(HNM)*HFM
           R7=CONJG(HFM)*HNP + CONJG(HFP)*HNM
           R8=CONJG(HFP)*HFM + CONJG(HNP)*HNM
           R9=CONJG(HNP)*HFM - CONJG(HNM)*HFP
           R10=CONJG(HNP)*HNM + CONJG(HFP)*HFM
           R11=CONJG(HN0)*HNM + CONJG(HF0)*HFM


           XILA=SQRT(2.*EPEPS*(EPEPS+1.))
           SIGE=-POL_ELEC*SQRT(2.*EPEPS*(1.-EPEPS))*SIN(EPPHI*PI/180.)*
     *     AIMAG(HN0*CONJG(HNM)+HF0*CONJG(HFM))

c
c           SIGPT=PX*(-XILA*SIN(EPPHI*PI/180.)*AIMAG(X1)-EPEPS*
c     *         SIN(EPPHI*PI/90.)*AIMAG(X2))
c     *       -PY*(AIMAG(Y1)+EPEPS*COS(EPPHI*PI/90.)*AIMAG(Y2)
c     *          +2.*EPEPS*AIMAG(Y3)+XILA*COS(EPPHI*PI/180.)*AIMAG(Y4))
c     *         +PZ*(EPEPS*SIN(EPPHI*PI/90.)*AIMAG(Z2)+
c     *              XILA*SIN(EPPHI*PI/180.)*AIMAG(Z1))
c            SIGET=-POL_ELEC*(-PX*(SQRT(2.*EPEPS*(1-EPEPS))*COS(EPPHI*PI/180.)
c     *               *REAL(X1)+SQRT(1-EPEPS**2)*REAL(X2))
c     *       +PY*SQRT(2.*EPEPS*(1-EPEPS))*SIN(EPPHI*PI/180.)*REAL(Y4)
c     *       +PZ*(SQRT(1-EPEPS**2)*REAL(Z2)+SQRT(2*EPEPS*(1-EPEPS))*
c     *        COS(EPPHI*PI/180.)*REAL(Z1)))
c           SIGMAO = SIGMAO + SIGE + SIGPT + SIGET
c           SIGMAO=SIGMAO*FKT
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C         The target is polarized along the beam axis.
c                                       PS = -1 anti-parallel to beam
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
           IF(EPQ2.LT.1.E-06.OR.EPEPS.LT.1.E-06)THEN
                THETA_GAM = 0.
           ELSE
           ENUE = (EPW**2 - XMPROT**2 + EPQ2)/2./XMPROT
           THETA_ELEC = 2.*ATAN(SQRT(EPQ2/2./(EPQ2+ENUE**2)*
     *                  (1.-EPEPS)/EPEPS))
           E_ELEC_SCATT = -ENUE/2 + SQRT((ENUE/2.)**2 +
     *                   EPQ2/4./SIN(THETA_ELEC/2.)**2)
           E_ELEC_PRIM = E_ELEC_SCATT + ENUE
           THETA_GAM = ATAN(SIN(THETA_ELEC)/
     *                 (E_ELEC_PRIM/E_ELEC_SCATT - COS(THETA_ELEC)))
           THETA_GAM = THETA_GAM*180./PI
           ENDIF
111        format (5f10.4)
           POL_X = POL_TARG*SIND(THETA_GAM)*COSD(EPPHI)
           POL_Y = POL_TARG*SIND(THETA_GAM)*SIND(EPPHI)
           POL_Z = POL_TARG*COSD(THETA_GAM)

           SIGPT = - POL_X*(SIND(EPPHI)*
     1  SQRT(2.*EPEPS*(1. + EPEPS))*AIMAG(X1) +
     2  EPEPS*SIND(2.*EPPHI)*AIMAG(X2)) -
     3  POL_Y*(AIMAG(Y1) +
     4  EPEPS*COSD(2.*EPPHI)*AIMAG(Y2) +
     5  2.*EPEPS*AIMAG(Y3) +
     6  SQRT(2.*EPEPS*(1.+EPEPS))*COSD(EPPHI)*AIMAG(Y4)) +
     7  POL_Z*(EPEPS*SIND(2.*EPPHI)*AIMAG(Z2) +
     8  SQRT(2.*EPEPS*(1. + EPEPS))*SIND(EPPHI)*AIMAG(Z1))

            SIGET=-POL_ELEC*(-POL_X*(SQRT(2.*EPEPS*(1.-EPEPS))*
     *    COSD(EPPHI)*REAL(X1)+SQRT(1.-EPEPS**2)*REAL(X2))
     *  +POL_Y*SQRT(2.*EPEPS*(1.-EPEPS))*SIND(EPPHI)*REAL(Y4)
     *  +POL_Z*(SQRT(1.-EPEPS**2)*REAL(Z2)+SQRT(2.*EPEPS*(1.-EPEPS))*
     *   COSD(EPPHI)*REAL(Z1)))

ccc        SIGET = POL_ELEC*POL_TARG*(SIND(THETA_GAM)*(COSD(EPPHI))**2*
ccc     *  (REAL(Y4)*sqrt(2.*EPEPS*(1.-EPEPS) +
CCC     *  SQRT(2.*EPEPS*(1.-EPEPS))*REAL(X1))+
ccc     *  COSD(EPPHI)*(SIND(THETA_GAM)*SQRT(1.-EPEPS**2)*REAL(X2) -
ccc     *  COSD(THETA_GAM)*SQRT(2.*EPEPS*(1.-EPEPS))*REAL(Z1)) -
ccc     *  SIND(THETA_GAM)*REAL(Y4) - COSD(THETA_GAM)*
ccc     *  SQRT(1.-EPEPS**2)*REAL(Z2))

           SIGMAO = SIGMAO + SIGE + SIGPT + SIGET
           SIGMAO = SIGMAO*FKT

           RETURN
           END





      
*****************************************************************************
*****************************************************************************

      function sind(x)
      data rad/.01745329/

      sind=sin(rad*x)
      return
      end
*****************************************************************************
*****************************************************************************


*****************************************************************************
*****************************************************************************

      function cosd(x)
      data rad/.01745329/

      cosd=cos(rad*x)
      return
      end
*****************************************************************************
*****************************************************************************



*****************************************************************************
*****************************************************************************

      function tand(x)
      DATA PI  /3.1415926  /
      tand=tan(pi*x/180)
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************
      
      function asind(x)
      DATA PI /3.1415926  /
      asind=180*asin(x)/pi
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************
      
      function acosd(x)
      DATA PI /3.1415926  /
      acosd=180*acos(x)/pi
      return
      end
*****************************************************************************
*****************************************************************************

*****************************************************************************
*****************************************************************************
      
      function atand(x)
      DATA PI /3.1415926  /
      atand=180.*atan(x)/pi
      return
      end
*****************************************************************************
*****************************************************************************
