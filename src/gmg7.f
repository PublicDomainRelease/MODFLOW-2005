      MODULE GMGMODULE
        INTEGER,SAVE,POINTER  ::IITER,IADAMPGMG,ISM,ISC,IOUTGMG
        INTEGER,SAVE,POINTER  ::ISIZ,IPREC,IIOUT
        INTEGER,SAVE,POINTER  ::SITER,TSITER
        INTEGER,SAVE,POINTER  ::GMGID
        REAL   ,SAVE,POINTER  ::HCLOSEGMG,RCLOSEGMG,DAMPGMG
        DOUBLE PRECISION,SAVE,POINTER  :: RELAXGMG
      TYPE GMGTYPE
        INTEGER,POINTER  ::IITER,IADAMPGMG,ISM,ISC,IOUTGMG
        INTEGER,POINTER  ::ISIZ,IPREC,IIOUT
        INTEGER,POINTER  ::SITER,TSITER
        INTEGER,POINTER  ::GMGID
        REAL   ,POINTER  ::HCLOSEGMG,RCLOSEGMG,DAMPGMG
        DOUBLE PRECISION,POINTER  :: RELAXGMG
      END TYPE
      TYPE(GMGTYPE), SAVE ::GMGDAT(10)
      END MODULE GMGMODULE
C
      SUBROUTINE GMG7AR(IN,MXITER,IGRID)
C--------------------------------------------------------------------
C     READS INPUT FROM FILE TYPE GMG SPECIFIED IN NAME FILE
C     ALLOCATES GMG SOLVER
C     EXPLICIT DECLERATIONS
C--------------------------------------------------------------------
      USE GLOBAL,   ONLY:IOUT,NCOL,NROW,NLAY
      USE GMGMODULE,ONLY:IITER,IADAMPGMG,ISM,ISC,IOUTGMG,ISIZ,
     1                   IPREC,IIOUT,SITER,TSITER,GMGID,HCLOSEGMG,
     2                   RCLOSEGMG,DAMPGMG,RELAXGMG
      IMPLICIT NONE
      CHARACTER*200 LINE
      INTEGER IN,MXITER,IGRID,IERR
C
C--------------------------------------------------------------------
C     ALLOCATE POINTERS 
C--------------------------------------------------------------------
      ALLOCATE(IITER,IADAMPGMG,ISM,ISC,IOUTGMG,ISIZ,IPREC,IIOUT,
     1         SITER,TSITER,GMGID)
      ALLOCATE(HCLOSEGMG,RCLOSEGMG,DAMPGMG,RELAXGMG)
C
C--------------------------------------------------------------------
C     READ AND PRINT COMMENTS
C--------------------------------------------------------------------
      CALL URDCOM(IN,IOUT,LINE)
C
C--------------------------------------------------------------------
C     READ INPUT FILE
C--------------------------------------------------------------------
      READ(LINE,*) RCLOSEGMG,IITER,HCLOSEGMG,MXITER
      CALL URDCOM(IN,IOUT,LINE)
      READ(LINE,*) DAMPGMG,IADAMPGMG,IOUTGMG
      CALL URDCOM(IN,IOUT,LINE)
      READ(LINE,*) ISM,ISC
C
      SITER=0
      TSITER=0
      RELAXGMG=0.0D0
      IF(ISC .EQ. 4) THEN
        CALL URDCOM(IN,IOUT,LINE)
        READ(LINE,*) RELAXGMG
      END IF
C
      IF(DAMPGMG .LE. 0.0 .OR. DAMPGMG .GT. 1.0) DAMPGMG=1.0
      IIOUT=IOUT
      IF(IOUTGMG .GT. 2) IIOUT=6
C
C--------------------------------------------------------------------
C     ALLOCATE
C--------------------------------------------------------------------
C
C---- CHECK FOR FORCED DOUBLE PRECISION
C
      IPREC=0
      IF(KIND(DAMPGMG) .EQ. 8) IPREC=1
C
      CALL MF2KGMG_ALLOCATE(GMGID,NCOL,NROW,NLAY,IPREC,ISM,ISC,
     &                      RELAXGMG,ISIZ,IERR)
      IF(IERR .NE. 0) THEN
        CALL USTOP('ALLOCATION ERROR IN SUBROUTINE GMG1ALG')
      END IF
C
      WRITE(IIOUT,500) RCLOSEGMG,IITER,HCLOSEGMG,MXITER,
     &                 DAMPGMG,IADAMPGMG,IOUTGMG,
     &                 ISM,ISC,RELAXGMG
C
      IF(IADAMPGMG .NE. 0) WRITE(IIOUT,510)
      IF(ISM .EQ. 0) WRITE(IIOUT,520)
      IF(ISM .EQ. 1) WRITE(IIOUT,525)
      IF(ISC .EQ. 0) WRITE(IIOUT,530)
      IF(ISC .EQ. 1) WRITE(IIOUT,531)
      IF(ISC .EQ. 2) WRITE(IIOUT,532)
      IF(ISC .EQ. 3) WRITE(IIOUT,533)
      IF(ISC .EQ. 4) WRITE(IIOUT,534)
C
      WRITE(IIOUT,540) ISIZ
C
C--------------------------------------------------------------------
C     FORMAT STATEMENTS
C--------------------------------------------------------------------
  500 FORMAT(1X,'-------------------------------------------------',/,
     &       1X,'GMG -- PCG GEOMETRIC MULTI-GRID SOLUTION PACKAGE:',/,
     &       1X,'-------------------------------------------------',/,
     &       1X,'RCLOSE  = ',1P,E8.2,'; INNER CONVERGENCE CRITERION',/,
     &       1X,'IITER   = ',I8,'; MAX INNER ITERATIONS            ',/,
     &       1X,'HCLOSE  = ',1P,E8.2,'; OUTER CONVERGENCE CRITERION',/,
     &       1X,'MXIITER = ',I8,'; MAX OUTER ITERATIONS            ',/,
     &       1X,'DAMP    = ',1P,E8.2,'; DAMPING PARAMETER          ',/,
     &       1X,'IADAMP  = ',I8,'; ADAPTIVE DAMPING FLAG           ',/,
     &       1X,'IOUTGMG = ',I8,'; OUTPUT CONTROL FLAG             ',/,
     &       1X,'ISM     = ',I8,'; SMOOTHER FLAG                   ',/,
     &       1X,'ISC     = ',I8,'; COARSENING FLAG                 ',/,
     &       1X,'RELAX   = ',1P,E8.2,'; RELAXATION PARAMETER       ',/,
     &       1X,"-------------------------------------------------")
C
  510 FORMAT(1X,"COOLEY'S ADAPTIVE DAMPING METHOD IMPLEMENTED")
  520 FORMAT(1X,'ILU SMOOTHING IMPLEMENTED')
  525 FORMAT(1X,'SGS SMOOTHING IMPLEMENTED')
C
  530 FORMAT(1X,'FULL COARSENING')
  531 FORMAT(1X,'COARSENING ALONG COLUMNS AND ROWS ONLY')
  532 FORMAT(1X,'COARSENING ALONG ROWS AND LAYERS ONLY')
  533 FORMAT(1X,'COARSENING ALONG COLUMNS AND LAYERS ONLY')
  534 FORMAT(1X,'NO COARSENING')
C
  540 FORMAT(1X,'-------------------------------------------------',/,
     &       1X,I4,' MEGABYTES OF MEMORY ALLOCATED BY GMG',/,
     &       1X,'-------------------------------------------------',/)
C
      CALL GMG7PSV(IGRID)
      RETURN
      END
C***********************************************************************
      SUBROUTINE GMG7AP(HNEW,RHS,CR,CC,CV,HCOF,HNOFLO,IBOUND,
     &                  IITER,MXITER,RCLOSE,HCLOSE,
     &                  KITER,KSTP,KPER,
     &                  ICNVG,SITER,TSITER,DAMP,IADAMP,
     &                  IOUTGMG,IOUT,GMGID)
C***********************************************************************
C     GMG7AP CALLS THE FOLLOWING FUNCTIONS FROM THE GMG LIBRARY:
C
C    MF2KMG_ASSEMBLE:
C      -- ASSEMBLES CELL-CENTERED FINITE DIFFERENCE
C         MATRIX AND MULTIGRID PRECONDITIONED CONJUGATE
C         GRADIENT SOLVER.  RETURNS L2-NORM OF RESIDUAL BIGR0
C
C    MF2KMG_EVAL:
C      -- APPROXIMATES THE HEAD CHANGE E FOR A*E=R WHERE A IS THE
C         CCFD MATRIX, AND R IS THE INITIAL RESIDUAL.
C         RETURNS THE L2-NORM OF THE RESIDUAL R-A*E.
C
C    MF2KGMG_BIGH:
C      -- COMPUTES MAX HEAD CHANGE BIGH AND RETURNS LOCATION
C         (COL,ROW,LAY) OF MAX HEAD CHANGE.  ABSOLUTE VALUE
C         OF BIGH IS MAX-NORM OF HEAD CHANGE.
C
C    MF2KMG_UPDATE:
C      -- ADDS THE CORRECTION TO THE HEADS HNEW=HNEW+DAMP*E.
C--------------------------------------------------------------------
      IMPLICIT NONE
      REAL RHS(*),CR(*),CC(*),CV(*),HCOF(*)
      REAL HNOFLO,RCLOSE,HCLOSE,DAMP
      DOUBLE PRECISION HNEW(*)
      INTEGER IBOUND(*)
      INTEGER MXITER,IITER,KITER,KSTP,KPER,ICNVG,IOUTGMG,IOUT 
      INTEGER SITER,TSITER
      INTEGER IADAMP
C
      INTEGER IIOUT,GMGID
      DOUBLE PRECISION DRCLOSE
C
      INTEGER ITER,IERR
      INTEGER JBIGH,IBIGH,KBIGH
C
      DOUBLE PRECISION BIGH,BIGR,BIGR0
      DOUBLE PRECISION S,DH,DAMP0
      DOUBLE PRECISION, SAVE ::BIGH0
      DOUBLE PRECISION, SAVE ::DDAMP
C--------------------------------------------------------------------
C
C---- INITIALIZE VARIABLES
C
      ICNVG=0
      IIOUT=IOUT
      IF(IOUTGMG .GT. 2) IIOUT=6
      IF(KITER .EQ. 1) DDAMP=DAMP
      DAMP0=DDAMP
C
C--------------------------------------------------------------------
C     ASSEMBLE SOLVER
C--------------------------------------------------------------------
      CALL MF2KGMG_ASSEMBLE(GMGID,BIGR0,CR,CC,CV,HCOF,HNEW,RHS,HNOFLO,
     &                      IBOUND,IERR)
      IF(IERR .NE. 0) THEN
        CALL USTOP('GMG ASSEMBLY ERROR IN SUBROUTINE GMG1AP')
      END IF
C
C--------------------------------------------------------------------
C     SCALE CLOSURE CRITERION FOR INNER ITERATION BASED ON CURRENT
C     VALUE OF DAMPING AND INITIAL RESIDUAL.
C--------------------------------------------------------------------
      DRCLOSE=DDAMP*RCLOSE+(1.0D0-DDAMP)*BIGR0
C
C--------------------------------------------------------------------
C     COMPUTE HEAD CHANGE
C--------------------------------------------------------------------
      CALL MF2KGMG_EVAL(GMGID,ITER,BIGR,DRCLOSE,IITER,IOUTGMG,IIOUT)
      SITER=SITER+ITER
C
C--------------------------------------------------------------------
C     COMPUTE MAX HEAD CHANGE
C--------------------------------------------------------------------
      CALL MF2KGMG_BIGH(GMGID,BIGH,JBIGH,IBIGH,KBIGH)
C
C--------------------------------------------------------------------
C     CHECK FOR CLOSURE
C--------------------------------------------------------------------
       IF(MXITER .EQ. 1 .AND. BIGR .LE. RCLOSE) THEN
         DDAMP=1.0;
         ICNVG=1
         GOTO 100
       END IF
C
       IF(ABS(BIGH) .LE. HCLOSE .AND. BIGR .LE. RCLOSE) THEN
         DDAMP=1.0;
         ICNVG=1
         GOTO 100
       END IF
C
C--------------------------------------------------------------------------
C     ADJUST DAMPING PARAMETER USING COOLEY'S METHOD
C--------------------------------------------------------------------------
      IF(IADAMP .NE. 0) THEN
        IF(KITER .GT. 1) THEN
          DH=BIGH/BIGH0
          S=DH/DDAMP
          IF(S .GE. -1.0D0) THEN
            DDAMP=(3.0D0+S)/(3.0D0+ABS(S))
          ELSE
            DDAMP=0.5D0/ABS(S)
          END IF
          IF(DDAMP .LT. DAMP) DDAMP=DAMP
        END IF
      END IF
C
C--------------------------------------------------------------------
C     ADD CORECTION
C--------------------------------------------------------------------
  100 CONTINUE
      CALL MF2KGMG_UPDATE(GMGID,HNEW,DDAMP)
      BIGH0=BIGH
C
      IF(IOUTGMG .NE. 0) THEN
        WRITE(IIOUT,510) ITER,DDAMP,BIGR,
     &                   ABS(BIGH),JBIGH,IBIGH,KBIGH
        IF(ICNVG .EQ. 1 ) THEN
          TSITER=TSITER+SITER
          WRITE(IIOUT,500) KSTP,KPER,KITER,SITER,TSITER
          SITER=0
        END IF
      END IF
C
C--------------------------------------------------------------------
C     FORMAT STATEMENTS
C--------------------------------------------------------------------
C
  500 FORMAT(1X,'-------------------------------',/,
     &       1X,'TIME STEP            : ',I6,/,
     &       1X,'STRESS PERIOD        : ',I6,/,
     &       1X,'GMG CALLS            : ',I6,/,
     &       1X,'PCG ITERATIONS       : ',I6,/,
     &       1X,'-------------------------------',/,
     &       1X,'TOTAL PCG ITERATIONS : ',I6,/,
     &       1X,'-------------------------------',/)
C
  510 FORMAT(1X,'-------------------------------------',
     &          '--------------------',/,
     &       1X,'PCG ITERATIONS                    : ',I4,/,
     &       1X,'DAMPING                           : ',0P,F5.3,/,
     &       1X,'L2-NORM OF RESIDUAL               : ',1P,E10.4,/,
     &       1X,'MAX HEAD CHANGE                   : ',1P,E10.4,/,
     &       1X,'MAX HEAD CHANGE AT (COL,ROW,LAY)  : (',
     &           I6,',',I6,',',I6,')',/,
     &       1X,'-------------------------------------',
     &          '--------------------',/)
C
      RETURN
      END

C***********************************************************************
C     SUBROUTINE RESPRINT IS CALLED FROM THE SOLVER FOR OUTPUT
C     OF REDUCTION HISTORY.
C
C     THE ITERATION (I), THE RESIDUAL (RES), AND THE
C     CONVERGENCE FACTOR (CFAC) ARE PRINTED.
C***********************************************************************
      SUBROUTINE RESPRINT(IOUT,I,RES,CFAC)
      IMPLICIT NONE
      INTEGER IOUT,I
      DOUBLEPRECISION RES,CFAC
C
C---- PRINT RESIDUALS
C
      WRITE(IOUT,100) I,RES,CFAC
C
C--------------------------------------------------------------------
C     FORMAT STATEMENTS
C--------------------------------------------------------------------
  100 FORMAT(1X,'ITER:',I4,
     &       2X,'RES: ',1P,E10.4,
     &       2X,'CFAC: ',0P,F5.3)
C
      RETURN
      END
C
      SUBROUTINE GMG7DA(IGRID)
C  Deallocate GMG data 
      USE GMGMODULE
      CALL GMG7PNT(IGRID)
      CALL MF2KGMG_FREE(GMGID)
      DEALLOCATE(IITER,IADAMPGMG,ISM,ISC,IOUTGMG,ISIZ,IPREC,IIOUT,
     1           SITER,TSITER,GMGID)
      DEALLOCATE(HCLOSEGMG,RCLOSEGMG,DAMPGMG,RELAXGMG)
C
      RETURN
      END
C
      SUBROUTINE GMG7PNT(IGRID)
C  Set pointers to GMG data for a grid
      USE GMGMODULE
C
      IITER=>GMGDAT(IGRID)%IITER
      IADAMPGMG=>GMGDAT(IGRID)%IADAMPGMG
      ISM=>GMGDAT(IGRID)%ISM
      ISC=>GMGDAT(IGRID)%ISC
      IOUTGMG=>GMGDAT(IGRID)%IOUTGMG
      ISIZ=>GMGDAT(IGRID)%ISIZ
      IPREC=>GMGDAT(IGRID)%IPREC
      IIOUT=>GMGDAT(IGRID)%IIOUT
      SITER=>GMGDAT(IGRID)%SITER
      TSITER=>GMGDAT(IGRID)%TSITER
      GMGID=>GMGDAT(IGRID)%GMGID
      HCLOSEGMG=>GMGDAT(IGRID)%HCLOSEGMG
      RCLOSEGMG=>GMGDAT(IGRID)%RCLOSEGMG
      DAMPGMG=>GMGDAT(IGRID)%DAMPGMG
      RELAXGMG=>GMGDAT(IGRID)%RELAXGMG
C
      RETURN
      END
C
      SUBROUTINE GMG7PSV(IGRID)
C  Save pointers to GMG data
      USE GMGMODULE
C
      GMGDAT(IGRID)%IITER=>IITER
      GMGDAT(IGRID)%IADAMPGMG=>IADAMPGMG
      GMGDAT(IGRID)%ISM=>ISM
      GMGDAT(IGRID)%ISC=>ISC
      GMGDAT(IGRID)%ISIZ=>ISIZ
      GMGDAT(IGRID)%IOUTGMG=>IOUTGMG
      GMGDAT(IGRID)%IPREC=>IPREC
      GMGDAT(IGRID)%IIOUT=>IIOUT
      GMGDAT(IGRID)%SITER=>SITER
      GMGDAT(IGRID)%TSITER=>TSITER
      GMGDAT(IGRID)%GMGID=>GMGID
      GMGDAT(IGRID)%HCLOSEGMG=>HCLOSEGMG
      GMGDAT(IGRID)%RCLOSEGMG=>RCLOSEGMG
      GMGDAT(IGRID)%DAMPGMG=>DAMPGMG
      GMGDAT(IGRID)%RELAXGMG=>RELAXGMG
C
      RETURN
      END

