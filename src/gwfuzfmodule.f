      MODULE GWFUZFMODULE
        CHARACTER(LEN=64) :: Version_uzf
        REAL,PARAMETER :: CLOSEZERO=1.0E-15
        DOUBLE PRECISION,PARAMETER :: NEARZERO=1.0D-30
        DOUBLE PRECISION,PARAMETER :: ZEROD15=1.0D-15, ZEROD9=1.0D-09
        DOUBLE PRECISION,PARAMETER :: ZEROD6=1.0D-06
        DOUBLE PRECISION,PARAMETER :: ZEROD7=1.0D-07
        INTEGER         ,PARAMETER :: IRUNBIG = 10000
        INTEGER,SAVE,POINTER :: NUMCELLS, TOTCELLS, Iseepsupress, IPRCNT
        INTEGER,SAVE,POINTER :: ITHTIFLG, ITHTRFLG
        DOUBLE PRECISION,SAVE :: THETAB, FLUXB, FLUXHLD2
        DOUBLE PRECISION,SAVE,DIMENSION(:),POINTER :: CHECKTIME
        INTEGER,SAVE,DIMENSION(:),POINTER :: MORE
        INTEGER,SAVE,DIMENSION(:,:),POINTER :: LAYNUM
        INTEGER,SAVE,POINTER   ::NUZTOP, IUZFOPT, IRUNFLG, IETFLG, IUZM
        INTEGER,SAVE,POINTER   ::IUZFCB1, IUZFCB2, NTRAIL, NWAV, NSETS
        INTEGER,SAVE,POINTER   ::IUZFB22, IUZFB11
        INTEGER,SAVE,POINTER   ::NUZGAG, NUZGAGAR, NUZCL, NUZRW, IGSFLOW
        INTEGER,SAVE,POINTER   ::RTSOLUTE
        INTEGER,SAVE,  DIMENSION(:),    POINTER :: ITRLSTH
        INTEGER,SAVE,  DIMENSION(:,:),  POINTER :: IRUNBND, IUZFBND
        INTEGER,SAVE,  DIMENSION(:,:),  POINTER :: IUZLIST, NWAVST
        INTEGER,SAVE,  DIMENSION(:,:),  POINTER :: IUZHOLD
        INTEGER,SAVE,  DIMENSION(:,:),  POINTER :: LTRLIT, LTRLST
        INTEGER,SAVE,  DIMENSION(:,:),  POINTER :: ITRLIT, ITRLST
        REAL,   SAVE,POINTER   ::TOTRUNOFF, SURFDEP
        REAL,   SAVE,  DIMENSION(:),    POINTER :: FBINS
        REAL,   SAVE,  DIMENSION(:,:),  POINTER :: SEEPOUT, EXCESPP, VKS
        REAL,   SAVE,  DIMENSION(:,:),  POINTER :: AIR_ENTRY, H_ROOT
        REAL,   SAVE,  DIMENSION(:,:),  POINTER :: REJ_INF
        REAL,   SAVE,  DIMENSION(:,:),  POINTER :: TO_CFP
        REAL,   SAVE,  DIMENSION(:,:),  POINTER :: EPS, THTS, THTI, THTR
        REAL,   SAVE,  DIMENSION(:,:),  POINTER :: PETRATE, ROOTDPTH
        REAL,   SAVE,  DIMENSION(:,:),  POINTER :: WCWILT, FINF
        REAL,   SAVE,  DIMENSION(:,:),  POINTER :: UZFETOUT, GWET
        REAL,   SAVE,  DIMENSION(:,:),  POINTER :: FNETEXFIL
        DOUBLE PRECISION, SAVE, DIMENSION(:),  POINTER :: CUMUZVOL 
        DOUBLE PRECISION, SAVE, DIMENSION(:),  POINTER :: UZTSRAT
        DOUBLE PRECISION, SAVE, DIMENSION(:,:),POINTER :: UZFLWT, UZSTOR
        DOUBLE PRECISION, SAVE, DIMENSION(:,:),POINTER :: UZDPIT, UZDPST
        DOUBLE PRECISION, SAVE, DIMENSION(:,:),POINTER :: UZTHIT, UZTHST
        DOUBLE PRECISION, SAVE, DIMENSION(:,:),POINTER :: UZSPIT, UZSPST
        DOUBLE PRECISION, SAVE, DIMENSION(:,:),POINTER :: UZFLIT, UZFLST
        DOUBLE PRECISION, SAVE, DIMENSION(:,:),POINTER :: DELSTOR
        DOUBLE PRECISION, SAVE, DIMENSION(:,:),POINTER :: UZOLSFLX
        DOUBLE PRECISION, SAVE, DIMENSION(:,:),POINTER :: HLDUZF
        DOUBLE PRECISION, SAVE, DIMENSION(:,:,:),POINTER :: UZTOTBAL
        DOUBLE PRECISION, SAVE, DIMENSION(:,:,:),POINTER :: GRIDSTOR
        DOUBLE PRECISION, SAVE, DIMENSION(:,:,:),POINTER :: GRIDET
      TYPE GWFUZFTYPE
        INTEGER,     POINTER   ::NUZTOP, IUZFOPT, IRUNFLG, IETFLG, IUZM
        INTEGER,     POINTER   ::IUZFCB1, IUZFCB2, NTRAIL, NWAV, NSETS
        INTEGER,     POINTER   ::IUZFB22, IUZFB11
        INTEGER,     POINTER   ::NUZGAG, NUZGAGAR, NUZCL, NUZRW, IGSFLOW
        INTEGER,     POINTER   ::RTSOLUTE
        INTEGER,     POINTER   ::NUMCELLS, TOTCELLS, Iseepsupress,IPRCNT
        INTEGER,     POINTER   ::ITHTIFLG, ITHTRFLG
        DOUBLE PRECISION,DIMENSION(:),  POINTER :: CHECKTIME
        INTEGER,DIMENSION(:),POINTER :: MORE
        INTEGER,DIMENSION(:,:),POINTER :: LAYNUM
        INTEGER,       DIMENSION(:),    POINTER :: ITRLSTH
        INTEGER,       DIMENSION(:,:),  POINTER :: IRUNBND, IUZFBND
        INTEGER,       DIMENSION(:,:),  POINTER :: IUZLIST, NWAVST
        INTEGER,       DIMENSION(:,:),  POINTER :: IUZHOLD
        INTEGER,       DIMENSION(:,:),  POINTER :: LTRLIT, LTRLST
        INTEGER,       DIMENSION(:,:),  POINTER :: ITRLIT, ITRLST
        REAL,          POINTER            ::TOTRUNOFF, SURFDEP
        REAL,          DIMENSION(:),    POINTER :: FBINS
        REAL,          DIMENSION(:,:),  POINTER :: SEEPOUT, EXCESPP, VKS
        REAL,          DIMENSION(:,:),  POINTER :: AIR_ENTRY, H_ROOT
        REAL,          DIMENSION(:,:),  POINTER :: REJ_INF
        REAL,          DIMENSION(:,:),  POINTER :: TO_CFP
        REAL,          DIMENSION(:,:),  POINTER :: EPS, THTS, THTI, THTR
        REAL,          DIMENSION(:,:),  POINTER :: PETRATE, ROOTDPTH
        REAL,          DIMENSION(:,:),  POINTER :: WCWILT, FINF
        REAL,          DIMENSION(:,:),  POINTER :: UZFETOUT, GWET
        REAL,          DIMENSION(:,:),  POINTER :: FNETEXFIL
        DOUBLE PRECISION,       DIMENSION(:),  POINTER :: CUMUZVOL
        DOUBLE PRECISION,       DIMENSION(:),  POINTER :: UZTSRAT
        DOUBLE PRECISION,       DIMENSION(:,:),POINTER :: UZFLWT, UZSTOR
        DOUBLE PRECISION,       DIMENSION(:,:),POINTER :: UZDPIT, UZDPST
        DOUBLE PRECISION,       DIMENSION(:,:),POINTER :: UZTHIT, UZTHST
        DOUBLE PRECISION,       DIMENSION(:,:),POINTER :: UZSPIT, UZSPST
        DOUBLE PRECISION,       DIMENSION(:,:),POINTER :: UZFLIT, UZFLST
        DOUBLE PRECISION,       DIMENSION(:,:),POINTER :: DELSTOR
        DOUBLE PRECISION,       DIMENSION(:,:),POINTER :: UZOLSFLX
        DOUBLE PRECISION,       DIMENSION(:,:),POINTER :: HLDUZF
        DOUBLE PRECISION,       DIMENSION(:,:,:),POINTER :: UZTOTBAL
        DOUBLE PRECISION,       DIMENSION(:,:,:),POINTER :: GRIDSTOR
        DOUBLE PRECISION,       DIMENSION(:,:,:),POINTER :: GRIDET
      END TYPE
      TYPE(GWFUZFTYPE), SAVE:: GWFUZFDAT(10)
      END MODULE GWFUZFMODULE
