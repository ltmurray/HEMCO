!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_ch4_isotopologues_mod.F90
!
! !DESCRIPTION: Partitions hierarchy-resolved, sector-specific bulk methane
!  emissions into stable methane isotopologue emissions. Carbon-13 and
!  deuterium signatures are supplied as delta values relative to VPDB and
!  VSMOW. Molecular weights come from the HEMCO species registry so all
!  mass-to-mole conversions remain consistent with species_database.yml.
!\
!\
! !INTERFACE:
!
MODULE HCOX_CH4_Isotopologues_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Calc_Mod,    ONLY : HCO_EvalFld
  USE HCO_State_Mod,   ONLY : HCO_State, HCO_GetHcoID
  USE HCO_Types_Mod,   ONLY : FileData, HCO_UFLAG_ALWAYS, HCO_UFLAG_3HR
  USE HCOX_State_Mod,  ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: HCOX_CH4_Isotopologues_Init
  PUBLIC :: HCOX_CH4_Isotopologues_Run
  PUBLIC :: HCOX_CH4_Isotopologues_Final

  INTEGER, PARAMETER :: nSector       = 16
  INTEGER, PARAMETER :: nIsotopologue = 4

  CHARACTER(LEN=6), PARAMETER :: Sector(nSector) = [                       &
       'OIL   ', 'GAS   ', 'COL   ', 'LIV_EF', 'LDF   ', 'WST   ',        &
       'RIC   ', 'OTA   ', 'BBN   ', 'WTL   ', 'SEE   ', 'LAK   ',        &
       'TER   ', 'OCN   ', 'RSV   ', 'LIV_MM' ]

  INTEGER, PARAMETER :: Category(nSector) = [                              &
       1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 ]

  CHARACTER(LEN=7), PARAMETER :: Isotopologue(nIsotopologue) = [           &
       'C12H4  ', 'C13H4  ', 'C12H3D ', 'C13H3D ' ]

  ! Defaults overridden by extension 115 options in HEMCO_Config.rc
  REAL(hp) :: R13_VPDB = 0.0111802_hp
  REAL(hp) :: RD_VSMOW = 1.5576e-4_hp

  ! Cache the sector-specific molecular fractions.  Signature fields are
  ! typically annual, whereas HEMCO emissions are evaluated every emission
  ! timestep.  Recalculate a sector only when either signature field changes.
  REAL(hp), ALLOCATABLE :: D13C_Cache(:,:,:)
  REAL(hp), ALLOCATABLE :: DD_Cache(:,:,:)
  REAL(hp), ALLOCATABLE :: Fraction_Cache(:,:,:,:)
  LOGICAL               :: CacheValid(nSector) = .FALSE.
  INTEGER               :: UpdateInterval(nSector) = 0
  INTEGER               :: LastYear(nSector)       = -1
  INTEGER               :: LastMonth(nSector)      = -1
  INTEGER               :: LastDay(nSector)        = -1
  INTEGER               :: LastHour(nSector)       = -1
  INTEGER               :: LastMinute(nSector)     = -1

CONTAINS
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_CH4_Isotopologues_Init
!
! !DESCRIPTION: Initializes the stable methane isotopologue extension.
!\
!\
! !INTERFACE:
!
  SUBROUTINE HCOX_CH4_Isotopologues_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod, ONLY : GetExtNr, GetExtOpt
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState
    CHARACTER(LEN=*), INTENT(IN)   :: ExtName
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER       :: ExtState
    INTEGER, INTENT(INOUT)         :: RC
!EOP
!------------------------------------------------------------------------------
!BOC

    INTEGER            :: AllocStat
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    RC      = HCO_SUCCESS
    ThisLoc = ' -> at HCOX_CH4_Isotopologues_Init'
    ExtState%CH4_Isotopologues = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )

    CALL GetExtOpt( HcoState%Config, ExtState%CH4_Isotopologues,             &
                    'R13_VPDB', OptValHp=R13_VPDB, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error reading the R13_VPDB extension option'
       CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    CALL GetExtOpt( HcoState%Config, ExtState%CH4_Isotopologues,             &
                    'RD_VSMOW', OptValHp=RD_VSMOW, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error reading the RD_VSMOW extension option'
       CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( R13_VPDB <= 0.0_hp .OR. RD_VSMOW <= 0.0_hp ) THEN
       ErrMsg = 'Reference isotope ratios must both be positive'
       CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( ALLOCATED(D13C_Cache)     ) DEALLOCATE(D13C_Cache)
    IF ( ALLOCATED(DD_Cache)       ) DEALLOCATE(DD_Cache)
    IF ( ALLOCATED(Fraction_Cache) ) DEALLOCATE(Fraction_Cache)
    IF ( ALLOCATED(HcoState%CategoryFlux) ) THEN
       DEALLOCATE(HcoState%CategoryFlux)
    ENDIF
    ALLOCATE( D13C_Cache(HcoState%NX,HcoState%NY,nSector),              &
              DD_Cache(HcoState%NX,HcoState%NY,nSector),                &
              Fraction_Cache(HcoState%NX,HcoState%NY,nIsotopologue,     &
                             nSector),                                  &
              HcoState%CategoryFlux(HcoState%NX,HcoState%NY,            &
                                    HcoState%NZ,nSector),                &
              STAT=AllocStat )
    IF ( AllocStat /= 0 ) THEN
       ErrMsg = 'Could not allocate CH4 isotopologue fraction cache'
       CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    CacheValid = .FALSE.
    HcoState%CategoryFluxHcoID = HCO_GetHcoID( 'CH4', HcoState )
    UpdateInterval = 0
    LastYear   = -1
    LastMonth  = -1
    LastDay    = -1
    LastHour   = -1
    LastMinute = -1

  END SUBROUTINE HCOX_CH4_Isotopologues_Init
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_CH4_Isotopologues_Run
!
! !DESCRIPTION: Reassembles bulk CH4 emissions one sector at a time, converts
!  source signatures to molecular isotopologue fractions, and adds the
!  corresponding mass flux to every isotopologue present in the model.
!  Bulk CH4 emissions are left unchanged.
!\
!\
! !INTERFACE:
!
  SUBROUTINE HCOX_CH4_Isotopologues_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_Clock_Mod,    ONLY : HcoClock_Get
    USE HCO_DataCont_Mod, ONLY : ListCont_Find
    USE HCO_FluxArr_Mod, ONLY : HCO_EmisAdd
    USE HCO_Types_Mod,   ONLY : ListCont
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER       :: ExtState
    TYPE(HCO_State), POINTER       :: HcoState
    INTEGER, INTENT(INOUT)         :: RC
!EOP
!------------------------------------------------------------------------------
!BOC

    INTEGER            :: idCH4, idIso(nIsotopologue)
    INTEGER            :: I, J, L, N, S, SigIdx, ThisCat
    INTEGER            :: YYYY, MM, DDnow, HH, MN
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    CHARACTER(LEN=63)  :: FieldName
    LOGICAL            :: Found
    LOGICAL            :: RefreshCache, SignatureLoaded(nSector)
    TYPE(ListCont), POINTER :: D13C_Lct, DD_Lct
    REAL(hp)           :: MW_CH4, MW_Iso(nIsotopologue)
    REAL(hp)           :: R13, RD, F12, F13, FH, FD
    REAL(hp)           :: D13C(HcoState%NX,HcoState%NY)
    REAL(hp)           :: DD(HcoState%NX,HcoState%NY)
    REAL(hp)           :: BulkFlux(HcoState%NX,HcoState%NY,HcoState%NZ)
    REAL(hp)           :: AssembledFlux(HcoState%NX,HcoState%NY,             &
                                        HcoState%NZ)
    REAL(hp)           :: IsoFlux(HcoState%NX,HcoState%NY,HcoState%NZ,      &
                                  nIsotopologue)

    RC      = HCO_SUCCESS
    ThisLoc = ' -> at HCOX_CH4_Isotopologues_Run'
    D13C_Lct => NULL()
    DD_Lct   => NULL()

    IF ( ExtState%CH4_Isotopologues <= 0 ) RETURN

    idCH4 = HCO_GetHcoID( 'CH4', HcoState )
    IF ( idCH4 < 1 ) THEN
       ErrMsg = 'CH4 isotope extension requires bulk CH4'
       CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    idIso  = -1
    MW_Iso = 0.0_hp
    DO N = 1, nIsotopologue
       idIso(N) = HCO_GetHcoID( TRIM(Isotopologue(N)), HcoState )
       IF ( idIso(N) > 0 ) MW_Iso(N) = HcoState%Spc(idIso(N))%MW_g
    ENDDO

    ! The extension is intentionally a no-op if no isotopologue is active.
    IF ( ALL(idIso < 1) ) RETURN

    MW_CH4 = HcoState%Spc(idCH4)%MW_g
    IF ( MW_CH4 <= 0.0_hp .OR.                                     &
         ANY(MW_Iso <= 0.0_hp .AND. idIso > 0) ) THEN
       ErrMsg = 'Invalid methane molecular weight in species database'
       CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    CALL HcoClock_Get( HcoState%Clock, cYYYY=YYYY, cMM=MM, cDD=DDnow,  &
                       cH=HH, cM=MN, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Could not read the HEMCO clock'
       CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    AssembledFlux = 0.0_hp
    SignatureLoaded = .FALSE.
    DO S = 1, nSector + 1

       IF ( S <= nSector ) THEN
          SigIdx  = S
          ThisCat = Category(S)
          BulkFlux = HcoState%CategoryFlux(:,:,:,S)
          AssembledFlux = AssembledFlux + BulkFlux
       ELSE
          ! CH4 created by preceding extensions is not part of the base
          ! emission list.  The only such source in this configuration is
          ! biomass burning from GFED or GFAS, so partition the residual with
          ! the biomass-burning signature.
          SigIdx  = 9
          ThisCat = Category(SigIdx)
          BulkFlux = MAX( 0.0_hp,                                       &
                          HcoState%Spc(idCH4)%Emis%Val - AssembledFlux )
       ENDIF

       IF ( .NOT. SignatureLoaded(SigIdx) ) THEN
          IF ( UpdateInterval(SigIdx) == 0 ) THEN
             FieldName = 'CH4ISO_D13C_' // TRIM(Sector(SigIdx))
             Found = .FALSE.
             D13C_Lct => NULL()
             CALL ListCont_Find( HcoState%EmisList, TRIM(FieldName),       &
                                 Found, D13C_Lct )
             IF ( .NOT. Found ) THEN
                ErrMsg = 'Missing CH4 isotope signature field: ' //       &
                         TRIM(FieldName)
                CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
                EXIT
             ENDIF
             IF ( .NOT. ASSOCIATED(D13C_Lct%Dct) ) THEN
                ErrMsg = 'Invalid CH4 isotope signature field: ' //       &
                         TRIM(FieldName)
                CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
                EXIT
             ENDIF
             IF ( .NOT. ASSOCIATED(D13C_Lct%Dct%Dta) ) THEN
                ErrMsg = 'Invalid CH4 isotope signature data: ' //        &
                         TRIM(FieldName)
                CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
                EXIT
             ENDIF

             FieldName = 'CH4ISO_DD_' // TRIM(Sector(SigIdx))
             Found = .FALSE.
             DD_Lct => NULL()
             CALL ListCont_Find( HcoState%EmisList, TRIM(FieldName),       &
                                 Found, DD_Lct )
             IF ( .NOT. Found ) THEN
                ErrMsg = 'Missing CH4 isotope signature field: ' //       &
                         TRIM(FieldName)
                CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
                EXIT
             ENDIF
             IF ( .NOT. ASSOCIATED(DD_Lct%Dct) ) THEN
                ErrMsg = 'Invalid CH4 isotope signature field: ' //       &
                         TRIM(FieldName)
                CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
                EXIT
             ENDIF
             IF ( .NOT. ASSOCIATED(DD_Lct%Dct%Dta) ) THEN
                ErrMsg = 'Invalid CH4 isotope signature data: ' //        &
                         TRIM(FieldName)
                CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
                EXIT
             ENDIF

             UpdateInterval(SigIdx) = MAX(                              &
                  SignatureUpdateInterval(D13C_Lct%Dct%Dta),             &
                  SignatureUpdateInterval(DD_Lct%Dct%Dta) )
          ENDIF

          RefreshCache = SignatureCheckIsDue( SigIdx, YYYY, MM, DDnow,    &
                                              HH, MN )
          IF ( .NOT. RefreshCache ) THEN
             SignatureLoaded(SigIdx) = .TRUE.
          ELSE

             FieldName = 'CH4ISO_D13C_' // TRIM(Sector(SigIdx))
             Found = .FALSE.
             CALL HCO_EvalFld( HcoState, TRIM(FieldName), D13C, RC,       &
                               Found=Found )
             IF ( RC /= HCO_SUCCESS .OR. .NOT. Found ) THEN
                ErrMsg = 'Missing CH4 isotope signature field: ' //      &
                         TRIM(FieldName)
                CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
                EXIT
             ENDIF

             FieldName = 'CH4ISO_DD_' // TRIM(Sector(SigIdx))
             Found = .FALSE.
             CALL HCO_EvalFld( HcoState, TRIM(FieldName), DD, RC,         &
                               Found=Found )
             IF ( RC /= HCO_SUCCESS .OR. .NOT. Found ) THEN
                ErrMsg = 'Missing CH4 isotope signature field: ' //      &
                         TRIM(FieldName)
                CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
                EXIT
             ENDIF

             RefreshCache = .NOT. CacheValid(SigIdx)
             IF ( CacheValid(SigIdx) ) THEN
                RefreshCache =                                           &
                     ANY( D13C /= D13C_Cache(:,:,SigIdx) ) .OR.          &
                     ANY( DD   /= DD_Cache(:,:,SigIdx)   )
             ENDIF

             IF ( RefreshCache ) THEN
                DO J = 1, HcoState%NY
                DO I = 1, HcoState%NX

                   R13 = R13_VPDB *                                      &
                         ( 1.0_hp + ( D13C(I,J) / 1000.0_hp ) )
                   RD  = RD_VSMOW *                                      &
                         ( 1.0_hp + ( DD(I,J) / 1000.0_hp ) )

                   IF ( R13 < 0.0_hp .OR. RD < 0.0_hp .OR.               &
                        RD >= ( 1.0_hp / 3.0_hp ) ) THEN
                      WRITE(ErrMsg,*) 'Invalid isotope signature for ',   &
                           'sector ', TRIM(Sector(SigIdx)),                &
                           ' at grid box ', I, J
                      CALL HCO_ERROR( TRIM(ErrMsg), RC, ThisLoc )
                      EXIT
                   ENDIF

                   F13 = R13 / ( 1.0_hp + R13 )
                   F12 = 1.0_hp - F13
                   FD  = RD  / ( 1.0_hp + RD  )
                   FH  = 1.0_hp - FD

                   ! Retain only methane containing zero or one D atom.
                   ! The sum of these four fractions is slightly less
                   ! than one; the residual is multiply deuterated CH4.
                   Fraction_Cache(I,J,1,SigIdx) =       F12 * FH**4
                   Fraction_Cache(I,J,2,SigIdx) =       F13 * FH**4
                   Fraction_Cache(I,J,3,SigIdx) =                      &
                        4.0_hp * F12 * FH**3 * FD
                   Fraction_Cache(I,J,4,SigIdx) =                      &
                        4.0_hp * F13 * FH**3 * FD
                ENDDO
                IF ( RC /= HCO_SUCCESS ) EXIT
                ENDDO
                IF ( RC /= HCO_SUCCESS ) EXIT

                D13C_Cache(:,:,SigIdx) = D13C
                DD_Cache(:,:,SigIdx)   = DD
                CacheValid(SigIdx)     = .TRUE.
             ENDIF
             LastYear(SigIdx)   = YYYY
             LastMonth(SigIdx)  = MM
             LastDay(SigIdx)    = DDnow
             LastHour(SigIdx)   = HH
             LastMinute(SigIdx) = MN
          ENDIF
          SignatureLoaded(SigIdx) = .TRUE.
       ENDIF

       DO N = 1, nIsotopologue
          IF ( idIso(N) < 1 ) CYCLE
          DO L = 1, HcoState%NZ
          DO J = 1, HcoState%NY
          DO I = 1, HcoState%NX
             IsoFlux(I,J,L,N) = BulkFlux(I,J,L) *                       &
                                Fraction_Cache(I,J,N,SigIdx) *           &
                                MW_Iso(N) / MW_CH4
          ENDDO
          ENDDO
          ENDDO
       ENDDO

       DO N = 1, nIsotopologue
          IF ( idIso(N) < 1 ) CYCLE
          CALL HCO_EmisAdd( HcoState, IsoFlux(:,:,:,N), idIso(N), RC,      &
                            ExtNr=ExtState%CH4_Isotopologues,               &
                            Cat=ThisCat )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Could not add emissions for ' //                    &
                      TRIM(Isotopologue(N))
             CALL HCO_ERROR( ErrMsg, RC, ThisLoc )
             EXIT
          ENDIF
       ENDDO
       IF ( RC /= HCO_SUCCESS ) EXIT
    ENDDO

  END SUBROUTINE HCOX_CH4_Isotopologues_Run
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SignatureUpdateInterval
!
! !DESCRIPTION: Returns how often a HEMCO signature field can change.
!\
!\
! !INTERFACE:
!
  FUNCTION SignatureUpdateInterval( Dta ) RESULT( Interval )
!
! !INPUT PARAMETERS:
!
    TYPE(FileData), POINTER :: Dta
!
! !RETURN VALUE:
!
    INTEGER :: Interval
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( Dta%UpdtFlag == HCO_UFLAG_ALWAYS ) THEN
       Interval = 7
    ELSEIF ( Dta%UpdtFlag == HCO_UFLAG_3HR ) THEN
       Interval = 6
    ELSEIF ( Dta%ncHrs(1) /= Dta%ncHrs(2) ) THEN
       Interval = 5
    ELSEIF ( Dta%ncDys(1) /= Dta%ncDys(2) ) THEN
       Interval = 4
    ELSEIF ( Dta%ncMts(1) /= Dta%ncMts(2) ) THEN
       Interval = 3
    ELSEIF ( Dta%ncYrs(1) /= Dta%ncYrs(2) ) THEN
       Interval = 2
    ELSE
       Interval = 1
    ENDIF

  END FUNCTION SignatureUpdateInterval
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SignatureCheckIsDue
!
! !DESCRIPTION: Returns true when a signature field may have advanced to a
!  new HEMCO time slice.  Values are compared before fractions are rebuilt.
!\
!\
! !INTERFACE:
!
  FUNCTION SignatureCheckIsDue( S, YYYY, MM, DD, HH, MN ) RESULT( IsDue )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: S, YYYY, MM, DD, HH, MN
!
! !RETURN VALUE:
!
    LOGICAL :: IsDue
!EOP
!------------------------------------------------------------------------------
!BOC

    IsDue = .NOT. CacheValid(S)
    IF ( IsDue ) RETURN

    SELECT CASE ( UpdateInterval(S) )
       CASE ( 1 )
          IsDue = .FALSE.
       CASE ( 2 )
          IsDue = YYYY /= LastYear(S)
       CASE ( 3 )
          IsDue = YYYY /= LastYear(S) .OR. MM /= LastMonth(S)
       CASE ( 4 )
          IsDue = YYYY /= LastYear(S)  .OR. MM /= LastMonth(S) .OR.       &
                  DD   /= LastDay(S)
       CASE ( 5 )
          IsDue = YYYY /= LastYear(S)  .OR. MM /= LastMonth(S) .OR.       &
                  DD   /= LastDay(S)   .OR. HH /= LastHour(S)
       CASE ( 6 )
          IsDue = YYYY /= LastYear(S)  .OR. MM /= LastMonth(S) .OR.       &
                  DD   /= LastDay(S)   .OR.                              &
                  HH / 3 /= LastHour(S) / 3
       CASE DEFAULT
          IsDue = YYYY /= LastYear(S)  .OR. MM /= LastMonth(S) .OR.       &
                  DD   /= LastDay(S)   .OR. HH /= LastHour(S)  .OR.       &
                  MN   /= LastMinute(S)
    END SELECT

  END FUNCTION SignatureCheckIsDue
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_CH4_Isotopologues_Final
!
! !DESCRIPTION: Finalizes the stable methane isotopologue extension.
!\
!\
! !INTERFACE:
!
  SUBROUTINE HCOX_CH4_Isotopologues_Final( ExtState )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER :: ExtState
!EOP
!------------------------------------------------------------------------------
!BOC

    ExtState%CH4_Isotopologues = -1
    IF ( ALLOCATED(D13C_Cache)     ) DEALLOCATE(D13C_Cache)
    IF ( ALLOCATED(DD_Cache)       ) DEALLOCATE(DD_Cache)
    IF ( ALLOCATED(Fraction_Cache) ) DEALLOCATE(Fraction_Cache)
    CacheValid = .FALSE.
    UpdateInterval = 0
    LastYear   = -1
    LastMonth  = -1
    LastDay    = -1
    LastHour   = -1
    LastMinute = -1

  END SUBROUTINE HCOX_CH4_Isotopologues_Final
!EOC

END MODULE HCOX_CH4_Isotopologues_Mod
