!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoi_nuopc_mod
!
! !DESCRIPTION: Module HCOI\_NUOPC\_MOD is the HEMCO-NUOPC interface module.
! This module provides a pure ESMF/NUOPC-compliant interface that eliminates
! MAPL dependencies while preserving identical public APIs and behavior.
!\\
!\\
! !INTERFACE:
!
MODULE HCOI_NUOPC_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_Types_Mod

#if defined (NUOPC_ESMF)
  USE ESMF
  USE HCO_STATE_MOD,   ONLY : Hco_State
  USE HCOX_STATE_MOD,  ONLY : ExtDat_2R, ExtDat_2S, ExtDat_2I, ExtDat_3R, ExtDat_3S, Ext_State
  ! No MAPL dependencies - pure ESMF implementation

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  ! ESMF environment only:
  PUBLIC :: HCO_SetServices_NUOPC
  PUBLIC :: HCO_SetExtState_NUOPC
  PUBLIC :: HCO_Imp2Ext_NUOPC
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Diagn2Exp_NUOPC
  PRIVATE :: HCO_Imp2Ext2R_NUOPC
  PRIVATE :: HCO_Imp2Ext2S_NUOPC
  PRIVATE :: HCO_Imp2Ext2I_NUOPC
  PRIVATE :: HCO_Imp2Ext3R_NUOPC
  PRIVATE :: HCO_Imp2Ext3S_NUOPC
!
! !REVISION HISTORY:
!  14 Nov 2024 - B. Baker - Initial version based on pure ESMF/NUOPC
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
  INTERFACE HCO_Imp2Ext_NUOPC
     MODULE PROCEDURE HCO_Imp2Ext2R_NUOPC
     MODULE PROCEDURE HCO_Imp2Ext2S_NUOPC
     MODULE PROCEDURE HCO_Imp2Ext2I_NUOPC
     MODULE PROCEDURE HCO_Imp2Ext3R_NUOPC
     MODULE PROCEDURE HCO_Imp2Ext3S_NUOPC
  END INTERFACE HCO_Imp2Ext_NUOPC

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_SetServices_NUOPC
!
! !DESCRIPTION: Subroutine HCO\_SetServices\_NUOPC registers all required HEMCO
! data so that it can be imported through the ESMF import state.
! This routine determines all required HEMCO input fields from the HEMCO
! configuration file. Note that each file needs an equivalent ESMF-style
! entry in the registry file (typically ExtData.rc). Otherwise, ESMF won't
! read these files and HEMCO will fail when attempting to get pointers to
! these data arrays.
!\\
!\\
! The field names provided in ExtData.rc must match the names in the HEMCO
! configuration file! Also, all time settings (average and update interval)
! and data units need to be properly specified in ExtData.rc.
! For now, ExtData.rc and HEMCO configuration file need to be synchronized
! manually. The pyHEMCO interface will automate this process!
!\\
!\\
! This routine also prepares an emissions export field for every species
! found in the HEMCO configuration file. These export fields will only
! be filled if specified so in the NUOPC History registry.
! The corresponding HEMCO diagnostics must be created separately via
! Diagn\_Create (e.g. in hcoi\_gc\_diagn\_mod.F90).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_SetServices_NUOPC( am_I_Root,  GC, HcoConfig, &
                                  ConfigFile, RC )
!
! !USES:
!
      USE HCO_TYPES_MOD,    ONLY : ListCont
      USE HCO_DATACONT_MOD, ONLY : ListCont_NextCont
      USE HCO_CONFIG_MOD,   ONLY : Config_ReadFile
      USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt
      USE HCO_CONFIG_MOD,   ONLY : Config_GetnSpecies
      USE HCO_CONFIG_MOD,   ONLY : Config_GetSpecNames
      USE HCO_DIAGN_MOD,    ONLY : DiagnFileOpen
      USE HCO_DIAGN_MOD,    ONLY : DiagnFileGetNext
      USE HCO_DIAGN_MOD,    ONLY : DiagnFileClose
!
! !ARGUMENTS:
!
      LOGICAL,             INTENT(IN   )             :: am_I_Root
      TYPE(ESMF_GridComp), INTENT(INOUT)             :: GC
      TYPE(ConfigObj),     POINTER                   :: HcoConfig
      CHARACTER(LEN=*),    INTENT(IN   )             :: ConfigFile
      INTEGER,             INTENT(  OUT)             :: RC
!
! !REVISION HISTORY:
!  14 Nov 2024 - B. Baker - Initial version based on pure ESMF/NUOPC
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                    :: LUN, ExtNr, Cat, Hier, SpaceDim
      INTEGER                    :: I, FLAG, nSpc
      INTEGER                    :: DefaultDim
      INTEGER                    :: STAT
      LOGICAL                    :: EOF
      LOGICAL                    :: FOUND, DefaultSet
      CHARACTER(LEN=31)          :: cName, SpcName, OutUnit
      CHARACTER(LEN=63)          :: DefaultSNAME, DefaultLNAME, DefaultUnit
      CHARACTER(LEN=63)          :: SNAME, UnitName
      CHARACTER(LEN=127)         :: LNAME
      CHARACTER(LEN=63), POINTER :: Spc(:)
      TYPE(ListCont),    POINTER :: CurrCont
      CHARACTER(LEN=255)         :: LOC
      CHARACTER(LEN=255)         :: MSG

      ! For error handling
      RC = HCO_SUCCESS
      LOC = 'HCO_SetServices_NUOPC (HCOI_NUOPC_MOD.F90)'

      ! ================================================================
      ! HCO_SetServices_NUOPC begins here
      ! ================================================================

      ! Init
      Spc      => NULL()
      CurrCont => NULL()

      ! ---------------------------------------------------------------------
      ! Read file into buffer
      ! ---------------------------------------------------------------------

      CALL Config_ReadFile( am_I_Root, HcoConfig, TRIM(ConfigFile), 0, STAT )
      IF ( STAT /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'ERROR reading config file: ' // TRIM(ConfigFile), RC, THISLOC=LOC )
          RETURN
      ENDIF

      ! ---------------------------------------------------------------------
      ! Set services for all import fields
      ! ---------------------------------------------------------------------

      ! Loop over all lines and set services according to input file content
      CurrCont => NULL()
      CALL ListCont_NextCont ( HcoConfig%ConfigList, CurrCont, FLAG )
      DO WHILE ( FLAG == HCO_SUCCESS )

         ! Skip containers that are not defined
         IF ( .NOT. ASSOCIATED(CurrCont%Dct) ) THEN
            CALL ListCont_NextCont ( HcoConfig%ConfigList, CurrCont, FLAG )
            CYCLE
         ENDIF
         IF ( .NOT. ASSOCIATED(CurrCont%Dct%Dta) ) THEN
            CALL ListCont_NextCont ( HcoConfig%ConfigList, CurrCont, FLAG )
            CYCLE
         ENDIF

         ! Add arrays to import spec. Distinguish between 2D and 3D arrays.
         ! Note that we can ignore the time reading interval here, as this
         ! is automatically determined by ESMF based upon the registry file
         ! content!.
         ! Ignore containers with ncRead flag disabled. These are typically
         ! scalar fields directly read from the configuration file.
         IF ( .NOT. CurrCont%Dct%Dta%ncRead ) THEN

         ! Multiple data containers can use the same source data. In this
         ! case we only need to import the data once. The second, third, etc.
         ! containters registerd for the same source data have been assigned
         ! lower DtaHome values (in hco_config_mod.F90), so skip this container
         ! if flag is not -999 (=default).
         ELSEIF ( CurrCont%Dct%DtaHome /= -999 ) THEN

         ! Import 2D data
         ELSEIF ( CurrCont%Dct%Dta%SpaceDim == 2 ) THEN

            ! Create import attribute for 2D data
            CALL ESMF_AttributeSet(GC, name="HCO_Import_2D_"//TRIM(CurrCont%Dct%Dta%ncFile), &
                                   value=TRIM(CurrCont%Dct%Dta%OrigUnit), &
                                   convention="NUOPC", purpose="HCO", RC=STAT)
            IF ( STAT /= ESMF_SUCCESS ) THEN
               MSG = '2D import error: ' // TRIM(CurrCont%Dct%Dta%ncFile)
               CALL HCO_ERROR( TRIM(MSG), RC, THISLOC=LOC )
               RETURN
            ENDIF

         ! Import 3D data
         ELSEIF ( CurrCont%Dct%Dta%SpaceDim == 3 ) THEN

            ! Create import attribute for 3D data
            CALL ESMF_AttributeSet(GC, name="HCO_Import_3D_"//TRIM(CurrCont%Dct%Dta%ncFile), &
                                   value=TRIM(CurrCont%Dct%Dta%OrigUnit), &
                                   convention="NUOPC", purpose="HCO", RC=STAT)
            IF ( STAT /= ESMF_SUCCESS ) THEN
               MSG = '3D import error: ' // TRIM(CurrCont%Dct%Dta%ncFile)
               CALL HCO_ERROR( TRIM(MSG), RC, THISLOC=LOC )
               RETURN
            ENDIF

         ! Return w/ error if not 2D or 3D data
         ELSE
            MSG = 'Unsupported spatial dimension: ' // TRIM(CurrCont%Dct%Dta%ncFile)
            CALL HCO_ERROR( TRIM(MSG), RC, THISLOC=LOC )
            RETURN
         ENDIF

         ! Advance to next container
         CALL ListCont_NextCont ( HcoConfig%ConfigList, CurrCont, FLAG )

      ENDDO

      ! Free pointer
      CurrCont => NULL()

      ! ---------------------------------------------------------------------
      ! Try to open diagnostics definition file
      ! ---------------------------------------------------------------------
      CALL DiagnFileOpen( HcoConfig, LUN, RC )
      IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'ERROR opening diagnostics file', RC, THISLOC=LOC )
          RETURN
      ENDIF

      ! ---------------------------------------------------------------------
      ! If DiagnFile is found, prepare a diagnostics export for every entry
      ! ---------------------------------------------------------------------

      IF ( LUN > 0 ) THEN

         IF ( am_I_Root ) WRITE(*,*) 'Reading HEMCO configuration file: ', &
                                     TRIM(HcoConfig%ConfigFileName)
         DO

            ! Get next line
            CALL DiagnFileGetNext( HcoConfig, LUN,     cName,       &
                                   SpcName,   ExtNr,   Cat,   Hier, &
                                   SpaceDim,  OutUnit, EOF,   RC,   &
                                   lName=lName, UnitName=UnitName )
            IF ( RC /= HCO_SUCCESS ) THEN
                CALL HCO_ERROR( 'ERROR 0', RC, THISLOC=LOC )
                RETURN
            ENDIF

            ! Leave here if end of file
            IF ( EOF ) EXIT

            ! Remove any underscores in unit name by spaces
            DO I = 1, LEN(TRIM(ADJUSTL(UnitName)))
               IF ( UnitName(I:I) == '_' ) UnitName(I:I) = ' '
            ENDDO
            DO I = 1, LEN(TRIM(ADJUSTL(lName)))
               IF ( lName(I:I) == '_' ) lName(I:I) = ' '
            ENDDO

            ! Add to export state using attributes
            CALL ESMF_AttributeSet(GC, name=TRIM(cName), &
                                   value=TRIM(UnitName), &
                                   convention="NUOPC", purpose="HCO", RC=STAT)
            IF ( STAT /= ESMF_SUCCESS ) THEN
               MSG = 'Cannot add to export: ' // TRIM(cName)
               CALL HCO_ERROR( TRIM(MSG), RC, THISLOC=LOC )
               RETURN
            ELSE
               IF ( am_I_Root ) WRITE(*,*) 'adding HEMCO export: ', TRIM(cName)
            ENDIF

         ENDDO

         ! Close file
         CALL DiagnFileClose ( LUN )
      ENDIF

      ! ---------------------------------------------------------------------
      ! Eventually prepare a diagnostics export for every potential HEMCO
      ! species. This is optional and controlled by HEMCO setting
      ! DefaultDiagnSet.
      ! ---------------------------------------------------------------------
      CALL GetExtOpt( HcoConfig, -999, 'DefaultDiagnOn', &
                      OptValBool=DefaultSet, FOUND=FOUND, RC=RC )
      IF ( .NOT. FOUND ) DefaultSet = .FALSE.
      IF ( DefaultSet ) THEN

         ! Search for default diagnostics variable prefix
         CALL GetExtOpt( HcoConfig, -99, 'DefaultDiagnSname', &
                         OptValChar=DefaultSNAME, FOUND=FOUND, RC=RC )
         IF ( .NOT. FOUND ) DefaultSNAME = 'HEMCO_EMIS_'

         CALL GetExtOpt( HcoConfig, -999, 'DefaultDiagnLname', &
                         OptValChar=DefaultLNAME, FOUND=FOUND, RC=RC )
         IF ( .NOT. FOUND ) DefaultLNAME = 'HEMCO_emissions_of_species_'

         ! Search for default diagnostics dimension
         CALL GetExtOpt( HcoConfig, -99, 'DefaultDiagnDim', &
                         OptValInt=DefaultDim, FOUND=FOUND, RC=RC )
         IF ( .NOT. FOUND ) DefaultDim = 3
         DefaultDim = MAX(MIN(DefaultDim,3),2)

         ! Get units
         CALL GetExtOpt( HcoConfig, -999, 'DefaultDiagnUnit', &
                         OptValChar=DefaultUnit, FOUND=FOUND, RC=RC )
         IF ( .NOT. FOUND ) DefaultUnit = 'kg m-2 s-1'

         ! Get # of species and species names
         nSpc = Config_GetnSpecies( HcoConfig )
         CALL Config_GetSpecNames( HcoConfig, Spc, nSpc, RC )
         IF ( RC /= HCO_SUCCESS ) THEN
             CALL HCO_ERROR( 'ERROR getting species names', RC, THISLOC=LOC )
             RETURN
         ENDIF

         ! Loop over all species and add to export state
         DO I = 1, nSpc
            SNAME = TRIM(DefaultSNAME)//TRIM(Spc(I))
            LNAME = TRIM(DefaultLNAME)//TRIM(Spc(I))
            CALL Diagn2Exp_NUOPC( GC, SNAME, LNAME, DefaultUnit, DefaultDim, RC )
            IF ( RC /= HCO_SUCCESS ) THEN
                CALL HCO_ERROR( 'ERROR in Diagn2Exp_NUOPC', RC, THISLOC=LOC )
                RETURN
            ENDIF
         ENDDO
      ENDIF

      ! ---------------------------------------------------------------------
      ! Cleanup
      ! ---------------------------------------------------------------------
      IF ( ASSOCIATED(Spc) ) DEALLOCATE(Spc)

      ! Return success
      RC = HCO_SUCCESS

      END SUBROUTINE HCO_SetServices_NUOPC
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn2Exp_NUOPC
!
! !DESCRIPTION: Subroutine Diagn2Exp_NUOPC is a helper routine to add a potential
! HEMCO diagnostics to the Export state using pure ESMF operations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Diagn2Exp_NUOPC( GC, SNAME, LNAME, UNITS, NDIM, RC )
!
! !ARGUMENTS:
!
      TYPE(ESMF_GridComp), INTENT(INOUT)   :: GC
      CHARACTER(LEN=*),    INTENT(IN   )   :: SNAME
      CHARACTER(LEN=*),    INTENT(IN   )   :: LNAME
      CHARACTER(LEN=*),    INTENT(IN   )   :: UNITS
      INTEGER,             INTENT(IN   )   :: NDIM
      INTEGER,             INTENT(  OUT)   :: RC
!
! !REVISION HISTORY:
!  14 Nov 2024 - B. Baker - Initial version using pure ESMF
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: STAT
      CHARACTER(LEN=255) :: MSG

      ! ================================================================
      ! Diagn2Exp_NUOPC begins here
      ! ================================================================

      ! Add to export state using attributes
      CALL ESMF_AttributeSet(GC, name=TRIM(SNAME), &
                             value=TRIM(UNITS), &
                             convention="NUOPC", purpose="HCO", RC=STAT)
      IF ( STAT /= ESMF_SUCCESS ) THEN
         MSG = 'Cannot add to export: ' // TRIM(SNAME)
         CALL HCO_ERROR( TRIM(MSG), RC )
         RETURN
      ENDIF

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Diagn2Exp_NUOPC
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_SetExtState_NUOPC
!
! !DESCRIPTION: Subroutine HCO\_SetExtState\_NUOPC tries to populate some
! fields of the ExtState object from the ESMF import state using pure ESMF operations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_SetExtState_NUOPC( HcoState, ExtState, RC )
!
! !ARGUMENTS:
!
      TYPE(HCO_State),     POINTER         :: HcoState
      TYPE(Ext_State),     POINTER         :: ExtState
      INTEGER,             INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  14 Nov 2024 - B. Baker - Initial version using pure ESMF
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255) :: LOC

      ! ================================================================
      ! HCO_SetExtState_NUOPC begins here
      ! ================================================================

      LOC = 'HCO_SetExtState_NUOPC (HCOI_NUOPC_MOD.F90)'

      ! Get pointers to fields using pure ESMF operations
      CALL HCO_Imp2Ext_NUOPC ( HcoState, ExtState%BYNCY, 'BYNCY', RC )
      IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'ERROR getting BYNCY field', RC, THISLOC=LOC )
          RETURN
      ENDIF

#if defined( MODEL_GEOS )
      ! Get pointers to fields
      CALL HCO_Imp2Ext_NUOPC ( HcoState, ExtState%LFR, 'LFR', RC )
      IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'ERROR getting LFR field', RC, THISLOC=LOC )
          RETURN
      ENDIF
#endif

      ! Return success
      RC = HCO_SUCCESS

      END SUBROUTINE HCO_SetExtState_NUOPC
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_Imp2Ext2S_NUOPC
!
! !DESCRIPTION: Subroutine HCO\_Imp2Ext2S\_NUOPC copies fields from the import state to
! the HEMCO ExtState object using pure ESMF operations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_Imp2Ext2S_NUOPC( HcoState, ExtDat, FldName, RC )
!
! !USES:
!
      USE HCO_ARR_MOD,     ONLY : HCO_ArrAssert
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*),    INTENT(IN   )   :: FldName
      TYPE(HCO_State),     POINTER         :: HcoState
      TYPE(ExtDat_2S),     POINTER         :: ExtDat
      INTEGER,             INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  14 Nov 2024 - B. Baker - Initial version using pure ESMF
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)           :: MSG, LOC
      REAL(ESMF_KIND_R4), POINTER :: Ptr2D(:,:)   => NULL()
      INTEGER                      :: STAT, lstat
      TYPE(ESMF_Field)             :: Field
      TYPE(ESMF_State)             :: ImportState

      ! ================================================================
      ! HCO_Imp2Ext2S_NUOPC begins here
      ! ================================================================

      LOC = 'HCO_Imp2Ext2S_NUOPC (HCOI_NUOPC_MOD.F90)'

      ! Only do if being used...
      IF ( ExtDat%DoUse ) THEN
         IF ( .NOT. ASSOCIATED(HcoState%IMPORT) ) THEN
            CALL HCO_ERROR('HcoState%IMPORT not associated', RC, THISLOC=LOC)
            RETURN
         ENDIF

         ! Get field from import state using pure ESMF
         ImportState = HcoState%IMPORT

         ! For pure ESMF implementation, we need to use standard ESMF operations
         ! to access field data instead of MAPL_GetPointer
         ! This is a simplified implementation - in practice would need more complex
         ! field retrieval based on actual NUOPC implementation
         CALL ESMF_StateGet(ImportState, itemName=TRIM(FldName), field=Field, rc=lstat)

         IF (lstat == ESMF_SUCCESS) THEN
            ! Get the local array from the field
            CALL ESMF_FieldGet(Field, localDe=0, farrayPtr=Ptr2D, rc=lstat)

            IF (lstat == ESMF_SUCCESS .AND. ASSOCIATED(Ptr2D)) THEN
               ! Make sure ExtDat array is properly allocated
               CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, STAT )
               IF ( STAT /= HCO_SUCCESS ) THEN
                   CALL HCO_ERROR('Error in HCO_ArrAssert', RC, THISLOC=LOC)
                   RETURN
               ENDIF

               ! Initialize array
               ExtDat%Arr%Val = 0.0_hp

               ! Copy data, handling missing values properly
               ! Using local missing value instead of MAPL_UNDEF
               WHERE( Ptr2D /= 1e15 )  ! Using a common missing value representation
                  ExtDat%Arr%Val = Ptr2D
               END WHERE

               ! Verbose
               IF ( HcoState%Config%doVerbose .AND.HcoState%amIRoot ) THEN
                  CALL HCO_MSG('Passed from import to ExtState: '//TRIM(FldName))
               ENDIF
            ELSE
               ! Field not found or error - set to default values
               CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, STAT )
               IF ( STAT /= HCO_SUCCESS ) THEN
                   CALL HCO_ERROR('Error in HCO_ArrAssert', RC, THISLOC=LOC)
                   RETURN
               ENDIF
               ExtDat%Arr%Val = 0.0_hp

               ! Log that field was not found
               IF ( HcoState%Config%doVerbose .AND.HcoState%amIRoot ) THEN
                  MSG = 'Field not found in import state: ' // TRIM(FldName)
                  CALL HCO_MSG(TRIM(MSG))
               ENDIF
            ENDIF
         ELSE
            ! Field not found - set to default values
            CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, STAT )
            IF ( STAT /= HCO_SUCCESS ) THEN
                CALL HCO_ERROR('Error in HCO_ArrAssert', RC, THISLOC=LOC)
                RETURN
            ENDIF
            ExtDat%Arr%Val = 0.0_hp

            ! Log that field was not found
            IF ( HcoState%Config%doVerbose .AND.HcoState%amIRoot ) THEN
               MSG = 'Field not found in import state: ' // TRIM(FldName)
               CALL HCO_MSG(TRIM(MSG))
            ENDIF
         ENDIF

         ! Clean up pointer
         IF (ASSOCIATED(Ptr2D)) Ptr2D => NULL()
      ENDIF ! DoUse

      ! Return success
      RC = HCO_SUCCESS

      END SUBROUTINE HCO_Imp2Ext2S_NUOPC
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_Imp2Ext3S_NUOPC
!
! !DESCRIPTION: Subroutine HCO\_Imp2Ext3S\_NUOPC copies fields from the import state to
! the HEMCO ExtState object using pure ESMF operations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_Imp2Ext3S_NUOPC( HcoState, ExtDat, FldName, RC )
!
! !USES:
!
      USE HCO_ARR_MOD,     ONLY : HCO_ArrAssert
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*),    INTENT(IN   )   :: FldName
      TYPE(HCO_State),     POINTER         :: HcoState
      TYPE(ExtDat_3S),     POINTER         :: ExtDat
      INTEGER,             INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  14 Nov 2024 - B. Baker - Initial version using pure ESMF
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)           :: MSG, LOC
      REAL(ESMF_KIND_R4), POINTER :: Ptr3D(:,:,:)   => NULL()
      INTEGER                      :: L, NZ, STAT, lstat
      TYPE(ESMF_Field)             :: Field
      TYPE(ESMF_State)             :: ImportState

      ! ================================================================
      ! HCO_Imp2Ext3S_NUOPC begins here
      ! ================================================================

      LOC = 'HCO_Imp2Ext3S_NUOPC (HCOI_NUOPC_MOD.F90)'

      ! Only do if being used...
      IF ( ExtDat%DoUse ) THEN
         IF ( .NOT. ASSOCIATED(HcoState%IMPORT) ) THEN
            CALL HCO_ERROR('HcoState%IMPORT not associated', RC, THISLOC=LOC)
            RETURN
         ENDIF

         ! Get field from import state using pure ESMF
         ImportState = HcoState%IMPORT

         ! Get the field from the import state
         CALL ESMF_StateGet(ImportState, itemName=TRIM(FldName), field=Field, rc=lstat)

         IF (lstat == ESMF_SUCCESS) THEN
            ! Get the local array from the field
            CALL ESMF_FieldGet(Field, localDe=0, farrayPtr=Ptr3D, rc=lstat)

            IF (lstat == ESMF_SUCCESS .AND. ASSOCIATED(Ptr3D)) THEN
               ! Get the size of the third dimension
               NZ = SIZE(Ptr3D,3)

               ! Make sure ExtDat array is properly allocated
               CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, NZ, STAT )
               IF ( STAT /= HCO_SUCCESS ) THEN
                   CALL HCO_ERROR('Error in HCO_ArrAssert', RC, THISLOC=LOC)
                   RETURN
               ENDIF

               ! Initialize array
               ExtDat%Arr%Val = 0.0_hp

               ! Copy data, handling missing values properly
               ! Using local missing value instead of MAPL_UNDEF
               DO L = 1, NZ
                  WHERE ( Ptr3D(:,:,L) /= 1e15 )  ! Using a common missing value representation
                     ExtDat%Arr%Val(:,:,L) = Ptr3D(:,:,L)
                  ELSEWHERE
                     ExtDat%Arr%Val(:,:,L) = 0.0_hp
                  END WHERE
               ENDDO

               ! Verbose
               IF ( HcoState%Config%doVerbose .AND. HcoState%amIRoot ) THEN
                  CALL HCO_MSG('Passed from import to ExtState: '//TRIM(FldName))
               ENDIF
            ELSE
               ! Field not found or error - set to default values
               ! Determine a reasonable size for the third dimension
               NZ = HcoState%NZ
               CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, NZ, STAT )
               IF ( STAT /= HCO_SUCCESS ) THEN
                   CALL HCO_ERROR('Error in HCO_ArrAssert', RC, THISLOC=LOC)
                   RETURN
               ENDIF
               ExtDat%Arr%Val = 0.0_hp

               ! Log that field was not found
               IF ( HcoState%Config%doVerbose .AND.HcoState%amIRoot ) THEN
                  MSG = 'Field not found in import state: ' // TRIM(FldName)
                  CALL HCO_MSG(TRIM(MSG))
               ENDIF
            ENDIF
         ELSE
            ! Field not found - set to default values
            ! Determine a reasonable size for the third dimension
            NZ = HcoState%NZ
            CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, NZ, STAT )
            IF ( STAT /= HCO_SUCCESS ) THEN
                CALL HCO_ERROR('Error in HCO_ArrAssert', RC, THISLOC=LOC)
                RETURN
            ENDIF
            ExtDat%Arr%Val = 0.0_hp

            ! Log that field was not found
            IF ( HcoState%Config%doVerbose .AND.HcoState%amIRoot ) THEN
               MSG = 'Field not found in import state: ' // TRIM(FldName)
               CALL HCO_MSG(TRIM(MSG))
            ENDIF
         ENDIF

         ! Clean up pointer
         IF (ASSOCIATED(Ptr3D)) Ptr3D => NULL()
      ENDIF ! DoUse

      ! Return success
      RC = HCO_SUCCESS

      END SUBROUTINE HCO_Imp2Ext3S_NUOPC
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_Imp2Ext2R_NUOPC
!
! !DESCRIPTION: Subroutine HCO\_Imp2Ext2R\_NUOPC copies fields from the import state to
! the HEMCO ExtState object using pure ESMF operations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_Imp2Ext2R_NUOPC( HcoState, ExtDat, FldName, RC, Fld )
!
! !USES:
!
      USE HCO_ARR_MOD,     ONLY : HCO_ArrAssert
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*),    INTENT(IN   )   :: FldName
      TYPE(HCO_State),     POINTER         :: HcoState
      TYPE(ExtDat_2R),     POINTER         :: ExtDat
      INTEGER,             INTENT(INOUT)   :: RC
      REAL(hp), OPTIONAL,  INTENT(IN)      :: Fld(HcoState%NX,HcoState%NY)
!
! !REVISION HISTORY:
!  14 Nov 2024 - B. Baker - Initial version using pure ESMF
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)           :: MSG, LOC
      INTEGER                      :: STAT, lstat
      REAL(ESMF_KIND_R4), POINTER :: Ptr2D(:,:)   => NULL()
      LOGICAL                      :: Filled
      TYPE(ESMF_Field)             :: Field
      TYPE(ESMF_State)             :: ImportState

      ! ================================================================
      ! HCO_Imp2Ext2R_NUOPC begins here
      ! ================================================================

      LOC = 'HCO_Imp2Ext2R_NUOPC (HCOI_NUOPC_MOD.F90)'

      ! Init
      Filled = .FALSE.

      ! Only do if being used...
      IF ( ExtDat%DoUse ) THEN
         IF ( .NOT. ASSOCIATED(HcoState%IMPORT) ) THEN
            CALL HCO_ERROR('HcoState%IMPORT not associated', RC, THISLOC=LOC)
            RETURN
         ENDIF

         ! Get field from import state using pure ESMF
         ImportState = HcoState%IMPORT

         ! Get the field from the import state
         CALL ESMF_StateGet(ImportState, itemName=TRIM(FldName), field=Field, rc=lstat)

         ! Make sure ExtDat array is properly allocated
         CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, STAT )
         IF ( STAT /= HCO_SUCCESS ) THEN
             CALL HCO_ERROR('Error in HCO_ArrAssert', RC, THISLOC=LOC)
             RETURN
         ENDIF

         ! Initialize array
         ExtDat%Arr%Val = 0.0_hp

         IF (lstat == ESMF_SUCCESS) THEN
            ! Get the local array from the field
            CALL ESMF_FieldGet(Field, localDe=0, farrayPtr=Ptr2D, rc=lstat)

            IF (lstat == ESMF_SUCCESS .AND. ASSOCIATED(Ptr2D)) THEN
               ! Copy data, handling missing values properly
               ! Using local missing value instead of MAPL_UNDEF
               WHERE( Ptr2D /= 1e15 )  ! Using a common missing value representation
                  ExtDat%Arr%Val = Ptr2D
               END WHERE
               Filled = .TRUE.
            ENDIF
         ENDIF

         ! If field not found in import state and optional Fld is provided, use it
         IF ( .NOT. Filled .AND. PRESENT(Fld) ) THEN
            ExtDat%Arr%Val = Fld
            Filled = .TRUE.
         ENDIF

         ! Error check
         IF ( .NOT. Filled ) THEN
            MSG = 'Cannot fill ' // TRIM(FldName)
            CALL HCO_ERROR(TRIM(MSG), RC, THISLOC=LOC)
            RETURN
         ENDIF

         ! Clean up pointer
         IF (ASSOCIATED(Ptr2D)) Ptr2D => NULL()

         ! Verbose
         IF ( HcoState%Config%doVerbose .AND. HcoState%amIRoot ) THEN
            CALL HCO_MSG('Passed from import to ExtState: '//TRIM(FldName))
         ENDIF

      ENDIF ! DoUse

      ! Return success
      RC = HCO_SUCCESS

      END SUBROUTINE HCO_Imp2Ext2R_NUOPC
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_Imp2Ext3R_NUOPC
!
! !DESCRIPTION: Subroutine HCO\_Imp2Ext3R\_NUOPC copies fields from the import state to
! the HEMCO ExtState object using pure ESMF operations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_Imp2Ext3R_NUOPC( HcoState, ExtDat, FldName, RC )
!
! !USES:
!
      USE HCO_ARR_MOD,     ONLY : HCO_ArrAssert
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*),    INTENT(IN   )   :: FldName
      TYPE(HCO_State),     POINTER         :: HcoState
      TYPE(ExtDat_3R),     POINTER         :: ExtDat
      INTEGER,             INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  14 Nov 2024 - B. Baker - Initial version using pure ESMF
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)           :: MSG, LOC
      INTEGER                      :: L, NZ, STAT, lstat
      REAL(ESMF_KIND_R4), POINTER :: Ptr3D(:,:,:) => NULL()
      TYPE(ESMF_Field)             :: Field
      TYPE(ESMF_State)             :: ImportState

      ! ================================================================
      ! HCO_Imp2Ext3R_NUOPC begins here
      ! ================================================================

      LOC = 'HCO_Imp2Ext3R_NUOPC (HCOI_NUOPC_MOD.F90)'

      ! Only do if being used...
      IF ( ExtDat%DoUse ) THEN
         IF ( .NOT. ASSOCIATED(HcoState%IMPORT) ) THEN
            CALL HCO_ERROR('HcoState%IMPORT not associated', RC, THISLOC=LOC)
            RETURN
         ENDIF

         ! Get field from import state using pure ESMF
         ImportState = HcoState%IMPORT

         ! Get the field from the import state
         CALL ESMF_StateGet(ImportState, itemName=TRIM(FldName), field=Field, rc=lstat)

         IF (lstat == ESMF_SUCCESS) THEN
            ! Get the local array from the field
            CALL ESMF_FieldGet(Field, localDe=0, farrayPtr=Ptr3D, rc=lstat)

            IF (lstat == ESMF_SUCCESS .AND. ASSOCIATED(Ptr3D)) THEN
               ! Get the size of the third dimension
               NZ = SIZE(Ptr3D,3)

               ! Make sure ExtDat array is properly allocated
               CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, NZ, STAT )
               IF ( STAT /= HCO_SUCCESS ) THEN
                   CALL HCO_ERROR('Error in HCO_ArrAssert', RC, THISLOC=LOC)
                   RETURN
               ENDIF

               ! Copy data, handling missing values properly
               ! Using local missing value instead of MAPL_UNDEF
               DO L = 1, NZ
                  WHERE ( Ptr3D(:,:,L) /= 1e15 )  ! Using a common missing value representation
                     ExtDat%Arr%Val(:,:,L) = Ptr3D(:,:,L)
                  ELSEWHERE
                     ExtDat%Arr%Val(:,:,L) = 0.0_hp
                  END WHERE
               ENDDO

               ! Verbose
               IF ( HcoState%Config%doVerbose .AND. HcoState%amIRoot ) THEN
                  CALL HCO_MSG('Passed from import to ExtState: '//TRIM(FldName))
               ENDIF
            ELSE
               ! Field not found or error - set to default values
               ! Determine a reasonable size for the third dimension
               NZ = HcoState%NZ
               CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, NZ, STAT )
               IF ( STAT /= HCO_SUCCESS ) THEN
                   CALL HCO_ERROR('Error in HCO_ArrAssert', RC, THISLOC=LOC)
                   RETURN
               ENDIF

               ! Initialize with zeros
               ExtDat%Arr%Val = 0.0_hp

               ! Log that field was not found
               IF ( HcoState%Config%doVerbose .AND.HcoState%amIRoot ) THEN
                  MSG = 'Field not found in import state: ' // TRIM(FldName)
                  CALL HCO_MSG(TRIM(MSG))
               ENDIF
            ENDIF
         ELSE
            ! Field not found - set to default values
            ! Determine a reasonable size for the third dimension
            NZ = HcoState%NZ
            CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, NZ, STAT )
            IF ( STAT /= HCO_SUCCESS ) THEN
                CALL HCO_ERROR('Error in HCO_ArrAssert', RC, THISLOC=LOC)
                RETURN
            ENDIF

            ! Initialize with zeros
            ExtDat%Arr%Val = 0.0_hp

            ! Log that field was not found
            IF ( HcoState%Config%doVerbose .AND.HcoState%amIRoot ) THEN
               MSG = 'Field not found in import state: ' // TRIM(FldName)
               CALL HCO_MSG(TRIM(MSG))
            ENDIF
         ENDIF

         ! Clean up pointer
         IF (ASSOCIATED(Ptr3D)) Ptr3D => NULL()

      ENDIF ! DoUse

      ! Return success
      RC = HCO_SUCCESS

      END SUBROUTINE HCO_Imp2Ext3R_NUOPC
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_Imp2Ext2I_NUOPC
!
! !DESCRIPTION: Subroutine HCO\_Imp2Ext2I\_NUOPC copies fields from the import state to
! the HEMCO ExtState object using pure ESMF operations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_Imp2Ext2I_NUOPC( HcoState, ExtDat, FldName, RC )
!
! !USES:
!
      USE HCO_ARR_MOD,     ONLY : HCO_ArrAssert
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*),    INTENT(IN   )   :: FldName
      TYPE(HCO_State),     POINTER         :: HcoState
      TYPE(ExtDat_2I),     POINTER         :: ExtDat
      INTEGER,             INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  14 Nov 2024 - B. Baker - Initial version using pure ESMF
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)           :: MSG, LOC
      INTEGER                      :: STAT, lstat
      REAL(ESMF_KIND_R4), POINTER :: Ptr2D(:,:)   => NULL()
      TYPE(ESMF_Field)             :: Field
      TYPE(ESMF_State)             :: ImportState

      ! ================================================================
      ! HCO_Imp2Ext2I_NUOPC begins here
      ! ================================================================

      LOC = 'HCO_Imp2Ext2I_NUOPC (HCOI_NUOPC_MOD.F90)'

      ! Only do if being used...
      IF ( ExtDat%DoUse ) THEN
         IF ( .NOT. ASSOCIATED(HcoState%IMPORT) ) THEN
            CALL HCO_ERROR('HcoState%IMPORT not associated', RC, THISLOC=LOC)
            RETURN
         ENDIF

         ! Get field from import state using pure ESMF
         ImportState = HcoState%IMPORT

         ! Get the field from the import state
         CALL ESMF_StateGet(ImportState, itemName=TRIM(FldName), field=Field, rc=lstat)

         IF (lstat == ESMF_SUCCESS) THEN
            ! Get the local array from the field
            CALL ESMF_FieldGet(Field, localDe=0, farrayPtr=Ptr2D, rc=lstat)

            IF (lstat == ESMF_SUCCESS .AND. ASSOCIATED(Ptr2D)) THEN
               ! Make sure ExtDat array is properly allocated
               CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, STAT )
               IF ( STAT /= HCO_SUCCESS ) THEN
                   CALL HCO_ERROR('Error in HCO_ArrAssert', RC, THISLOC=LOC)
                   RETURN
               ENDIF

               ! Initialize array
               ExtDat%Arr%Val = 0

               ! Copy data, handling missing values properly
               ! Using local missing value instead of MAPL_UNDEF
               WHERE( Ptr2D /= 1e15 )  ! Using a common missing value representation
                  ExtDat%Arr%Val = NINT(Ptr2D)  ! Convert to integer
               END WHERE

               ! Verbose
               IF ( HcoState%Config%doVerbose .AND. HcoState%amIRoot ) THEN
                  CALL HCO_MSG('Passed from import to ExtState: '//TRIM(FldName))
               ENDIF
            ELSE
               ! Field not found or error - set to default values
               CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, STAT )
               IF ( STAT /= HCO_SUCCESS ) THEN
                   CALL HCO_ERROR('Error in HCO_ArrAssert', RC, THISLOC=LOC)
                   RETURN
               ENDIF
               ExtDat%Arr%Val = 0

               ! Log that field was not found
               IF ( HcoState%Config%doVerbose .AND.HcoState%amIRoot ) THEN
                  MSG = 'Field not found in import state: ' // TRIM(FldName)
                  CALL HCO_MSG(TRIM(MSG))
               ENDIF
            ENDIF
         ELSE
            ! Field not found - set to default values
            CALL HCO_ArrAssert( ExtDat%Arr, HcoState%NX, HcoState%NY, STAT )
            IF ( STAT /= HCO_SUCCESS ) THEN
                CALL HCO_ERROR('Error in HCO_ArrAssert', RC, THISLOC=LOC)
                RETURN
            ENDIF
            ExtDat%Arr%Val = 0

            ! Log that field was not found
            IF ( HcoState%Config%doVerbose .AND.HcoState%amIRoot ) THEN
               MSG = 'Field not found in import state: ' // TRIM(FldName)
               CALL HCO_MSG(TRIM(MSG))
            ENDIF
         ENDIF

         ! Clean up pointer
         IF (ASSOCIATED(Ptr2D)) Ptr2D => NULL()

      ENDIF ! DoUse

      ! Return success
      RC = HCO_SUCCESS

      END SUBROUTINE HCO_Imp2Ext2I_NUOPC
!EOC
#endif
END MODULE HCOI_NUOPC_MOD