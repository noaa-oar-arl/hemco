!BOC
#if defined ( ESMF_ ) && !defined( HEMCO_STANDALONE )
! The 'standard' HEMCO I/O module is used for:
! - Pure ESMF applications without MAPL dependencies
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoio_read_esmf_mod.F90
!
! !DESCRIPTION: Module HCOIO\_Read\_ESMF\_Mod is the HEMCO interface for
!  data reading within the pure ESMF environment (without MAPL dependencies).
!
!  This module implements the pure ESMF environment.
!\\
!\\
! !INTERFACE:
!
MODULE HCOIO_Read_Mod
!
! !USES:
!
  USE HCO_Types_Mod
  USE HCO_Error_Mod
  USE HCO_State_Mod,       ONLY : Hco_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOIO_Read
  PUBLIC  :: HCOIO_CloseAll
!
! !REVISION HISTORY:
!  14 Nov 2024 - B. Baker - Initial version based on pure ESMF
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_Read (ESMF version without MAPL)
!
! !DESCRIPTION: Interface routine between ESMF and HEMCO to obtain
! the data array for a given HEMCO data container. The data is obtained
! through the ESMF State interface. The HEMCO source file attribute is taken
! to identify the ESMF field name.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_Read( HcoState, Lct, RC )
!
! !USES:
!
    USE ESMF
#ifdef ESMF_8
    USE ESMF_FieldGetMod, ONLY : ESMF_FieldGet
    USE ESMF_StateGetMod, ONLY : ESMF_StateGet
#endif
    USE HCO_FILEDATA_MOD, ONLY : FileData_ArrInit

!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState
    TYPE(ListCont),   POINTER        :: Lct
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  14 Nov 2024 - B. Baker - Initial version based on pure ESMF
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                    :: II, JJ, LL, TT
    INTEGER                    :: I, J, L, T
    INTEGER                    :: STAT, lstat
    REAL,             POINTER  :: Ptr3D(:,:,:)
    REAL,             POINTER  :: Ptr2D(:,:)
    TYPE(ESMF_State), POINTER  :: IMPORT
    TYPE(ESMF_Field)           :: Field
    CHARACTER(LEN=255)         :: MSG
    CHARACTER(LEN=255), PARAMETER :: LOC = 'HCOIO_READ (hcoio_read_esmf_mod.F90)'
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam

    !=================================================================
    ! HCOIO_READ begins here
    !=================================================================

    ! For error handling
    Iam = LOC
    CALL HCO_ENTER( HcoState%Config%Err,  LOC, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
        CALL HCO_ERROR( 'ERROR 0', RC, THISLOC=LOC )
        RETURN
    ENDIF

    ! Point to ESMF IMPORT object
    IMPORT => HcoState%IMPORT
    IF (.NOT. ASSOCIATED(IMPORT)) THEN
        CALL HCO_ERROR('HcoState%IMPORT not associated', RC, THISLOC=LOC)
        RETURN
    ENDIF

    ! Init pointers
    Ptr3D => NULL()
    Ptr2D => NULL()

    ! Verbose?
    IF ( HcoState%Config%doVerbose ) THEN
       MSG = 'Reading from ESMF State: ' // TRIM(Lct%Dct%Dta%ncFile)
       CALL HCO_MSG(MSG,LUN=HcoState%Config%hcoLogLUN)
    ENDIF

    !-----------------------------------------------------------------
    ! Read 3D data from ESMF
    !-----------------------------------------------------------------
    IF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN
       ! Get field from import state using the canonical field name
       CALL ESMF_StateGet(IMPORT, itemName=TRIM(Lct%Dct%Dta%ncFile), field=Field, rc=lstat)

       IF (lstat == ESMF_SUCCESS) THEN
          ! Get the local array from the field
          CALL ESMF_FieldGet(Field, localDe=0, farrayPtr=Ptr3D, rc=lstat)

          IF (lstat == ESMF_SUCCESS .AND. ASSOCIATED(Ptr3D)) THEN
             ! Get array dimensions
             II = SIZE(Ptr3D,1)
             JJ = SIZE(Ptr3D,2)
             LL = SIZE(Ptr3D,3)
             TT = 1

             ! Define HEMCO array if not yet defined.
             IF ( .NOT. ASSOCIATED(Lct%Dct%Dta%V3) ) THEN
                ! Use pointer if types match
                CALL FileData_ArrInit( Lct%Dct%Dta, TT, 0, 0, 0, RC )
                IF ( RC /= HCO_SUCCESS ) THEN
                    CALL HCO_ERROR( 'ERROR 1', RC, THISLOC=LOC )
                    RETURN
                ENDIF
             ENDIF

             ! Pointer to data. HEMCO expects data to have surface level at
             ! index 1 ('up').
             Lct%Dct%Dta%V3(1)%Val => Ptr3D(:,:,LL:1:-1)

             ! Verbose
             IF ( HcoState%Config%doVerbose .AND. HcoState%amIRoot ) THEN
                MSG = 'HEMCO: array pointer vertically flipped relative to ESMF Import ' // TRIM(Lct%Dct%Dta%ncFile)
                CALL HCO_MSG(MSG)
             ENDIF
          ELSE
             MSG = 'Cannot get 3D pointer: ' // TRIM(Lct%Dct%Dta%ncFile)
             CALL HCO_ERROR( MSG, RC )
             RETURN
          ENDIF
       ELSE
          MSG = 'Cannot find 3D field: ' // TRIM(Lct%Dct%Dta%ncFile)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

    !-----------------------------------------------------------------
    ! Read 2D data from ESMF
    !-----------------------------------------------------------------
    ELSEIF ( Lct%Dct%Dta%SpaceDim == 2 ) THEN

       ! Get field from import state using the canonical field name
       CALL ESMF_StateGet(IMPORT, itemName=TRIM(Lct%Dct%Dta%ncFile), field=Field, rc=lstat)

       IF (lstat == ESMF_SUCCESS) THEN
          ! Get the local array from the field
          CALL ESMF_FieldGet(Field, localDe=0, farrayPtr=Ptr2D, rc=lstat)

          IF (lstat == ESMF_SUCCESS .AND. ASSOCIATED(Ptr2D)) THEN
             ! Get array dimensions
             II = SIZE(Ptr2D,1)
             JJ = SIZE(Ptr2D,2)
             LL = 1
             TT = 1

             ! Define HEMCO array pointer if not yet defined
             IF ( .NOT. ASSOCIATED(Lct%Dct%Dta%V2) ) THEN
                CALL FileData_ArrInit( Lct%Dct%Dta, TT, 0, 0, RC )
                IF ( RC /= HCO_SUCCESS ) THEN
                    CALL HCO_ERROR( 'ERROR 2', RC, THISLOC=LOC )
                    RETURN
                ENDIF
             ENDIF

             ! Pointer to data
             Lct%Dct%Dta%V2(1)%Val => Ptr2D
          ELSE
             MSG = 'Cannot get 2D pointer: ' // TRIM(Lct%Dct%Dta%ncFile)
             CALL HCO_ERROR( MSG, RC )
             RETURN
          ENDIF
       ELSE
          MSG = 'Cannot find 2D field: ' // TRIM(Lct%Dct%Dta%ncFile)
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF

    ENDIF

    !-----------------------------------------------------------------
    ! Cleanup and leave
    !-----------------------------------------------------------------
    Ptr3D  => NULL()
    Ptr2D  => NULL()
    IMPORT => NULL()

    ! Return w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err,  RC )

  END SUBROUTINE HCOIO_Read
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_CloseAll
!
! !DESCRIPTION: Subroutine HCOIO\_CloseAll makes sure that there is no open
! netCDF file left in the stream. This is a stub as there is no such handling
! within HEMCO for the ESMF environment, it is performed by ESMF.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_CloseAll( HcoState, RC )
!
! !INPUT PARAMTERS:
!
    TYPE(HCO_State), POINTER          :: HcoState    ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  14 Nov 2024 - B. Baker - Initial version based on pure ESMF
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !======================================================================
    ! HCOIO_CloseAll begins here
    !======================================================================

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_CloseAll
!EOC
END MODULE HCOIO_Read_Mod
#endif