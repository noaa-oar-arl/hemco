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
! !MODULE: hcoio_write_esmf_mod.F90
!
! !DESCRIPTION: Module HCOIO\_Write\_ESMF\_Mod is the HEMCO output
! interface for the pure ESMF environment (without MAPL dependencies).
! In a pure ESMF environment, the HEMCO diagnostics are not directly
! written to disk but passed to the gridded component export state, where
! they can be picked up by the ESMF history component.
!\\
!\\
! !INTERFACE:
!
MODULE HCOIO_Write_Mod
!
! !USES:
!
  USE HCO_ERROR_MOD
 USE HCO_DIAGN_MOD

  IMPLICIT NONE
 PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOIO_Write
!
! !REMARKS:
!  HEMCO diagnostics are still in testing mode. We will fully activate them
!  at a later time.  They will be turned on when debugging & unit testing.
!
! !REVISION HISTORY:
!  14 Nov 2024 - B. Baker - Initial version based on pure ESMF
!  See https://github.com/geoschem/hemco for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                   Harmonized Emissions Component (HEMCO)                    !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_Write
!
! !DESCRIPTION: Subroutine HCOIO\_Write is the interface routine to
! link the HEMCO diagnostics arrays to the corresponding data pointers of the
! ESMF history component using pure ESMF operations.
!\\
!\\
! Since the history component internally organizes many diagnostics tasks such
! as output scheduling, file writing, and data averaging, all HEMCO diagnostics
! are made available to the history component on every time step, e.g. the
! entire content of the HEMCO diagnostics list is 'flushed' every time this
! subroutine is called.
!\\
!\\
! For now, all diagnostics data is copied to the corresponding ESMF data
! pointer so that this routine works for cases where the HEMCO precision is
! not equal to the ESMF precision.
!\\
!\\
! Once the HEMCO precision is pegged to the ESMF precision, we can just
! establish pointers between the export arrays and the diagnostics the first
! time this routine is called.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_WRITE ( HcoState, RC, OnlyIfFirst, COL )
!
! !USES:
!
    USE ESMF
#ifdef ESMF_8
    USE ESMF_FieldGetMod, ONLY : ESMF_FieldGet
    USE ESMF_StateGetMod, ONLY : ESMF_StateGet
#endif
    USE HCO_Types_Mod, ONLY : DiagnCont
    USE HCO_State_Mod, ONLY : HCO_State

!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState    ! HEMCO state object
    LOGICAL,          OPTIONAL, INTENT(IN   ) :: OnlyIfFirst !
    INTEGER,          OPTIONAL, INTENT(IN   ) :: COL         ! Collection Nr.
!
! !INPUT/OUTPUT PARAMETERS:
!

    INTEGER,                    INTENT(INOUT) :: RC          ! Failure or success
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
    TYPE(DiagnCont), POINTER  :: ThisDiagn
    INTEGER                   :: PS, FLAG, STAT, lstat
    CHARACTER(LEN=255)        :: MSG
    LOGICAL                   :: EOI
    REAL, POINTER             :: Ptr2D(:,:)
    REAL, POINTER             :: Ptr3D(:,:,:)
    TYPE(ESMF_Field)          :: Field
    TYPE(ESMF_State)          :: ExportState

    CHARACTER(LEN=255), PARAMETER :: LOC = 'HCOIO_WRITE_ESMF (hcoio_write_esmf_mod.F90)'

    !=================================================================
    ! HCOIO_WRITE_ESMF begins here!
    !=================================================================

    ! Assume success until otherwise
    RC  = HCO_SUCCESS

    ! Init
    ThisDiagn => NULL()
    Ptr2D     => NULL()
    Ptr3D     => NULL()

    ! Collection number
    PS = HcoState%Diagn%HcoDiagnIDDefault
    IF ( PRESENT(COL) ) PS = COL

    ! In an ESMF environment, always get all diagnostics since output
    ! is scheduled through ESMF History!
    EOI = .FALSE.

    !-----------------------------------------------------------------
    ! Connect diagnostics to export state.
    !-----------------------------------------------------------------

    ! Loop over all diagnostics in diagnostics list
    ThisDiagn => NULL()
    DO WHILE ( .TRUE. )

       ! Get next diagnostics in list. This will return the next
       ! diagnostics container that contains content to be written
       ! out on this time step.
       CALL Diagn_Get ( HcoState, EOI, ThisDiagn, FLAG, RC, COL=PS )

       IF ( RC /= HCO_SUCCESS ) THEN
           CALL HCO_ERROR( 'ERROR 0', RC, THISLOC=LOC )
           RETURN
       ENDIF
       IF ( FLAG /= HCO_SUCCESS ) EXIT

       ! Only write diagnostics if this is the first Diagn_Get call for
       ! this container and time step.
       IF ( PRESENT(OnlyIfFirst) ) THEN
          IF ( OnlyIfFirst .AND. ThisDiagn%nnGetCalls > 1 ) CYCLE
       ENDIF

       ! Get pointer to ESMF EXPORT field and pass data to it (if found):

       ! 2D...
       IF ( ThisDiagn%SpaceDim == 2 ) THEN
          ! Use pure ESMF operations to get field from export state
          ExportState = HcoState%EXPORT

          ! Try to get the field from export state
          CALL ESMF_StateGet(ExportState, itemName=TRIM(ThisDiagn%cName), field=Field, rc=lstat)
          IF (lstat == ESMF_SUCCESS) THEN
             ! Get the local array from the field
             CALL ESMF_FieldGet(Field, localDe=0, farrayPtr=Ptr2D, rc=lstat)

             IF (lstat == ESMF_SUCCESS .AND. ASSOCIATED(Ptr2D)) THEN
                IF ( ASSOCIATED(ThisDiagn%Arr2D) ) THEN
                   Ptr2D = ThisDiagn%Arr2D%Val
                ENDIF
             ENDIF
          ENDIF

       ! ... or 3D
       ELSEIF ( ThisDiagn%SpaceDim == 3 ) THEN
          ! Use pure ESMF operations to get field from export state
          ExportState = HcoState%EXPORT

          ! Try to get the field from export state
          CALL ESMF_StateGet(ExportState, itemName=TRIM(ThisDiagn%cName), field=Field, rc=lstat)
          IF (lstat == ESMF_SUCCESS) THEN
             ! Get the local array from the field
             CALL ESMF_FieldGet(Field, localDe=0, farrayPtr=Ptr3D, rc=lstat)

             IF (lstat == ESMF_SUCCESS .AND. ASSOCIATED(Ptr3D)) THEN
                IF ( ASSOCIATED(ThisDiagn%Arr3D) ) THEN
                   Ptr3D(:,:,:) = ThisDiagn%Arr3D%Val(:,:,HcoState%NZ:1:-1)
                ENDIF
             ENDIF
          ENDIF
       ENDIF

       ! Free pointer
       Ptr2D => NULL()
       Ptr3D => NULL()
    ENDDO

    ! Cleanup
    ThisDiagn => NULL()

    ! Return
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_Write
!EOC
END MODULE HCOIO_Write_Mod
#endif