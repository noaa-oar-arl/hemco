!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hemcox_dustginoux_mod.F90
!
! !DESCRIPTION: Paul GINOUX dust source function.  This subroutine updates
!  the surface mixing ratio of dust aerosols for NDSTBIN size bins.  The
!  uplifting of dust depends in space on the source function, and in time
!  and space on the soil moisture and surface wind speed (10 meters).  Dust
!  is uplifted if the wind speed is greater than a threshold velocity which
!  is calculated with the formula of Marticorena et al.  (JGR, v.102,
!  pp 23277-23287, 1997).  To run this subroutine you need the source
!  function which can be obtained by contacting Paul Ginoux at
!  ginoux@rondo.gsfc.nasa.gov/  If you are not using GEOS DAS met fields,
!  you will most likely need to adapt the adjusting parameter.
!\\
!\\
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities.
!\\
!\\
! References:
!
! \begin{enumerate}
! \item Ginoux, P., M. Chin, I. Tegen, J. Prospero, B. Hoben, O. Dubovik,
!        and S.-J. Lin, "Sources and distributions of dust aerosols simulated
!        with the GOCART model", J. Geophys. Res., 2001
! \item Chin, M., P. Ginoux, S. Kinne, B. Holben, B. Duncan, R. Martin,
!        J. Logan, A. Higurashi, and T. Nakajima, "Tropospheric aerosol
!        optical thickness from the GOCART model and comparisons with
!        satellite and sunphotometers measurements", J. Atmos Sci., 2001.
! \end{enumerate}
!
! !AUTHOR:
!  Barry D. Baker (barry.baker@noaa.gov)
!
! !INTERFACE:
!
MODULE HCOX_DustFengsha_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE HCO_State_Mod,  ONLY : HCO_State
  USE HCOX_State_Mod, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HcoX_DustFengsha_Run
  PUBLIC :: HcoX_DustFengsha_Init
  PUBLIC :: HcoX_DustFengsha_Final
  PUBLIC :: HcoX_DustFengsha_GetChDust
!
! !REVISION HISTORY:
!  08 Apr 2004 - T. D. Fairlie - Initial version
!  (1 ) Added OpenMP parallelization (bmy, 4/8/04)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  26 Feb 2013 - R. Yantosca - Now accept Input_Opt via the arg list
!  11 Dec 2013 - C. Keller   - Now a HEMCO extension.
!  29 Sep 2014 - R. Yantosca - Now make NBINS a variable and not a parameter
!  29 Sep 2014 - R. Yantosca - Now use F90 free-format indentation
!  08 Jul 2015 - M. Sulprizio- Now include dust alkalinity source (tdf 04/10/08)
!  14 Oct 2016 - C. Keller   - Now use HCO_EvalFld instead of HCO_GetPtr.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Quantities related to dust bins
  INTEGER, PARAMETER   :: NBINS = 4
  INTEGER              :: ExtNr    = -1     ! Extension number  for DustFengsha
  INTEGER              :: ExtNrAlk = -1     ! Extension number  for DustAlk
  INTEGER, ALLOCATABLE :: HcoIDs    (:)     ! HEMCO species IDs for DustFengsha
  INTEGER, ALLOCATABLE :: HcoIDsAlk (:)     ! HEMCO species IDs for DustAlk
  INTEGER, ALLOCATABLE :: IPOINT    (:)     ! 1=sand, 2=silt, 3=clay
  REAL,    ALLOCATABLE :: FRAC_S    (:)     !
  REAL,    ALLOCATABLE :: DUSTDEN   (:)     ! dust density     [kg/m3]
  REAL,    ALLOCATABLE :: DUSTREFF  (:)     ! effective radius [um]

  ! Source functions (get from HEMCO core)
  REAL(hp), POINTER    :: SRCE_SAND(:,:) => NULL()
  REAL(hp), POINTER    :: SRCE_SILT(:,:) => NULL()
  REAL(hp), POINTER    :: SRCE_CLAY(:,:) => NULL()
  REAL(hp), POINTER    :: SRCE_TYP (:,:) => NULL()
  REAL(hp), POINTER    :: SRCE_SSM (:,:) => NULL()

  ! Transfer coefficient (grid-dependent)
  REAL(dp)             :: CH_DUST

  ! Fundamental physical constants
  REAL*8,  PARAMETER   :: GAS_CST_UNV      = 8.3144598d0
  REAL*8,  PARAMETER   :: MMW_H2O          = 1.8015259d-02
  REAL*8,  PARAMETER   :: MMW_DRY_AIR      = 28.97d-3
  REAL*8,  PARAMETER   :: CST_VON_KRM      = 0.4d0
  REAL*8,  PARAMETER   :: GRV_SFC          = 9.80665d0
  REAL*8,  PARAMETER   :: GAS_CST_DRY_AIR  = 287.0d0
  REAL*8,  PARAMETER   :: RDS_EARTH        = 6.37122d+6
  REAL*8,  PARAMETER   :: GAS_CST_H2O      = 461.65D0
  REAL*8,  PARAMETER   :: SPC_HEAT_DRY_AIR = 1005.0d0
  REAL*8,  PARAMETER   :: TPT_FRZ_PNT      = 273.15d0


CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustFengsha_Run
!
! !DESCRIPTION: Subroutine HcoX\_DustFengsha\_Run is the driver routine
! for the Paul Fengsha dust source function HEMCO extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoX_DustFengsha_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_Calc_Mod,     ONLY : HCO_EvalFld
    USE HCO_EmisList_Mod, ONLY : HCO_GetPtr
    USE HCO_FluxArr_Mod,  ONLY : HCO_EmisAdd
    USE HCO_Clock_Mod,    ONLY : HcoClock_First
    USE HCO_GEOTOOLS_MOD, ONLY : HCO_LANDTYPE
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root   ! Are we on the root CPU?
    TYPE(Ext_State), POINTER        :: ExtState    ! Options for this ext
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState    ! HEMCO state object
    INTEGER,         INTENT(INOUT)  :: RC          ! Success or failure?
!
! !REMARKS:
!    SRCE_FUNK Source function                               (-)
!              for 1: Sand, 2: Silt, 3: Clay
!                                                                             .
!    DUSTDEN   Dust density                                  (kg/m3)
!    DUSTREFF  Effective radius                              (um)
!    AD        Air mass for each grid box                    (kg)
!    NTDT      Time step                                     (s)
!    W10m      Velocity at the anemometer level (10meters)   (m/s)
!    GWET      Surface wetness                               (-)
!                                                                             .
!  Dust properties used in GOCART
!                                                                             .
!  Size classes: 01-1, 1-1.8, 1.8-3, 3-6 (um)
!  Radius: 0.7, 1.5, 2.5, 4  (um)
!  Density: 2500, 2650, 2650, 2650 (kg/m3)
!
! !REVISION HISTORY:

!EOP
!------------------------------------------------------------------------------
!BOC

!
! !LOCAL VARIABLES:
!
    ! SAVED scalars
!    LOGICAL, SAVE     :: FIRST = .TRUE.
    ! Parameters
    ! Scalars
    INTEGER           :: I, J, N, M, tmpID, stype
    LOGICAL           :: ERR
    REAL*8            :: W10M     ! 10 m Wind Speed
    REAL*8            :: DEN      ! Air Density at lowest model level
    REAL*8            :: DIAM
    REAL*8            :: SNOW     ! Snow height
    REAL*8            :: U_TS0    ! Dry Threshold Velocity
    REAL*8            :: U_TS     ! Wet threshold Velocity (u_ts0 * H / R)
    REAL*8            :: UST      ! Threshold Velocity
    REAL*8            :: t        ! 2m temperature
    REAL*8            :: P        ! sfc pressure
    REAL*8            :: H        ! Fecan Soil Moisture Correction
    REAL*8            :: Z0       ! roughness lengh
    REAL*8            :: R        ! Drag Partition
    REAL*8            :: Q        ! Horizontal Flux
    REAL*8            :: SSM      ! Sediment Supply Map
    REAL*8            :: BULKD    ! Soil Bulk Density [kg/m3]
    REAL*8            :: RHOA     ! Air Density
    REAL*8            :: VSOILM   ! Volumetric Soil Moisture [m3/m3]
    REAL*8            :: AREA_M2  ! grid area [m2]
    REAL*8            :: ORO      ! Orography
    REAL*8            :: CLAY, SAND, SILT ! Clay sand and silt fraction
    REAL*8            :: SRCE_P, REYNOL, ALPHA,  BETA
    REAL*8            :: GAMMA,  CW,     DTSRCE, G
    REAL              :: DSRC
    CHARACTER(LEN=63) :: MSG

    ! Arrays
    REAL*8, TARGET  :: DUST_EMI_TOTAL(HcoState%NX, HcoState%NY)
    REAL(hp), TARGET  :: FLUX(HcoState%NX,HcoState%NY,NBINS)
    REAL(hp), TARGET  :: FLUX_ALK(HcoState%NX,HcoState%NY,NBINS)
!    REAL(hp), TARGET  :: DTHRES(HcoState%NX, HcoState%NY)

    ! Pointers
    REAL(hp), POINTER :: Arr2D(:,:)

    !=======================================================================
    ! HCOX_DUSTFENGSHA_RUN begins here!
    !=======================================================================

    ! Return if extension is disabled
    IF ( .NOT. ExtState%DustFengsha ) RETURN

    ! Enter
    CALL HCO_ENTER(HcoState%Config%Err,'HCOX_DustFengsha_Run (hcox_dustfengsha_mod.F90)',RC)
!    write(*,*) 'HERE ->>>>>>>>>>>>>>>>>>>>>'
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set gravity at earth surface (cm/s^2)
    G  = HcoState%Phys%g0 * 1.0d2

    ! Emission timestep [s]
    DTSRCE  = HcoState%TS_EMIS

    ! Initialize total dust emissions array [kg/m2/s]
    DUST_EMI_TOTAL = 0 ! 0.0_hp
!    DTHRES = 0.0_hp

    ! Error check
    ERR     = .FALSE.

    ! Init
    FLUX     = 0.0_hp
    FLUX_ALK = 0.0_hp
    Arr2D    => NULL()

    !=================================================================
    ! Point to DUST source functions
    !=================================================================
!    IF ( HcoClock_First(HcoState%Clock,.TRUE.) ) THEN

       ! Sand
       CALL HCO_EvalFld ( am_I_Root, HcoState, 'SAND', SRCE_SAND, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Silt
       CALL HCO_EvalFld ( am_I_Root, HcoState, 'SILT', SRCE_SILT, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Clay
       CALL HCO_EvalFld ( am_I_Root, HcoState, 'CLAY', SRCE_CLAY, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Soil Bulk Density
       CALL HCO_EvalFld ( am_I_Root, HcoState, 'SOILT', SRCE_TYP, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Sediment Supply Map
       CALL HCO_EvalFld ( am_I_Root, HcoState, 'SSM', SRCE_SSM, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN


!    ENDIF

    !=================================================================
    ! Compute dust emisisons
    !=================================================================
       DUST_EMI_TOTAL = 0.0d0

!$OMP PARALLEL DO                                             &
!$OMP DEFAULT( SHARED )                                       &
!$OMP PRIVATE( I,      J,     M,      N,      DEN,   DIAM   ) &
!$OMP PRIVATE( REYNOL, ALPHA, BETA,   GAMMA,  U_TS0, U_TS   ) &
!$OMP PRIVATE( CW,     W10M,  SRCE_P, RC                    ) &
!$OMP SCHEDULE( DYNAMIC )
       ! Loop over grid boxes
       DO J = 1, HcoState%NY
          DO I = 1, HcoState%NX
             ! orography
             ! Ocean is 0 - land is 1 - ice is 2
             ! write(*,*) 'CALCULATE ORO'
             ! oro = ExtState%WLI%Arr%Val(I,J) !HCO_LANDTYPE(ExtState%WLI%Arr%Val(I,J), ExtState%ALBD%Arr%Val(I,J))

             if (oro .eq. 0 ) then
                cycle
             endif
             ! write(*,*) 'ORO1 cycle <<<<<<<<<<<<<<<<'
             if (oro.eq.2) then
                cycle
             endif
             ! write(*,*) 'ORO2 Cycle <<<<<<<<<<<<<<<<'
             ! fractional clay sand silt (0 -> 1)
             clay = SRCE_CLAY(i,j)
             sand = SRCE_SAND(i,j)
!             silt = SRCE_SILT(i,j)

             ! Z0
             ! write(*,*) 'GET z0'
             z0 = ExtState%Z0%Arr%Val(I,J)
             if (z0 > 0.2) then
                cycle ! write(*,*) 'Z0 > 0.25'
             end if
             ! write(*,*) 'Z0 Cycle <<<<<<<<<<<<<<<<<<<<<'

             if (sand .gt. 110) then
                cycle
             endif
             !write(*,*) 'clay sand cycle <<<<<<<<<<<<<<<<<<<<'
             !====================================================================
             ! FENGSHA uses the threshold values from Dale Gillettes many papers
             ! that is based on soil type
             ! Read from the namelist
             !====================================================================
             !             stype = NINT(SRCE_TYP(I,J))
             !             write(*,*) stype,SRCE_TYP(I,J), '<<<<<<<<<<<<<STYPE'
             !if (stype .lt. 1) then
             !cycle
             !end if
             
!             write(*,*) stype, 'stype  <<<<<<<<<<<<<<'
             call HcoX_DustFengsha_styp(clay,sand,MIN(100., 100. - clay - sand), stype)
             !             call HcoX_DustFengsha_utst(stype,u_ts0)

             ! USTAR friction velocity from the model [m/s]
             ! write(*,*) 'GET USTAR'
             ust = ExtState%USTAR%Arr%Val(I,J)

             ! get the drag partition
             ! write(*,*) 'CALL DRAG PARTITION'
             call HcoX_DustFengsha_drag(z0,R)

             ! get the soil moisture correction
             ! write(*,*) 'GET SOIL MOISTURE'
             VSOILM = ExtState%GWETTOP%Arr%Val(I,J)

             ! write(*,*) 'CALL GET Moisture Correction'
             call HcoX_DustFengsha_moisture(VSOILM, clay, sand, H)
             !call HcoX_DustFengsha_Shao_Moisture(VSOILM,H)

             ! air density [kg / m3]
             ! write(*,*) 'GET 10m Temperature'
             ! write(*,*) 'TEMP: ', ExtState%T2M%Arr%Val(I,J)
             t = ExtState%T2M%Arr%Val(I,J)
             !t = t + 273.15d0

             !             write(*,*) 'GET PRES_SFC'
!             write(*,*) 'P: ', ExtState%PRES_SFC%Arr%Val(I,J)
!             P = ExtState%PRES_SFC%Arr%Val(I,J) !* 100.0d0
             !             p = p * 100
             p = 100000.
             rhoa = p / t / 287.05d0
             ! snow [ m H20 ]
             snow = ExtState%SNOWHGT%Arr%Val(I,J) !/ 1000.d0

             ! Calculate total vertical mass flux (note beta has units of m^-1)
             ! Beta acts to tone down dust in areas with so few dust-sized particles that the
             ! lofting efficiency decreases.  Otherwise, super sandy zones would be huge dust
             ! producers, which is generally not the case.  Equation derived from wind-tunnel
             ! experiments (see MB95).
             if (clay <= 0.2) then
                beta = 10**(13.4*clay-6.0)
             else
                beta = 2.E-4
             endif

             call HcoX_DustFengsha_hflux(ust, u_ts0 * H / R, Q)

             AREA_M2 = HcoState%Grid%AREA_M2%Val( I, J )
             SSM = SRCE_SSM(I,J)
             DUST_EMI_TOTAL(I,J) = beta * Q * SSM**2 * rhoa / GRV_SFC
             !         write(*,*) ust, snow, ssm, z0
             ! Increment total dust emissions [kg/m2/s] (L. Zhang, 6/26/15)
             DO N=1, NBINS
                SELECT CASE( N )
                CASE( 1 )
                   FLUX(I,J,N) = clay ! DUST_EMI_TOTAL(I,J) * 0.0766d0
                CASE( 2 )
                   FLUX(I,J,N) = sand ! DUST_EMI_TOTAL(I,J) * 0.1924d0
                CASE( 3 )
                   FLUX(I,J,N) = stype ! DUST_EMI_TOTAL(I,J) * 0.3491d0
                CASE( 4 )
                   FLUX(I,J,N) = DUST_EMI_TOTAL(I,J) ! DUST_EMI_TOTAL(I,J) * 0.3819d0
                END SELECT
                IF ( ExtNrAlk > 0 ) THEN
                   FLUX_ALK(I,J,N) = 0.04 * FLUX(I,J,N)
                END IF
             END DO ! N

          ENDDO ! NX
       ENDDO ! NY

!$OMP END PARALLEL DO

    ! Error check
    IF ( ERR ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    !=======================================================================
    ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS
    !=======================================================================
    DO N = 1, NBINS
       IF ( HcoIDs(N) > 0 ) THEN

          ! Add flux to emission array
          CALL HCO_EmisAdd( am_I_Root, HcoState, FLUX(:,:,N), &
                            HcoIDs(N), RC,       ExtNr=ExtNr   )
          IF ( RC /= HCO_SUCCESS ) THEN
             WRITE(MSG,*) 'HCO_EmisAdd error: dust bin ', N
             CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
             RETURN
          ENDIF

       ENDIF

       ! IF ( HcoIDsAlk(N) > 0 ) THEN

       !    ! Add flux to emission array
       !    CALL HCO_EmisAdd( am_I_Root,    HcoState, FLUX_Alk(:,:,N), &
       !                      HcoIDsAlk(N), RC,       ExtNr=ExtNrAlk   )
       !    IF ( RC /= HCO_SUCCESS ) THEN
       !       WRITE(MSG,*) 'HCO_EmisAdd error: dust alkalinity bin ', N
       !       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       !       RETURN
       !    ENDIF

       ! ENDIF

    ENDDO

    ! CALL HCO_EmisAdd( am_I_Root, HcoState, DTHRES, 5, RC, ExtNr=ExtNr )
    ! IF ( RC /= HCO_SUCCESS ) THEN
    !    WRITE(MSG,*) 'HCO_EmisAdd error: DTHRES ',5
    !    CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
    !    return
    ! END IF

    ! Leave w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HcoX_DustFengsha_Run


  subroutine HcoX_DustFengsha_utst(styp,ut)
    integer,                            intent(in)  :: styp
    real*8,                 intent(out) :: ut
!    ut = uth(styp)
    real *8 :: uth(13) = &
         (/ 0.07,   & ! Sand          - 1
         0.15,    & ! Loamy Sand      - 2
         0.27,    & ! Sandy Loam      - 3
         0.40,    & ! Silt Loam       - 4
         0.45,    & ! Silt            - 5
         0.45,    & ! Loam            - 6
         0.55,    & ! Sandy Clay Loam - 7
         0.45,    & ! Silty Clay Loam - 8
         0.50,    & ! Clay Loam       - 9
         0.45,    & ! Sandy Clay      - 10
         0.50,    & ! Silty Clay      - 11
         0.45,    & ! Clay            - 12
         9.999 /)   ! Other           - 13
    ut = uth(styp)
    return
  end subroutine HcoX_DustFengsha_utst

!   subroutine HcoX_DustFengsha_moisture(VSOILM, clay, soil_type, H)
!     ! INPUTS
!     real*8, intent(in)  :: VSOILM ! volumetric soil moisture [m3/m3]
!     real*8, intent(in)  :: clay   ! clay fraction [0-1]
!     integer, intent(in) :: soil_type
! !    real*8, intent(in)  :: bulkd
!     ! OUTPUTS
!     real*8, intent(out) :: H ! soil moisture correction - Fecan 1999

!     ! Local
!     real*8              :: gmoist ! gravimetric soil moisture
!     real*8              :: bulk_dens_dry
!     real*8              :: limit

!     ! parameters
!     real*8, parameter   :: bulk_dens = 2650.0d0
!     real*8, parameter   :: h20_dens = 1000.0d0
!     real*8              :: soilml1( 13 ) = (/ 0.395,    & ! Sand
!          0.410,    & ! Loamy Sand
!          0.435,    & ! Sandy Loam
!          0.485,    & ! Silt Loam
!          0.476,    & ! Silt
!          0.451,    & ! Loam
!          0.420,    & ! Sandy Clay Loam
!          0.477,    & ! Silty Clay Loam
!          0.476,    & ! Clay Loam
!          0.426,    & ! Sandy Clay
!          0.482,    & ! Silty Clay
!          0.482,    & ! Clay
!          0.482 /)    ! Other
!     ! Bulk density of dry surface soil  [kg m-3]
!     bulk_dens_dry = bulk_dens * ( 1.0d0 - VSOILM)
!     ! Gravimetric water content [ kg kg-1]
!     gmoist = 100 * VSOILM * h20_dens / bulk_dens_dry
!     if (gmoist.ge.1e10) then
!        gmoist = 0.
!     endif

!     ! gravimetric water threshold
!     limit = .0014* clay ** 2 + .17 * clay
!     if (gmoist > limit) then
!        H = sqrt( 1.0d0 + 1.21D0 * (gmoist - soilml1(soil_type)) ** 0.68D0 )
!     else
!        H = 1.0d0
!     endif
!     return
!   end subroutine HcoX_DustFengsha_moisture
  subroutine HcoX_DustFengsha_moisture(VSOILM, clay, sand, H)
    ! INPUTS
    real*8, intent(in)  :: VSOILM ! volumetric soil moisture [m3/m3]
    real*8, intent(in)  :: clay   ! clay fraction [0-1]
    real*8, intent(in) :: sand
!    real*8, intent(in)  :: bulkd
    ! OUTPUTS
    real*8, intent(out) :: H ! soil moisture correction - Fecan 1999

    ! Local
    real*8              :: gwc ! gravimetric soil moisture
    real*8              :: bulk_dens_dry
    real*8              :: limit
    real*8              :: wsat
    real*8              :: wdry
    real*8              :: wopt
    real*8              :: mpot
    real*8              :: b      !exponent b

    ! parameters
    real*8, parameter   :: bulk_dens = 2650.0d0
    real*8, parameter   :: h20_dens = 1000.0d0

    ! exponent for soil maptric potential [ unitless ]
    b = 2.91 + 0.159 * clay !* 100.

    ! saturated soil matric potential [ mm H2O ]
    mpot = 10 * (10 ** (1.88 - 0.0131 * sand  )) ! * 100.))

    ! saturated volumentric water content [ m3 m-3 ]
    wsat = 0.489 - 0.00126 * sand !* 100.

    ! volumetric water content of dry soil [ m3 m-3 ]
    ! wdry = wsat * (316230 / mpot ** (1/ b))
    
    ! Bulk density of dry surface soil  [kg m-3]
    bulk_dens_dry = bulk_dens * ( 1.0d0 - wsat)
    ! Gravimetric water content [ kg kg-1]
    gwc = VSOILM * h20_dens / bulk_dens_dry
    if (gwc.ge.1e10) then
       gwc = 0.
    endif

    ! gravimetric water threshold
    limit = .0014* clay ** 2 + .17 * clay
    if (gwc > limit) then
       H = sqrt( 1. + 1.21 * (100 * (gwc - limit) ) ** 0.68)
    else
       H = 1.0d0
    endif
    return
  end subroutine HcoX_DustFengsha_moisture

  subroutine HcoX_DustFengsha_Shao_Moisture(w, H)
    ! INPUTS
    real*8, intent(in)  :: w ! volumetric soil moisture [m3/m3]
    real*8, intent(out) :: H ! soil moisture correction

    if (w > 0.03) then
       H = exp(95.3 * w - 2.029)
    else
       H = exp(22.7 * w)
    endif
    return
  end subroutine HcoX_DustFengsha_Shao_Moisture

  subroutine HcoX_DustFengsha_styp(clay, sand, silt, type)
    !---------------------------------------------------------------
    ! Function: calculate soil type based on USDA definition.
    ! Source: USDA soil texture calculator
    !
    !---------------------------------------------------------------
    REAL*8, intent(in) ::  clay, sand, silt
    integer, intent(out) ::  type
 !    real*8 :: cly, snd, slt
    type = 0

    if (silt+1.5*clay .lt. 15)                                                                type = 1      ! sand
    if (silt+1.5*clay .ge. 15 .and.silt+1.5*clay .lt. 30)                                       type = 2      ! loamy sand
    if (clay .ge. 7 .and. clay .lt. 20 .and. sand .gt. 52 .and. silt+2*clay .ge. 30)             type = 3      ! sandy loam (cond 1)
    if (clay .lt. 7 .and. silt .lt. 50 .and. silt+2*clay .ge. 30)                               type = 3      ! sandy loam (cond 2)
    if (silt .ge. 50 .and. clay .ge. 12 .and.clay .lt. 27 )                                    type = 4      ! silt loam (cond 1)
    if (silt .ge. 50 .and. silt .lt. 80 .and.clay .lt. 12)                                     type = 4      ! silt loam (cond 2)
    if (silt .ge. 80 .and. clay .lt. 12)                                                      type = 5      ! silt
    if (clay .ge. 7  .and. clay .lt. 27 .and.silt .ge. 28 .and. silt .lt. 50 .and.sand .le. 52)  type = 6      ! loam
    if (clay .ge. 20 .and. clay .lt. 35 .and.silt .lt. 28 .and. sand .gt. 45)                   type = 7      ! sandy clay loam
    if (clay .ge. 27 .and. clay .lt. 40 .and.sand .lt. 20)                                     type = 8      ! silt clay loam
    if (clay .ge. 27 .and. clay .lt. 40 .and.sand .ge. 20 .and. sand .le. 45)                   type = 9      ! clay loam
    if (clay .ge. 35 .and. sand .gt. 45)                                                      type = 10     ! sandy clay
    if (clay .ge. 40 .and. silt .ge. 40)                                                      type = 11     ! silty clay
    if (clay .ge. 40 .and. sand .le. 45 .and.silt .lt. 40)                                     type = 12     ! clay
    if ( type .eq. 0) type = 13
    return
  end subroutine HcoX_DustFengsha_styp

  subroutine HcoX_DustFengsha_drag(z0,R)
    real*8, intent(in) :: z0
    real*8, intent(out) :: R
    real, parameter :: z0s = 1.0e-04 !Surface roughness for ideal bare surface [m]
    ! ------------------------------------------------------------------------
    ! Function: Calculates the MacKinnon et al. 2004 Drag Partition Correction
    !
    !   R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)
    !
    !--------------------------------------------------------------------------
    ! Drag partition correction. See MacKinnon et al. (2004),
    !     doi:10.1016/j.geomorph.2004.03.009
    R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)
    ! Drag partition correction. See Marticorena et al. (1997),
    !     doi:10.1029/96JD02964
    !R = 1.0 - log(z0 / z0s) / log( 0.7 * (10./z0s) ** 0.8)
    return
  end subroutine HcoX_DustFengsha_drag

  subroutine HcoX_DustFengsha_hflux(ust,utst, Q)
    !---------------------------------------------------------------------
    ! Function: Calculates the Horizontal Saltation Flux, Q, and then
    !           calculates the vertical flux.
    !
    ! formula of Draxler & Gillette (2001) Atmos. Environ.
    ! Q = U* ( U*^2 - Ut*^2 ) = u*^3 ( 1 - (Ut* / u*) ^ 2 )
    !
    ! where:
    !     F   = vertical emission flux  [g/m**2-s]
    !     K   = constant 2.0E-04                      [1/m]
    !     A   = 0~3.5  mean = 2.8  (fudge factor)
    !     U*  = friction velocity                     [m/s]
    !     Ut* = threshold friction velocity           [m/s]
    !
    !--------------------------------------------------------------------
    real*8, intent(in) :: ust    ! friction velocity
    real*8, intent(in) :: utst   ! modified threshold velocity
    real*8, intent(out) :: Q
    !    write(*,*) 'UST:',ust, 'UTST:', utst
    if (ust < utst) then
       Q = 0.
    else
       Q = ust * (ust * ust - utst * utst)
    end if
    return
  end subroutine HcoX_DustFengsha_hflux
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustFengsha_Init
!
! !DESCRIPTION: Subroutine HcoX\_DustFengsha\_distribution calculates the Kok
! particle distribution
!\\
!\\
! !INTERFACE:
!   SUBROUTINE HcoX_DustFengsha_distribution(flux)

! !
! ! !USES:
! !

! !
! ! !INPUT PARAMETERS:
! !
!     real*8           :: flux
! !
! ! !INPUT/OUTPUT PARAMETERS:
! !
!       ! -- dust parameters
!     integer,            dimension(ndust), parameter :: ipoint    = (/ 3, 2, 2, 2, 2 /)
!     real*8, dimension(ndust), parameter :: den_dust  = (/   2500.,  2650.,  2650.,  2650.,  2650. /)
!     real*8, dimension(ndust), parameter :: reff_dust = (/ 0.73D-6, 1.4D-6, 2.4D-6, 4.5D-6, 8.0D-6 /)
!     real*8, dimension(ndust), parameter :: frac_s    = (/     0.1,   0.25,   0.25,   0.25,   0.25 /)
!     real*8, dimension(ndust), parameter :: lo_dust   = (/  0.1D-6, 1.0D-6, 1.8D-6, 3.0D-6, 6.0D-6 /)
!     real*8, dimension(ndust), parameter :: up_dust   = (/  1.0D-6, 1.8D-6, 3.0D-6, 6.0D-6,10.0D-6 /)
!     real*8, dimension(ndust, 12)        :: ch_dust   = 0.8e-09_CHEM_KIND_R8
!     dvol_tot=0.
!     DO n=1,nmx
!        dlndp(n)=LOG(up_dust(n)/lo_dust(n))
!        dvol(n)=(2.0*reff_dust(n)/cv)*(1.+ERF(LOG(2.0*reff_dust(n)/mmd_dust)/(SQRT(2.)*LOG(gsd_dust))))*&
!             EXP(-(2.0*reff_dust(n)/lambda)**3.0)*dlndp(n)
!        dvol_tot=dvol_tot+dvol(n)
!        ! Convert mass flux to volume flux
!        emit_vol=emit/den_dust(n) ! (m s^-1)
!     END DO
!     DO n=1,nmx
!        distr_dust(n)=dvol(n)/dvol_tot
!        !print *,"distr_dust(",n,")=",distr_dust(n)
!     END DO
!     ! Now distribute total vertical emission into dust bins and update concentration.
!     DO n=1,nmx
!        DO i=1,imx
!           DO j=1,jmx
!              ! Calculate total mass emitted
!              dsrc = emit_vol*den_dust(n)*distr_dust(n)*dxy(j)*dt1  ! (kg)
!              IF (dsrc < 0.0) dsrc = 0.0
!              ! Update dust mixing ratio at first model level.
!              tc(i,j,1,n) = tc(i,j,1,n) + dsrc / airmas(i,j,1) ! (kg/kg)
!              !   bems(i,j,n) = dsrc  ! diagnostic
!              bems(i,j,n) = 1000.*dsrc/(dxy(j)*dt1) ! diagnostic (g/m2/s)
!           END DO
!        END DO
!     END DO
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustFengsha_Init
!
! !DESCRIPTION: Subroutine HcoX\_DustFengsha\_Init initializes the HEMCO
! DUSTFENGSHA extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoX_DustFengsha_Init( am_I_Root, HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod, ONLY : GetExtNr, GetExtOpt
    USE HCO_State_Mod,   ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO State object
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName    ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState   ! Extension options
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  11 Dec 2013 - C. Keller   - Now a HEMCO extension
!  26 Sep 2014 - R. Yantosca - Updated for TOMAS
!  29 Sep 2014 - R. Yantosca - Now initialize NBINS from HcoState%N_DUST_BINS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                        :: N, AS, nSpc, nSpcAlk
    CHARACTER(LEN=255)             :: MSG
    REAL(dp)                       :: Mp, Rp, TmpScal
    LOGICAL                        :: FOUND

    ! Arrays
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNamesAlk(:)

    !=======================================================================
    ! HCOX_DUSTFENGSHA_INIT begins here!
    !=======================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Check for dust alkalinity option
    ExtNrAlk = GetExtNr( HcoState%Config%ExtList, 'DustAlk' )

    ! Enter
    CALL HCO_ENTER(HcoState%Config%Err,'HCOX_DustFengsha_Init (hcox_dustfengsha_mod.F90)',RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get the expected number of dust species
!    NBINS = HcoState%nDust

    ! Get the actual number of dust species defined for DustFengsha extension
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get the dust alkalinity species defined for DustAlk option
    IF ( ExtNrAlk > 0 ) THEN
       CALL HCO_GetExtHcoID( HcoState,    ExtNrAlk, HcoIDsAlk, &
                             SpcNamesAlk, nSpcAlk,  RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Make sure the # of dust species is as expected
    IF ( nSpc /= NBINS ) THEN
       WRITE( MSG, 100 ) NBINS, nSpc
 100   FORMAT( 'Expected ', i3, ' DustFengsha species but only found ', i3, &
               ' in the HEMCO configuration file!  Exiting...' )
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! Set scale factor: first try to read from configuration file. If
    ! not specified, call wrapper function which sets teh scale factor
    ! based upon compiler switches.
    CALL GetExtOpt( HcoState%Config, ExtNr, 'Mass tuning factor', &
                     OptValDp=TmpScal, Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set parameter FLX_MSS_FDG_FCT to specified tuning factor. Get from
    ! wrapper routine if not defined in configuration file
    IF ( FOUND ) THEN
       CH_DUST = TmpScal
    ELSE
       ! Get global mass flux tuning factor
       CH_DUST = HcoX_DustFengsha_GetCHDust()
       IF ( CH_DUST < 0.0_dp ) THEN
          RC = HCO_FAIL
          RETURN
       ENDIF
    ENDIF

    ! Verbose mode
    IF ( am_I_Root ) THEN
       MSG = 'Use Fengsha dust emissions (extension module)'
       CALL HCO_MSG(HcoState%Config%Err,MSG )

       IF ( ExtNrAlk > 0 ) THEN
          MSG = 'Use dust alkalinity option'
          CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
       ENDIF

       MSG = 'Use the following species (Name: HcoID):'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       DO N = 1, nSpc
          WRITE(MSG,*) TRIM(SpcNames(N)), ':', HcoIDs(N)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDDO
       IF ( ExtNrAlk > 0 ) THEN
          DO N = 1, nSpcAlk
             WRITE(MSG,*) TRIM(SpcNamesAlk(N)), ':', HcoIDsAlk(N)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDDO
       ENDIF

       WRITE(MSG,*) 'Global mass flux tuning factor: ', CH_DUST
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP2='-')
    ENDIF

    ! Allocate vectors holding bin-specific informations
    ALLOCATE ( IPOINT  (NBINS) )
    ALLOCATE ( FRAC_S  (NBINS) )
    ALLOCATE ( DUSTDEN (NBINS) )
    ALLOCATE ( DUSTREFF(NBINS) )

    ! Allocate arrays
    ALLOCATE ( SRCE_SAND ( HcoState%NX, HcoState%NY ), &
               SRCE_SILT ( HcoState%NX, HcoState%NY ), &
               SRCE_CLAY ( HcoState%NX, HcoState%NY ), &
               SRCE_TYP  ( HcoState%NX, HcoState%NY ), &
               SRCE_SSM  ( HcoState%NX, HcoState%NY ), &
               STAT = AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR(HcoState%Config%Err,'Allocation error', RC )
       RETURN
    ENDIF
    SRCE_SAND = 0.0_hp
    SRCE_SILT = 0.0_hp
    SRCE_CLAY = 0.0_hp
    SRCE_TYP  = 0.0_hp
    SRCE_SSM  = 0.0_hp

    !=======================================================================
    ! Setup for simulations that use 4 dust bins (w/ or w/o TOMAS)
    !=======================================================================

    ! Fill bin-specific information
    IF ( NBINS == 4 ) THEN

       IPOINT  (1:NBINS) = (/ 3,       2,       2,       2       /)
       FRAC_S  (1:NBINS) = (/ 0.095d0, 0.3d0,   0.3d0,   0.3d0   /)
       DUSTDEN (1:NBINS) = (/ 2500.d0, 2650.d0, 2650.d0, 2650.d0 /)
       DUSTREFF(1:NBINS) = (/ 0.73d-6, 1.4d-6,  2.4d-6,  4.5d-6  /)

    ELSE

#if !defined( TOMAS )
       MSG = 'Cannot have > 4 FENGSHA dust bins unless you are using TOMAS!'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
#endif

    ENDIF

#if defined( TOMAS )

    !=======================================================================
    ! Setup for TOMAS simulations using more than 4 dust bins
    !
    ! from Fengsha:
    ! The U.S. Department of Agriculture (USDA) defines particles
    ! with a radius between 1 um and 25 um as silt, and below 1 um
    ! as clay [Hillel, 1982]. Mineralogical silt particles are mainly
    ! composed of quartz, but they are often coated with strongly
    ! adherent clay such that their physicochemical properties are
    ! similar to clay [Hillel, 1982].
    !
    ! SRCE_FUNC Source function
    ! for 1: Sand, 2: Silt, 3: Clay
    !=======================================================================
    IF ( NBINS == HcoState%MicroPhys%nBins ) THEN

       !--------------------------------------------------------------------
       ! Define the IPOINT array based on particle size
       !--------------------------------------------------------------------

       ! Loop over # of TOMAS bins
       DO N = 1, HcoState%MicroPhys%nBins

          ! Compute particle mass and radius
          Mp = 1.4 * HcoState%MicroPhys%BinBound(N)
          Rp = ( ( Mp /2500. ) * (3./(4.*HcoState%Phys%PI)))**(0.333)

          ! Pick the source function based on particle size
          IF ( Rp < 1.d-6 ) THEN
             IPOINT(N) = 3
          ELSE
             IPOINT(N) = 2
          END IF
       END DO

       !--------------------------------------------------------------------
       ! Set up dust density (DUSTDEN) array
       !--------------------------------------------------------------------
       DO N = 1, HcoState%MicroPhys%nBins
          IF ( HcoState%MicroPhys%BinBound(N) < 4.0D-15 ) THEN
             DUSTDEN(N)  = 2500.d0
          ELSE
             DUSTDEN(N)  = 2650.d0
          ENDIF
       ENDDO

       !--------------------------------------------------------------------
       ! Set up dust density (DUSTDEN) array
       !--------------------------------------------------------------------
       DO N = 1, HcoState%MicroPhys%nBins
          DUSTREFF(N) = 0.5d0                                              &
                      * ( SQRT( HcoState%MicroPhys%BinBound(N) *      &
                                HcoState%MicroPhys%BinBound(N+1) )    &
                      /   DUSTDEN(N) * 6.d0/HcoState%Phys%PI )**( 0.333d0 )
       ENDDO

       !--------------------------------------------------------------------
       ! Set up the FRAC_S array
       !--------------------------------------------------------------------

       ! Initialize
       FRAC_S( 1:HcoState%MicroPhys%nBins )           = 0d0

# if  defined( TOMAS12 ) || defined( TOMAS15 )

       !---------------------------------------------------
       ! TOMAS simulations with 12 or 15 size bins
       !---------------------------------------------------
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 1  )  = 7.33E-10
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 2  )  = 2.032E-08
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 3  )  = 3.849E-07
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 4  )  = 5.01E-06
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 5  )  = 4.45E-05
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 6  )  = 2.714E-04
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 7  )  = 1.133E-03
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 8  )  = 3.27E-03
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 9  )  = 6.81E-03
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 10 )  = 1.276E-02
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 11 )  = 2.155E-01
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 12 )  = 6.085E-01

# else

       !---------------------------------------------------
       ! TOMAS simulations with 30 or 40 size bins
       !---------------------------------------------------
       FRAC_S( HcoState%MicroPhys%nActiveModeBins +  1 )  = 1.05d-10
       FRAC_S( HcoState%MicroPhys%nActiveModeBins +  2 )  = 6.28d-10
       FRAC_S( HcoState%MicroPhys%nActiveModeBins +  3 )  = 3.42d-09
       FRAC_S( HcoState%MicroPhys%nActiveModeBins +  4 )  = 1.69d-08
       FRAC_S( HcoState%MicroPhys%nActiveModeBins +  5 )  = 7.59d-08
       FRAC_S( HcoState%MicroPhys%nActiveModeBins +  6 )  = 3.09d-07
       FRAC_S( HcoState%MicroPhys%nActiveModeBins +  7 )  = 1.15d-06
       FRAC_S( HcoState%MicroPhys%nActiveModeBins +  8 )  = 3.86d-06
       FRAC_S( HcoState%MicroPhys%nActiveModeBins +  9 )  = 1.18d-05
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 10 )  = 3.27d-05
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 11 )  = 8.24d-05
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 12 )  = 1.89d-04
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 13 )  = 3.92d-04
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 14 )  = 7.41d-04
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 15 )  = 1.27d-03
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 16 )  = 2.00d-03
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 17 )  = 2.89d-03
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 18 )  = 3.92d-03
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 19 )  = 5.26d-03
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 20 )  = 7.50d-03
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 21 )  = 1.20d-02
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 22 )  = 2.08d-02
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 23 )  = 3.62d-02
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 24 )  = 5.91d-02
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 25 )  = 8.74d-02
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 26 )  = 1.15d-01
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 27 )  = 1.34d-01
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 28 )  = 1.37d-01
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 29 )  = 1.24d-01
       FRAC_S( HcoState%MicroPhys%nActiveModeBins + 30 )  = 9.85d-02

# endif

    ELSE

       ! Stop w/ error message
       CALL HCO_ERROR( HcoState%Config%Err, 'Wrong number of TOMAS dust bins!', RC )

    ENDIF

#endif

    !=====================================================================
    ! Activate fields in ExtState used by Fengsha dust
    !=====================================================================

    ! Activate met. fields required by this module
    ExtState%USTAR%DoUse   = .TRUE.
    ExtState%Z0%DoUse      = .TRUE.
    ExtState%PRES_SFC%DoUse     = .TRUE.
    ExtState%ALBD%DoUse    = .TRUE.
    ExtState%WLI%DoUse     = .TRUE.
    ExtState%SNOWHGT%DoUse = .TRUE.
    ExtState%T2M%DoUse     = .TRUE.
    ExtState%GWETTOP%DoUse = .TRUE.

    ! Activate this module
    ExtState%DustFengsha    = .TRUE.

    ! Leave w/ success
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HcoX_DustFengsha_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustFengsha_Final
!
! !DESCRIPTION: Subroutine HcoX\_DustFengsha\_Final finalizes the HEMCO
! DUSTFENGSHA extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoX_DustFengsha_Final()
!
! !REVISION HISTORY:
!  11 Dec 2013 - C. Keller - Now a HEMCO extension
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! HCOX_DUSTFENGSHA_FINAL begins here!
    !=======================================================================

    ! Free pointer
    IF ( ASSOCIATED( SRCE_SAND ) ) DEALLOCATE( SRCE_SAND )
    IF ( ASSOCIATED( SRCE_SILT ) ) DEALLOCATE( SRCE_SILT )
    IF ( ASSOCIATED( SRCE_CLAY ) ) DEALLOCATE( SRCE_CLAY )
    IF ( ASSOCIATED( SRCE_TYP  ) ) DEALLOCATE( SRCE_TYP  )
    IF ( ASSOCIATED( SRCE_SSM  ) ) DEALLOCATE( SRCE_SSM  )

    ! Cleanup option object
    IF ( ALLOCATED( IPOINT    ) ) DEALLOCATE( IPOINT    )
    IF ( ALLOCATED( FRAC_S    ) ) DEALLOCATE( FRAC_S    )
    IF ( ALLOCATED( DUSTDEN   ) ) DEALLOCATE( DUSTDEN   )
    IF ( ALLOCATED( DUSTREFF  ) ) DEALLOCATE( DUSTREFF  )
    IF ( ALLOCATED( HcoIDs    ) ) DEALLOCATE( HcoIDs    )
    IF ( ALLOCATED( HcoIDsALK ) ) DEALLOCATE( HcoIDsALK )

  END SUBROUTINE HcoX_DustFengsha_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustFengsha_GetChDust
!
! !DESCRIPTION: Function HCOX\_DustFengsha\_GetChDust returns the CH\_DUST
! parameter for the current simulation type.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCOX_DustFengsha_GetChDust() RESULT( CH_DUST )
!
! !RETURN VALUE:
!
    REAL*8 :: CH_DUST
!
! !REMARKS:
!  The logic in the #ifdefs may need to be cleaned up later on.  We have
!  just replicated the existing code in pre-HEMCO versions of dust_mod.F.
!
! !REVISION HISTORY:
!  11 Dec 2013 - C. Keller   - Initial version
!  25 Sep 2014 - R. Yantosca - Updated for TOMAS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Transfer coeff for type natural source  (kg*s2/m5)
    ! Emission reduction factor for China-nested grid domain (win, 4/27/08)

#if defined( GRID4x5 )

    !-----------------------------------------------------------------------
    ! All 4x5 simulations (including TOMAS)
    !-----------------------------------------------------------------------
    CH_DUST  = 9.375d-10

#elif defined( GRID1x1 ) && defined( NESTED_CH )

    !-----------------------------------------------------------------------
    ! Note: monthly emission over the China nested-grid domain is about
    !       2 times higher compared to the same domain in 4x5 resolution
    !       Thus applying 1/2  factor to correct the emission.
    !
    !%%% NOTE: MAY NEED TO UPDATE THIS STATEMENT FOR HIGHER RESOLUTION
    !%%% NESTED GRIDS.  THIS WAS ORIGINALLY DONE FOR THE GEOS-3 1x1
    !%%% NESTED GRID.  LOOK INTO THIS LATER.  (bmy, 9/25/14)
    !-----------------------------------------------------------------------
    CH_DUST  = 9.375d-10 * 0.5d0

#else

    !-----------------------------------------------------------------------
    ! All other resolutions
    !-----------------------------------------------------------------------

    ! Start w/ same value as for 4x5
    CH_DUST  = 9.375d-10

#if defined( TOMAS )
    ! KLUDGE: For TOMAS simulations at grids higher than 4x5 (e.g. 2x25),
    ! then multiplyCH_DUST by 0.75.  (Sal Farina)
    CH_DUST  = CH_DUST * 0.75d0
#endif

#endif

  END FUNCTION HCOX_DustFengsha_GetChDust

!EOC
END MODULE HCOX_DustFengsha_Mod
