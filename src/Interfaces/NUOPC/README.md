# HEMCO NUOPC Interface

This directory contains the HEMCO NUOPC (NASA Unified Observation Processing and Communications) interface module that provides pure ESMF functionality without MAPL dependencies.

## Overview

The HEMCO NUOPC interface (`hcoi_nuopc_mod.F90`) provides a clean, pure ESMF implementation that can be used as an alternative to the MAPL-based interface. This implementation:

- Eliminates MAPL dependencies
- Maintains identical public APIs to the MAPL interface
- Supports all the same functionality as the original interface
- Uses standard ESMF operations instead of MAPL-specific functions

## Public Interface

The module provides the following public functions:

- `HCO_SetServices_NUOPC` - Registers HEMCO data with ESMF import state
- `HCO_SetExtState_NUOPC` - Populates ExtState from ESMF import state
- `HCO_Imp2Ext_NUOPC` - Copies fields from import to ExtState (with 2D/3D, real/single/integer variants)

## Build Configuration

To enable the NUOPC interface, define the `NUOPC_ESMF` preprocessor flag during compilation:

```cmake
target_compile_definitions(HCOI_NUOPC
    PRIVATE
        ESMF_             # Enable ESMF functionality
        NUOPC_ESMF        # Enable NUOPC-specific code paths
)
```

The build system will automatically link against the ESMF library when the NUOPC interface is enabled.

## Conditional Compilation

The implementation uses the following conditional compilation patterns:

```fortran
#if defined(NUOPC_ESMF)
  ! Pure ESMF implementation
#elif defined(MAPL_ESMF)
  ! MAPL-specific implementation
#else
  ! Standalone implementation
#endif
```

## Core Module Integration

The following core modules have been updated to support both MAPL and NUOPC interfaces:

- `hco_error_mod.F90` - Handles missing value definitions and error handling
- `hco_restart_mod.F90` - Supports restart operations with both interfaces
- `hco_state_mod.F90` - Maintains ESMF state objects

## Backward Compatibility

The implementation maintains full backward compatibility with existing MAPL-based implementations. When `MAPL_ESMF` is defined, the original behavior is preserved.

## Usage

To use the NUOPC interface, simply define `NUOPC_ESMF` instead of `MAPL_ESMF` when building HEMCO. All existing APIs remain identical, ensuring seamless transition between interfaces.