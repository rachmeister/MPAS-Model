! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_mixed_layer_eddies
!
!> \brief MPAS ocean initialize case -- Mixed Layer Eddies
!> \author Luke Van Roekel
!> \date   04/15/2019
!> \details
!>  This module contains the routines for initializing the
!>  the mixed layer eddies test case (Fox Kemper et al 2008)
!
!-----------------------------------------------------------------------

module ocn_init_mixed_layer_eddies

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar

   use ocn_constants
   use ocn_init_vertical_grids
   use ocn_init_cell_markers

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_init_setup_mixed_layer_eddies, &
             ocn_init_validate_mixed_layer_eddies

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_mixed_layer_eddies
!
!> \brief   Setup for mixed layer eddies test case
!> \author  Luke Van Roekel
!> \date    04/15/2019
!> \details
!>  This routine sets up the initial conditions for the mixed layer eddies test case.
!>  It should also ensure the mesh that was input is valid for the configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_mixed_layer_eddies(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr
      real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal
      real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal
      real (kind=RKIND) :: temperature, yOffset, xPerturbationMin, xPerturbationMax
      real (kind=RKIND) :: perturbationWidth, stratification, buoyancy

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: forcingPool

      integer :: iCell, k, idx

      ! Define config variable pointers
      character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid
      logical, pointer :: config_mixed_layer_eddies_use_shortwave
      real (kind=RKIND), pointer :: config_mixed_layer_eddies_background_stratification,  &
                                    config_mixed_layer_eddies_ML_stratification,          &
                                    config_mixed_layer_eddies_horizontal_buoyancy_grad,   &
                                    config_mixed_layer_eddies_frontal_width,              &
                                    config_mixed_layer_eddies_mixed_layer_depth,          &
                                    config_mixed_layer_eddies_surface_cooling,            &
                                    config_mixed_layer_eddies_background_buoyancy,        &
                                    config_mixed_layer_eddies_coriolis_parameter,         &
                                    config_mixed_layer_eddies_bottom_depth,               &
                                    config_mixed_layer_eddies_density_0,                  &
                                    config_mixed_layer_eddies_alpha,                      &
                                    config_mixed_layer_eddies_beta,                       &
                                    config_mixed_layer_eddies_tref,                       &
                                    config_mixed_layer_eddies_sref,                       &
                                    config_mixed_layer_eddies_windStressZonal,            &
                                    config_mixed_layer_eddies_windStressMeridional

      ! Define dimension pointers
      integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
      integer, pointer :: index_temperature, index_salinity, index_tracer1

      ! Define variable pointers
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell,refBottomDepth, refZMid, &
                                                  vertCoordMovementWeights, bottomDepth, &
                                                  fCell, fEdge, fVertex, dcEdge,         &
                                                  sensibleHeatFlux, windStressZonal,     &
                                                  windStressMeridional, shortWaveHeatFlux

      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, debugTracers
      real (kind=RKIND), parameter :: salinity = 35.0 !for now just have constant salinity

      ! Define local interfaceLocations variable
      real (kind=RKIND), dimension(:), pointer :: interfaceLocations

      logical, pointer :: on_a_sphere

      iErr = 0

      call mpas_pool_get_config(ocnConfigs, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('baroclinic_channel')) return

      call mpas_pool_get_config(ocnConfigs, 'config_vertical_grid', config_vertical_grid)

      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_background_stratification', &
                                config_mixed_layer_eddies_background_stratification)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_ML_stratification', &
                                config_mixed_layer_eddies_ML_stratification)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_horizontal_buoyancy_grad', &
                                config_mixed_layer_eddies_horizontal_buoyancy_grad)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_bottom_depth',              &
                                config_mixed_layer_eddies_bottom_depth)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_frontal_width', &
                                config_mixed_layer_eddies_frontal_width)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_mixed_layer_depth', &
                                config_mixed_layer_eddies_mixed_layer_depth)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_surface_cooling', &
                                config_mixed_layer_eddies_surface_cooling)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_background_buoyancy',     &
                                config_mixed_layer_eddies_background_buoyancy)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_coriolis_parameter', &
                                config_mixed_layer_eddies_coriolis_parameter)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_use_shortwave', &
                                config_mixed_layer_eddies_use_shortwave)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_windStressZonal', &
                                config_mixed_layer_eddies_windStressZonal)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_windStressMeridional', &
                                config_mixed_layer_eddies_windStressMeridional)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_density_0', &
                                config_mixed_layer_eddies_density_0)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_alpha', &
                                config_mixed_layer_eddies_alpha)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_beta', &
                                config_mixed_layer_eddies_beta)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_tref', &
                                config_mixed_layer_eddies_tref)
      call mpas_pool_get_config(ocnConfigs, 'config_mixed_layer_eddies_sref', &
                                config_mixed_layer_eddies_sref)

      ! Determine vertical grid for configuration
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      if ( on_a_sphere ) call mpas_log_write('The baroclinic channel configuration ' &
              // 'can only be applied to a planar mesh. Exiting...', MPAS_LOG_CRIT)

      allocate(interfaceLocations(nVertLevelsP1))
      call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations )

      ! Initalize min/max values to large positive and negative values
      yMin = 1.0E10_RKIND
      yMax = -1.0E10_RKIND
      xMin = 1.0E10_RKIND
      xMax = -1.0E10_RKIND
      dcEdgeMin = 1.0E10_RKIND

      ! Determine local min and max values.
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

        yMin = min( yMin, minval(yCell(1:nCellsSolve)))
        yMax = max( yMax, maxval(yCell(1:nCellsSolve)))
        xMin = min( xMin, minval(xCell(1:nCellsSolve)))
        xMax = max( xMax, maxval(xCell(1:nCellsSolve)))
        dcEdgeMin = min( dcEdgeMin, minval(dcEdge(1:nEdgesSolve)))

        block_ptr => block_ptr % next
      end do

      ! Determine global min and max values.
      call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
      call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
      call mpas_dmpar_min_real(domain % dminfo, xMin, xMinGlobal)
      call mpas_dmpar_max_real(domain % dminfo, xMax, xMaxGlobal)
      call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

      yMidGlobal = (yMinGlobal + yMaxGlobal) * 0.5_RKIND

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
        call mpas_pool_get_dimension(tracersPool, 'index_tracer1', index_tracer1)

        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'fCell', fCell)
        call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
        call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

        call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
        call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal, 1)
        call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional, 1)
        call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)

        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
        call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

        call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
        call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

        ! Set refBottomDepth and refZMid
        do k = 1, nVertLevels
            refBottomDepth(k) = config_mixed_layer_eddies_bottom_depth * interfaceLocations(k+1)
            refZMid(k) = - 0.5_RKIND * (interfaceLocations(k+1) + interfaceLocations(k)) * config_mixed_layer_eddies_bottom_depth
        end do

        ! Set vertCoordMovementWeights
        vertCoordMovementWeights(:) = 1.0_RKIND

        do iCell = 1, nCellsSolve
           ! Set debug tracer
           if ( associated(debugTracers) ) then
              idx = index_tracer1
              do k = 1, nVertLevels
                 debugTracers(idx, k, iCell) = 1.0_RKIND
              enddo
           end if

           ! Set stratification based on northern half of domain temperature
           if ( associated(activeTracers) ) then
              idx = index_temperature
              do k = nVertLevels, 1, -1

                 if ( refBottomDepth(k) < config_mixed_layer_eddies_mixed_layer_depth ) then
                   stratification = config_mixed_layer_eddies_ML_stratification * &
                                    config_mixed_layer_eddies_coriolis_parameter
                 else
                   stratification = config_mixed_layer_eddies_background_stratification * &
                                    config_mixed_layer_eddies_coriolis_parameter
                 endif

                 buoyancy = stratification*(config_mixed_layer_eddies_mixed_layer_depth - refBottomDepth(k)) + &
                          config_mixed_layer_eddies_frontal_width * config_mixed_layer_eddies_horizontal_buoyancy_grad * &
                          config_mixed_layer_eddies_coriolis_parameter * 0.5_RKIND * tanh(2.0_RKIND*(yCell(iCell) - yMidGlobal) / config_mixed_layer_eddies_frontal_width) &
                          + config_mixed_layer_eddies_background_buoyancy

                 temperature = config_mixed_layer_eddies_tref - buoyancy / (9.806_RKIND * config_mixed_layer_eddies_alpha)
                 activeTracers(idx, k, iCell) = temperature
              end do

                         ! Set salinity
              idx = index_salinity
              activeTracers(idx, :, iCell) = salinity
           end if

           ! Set layerThickness and restingThickness
           do k = 1, nVertLevels
              layerThickness(k, iCell) = config_mixed_layer_eddies_bottom_depth * ( interfaceLocations(k+1) &
                                       - interfaceLocations(k) )
              restingThickness(k, iCell) = config_mixed_layer_eddies_bottom_depth * ( interfaceLocations(k+1) &
                                         - interfaceLocations(k) )

           end do

           ! Set bottomDepth
           bottomDepth(iCell) = config_mixed_layer_eddies_bottom_depth

           ! Set maxLevelCell
           maxLevelCell(iCell) = nVertLevels
        end do

        ! add forcing, wind stress, sensible heat, and diurnal max
        sensibleHeatFlux(:) = config_mixed_layer_eddies_surface_cooling 
        windStressZonal(:) = config_mixed_layer_eddies_windStressZonal
        windStressMeridional(:) = config_mixed_layer_eddies_windStressMeridional

        if ( config_mixed_layer_eddies_use_shortwave ) THEN
          shortwaveHeatFlux(:) = -config_mixed_layer_eddies_surface_cooling / (3.14159_RKIND * &
                                         config_mixed_layer_eddies_alpha)
        endif

        ! Set Coriolis parameters
        fCell(:) = config_mixed_layer_eddies_coriolis_parameter
        fEdge(:) = config_mixed_layer_eddies_coriolis_parameter
        fVertex(:) = config_mixed_layer_eddies_coriolis_parameter

        block_ptr => block_ptr % next
      end do

      deallocate(interfaceLocations)

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_mixed_layer_eddies!}}}

!***********************************************************************
!
!  routine ocn_init_validate_baroclinic_channel
!
!> \brief   Validation for baroclinic channel test case
!> \author  Doug Jacobsen
!> \date    02/20/2014
!> \details
!>  This routine validates the configuration options for the baroclinic channel test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_mixed_layer_eddies(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_mixed_layer_eddies_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('baroclinic_channel')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddies_vert_levels', config_mixed_layer_eddies_vert_levels)

      if(config_vert_levels <= 0 .and. config_mixed_layer_eddies_vert_levels > 0) then
         config_vert_levels = config_mixed_layer_eddies_vert_levels
      else if (config_vert_levels <= 0) then
        call mpas_log_write( 'Validation failed for mixed layer eddies test case. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_mixed_layer_eddies!}}}

!***********************************************************************

end module ocn_init_mixed_layer_eddies

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker

