program FatesTestPatch

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : itrue, ifalse, g_per_kg 
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use FatesArgumentUtils,          only : command_line_arg
  use FatesCohortMod,              only : fates_cohort_type
  use FatesPatchMod,               only : fates_patch_type
  use FatesFactoryMod,             only : InitializeGlobals, GetSyntheticPatch
  use SyntheticPatchTypes,         only : synthetic_patch_array_type
  use EDCanopyStructureMod,        only : UpdatePatchLAI, UpdateCohortLAI
  use PRTParametersMod,            only : prt_params
  use EDParamsMod,                 only : nclmax, dinc_vai, dlower_vai
  use PRTGenericMod,               only : leaf_organ, carbon12_element
  use FatesInterfaceTypesMod,      only : hlm_use_sp !!  Flag to use FATES satellite phenology (LAI) mode 1 = TRUE, 0 = FALSE
  use FatesAllometryMod,           only : tree_lai_sai, decay_coeff_vcmax, bleaf, tree_lai


  
  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader)             :: param_reader ! param reader instance
  type(synthetic_patch_array_type)               :: patch_data   ! array of synthetic patches
  character(len=:),                  allocatable :: param_file   ! input parameter file
  type(fates_patch_type),            pointer     :: patch        ! patch
  type(fates_cohort_type),           pointer     :: cohort       ! cohort
  integer                                        :: i            ! patch array location
  integer                                        :: cl           ! Canopy layer index
  integer                                        :: ft           ! Plant functional type index
  real(r8)                                       :: leaf_c       ! leaf carbon [kg]
  real(r8)                                       :: treesai      ! stem area index within crown m2/m2
  real(r8)                                       :: treelai           ! the in-crown leaf area index for the plant [m2 leaf/m2 crown footprint]
  real(r8)                                       :: leafc_per_unitarea ! KgC of leaf per m2 area of ground.
  real(r8)                                       :: slat               ! the sla of the top leaf layer. m2/kgC
  real(r8)                                       :: canopy_lai_above   ! total LAI of canopy layer overlying this tree
  real(r8)                                       :: kn                 ! coefficient for exponential decay of 1/sla and vcmax with canopy depth
  real(r8)                                       :: sla_max            ! Observational constraint on how large sla (m2/gC) can become
  real(r8)                                       :: leafc_slamax       ! Leafc_per_unitarea at which sla_max is reached
  real(r8)                                       :: clim               ! Upper limit for leafc_per_unitarea in exponential tree_lai function
  real(r8)                                       :: target_lai
  real(r8)                                       :: target_bleaf


  

  ! CONSTANTS:
  integer,  parameter :: num_levsoil = 10      ! number of soil layers
  real(r8), parameter :: step_size = 1800.0_r8 ! step-size [s]

  !read in parameter file name from command line
  param_file = command_line_arg(1)
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  
  ! initialize some global data we need
  call InitializeGlobals(step_size)
  
  ! get all the patch data
  call patch_data%GetSyntheticPatchData()
  
  i = patch_data%PatchDataPosition(patch_name='3layers')
  call GetSyntheticPatch(patch_data%patches(i), num_levsoil, patch)
  
  ! update the patch LAI
  !! calculate total_canopy_area
  cohort => patch%shortest
  patch%total_canopy_area = 0.0_r8
  do while (associated(cohort))
  if (Cohort%canopy_layer==1)then
      patch%total_canopy_area = patch%total_canopy_area + cohort%c_area
      if( prt_params%woody(cohort%pft) == 1)then
            patch%total_tree_area = patch%total_tree_area + cohort%c_area
      endif
  endif
  cohort => cohort%taller
  end do

  !! update the cohort LAI
  patch%canopy_layer_tlai(:) = 0.0_r8
  do cl = 1, nclmax
    cohort => patch%tallest
    do while (associated(cohort))
  
      if (cohort%canopy_layer .eq. cl) then
        ft = cohort%pft
        !call UpdateCohortLAI(cohort, patch%canopy_layer_tlai(cl), patch%total_canopy_area)
        leaf_c = cohort%prt%GetState(leaf_organ,carbon12_element)
        !write(*,*) prt_params%allom_lmode(ft)

  
        !call  tree_lai_sai(leaf_c, cohort%pft, cohort%c_area, cohort%n,           &
        !  cohort%canopy_layer, patch%canopy_layer_tlai, cohort%vcmax25top, cohort%dbh, cohort%crowndamage,          &
        !  cohort%canopy_trim, cohort%efstem_coh, 4, cohort%treelai, treesai )
        !!! treelai calculation
        slat = g_per_kg * prt_params%slatop(ft) ! m2/g to m2/kg
        leafc_per_unitarea = leaf_c/(cohort%c_area/cohort%n) !KgC/m2
        
        !write(*,*) leafc_per_unitarea
        if(leafc_per_unitarea > 0.0_r8)then
          if (cohort%canopy_layer==1) then ! if in we are in the canopy (top) layer)
            canopy_lai_above = 0._r8
          else
            canopy_lai_above = sum(patch%canopy_layer_tlai(1:cl-1))
          end if
          ! Coefficient for exponential decay of 1/sla with canopy depth:
          kn = decay_coeff_vcmax(cohort%vcmax25top, &
                              prt_params%leafn_vert_scaler_coeff1(ft), &
                              prt_params%leafn_vert_scaler_coeff2(ft))
          !write(*,*)ft, leaf_c
          ! take PFT-level maximum SLA value, even if under a thick canopy (which has units of m2/gC),
          ! and put into units of m2/kgC
          sla_max = g_per_kg*prt_params%slamax(ft)
          ! Leafc_per_unitarea at which sla_max is reached due to exponential sla profile in canopy:
          leafc_slamax = (slat - sla_max * exp(-1.0_r8 * kn * canopy_lai_above)) / &
                (-1.0_r8 * kn * slat * sla_max)
          if(leafc_slamax < 0.0_r8)then
              leafc_slamax = 0.0_r8
          endif
          
          if (leafc_per_unitarea <= leafc_slamax)then
            treelai = (log(exp(-1.0_r8 * kn * canopy_lai_above) - &
              kn * slat * leafc_per_unitarea) + &
              (kn * canopy_lai_above)) / (-1.0_r8 * kn)
            clim = (exp(-1.0_r8 * kn * canopy_lai_above)) / (kn * slat)
            !write(*,*) 1, treelai
          else if(leafc_per_unitarea > leafc_slamax)then
            treelai = ((log(exp(-1.0_r8 * kn * canopy_lai_above) - kn * slat * leafc_slamax) + &
              (kn * canopy_lai_above)) / (-1.0_r8 * kn)) + &
              (leafc_per_unitarea - leafc_slamax) * sla_max
            clim = (exp(-1.0_r8 * kn * canopy_lai_above)) / (kn * slat)
            !write(*,*) 2, treelai
          end if
        else
          treelai = 0.0_r8
        end if 
        
        !!! treesai calculation
        call bleaf(cohort%dbh, ft, cohort%crowndamage, cohort%canopy_trim, 1.0_r8, target_bleaf)
        target_lai = tree_lai(target_bleaf, ft, cohort%c_area, cohort%n, cl, patch%canopy_layer_tlai, cohort%vcmax25top) 
        treesai    = cohort%efstem_coh * prt_params%allom_sai_scaler(ft) * target_lai
        !write(*,*) target_lai, treelai
        !write(*,*) cohort%efstem_coh, prt_params%allom_sai_scaler(ft), target_lai

        cohort%treelai = treelai
        if (hlm_use_sp .eq. ifalse) then !hlm_use_sp = ifalse
            cohort%treesai = treesai
        end if
   
        ! Number of actual vegetation layers in this cohort's crown
        cohort%nv =  count((cohort%treelai+cohort%treesai) .gt. dlower_vai(:)) + 1
        !write(*,*) cohort%treelai, cohort%treesai
        !write(*,*) cohort%treelai, cohort%treesai, cohort%nv

        patch%nleaf(cl,ft) = max(patch%nleaf(cl,ft),cohort%NV)
        patch%canopy_layer_tlai(cl) = patch%canopy_layer_tlai(cl) + cohort%treelai *cohort%c_area/patch%total_canopy_area
      end if   

    cohort => cohort%shorter
    end do
  end do


  !write(*,*) patch%canopy_layer_tlai!patch%total_canopy_area,patch%canopy_layer_tlai(1),patch%canopy_layer_tlai(2)

  ! print out list in ascending order
  cohort => patch%shortest
  write(*,*) 'Updated Patch LAI: '
  do while (associated(cohort))
    !write(*,*) cohort%pft,cohort%height,cohort%treelai,cohort%c_area,cohort%canopy_layer
    cohort => cohort%taller
  end do
  write(*,*) ' '

end program FatesTestPatch

