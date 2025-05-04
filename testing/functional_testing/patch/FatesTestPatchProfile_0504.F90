program FatesTestPatch

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : itrue, nearzero
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader
  use FatesArgumentUtils,          only : command_line_arg
  use FatesCohortMod,              only : fates_cohort_type
  use FatesPatchMod,               only : fates_patch_type
  use FatesAllometryMod,           only : CrownDepth
  use FatesFactoryMod,             only : InitializeGlobals, GetSyntheticPatch
  use SyntheticPatchTypes,         only : synthetic_patch_array_type
  use EDCanopyStructureMod,        only : UpdatePatchLAI, UpdateCohortLAI
  use PRTParametersMod,            only : prt_params
  use EDParamsMod,                 only : nclmax, maxpft, nlevleaf, dinc_vai, dlower_vai
  use FatesRadiationMemMod,        only : num_rad_stream_types

  
  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader)             :: param_reader ! param reader instance
  type(synthetic_patch_array_type)               :: patch_data   ! array of synthetic patches
  character(len=:),                  allocatable :: param_file   ! input parameter file
  type(fates_patch_type),            pointer     :: patch        ! patch
  type(fates_cohort_type),           pointer     :: cohort       ! cohort
  integer                                        :: i            ! patch array location
  integer                                        :: iv           ! Vertical leaf layer index
  integer                                        :: cl           ! Canopy layer index
  integer                                        :: ft           ! Plant functional type index
  real(r8)                                       :: fleaf        ! fraction of cohort incepting area that is leaves.
  real(r8)                                       :: crown_depth  ! Current cohort's crown depth
  real(r8)                                       :: remainder    !Thickness of layer at bottom of canopy.
  real(r8)                                       :: fraction_exposed   ! how much of this layer is not covered by snow?
  real(r8)                                       :: elai_layer,tlai_layer    ! leaf area per canopy area
  real(r8)                                       :: esai_layer,tsai_layer    ! stem area per canopy area
  real(r8)                                       :: lai 
  real(r8)                                       :: sai 
  real(r8)                                       :: layer_top_height
  real(r8)                                       :: layer_bottom_height
  real(r8)                                       :: sum_tlai
  real(r8),                          allocatable :: avg_tlai_profile(:), height_bottom(:,:,:), height_top(:,:,:), patch_total_lai(:,:), patch_total_canopy_area(:,:) 
  integer                                        :: n, m, o, count
  logical                                        :: valid
  logical                                        :: re_allocate              ! Should we re-allocate the patch arrays?
  integer                                        :: prev_nveg                ! Previous number of vegetation layers
  integer                                        :: nveg                     ! Number of vegetation layers
  integer                                        :: ncan                     ! Number of canopy layers
  integer                                        :: prev_ncan                ! Number of canopy layers previously
  integer                                        :: npft  
    
  ! CONSTANTS:
  character(len=*), parameter :: out_file = 'patchprofile_out.nc'    ! output file
  integer,  parameter :: num_levsoil = 10      ! number of soil layers
  real(r8), parameter :: step_size = 1800.0_r8 ! step-size [s]
  logical, parameter  :: preserve_b4b = .true.


  interface

    subroutine WritePatchProfileData(out_file, numcl, numpft, numnleaf, tlai_profile, height_bottom, height_top, total_lai, dbh, idspft)
      use FatesConstantsMod,  only : r8 => fates_r8
      use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
      use FatesUnitTestIOMod, only : WriteVar
      use FatesUnitTestIOMod, only : RegisterVar
      use FatesUnitTestIOMod, only : EndNCDef
      use FatesUnitTestIOMod, only : type_double, type_int

      implicit none
      character(len=*), intent(in) :: out_file                 
      integer,          intent(in) :: numcl                  
      integer,          intent(in) :: numpft              
      integer,          intent(in) :: numnleaf                 
      real(r8),         intent(in) :: tlai_profile(:,:,:)  
      real(r8),         intent(in) :: height_bottom(:,:,:)
      real(r8),         intent(in) :: height_top(:,:,:)  
      real(r8),         intent(in) :: total_lai(:,:)  
      real(r8),         intent(in) :: dbh(:)
      integer,         intent(in)  :: idspft(:)

    end subroutine WritePatchProfileData

  end interface


  !read in parameter file name from command line
  param_file = command_line_arg(1)
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  
  ! initialize some global data we need
  call InitializeGlobals(step_size)
  
  ! get all the patch data
  call patch_data%GetSyntheticPatchData()
  
  i = patch_data%PatchDataPosition(patch_name='cohort8')
  call GetSyntheticPatch(patch_data%patches(i), num_levsoil, patch)


  ! ----------------------------------------------------------------------------------------
  ! calculate total canopy area for the patch -> canopy_summarization
  !write(*,*) patch%total_canopy_area,patch%total_tree_area
  patch%total_canopy_area = 0.0_r8
  patch%total_tree_area = 0.0_r8

  cohort => patch%shortest
  do while(associated(cohort))
      ft = cohort%pft
      if(cohort%canopy_layer==1)then
        patch%total_canopy_area = patch%total_canopy_area + cohort%c_area
        if( prt_params%woody(ft) == itrue)then
          patch%total_tree_area = patch%total_tree_area + cohort%c_area
        endif
      endif
      cohort => cohort%taller
  enddo 
  !write(*,*) patch%total_canopy_area,patch%total_tree_area

  ! ----------------------------------------------------------------------------------------
  ! calculate leaf layer and total leaf area index for the patch
  !write(*,*) patch%nleaf,patch%canopy_layer_tlai
  patch%nleaf(:,:) = 0
  patch%canopy_layer_tlai(:) = 0.0_r8
  
  call UpdatePatchLAI(patch)
  !write(*,*)  size(patch%nleaf, 1), size(patch%nleaf, 2)!2 layers, 16 pft
  !write(*,*) patch%nleaf,patch%canopy_layer_tlai, patch%total_canopy_area

  !allocate(patch%elai_profile(nclmax,maxpft,nlevleaf)) 
  !allocate(patch%tlai_profile(nclmax,maxpft,nlevleaf)) 
  !allocate(patch%esai_profile(nclmax,maxpft,nlevleaf)) 
  !allocate(patch%tsai_profile(nclmax,maxpft,nlevleaf)) 
  !allocate(patch%canopy_area_profile(nclmax,maxpft,nlevleaf)) 

  !write(*,*) patch%ncl_p,patch%nleaf(:,:),maxpft
  patch%ncl_p = 4 !2
  npft = maxpft
  !write(*,*) patch_data%patches(i)%dbhs(:), patch_data%patches(i)%pft_ids(:)
  !!! call patch%ReAllocateDynamics()
  ncan = patch%ncl_p
  nveg = maxval(patch%nleaf(:,:))
  re_allocate = .true.
  if(re_allocate) then
    nveg = nveg + 1 
    allocate(patch%tlai_profile(ncan,npft,nveg))
    allocate(patch%tsai_profile(ncan,npft,nveg))
    allocate(patch%elai_profile(ncan,npft,nveg))
    allocate(patch%esai_profile(ncan,npft,nveg))
    allocate(patch%canopy_area_profile(ncan,npft,nveg))
    allocate(patch%f_sun(ncan,npft,nveg))
    allocate(patch%fabd_sun_z(ncan,npft,nveg))
    allocate(patch%fabd_sha_z(ncan,npft,nveg))
    allocate(patch%fabi_sun_z(ncan,npft,nveg))
    allocate(patch%fabi_sha_z(ncan,npft,nveg))
    allocate(patch%nrmlzd_parprof_pft_dir_z(num_rad_stream_types,ncan,npft,nveg))
    allocate(patch%nrmlzd_parprof_pft_dif_z(num_rad_stream_types,ncan,npft,nveg))
    allocate(patch%ed_parsun_z(ncan,npft,nveg))
    allocate(patch%ed_parsha_z(ncan,npft,nveg))
    allocate(patch%ed_laisun_z(ncan,npft,nveg))
    allocate(patch%ed_laisha_z(ncan,npft,nveg))
    allocate(patch%parprof_pft_dir_z(ncan,npft,nveg))
    allocate(patch%parprof_pft_dif_z(ncan,npft,nveg))
    allocate(height_bottom(ncan,npft,nveg))
    allocate(height_top(ncan,npft,nveg))
    allocate(patch_total_lai(ncan,nveg))
    allocate(patch_total_canopy_area(ncan,nveg))
  end if
  !!! call patch%ReAllocateDynamics()
  call patch%NanDynamics()
  call patch%ZeroDynamics()


  if_any_canopy_area: if (patch%total_canopy_area > nearzero ) then
    cohort => patch%shortest
    do while(associated(cohort))
      ft = cohort%pft
      cl = cohort%canopy_layer
      if_preserve_b4b: if(preserve_b4b) then
        lai = cohort%treelai * cohort%c_area/patch%total_canopy_area
        sai = cohort%treesai * cohort%c_area/patch%total_canopy_area
        if( (cohort%treelai+cohort%treesai) > nearzero)then
          fleaf = lai / (lai+sai)
        else
          fleaf = 0._r8
        endif

        call CrownDepth(cohort%height,cohort%pft,crown_depth)

        write(*,*) cohort%height,cohort%dbh,crown_depth
        do iv = 1,cohort%NV
          layer_top_height = cohort%height - ( real(iv-1,r8)/cohort%NV * crown_depth )
          layer_bottom_height = cohort%height - ( real(iv,r8)/cohort%NV * crown_depth )
          !!write(*,*) 'layer_top_height', layer_top_height,iv
          write(*,*) 'layer_bottom_height', layer_bottom_height
          height_bottom(cl,ft,iv) = layer_bottom_height
          height_top(cl,ft,iv)    = layer_top_height


          fraction_exposed = 1.0_r8 !!! snow depth of the site is unknown, set to 1.0
          if(iv==cohort%NV) then
            remainder = (cohort%treelai + cohort%treesai) - (dlower_vai(iv) - dinc_vai(iv))
            if(remainder > dinc_vai(iv) )then
                write(*, *)'ED: issue with remainder', &
                    cohort%treelai,cohort%treesai,dinc_vai(iv), & 
                    cohort%NV,remainder
            endif
          else
            remainder = dinc_vai(iv)
          end if
          
          patch%tlai_profile(cl,ft,iv) = patch%tlai_profile(cl,ft,iv) + &
              remainder * fleaf * cohort%c_area/patch%total_canopy_area

          patch_total_lai(cl,iv) = patch_total_lai(cl,iv) + &
              remainder * fleaf * cohort%c_area/patch%total_canopy_area
          
          patch%elai_profile(cl,ft,iv) = patch%elai_profile(cl,ft,iv) + &
              remainder * fleaf * cohort%c_area/patch%total_canopy_area * &
              fraction_exposed

          patch%tsai_profile(cl,ft,iv) = patch%tsai_profile(cl,ft,iv) + &
              remainder * (1._r8 - fleaf) * cohort%c_area/patch%total_canopy_area

          patch%esai_profile(cl,ft,iv) = patch%esai_profile(cl,ft,iv) + &
              remainder * (1._r8 - fleaf) * cohort%c_area/patch%total_canopy_area * &
              fraction_exposed

          patch%canopy_area_profile(cl,ft,iv) = patch%canopy_area_profile(cl,ft,iv) + &
              cohort%c_area/patch%total_canopy_area

          patch_total_canopy_area(cl,iv) = patch_total_canopy_area(cl,iv) + &
              cohort%c_area/patch%total_canopy_area
          
        end do

      end if if_preserve_b4b

      cohort => cohort%taller
    enddo !cohort

    
    do cl = 1,patch%NCL_p
      do ft = 1,npft
        do iv = 1,patch%nleaf(cl,ft)
          if(sum(patch%canopy_area_profile(cl,:,iv)) > 1.0001_r8 ) then
            write(*,*) 'ERROR1'

            cohort => patch%shortest
            do while(associated(cohort))
              if (cohort%canopy_layer==cl)then
                write(*,*) 'ERROR2'
              end if
              cohort => cohort%taller
            end do

          end if
        end do
      end do
    end do

    do cl = 1,patch%NCL_p    
      do ft = 1,npft
        do iv = 1,patch%nleaf(cl,ft)
          if( patch%canopy_area_profile(cl,ft,iv) > nearzero ) then
              !write(*,*) cl, ft, iv
              !write(*,*) patch%tlai_profile(cl,ft,iv), patch%canopy_area_profile(cl,ft,iv)
              patch%tlai_profile(cl,ft,iv) = patch%tlai_profile(cl,ft,iv) / &
                    patch%canopy_area_profile(cl,ft,iv)
              
              patch_total_lai(cl,iv) = patch_total_lai(cl,iv) / &
                    patch_total_canopy_area(cl,iv)

              patch%tsai_profile(cl,ft,iv) = patch%tsai_profile(cl,ft,iv) / &
                    patch%canopy_area_profile(cl,ft,iv)

              patch%elai_profile(cl,ft,iv) = patch%elai_profile(cl,ft,iv) / &
                    patch%canopy_area_profile(cl,ft,iv)

              patch%esai_profile(cl,ft,iv) = patch%esai_profile(cl,ft,iv) / &
                    patch%canopy_area_profile(cl,ft,iv)
          end if
        enddo
      enddo
    enddo

  end if if_any_canopy_area

  !write(*,*)  size(patch%tlai_profile, 1), size(patch%tlai_profile, 2), size(patch%tlai_profile, 3)
  !write(*,*)  nveg, patch%nleaf
  write(*,*) patch_total_lai
  ! write out data to netcdf file
  call WritePatchProfileData(out_file, patch%NCL_p , npft, nveg, patch%tlai_profile, height_bottom, height_top, patch_total_lai, patch_data%patches(i)%dbhs(:), patch_data%patches(i)%pft_ids(:))

  ! ----------------------------------------------------------------------------------------
  ! Calculate the total LAI profile for the patch, may be wrong? not average but sum
  count = 0
  do cl = 1,size(patch%tlai_profile, 1)
    do iv = 1,size(patch%tlai_profile, 3)
      if (sum(patch%tlai_profile(cl,:,iv)) > nearzero) then
        count = count + 1
      end if
    end do
  end do
  write(*,*) count


  allocate(avg_tlai_profile(count))
  avg_tlai_profile = 0.0_r8
  count = 0
  do cl = 1,size(patch%tlai_profile, 1)
    do iv = 1, size(patch%tlai_profile, 3)

      o = 0
      sum_tlai = 0.0_r8
      do ft = 1, maxpft
        if (patch%tlai_profile(cl,ft,iv) > nearzero) then
          sum_tlai = sum_tlai + patch%tlai_profile(cl,ft,iv)
          o = o + 1
        end if
      end do

      if (o > 0) then
        count = count + 1
        avg_tlai_profile(count) = sum_tlai/o     
      end if

    end do
  end do

  !write(*,*) avg_tlai_profile
end program FatesTestPatch

! ----------------------------------------------------------------------------------------

subroutine WritePatchProfileData(out_file, numcl, numpft, numnleaf, tlai_profile, height_bottom, height_top, total_lai, dbh, idspft)
  !
  ! DESCRIPTION:
  ! Writes out data from the patch profile test
  !
  use FatesConstantsMod,  only : r8 => fates_r8
  use FatesUnitTestIOMod, only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod, only : WriteVar
  use FatesUnitTestIOMod, only : RegisterVar
  use FatesUnitTestIOMod, only : EndNCDef
  use FatesUnitTestIOMod, only : type_double, type_int

  implicit none
  ! ARGUMENTS:
  character(len=*), intent(in) :: out_file                 ! output file name
  integer,          intent(in) :: numcl                    ! number of canopy layers
  integer,          intent(in) :: numpft                   ! number of pfts
  integer,          intent(in) :: numnleaf                 ! number of leaf layers
  !real(r8),         intent(in) :: avg_tlai_profile(:)     ! lai profile [m2/m2] on non-zero layers
  real(r8),         intent(in) :: tlai_profile(:,:,:)      ! lai profile [m2/m2]
  real(r8),         intent(in) :: height_bottom(:,:,:)    ! bottom height of leaf layer [m]
  real(r8),         intent(in) :: height_top(:,:,:)       ! top height of leaf layer [m]
  real(r8),         intent(in) :: total_lai(:,:)          ! total lai profile [m2/m2]
  real(r8),         intent(in) :: dbh(:)                   ! dbh [cm]
  integer,         intent(in)  :: idspft(:)                 ! pft index

  ! LOCALS:
  integer, allocatable :: pft_indices(:) ! array of pft indices to write out
  integer, allocatable :: cl_indices(:)  ! array of canopy layer indices to write out
  integer, allocatable :: nleaf_indices(:) ! array of leaf layer indices to write out
  integer              :: i              ! looping index
  integer              :: ncid           ! netcdf file id
  character(len=8)     :: dim_names(3)   ! dimension names
  integer              :: dimIDs(3)      ! dimension IDs
  integer              :: clID, pftID, nleafID   ! variable IDs for dimensions
  integer              :: tlaiID, totallaiID     ! variable ID for tlai
  integer              :: height_bottomID, height_topID ! variable ID for height
  integer              :: dbhID, idspftID        ! variable ID for tlai

  ! create pft indices
  allocate(pft_indices(numpft))
  do i = 1, numpft
    pft_indices(i) = i
  end do

  ! create cl indices
  allocate(cl_indices(numcl))
  do i = 1, numcl
    cl_indices(i) = i
  end do

  ! create nleaf indices
  allocate(nleaf_indices(numnleaf))
  do i = 1, numnleaf
    nleaf_indices(i) = i
  end do

  ! dimension names
  dim_names = [character(len=12) :: 'canopylayer', 'pft', 'nleaf']

  ! open file
  call OpenNCFile(trim(out_file), ncid, 'readwrite')

  ! register dimensions
  call RegisterNCDims(ncid, dim_names, (/numcl, numpft, numnleaf/), 3, dimIDs)

  ! register canopylayer
  call RegisterVar(ncid, dim_names(1), dimIDs(1:1), type_int,         &
    [character(len=20)  :: 'units', 'long_name'],                        &
    [character(len=150) :: '', 'canopylayer'], 2, clID)

  ! register pft
  call RegisterVar(ncid, dim_names(2), dimIDs(2:2), type_int,      &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'plant functional type'], 2, pftID)

  ! register nleaf
  call RegisterVar(ncid, dim_names(3), dimIDs(3:3), type_int,      &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'number of leaf layers'], 2, nleafID)

  ! register dbh
  call RegisterVar(ncid, 'dbh', dimIDs(2:2), type_double,      &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: 'cm', 'diameter at breast height'], 2, dbhID)

  ! register pft index
  call RegisterVar(ncid, 'idspft', dimIDs(2:2), type_int,      &
    [character(len=20)  :: 'units', 'long_name'],                  &
    [character(len=150) :: '', 'index of plant functional type'], 2, idspftID)

  ! register tlai_profile
  call RegisterVar(ncid, 'tlai_profile', dimIDs(1:3), type_double, &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],            &
    [character(len=150) :: 'cl pft nleaf', 'm2/m2', 'total leaf area index profile'], &
    3, tlaiID)

  ! register height_bottom
  call RegisterVar(ncid, 'height_bottom', dimIDs(1:3), type_double, &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],            &
    [character(len=150) :: 'cl pft nleaf', 'm', 'bottom height of leaf layer'], &
    3, height_bottomID)

  ! register height_top
  call RegisterVar(ncid, 'height_top', dimIDs(1:3), type_double, &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],            &
    [character(len=150) :: 'cl pft nleaf', 'm', 'top height of leaf layer'], &
    3, height_topID)
  
  ! reguster total_lai
  call RegisterVar(ncid, 'total_lai', [dimIDs(1), dimIDs(3)], type_double, &
    [character(len=20)  :: 'coordinates', 'units', 'long_name'],            &
    [character(len=150) :: 'cl nleaf', 'm2/m2', 'total leaf area index profile (no pft)'], &
    2, totallaiID)

  ! finish defining variables
  call EndNCDef(ncid)

  ! write out data
  call WriteVar(ncid, clID, cl_indices(:))
  call WriteVar(ncid, pftID, pft_indices(:))
  call WriteVar(ncid, nleafID, nleaf_indices(:))
  call WriteVar(ncid, dbhID, dbh(:))
  call WriteVar(ncid, idspftID, idspft(:))
  call WriteVar(ncid, tlaiID, tlai_profile(:,:,:))
  call WriteVar(ncid, height_bottomID, height_bottom(:,:,:))
  call WriteVar(ncid, height_topID, height_top(:,:,:))
  call WriteVar(ncid, totallaiID, total_lai(:,:))

  ! close the file
  call CloseNCFile(ncid)

end subroutine WritePatchProfileData
