!==========================================================================================!
!==========================================================================================!
! This subroutine reads ED-2.0 history files (text files), with  more allometric info.
! You can set allometric variables to 0, in which case they will be auto-calculated.
!------------------------------------------------------------------------------------------!
subroutine read_ext_ed20_history_file


   use ed_max_dims         , only : n_pft                       & ! intent(in)
                                  , huge_patch                  & ! intent(in)
                                  , huge_cohort                 & ! intent(in)
                                  , max_water                   & ! intent(in)
                                  , str_len                     & ! intent(in)
                                  , maxfiles                    & ! intent(in)
                                  , maxlist                     ! ! intent(in)
   use pft_coms            , only : q                           & ! intent(in)
                                  , qsw                         & ! intent(in)
                                  , qbark                       & ! intent(in)
                                  , SLA                         & ! intent(in)
                                  , min_dbh                     & ! intent(in)
                                  , min_bdead                   & ! intent(in)
                                  , is_grass                    & ! intent(in)
                                  , include_pft                 & ! intent(in)
                                  , include_pft_ag              & ! intent(in)
                                  , pft_1st_check               & ! intent(in)
                                  , agf_bs                      & ! intent(in)
                                  , f_bstorage_init             & ! intent(in)
                                  , include_these_pft           ! ! intent(in)
   use ed_misc_coms        , only : sfilin                      & ! intent(in)
                                  , ied_init_mode               ! ! intent(in)
   use consts_coms         , only : pio180                      & ! intent(in)
                                  , pio4                        & ! intent(in)
                                  , almost_zero                 ! ! intent(in)
   use ed_misc_coms        , only : use_target_year             & ! intent(in)
                                  , restart_target_year         ! ! intent(in)
   use ed_state_vars       , only : polygontype                 & ! variable type
                                  , sitetype                    & ! variable type
                                  , patchtype                   & ! variable type
                                  , edtype                      & ! variable type
                                  , edgrid_g                    & ! variable type
                                  , allocate_sitetype           & ! subroutine
                                  , allocate_patchtype          ! ! subroutine
   use grid_coms           , only : ngrids                      & ! intent(in)
                                  , nzg                         ! ! intent(in)
   use allometry           , only : bd2dbh                      & ! function
                                  , dbh2h                       & ! function
                                  , size2bd                     & ! function
                                  , size2bl                     & ! function
                                  , size2bt                     & ! function
                                  , size2xb                     & ! function
                                  , ed_balive                   & ! function
                                  , ed_biomass                  & ! function
                                  , area_indices                ! ! subroutine
   use fuse_fiss_utils     , only : sort_cohorts                & ! subroutine
                                  , sort_patches                ! ! subroutine
   use disturb_coms        , only : ianth_disturb               ! ! intent(in)
   use decomp_coms         , only : decomp_scheme               & ! intent(in)
                                  , agf_fsc                     & ! intent(in)
                                  , agf_stsc                    & ! intent(in)
                                  , f0_msc                      & ! intent(in)
                                  , f0_ssc                      & ! intent(in)
                                  , f0_psc                      & ! intent(in)
                                  , c2n_structural              ! ! intent(in)
   use physiology_coms     , only : iddmort_scheme              & ! intent(in)
                                  , trait_plasticity_scheme     ! ! intent(in)
   use update_derived_utils, only : update_cohort_plastic_trait ! ! subroutine
   use ed_type_init        , only : init_ed_cohort_vars         & ! subroutine
                                  , init_ed_patch_vars          & ! subroutine
                                  , init_ed_site_vars           & ! subroutine
                                  , init_ed_poly_vars           ! ! subroutine
   implicit none

   !----- Local constants. ----------------------------------------------------------------!
   real(kind=8), parameter :: min_area = 1.d-7           ! Minimum acceptable area.
   real(kind=8), parameter :: min_ok   = 1.d-20          ! Minimum acceptable value for
                                                         !    any restart variable.
   !----- Local variables. ----------------------------------------------------------------!
   type(edtype)          , pointer                        :: cgrid
   type(polygontype)     , pointer                        :: cpoly
   type(sitetype)        , pointer                        :: csite
   type(patchtype)       , pointer                        :: cpatch
   character(len=str_len), dimension(maxlist)             :: full_list
   character(len=str_len), dimension(maxfiles)            :: site_list
   character(len=str_len), dimension(maxfiles)            :: pss_list
   character(len=str_len), dimension(maxfiles)            :: css_list
   character(len=str_len), dimension(huge_patch)          :: pname
   character(len=str_len)                                 :: pss_name
   character(len=str_len)                                 :: css_name
   character(len=str_len)                                 :: cdum
   character(len=str_len), dimension(huge_cohort)         :: cname
   character(len=str_len), dimension(huge_cohort)         :: cpname
   integer               , dimension(huge_patch)          :: trk
   integer               , dimension(huge_patch)          :: sitenum
   integer               , dimension(huge_cohort)         :: leaves_on
   integer               , dimension(huge_cohort)         :: ipft
   integer                                                :: year
   integer                                                :: igr
   integer                                                :: ipy
   integer                                                :: isi
   integer                                                :: ipa
   integer                                                :: ico
   integer                                                :: ip
   integer                                                :: ip2
   integer                                                :: ic
   integer                                                :: ic2
   integer                                                :: nwater
   integer                                                :: ierr
   integer                                                :: nf
   integer                                                :: nflist
   integer                                                :: nflsite
   integer                                                :: nflpss
   integer                                                :: nflcss
   integer                                                :: nclosest
   integer                                                :: ncohorts
   integer                                                :: npatchco
   integer                                                :: npatches
   integer                                                :: nsitepat
   integer                                                :: npatch2
   integer                                                :: nw
   logical              , dimension(huge_cohort)          :: add_this_cohort
   logical                                                :: renumber_pfts
   logical                                                :: site_match
   real                 , dimension(max_water)            :: depth
   real                 , dimension(huge_patch)           :: time
   real                 , dimension(huge_patch)           :: age
   real                 , dimension(huge_patch)           :: area
   real                 , dimension(huge_patch)           :: fsc
   real                 , dimension(huge_patch)           :: stsc
   real                 , dimension(huge_patch)           :: stsl
   real                 , dimension(huge_patch)           :: ssc
   real                 , dimension(huge_patch)           :: msn
   real                 , dimension(huge_patch)           :: fsn
   real                 , dimension(max_water,huge_patch) :: water
   real                 , dimension(12,huge_cohort)       :: cb
   real                 , dimension(12,huge_cohort)       :: cb_max
   real                 , dimension(huge_cohort)          :: balive
   real                 , dimension(huge_cohort)          :: avgRg
   real                 , dimension(huge_cohort)          :: nplant
   real                 , dimension(huge_cohort)          :: bleaf
   real                 , dimension(huge_cohort)          :: broot
   real                 , dimension(huge_cohort)          :: bsapwooda
   real                 , dimension(huge_cohort)          :: bsapwoodb
   real                 , dimension(huge_cohort)          :: bbarka
   real                 , dimension(huge_cohort)          :: bbarkb
   real                 , dimension(huge_cohort)          :: bdeada
   real                 , dimension(huge_cohort)          :: bdeadb
   real                 , dimension(huge_cohort)          :: hite
   real                 , dimension(huge_cohort)          :: dbh
   real                 , dimension(huge_cohort)          :: ctime
   real                 , dimension(maxfiles)             :: slon_list,slat_list
   real                 , dimension(maxfiles)             :: plon_list,plat_list
   real                 , dimension(maxfiles)             :: clon_list,clat_list
   real                 , dimension(maxfiles)             :: file_pdist,file_cdist
   real                                                   :: dummy
   real                                                   :: area_tot
   real                                                   :: area_sum
   real(kind=8)         , dimension(max_water)            :: dwater
   real(kind=8)                                           :: dage
   real(kind=8)                                           :: darea
   real(kind=8)                                           :: dfsc
   real(kind=8)                                           :: dstsc
   real(kind=8)                                           :: dstsl
   real(kind=8)                                           :: dssc
   real(kind=8)                                           :: dpsc
   real(kind=8)                                           :: dmsn
   real(kind=8)                                           :: dfsn
   real(kind=8)                                           :: bdead
   !----- External function. --------------------------------------------------------------!
   real                 , external                        :: sngloff
   real                 , external                        :: dist_gc
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !      Now we loop over all all grids and polygons, and fill them with patches and      !
   ! cohorts from the closest polygon.                                                     !
   !---------------------------------------------------------------------------------------!
   gridloop: do igr = 1,ngrids

      !----- Retrieve all files with the specified prefix. --------------------------------!
      call ed_filelist(full_list,sfilin(igr),nflist)

      !----- Retrieve LON/LAT information for sites ---------------------------------------!
      renumber_pfts = .false.

      !----- Retrieve LON/LAT information for patches and cohorts -------------------------!
      call ed1_fileinfo('.pss',nflist,full_list,nflpss,pss_list,plon_list,plat_list)
      call ed1_fileinfo('.css',nflist,full_list,nflcss,css_list,clon_list,clat_list)

      cgrid => edgrid_g(igr)

      polyloop: do ipy = 1,cgrid%npolygons

         cpoly => cgrid%polygon(ipy)


         !---------------------------------------------------------------------------------!
         !    Initialise the distances as very large numbers, so if we don't fill all the  !
         ! patches and cohorts, we will not going to take non-sense as a valid polygon.    !
         !---------------------------------------------------------------------------------!
         file_pdist(:) = 1.e20
         file_cdist(:) = 1.e20

         !---------------------------------------------------------------------------------!
         !    Compute the distances between every polygon in the restart files and the     !
         ! current polygon.                                                                !
         !---------------------------------------------------------------------------------!
         do nf=1,nflpss
            file_pdist(nf) = dist_gc(cgrid%lon(ipy),plon_list(nf)                          &
                                    ,cgrid%lat(ipy),plat_list(nf) )
         end do
         do nf=1,nflcss
            file_cdist(nf) = dist_gc(cgrid%lon(ipy),clon_list(nf)                          &
                                    ,cgrid%lat(ipy),clat_list(nf) )
         end do
         !---------------------------------------------------------------------------------!


         find_nonwater: do nf=1,nflpss
            !------------------------------------------------------------------------------!
            !     Find the file that is the closest to the current polygon, based on the   !
            ! distance vector.                                                             !
            !------------------------------------------------------------------------------!
            nclosest = minloc(file_pdist,dim=1)
            pss_name = trim(pss_list(nclosest))
            write (unit=*,fmt='(2a)') 'Using patch file: ',trim(pss_name)

            !------------------------------------------------------------------------------!
            !    Open the patch file and read in all patches.                              !
            !------------------------------------------------------------------------------!
            open(unit=12,file=trim(pss_name),form='formatted',status='old',action='read')
            read(unit=12,fmt='(a4)')  cdum

            !----- Read the other information from the header (if there is any...). -------!
            nwater = 1
            if (ied_init_mode == 1 ) then
               read (unit=12,fmt=*) cdum,nwater
               read (unit=12,fmt=*) cdum,depth(1:nwater)
               read (unit=12,fmt=*)
            end if

            !------------------------------------------------------------------------------!
            !     Now we loop over all patches and decide whether they should be included  !
            ! or not.                                                                      !
            !------------------------------------------------------------------------------!
            ip      = 1
            sitenum = 0
            count_patches: do

               !---------------------------------------------------------------------------!
               !     We must check whether we are not exceeding the maximum number of      !
               ! patches that we can read.                                                 !
               !---------------------------------------------------------------------------!
               if (ip > huge_patch) then
                  write (unit=*,fmt='(a,1x,a)')  ' In file:',trim(pss_name)
                  write (unit=*,fmt='(a)')       ' Number of patches is > HUGE_PATCH...'
                  write (unit=*,fmt='(a,1x,i7)') ' HUGE_PATCH:',huge_patch
                  write (unit=*,fmt='(a)')       ' Increase HUGE_PATCH to read this...'
                  call fatal_error('Too many patches to be read...'                        &
                                  ,'read_ext_ed20_history_file','ed_history_io.f90')
               end if

               !---------------------------------------------------------------------------!
               !    If we can still add new patches, read the next one, according to the   !
               ! input format type.                                                        !
               !---------------------------------------------------------------------------!
               !----- Standard ED-2.0 file. --------------------------------------------!
               read(unit=12,fmt=*,iostat=ierr) time(ip),pname(ip),trk(ip),dage,darea    &
                                              ,dwater(1),dfsc,dstsc,dstsl,dssc,dummy    &
                                              ,dmsn,dfsn

               !------------------------------------------------------------------------!
               !     Check whether the file has hit the end, and if so, leave the loop. !
               !------------------------------------------------------------------------!
               if(ierr /= 0)exit count_patches

               !----- Copy the double-precision scratch variables to the arrays. -------!
               area   (ip) = sngloff(darea    ,min_area)
               age    (ip) = sngloff(dage     ,min_ok  )
               fsc    (ip) = sngloff(dfsc     ,min_ok  )
               stsc   (ip) = sngloff(dstsc    ,min_ok  )
               stsl   (ip) = sngloff(dstsl    ,min_ok  )
               ssc    (ip) = sngloff(dssc     ,min_ok  )
               msn    (ip) = sngloff(dmsn     ,min_ok  )
               fsn    (ip) = sngloff(dfsn     ,min_ok  )
               water(1,ip) = sngloff(dwater(1),min_ok  )

               !---------------------------------------------------------------------------!



               !---------------------------------------------------------------------------!
               !      By adding ip only if the area is above minimum, we avoid including   !
               ! patches that are tiny since they will be overwritten by the next patch.   !
               !---------------------------------------------------------------------------!
               if (area(ip) > min_area) ip = ip + 1
               !---------------------------------------------------------------------------!

            end do count_patches
            !------------------------------------------------------------------------------!



            !------------------------------------------------------------------------------!
            !      Here we determine the number of patches.  Unless there are no patches,  !
            ! (this may happen if the file used was all water.  In that case, then we      !
            ! should use the next closest file.                                            !
            !------------------------------------------------------------------------------!
            npatches = max(ip-1,0)

            close(unit=12,status='keep')

            if (npatches > 0) then

               !---------------------------------------------------------------------------!
               !      We have found a suitable file, we will leave the loop, after we find !
               ! the name of the corresponding cohort file.                                !
               !---------------------------------------------------------------------------!
               nclosest = minloc(abs(file_pdist(nclosest)-file_cdist),dim=1)
               css_name = trim(css_list(nclosest))
               write (unit=*,fmt='(2a)') 'Using cohort file: ',trim(css_name)
               exit find_nonwater
            else
               !----- The closest file was no good, so we make it far away for now. -------!
               file_pdist(nclosest) = 1.e20
            end if
         end do find_nonwater
         !---------------------------------------------------------------------------------!


         !---------------------------------------------------------------------------------!
         !     Allocate the patches.                                                       !
         !---------------------------------------------------------------------------------!
         !------------------------------------------------------------------------------!
         !      We don't have information on sites, assume all sites to be the same.    !
         !------------------------------------------------------------------------------!
         do isi = 1, cpoly%nsites
            csite => cpoly%site(isi)

            !----- Allocate the patches in this site. ----------------------------------!
            call allocate_sitetype(csite,npatches)

            do ip=1,npatches

               !----- Translate the land use type. -------------------------------------!
               select case (trk(ip))
               case (-1)
                  !----- Forest plantation. --------------------------------------------!
                  csite%dist_type (ip) = 2
                  !---------------------------------------------------------------------!
               case (0)
                  !----- Agriculture. --------------------------------------------------!
                  csite%dist_type (ip) = 1
                  !---------------------------------------------------------------------!
               case (1)
                  !---------------------------------------------------------------------!
                  !     Secondary forest.  This was ambiguous in previous versions,     !
                  ! so we decide which type to point this patch based on whether        !
                  ! anthropogenic disturbance is turned on or not.  Using flags 5 or    !
                  ! 6 is unambiguous.                                                   !
                  !---------------------------------------------------------------------!
                  select case (ianth_disturb)
                  case (0)
                     !----- No anthropogenic, assume abandoned lands. ------------------!
                     csite%dist_type (ip) = 5
                     !------------------------------------------------------------------!
                  case (1,2)
                     !----- Anthropogenic, assume logging. -----------------------------!
                     csite%dist_type (ip) = 6
                     !------------------------------------------------------------------!
                  end select
                  !---------------------------------------------------------------------!
               case (2)
                  !----- Primary forest. -----------------------------------------------!
                  csite%dist_type (ip) = 3
                  !---------------------------------------------------------------------!
               case (3)
                  !----- Burnt patch. --------------------------------------------------!
                  csite%dist_type (ip) = 4
                  !---------------------------------------------------------------------!
               case (4)
                  !----- Abandoned land (secondary growth). ----------------------------!
                  csite%dist_type (ip) = 5
                  !---------------------------------------------------------------------!
               case (5)
                  !----- Logged forest. ------------------------------------------------!
                  csite%dist_type (ip) = 6
                  !---------------------------------------------------------------------!
               end select
               !------------------------------------------------------------------------!

               csite%fbeam             (ip) = 1.0
               csite%light_type        (ip) = 1
               csite%age               (ip) = age (ip)
               csite%area              (ip) = area(ip)
               csite%fast_grnd_C       (ip) =        agf_fsc  * fsc (ip)
               csite%fast_soil_C       (ip) = (1.0 - agf_fsc) * fsc (ip)
               csite%slow_soil_C       (ip) = ssc (ip)
               csite%structural_grnd_C (ip) =        agf_stsc  * stsc(ip)
               csite%structural_soil_C (ip) = (1.0 - agf_stsc) * stsc(ip)
               csite%structural_grnd_L (ip) =        agf_stsc  * stsl(ip)
               csite%structural_soil_L (ip) = (1.0 - agf_stsc) * stsl(ip)
               csite%mineralized_soil_N(ip) = msn (ip)
               csite%fast_grnd_N       (ip) =        agf_fsc  * fsn (ip)
               csite%fast_soil_N       (ip) = (1.0 - agf_fsc) * fsn (ip)
               csite%structural_grnd_N (ip) = csite%structural_grnd_C (ip)              &
                                            / c2n_structural
               csite%structural_soil_N (ip) = csite%structural_soil_C (ip)              &
                                            / c2n_structural
               csite%pname             (ip) = trim(pname(ip))
               csite%sum_dgd           (ip) = 0.0
               csite%sum_chd           (ip) = 0.0
               csite%cohort_count      (ip) = 0


               !------------------------------------------------------------------------!
               !     Check decomposition scheme before assigning microbial carbon.      !
               !------------------------------------------------------------------------!
               select case (decomp_scheme)
               case (5)
                  csite%microbial_soil_C(ip) = f0_msc * ssc(ip)
                  csite%slow_soil_C     (ip) = f0_ssc * ssc(ip)
                  csite%passive_soil_C  (ip) = f0_psc * ssc(ip)
               case default
                  csite%microbial_soil_C(ip) = 0.0
                  csite%slow_soil_C     (ip) = ssc(ip)
                  csite%passive_soil_C  (ip) = 0.0
               end select
               !------------------------------------------------------------------------!
            end do
            !---------------------------------------------------------------------------!

            !---- Initialize the cohort counts per patch. ------------------------------!
            csite%cohort_count(:) = 0
            !---------------------------------------------------------------------------!
         end do
         !------------------------------------------------------------------------------!

         close(unit=12,status='keep')

         !---------------------------------------------------------------------------------!
         !    Open the cohort file and read in all cohorts.                                !
         !---------------------------------------------------------------------------------!
         open(unit=12,file=trim(css_name),form='formatted',status='old')
         read(unit=12,fmt='(a4)')  cdum

         if (ied_init_mode == 1) then
            read(unit=12,fmt=*) !---- Skip second line. -----------------------------------!
         end if
         ic = 0

         !---------------------------------------------------------------------------------!
         !     Now we loop over all patches and decide whether they should be included     !
         ! or not.                                                                         !
         !---------------------------------------------------------------------------------!
         read_cohorts: do

            ic = ic + 1
            add_this_cohort(ic) = .true.

            !------------------------------------------------------------------------------!
            !     We must check whether we are not exceeding the maximum number of patches !
            ! that we can read.                                                            !
            !------------------------------------------------------------------------------!
            if (ic > huge_cohort) then
               write (unit=*,fmt='(a,1x,a)')  ' In file:',trim(css_name)
               write (unit=*,fmt='(a)')       ' Number of cohorts is > HUGE_COHORT...'
               write (unit=*,fmt='(a,1x,i7)') ' HUGE_COHORT:',huge_cohort
               write (unit=*,fmt='(a)')       ' Increase HUGE_COHORT to read this...'
               call fatal_error('Too many cohorts to be read...'                           &
                               ,'read_ext_ed20_history_file','ed_history_io.f90')
            end if

            !------------------------------------------------------------------------------!
            !    If we can still add new cohorts, read the next one, according to the      !
            ! input format type.                                                           !
            !------------------------------------------------------------------------------!
            read(unit=12,fmt=*,iostat=ierr) ctime(ic),cpname(ic),cname(ic),dbh(ic)      &
                                           ,hite(ic),ipft(ic),nplant(ic),bleaf(ic)      &
                                           ,broot(ic),bsapwooda(ic),bsapwoodb(ic)       &
                                           ,bbarka(ic),bbarkb(ic),bdeada(ic),bdeadb(ic) &
            !---------------------------------------------------------------------------!
            !     Check whether the file has hit the end, and if so, leave the loop.    !
            !---------------------------------------------------------------------------!
            if(ierr /= 0)exit read_cohorts


            !----- No carbon balance information.  Assign 1. ---------------------------!
            cb(1:12,ic)     = 1.0
            cb_max(1:12,ic) = 1.0


            !------------------------------------------------------------------------------!
            !     The PFT classes has changed between different ED versions, here we       !
            ! standardise this.                                                            !
            !------------------------------------------------------------------------------!
            if(renumber_pfts) then
               if (ipft(ic) < 100) then
                  ipft(ic) = ipft(ic) + 1
                  if(ipft(ic) >= 5) ipft(ic) = ipft(ic) - 3
               else
                  ipft(ic) = ipft(ic) - 100
               end if
            end if

            !----- Check if the year matches.  If not, we will ignore this cohort. --------!
            year = int(ctime(ic))
            if(use_target_year == 1 .and. year /= restart_target_year) then
               add_this_cohort(ic) = .false.
            end if

            !----- Remove cohort in case nplant > 0. --------------------------------------!
            if(nplant(ic) < tiny(1.0)) add_this_cohort(ic) = .false.


            !------------------------------------------------------------------------------!
            !     Find site and patch to which this cohort belong, and start counting how  !
            ! many we need to allocate.                                                    !
            !------------------------------------------------------------------------------!
            put_cohort:do isi=1,cpoly%nsites
               csite => cpoly%site(isi)
               do ipa=1,csite%npatches

                  !------------------------------------------------------------------------!
                  !    Here we test whether the PFT of this cohort is expected to be       !
                  ! included.                                                              !
                  !------------------------------------------------------------------------!
                  if (.not. include_pft(ipft(ic))) then
                     !----- This PFT wasn't expected... -----------------------------------!
                     select case (pft_1st_check)
                     case (0)
                        !----- Stop the run. ----------------------------------------------!
                        write (unit=*,fmt='(a,1x,i5,1x,a)')                                &
                             'I found a cohort with PFT=',ipft(ic)                         &
                            ,' and it is not in your include_these_pft...'
                        call fatal_error('Invalid PFT in history file'                     &
                                        ,'read_ext_ed20_history_file','ed_history_io.f90')

                     case (1)
                        !----- Add the unexpected PFT to the list of possible PFTs. -------!
                        write (unit=*,fmt='(a,1x,i5,1x,a)')                                &
                             'I found a cohort with PFT=',ipft(ic)                         &
                            ,'... Including this PFT in your include_these_pft...'
                        include_pft(ipft(ic))                 = .true.
                        include_these_pft(count(include_pft)) = ipft(ic)
                        call sort_up(include_these_pft,n_pft)
                        if (is_grass(ipft(ic))) include_pft_ag(ipft(ic)) = .true.

                     case (2)
                        !----- Ignore the cohort. -----------------------------------------!
                        write (unit=*,fmt='(a,1x,i5,1x,a)')                                &
                             'I found a cohort with PFT=',ipft(ic),'... Ignoring it...'
                        add_this_cohort(ic) = .false.
                     end select
                  end if

                  if (trim(csite%pname(ipa)) == trim(cpname(ic) ) .and.                    &
                      add_this_cohort(ic)                              ) then
                     csite%cohort_count(ipa) = csite%cohort_count(ipa) + 1
                     exit put_cohort
                  end if
               end do
            end do put_cohort
         end do read_cohorts

         !----- Find the total number of cohorts. -----------------------------------------!
         ncohorts = max(ic-1,0)


         close (unit=12,status='keep')

         loop_sites: do isi=1,cpoly%nsites
            csite => cpoly%site(isi)

            loop_patches: do ipa=1,csite%npatches
               cpatch => csite%patch(ipa)

               if (csite%cohort_count(ipa) /= 0) then

                  call allocate_patchtype(cpatch,csite%cohort_count(ipa))
                  csite%plant_ag_biomass(ipa) = 0.
                  ic2 = 0
                  do ic = 1,ncohorts

                     if (trim(csite%pname(ipa)) == trim(cpname(ic) ) .and.                 &
                         add_this_cohort(ic)                              ) then

                        ic2 = ic2 + 1
                        cpatch%pft(ic2) = ipft(ic)
                        cpatch%nplant(ic2) = nplant(ic)


                        !------------------------------------------------------------------!
                        !     Select which history file we are using.  We use DBH only     !
                        ! when ied_init_mode is 6 (inventory initialisation, when DBH is   !
                        ! most likely to be the actual measured variable), otherwise we    !
                        ! use BDEAD instead.                                               !
                        !------------------------------------------------------------------!
                        !----- Inventory.  Read DBH and find the other stuff. ----------!
                        cpatch%dbh(ic2)    = max(dbh(ic),min_dbh(ipft(ic)))
                        cpatch%hite(ic2)   = dbh2h(cpatch%pft(ic2),cpatch%dbh(ic2))

                        !----- Respect BDEAD fractions if they have been specified -----!
                        bdead              = size2bd(cpatch%dbh(ic2),cpatch%hite(ic2)   &
                                                       ,ipft(ic))
                        if (bdeada(ic) .gt. 0) then
                          cpatch%bdeada(ic2) = bdeada(ic)
                        else
                          cpatch%bdeada(ic2) = agf_bs(ipft(ic))  * bdead
                        end if

                        if (bdeadb(ic) .gt. 0) then
                          cpatch%bdeadb(ic2) = bdeadb(ic)
                        else
                          cpatch%bdeadb(ic2) = (1.0 - agf_bs(ipft(ic))) * bdead
                        end if
                        !------------------------------------------------------------------!

                        !------------------------------------------------------------------!
                        !     Initialise SLA with the look-up table value, this may be     !
                        ! updated during phenology initialisation, but an initial assign-  !
                        ! ment is needed to obtain area indices.                           !
                        !------------------------------------------------------------------!
                        cpatch%sla(ic2) = SLA(ipft(ic))
                        !------------------------------------------------------------------!


                        !------------------------------------------------------------------!
                        ! Respect the allometry fractions that have been specified;
                        ! otherwise use the allometry functions to define leaf and the
                        ! other live biomass pools.                                                           !
                        !------------------------------------------------------------------!
                        !----- BLEAF ------------------------------------------------------!
                        if (bleaf(ic) .gt. 0) then
                          cpatch%bleaf(ic2) = bleaf(ic)
                        else
                          cpatch%bleaf(ic2)     = size2bl(cpatch%dbh(ic2),cpatch%hite(ic2)   &
                                                       ,ipft(ic))
                        end if

                        !----- BROOT ------------------------------------------------------!
                        if (broot(ic) .gt. 0) then
                          cpatch%broot(ic2) = broot(ic)
                        else
                          cpatch%broot(ic2)     = cpatch%bleaf(ic2) * q(ipft(ic))
                        end if

                        !----- BSAPWOODA---------------------------------------------------!
                        if (bsapwooda(ic) .gt. 0) then
                          cpatch%bsapwooda(ic2) = bsapwooda(ic)
                        else
                          cpatch%bsapwooda(ic2) = agf_bs(ipft(ic))                           &
                                                * cpatch%bleaf(ic2)                          &
                                                * qsw(ipft(ic))   * cpatch%hite(ic2)
                        end if

                        !----- BSAPWOODB --------------------------------------------------!
                        if (bsapwoodb(ic) .gt. 0) then
                          cpatch%bsapwoodb(ic2) = bsapwoodb(ic)
                        else
                          cpatch%bsapwoodb(ic2) = (1.-agf_bs(ipft(ic)))                      &
                                                * cpatch%bleaf(ic2)                          &
                                                * qsw(ipft(ic))   * cpatch%hite(ic2)
                        end if

                        !----- BBARKA -----------------------------------------------------!
                        if (bbarka(ic) .gt. 0) then
                          cpatch%bbarka(ic2) = bbarka(ic)
                        else
                          cpatch%bbarka(ic2)    = agf_bs(ipft(ic))                           &
                                                * cpatch%bleaf(ic2)                          &
                                                * qbark(ipft(ic)) * cpatch%hite(ic2)
                        end if

                        !----- BBARKB -----------------------------------------------------!
                        if (bbarkb(ic) .gt. 0) then
                          cpatch%bbarkb(ic2) = bbarkb(ic)
                        else
                          cpatch%bbarkb(ic2)    = (1.-agf_bs(ipft(ic)))                      &
                                                * cpatch%bleaf(ic2)                          &
                                                * qbark(ipft(ic)) * cpatch%hite(ic2)
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Start plants with full phenology, we will take care of       !
                        ! phenology after this sub-routine.                                !
                        !------------------------------------------------------------------!
                        cpatch%phenology_status(ic2) = 0
                        !------------------------------------------------------------------!



                        !----- Assign biomass of living tissues. --------------------------!
                        cpatch%balive(ic2) = ed_balive(cpatch, ic2)
                        !------------------------------------------------------------------!


                        !----- Initialise storage biomass (after setting balive). ---------#
                        cpatch%bstorage(ic2)  = max( almost_zero                           &
                                                   , f_bstorage_init(ipft(ic)))            &
                                              * cpatch%balive(ic2)
                        !------------------------------------------------------------------!


                        !------------------------------------------------------------------!
                        !     In case we are representing trait plasticity, update traits  !
                        ! (SLA, Vm0).  This must be done before calculating LAI.           !
                        !------------------------------------------------------------------!
                        select case (trait_plasticity_scheme)
                        case (0)
                           continue
                        case default
                           call update_cohort_plastic_trait(cpatch,ic2)
                        end select
                        !------------------------------------------------------------------!



                        !----- Assign LAI, WAI, and CAI -----------------------------------!
                        call area_indices(cpatch, ic2)
                        !------------------------------------------------------------------!



                        !------------------------------------------------------------------!
                        !     Initialise the carbon balance.  We ignore the carbon balance !
                        ! even for ED-1.0, the models are so different that there is no    !
                        ! reason to use the stored value.  For initial conditions, we      !
                        ! always assume storage biomass for the previous months so         !
                        ! the scale is correct (carbon balance is given in kgC/pl).        !
                        ! similar to what is typically solved in the model.  The current   !
                        ! month carbon balance must be initialised consistently with the   !
                        ! iddmort_scheme we are using.                                     !
                        !------------------------------------------------------------------!
                        cpatch%cb         (1:12,ic2) = cpatch%bstorage(ic2)
                        cpatch%cb_lightmax(1:12,ic2) = cpatch%bstorage(ic2)
                        cpatch%cb_moistmax(1:12,ic2) = cpatch%bstorage(ic2)
                        cpatch%cb_mlmax   (1:12,ic2) = cpatch%bstorage(ic2)
                        select case (iddmort_scheme)
                        case (0)
                           !------ Storage is not accounted. ------------------------------!
                           cpatch%cb         (13,ic2) = 0.0
                           cpatch%cb_lightmax(13,ic2) = 0.0
                           cpatch%cb_moistmax(13,ic2) = 0.0
                           cpatch%cb_mlmax   (13,ic2) = 0.0
                           !---------------------------------------------------------------!
                        case (1)
                           !------ Storage is accounted. ----------------------------------!
                           cpatch%cb         (13,ic2) = cpatch%bstorage(ic2)
                           cpatch%cb_lightmax(13,ic2) = cpatch%bstorage(ic2)
                           cpatch%cb_moistmax(13,ic2) = cpatch%bstorage(ic2)
                           cpatch%cb_mlmax   (13,ic2) = cpatch%bstorage(ic2)
                           !---------------------------------------------------------------!
                        end select
                        cpatch%cbr_bar          (ic2) = 1.0
                        !------------------------------------------------------------------!



                        !----- Above ground biomass, use the allometry. -------------------!
                        cpatch%agb    (ic2) = ed_biomass(cpatch, ic2)
                        cpatch%basarea(ic2) = pio4 * cpatch%dbh(ic2) * cpatch%dbh(ic2)
                        cpatch%btimber(ic2) = size2bt( cpatch%dbh       (ic2)              &
                                                     , cpatch%hite      (ic2)              &
                                                     , cpatch%bdeada    (ic2)              &
                                                     , cpatch%bsapwooda (ic2)              &
                                                     , cpatch%bbarka    (ic2)              &
                                                     , cpatch%pft       (ic2) )
                        cpatch%thbark (ic2) = size2xb( cpatch%dbh       (ic2)              &
                                                     , cpatch%hite      (ic2)              &
                                                     , cpatch%bbarka    (ic2)              &
                                                     , cpatch%bbarkb    (ic2)              &
                                                     , cpatch%pft       (ic2) )

                        !----- Growth rates, start with zero. -----------------------------!
                        cpatch%dagb_dt  (ic2)  = 0.
                        cpatch%dlnagb_dt(ic2)  = 0.
                        cpatch%dba_dt   (ic2)  = 0.
                        cpatch%dlnba_dt (ic2)  = 0.
                        cpatch%ddbh_dt  (ic2)  = 0.
                        cpatch%dlndbh_dt(ic2)  = 0.

                        !------------------------------------------------------------------!
                        !      Initialise other cohort variables.  Some of them won't be   !
                        ! updated unless the lai exceeds lai_min.                          !
                        !------------------------------------------------------------------!
                        cpatch%fsw(ic2)   = 1.0
                        cpatch%gpp(ic2)   = 0.0
                        cpatch%par_l(ic2) = 0.0
                        !------------------------------------------------------------------!


                        !----- Update the patch level above-ground biomass. ---------------!
                        csite%plant_ag_biomass(ipa) = csite%plant_ag_biomass(ipa)          &
                                                    + cpatch%agb(ic2) * cpatch%nplant(ic2)
                         !------------------------------------------------------------------!
                    end if
                  end do
               end if
            end do loop_patches
         end do loop_sites

         close (unit=12,status='keep')

         !----- Initialise all the other site-, patch-, and cohort-level variables. -------!
         do isi = 1,cpoly%nsites

            area_sum = 0.0
            ncohorts = 0


            !----- Make sure that the total patch area is 1. ------------------------------!
            csite => cpoly%site(isi)
            area_tot      = sum(csite%area(1:csite%npatches))
            csite%area(:) = csite%area(:)/area_tot

            !----- Find the patch-level LAI, WAI, and CAI. --------------------------------!
            do ipa=1,csite%npatches
               area_sum        = area_sum + csite%area(ipa)

               cpatch => csite%patch(ipa)
               do ico = 1,cpatch%ncohorts
                  ncohorts        = ncohorts + 1
               end do
            end do

            !----- Initialise the cohort variables, then sort them by size. ---------------!
            do ipa = 1,csite%npatches
               cpatch => csite%patch(ipa)
               do ico = 1,cpatch%ncohorts
                  call init_ed_cohort_vars(cpatch,ico,cpoly%lsl(isi),nzg                   &
                                          ,cpoly%ntext_soil(:,isi))
               end do

               !----- Make sure that cohorts are organised from tallest to shortest. ------!
               call sort_cohorts(cpatch)
            end do

            !----- Initialise the patch-level variables. ----------------------------------!
            call init_ed_patch_vars(csite,1,csite%npatches,cpoly%lsl(isi))

            !----- Make sure that patches are organised from oldest to youngest. ----------!
            call sort_patches(csite)
         end do

         !----- Initialise site-level variables. ------------------------------------------!
         call init_ed_site_vars(cpoly)


         !----- Get a diagnostic on the polygon's vegetation. -----------------------------!
         ncohorts = 0

         do isi = 1,cpoly%nsites
            nsitepat = 0
            csite => cpoly%site(isi)

            do ipa = 1,csite%npatches
               npatchco        = 0

               cpatch => csite%patch(ipa)
               do ico = 1,cpatch%ncohorts
                  ncohorts        = ncohorts+1
                  npatchco        = npatchco+1
               end do
               csite%cohort_count(ipa) = npatchco
               nsitepat                = nsitepat + 1
            end do

            cpoly%patch_count(isi) = nsitepat
         end do
      end do polyloop

      !----- Initialise the polygon-level variables. --------------------------------------!
      call init_ed_poly_vars(cgrid)
   end do gridloop

   return
end subroutine read_ext_ed20_history_file
!==========================================================================================!
!==========================================================================================!
