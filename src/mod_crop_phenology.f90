module mod_crop_phenology
    use mod_grid
    use mod_utility, only: round, round_2darray
    implicit none

    type file_phenology_r!
        integer::unit                           ! unit associated to the file
        real(dp),dimension(:,:),pointer::tab    ! daily crop parameters for each land uses
    end type file_phenology_r!
    !
    type file_phenology_i!
        integer::unit                           ! unit associated to the file
        integer,dimension(:,:),pointer::tab     ! daily crop parameters for each land uses
    end type file_phenology_i!
    
    type k_cb_matrices
        real(dp), dimension(:,:), pointer::low  ! minimum kcb
        real(dp), dimension(:,:), pointer::high ! maximum kcb
        real(dp), dimension(:,:), pointer::mid  ! kcb between the 1st and 2nd stage
    end type k_cb_matrices

    type crop_pheno_info!
        ! include all crop parameters for each calculation cell
        integer,dimension(:,:),pointer::ii0                 ! start growing day of the crop [doy]
        integer,dimension(:,:),pointer::iid                 ! duration of the growing period [day]
        integer,dimension(:,:),pointer::iie                 ! harvest day [doy]
        integer,dimension(:),pointer::n_crops_by_year       ! number of crop in a year
        type(file_phenology_r)::k_cb                        ! base crop coefficient
        type(file_phenology_r)::h                           ! crop height
        type(file_phenology_r)::z_r                         ! root depth
        type(file_phenology_r)::lai                         ! leaf area index
        type(file_phenology_r)::f_c                         ! cover fraction
        type(file_phenology_r)::r_stress                    ! plant resistance to water stress
        type(file_phenology_i)::cn_day                      ! %AB% soil moisture adjusted cn value 
        type(k_cb_matrices)::kcb_phases                     ! k_cb change points during phenology
        real(dp),dimension(:,:),pointer::p_raw_const        ! readily available water factor [-]
        real(dp),dimension(:,:),pointer::a                  ! interception coefficient according to the Von Hoyningen-Hune & Braden method [-]
        real(dp),dimension(:,:),pointer::d_r_max            ! maximum root depth [m]
        real(dp),dimension(:,:),pointer::max_RF_t           ! maximum fraction of active roots in transpirative layer [-]
        real(dp),dimension(:,:),pointer::T_lim              ! limit temperature threshold - thermic stress [°C]
        real(dp),dimension(:,:),pointer::T_crit             ! critical temperature threshold - thermic stress [°C]
        real(dp),dimension(:,:),pointer::HI                 ! harvest index
        real(dp),dimension(:,:),pointer::Ky_tot             ! water stress coefficient - overall
        real(dp),dimension(:,:,:),pointer::Ky_pheno         ! water stress coefficient - phases
        integer,dimension(:,:),pointer::irrigation_class    ! 1 = irrigated, 0 = not irrigated
        integer,dimension(:,:),pointer::cn_class            ! CN class
        real(dp),dimension(:,:,:),pointer::wp_adj           ! normalized biomass water productivity
         
    end type crop_pheno_info   !es: valday(:)%kcb%unit; valday(:)%kcb%tab(:,:)!

    type crop_pars_matrices!
        ! store the crop parameters for each calculation cells
        ! see crop_pheno_info for details
        real(dp),dimension(:,:),pointer::k_cb
        real(dp),dimension(:,:),pointer::h
        real(dp),dimension(:,:),pointer::d_r
        real(dp),dimension(:,:),pointer::lai
        integer,dimension(:,:),pointer::cn_day
        real(dp),dimension(:,:),pointer::f_c
        integer,dimension(:,:),pointer::irrigation_class
        integer,dimension(:,:),pointer::cn_class
        real(dp),dimension(:,:),pointer::p
        real(dp),dimension(:,:),pointer::a
        real(dp),dimension(:,:),pointer::d_t_max
        real(dp),dimension(:,:),pointer::RF_t_max
        real(dp),dimension(:,:),pointer::RF_e
        real(dp),dimension(:,:),pointer::RF_t
        real(dp),dimension(:,:),pointer::T_lim
        real(dp),dimension(:,:),pointer::T_crit
        real(dp),dimension(:,:),pointer::HI
        real(dp),dimension(:,:),pointer::Ky_tot
        real(dp),dimension(:,:,:),pointer::Ky_pheno
        real(dp),dimension(:,:),pointer::k_cb_low
        real(dp),dimension(:,:),pointer::k_cb_mid
        real(dp),dimension(:,:),pointer::k_cb_high
        real(dp),dimension(:,:),pointer::wp_adj
        real(dp),dimension(:,:),pointer::p_day
        real(dp),dimension(:,:),pointer::k_cb_old           ! k_cb of previous day
        integer,dimension(:,:),pointer::n_crop_in_year
        integer,dimension(:,:),pointer::pheno_idx           ! phenological stage index
        real(dp),dimension(:,:),pointer::r_stress           ! plant resistance to (water) stress
    end type crop_pars_matrices!

    type crop_matrices
    ! crop_pheno_info for details    
        integer,dimension(:,:,:),pointer::ii0
        integer,dimension(:,:,:),pointer::iie
        real(dp),dimension(:,:,:),pointer::iid
        integer,dimension(:,:,:),pointer::ii0_ref           ! harvest date in meteorological reference series
        integer,dimension(:,:,:),pointer::iie_ref           ! harvest date in meteorological reference series
        real(dp),dimension(:,:,:),pointer::iid_ref          ! crop cycle length in meteorological reference series
        real(dp),dimension(:,:,:),pointer::dij              ! coefficient of crop cycle expansion in relation to meteorological reference series
        integer,dimension(:,:,:),pointer::TSP_high
        integer,dimension(:,:,:),pointer::TSP_low
        real(dp),dimension(:,:,:),pointer::wp_adj
        real(dp),dimension(:,:,:),pointer::HI
        real(dp),dimension(:,:,:),pointer::Ky_tot
        real(dp),dimension(:,:,:,:),pointer::Ky_pheno
        real(dp),dimension(:,:,:),pointer::T_crit
        real(dp),dimension(:,:,:),pointer::T_lim
        real(dp),dimension(:,:,:),pointer::k_cb_min
        real(dp),dimension(:,:,:),pointer::k_cb_mid
        real(dp),dimension(:,:,:),pointer::k_cb_max
    end type crop_matrices

    contains

    subroutine make_random_emergence(info_pheno,meteo_weight,dir_meteo,domain,soiluse,crop_mat,irandom,year_length)!
        !randomization of emergence date from phenological series at meteorological stations

        implicit none!
        real(dp),dimension(:,:,:),intent(in)::meteo_weight!
        integer,dimension(:,:,:),intent(in)::dir_meteo!
        type(grid_i),intent(in)::domain!
        integer,dimension(:,:),intent(in)::soiluse!
        type(crop_pheno_info),dimension(:),intent(in)::info_pheno!
        type(crop_matrices),intent(out)::crop_mat
        integer,dimension(:,:),intent(in)::irandom
        integer,intent(in)::year_length
        
        real(dp),dimension(size(crop_mat%ii0,1),size(crop_mat%ii0,2),size(crop_mat%ii0,3))::ii0_r
        real(dp),dimension(size(crop_mat%iie,1),size(crop_mat%iie,2),size(crop_mat%iie,3))::iie_r
        integer::i,j,k,z
        
        ! initialization
        ii0_r = 0.d0
        iie_r = 0.d0
        
        ! spatial distribution of emergence date and phenological cycle duration by using meteorological weights
        do k=1,size(meteo_weight,3)!
            do j=1,size(domain%mat,2)!
                do i=1,size(domain%mat,1)!
                    if(domain%mat(i,j)/=domain%header%nan)then!
                        ii0_r(i,j,:)   =info_pheno(dir_meteo(i,j,k))%ii0(soiluse(i,j),:)*meteo_weight(i,j,k) + ii0_r(i,j,:)
                        iie_r(i,j,:)   =info_pheno(dir_meteo(i,j,k))%iie(soiluse(i,j),:)*meteo_weight(i,j,k) + iie_r(i,j,:)
                        crop_mat%iid(i,j,:)=info_pheno(dir_meteo(i,j,k))%iid(soiluse(i,j),:)*meteo_weight(i,j,k) &
                                             + crop_mat%iid(i,j,:)  ! crop cycle length (real)
                    end if!
                end do!                
            end do!
        end do!

        do z=1,size(crop_mat%ii0,3)
            crop_mat%ii0(:,:,z)=merge(nint(ii0_r(:,:,z)),domain%header%nan,domain%mat/=domain%header%nan)   ! emergence date
            crop_mat%iie(:,:,z)=merge(nint(iie_r(:,:,z)),domain%header%nan,domain%mat/=domain%header%nan)   ! harvest date
        end do
        
        do j=1,size(domain%mat,2)!
            do i=1,size(domain%mat,1)!
                if(domain%mat(i,j)/=domain%header%nan)then!
                    do z=1,size(crop_mat%ii0,3)
                        crop_mat%ii0_ref(i,j,:) = info_pheno(dir_meteo(i,j,1))%ii0(soiluse(i,j),:)   ! emergence date in meteorological reference series
                        crop_mat%iie_ref(i,j,:) = info_pheno(dir_meteo(i,j,1))%iie(soiluse(i,j),:)   ! harvest date in meteorological reference series
                        crop_mat%iid_ref(i,j,:) = info_pheno(dir_meteo(i,j,1))%iid(soiluse(i,j),:)   ! crop cycle length in meteorological reference series
                        crop_mat%dij(i,j,:)     = crop_mat%iid_ref(i,j,:)/crop_mat%iid(i,j,:)        ! coefficient of crop cycle expansion in relation to meteorological reference series
                    end do
                end if
            end do
        end do
        
        do z=1, size(crop_mat%ii0,3)
            crop_mat%ii0_ref(:,:,z) = merge(crop_mat%ii0_ref(:,:,z),domain%header%nan,domain%mat/=domain%header%nan)
            crop_mat%ii0_ref(:,:,z) = merge(crop_mat%ii0_ref(:,:,z),domain%header%nan,crop_mat%ii0_ref(:,:,z)/=0)
            crop_mat%iie_ref(:,:,z) = merge(crop_mat%iie_ref(:,:,z),domain%header%nan,domain%mat/=domain%header%nan)
            crop_mat%iie_ref(:,:,z) = merge(crop_mat%iie_ref(:,:,z),domain%header%nan,crop_mat%iie_ref(:,:,z)/=0)
            crop_mat%iid_ref(:,:,z) = merge(crop_mat%iid_ref(:,:,z),dble(domain%header%nan),domain%mat/=domain%header%nan)
            crop_mat%iid_ref(:,:,z) = merge(crop_mat%iid_ref(:,:,z),dble(domain%header%nan),crop_mat%iid_ref(:,:,z)/=0)
            crop_mat%dij(:,:,z)     = merge(crop_mat%dij(:,:,z),dble(domain%header%nan),domain%mat/=domain%header%nan)
            crop_mat%dij(:,:,z)     = merge(crop_mat%dij(:,:,z),dble(domain%header%nan),crop_mat%iid(:,:,z)/=0)
        end do
        
        do z=1,size(crop_mat%ii0,3)
            crop_mat%TSP_high(:,:,z) = nint(0.75 * crop_mat%iid(:,:,z)) + crop_mat%ii0(:,:,z) + irandom
            crop_mat%TSP_low(:,:,z)  = nint(0.45 * crop_mat%iid(:,:,z)) + crop_mat%ii0(:,:,z) + irandom
        end do
          
        ! values adjustment for double-years crops
        crop_mat%TSP_high = merge(crop_mat%TSP_high-year_length, crop_mat%TSP_high, crop_mat%TSP_high>year_length)
        crop_mat%TSP_low  = merge(crop_mat%TSP_low-year_length,  crop_mat%TSP_low,  crop_mat%TSP_low>year_length)
        
    end subroutine make_random_emergence!

    subroutine populate_crop_yield_matrices(info_pheno,dir_phenofases,domain,soil_use,crop_mat,year)
        ! populate of crop matrices with crop yield parameters
        implicit none!
        type(crop_pheno_info),dimension(:),intent(in)::info_pheno
        integer,dimension(:,:),intent(in)::dir_phenofases!
        type(grid_i),intent(in)::domain
        integer,dimension(:,:),intent(in)::soil_use
        type(crop_matrices),intent(inout)::crop_mat
        integer,intent(in)::year
        integer::i,j
        
        do j=1,size(domain%mat,2)!
            do i=1,size(domain%mat,1)!
                if(domain%mat(i,j)/=domain%header%nan)then!
                    crop_mat%wp_adj(i,j,:) = info_pheno(dir_phenofases(i,j))%wp_adj(soil_use(i,j),:,year)
                    crop_mat%HI(i,j,:) = info_pheno(dir_phenofases(i,j))%HI(soil_use(i,j),:)
                    crop_mat%Ky_tot(i,j,:) = info_pheno(dir_phenofases(i,j))%Ky_tot(soil_use(i,j),:)
                    crop_mat%Ky_pheno(i,j,:,:) = info_pheno(dir_phenofases(i,j))%Ky_pheno(soil_use(i,j),:,:)
                    crop_mat%T_crit(i,j,:) = info_pheno(dir_phenofases(i,j))%T_crit(soil_use(i,j),:)
                    crop_mat%T_lim(i,j,:) = info_pheno(dir_phenofases(i,j))%T_lim(soil_use(i,j),:)
                    crop_mat%k_cb_min(i,j,:) = info_pheno(dir_phenofases(i,j))%kcb_phases%low(soil_use(i,j),:)
                    crop_mat%k_cb_mid(i,j,:) = info_pheno(dir_phenofases(i,j))%kcb_phases%mid(soil_use(i,j),:)
                    crop_mat%k_cb_max(i,j,:) = info_pheno(dir_phenofases(i,j))%kcb_phases%high(soil_use(i,j),:)
                end if
            end do
        end do

    end subroutine populate_crop_yield_matrices

    subroutine compute_regrow_window(info_pheno, regrow_start, regrow_end, season_end)
        ! FORCED CUTS support.
        ! For every (weather station, crop) pair, locate in the REFERENCE phenological
        ! series one complete regrowth cycle, i.e. the interval that goes from a cut
        ! (Kcb back to its minimum) up to the day before the following cut.
        ! That window is later used to replay the original regrowth shape starting
        ! from an externally imposed cut date, so the curve computed by cropcoef is
        ! preserved and only re-anchored in time.
        implicit none
        type(crop_pheno_info),dimension(:),intent(in)::info_pheno
        integer,dimension(:,:),allocatable,intent(out)::regrow_start, regrow_end
        ! season_end: last day of active growth in the reference series. After that day the
        ! cell must go back to the ORIGINAL curve, so that autumn senescence and winter
        ! dormancy are the real ones and the crop does not stay green until December.
        integer,dimension(:,:),allocatable,intent(out)::season_end
        integer::n_ws, n_crop, n_days, i, c, d, d_start, d_end, cut1, cut2
        real(dp)::peak, base, thr, tol, mature

        n_ws   = size(info_pheno)
        ! La tabella pheno prodotta da cropcoef contiene tutti gli anni della
        ! simulazione uno dopo l'altro (es. 1096 righe per 2012-2014), ma il
        ! puntatore fenologico viene sempre limitato a year_length: le righe oltre
        ! il primo anno non vengono mai lette. Cercare qui la stagione su tutta la
        ! tabella darebbe un season_end di ~1000 e renderebbe inutile il controllo
        ! di fine stagione. Si guarda quindi solo il primo anno.
        n_days = min(size(info_pheno(1)%k_cb%tab,1), 366)
        n_crop = size(info_pheno(1)%k_cb%tab,2)
        allocate(regrow_start(n_ws,n_crop))
        allocate(regrow_end(n_ws,n_crop))
        allocate(season_end(n_ws,n_crop))
        regrow_start = 0
        regrow_end   = 0
        season_end   = 0
        tol = 1.0d-6

        do i=1,n_ws
            do c=1,n_crop
                ! ---- growing season: Kcb is zero during winter dormancy ----
                d_start = 0
                do d=1,n_days
                    if (info_pheno(i)%k_cb%tab(d,c) > tol) then
                        d_start = d
                        exit
                    end if
                end do
                if (d_start == 0) cycle          ! no crop at all on this land use

                d_end = 0
                do d=n_days,1,-1
                    if (info_pheno(i)%k_cb%tab(d,c) > tol) then
                        d_end = d
                        exit
                    end if
                end do
                season_end(i,c) = d_end

                ! ---- cuts are SHARP DROPS inside the season, NOT returns to zero ----
                ! Between two cuts Kcb falls back to its "cut" level (e.g. 0.30 for alfalfa),
                ! while zero only means winter dormancy. Looking for the absolute minimum
                ! would therefore find the end of the season instead of the first cut.
                peak = maxval(info_pheno(i)%k_cb%tab(d_start:d_end,c))
                base = minval(info_pheno(i)%k_cb%tab(d_start:d_end,c), &
                    & mask=info_pheno(i)%k_cb%tab(d_start:d_end,c) > tol)
                thr = 0.2d0 * (peak - base)
                if (thr <= tol) cycle             ! flat curve: no cut to replay

                ! ---- a cut always starts from a MATURE canopy ----
                ! Requiring the drop to begin above half way between the post-cut level
                ! and the seasonal peak rules out the drops that are not cuts: the short
                ! zero block cropcoef writes at the beginning of the series, and the
                ! senescence steps. Without this test, when HarvestDate_max is 366 (so
                ! that Kcb never returns to zero in winter) the scan starts on 1 January
                ! and takes that initial artefact as the first cut: the replayed
                ! "regrowth" is then a flat line at the post-cut level and the crop
                ! never grows again after a forced cut.
                mature = base + 0.5d0 * (peak - base)

                cut1 = 0
                cut2 = 0
                do d=d_start+1,d_end
                    if (info_pheno(i)%k_cb%tab(d-1,c) - info_pheno(i)%k_cb%tab(d,c) > thr &
                        & .and. info_pheno(i)%k_cb%tab(d-1,c) > mature) then
                        if (cut1 == 0) then
                            cut1 = d
                        else
                            cut2 = d
                            exit
                        end if
                    end if
                end do
                if (cut1 == 0) cycle              ! crop without cuts: leave it to the standard engine

                regrow_start(i,c) = cut1          ! day of the cut: Kcb is at its post-cut level
                if (cut2 > 0) then
                    regrow_end(i,c) = cut2 - 1    ! regrowth ends the day before the next cut
                else
                    regrow_end(i,c) = d_end
                end if
            end do
        end do

    end subroutine compute_regrow_window

    subroutine populate_crop_pars_matrices(crop_pars_mat,info_pheno,irandom,doy,ws_idx,&
                                           & domain,soil_use,y, year_length, crop_mat, &
                                           & fc_last_cut, fc_regrow_start, fc_regrow_end, fc_season_end, &
                                           & fc_first_cut)!
        ! populate crop parameters matrices from weather stations time series
        integer,intent(in)::doy,year_length,y!
        type(grid_i),intent(in)::domain,soil_use!
        integer,dimension(:,:),intent(in)::ws_idx                   !index in the list of weather stations
        type(crop_pheno_info),dimension(:),intent(in)::info_pheno
        type(crop_pars_matrices),intent(inout)::crop_pars_mat
        integer,dimension(:,:),intent(inout)::irandom               ! pseudorandom parameter that shifts crop cycle
        type(crop_matrices),intent(in)::crop_mat
        ! ---- FORCED CUTS (optional): externally imposed harvest dates per cell ----
        ! fc_last_cut(i,j)   : doy of the most recent forced cut for that cell (0 = none yet this year)
        ! fc_regrow_start/end: for each (weather station, crop), the window in the reference
        !                      phenological series describing one regrowth cycle.
        ! When a cell has had a forced cut, the phenological pointer (doy_s) is re-anchored to
        ! the start of the reference regrowth, so the ORIGINAL curve shape (from cropcoef) is
        ! preserved but restarted on the imposed date. All arguments are optional: when absent
        ! the routine behaves exactly like the standard version.
        integer,dimension(:,:),intent(in),optional::fc_last_cut
        integer,dimension(:,:),intent(in),optional::fc_regrow_start
        integer,dimension(:,:),intent(in),optional::fc_regrow_end
        integer,dimension(:,:),intent(in),optional::fc_season_end
        ! fc_first_cut(i,j): doy of the FIRST imposed cut of the year for that cell (0 = none).
        ! Needed to prevent the GDD calendar from firing its own cut before the satellite one.
        integer,dimension(:,:),intent(in),optional::fc_first_cut

        integer::i,j
        integer::doy_s ! shifted day of the year
        integer::doy_std ! standard (not re-anchored) pointer, kept for the end of season
        logical::use_forced_cuts
        integer::rs, re, se
        integer::d_end_season, overshoot

        use_forced_cuts = present(fc_last_cut) .and. present(fc_regrow_start) &
            & .and. present(fc_regrow_end) .and. present(fc_season_end) .and. present(fc_first_cut)

        do j=1,size(domain%mat,2)!
            do i=1,size(domain%mat,1)!
                scans_domain: if(domain%mat(i,j) /= domain%header%nan)then!
                    if (crop_mat%ii0(i,j,crop_pars_mat%n_crop_in_year(i,j)) == 0 &
                        & .and. crop_mat%iie(i,j,crop_pars_mat%n_crop_in_year(i,j)) == 0 ) then  ! no crop
                        doy_s = doy - irandom(i,j)
                    else if (crop_mat%ii0(i,j,crop_pars_mat%n_crop_in_year(i,j)) < &
                        crop_mat%iie(i,j,crop_pars_mat%n_crop_in_year(i,j))) then            ! annuals or perennials
                        ! emergence date is shifted as ii0(i,j,cs)-irandom(i,j,cs)
                        ! nint((gg-ii0(i,j))*dij(i,j)) contracts/expands the series
                        ! randomization of emergence date (ii0/irandom) and factor of dilatation (dij) are used to calculate gg1
                        doy_s = crop_mat%ii0_ref(i,j,crop_pars_mat%n_crop_in_year(i,j)) - irandom(i,j) + &
                            & nint((doy-crop_mat%ii0(i,j,crop_pars_mat%n_crop_in_year(i,j))) * &
                            & crop_mat%dij(i,j,crop_pars_mat%n_crop_in_year(i,j)))
                    else                                                                                                ! biennals
                        if (doy <= crop_mat%iie(i,j,crop_pars_mat%n_crop_in_year(i,j)) + irandom(i,j)) then !%EAC% add equal condition
                            ! from 1/1 to harvest date, only the contraction/expansion of crop cycle is taken into account
                            ! in the first part of the year, the limits are 1/1 and harvest date (only contraction/expansion)
                            doy_s = doy * (crop_mat%iie_ref(i,j,crop_pars_mat%n_crop_in_year(i,j))) / &
								& crop_mat%iie(i,j,crop_pars_mat%n_crop_in_year(i,j)) - irandom(i,j)
                        else if (doy >= crop_mat%ii0(i,j,crop_pars_mat%n_crop_in_year(i,j)) + irandom(i,j)) then
                            
                            ! from emergence to 31/12, both parameters are taken into account: the randomization of emergence date and the contraction/expansion of crop cycle
                            ! in the second part of the year, the limits are the randomized emergence date and 31/12 (day 365/366)
                            doy_s = (doy - crop_mat%ii0(i,j,crop_pars_mat%n_crop_in_year(i,j))) * &
								& (year_length - crop_mat%ii0_ref(i,j,crop_pars_mat%n_crop_in_year(i,j))) / &
                                & (year_length - crop_mat%ii0(i,j,crop_pars_mat%n_crop_in_year(i,j))) + &
                                & crop_mat%ii0_ref(i,j,crop_pars_mat%n_crop_in_year(i,j)) - irandom (i,j)
                        else
                            doy_s = crop_mat%iie_ref(i,j,crop_pars_mat%n_crop_in_year(i,j)) +1    ! it points to a null Kcb
                        end if
                    end if
                    
                    ! adjust pointing if gg1 doesn't belong to [1,year_length] (it takes into account irandom and not agricultural soil uses)
                    if (doy_s < 1) doy_s = 1
                    if (doy_s > year_length) doy_s = year_length

                    ! %EAC% fix back shifting
                    ! if the second crop, then check if new doy overlap the series of the previous crop
                    ! if so, delete the overlap
                    if (crop_pars_mat%n_crop_in_year(i,j) == 2) then !
                        if (doy_s < crop_mat%iie(i,j,1)) then ! compare with the first crop
                                doy_s = doy
                        end if
                    end if
                    
                    ! ---- FORCED CUTS: re-anchor the phenological pointer ----------------
                    ! If this cell has already received an externally imposed cut this year,
                    ! the pointer is moved to the beginning of the reference regrowth cycle
                    ! and advances one day per day. The shape of the curve is therefore the
                    ! original one computed by cropcoef, only restarted on the imposed date.
                    ! Before the first forced cut of the year the standard pointer is kept,
                    ! so spring green-up is left untouched.
                    if (use_forced_cuts) then
                        rs = fc_regrow_start(ws_idx(i,j), soil_use%mat(i,j))
                        re = fc_regrow_end(ws_idx(i,j), soil_use%mat(i,j))
                        se = fc_season_end(ws_idx(i,j), soil_use%mat(i,j))
                        doy_std = doy_s                     ! standard pointer, kept as reference
                        if (rs > 0 .and. re >= rs .and. doy_std > se .and. fc_last_cut(i,j) > 0) then
                            ! (c) STAGIONE DI RIFERIMENTO FINITA, ma la cella ha ricevuto
                            ! tagli da satellite. La curva di cropcoef qui vale zero, perche'
                            ! HarvestDate_max chiude la stagione: seguirla farebbe crollare
                            ! l'ET a picco, cosa che le misure non mostrano. La coltura viene
                            ! invece tenuta nello stato POST-TAGLIO (rs), cioe' canopy bassa
                            ! ma viva: Kcb, LAI, altezza e profondita' radicale restano quelli
                            ! del giorno del taglio nella serie di riferimento, e l'ET scende
                            ! da sola perche' scende ET0. Le celle senza tagli imposti non
                            ! passano di qui e mantengono il comportamento originale.
                            ! Il rientro non e' istantaneo: circa il 10% delle celle ha
                            ! l'ultimo taglio prima del giorno 230 e a fine stagione si trova
                            ! ancora sul plateau, quindi un salto secco allo stato post-taglio
                            ! sarebbe l'ennesimo gradino verticale. Il puntatore ripercorre
                            ! invece all'indietro la finestra di ricrescita, un giorno al giorno:
                            ! la canopy cala seguendo esattamente la forma con cui era cresciuta,
                            ! i salti giornalieri sono quelli della rampa di crescita (che il
                            ! modello gia' accetta) e il raccordo con l'ultimo giorno di stagione
                            ! e' continuo. Arrivata allo stato post-taglio la coltura ci resta.
                            ! overshoot: di quanto il puntatore standard ha superato la fine
                            ! della stagione di riferimento. Si misura sul puntatore, non sul
                            ! calendario, perche' i due non avanzano alla stessa velocita'
                            ! (dij contrae o dilata il ciclo). Cosi' il raccordo con l'ultimo
                            ! giorno di stagione e' esatto qualunque sia dij.
                            overshoot = doy_std - se
                            d_end_season = rs + ((doy - overshoot) - fc_last_cut(i,j))
                            if (d_end_season > re) d_end_season = re
                            if (d_end_season < rs) d_end_season = rs
                            doy_s = d_end_season - overshoot
                            if (doy_s < rs) doy_s = rs
                            if (doy_s > re) doy_s = re
                        else if (rs > 0 .and. re >= rs .and. doy_std <= se) then
                            if (fc_last_cut(i,j) > 0) then
                                ! (a) after an imposed cut: replay the reference regrowth from it
                                doy_s = rs + (doy - fc_last_cut(i,j))
                                ! hold on the mature plateau if the next cut is late: the crop
                                ! must not run into the next cut of the reference series
                                if (doy_s > re) doy_s = re
                                if (doy_s < 1) doy_s = 1
                                if (doy_s > year_length) doy_s = year_length
                            else if (fc_first_cut(i,j) > 0) then
                                ! (b) BEFORE the first imposed cut of the year: spring green-up is
                                ! left untouched, but the pointer is not allowed to reach the cut
                                ! of the reference curve. Without this the GDD calendar would fire
                                ! a spurious cut of its own before the satellite one.
                                if (rs > 1 .and. doy_s > rs - 1) doy_s = rs - 1
                            end if
                        end if
                    end if
                    ! ---------------------------------------------------------------------

                    ! TODO: parameters overloading can be moved to a specific subroutine

                    ! conveniently updates phenological data from its series - update occurs only if Kcb varies
                    crop_pars_mat%k_cb(i,j)=info_pheno(ws_idx(i,j))%k_cb%tab(doy_s,soil_use%mat(i,j))
                    
                    if (crop_pars_mat%k_cb(i,j) /= crop_pars_mat%k_cb_low(i,j) .or. &
                        crop_pars_mat%k_cb_old(i,j) /= crop_pars_mat%k_cb_low(i,j)) then
                        if (crop_pars_mat%k_cb_old(i,j) > crop_pars_mat%k_cb_low(i,j) .and. &
                            crop_pars_mat%k_cb(i,j) == crop_pars_mat%k_cb_low(i,j) &
                            & .and. info_pheno(ws_idx(i,j))%n_crops_by_year(soil_use%mat(i,j))>1) then
                            if (crop_pars_mat%n_crop_in_year(i,j) < &
                                info_pheno(ws_idx(i,j))%n_crops_by_year(soil_use%mat(i,j))) then       ! cult_switch cycle
                                crop_pars_mat%n_crop_in_year(i, j) = crop_pars_mat%n_crop_in_year(i, j) + 1
                                crop_pars_mat%pheno_idx(i,j) = 1
                            else if (crop_pars_mat%n_crop_in_year(i,j) == &
                                info_pheno(ws_idx(i,j))%n_crops_by_year(soil_use%mat(i,j))) then ! if cult_switch cycle ends, switch is set back to 1
                                crop_pars_mat%n_crop_in_year(i, j) = 1
                                crop_pars_mat%pheno_idx(i,j) = 1
                            end if
                        end if
                        crop_pars_mat%h(i,j)=info_pheno(ws_idx(i,j))%h%tab(doy_s,soil_use%mat(i,j))
                        crop_pars_mat%d_r(i,j)=info_pheno(ws_idx(i,j))%z_r%tab(doy_s,soil_use%mat(i,j))
                        crop_pars_mat%lai(i,j)=info_pheno(ws_idx(i,j))%lai%tab(doy_s,soil_use%mat(i,j))
                        crop_pars_mat%cn_day(i,j)=info_pheno(ws_idx(i,j))%cn_day%tab(doy_s,soil_use%mat(i,j))
                        crop_pars_mat%f_c(i,j)=info_pheno(ws_idx(i,j))%f_c%tab(doy_s,soil_use%mat(i,j))
                        crop_pars_mat%r_stress(i,j)=info_pheno(ws_idx(i,j))%r_stress%tab(doy_s,soil_use%mat(i,j))
                                               
                        crop_pars_mat%irrigation_class(i,j)=&
                            info_pheno(ws_idx(i,j))%irrigation_class(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j))
                        crop_pars_mat%cn_class(i,j)= &
                            info_pheno(ws_idx(i,j))%cn_class(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j))
                        ! TEST: replace here with variable p value
                        !crop_pars_mat%p(i,j)=info_pheno(ws_idx(i,j))%p_raw%tab(doy_s,soil_use%mat(i,j))
                        crop_pars_mat%p(i,j)= &
                            info_pheno(ws_idx(i,j))%p_raw_const(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j))

                        crop_pars_mat%a(i,j)=&
                            info_pheno(ws_idx(i,j))%a(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j))
                        crop_pars_mat%d_t_max(i,j)= &
                            info_pheno(ws_idx(i,j))%d_r_max(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j))
                        crop_pars_mat%RF_t_max(i,j)= &
                            info_pheno(ws_idx(i,j))%max_RF_t(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j))
                        crop_pars_mat%T_lim(i,j)= &
                            info_pheno(ws_idx(i,j))%T_lim(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j))
                        crop_pars_mat%T_crit(i,j)= &
                            info_pheno(ws_idx(i,j))%T_crit(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j))
                        crop_pars_mat%HI(i,j)= &
                            info_pheno(ws_idx(i,j))%HI(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j))
                        crop_pars_mat%Ky_tot(i,j)= &
                            info_pheno(ws_idx(i,j))%Ky_tot(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j))
                        crop_pars_mat%Ky_pheno(i,j,:)= &
                            info_pheno(ws_idx(i,j))%Ky_pheno(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j),:)
                        crop_pars_mat%k_cb_low(i,j)= &
                            info_pheno(ws_idx(i,j))%kcb_phases%low(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j))
                        crop_pars_mat%k_cb_mid(i,j)= &
                            info_pheno(ws_idx(i,j))%kcb_phases%mid(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j))
                        crop_pars_mat%k_cb_high(i,j)= &
                            info_pheno(ws_idx(i,j))%kcb_phases%high(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j))
                        crop_pars_mat%wp_adj(i,j) = &
                            info_pheno(ws_idx(i,j))%wp_adj(soil_use%mat(i,j),crop_pars_mat%n_crop_in_year(i,j),y)
                    end if
                end if scans_domain
            end do!
        end do!
    end subroutine populate_crop_pars_matrices!
    
    subroutine calculate_RF_t(d_t, crop_par_mat,domain)!
        ! calculate root fraction in both evaporative and transpirative layer 
        real(dp), dimension(:,:), intent(in)::d_t
        type(grid_i),intent(in)::domain!
        type(crop_pars_matrices),intent(inout)::crop_par_mat
        
        ! populate pheno%RF_t & pheno%RF_e
        crop_par_mat%RF_t = merge(crop_par_mat%RF_t_max * (d_t/crop_par_mat%d_t_max)*&
                                 (1.0D0 / (crop_par_mat%RF_t_max*(d_t/crop_par_mat%d_t_max)+(1.0D0-crop_par_mat%RF_t_max))), &
                                  crop_par_mat%RF_t, crop_par_mat%RF_t_max /= domain%header%nan)
        crop_par_mat%RF_t = round_2darray(crop_par_mat%RF_t,6)
        crop_par_mat%RF_e = merge(1-  crop_par_mat%RF_t, crop_par_mat%RF_e, crop_par_mat%RF_t_max /= domain%header%nan)

    end subroutine calculate_RF_t!
    

end module mod_crop_phenology
