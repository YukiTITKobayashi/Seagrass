! **********************************************************************
! *                                                                    *
! *   Test program of mod_reef_ecosys                                  *
! *                                                                    *
! **********************************************************************
!現在の問題：要修正で検索
!なぜかinputの値が10行ほどずれている，updateのsubroutineはreef_ecosysにはない


#include "cppdefs.h"


    PROGRAM kob_test

        USE mod_param
        USE mod_input
        USE mod_output
        USE mod_seagrass
        USE mod_seagrass_flux


        implicit none

        !real(8), parameter :: dt = 0.1    ![hour] 6min
        real(8), parameter :: dt = 0.01    ![hour] 36sec

        integer :: istep, iprint, iend
        integer :: kSetting

        


!  For Output      
        !real(8), parameter :: OUTPUT_INTERVAL = 5.0d0     ! Output interval (min)


!----- Open output files -------------------------

        CALL files_open





!----- Set initial conditions -------------------------

        time = 0.0d0             !day (from 2/22 10:00 --> 2/23 10:00: time = 24.0d0/24.0d0)
        start_time = 10./24.
        start_day = 31 + 22      !day (Feb. 22)

        end_time = 10./24.
        end_day = 31 + 29        !day (Feb. 23)
        !end_time = 11./24.
        !end_day = 31 + 22        !day (Feb. 22)

        term = (end_day + end_time) - (start_day + start_time)          ![day]

!Setting of condition (kSetting)
!
!       kSetting = 1: normal
        kSetting = 1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !----- Import data -----------------------------------
        CALL read_timeseries_tide
        CALL read_timeseries_others    
        !-----------------------------------------------------    

        istep=0
        iprint=0

        CALL setdata(kSetting)

        CALL initialize_seagrass
        CALL initialize_seagrass_flux


!----- Write data labels -------------------------
        CALL write_env_label(11)

!----- Main loop -------------------------------------------


!      do istep=1, int(60.*60./dt)+1       ! 1 hour
!      do istep=1, int(24.*60.*60./dt) * 14 +1      ! 5 days
!      do istep=1, int(24.*60.*60./dt) * 7 +1      ! 7 days
!      do istep=1, int(24.*60.*60./dt) * 9 +1      ! 9 days
!      do istep=1, int(24.*60.*60./dt) * 30 +1      ! 30 days
!      do istep=1, int(24.*60.*60./dt) *365*3 +1      ! 3 year
        iend = int((24./dt) * term) + 1                                      ! 1 day
        do istep=1, iend                                                ! 1 day
!        do istep=1, int(1./dt) * 1 +1         ! 1 hour
!       do istep=1, iend
            !time = 0.3125

!------ Set environmental parameters ----------------------------

            CALL setdata(kSetting)

!----- Ecosystem model ----------------------------------------

            CALL seagrass_flux

!----- light attenuation -----
            CALL seagrass(dt)
            !Pg = 0.0d0
            !I_d_ave = 0.0d0
            !do leaf_pos = 1, part_leaf
                !print *, SGRASS(Nsg)%Pg
                !I_d = PFDsurf * exp( - attenuation_coeff * (-bathym + tide - SGRASS(Nsg)%Leaf%length) )
                !if (-bathym + tide - leaf_pos * div_leaf >= 0) then
                        !I_d(leaf_pos) = PFDsurf * exp( - attenuation_coeff * (-bathym + tide - leaf_pos * div_leaf ) )
                !else
                        !I_d(leaf_pos) = PFDsurf
                !end if
                !I_d_ave = I_d_ave + I_d(leaf_pos)/part_leaf
                !Pg_part(leaf_pos) = pmax*tanh(I_d(leaf_pos)/pIk)/3600.d0   !Light response curve [mmolC/m2/s]
                !Pg = Pg + Pg_part(leaf_pos)/part_leaf
            !end do
            
            !R  = p0/3600.d0   !Constant [mmolC/m2/s]

            !print *, Pg_part(60)
            !DICuptake= Pg - R
            !DOuptake = R - Pg

            !IF(DICamb<=0.d0 .or. DINamb <=0.d0 .or. DIPamb <= 0.d0) THEN !-----For Error handling
            !    Pg = 0.d0
            !ENDIF
            !IF(DOamb<=0.d0) THEN !-----For Error handling
            !    R = 0.d0
            !ENDIF

!  -- Print --------------------------------------------
            
            !print *, tide
            
            !print '(f5.2)', salinity_wc
            !print *, salinity_pw

            !write(*,fmt='(A)',advance='no') "I_d: "
            !print *, I_d
            
            
            CALL write_env_data(11)             !初期条件なども踏まえてアップデートの前の位置に置く方が適する


            CALL update_seagrass_flux(dt)       !dV_leaf, dV_rhizome, osm_savのアップデート

            time=time+dt/24.


            !CALL write_env_data(11)            !reef_ecosysに倣うならこの位置
            

        end do

        !print *, tide_data
        !print *, salinity_data_pw
        !print *, salinity_pw
        print *, dV_sum
        
        CALL files_close
 

    END PROGRAM kob_test


