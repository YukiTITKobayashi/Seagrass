#include "cppdefs.h"

    MODULE mod_output

        implicit none

    CONTAINS
    
! **********************************************************************
!  Open files
! **********************************************************************

        SUBROUTINE files_open
        
            implicit none

#if defined SGFLUX_HOFLER_THODAY
            open(11,file='./output/kob_test_36sec_HT.csv')
#elif defined SGFLUX_NATURAL
            open(11,file='./output/kob_test_36sec_NT.csv')
#elif defined SGFLUX_SIMPLE
            open(11,file='./output/kob_test_36sec_SP.csv')
#endif

            RETURN

        END SUBROUTINE files_open

        SUBROUTINE files_close
        
            implicit none

            close(11)

            RETURN

        END SUBROUTINE files_close

        SUBROUTINE write_env_label(fid)

            implicit none

            integer, intent(in) :: fid

            ! write(fid,*) 'time, ', 'short_radi, ', 'salinity_wc, ', 'salinity_pw, ', 'osm_sav, ', 'tide'

            !導出順的にはosm_savはJv_leafなどから導かれた値であるため最後．(時間差分法的にそうせざるを得ない)

            write(fid,*) 'time,      ' &
#if defined SGFLUX_HOFLER_THODAY
                    &  , 'day,                        ', 'short_radi,                 ', 'salinity_wc,                ' &
                    &  , 'salinity_pw,                ', 'sal_equiv_sav_total,        ', 'osm_sav,                    ' &
                    &  , 'tide,                       ' &
                    &  , 'psi_wc,                     ', 'psi_pw,                     ', 'psi_sav_total,              ' &
                    &  , 'theta_sav,                  ', 'psi_sav_p,                  ', 'psi_sav_pi,                 ' &
                    &  , 'Jv_leaf,                    ', 'Jv_rhizome,                 ', 'dV_sum,                     ' &
#elif defined SGFLUX_NATURAL
                    &  , 'day,                        ', 'short_radi,                 ', 'salinity_wc,                ' &
                    &  , 'salinity_pw,                ', 'sal_equiv_sav_total,        ', 'osm_sav,                    ' &
                    &  , 'tide,                       ' &
                    &  , 'psi_wc,                     ', 'psi_pw,                     ', 'psi_sav_total,              ' &
                    &  , 'theta_sav,                  ', 'psi_sav_p,                  ', 'psi_sav_pi,                 ' &
                    &  , 'Jv_leaf,                    ', 'Jv_rhizome,                 ', 'dV_sum,                     ' &
#elif defined SGFLUX_SIMPLE
                    &  , 'day,                        ', 'short_radi,                 ', 'salinity_wc,                ' &
                    &  , 'salinity_pw,                ', 'tide,                       ', 'psi_wc,                     ' &
                    &  , 'psi_pw,                     ', 'Jv_upward,                  ' &
#endif
                    &  , 'Pg,                         ', 'R,                          ' &
                    &  , 'DIC_sav,                    ', 'DIN_sav,                    ', 'DO_sav'




            RETURN

        END SUBROUTINE write_env_label



        SUBROUTINE write_env_data(fid)
            
            USE mod_param
            USE mod_seagrass_flux
            USE mod_seagrass

            implicit none

            integer, intent(in) :: fid
            character(len=10) :: disp
            !write(*,fmt='(A)',advance='no') "time: "
            !write(*,fmt='(A)',advance='no') disp
            !write(*,fmt='(A)') "       "

            disp = time_display(time,start_time)

            write(fid,*) disp,',', time,',', rad,',',salinity_wc,',',salinity_pw &
#if defined SGFLUX_HOFLER_THODAY
            & ,',', sal_equiv_sav_total,',',osm_sav,',', tide &
            & ,',', psi_wc,',', psi_pw,',', psi_sav_total,',', theta_sav,',', psi_sav_p,',', psi_sav_pi &
            & ,',', Jv_leaf,',', Jv_rhizome,',', dV_sum &
#elif defined SGFLUX_NATURAL
            & ,',', sal_equiv_sav_total,',',osm_sav,',', tide &
            & ,',', psi_wc,',', psi_pw,',', psi_sav_total,',', theta_sav,',', psi_sav_p,',', psi_sav_pi &
            & ,',', Jv_leaf,',', Jv_rhizome,',', dV_sum &
#elif defined SGFLUX_SIMPLE
            & ,',', tide      &
            & ,',', psi_wc,',', psi_pw &
            & ,',', Jv_upward &
#endif
            ,',', Pg,',', R,',', DIC_sav,',', DIN_sav,',', DO_sav

            RETURN
    
        END SUBROUTINE write_env_data


        CHARACTER(len=10) FUNCTION time_display(time,start_time)

            implicit none

            real(8), intent(in) :: time
            real(8), intent(in) :: start_time

            real(8) now
            integer day
            integer hour
            integer minute
            integer second

            now = start_time + time
            day    = int(now)
            hour   = int((now - day) * 24.e0)
            minute = int((now - day - hour/24.e0) * 24.e0 * 60.e0)
            second = int((now - day - hour/24.e0 - minute/24.e0/60.e0) * 24.e0 * 60.e0 * 60.e0)

            write(time_display(1:1),'(I1.1)') day
            write(time_display(2:2),'(A)') "D"
            write(time_display(3:4),'(I2.2)') hour
            write(time_display(5:5),'(A)') ":"
            write(time_display(6:7),'(I2.2)') minute
            write(time_display(8:8),'(A)') ":"
            write(time_display(9:10),'(I2.2)') second
            !print *, time_display

        END FUNCTION



        SUBROUTINE write_env_data_old(fid,disp_,rad_,salinity_wc_,salinity_pw_,osm_sav_,tide_, &
            & psi_wc_,psi_pw_,psi_sav_,theta_sav_,psi_sav_p_,psi_sav_pi_,Jv_leaf_,Jv_rhizome_)
            
            USE mod_param

            implicit none

            integer, intent(in) :: fid
            character(len=10), intent(in) :: disp_
            real(8), intent(in) :: rad_
            real(8), intent(in) :: salinity_wc_
            real(8), intent(in) :: salinity_pw_
            real(8), intent(in) :: osm_sav_
            real(8), intent(in) :: tide_
            real(8), intent(in) :: psi_wc_
            real(8), intent(in) :: psi_pw_
            real(8), intent(in) :: psi_sav_
            real(8), intent(in) :: theta_sav_
            real(8), intent(in) :: psi_sav_p_
            real(8), intent(in) :: psi_sav_pi_
            real(8), intent(in) :: Jv_leaf_
            real(8), intent(in) :: Jv_rhizome_


            

            write(11,*) disp_,',', rad_,',',salinity_wc_,',',salinity_pw_,',', osm_sav_,',', tide_      &
            & ,',', psi_wc_,',', psi_pw_,',', psi_sav_,',', theta_sav_,',', psi_sav_p_,',', psi_sav_pi_ &
            & ,',', Jv_leaf_,',', Jv_rhizome_

            RETURN
    
        END SUBROUTINE write_env_data_old



    END MODULE mod_output