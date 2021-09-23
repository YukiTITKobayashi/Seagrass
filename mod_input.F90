!--------------------------------------------------------------------------------
!
!              Input module
!mod_input.F90だけ単位を時間にしてしまっているので後で直す．(lin_interpol3も含めて)
!--------------------------------------------------------------------------------

#include "cppdefs.h"


    MODULE mod_input


    CONTAINS

        SUBROUTINE read_timeseries_tide

            USE mod_param
            USE mod_fileio
            
            implicit none

            integer :: nm_data, ios
            integer :: i,j

!
! ----- READ Tide data (m) ---------------------------------------------
!
#if defined TIDE_DATA_FILE

            open(77,file='./input/tide_data.dat')   !Any row * 1 column array per dt_tide hours
            dt_tide=1.0d0  ![hour] interval

!  ---- count data number -----
            nm_data = 0
            do
                read(77,*,iostat=ios)
                if(ios==-1) exit
                nm_data = nm_data+1
            end do
            allocate( tide_data(nm_data) )
            rewind(77)
! ---- read data -----
            do i =1, nm_data
                read(77,*) tide_data(i)
            end do
            close(77)

#elif defined TIDE_CSV_FILE

            CALL csv_read_char(77,"./input/naotide1_append.csv",tide_csv,1,nm_data,column_tide)     !From module: mod_fileio
            dt_tide=1.0d0  ![hour] interval
            nm_tide = nm_data

            allocate( tide_data(nm_data) )

            read(tide_csv(:,2),*) tide_data(:)   !文字列型2次元配列から実数型1次元配列を抽出
            tide_data(:) = tide_data(:) / 100d0
                
#elif defined TIDE_FORMULA1
            dt_tide = 0.1d0   !6min
            day_ex  = term    ![day] period of simulation
            h_max   = 0.6d0     ![m] maximum tide 
            ave_tide_target = - 0.3d0     ![m] INITIAL & average tide (should be negative)
            nm_data = day_ex * 24.0d0 / dt_tide + 1.0d0
            allocate( tide_data(nm_data) )
            nm_tide = nm_data
            do i =1, nm_data
                tide_data(i) = h_max * sin( 2.0d0 * pi * (i * dt_tide / 24.0d0) ) + ave_tide_target
            end do

#elif defined TIDE_FORMULA2
            dt_tide = 0.1d0   !6min
            day_ex  = term    ![day] period of simulation
            h_max   = 0.6d0     ![m] maximum tide 
            ave_tide_target = - 0.0d0     ![m] INITIAL & average tide (should be negative)
            arg_tide = 0.5d0        ![rad/pi] argument of tide
            nm_data = day_ex * 24.0d0 / dt_tide + 1.0d0
            allocate( tide_data(nm_data) )
            nm_tide = nm_data
            do i =1, nm_data
                tide_data(i) = h_max * sin( pi * (4.0d0 * (i * dt_tide / 24.0d0) + arg_tide) ) + ave_tide_target
                !tide_data(i) = h_max * sin( pi * (4.0d0 * (i * dt_tide / 24.0d0) + 0.5d0) ) + ave_tide_target
            end do            
#endif

        END SUBROUTINE read_timeseries_tide


        SUBROUTINE read_timeseries_others

            USE mod_param
            USE mod_fileio

            implicit none

            integer :: nm_data, ios
            integer :: i,j

!
! ----- READ ambient water column & pore water salinity data (psu) ---------------------------------------------
!

#if defined SAL_CSV_FILE
            ! Busuanga_salinity_point_append.csvは2/22 10:00のsalinity_wc,salinity_pwに11:00の平均値を追加している
            CALL csv_read_char(77,"./input/Busuanga_salinity_point_append.csv",salinity_csv,1,nm_salinity,column_salinity)     !From module: mod_fileio
            dt_salinity=1.0d0  ![hour] interval

            allocate( salinity_data_wc(nm_salinity), salinity_data_pw(nm_salinity) )
            read(salinity_csv(:,2),*) salinity_data_wc(:)   !water column
            read(salinity_csv(:,3),*) salinity_data_pw(:)   !pore water

#elif defined SAL_FORMULA_FROM_TIDE
            A_S_wc = 7.5d0 / h_max      !12.5 (h_max=0.6)
            A_S_pw = 1.5d0 / h_max      !2.5 (h_max=0.6)
            ave_sal_wc_target = 27.5d0
            ave_sal_pw_target = 32.0d0

            dt_salinity = dt_tide
            nm_salinity = nm_tide
            allocate( salinity_data_wc(nm_salinity), salinity_data_pw(nm_salinity) )
            salinity_data_wc(:) = A_S_wc * tide_data(:) + ave_sal_wc_target
            salinity_data_pw(:) = A_S_pw * tide_data(:) + ave_sal_pw_target
#endif

        END SUBROUTINE read_timeseries_others


!
! **********************************************************************
!  Set environmental condition
! **********************************************************************
!
        SUBROUTINE setdata(kSetting_)
!     Setting of condition (kSetting)
!
!     nSetting = 1: normal

            USE mod_param

            implicit none
                    
            integer, intent(in) :: kSetting_

            IF (kSetting_ .eq. 1) then
!  -- Set Tide data ---------------------------------------------
                tide        = lin_interpol3(time,tide_data,dt_tide,nm_tide)

!  -- Set other data ---------------------------------------------
!  -- Salinity -----------  ----------------------------------
                salinity_wc = lin_interpol3(time,salinity_data_wc,dt_salinity,nm_salinity)
                salinity_pw = lin_interpol3(time,salinity_data_pw,dt_salinity,nm_salinity)

!  -- Solar radiation --------------------------------------------
                rad = short_radi(time,start_time,start_day,lat,temp,Hum,2)
                !write(*,fmt='(A)',advance='no') "short_radi: "
                !print *, rad
                !write(*,fmt='(A)',advance='yes') " "
                
                PFDsurf = rad * 4.57   ![umol m-2 s-1] <-- [W m-2] (under sunny env.)
                !write(*,fmt='(A)',advance='no') "PFDsurf: "
                !print *, PFDsurf

                !Testing!!!!!
                DO_wc = PFDsurf / 150.0d0 + 1
            
            END IF

            return

        END SUBROUTINE setdata
            


        REAL(8) FUNCTION lin_interpol3(time,dataset,dt_data,datamax)
        ! **********************************************************************
        ! *                                                                    *
        ! * FUNCTION    :  linear interpolation between input data.            *
        ! *              time:     progress time (day)                         *
        ! *              dataset(i): data set                                  *
        ! *              dt_data: data interval (hour)                         *
        ! * !!!!Editted to make a linear interpolation by the unit of hours    *
        ! **********************************************************************
        !
            implicit none
        
            integer i
        !      real lin_interpol
            real(8) dataset(datamax)
            real(8) time,dt_data
            integer datamax
        
            i=int(time*24.e0/dt_data)+1     !Accurate decimal expression is necessary with rounding down
            !timestep = time*24.e0/dt
            !datastep = dt_data/dt
            !i = int(timestep/datastep) + 1
              
            !print *, time

            if (i.ge.datamax) i=datamax-1
        
            lin_interpol3 =dataset(i)+(dataset(i+1)-dataset(i))      &
           &       *(time*24.e0/dt_data-real(i-1))

            return
        
        END FUNCTION lin_interpol3
            
! *****function solar_radi (not used)***********************************


        

        real(8) function solar_radi(time,PFD,start_time)
        ! **********************************************************************
        ! *                                                                    *
        ! * FUNCTION    :  caliculate solar radiation (W/m2).                  *
        ! *              time:       progress time (day)                       *
        ! *              start_time: start time (day)                          *
        ! *                 ex.) 0:00-> start_time=0; 12:00-> start_time=0.5   *
        ! *                                                                    *
        ! **********************************************************************
    
            implicit none
            
            real(8), parameter :: pi = 3.141592654d0
    
            real(8) time, PFD,start_time
    
        !      solar_radi=PFD/4.57*sin(2.*pi*(time+start_time-0.25))
            solar_radi=PFD/1.82d0*sin(2.*pi*(time+start_time-0.25))
                    !convert photon flux density (umol m-2 s-1) to solar radiation (W m-2)
    
            if(solar_radi.lt.0.) then
            solar_radi=0.
            endif
    
            return
    
        end function solar_radi

! **********************************************************************
        
        
        real(8) function short_radi(time,start_time,start_day, lat   &
            &                         ,temp,Hum, mode )
        ! **********************************************************************
        ! *                                                                    *
        ! * FUNCTION    :  Caliculate short wave radiation (W/m2).             *
        ! *                 by Zillman equation                                *
        ! *                                                                    *
        ! *              time:       progress time (day)                       *
        ! *              start_time: start time (day)                          *
        ! *                 ex.) 0:00-> start_time=0; 12:00-> start_time=0.5   *
        ! *              start_day:  start day (day)                           *
        ! *                 ex.) Jan. 1-> start_day=1; Jun 22-> start_day=174  *
        ! *                                                                    *
        ! *              lat: latitude (degree)                                *
        ! *              Hum: rerative humidity (%)                            *
        ! *              temp: temperature (degree C)                          *
        ! *                                                                    *
        ! *              Hum: rerative humidity (%)                            *
        ! *                                                                    *
        ! *              mode: mode =1 -> day time progress                    *
        ! *                    mode =2 -> loop same day condition              *
        ! * (Added start_time)                                                 *
        ! **********************************************************************
        !
            implicit none
            
            real(8), parameter :: pi = 3.141592654d0
            real(8), parameter :: Sc = 1353.0d0        !Solar constant (W m-2)
    
            real(8) :: time, start_time,lat,temp, Hum
            integer :: mode
            real(8) :: start_day
            real(8) :: yDay !day of year
            real(8) :: cosZ
            real(8) :: DA    !DA=2.*pi*yDay/365.
            real(8) :: e_vep !vapor pressure (Pa)
            real(8) :: dec   !declination (radian)
            real(8) :: HA    !hour angle (radian)
            real(8) :: lati
            integer :: i
    
            e_vep=611.*10.**(7.5*(temp/(temp+273.15-35.86))) * Hum/100.
    
            HA=(0.5-time-start_time)*2.*pi
    
            if(mode .eq. 1) then
            yDay=start_day+Int(start_time + time) ! day time progress
            else
            yDay=start_day            ! light condition releated same day.
            endif
            
            yDay=mod(yDay-1.,365.)+1.
    
            dec=23.44*cos((172.-yDay)*2.*pi/365.) *pi/180.
    !      DA=2.*pi*yDay/365.
    !      dec=0.006918-0.399912*cos(DA)   +0.070257*sin(DA)
    !     &            -0.006758*cos(2.*DA)+0.000907*sin(2.*DA)
    !     &            -0.002697*cos(3.*DA)+0.001480*sin(3.*DA)
    !
            lati=lat *pi/180. !degree -> radian
    
            cosZ=sin(lati)*sin(dec)+cos(lati)*cos(dec)*cos(HA)
    
    
    
            if(cosZ.lt.0.) then
            short_radi=0.
            else
            short_radi=Sc*cosZ**2./((cosZ+2.7)*e_vep*1.e-5 + 1.085*cosZ +0.10)
            endif
    
    !----for debug---------------------------------------------------------c
    !        if(i.eq.7200) then
    !!        if(iprint.eq.1000) then
    !!        if(iprint.eq.1) then
    !         write(52,*) time,HA,yDay, dec,lati,cosZ,short_radi,e_vep
    !        endif
    !------------------------------------------------------------------c        
    
            return
    
        end function short_radi 

    END MODULE mod_input