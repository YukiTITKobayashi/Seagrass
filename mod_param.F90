!--------------------------------------------------------------------------------
!
!  Parameter module
!
!--------------------------------------------------------------------------------

#include "cppdefs.h"

    MODULE mod_param

        implicit none


! ----- mod_input.F90 --------------------------------------------------
        
        real(8) :: time                         ![day]

        real(8) :: start_time                   ![day]
        real(8) :: start_day                    ![day]

        real(8) :: end_time
        real(8) :: end_day

        real(8) :: term                         ![day]

        real(8), allocatable :: tide_data(:)    !imported data with any row * 1 column array per dt_tide hours
        real(8) :: tide                         !
        real(8) :: dt_tide                      ![hour] interval to import tide_data
        integer :: nm_tide                      !total row of tide_data     


#if defined TIDE_CSV_FILE
        character(50), allocatable :: tide_csv(:,:)
        integer :: column_tide

#elif defined TIDE_FORMULA1
        real(8) :: day_ex                       ![day] period of simulation
        real(8) :: h_max                        ![m] maximum tide
        real(8) :: ave_tide_target              ![m] average tide which is searched as target(should be negative)


#elif defined TIDE_FORMULA2
        real(8) :: day_ex                       ![day] period of simulation
        real(8) :: h_max                        ![m] maximum tide
        real(8) :: ave_tide_target              ![m] average tide which is searched as target(should be negative)
        real(8) :: arg_tide                     ![rad/pi] argument of tide

#endif
        !Salinity of water column
        real(8), allocatable :: salinity_data_wc(:)        !imported data with any row * 1 column array per dt_salinity hours
        real(8) :: salinity_wc

        !Salinity of pore water
        real(8), allocatable :: salinity_data_pw(:)        !imported data with any row * 1 column array per dt_salinity hours
        real(8) :: salinity_pw

        real(8) :: dt_salinity
        integer :: nm_salinity

#if defined SAL_CSV_FILE
        character(50), allocatable :: salinity_csv(:,:)
        integer :: column_salinity
#elif defined SAL_FORMULA_FROM_TIDE
        real(8) :: A_S_wc
        real(8) :: ave_sal_wc_target
        real(8) :: A_S_pw
        real(8) :: ave_sal_pw_target
#endif

        real(8) :: rad                  ![W m-2]
        real(8) :: PFDsurf              ![umol m-2 s-1]
        real(8) :: lat = 12.06

        real(8) :: temp = 27
        real(8) :: Hum = 50


! ----- mod_test.F90 --------------------------------------------------
        integer :: ng = 1
        real(8), parameter :: attenuation_coeff = 1.2d0         ![m-1]
        real(8) :: bathym = -1.5d0       ![m]


!-----Other scientific parameters-----
        real(8), parameter :: pi  = 3.14159265359d0             ! Circle ratio
        real(8), parameter :: sigma_water_air = 7.28d-2         ![Nm-1]Surface tension of water on air at 20degC
        real(8), parameter :: Cp_water = 4182.0d0               ![J kg-1 K-1] Heat capacity of water
        real(8), parameter :: lamda_water = 2.454d6             ![J kg-1] Vaporization heat of water at 20degC    
        real(8), parameter :: V_water = 18.05d-6                ![m-3 mol-1] Volume of water at 20degC
        real(8), parameter :: M_sea = 58.4*0.779 + 95.2*0.096 + 120.4*0.061 + 136.1*0.04 + 74.6*0.021    ![g/mol] Average molar mass of seawater
        real(8), parameter :: R_gas = 8.31446                   ![J K-1 mol-1] Gas constant

!!!For testing
        real(8) :: test
        real(8) :: DO_wc


    CONTAINS


    END MODULE mod_param

        