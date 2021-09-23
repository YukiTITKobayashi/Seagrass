!--------------------------------------------------------------------------------
!
!  Seagrass module ver. 1
!
!--------------------------------------------------------------------------------

#include "cppdefs.h"

    MODULE mod_seagrass

        USE mod_seagrass_flux
        
        implicit none

        

        !!-!!--!-----For metabolic processes-----
        
        !integer, parameter :: Nsg = 1    !! Number of seagrass groups

        !TYPE T_TISSUE_S
            
            !real(8), pointer :: length
            !real(8), pointer :: surface_area
            !real(8), pointer :: volume
            !real(8), pointer :: xylem_radius
            !real(8), pointer :: elasticity   !of cell wall
            !real(8), pointer :: Lp           !water flow conductance [m s-1 Pa-1]
            !real(8), pointer :: Jv           !Volume flux density [m3 m-2 s-1]

            !Element is (time) --> need much memory thus not used
            !These values differ by the time
            !real(8), pointer :: psi        !water potential[MPa]
            !real(8), pointer :: psi_P      !pressure potential[MPa]
            !real(8), pointer :: psi_pi     !osmotic potential[MPa] (including matrix potential)
            !real(8), pointer :: turgor     !turgor pressure[MPa]
            !real(8), pointer :: osmolality_total
            !real(8), pointer :: osmolality_non_vacuolar
            !real(8), pointer :: theta      !moisture content [%]
            !real(8), pointer :: W          !moisture content [m3]
            !real(8), pointer :: nut
            !real(8), pointer :: growth  ! growth rate

        !END TYPE T_TISSUE_S

        !TYPE T_SGRASS1
            !Elements are (time), which are different from mod_seagrass and mod_coral
            !real(8), pointer :: Pg
            !real(8), pointer :: R(:)
            !real(8), pointer :: QC(:)
            !integer, pointer :: num_leaves(:)
            !TYPE (T_TISSUE_S), pointer :: Leaf    !(number of leaves)
            !TYPE (T_TISSUE_S), pointer :: Sheath
            !TYPE (T_TISSUE_S), pointer :: Rhizome
            !TYPE (T_TISSUE_S), pointer :: rhizome
        !END TYPE T_SGRASS1

        !TYPE (T_SGRASS1), allocatable :: SGRASS(:)  !(number of species)
        !!-!--!!

        
        real(8) :: Leaf_length      ![m]
        integer :: div_leaf         ![m]
        integer :: part_leaf
        real(8), allocatable :: I_d(:)              ![umol m-2 s-1]
        real(8) :: I_d_ave
        real(8), allocatable :: Pg_part(:)          ![mmolC/m2/s]
        real(8) :: Pg               ![mmolC/m2/s]
        real(8) :: R                                ![mmolC/m2/s]

        real(8) :: R_dark
        real(8) :: R_light

        real(8), parameter :: g_R_per_Pg = 0.2d0    !coefficient of activeness for light-attenuated respiration


        real(8) :: PFD
        real(8) :: rho_sw = 1000.0d0        !kg/m3
        real(8) :: DIC_wc = 1900.0d0        !umol/kg
        !real(8) :: DO_wc                  !umol/L (from mod_input)
        real(8) :: DIN_wc = 10.0d0          !umol/L
        !real(8) :: DIP_wc = 1.0d0

        real(8) :: DIC_pw = 1900.0d0
        real(8) :: DO_pw  = 0.5d0
        real(8) :: DIN_pw = 40.0d0
        !real(8) :: DIP_pw


  ! output parameters
        real(8) :: DICuptake
        real(8) :: DOuptake
       
        real(8) :: DINuptake
        real(8) :: DIPuptake

        !!!!!!!!!!!!!!!!!!
        real(8) :: DIC_sav
        real(8) :: DIN_sav
        real(8) :: DO_sav



        ! --- C:N:P ratio of seagrass ---
      real(8), parameter :: nc=27./599.d0 !M.J.Atkinson and SV Smith(1983)
      real(8), parameter :: pc=1./599.d0
! --- Photosynthesis  Parameters ---
!      real(8), parameter :: pmax =  51.3d0 ! Watanabe et al. 2013
!      real(8), parameter :: pIk  = 589.65d0
!      real(8), parameter :: p0   =  15.05d0
      real(8), parameter :: pmax =  55.81d0  ! Nakamura & Nakamori 2009
      real(8), parameter :: pIk  = 671.8d0   !  Model skill = 0.990
      real(8), parameter :: p0   =  21.62d0  !

      !Thermal dependence
      real(8), parameter :: dH = 2803.0d3/6.0d0     ![J mol-1] 化学式から求めた反応熱
      real(8)            :: dS = dH/340.0d0         !エントロピー変化量 !340Kで活性が1/(1+A)倍になるためのエントロピー
      real(8), parameter :: Ea = 51.0d3             ![J mol-1] 活性化エネルギー（水？）
      real(8)            :: A_Arrhenius = 1/exp(-Ea/(R_gas*303.0d0))            
      !ArrheniusのAmplitude !25degCで1になるようなA(本当はArrhenius plotで導出)


      real(8) :: dp_I_d     !light dependence
      real(8) :: dp_temp    !thermal depencence
      real(8) :: dp_nutrient
      real(8) :: dp_DO

      
      real(8), parameter :: K_DIN = 15.0d0             ![mmol m-3]=[uM] (Touchette, Burkholder (2000)からseagrassについて大体の値)
      real(8) :: K_DIC = K_DIN / nc
      real(8) :: K_DO  = 3.0d0                        !適当な値

      integer :: leaf_pos



    CONTAINS

        !SUBROUTINE initialize_seagrass(Nsg)
        SUBROUTINE  initialize_seagrass

            implicit none

            !integer, intent(in) :: Nsg

            !allocate (SGRASS(Nsg))
            !SGRASS(Nsg)%Pg = 1000
            !SGRASS(Nsg)%Leaf%length = 0.8d0    ![m]
            !SGRASS(Nsg)%rhizome%length = 0.5d0    ![m]


            ! -----光合成・呼吸関連（書き途中）-----

            Leaf_length = 0.8d0
            part_leaf = 8
            div_leaf = Leaf_length/part_leaf

            DIC_sav = 2000.0d0  !(umol/kg)
            DIN_sav = 15.0d0    !(umol/L)
            DO_sav  = 4.0d0    !(umol/L)

            allocate(I_d(part_leaf))
            allocate(Pg_part(part_leaf))
  

        END SUBROUTINE initialize_seagrass

        SUBROUTINE seagrass(dt_)

            USE mod_param

            implicit none

            real(8), intent(in) :: dt_


            !----- light attenuation -----
            Pg = 0.0d0
            I_d_ave = 0.0d0

            !要修正：reef_ecosysの(i,j,n)と構造体定義がよくわからず，leaf_posを使っている

            dp_temp = A_Arrhenius * exp(-Ea/(R_gas*(temp+273.15d0))) &
            &         / (1.0d0 + A_Arrhenius * exp(((temp+273.15d0)*dS-dH)/(R_gas*(temp+273.15d0))))
            !Thermal dependence (簡易的にT=const.にしているが，パラメータはイメージで調整できている．)
            dp_nutrient = min(DIC_sav/(K_DIC + DIC_sav), DIN_sav/(K_DIN + DIN_sav))
            !Nutrient dependence DICとDINの積で表すと値が小さくなりすぎるのがよくなかったのでリービッヒの最小律に基づく

            do leaf_pos = 1, part_leaf
                !print *, SGRASS(Nsg)%Pg
                !I_d = PFDsurf * exp( - attenuation_coeff * (-bathym + tide - SGRASS(Nsg)%Leaf%length) )
                if (-bathym + tide - leaf_pos * div_leaf >= 0) then
                        I_d(leaf_pos) = PFDsurf * exp( - attenuation_coeff * (-bathym + tide - leaf_pos * div_leaf ) )
                else
                        I_d(leaf_pos) = PFDsurf
                end if
                I_d_ave = I_d_ave + I_d(leaf_pos)/part_leaf
                !Pg_part(leaf_pos) = pmax*tanh(I_d(leaf_pos)/pIk)/3600.d0   !Light response curve [mmolC/m2/s]

                !DEPENDENCE
                !dp_I_d  = tanh(I_d(leaf_pos)/pIk)
                dp_I_d  = I_d(leaf_pos)/(0.4d0 * pIk + I_d(leaf_pos))    !Light dependence
                
                Pg_part(leaf_pos) = pmax * dp_I_d * dp_temp * dp_nutrient /3600.d0   !Light response curve [mmolC/m2/s]
                !(要修正：tanh(E/Ek)とE/(Ek+E)が近くなるように0.4d0)
                Pg = Pg + Pg_part(leaf_pos)/part_leaf
            end do

            !print *, Pg
            
            !R_dark  = p0/3600.d0   !Constant [mmolC/m2/s]
            R_dark  = p0/10000.d0   !Constant [mmolC/m2/s]

            R_light = g_R_per_Pg  * Pg
            !print *, R_light

            dp_DO = DO_sav / (K_DO + DO_sav)

            R = dp_DO * dp_temp * (R_dark + R_light)

            !print *, Pg_part(60)
            !DICuptake= Pg - R
            !DOuptake = R - Pg

            IF(DIC_wc<=0.d0 .or. DIN_wc <=0.d0 .or. DIC_pw <= 0.d0 .or. DIN_pw <= 0.d0) THEN !-----For Error handling
                Pg = 0.d0
            ENDIF
            IF(DO_wc <=0.d0 .or. DO_pw <= 0.d0) THEN !-----For Error handling
                R = 0.d0
            ENDIF

            !  -- Print --------------------------------------------
            
            !print *, tide
            
            !print '(f5.2)', salinity_wc
            !print *, salinity_pw

            !write(*,fmt='(A)',advance='no') "I_d: "
            !print *, I_d

            !  -- Update ------------------------------------------
            IF (Jv_leaf >= 0 .and. Jv_rhizome >= 0) THEN
                DICuptake = Jv_leaf * DIC_wc + Jv_rhizome * DIC_pw  ![m s-1 * umol kg-1]
                DINuptake = Jv_leaf * DIN_wc + Jv_rhizome * DIN_pw  ![umol m-2 s-1]=[m s-1 * umol m-3]
                DOuptake  = Jv_leaf * DO_wc  + Jv_rhizome * DO_pw
            ELSE IF (Jv_leaf >= 0 .and. Jv_rhizome < 0) THEN
                DICuptake = Jv_leaf * DIC_wc  
                DINuptake = Jv_leaf * DIN_wc
                DOuptake  = Jv_leaf * DO_wc
            ELSE IF (Jv_leaf < 0 .and. Jv_rhizome >= 0) THEN
                DICuptake = Jv_rhizome * DIC_pw 
                DINuptake = Jv_rhizome * DIN_pw
                DOuptake  = Jv_rhizome * DO_pw
            ELSE IF (Jv_leaf < 0 .and. Jv_rhizome < 0) THEN
                DICuptake = 0.0d0  
                DINuptake = 0.0d0
                DOuptake  = 0.0d0
            ENDIF


            !要修正：SPの場合は合っているが，HT,NTはモデルの流入部断面積が管のように一定ではないから
            !/(l_leaf+l_rhizome)ではなく，*S/(theta*V)する必要がある
            !*rho/1000が隠れている
            DIC_sav = DIC_sav + ( DICuptake/1000 - Pg + R )  /theta_sav/(l_leaf+l_rhizome) * dt_ *3600!umol kg-1
            DIN_sav = DIN_sav + ( DINuptake - nc * Pg ) /theta_sav/(l_leaf+l_rhizome) * dt_ *3600
            DO_sav  = DO_sav  + ( DOuptake  + Pg - R )  /theta_sav/(l_leaf+l_rhizome) * dt_ *3600

            !要修正：飽和状態に近い時，uptakeが穏やかになることを考えていない
            !要修正：CO2の飽和溶解度31.22d3 mmol m-3?合わないからDIC=2500 umol m-3未満にする
            if (DIC_sav > 2500.0d0) then
                DIC_sav = 2500.0d0
            end if
            if (DIN_sav > 50.0d0) then
                DIN_sav = 50.0d0
            end if
            if (DO_sav > 8.0d0) then
                DO_sav = 8.0d0
            end if

        END SUBROUTINE seagrass
  

    END MODULE mod_seagrass