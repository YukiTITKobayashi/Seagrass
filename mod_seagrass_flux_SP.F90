!--------------------------------------------------------------------------------
!  mod_seagrass_flux_SP.F90
!  Module about Seagrass flux
!  The version considering cyrindrical flow and NOT considering volume of seagarass. (SimPle)
!--------------------------------------------------------------------------------

!initialize_seagrass (初期化にのみHofler-Thoday diagram利用)
!
!



#include "cppdefs.h"

MODULE mod_seagrass_flux
    
    USE mod_param
    USE mod_input
    USE mod_seagrass

    implicit none

    ! -----流束密度関連-----
    real(8) :: V_sav         ![m3] Volume of seagrass(width*thickness*height*shoot)

    real(8) :: psi_wc        ![Pa] Water Potential of Water Column
    real(8) :: psi_pw        ![Pa] Water Potential of Pore Water
    real(8) :: psi_sav_total       ![Pa] Water Potential in seagrass
    real(8) :: psi_sav_p     ![Pa] Pressure Component of seagrass (Turgor Pressure)
    real(8) :: psi_sav_pi    ![Pa] Osmotic Potential in seagrass
    real(8) :: theta_sav         !(ratio) Water Content
    real(8) :: osm_sav           !(PSU) Concentration of osmolyte which was converted to PSU

    !!-----E. Acoroides (後で構造体に入れ込む)
    real(8) :: sigma = 0.1d0          !反発係数，どれだけ圧力勾配通りに水が細胞壁を通過するか（一部は通過しない）(仮の値)
    real(8) :: rv_leaf = 1.0d-5        !1.0d-4        !管流を仮定した葉部通水組織の半径[m](=10um)
    real(8) :: rv_rhizome = 1.0d-5        !1.0d-4        !管流を仮定した地下茎通水組織の半径[m](=10um)
    
    real(8) :: width_leaf = 1.5d-2
    real(8) :: thickness_leaf = 5.0d-3
    real(8) :: l_leaf = 1.0d0        !0.5           !維管束の長さ[m]
    real(8) :: l_rhizome = 0.5d0        !0.5
    real(8) :: R_rhizome = 7.5d-3                   !地下茎の半径[m](=7.5mm)

        !チャネルを通過するイオン類の割合
    real(8) :: ch_leaf_in     = 0.01d0
    real(8) :: ch_leaf_out    = 0.01d0
    real(8) :: ch_rhizome_in  = 0.01d0
    real(8) :: ch_rhizome_out = 0.01d0
    !!-----

    real(8) :: eta = 1.25d-3           !海水の粘性係数[kg m-2 s-1](or [Pa*s])(15degC)
    real(8) :: Lp_leaf                 !葉部通水コンダクタンス[m s-1 Pa-1]
    real(8) :: Lp_rhizome                 !葉部通水コンダクタンス[m s-1 Pa-1]
    real(8) :: Jv_leaf                 !流出入流束密度[m3 m-2 s-1]([m s-1])
    real(8) :: Jv_rhizome                 ![m3 m-2 s-1]([m s-1])
    real(8) :: dV_leaf_dt              ![m3 s-1] 流量=単位時間あたり流出入体積
    real(8) :: dV_rhizome_dt
    real(8) :: osm_mass                ![PSU*m3] オスモライト溶質の濃度*体積の相当量

    !old functions
    real(8) :: dV_leaf                 ![m3] dt後の流出入体積
    real(8) :: dV_rhizome

    !  For test
    real(8) :: dV_sum
    
    !SP
    real(8) :: Jv_upward               !Pore Water --> Water Column方向の流束密度[m3 m-2 s-1]([m s-1])
    real(8) :: Lp_SAV


CONTAINS

    SUBROUTINE initialize_seagrass_flux

        return

    END SUBROUTINE initialize_seagrass_flux


    SUBROUTINE seagrass_flux

        !!!!!-----流束密度関連-----
        
        !塩分から水ポテンシャルを求める
        psi_wc  = WaterPotentialFromSal(salinity_wc)
        psi_pw  = WaterPotentialFromSal(salinity_pw)
        
        ! 浸透圧差から通水コンダクタンス求める(etaがtimeで変わることも考えてこの位置)
        Lp_leaf = rv_leaf**2/(8.0d0*eta*l_leaf)       !通水コンダクタンスを導出(Hagen-Poseuille流れ)
        Lp_rhizome = rv_rhizome**2/(8.0d0*eta*l_rhizome)    
        Lp_sav = Lp_leaf * Lp_rhizome / (Lp_leaf + Lp_rhizome)
        !流束密度
        Jv_upward = Lp_sav * sigma * (psi_pw - psi_wc)

    END SUBROUTINE seagrass_flux


    SUBROUTINE update_seagrass_flux(dt_)

        real(8), intent(in) :: dt_

        return

    END SUBROUTINE update_seagrass_flux



    REAL(8) FUNCTION WaterPotentialFromSal(Sal_)

        USE mod_param
        implicit none

        real(8), intent(in) :: Sal_      !Salinity(PSU)


        real(8) :: rho_
        real(8) :: C
        real(8) :: i_vanthoff
        real(8) :: psi_pi_

        !M[g/mol],R_gas[J K-1 mol-1]
        ! モル濃度を求める
        rho_ = 997.0d0 + Sal_                    ![kg m-3]=[g L-1]
        C = Sal_/M_sea / ((1.0d3/rho_)*1.0d-3)   ![mol/m3]
        ! Van't-Hoffの式
        i_vanthoff = 1.9d0 !Van't Hoff Coefficient
        psi_pi_ = -C * i_vanthoff * R_gas * (temp+273.15d0)  ![Pa]
        WaterPotentialFromSal = psi_pi_
        !print *, psi_pi_

    END FUNCTION WaterPotentialFromSal


END MODULE mod_seagrass_flux