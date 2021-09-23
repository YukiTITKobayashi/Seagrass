!--------------------------------------------------------------------------------
!  mod_seagrass_flux_NT.F90
!  Module about Seagrass flux
!  The version NOT using Hofler-Thoday diagram. (NaTural)
!--------------------------------------------------------------------------------

!initialize_seagrass_flux (目安のため初期化にのみHofler-Thoday diagram利用)
!
!     S_wc == S_equiv_sav_total(t=0) --> psi_sav_total(t=0) --> theta_sav(t=0) 
!                                                          |--> psi_sav_pi(t=0) --> osm_sav(t=0) 
!                                                          |--> psi_sav_p(t=0)
!
!seagrass_flux
!
!     psi_sav_pi --> osm_sav
!     psi_sav_p  --> bulk_modulus -->
!
!update_seagrass_flux


#include "cppdefs.h"

MODULE mod_seagrass_flux
    
    USE mod_param
    USE mod_input
    !USE mod_seagrass

    implicit none

    ! -----流束密度関連-----
    real(8) :: V_sav         ![m3] Volume of seagrass(width*thickness*height*shoot)

    real(8) :: sal_equiv_sav_total

    real(8) :: psi_wc        ![Pa] Water Potential of Water Column
    real(8) :: psi_pw        ![Pa] Water Potential of Pore Water
    real(8) :: psi_sav_total       ![Pa] Water Potential in seagrass
    real(8) :: psi_sav_p     ![Pa] Pressure Component of seagrass (Turgor Pressure)
    real(8) :: psi_sav_pi    ![Pa] Osmotic Potential in seagrass
    real(8) :: theta_sav         !(ratio) Water Content
    real(8) :: osm_sav           !(PSU) Concentration of osmolyte which was converted to PSU

    !!-----E. Acoroides (後で構造体に入れ込む)
    real(8) :: sigma = 0.10d0          !反発係数，どれだけ圧力勾配通りに水が細胞壁を通過するか（一部は通過しない）(仮の値)
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

    real(8) :: Jv_upward               !Pore Water --> Water Column方向の流束密度[m3 m-2 s-1]([m s-1])

    ! For mod_seagrass_flux_NT
    real(8) :: a_HT3
    real(8) :: dpsi_sav_p_dt
    real(8) :: dpsi_sav_pi1_dt       !1がついているものは，オスモライト濃度に由来する水ポテンシャルのみを考慮した値
    real(8) :: tau_psi_sav_pi1       !psi_sav_piが(Jv_leaf*psi_wc+Jv_rhizome*psi_pw)/(Jv_leaf+Jv_rhizome)-psi_sav_p/sigma=0
                                    !となるように指数変化する時の時定数
    real(8) :: dpsi_sav_p
    real(8) :: dpsi_sav_pi1
    real(8) :: psi_sav_pi1
    real(8) :: osm_sav1
    

CONTAINS

    SUBROUTINE initialize_seagrass_flux

        !!!!!-----流束密度関連-----
        !V_sav = 0.015d0 * 0.01d0 * 1.0d0
        !E. Acoroides leef width: 1.0~1.5cm, rhizome diameter: 1.5cm, length = 0.3~1.5m
        V_sav = width_leaf * thickness_leaf * l_leaf + pi * R_rhizome**2 * l_rhizome
        V_sav = V_sav * 1.0d-2 !スケールのテスト  (V_savを小さくするとtheta_savが1に近くなった時の挙動を見れる)
        print *, V_sav
        
        sal_equiv_sav_total = salinity_wc                             !海草内の全水ポテンシャルを海草外のものと一致させて1時間使う(10:00~11:00)

        !!!!!-----流束密度関連-----
            
        !塩分から水ポテンシャルを求める
        psi_wc  = WaterPotentialFromSal(salinity_wc)
        psi_pw  = WaterPotentialFromSal(salinity_pw)
        psi_sav_total = WaterPotentialFromSal(sal_equiv_sav_total)     !上で求めた溶質濃度(PSU)から浸透ポテンシャルを求める
        !print *, sal_equiv_sav_total
        
        ! ---Hofler-Thoday Diagramまたは何らかのプロセスによって、海草の全水ポテンシャルから圧ポテンシャル、浸透ポテンシャル→海草内溶質濃度がわかる。---
        !---もしpsi_sav_pが負なら、負の膨圧（しおれ状態）を示す。thetaは含水率。---
        !CALL HoflerThoday1(psi_sav_total,theta_sav,psi_sav_p,psi_sav_pi,osm_sav)
        !CALL HoflerThoday2(psi_sav_total,theta_sav,psi_sav_p,psi_sav_pi,osm_sav)
        CALL HoflerThoday3(psi_sav_total,theta_sav,psi_sav_p,psi_sav_pi,osm_sav,a_HT3)
        !write(*,fmt='(A)',advance='no') "mod_seagrass_flux.F90 theta_sav=: "
        !print *, theta_sav
        !print *, osm_sav

        dV_sum = 0.0d0

    END SUBROUTINE initialize_seagrass_flux

    SUBROUTINE seagrass_flux

        !!!!!-----流束密度関連-----  
        !塩分から水ポテンシャルを求める
        psi_wc  = WaterPotentialFromSal(salinity_wc)
        psi_pw  = WaterPotentialFromSal(salinity_pw)
        !psi_sav_total = psi_sav_p + psi_sav_pi

        ! 浸透圧差から通水コンダクタンス求める(etaがtimeで変わることも考えてこの位置)
        Lp_leaf = rv_leaf**2/(8.0d0*eta*l_leaf)       !通水コンダクタンスを導出(Hagen-Poseuille流れ)
        Lp_rhizome = rv_rhizome**2/(8.0d0*eta*l_rhizome)    
        !流束密度
        Jv_leaf = Lp_leaf * (sigma * (psi_wc - psi_sav_pi) - psi_sav_p)   ![m3 m-2 s-1]([m s-1])
        Jv_rhizome = Lp_rhizome * (sigma * (psi_pw - psi_sav_pi) - psi_sav_p)

        !水の流入体積
        dV_leaf_dt = Jv_leaf * pi * rv_leaf**2      ![m3 s-1]
        dV_rhizome_dt = Jv_rhizome * pi * rv_leaf**2

        !psi_sav_p（膨圧）の変化
        dpsi_sav_p_dt = 8.0d0 * a_HT3 * (theta_sav - 0.5d0) * 1.0d6 * (dV_leaf_dt + dV_rhizome_dt) / V_sav

        !psi_sav_pi（浸透圧）の変化
        tau_psi_sav_pi1 = 24.0d0 * 3600   ![s] (1day)
        dpsi_sav_pi1_dt = ( (Jv_leaf*psi_wc+Jv_rhizome*psi_pw)/(Jv_leaf+Jv_rhizome) - psi_sav_p / sigma ) / tau_psi_sav_pi1
        !dpsi_sav_pi1_dt = ( (Jv_leaf*psi_wc+Jv_rhizome*psi_pw)/(Jv_leaf+Jv_rhizome) - psi_sav_p ) / tau_psi_sav_pi1
        !print *, dpsi_sav_pi1_dt

    END SUBROUTINE seagrass_flux

    SUBROUTINE update_seagrass_flux(dt_)

        !kob_test.F90からdtを持ってくる
        real(8), intent(in) :: dt_

        !水の流入体積
        dV_leaf = dV_leaf_dt * dt_ * 3600.0d0   ![m3]
        dV_rhizome = dV_rhizome_dt * dt_ * 3600.0d0

        theta_sav = (theta_sav * V_sav + dV_leaf + dV_rhizome) / V_sav
        if (theta_sav < 0.5d0 .or. theta_sav > 1.0d0) stop "mod_seagrass_flux.F90: theta_sav is out of range."

        !psi_sav_p（膨圧）の変化
        dpsi_sav_p = dpsi_sav_p_dt * dt_ * 3600.0d0 
        psi_sav_p = psi_sav_p + dpsi_sav_p
        !psi_sav_p = ( 4.0d0 * a_HT3 * (theta_sav - 0.5d0)**2 ) * 1.0d6

        !psi_sav_pi（浸透圧）の変化
        dpsi_sav_pi1 = dpsi_sav_pi1_dt * dt_ !* 3600.0d0 
        psi_sav_pi1 = psi_sav_pi + dpsi_sav_pi1

        !for test
        dV_sum = dV_sum + dV_leaf + dV_rhizome

        osm_sav1 = SalFromWaterPotential(psi_sav_pi1)
        !print *, osm_sav1

        !上式によってmod_seagrass_flux_HTのosm_savとは少し異なる
        if ( dV_leaf >= 0 .and. dV_rhizome >= 0 ) then
            osm_mass = &
            & theta_sav * V_sav * osm_sav1 + ch_leaf_in * dV_leaf * salinity_wc + ch_rhizome_in * dV_rhizome * salinity_pw
        else if ( dV_leaf >= 0 .and. dV_rhizome < 0 ) then
            osm_mass = &
            & theta_sav * V_sav * osm_sav1 + ch_leaf_in * dV_leaf * salinity_wc + ch_rhizome_out * dV_rhizome * osm_sav
        else if ( dV_leaf < 0 .and. dV_rhizome >= 0 ) then
            osm_mass = &
            & theta_sav * V_sav * osm_sav1 + ch_leaf_out * dV_leaf * osm_sav1 + ch_rhizome_out * dV_rhizome * salinity_pw
        else if ( dV_leaf < 0 .and. dV_rhizome < 0 ) then
            osm_mass = &
            & theta_sav * V_sav * osm_sav1 + ch_leaf_out * dV_leaf * osm_sav1 + ch_rhizome_out * dV_rhizome * osm_sav
        end if

        !溶質濃度（オスモライト）の塩分換算を求めて代入
        osm_sav = osm_mass / (theta_sav * V_sav + dV_leaf + dV_rhizome)
        psi_sav_pi = WaterPotentialFromSal(osm_sav)

        psi_sav_total = psi_sav_p + psi_sav_pi

        sal_equiv_sav_total = SalFromWaterPotential(psi_sav_total)


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


    REAL(8) FUNCTION SalFromWaterPotential(psi_)

        USE mod_param       !temp, R_gas, M_sea
        implicit none

        real(8), intent(in) :: psi_

        real(8) :: rho_
        real(8) :: C
        real(8) :: i_vanthoff

        i_vanthoff = 1.9d0 !Van't Hoff Coefficient
        C = - psi_ / i_vanthoff / R_gas / (temp+273.15d0)
        !逆算して溶質濃度の塩分換算を求める
        SalFromWaterPotential = (-997.0d0+sqrt(997.0d0**2.0d0 - 4.0d0*(-M_sea*C)))/2.0d0

    END FUNCTION SalFromWaterPotential



    SUBROUTINE HoflerThoday3(psi_total_,theta_,psi_p_,psi_pi_,osm_,a)
        !Extended version of HoflerThoday2 by changing parameters
        
        !If a:b = 1:4 --> HoflerThoday2 version
        !psi_pi = -4(b-a)(theta-1)**2 - a [MPa] (theta, psi) = (0.5, -b), (1, -a) axis theta = 1
        !psi_p = 4a(theta-0.5)**2       [MPa] (theta, psi) = (0.5, 0), (1, a)     axis theta = 0.5
        !psi_total = psi_pi + psi_p
        !psi_total = -4(b-2a)(theta- (2b-3a)/2(b-2a))**2 + a**2/(b-2a)[MPa]
        
        USE mod_param       !R_gas, M_sea
        implicit none

        real(8), intent(in)  :: psi_total_        ![Pa] 海草内の全水ポテンシャル
        real(8), intent(out) :: theta_            ![ratio] 含水率
        real(8), intent(out) :: psi_p_            ![Pa] 圧ポテンシャル（膨圧）
        real(8), intent(out) :: psi_pi_           ![Pa] 浸透ポテンシャル
        real(8) :: C
        real(8) :: i_vanthoff
        real(8), intent(out) :: osm_              !(PSU) 逆算して溶質濃度の塩分換算を求めた擬似的な塩分
        real(8), intent(inout) :: a
        real(8) :: b 

        a = 1.0d0
        b = 6.0d0

        !Hofler-Thoday Diagramまたは何らかのプロセスによって、
        !海草の全水ポテンシャルから，含水率，圧ポテンシャル、浸透ポテンシャル，海草内溶質濃度の塩分換算がわかる。
        theta_ = (2*b-3*a)/(2*(b-2*a)) - sqrt(  ( (a**2/(b-2*a) - psi_total_ /1.0d6) / (4*(b-2*a)) ) )

        if (theta_ < 0.5) then
            stop "mod_seagrass_flux.F90: Water content is under 50%. (undefined region)"
        end if

        psi_p_  = (4*a * (theta_ - 0.5d0)**2) * 1.0d6
        psi_pi_ = (-4*(b-a) * (theta_ - 1.0d0)**2 - a) * 1.0d6
        
        i_vanthoff = 1.9d0 !Van't Hoff Coefficient
        C = -psi_pi_ / i_vanthoff / R_gas / (temp+273.15d0)
        !逆算して溶質濃度の塩分換算を求める
        osm_ = (-997.0d0+sqrt(997.0d0**2.0d0 - 4.0d0*(-M_sea*C)))/2.0d0

    END SUBROUTINE HoflerThoday3

    


END MODULE mod_seagrass_flux