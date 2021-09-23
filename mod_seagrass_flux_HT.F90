!--------------------------------------------------------------------------------
!  mod_seagrass_flux_HT.F90
!  Module about Seagrass flux
!  The version using Hofler-Thoday diagram, however, gave up making due to mathematical difficulty.
!--------------------------------------------------------------------------------

#include "cppdefs.h"

    MODULE mod_seagrass_flux
        
        USE mod_param
        USE mod_input
        !USE mod_seagrass

        implicit none

        ! -----流束密度関連-----
        real(8) :: V_sav         ![m3] Volume of seagrass(width*thickness*height*shoot)

        real(8) :: psi_wc        ![Pa] Water Potential of Water Column
        real(8) :: psi_pw        ![Pa] Water Potential of Pore Water
        real(8) :: psi_sav_total       ![Pa] Water Potential in seagrass
        real(8) :: psi_sav_p     ![Pa] Pressure Component of seagrass (Turgor Pressure)
        real(8) :: psi_sav_pi    ![Pa] Osmotic Potential in seagrass
        real(8) :: theta_sav         !(ratio) Water Content
        real(8) :: osm_sav           !(PSU) Concentration of osmolyte and ions inside the seeagrass
        real(8) :: sal_equiv_sav_total    !(PSU)Imaginary salinity for seagrass 
                                          !which is equivalent to the total water potential and converted to PSU

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

        real(8) :: Jv_upward               !Pore Water --> Water Column方向の流束密度[m3 m-2 s-1]([m s-1])

    CONTAINS

        SUBROUTINE initialize_seagrass_flux

            !!!!!-----流束密度関連-----
            !V_sav = 0.015d0 * 0.01d0 * 1.0d0
            !E. Acoroides leef width: 1.0~1.5cm, rhizome diameter: 1.5cm, length = 0.3~1.5m
            V_sav = width_leaf * thickness_leaf * l_leaf + pi * R_rhizome**2 * l_rhizome
            V_sav = V_sav * 1.0d-2 !スケールのテスト
            print *, V_sav
            !水の流入体積
            !dV_leaf = 0.0d0
            !dV_rhizome = 0.0d0
            ! ---塩分から水流束密度を求める---
            !psi_wc  = WaterPotentialFromSal(salinity_wc)
            !psi_pw  = WaterPotentialFromSal(salinity_pw)
            sal_equiv_sav_total = salinity_wc                  !海草内の全水ポテンシャルを海草外のものと一致させて1時間使う(10:00~11:00)
            !osm_sav = salinity_wc                             !海草内の全水ポテンシャルを海草外のものと一致させて1時間使う(10:00~11:00)

            ! ---Hofler-Thoday Diagramまたは何らかのプロセスによって、海草の全水ポテンシャルから圧ポテンシャル、浸透ポテンシャル→海草内溶質濃度がわかる。---
            !---もしpsi_sav_pが負なら、負の膨圧（しおれ状態）を示す。thetaは含水率。---
            !CALL HoflerThoday1(psi_sav_total,theta_sav,psi_sav_p,psi_sav_pi,osm_sav)
            !write(*,fmt='(A)',advance='no') "mod_seagrass_flux.F90 theta_sav=: "
            !print *, theta_sav

            ! 浸透圧差から通水コンダクタンス求める
            !Lp_leaf = rv_leaf**2/(8.0d0*eta*l)       !通水コンダクタンスを導出(Hagen-Poseuille流れ)    
            !流束密度
            !Jv_leaf = Lp * (sigma * (psi_wc - psi_sav_pi) - psi_sav_p)   ![m3 m-2 s-1]([m s-1])
            !Jv_rhizome = Lp * (sigma * (psi_pw - psi_sav_pi) - psi_sav_p)

            !For print
            dV_sum = 0.0d0

        END SUBROUTINE initialize_seagrass_flux

        SUBROUTINE seagrass_flux

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
            CALL HoflerThoday3(psi_sav_total,theta_sav,psi_sav_p,psi_sav_pi,osm_sav)
            !write(*,fmt='(A)',advance='no') "mod_seagrass_flux.F90 theta_sav=: "
            !print *, theta_sav
            !print *, osm_sav

            ! 浸透圧差から通水コンダクタンス求める(etaがtimeで変わることも考えてこの位置)
            Lp_leaf = rv_leaf**2/(8.0d0*eta*l_leaf)       !通水コンダクタンスを導出(Hagen-Poseuille流れ)
            Lp_rhizome = rv_rhizome**2/(8.0d0*eta*l_rhizome)    
            !流束密度
            Jv_leaf = Lp_leaf * (sigma * (psi_wc - psi_sav_pi) - psi_sav_p)   ![m3 m-2 s-1]([m s-1])
            Jv_rhizome = Lp_rhizome * (sigma * (psi_pw - psi_sav_pi) - psi_sav_p)

            !水の流入体積
            dV_leaf_dt = Jv_leaf * pi * rv_leaf**2      ![m3 s-1]
            dV_rhizome_dt = Jv_rhizome * pi * rv_leaf**2

        END SUBROUTINE seagrass_flux

        SUBROUTINE update_seagrass_flux(dt_)

            !kob_test.F90からdtを持ってくる
            real(8), intent(in) :: dt_

            !水の流入体積
            dV_leaf = dV_leaf_dt * dt_ * 3600.0d0   ![m3]
            dV_rhizome = dV_rhizome_dt * dt_ * 3600.0d0

            !for test
            dV_sum = dV_sum + dV_leaf + dV_rhizome

            !V_sav = V_sav + dV_leaf + dV_rhizome
            theta_sav = (theta_sav * V_sav + dV_leaf + dV_rhizome) / V_sav

            !osm_mass = theta_sav * V_sav * osm_sav + dV_leaf * salinity_wc + dV_rhizome * salinity_pw

            !水の流出入によって海草内の
            if ( dV_leaf >= 0 .and. dV_rhizome >= 0 ) then
                osm_mass = &
                & theta_sav * V_sav * osm_sav + ch_leaf_in * dV_leaf * salinity_wc + ch_rhizome_in * dV_rhizome * salinity_pw
            else if ( dV_leaf >= 0 .and. dV_rhizome < 0 ) then
                osm_mass = theta_sav * V_sav * osm_sav + ch_leaf_in * dV_leaf * salinity_wc + ch_rhizome_out * dV_rhizome * osm_sav
            else if ( dV_leaf < 0 .and. dV_rhizome >= 0 ) then
                osm_mass = theta_sav * V_sav * osm_sav + ch_leaf_out * dV_leaf * osm_sav + ch_rhizome_out * dV_rhizome * salinity_pw
            else if ( dV_leaf < 0 .and. dV_rhizome < 0 ) then
                osm_mass = theta_sav * V_sav * osm_sav + ch_leaf_out * dV_leaf * osm_sav + ch_rhizome_out * dV_rhizome * osm_sav
            end if

            !溶質濃度（オスモライト）の塩分換算を求めて代入
            !osm_sav = ((V_sav * osm_sav) + (dV_leaf * salinity_wc) + (dV_rhizome * salinity_pw)) / (V_sav + dV_leaf + dV_rhizome)
            sal_equiv_sav_total = osm_mass / (theta_sav * V_sav + dV_leaf + dV_rhizome)


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

        SUBROUTINE HoflerThoday1(psi_total_,theta_,psi_p_,psi_pi_,osm_)
            !version of "Plants and Microclimate" p.89
            
            USE mod_param       !R_gas, M_sea
            implicit none

            real(8), intent(in)  :: psi_total_        ![Pa] 海草内の全水ポテンシャル
            real(8), intent(out) :: theta_            ![ratio] 含水率
            real(8), intent(out) :: psi_p_            ![Pa] 圧ポテンシャル（膨圧）
            real(8), intent(out) :: psi_pi_           ![Pa] 浸透ポテンシャル
            real(8) :: C
            real(8) :: i_vanthoff
            real(8), intent(out) :: osm_              !(PSU) 逆算して溶質濃度の塩分換算を求めた擬似的な塩分

            !Hofler-Thoday Diagramまたは何らかのプロセスによって、
            !海草の全水ポテンシャルから，含水率，圧ポテンシャル、浸透ポテンシャル，海草内溶質濃度の塩分換算がわかる。
            !以下は『植物と微気象』P89の図の場合
            !psi_pi = 2.5theta-4.5 [MPa]
            !psi_p = 8theta-6 (theta>=0.75), psi_p = 0 (theta<0.75)
            !psi_total = psi_pi + psi_p
            !しおれ点は(theta,psi)=(0.75,-2.625MPa)

            IF (psi_total_>=-2.625d6) THEN
                theta_ = (psi_total_/1.0d6 + 10.5d0) / 10.5d0      !psi_total=0でtheta=1
            ELSE
                theta_ = (psi_total_/1.0d6 + 4.5d0) / 2.5d0        !膨圧0より下の場合
            ENDIF

            IF (theta_>=0.75) THEN
                psi_p_  = (8.0d0 * theta_ - 6.0d0) * 1.0d6
            else
                psi_p_  = 0.0d0
            ENDIF

            psi_pi_ = (2.5d0 * theta_ - 4.5d0) * 1.0d6
            
            i_vanthoff = 1.9d0 !Van't Hoff Coefficient
            C = -psi_pi_ / i_vanthoff / R_gas / (temp+273.15d0)
            !逆算して溶質濃度の塩分換算を求める
            osm_ = (-997.0d0+sqrt(997.0d0**2.0d0 - 4.0d0*(-M_sea*C)))/2.0d0

        END SUBROUTINE HoflerThoday1
 

        SUBROUTINE HoflerThoday2(psi_total_,theta_,psi_p_,psi_pi_,osm_)
            !modified version of Touchette(2006) by changing parameters
            !If Amp = 1.0d0 --> Touchette's version
            !以下はTouchette(2006)の場合
            !psi_pi = -3.6(theta-1)**2 - 0.3 [MPa] (theta, psi) = (0.5, -1.2), (1, -0.3) axis theta = 1
            !psi_p = 1.2(theta-0.5)**2       [MPa] (theta, psi) = (0.5, 0), (1, 0.3)     axis theta = 0.5
            !psi_total = psi_pi + psi_p
            !psi_total = -2.4(theta-1.25)**2 + 0.15 = -1.2(2theta**2 -5theta + 3) = -2.4(theta-1)(theta-1.5)[MPa]

            !If Amp = 4.0d0 --> Modified version
            !以下はTouchette(2006)のHofler-Thoday diagramをpsi方向にAmp=4倍したもの
            !psi_pi = -14.4(theta-1)**2 - 1.2 [MPa] (theta, psi) = (0.5, -4.8), (1, -1.2) axis theta = 1
            !psi_p = 4.8(theta-0.5)**2       [MPa] (theta, psi) = (0.5, 0), (1, 1.2)     axis theta = 0.5
            !psi_total = psi_pi + psi_p
            !psi_total = -9.6(theta-1.25)**2 + 0.6 = -4.8(2theta**2 -5theta + 3) = -9.6(theta-1)(theta-1.5)[MPa]
            
            USE mod_param       !R_gas, M_sea
            implicit none

            real(8), intent(in)  :: psi_total_        ![Pa] 海草内の全水ポテンシャル
            real(8), intent(out) :: theta_            ![ratio] 含水率
            real(8), intent(out) :: psi_p_            ![Pa] 圧ポテンシャル（膨圧）
            real(8), intent(out) :: psi_pi_           ![Pa] 浸透ポテンシャル
            real(8) :: C
            real(8) :: i_vanthoff
            real(8), intent(out) :: osm_              !(PSU) 逆算して溶質濃度の塩分換算を求めた擬似的な塩分
            real(8) :: Amp = 4.0d0

            !Hofler-Thoday Diagramまたは何らかのプロセスによって、
            !海草の全水ポテンシャルから，含水率，圧ポテンシャル、浸透ポテンシャル，海草内溶質濃度の塩分換算がわかる。
            theta_ = 1.25d0 - sqrt(  (0.6d0 - psi_total_ /1.0d6) / (Amp*2.4d0) )

            if (theta_ < 0.5) then
                stop "mod_seagrass_flux.F90: Water content is under 50%. (undefind region)"
            end if

            psi_p_  = ((Amp*1.2d0) * (theta_ - 0.5d0)**2) * 1.0d6
            psi_pi_ = (-(Amp*3.6d0) * (theta_ - 1.0d0)**2 - (Amp*0.3d0)) * 1.0d6
            
            i_vanthoff = 1.9d0 !Van't Hoff Coefficient
            C = -psi_pi_ / i_vanthoff / R_gas / (temp+273.15d0)
            !逆算して溶質濃度の塩分換算を求める
            osm_ = (-997.0d0+sqrt(997.0d0**2.0d0 - 4.0d0*(-M_sea*C)))/2.0d0

        END SUBROUTINE HoflerThoday2


        SUBROUTINE HoflerThoday3(psi_total_,theta_,psi_p_,psi_pi_,osm_)
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
            real(8) :: a = 1.0d0
            real(8) :: b = 6.0d0

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


!!!!!old functions
!!!!!!!!Vの更新をどこでするか迷った以前のバージョン
        SUBROUTINE initialize_seagrass_flux_old

            !!!!!-----流束密度関連-----
            !V_sav = 0.015d0 * 0.01d0 * 1.0d0
            !E. Acoroides leef width: 1.0~1.5cm, rhizome diameter: 1.5cm, length = 0.3~1.5m
            V_sav = 1.0d-6
            !水の流入体積
            dV_leaf = 0.0d0
            dV_rhizome = 0.0d0
            ! ---塩分から水流束密度を求める---
            !psi_wc  = WaterPotentialFromSal(salinity_wc)
            !psi_pw  = WaterPotentialFromSal(salinity_pw)
            !psi_sav_total = psi_wc                            !海草内の全水ポテンシャルを海草外のものと一致させて1時間使う(10:00~11:00)
            osm_sav  = salinity_wc                             !海草内の全水ポテンシャルを海草外のものと一致させて1時間使う(10:00~11:00)

            ! ---Hofler-Thoday Diagramまたは何らかのプロセスによって、海草の全水ポテンシャルから圧ポテンシャル、浸透ポテンシャル→海草内溶質濃度がわかる。---
            !---もしpsi_sav_pが負なら、負の膨圧（しおれ状態）を示す。thetaは含水率。---
            !CALL HoflerThoday1(psi_sav_total,theta_sav,psi_sav_p,psi_sav_pi,osm_sav)
            !write(*,fmt='(A)',advance='no') "mod_seagrass_flux.F90 theta_sav=: "
            !print *, theta_sav

            ! 浸透圧差から通水コンダクタンス求める
            !Lp = r_xylem**2/(8.0d0*eta*l)       !通水コンダクタンスを導出(Hagen-Poseuille流れ)    
            !流束密度
            !Jv_leaf = Lp * (sigma * (psi_wc - psi_sav_pi) - psi_sav_p)   ![m3 m-2 s-1]([m s-1])
            !Jv_rhizome = Lp * (sigma * (psi_pw - psi_sav_pi) - psi_sav_p)

        END SUBROUTINE initialize_seagrass_flux_old

        SUBROUTINE seagrass_flux_old(dt_)

            !kob_test.F90からdtを持ってくる
            real(8), intent(in) :: dt_

            !!!!!-----流束密度関連-----

            !水の流入体積
            dV_leaf = Jv_leaf * pi * rv_leaf**2 * dt_ * 3600.0d0   ![m3]
            dV_rhizome = Jv_rhizome * pi * rv_rhizome**2 * dt_ * 3600.0d0

            osm_mass = (V_sav * osm_sav) + (dV_leaf * salinity_wc) + (dV_rhizome * salinity_pw)
            V_sav = V_sav + dV_leaf + dV_rhizome
            
            
            !溶質濃度（オスモライト）の塩分換算を求めて代入
            !osm_sav = ((V_sav * osm_sav) + (dV_leaf * salinity_wc) + (dV_rhizome * salinity_pw)) / (V_sav + dV_leaf + dV_rhizome)
            osm_sav = osm_mass / V_sav

            !塩分から水ポテンシャルを求める
            psi_wc  = WaterPotentialFromSal(salinity_wc)
            psi_pw  = WaterPotentialFromSal(salinity_pw)
            psi_sav_total = WaterPotentialFromSal(osm_sav)     !上で求めた溶質濃度(PSU)から浸透ポテンシャルを求める
            
            ! ---Hofler-Thoday Diagramまたは何らかのプロセスによって、海草の全水ポテンシャルから圧ポテンシャル、浸透ポテンシャル→海草内溶質濃度がわかる。---
            !---もしpsi_sav_pが負なら、負の膨圧（しおれ状態）を示す。thetaは含水率。---
            CALL HoflerThoday1(psi_sav_total,theta_sav,psi_sav_p,psi_sav_pi,osm_sav)
            !CALL HoflerThoday2(psi_sav_total,theta_sav,psi_sav_p,psi_sav_pi,osm_sav)
            !write(*,fmt='(A)',advance='no') "mod_seagrass_flux.F90 theta_sav=: "
            !print *, theta_sav

            ! 浸透圧差から通水コンダクタンス求める
            Lp_leaf = rv_leaf**2/(8.0d0*eta*l_leaf)       !通水コンダクタンスを導出(Hagen-Poseuille流れ)   
            Lp_rhizome = rv_rhizome**2/(8.0d0*eta*l_rhizome)  
            !流束密度
            Jv_leaf = Lp_leaf * (sigma * (psi_wc - psi_sav_pi) - psi_sav_p)   ![m3 m-2 s-1]([m s-1])
            Jv_rhizome = Lp_rhizome * (sigma * (psi_pw - psi_sav_pi) - psi_sav_p)

        END SUBROUTINE seagrass_flux_old

!!!!!!!!容積の増加を考慮していた以前のバージョン

        SUBROUTINE update_seagrass_flux_old(dt_)

            !kob_test.F90からdtを持ってくる
            real(8), intent(in) :: dt_

            !水の流入体積
            dV_leaf = dV_leaf_dt * dt_ * 3600.0d0   ![m3]
            dV_rhizome = dV_rhizome_dt * dt_ * 3600.0d0

            osm_mass = (V_sav * osm_sav) + (dV_leaf * salinity_wc) + (dV_rhizome * salinity_pw)
            V_sav = V_sav + dV_leaf + dV_rhizome
  
            !溶質濃度（オスモライト）の塩分換算を求めて代入
            !osm_sav = ((V_sav * osm_sav) + (dV_leaf * salinity_wc) + (dV_rhizome * salinity_pw)) / (V_sav + dV_leaf + dV_rhizome)
            osm_sav = osm_mass / V_sav

            !for test
            !dV_sum = dV_sum + dV_leaf + dV_rhizome

        END SUBROUTINE update_seagrass_flux_old

        SUBROUTINE HoflerThoday2_old(psi_total_,theta_,psi_p_,psi_pi_,osm_)
            !version of Touchette(2006)
            !Cannot use because the initial water potential given from psi_wc is -2.6MPa, which is out of range.
            
            USE mod_param       !R_gas, M_sea
            implicit none

            real(8), intent(in)  :: psi_total_        ![Pa] 海草内の全水ポテンシャル
            real(8), intent(out) :: theta_            ![ratio] 含水率
            real(8), intent(out) :: psi_p_            ![Pa] 圧ポテンシャル（膨圧）
            real(8), intent(out) :: psi_pi_           ![Pa] 浸透ポテンシャル
            real(8) :: C
            real(8) :: i_vanthoff
            real(8), intent(out) :: osm_              !(PSU) 逆算して溶質濃度の塩分換算を求めた擬似的な塩分

            !Hofler-Thoday Diagramまたは何らかのプロセスによって、
            !海草の全水ポテンシャルから，含水率，圧ポテンシャル、浸透ポテンシャル，海草内溶質濃度の塩分換算がわかる。
            !以下はTouchette(2006)の場合
            !psi_pi = -3.6(theta-1)**2 - 0.3 [MPa] (theta, psi) = (0.5, -1.2), (1, -0.3) axis theta = 1
            !psi_p = 1.2(theta-0.5)**2       [MPa] (theta, psi) = (0.5, 0), (1, 0.3)     axis theta = 0.5
            !psi_total = psi_pi + psi_p
            !psi_total = -2.4(theta-1.25)**2 + 0.15 = -1.2(2theta**2 -5theta + 3) = -2.4(theta-1)(theta-1.5)[MPa]
            theta_ = 1.25d0 - sqrt(  (0.15d0 - psi_total_ /1.0d6) / 2.4d0 )

            !if (theta_ < 0.5) then
            !    stop "mod_seagrass_flux.F90: Water content is under 50%. (undefind region)"
            !end if

            psi_p_  = (1.2d0 * (theta_-0.5d0)**2) * 1.0d6
            psi_pi_ = (-3.6d0 * (theta_-1.0d0)**2 - 0.3d0) * 1.0d6
            
            i_vanthoff = 1.9d0 !Van't Hoff Coefficient
            C = -psi_pi_ / i_vanthoff / R_gas / (temp+273.15d0)
            !逆算して溶質濃度の塩分換算を求める
            osm_ = (-997.0d0+sqrt(997.0d0**2.0d0 - 4.0d0*(-M_sea*C)))/2.0d0

        END SUBROUTINE HoflerThoday2_old


    END MODULE mod_seagrass_flux