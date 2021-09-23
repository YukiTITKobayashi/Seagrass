!#include "cppdefs.h"
    MODULE mod_calc

        implicit none

    CONTAINS

        SUBROUTINE gauss_jordan(a0, x, b, n)

            !Gauss-Jordan法によるn個の連立一次方程式Ax=bの求解
            !牛島省『数値計算のためのFortran90/95プログラミング入門』森北出版（編集後消す！！！）

            integer, intent(in) :: n        !連立方程式の数(rank)
            real(8), intent(in) :: a0(n,n)   !Aは正方行列
            real(8), intent(in) :: b(n)     !bはn次元のベクトル
            real(8), intent(out) :: x(n)    !Ax=bの解

            integer i,k

            real(8) ar, a(n,n)

            a(:,:) = a0(:,:)
            x(:) = b(:)

            do k = 1,n

                if ( a(k,k) == 0.0d0 ) stop 'pivot = 0'
                ar = 1.0d0 / a(k,k)             !arは対角成分の逆数
                a(k,k) = 1.0d0                  !対角成分に1を設定
                a(k,k+1:n) = a(k,k+1:n) * ar    !k行のk+1~n列に*ar
                x(k) = x(k) * ar                !k行の右辺にも*ar

                do i = 1,n
                    if ( i /= k ) then
                        a(i,k+1:n) = a(i,k+1:n) - a(i,k) * a(k,k+1:n)
                        x(i) = x(i) - a(i,k) * x(k)
                        a(i,k) = 0.0d0
                    end if
                end do

            end do

        END SUBROUTINE gauss_jordan

        SUBROUTINE gaussian(a,x,b,n)
            !ガウスの消去法Gaussian elimination
        END SUBROUTINE gaussian

        SUBROUTINE set_random_ab(a,x,b,n)

            !nを取得、a,b,xを割付、aとbに乱数を設定

            integer, intent(out) :: n
            real(8), allocatable, intent(out) :: a(:,:)
            real(8), allocatable, intent(out) :: x(:)
            real(8), allocatable, intent(out) :: b(:)

            write(*,'(A)',advance='no') 'input n: '
            read(*,*) n
            if ( n < 1 ) stop 'n must be positive integer.'

            allocate(a(n,n))
            allocate(x(n))
            allocate(b(n))

            call random_number(a)
            call random_number(b)

        END SUBROUTINE set_random_ab
        
    END MODULE mod_calc

    PROGRAM main

        USE mod_calc
        implicit none

        real(8), allocatable :: a(:,:)
        real(8), allocatable :: b(:)
        real(8), allocatable :: x(:)
        real(8), allocatable :: rest(:)     !残差ベクトル rest = b - Ax
        integer n

        CALL set_random_ab(a,x,b,n)
        CALL gauss_jordan(a,x,b,n)
        allocate(rest(n))
        rest(:) = b(:) - matmul(a,x)
        write(*,*) 'Gauss-Jordan error = ', dot_product(rest,rest)
        deallocate( a,b,x )

    END PROGRAM main