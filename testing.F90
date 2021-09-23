PROGRAM testing
    implicit none

    real(8) time
    real(8) start_time

    integer hour
        !!Region
    integer min
        !!MiddleRegion
    integer sec
        !!EndRegion

    character(len=10) :: disp

    time = 2./24. + 30./60./24.
    start_time = 11./24.

    disp = time_display(time,start_time,hour,min,sec)

    print *, disp
    !!-!!--! Start (setting.json -> explicitFolding.rules)
    print *, hour
    !!-!-!-! Middle (setting.json -> explicitFolding.rules)
    print *, min
    !!-!--!! End (setting.json -> explicitFolding.rules)
    print *, sec

contains
    CHARACTER(len=10) FUNCTION time_display(time,start_time,hour,minute,second)

        implicit none

        real(8), intent(in) :: time
        real(8), intent(in) :: start_time

        real(8) now
        integer day
        integer,intent(out) :: hour
        integer,intent(out) :: minute
        integer,intent(out) :: second

        now = start_time + time
        day    = int(now)
        hour   = int((now - day) * 24.e0)
        minute = int((now - day - hour/24.e0) * 24.e0 * 60.e0)
        second = int((now - day - hour/24.e0 - minute/24.e0/60.e0) * 24.e0 * 60.e0 * 60.e0)

        write(time_display(1:1),'(I1.1)') day
        !!Region Level 1
        write(time_display(2:2),'(A)') "D"
        !!Region Level 2
        write(time_display(3:4),'(I2.2)') hour
        !!EndRegion Level 2
        write(time_display(5:5),'(A)') ":"
        !!MiddleRegion Level 1
        write(time_display(6:7),'(I2.2)') minute
        !!EndRegion Level 1
        write(time_display(8:8),'(A)') ":"
        write(time_display(9:10),'(I2.2)') second

    END FUNCTION
END PROGRAM testing