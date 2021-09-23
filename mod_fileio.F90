MODULE mod_fileio

	!Trusted 2020/08/29

	implicit none
CONTAINS

	SUBROUTINE csv_read(fid,filename,real_alloc_2d_array,buf,nm_row,nm_column)

		!サブルーチンの使用後、メインプログラムでのdeallocate(real_alloc_2d_array)を忘れない！
		!filenameをパスとする行数・列数がわからないアスキーファイルを読み込み、二次元配列real_alloc_2d_arrayに格納する
		!Fortranはカンマを無視して数値を読み込むので、.txt/.datと.csvでコードを変える必要はない。
		!SUBROUTINE buffer(buf)を使用。buf行のヘッダーを無視

		!可変長の文字列のファイル間受け渡しはcharacter(:), allocatableではなくcharacter(len=*)とする。
		!追記：仮引数のみで可能？メインプログラムのcsv_read.F90でcharacter(len=*)とするとエラーになるので、逆にcharacter(:), allocatableが正しい(?)
		integer, intent(in) :: fid
		character(len=*), intent(in) :: filename
		real(8), allocatable, intent(out) :: real_alloc_2d_array(:,:)
		integer, intent(in) :: buf
		integer, intent(out) :: nm_row
		integer, intent(out) :: nm_column

		!読み込んだ配列を表示するか否か
		logical :: isDisplayed = .false.

		integer :: i
		integer :: j
		integer :: ios

		character(1000) :: line

		!integer :: len_filename
		!len_filename = len_trim(filename)

		open(fid,file=filename, status='old')

    	CALL buffer(fid,buf)
    
		!Count row number
		nm_row = 0
		do
			read(fid,*,iostat=ios)
			if (ios==-1) exit
			nm_row = nm_row + 1
		end do
		rewind(fid)


		!Count column number (=the number of ',' + 1 or ' ' + 1) (クオーテーション未配慮、.TXT,.DAT未配慮)
		nm_column = 1
		CALL buffer(fid,buf)
		read(fid,'(A)') line

		if (len_trim(line) /= 0) then

			if ( index(filename,".csv") /= 0 .or. index(filename,".CSV") /= 0) then	!FORTRANでも文字列は流石に大文字小文字区別あり

				do j = 1,len_trim(line)
					if (line(j:j) == ',') then      !文字列の一文字指定はcharacter(i)ではなくcharacter(i:i)
						nm_column = nm_column + 1
					end if
				end do

			else if ( index(filename,".txt") /= 0 .or. index(filename,".dat") /= 0 ) then

				do j = 1,len_trim(line)
					if (line(j:j) == ' ') then      !文字列の一文字指定はcharacter(i)ではなくcharacter(i:i)
						nm_column = nm_column + 1
					end if
				end do

			else
				stop "User error statement: The extension must be .csv, .txt, or .dat. Check the data file or modify the code."
			end if

		else
			stop "User error statement: The first line is empty."     !どちらも','0コの1列と0列を区別
		end if

		rewind(fid)

		!allocate(character(len_trim(line)) :: line) 
		!characterの動的割付がうまく行かない。allocateしていない配列に代入は不可?
		!結局文字列長を求める際に、一行目の文字数よりより大きな要素数の文字型変数を定義する必要があるので、allocateさせる意味がない。


		allocate( real_alloc_2d_array(nm_row,nm_column) )
		!write(*,*) "Row:", nm_row, ", Column:", nm_column

		!Read data
		CALL buffer(fid,buf)
		do i = 1, nm_row
			read(fid,*) ( real_alloc_2d_array(i,j), j=1,nm_column )
		end do
		

		close(fid)

		if (isDisplayed) then
			write(*,*) "The array could be successfully loaded."
			do i = 1,nm_row
				write(*,'(999f15.8)') ( real_alloc_2d_array(i,j), j=1,nm_column )	!nm_column > 999だとエラーになる
			end do
		end if
	END SUBROUTINE csv_read

	SUBROUTINE csv_read_char(fid,filename,char_alloc_2d_array,buf,nm_row,nm_column)

		!サブルーチンの使用後、メインプログラムでのdeallocate(real_alloc_2d_array)を忘れない！
		!filenameをパスとする行数・列数がわからないアスキーファイルを読み込み、二次元配列real_alloc_2d_arrayに格納する
		!Fortranはカンマを無視して数値を読み込むので、.txt/.datと.csvでコードを変える必要はない。
		!SUBROUTINE buffer(buf)を使用。buf行のヘッダーを無視

		!可変長の文字列のファイル間受け渡しはcharacter(:), allocatableではなくcharacter(len=*)とする。
		!追記：仮引数のみで可能？メインプログラムのcsv_read.F90でcharacter(len=*)とするとエラーになるので、逆にcharacter(:), allocatableが正しい(?)
		integer, intent(in) :: fid
		character(len=*), intent(in) :: filename
		character(len=50), allocatable, intent(out) :: char_alloc_2d_array(:,:)
		integer, intent(in) :: buf
		integer, intent(out) :: nm_row
		integer, intent(out) :: nm_column

		!読み込んだ配列を表示するか否か
		logical :: isDisplayed = .false.

		integer :: i
		integer :: j
		integer :: ios

		character(1000) :: line

		integer :: col
		integer :: k

		!integer :: len_filename
		!len_filename = len_trim(filename)

		open(fid,file=filename, status='old')

    	CALL buffer(fid,buf)
    
		!Count row number
		nm_row = 0
		do
			read(fid,*,iostat=ios)
			if (ios==-1) exit
			nm_row = nm_row + 1
		end do
		rewind(fid)


		!Count column number (=the number of ',' + 1 or ' ' + 1) (クオーテーション未配慮、.TXT,.DAT未配慮)
		nm_column = 1
		CALL buffer(fid,buf)
		read(fid,'(A)') line

		if (len_trim(line) /= 0) then

			if ( index(filename,".csv") /= 0 .or. index(filename,".CSV") /= 0) then	!FORTRANでも文字列は流石に大文字小文字区別あり

				do j = 1,len_trim(line)
					if (line(j:j) == ',') then      !文字列の一文字指定はcharacter(i)ではなくcharacter(i:i)
						nm_column = nm_column + 1
					end if
				end do

			else if ( index(filename,".txt") /= 0 .or. index(filename,".dat") /= 0 ) then

				do j = 1,len_trim(line)
					if (line(j:j) == ' ') then      !文字列の一文字指定はcharacter(i)ではなくcharacter(i:i)
						nm_column = nm_column + 1
					end if
				end do

			else
				stop "User error statement: The extension must be .csv, .txt, or .dat. Check the data file or modify the code."
			end if

		else
			stop "User error statement: The first line is empty."     !どちらも','0コの1列と0列を区別
		end if

		rewind(fid)

		!allocate(character(len_trim(line)) :: line) 
		!characterの動的割付がうまく行かない。allocateしていない配列に代入は不可?
		!結局文字列長を求める際に、一行目の文字数よりより大きな要素数の文字型変数を定義する必要があるので、allocateさせる意味がない。


		allocate( char_alloc_2d_array(nm_row,nm_column) )
		!write(*,*) "Row:", nm_row, ", Column:", nm_column

		!Read data
		CALL buffer(fid,buf)

		do i = 1, nm_row

			read(fid,'(A)') line
			col = 1
			k = 0

			if (len_trim(line) /= 0) then

				if ( index(filename,".csv") /= 0 .or. index(filename,".CSV") /= 0) then	!FORTRANでも文字列は流石に大文字小文字区別あり
	
					do j = 1,len_trim(line)
						if (line(j:j) == ',') then      !文字列の一文字指定はcharacter(i)ではなくcharacter(i:i)
							char_alloc_2d_array(i,col) = line(k+1:j-1)
							col = col + 1
							k = j
						end if

						if (col == nm_column) exit	
					end do
					char_alloc_2d_array(i,col) = trim(line(k+1:len_trim(line)))
	
				else if ( index(filename,".txt") /= 0 .or. index(filename,".dat") /= 0 ) then
	
					do j = 1,len_trim(line)
						if (line(j:j) == ' ') then      !文字列の一文字指定はcharacter(i)ではなくcharacter(i:i)
							char_alloc_2d_array(i,col) = line(k+1:j-1)
							col = col + 1
							k = j
						end if

						if (col == nm_column) exit	
					end do
					char_alloc_2d_array(i,col) = trim(line(k+1:len_trim(line)))

				end if
	
			else
				stop "User error statement: The line is empty."     !どちらも','0コの1列と0列を区別
			end if

		end do


		!do i = 1, nm_row
		!	read(fid,'(A)') ( char_alloc_2d_array(i,j), j=1,nm_column )
		!end do
		

		close(fid)

		if (isDisplayed) then
			write(*,*) "The array could be successfully loaded."
			do i = 1,nm_row
				write(*,'(A)') ( char_alloc_2d_array(i,j), j=1,nm_column )	!nm_column > 999だとエラーになる
			end do
		end if
	END SUBROUTINE csv_read_char

	
	SUBROUTINE csv_write(fid,filename,real_alloc_2d_array,fmt_type,cell_length,nm_row,nm_column)

		!任意の行列数(nm_row,nm_column)を持つ二次元配列real_alloc_2d_arrayの実数値を
		!fmt_type、cell_lengthの書式でfilenameのCSVに書き込む。

		!サブルーチンの使用後にdeallocateする

		integer, intent(in) :: fid 
		character(len=*), intent(in) :: filename
		real(8), allocatable, intent(in) :: real_alloc_2d_array(:,:)

		!CFMT='(fmt_type cell_length.decimal_digit)'
		character(len=1), intent(in) :: fmt_type	!e or f
		integer, intent(in) :: cell_length			!15 for 'f' -> integer_digit = 6
		integer :: decimal_digit = 8				!8 for 'f' (real(8))
		character(9) :: CFMT

		integer, intent(in) :: nm_row
		integer, intent(in) :: nm_column

		!成功したことを表示させるか否か
		logical :: message = .false.

		integer i, j


		open(fid,file=filename, status='replace')
		CALL write_fmt(CFMT,fmt_type,cell_length,decimal_digit)

		do i = 1,nm_row
    
            do j = 1,nm_column-1
                
                write(fid,CFMT,advance='no') real_alloc_2d_array(i,j) !改行なし
                write(fid,fmt='(A)',advance='no') ','
    
            end do

			write(fid,CFMT) real_alloc_2d_array(i,nm_column)

        end do
		
		close(fid)

		if (message) then
			print *, "The array has been successfully written to the file."
		end if

	END SUBROUTINE csv_write

	SUBROUTINE csv_write2(fid,filename,real_alloc_2d_array,fmt_type,cell_length,first_row,final_row,first_column,final_column)

		!二次元配列のインデックスが１始まりではない場合

		!任意の行列数(nm_row,nm_column)を持つ二次元配列real_alloc_2d_arrayの実数値を
		!fmt_type、cell_lengthの書式でfilenameのCSVに書き込む。

		!サブルーチンの使用後にdeallocateする

		integer, intent(in) :: fid 
		character(len=*), intent(in) :: filename
		real(8), allocatable, intent(in) :: real_alloc_2d_array(:,:)

		!CFMT='(fmt_type cell_length.decimal_digit)'
		character(len=1), intent(in) :: fmt_type	!e or f
		integer, intent(in) :: cell_length			!15 for 'f' -> integer_digit = 6
		integer :: decimal_digit = 8				!8 for 'f' (real(8))
		character(9) :: CFMT

		integer, intent(in) :: first_row
		integer, intent(in) :: final_row
		integer, intent(in) :: first_column
		integer, intent(in) :: final_column

		!成功したことを表示させるか否か
		logical :: message = .false.

		integer i, j


		open(fid,file=filename, status='replace')
		CALL write_fmt(CFMT,fmt_type,cell_length,decimal_digit)

		do i = first_row,final_row
    
            do j = first_column,final_column-1
                
                write(fid,CFMT,advance='no') real_alloc_2d_array(i,j) !改行なし
                write(fid,fmt='(A)',advance='no') ','
    
            end do

			write(fid,CFMT) real_alloc_2d_array(i,final_column)

        end do
		
		close(fid)

		if (message) then
			print *, "The array has been successfully written to the file."
		end if

	END SUBROUTINE csv_write2

	SUBROUTINE buffer(fid,buf)

        !buf行のヘッダの読み飛ばし
		integer, intent(in) :: fid
        integer, intent(in) :: buf

        integer :: i

        if (buf==0) then
            continue
        else
            do i = 1,buf
                read(fid, '()') !ヘッダの読み飛ばし
            end do
        end if

	END SUBROUTINE buffer


	SUBROUTINE del_spaces(s)
        character (*), intent (inout) :: s
        character (len=len(s)) tmp
        integer i, j
        j = 1
        do i = 1, len(s)
        if (s(i:i)==' ') cycle
        tmp(j:j) = s(i:i)
        j = j + 1
        end do
        s = tmp(1:j-1)
    END SUBROUTINE del_spaces


    SUBROUTINE write_fmt(CFMT,fmt_type,cell_length,decimal_digit)

        character(len=9), intent(out) :: CFMT
        character(1), intent(in) :: fmt_type
        integer, intent(in) :: cell_length
        integer, intent(in) :: decimal_digit

        integer :: integer_digit

        logical :: message = .false.

        if (fmt_type == 'f') then

            integer_digit = cell_length - decimal_digit - 1     !'.'

            if (integer_digit >= 2) then

                if (message) then
                    print *, "Integer digit can be written until:", integer_digit, &
                    & "        (If negative value exists, written until", integer_digit-1,  ")"
                end if

                write(CFMT,'( "(", (A), i2, ".", i2, ")" )') fmt_type, cell_length, decimal_digit

            else
                stop "User error statement: More cell_length is needed."
            end if

        else if (fmt_type == 'e') then

            !'0.', 'E+xx'
            if (cell_length - decimal_digit >= 7) then

                write(CFMT,'( "(", (A), i2, ".", i2, ")" )') fmt_type, cell_length, decimal_digit

            else
                stop "User error statement: More cell_length is needed."
            end if

        else
            stop "User error statement: fmt_type is wrong."
        end if

    END SUBROUTINE write_fmt
	

END MODULE mod_fileio