subroutine lowercase_conversion(input_string,output_string)

    implicit none
    character(len=*) :: input_string
    character(len=*) :: output_string
    integer :: i, offset

    offset = ichar('a') - ichar('A')

    ! Initialize output string
    output_string = input_string

    ! Convert each character
    do i = 1, len_trim(input_string)
        if (input_string(i:i) >= 'A' .and. input_string(i:i) <= 'Z') then
            ! Convert uppercase to lowercase
            output_string(i:i) = char(ichar(input_string(i:i)) + offset)
        else 
            output_string(i:i) = input_string(i:i)
        endif
    end do

end subroutine lowercase_conversion