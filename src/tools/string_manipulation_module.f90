!> String manipulation module
!!
!! This module contains the subroutines and functions for manipulating strings
!!
!! Date         : 2025
!!
module string_manipulation_module

  implicit none
  public lower,             &
         sweep_blanks,      &
         get_n_in_line,     &
         commented_line,    &
         move_to_character, &
         go_to_string,      &
         get_number_lines,  &
         substitute_string, &
         get_keyword_and_value, &
         validate_and_parse_char, &
         validate_and_parse_int,  &
         validate_and_parse_real, &
         validate_and_parse_bool

contains

   !> Function for lowering all characters of a string
   !!    Input  : str        - String
   !!    Output : str_lower  - Lower Case String
   function lower(str) result(str_lower)
   
      implicit none
   
      !input variables
      character :: str*(*)
      character(len=len(str)) :: str_lower

      !internal variables
      integer :: i
      
      str_lower = str
      do i = 1, len(str)
        select case(str(i:i))
          case("A":"Z")
            str_lower(i:i) = achar(iachar(str_lower(i:i))+32)
        end select
      end do  
       
   end function lower


   !> Function for removing spaces into a string
   !! Readapated from https://www.tek-tips.com/viewthread.cfm?qid=1629864
   !!    Input  : in_str
   !!    Output : sweep_blanks
   character(200) function sweep_blanks(in_str)
  
      !input variables
      character(*), intent(in) :: in_str
      character(200) :: out_str

      !internal variables
      character :: ch
      integer :: j
       
      out_str = " "
      do j=1, len_trim(in_str)
         !get j-th char
         ch = in_str(j:j)
         if (ch .ne. " ".and.ch.ne.char(9)) then
           out_str = trim(out_str) // ch
         endif
         sweep_blanks = out_str
      end do
  
   end function sweep_blanks


   !> Function for getting how many elements are in a line
   !!    Input  : line
   !!    Output : number_
   integer function get_n_in_line(line) result(number_)
  
      implicit none
  
      !input variables
      character(len=200), intent(in) :: line
      
      !internal variables
      character(len=200) :: line_nospaces
  
      line_nospaces = sweep_blanks(line)
      number_ = len_trim(line) - len_trim(line_nospaces) + 1
  
   end function get_n_in_line

   !> Function for understanding if a line is commented (with ! or #)
   !!    Input  : line
   !!    Output : commented (T or F)
   logical function commented_line(line) result(commented)
 
      implicit none
 
      !input variables
      character(len=200), intent(in) :: line

      !internal variables
      character(len=200) :: line_adj
 
      line_adj = trim(adjustl(line))
      if ( line_adj(1:1).eq.'!'.or.line_adj(1:1).eq.'#') then
          commented = .true.
      else
          commented = .false.
      endif
 
   end function commented_line


   !> Function for moving to a specific character 
   !!    Input  : string
   !!    Input  : string2  - the character
   !!    Output : point    - pointer (integer)
   integer function move_to_character(string,string_2) result(point)
   
      implicit none
   
      !input variables
      character(len=200), intent(inout) :: string
      character(len=*), intent(in)      :: string_2
       
      string = adjustl(string)
      point = 1
      do while (point .lt. 200)
         if (string(point:point) .eq. string_2) then
            exit
         else
            point = point + 1
            cycle
         endif
      enddo
   
   end function move_to_character


   !> Subroutine for moving the reading of a file to a specific string
   !!    Input  : IIn           - Input File unit
   !!    Input  : string        - the string to be found
   !!    Output : found_string  - if the string is found (T or F)
   !!    Input  : num_lines     - number of lines to be read
   !!    Output : num_string    - how many strings (optional)
   subroutine go_to_string(IIn,string,found_string,num_lines,num_string)
  
      implicit none
  
      !input variables
      integer, intent(in)            :: IIn !unit 
      character(len=*), intent(in)   :: string       !string to be found
      logical, intent(out)           :: found_string !is found?
      integer, intent(in)            :: num_lines !number of lines to be read
      integer, optional, intent(out) :: num_string!number of strings found
  
      !internal variables
      integer            :: i
      character(len=200) :: line
  
      found_string = .false.
      num_string   = 1

      do i = 1, num_lines
         read(IIn,'(a)') line
         line = lower(line)
         if (trim(adjustl(line)).eq.string) then
            found_string = .true.
            exit
         endif
         num_string = num_string + 1
      enddo
      
   end subroutine go_to_string


   !> Subroutine for getting the number of lines of a file
   !!    Input  : IIn       - Input File unit
   !!    Output : nlines    - number of lines in the file
   subroutine get_number_lines(IIn,nlines)
  
      implicit none
  
      !input variables
      integer, intent(in)  :: IIn
      integer, intent(out) :: nlines
  
      !internal variables
      integer :: io
        
      nlines = 0
      do
        read(IIn,*,iostat=io)
        if (io.ne.0) exit
        nlines = nlines + 1
      end do
    
   end subroutine get_number_lines


   !> Function for substituing a character in a string
   !! Adapted from https://stackoverflow.com/questions/58938347/
   !! how-do-i-replace-a-character-in-the-string-with-another-charater-in-fortran
   !!    Input  : string           - Input File unit
   !!    Input  : search           - character to be looked for
   !!    Input  : substitute       - the substitution
   !!    Output : modified_string  - the modified string
   pure recursive function substitute_string(string,search,substitute) &
                                      result(modified_string)
   
      implicit none
       
      !input variables
      character(len=*), intent(in)  :: string
      character(len=*), intent(in)  :: search
      character(len=*), intent(in)  :: substitute
      character(len=:), allocatable :: modified_string

      !integer variables
      integer :: i
      integer :: stringLen
      integer :: searchLen
       
      stringLen = len(string)
      searchLen = len(search)
      if (stringLen.eq.0 .or. searchLen.eq.0) then
         modified_string = ""
         return
      else if (stringLen.lt.searchLen) then
         modified_string = string
         return
      else 
         i = 1
         do
            if (string(i:i+searchLen-1).eq.search) then
               modified_string = string(1:i-1) // substitute // &
               substitute_string(string(i+searchLen:stringLen), &
                                 search,                        &
                                 substitute)
               exit
            end if
            if (i+searchLen.gt.stringLen) then
               modified_string = string
               exit
            end if
            i = i + 1
            cycle
         end do
      end if
       
   end function substitute_string


   !> Subroutine for extracting keyword and value from a line as strings
   !!    Input  : line            -- line to read
   !!    Output : keyword         -- keyword to read
   !!    Output : value_          -- the result of parsing 
   !!    Input  : out_unit        -- output unit (optional)
   subroutine get_keyword_and_value(out_unit, line, keyword, value_)

      implicit none

      !input/output variables
      integer, intent(in) :: out_unit
      character(len=*), intent(in) :: line
      character(len=*), intent(out) :: keyword
      character(len=*), intent(out) :: value_

      !internal variables
      integer :: pos

      pos = index(line, ":")
      if (pos .eq. 0) then
         write(out_unit, '(a)') "Line "//trim(line)//" has the wrong format.&
                                & This should be -- keyword: value."
         stop
      endif

      keyword = trim(adjustl(line(1:pos-1))) !before :
      value_  = trim(adjustl(line(pos+1:))) !after  :

   end subroutine get_keyword_and_value


   !> Subroutine for validating and parsing input lines for type = character
   !!    Input  : line            -- line to read
   !!    Input  : expected_string -- keyword to read
   !!    Output : value_          -- the result of parsing 
   !!    Input  : out_unit        -- output unit (optional)
   subroutine validate_and_parse_char(line, expected_string, value_, out_unit)

      implicit none

      !input/output variables
      character(len=*), intent(in) :: line
      character(len=*), intent(in) :: expected_string
      character(len=*), intent(out) :: value_
      integer, optional, intent(in) :: out_unit

      !internal variables
      integer :: out_unit_internal
      character(len=200) :: keyword

      if (present(out_unit)) then
          out_unit_internal = out_unit
      else
          out_unit_internal = 6
      endif

      call get_keyword_and_value(out_unit_internal,line,keyword,value_)

      !get the first key name
      if (trim(keyword).ne.trim(expected_string)) then
         write(out_unit, '(a)') "Keyword mismatch: expected " // &
                                 trim(expected_string)//", found "//&
                                 trim(keyword)
         stop
      endif

   end subroutine validate_and_parse_char


   !> Subroutine for validating and parsing input lines for type = integer
   !!    Input  : line            -- line to read
   !!    Input  : expected_string -- keyword to read
   !!    Output : value_          -- the result of parsing 
   !!    Input  : out_unit        -- output unit (optional)
   subroutine validate_and_parse_int(line, expected_string, value_, out_unit)

      implicit none

      !input/output variables
      character(len=*), intent(in) :: line
      character(len=*), intent(in) :: expected_string
      integer, intent(out) :: value_
      integer, optional, intent(in) :: out_unit

      !internal variables
      integer :: out_unit_internal
      integer :: ierr
      character(len=200) :: keyword
      character(len=200) :: str_value

      if (present(out_unit)) then
          out_unit_internal = out_unit
      else
          out_unit_internal = 6
      endif

      call get_keyword_and_value(out_unit_internal,line, keyword, str_value)

      !get the first key name
      if (trim(keyword).ne.trim(expected_string)) then
         write(out_unit_internal, '(a)') "Keyword mismatch: expected " // &
                                          trim(expected_string)//", found "//&
                                          trim(keyword)
         stop
      endif

      read(str_value, *, iostat=ierr) value_
      if (ierr .ne. 0) then
         write(out_unit_internal,'(a)') "Value is not an integer: "// &
                                        trim(str_value)//" for keyword "// &
                                        trim(keyword)
         stop
      endif

   end subroutine validate_and_parse_int


   !> Subroutine for validating and parsing input lines for type = integer
   !!    Input  : line            -- line to read
   !!    Input  : expected_string -- keyword to read
   !!    Output : value_          -- the result of parsing 
   !!    Input  : out_unit        -- output unit (optional)
   subroutine validate_and_parse_real(line, expected_string, value_, out_unit)

      implicit none

      !input/output variables
      character(len=*), intent(in) :: line
      character(len=*), intent(in) :: expected_string
      real(kind=8), intent(out) :: value_
      integer, optional, intent(in) :: out_unit

      !internal variables
      integer :: out_unit_internal
      integer :: ierr
      character(len=200) :: str_value
      character(len=200) :: keyword

      if (present(out_unit)) then
          out_unit_internal = out_unit
      else
          out_unit_internal = 6
      endif

      call get_keyword_and_value(out_unit_internal,line, keyword, str_value)

      !get the first key name
      if (trim(keyword).ne.trim(expected_string)) then
         write(out_unit_internal, '(a)') "Keyword mismatch: expected " // &
                                          trim(expected_string)//", found "//&
                                          trim(adjustl(keyword))
         stop
      endif

      read(str_value, *, iostat=ierr) value_
      if (ierr.ne.0) then
         write(out_unit_internal,'(a)') "Value is not a float: "// &
                                        trim(str_value)//" for keyword "// &
                                        trim(keyword)
         stop
      endif

   end subroutine validate_and_parse_real


   !> Subroutine for validating and parsing input lines for type = boolean
   !!    Input  : line            -- line to read
   !!    Input  : expected_string -- keyword to read
   !!    Output : value_          -- the result of parsing 
   !!    Input  : out_unit        -- output unit (optional)
   subroutine validate_and_parse_bool(line, expected_string, value_, out_unit)

      implicit none

      !input/output variables
      character(len=*), intent(in) :: line
      character(len=*), intent(in) :: expected_string
      logical, intent(out) :: value_
      integer, optional, intent(in) :: out_unit

      !internal variables
      integer :: out_unit_internal
      character(len=200) :: str_value
      character(len=200) :: keyword

      if (present(out_unit)) then
          out_unit_internal = out_unit
      else
          out_unit_internal = 6
      endif

      call get_keyword_and_value(out_unit_internal,line, keyword, str_value)

      !get the first key name
      if (trim(keyword).ne.trim(expected_string)) then
         write(out_unit_internal, '(a)') "Keyword mismatch: expected " // &
                                          trim(expected_string)//", found "//&
                                          trim(adjustl(keyword))
         stop
      endif

      ! Validate and parse the value as True or False
      select case (trim(adjustl(str_value)))
      case ("True")
          value_ = .true.
      case ("False")
          value_ = .false.
      case default
          write(out_unit_internal, '(a)') "Invalid boolean value: "// &
                                           trim(str_value)//&
                                           " for keyword "//trim(keyword)// &
                                           ". Must be True or False."
          stop
      end select

   end subroutine validate_and_parse_bool


   !> Function for getting the nth word in a string
   !!    Input  : str             -- string to read
   !!    Input  : n               -- n-th word 
   !!    Output : w               -- the extracted word
   function nth_word_in_string(str, n) result(w)

      implicit none

      !input/output variables
      character(len=*), intent(in) :: str
      integer, intent(in) :: n
      character(len=200) :: w

      !internal variables
      character(len=200) :: trimmed
      integer :: start, end_, count_
   
      trimmed = trim(adjustl(str))
      start = 1
      count_ = 0
   
      do while (start.le.len(trimmed))
          ! skip spaces
          do while (start.le.len(trimmed).and.trimmed(start:start).eq. " ")
              start = start + 1
          end do
   
          if (start.gt.len(trimmed)) exit  !End of string
   
          !find the end of the word
          end_ = start
          do while (end_.le.len(trimmed).and.trimmed(end_:end_).ne." ")
              end_ = end_ + 1
          end do
   
          count_ = count_ + 1
          if (count_.eq.n) then
              w = trimmed(start:end_-1)
              return
          end if
   
          start = end_ + 1
      end do
   
      w = ""  ! no word found

   end function nth_word_in_string

end module string_manipulation_module
