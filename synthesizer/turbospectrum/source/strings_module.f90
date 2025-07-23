MODULE strings_module

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: to_upper, to_lower
  
  CHARACTER(len=26), PARAMETER :: uca = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  CHARACTER(len=26), PARAMETER :: lca = 'abcdefghijklmnopqrstuvwxyz'
    
CONTAINS

  CHARACTER FUNCTION to_upper_ch(ch)
    CHARACTER, INTENT(in) :: ch
    INTEGER :: ix
    to_upper_ch = ch
    ix = INDEX(lca,ch)
    IF (ix /= 0) to_upper_ch = uca(ix:ix)
  END FUNCTION to_upper_ch

  FUNCTION to_upper(str) RESULT(rslt)
    CHARACTER(len=*), INTENT(in) :: str
    CHARACTER(len=len(str)) :: rslt
    INTEGER :: ix
    DO ix = 1, LEN(str)
       rslt(ix:ix) = to_upper_ch(str(ix:ix))
    END DO    
  END FUNCTION to_upper

  CHARACTER FUNCTION to_lower_ch(ch)
    CHARACTER, INTENT(in) :: ch
    INTEGER :: ix
    to_lower_ch = ch
    ix = INDEX(uca,ch)
    IF (ix /= 0) to_lower_ch = lca(ix:ix)
  END FUNCTION to_lower_ch

  FUNCTION to_lower(str) RESULT(rslt)
    CHARACTER(len=*), INTENT(in) :: str
    CHARACTER(len=len(str)) :: rslt
    INTEGER :: ix
    DO ix = 1, LEN(str)
       rslt(ix:ix) = to_lower_ch(str(ix:ix))
    END DO
  END FUNCTION to_lower

END MODULE strings_module
