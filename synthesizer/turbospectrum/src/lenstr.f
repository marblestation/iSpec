C
       FUNCTION LENSTR(STRING)
*
* Returns the length of a string not counting trailing blanks
*
       CHARACTER*(*)  STRING
*
       DO 10 I = LEN(STRING), 1, -1
         IF(STRING(I:I) .NE. ' ') THEN
           LENSTR = I
           RETURN
         ENDIF
   10  CONTINUE
       LENSTR = 0
       RETURN
       END
         
        
