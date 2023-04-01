Module Module1
    'This example illustrates how the Int and Fix functions return integer portions of numbers. 
    'In the case of a negative number argument, the Int function returns the first negative 
    'integer less than or equal to the number; the Fix function returns the first negative 
    'integer greater than or equal to the number. The following example requires you to specify 
    'Option Strict Off because implicit conversions from type Double to type Integer are 
    'not allowed under Option Strict On:
    ' 
    'VB

    '    ' This code requires Option Strict Off 
    '    Dim MyNumber As Integer
    'MyNumber = Int(99.8)   ' Returns 99.
    'MyNumber = Fix(99.8)   ' Returns 99.

    'MyNumber = Int(-99.8)  ' Returns -100.
    'MyNumber = Fix(-99.8)  ' Returns -99.

    'MyNumber = Int(-99.2)  ' Returns -100.
    'MyNumber = Fix(-99.2)  ' Returns -99.

    'You can use the CInt function to explicitly convert other data types to type Integer with Option Strict Off. However, CInt rounds to the nearest integer instead of truncating the fractional part of numbers. For example:
    'VB

    'MyNumber = CInt(99.8)    ' Returns 100.
    'MyNumber = CInt(-99.8)   ' Returns -100.
    'MyNumber = CInt(-99.2)   ' Returns -99.

    'You can use the CInt function on the result of a call to Fix or Int to perform explicit conversion to integer without rounding. For example:
    'VB

    'MyNumber = CInt(Fix(99.8))   ' Returns 99.
    'MyNumber = CInt(Int(99.8))   ' Returns 99.


End Module
