Public Class Form1

    ' Return a word representation of a dollar amount.
    Private Function MillionToString(ByVal money As Double, Optional ByVal Before_Amount As String = " Taka ", Optional ByVal after_Paisa As String = " and Paisa ") As String
        ' Get the dollar and cents parts.
        Dim dollars As Double = Int(money)
        Dim cents As Double = CInt(100 * (money - dollars))

        ' Convert the parts into words.
        Dim dollars_str As String = NumberToStringMillion(dollars)
        Dim cents_str As String = NumberToStringMillion(cents)

        Return Before_Amount & dollars_str & after_Paisa & cents_str
    End Function

    ' Return a word representation of the whole number value.
    Private Function NumberToStringMillion(ByVal num As Double) As String
        ' Remove any fractional part.
        num = Int(num)

        ' If the number is 0, return zero.
        If num = 0 Then Return "zero"

        Static groups() As String = {"", "thousand", "million", "billion", "trillion", "quadrillion", "?", "??", "???", "????"}
        Dim result As String = ""

        ' Process the groups, smallest first.
        Dim quotient As Double
        Dim remainder As Integer
        Dim group_num As Integer = 0
        Do While num > 0
            ' Get the next group of three digits.
            quotient = Int(num / 1000)
            remainder = CInt(num - quotient * 1000)
            num = quotient

            ' Convert the group into words.
            result = GroupToWordsMillion(remainder) & _
                " " & groups(group_num) & ", " & _
                result

            ' Get ready for the next group.
            group_num += 1
        Loop

        ' Remove the trailing ", ".
        If result.EndsWith(", ") Then
            result = result.Substring(0, result.Length - 2)
        End If

        Return result.Trim()
    End Function

    ' Convert a number between 0 and 999 into words.
    Private Function GroupToWordsMillion(ByVal num As Integer) As String
        Static one_to_nineteen() As String = {"zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen", "seventeen", "eightteen", "nineteen"}
        Static multiples_of_ten() As String = {"twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety"}

        ' If the number is 0, return an empty string.
        If num = 0 Then Return ""

        ' Handle the hundreds digit.
        Dim digit As Integer
        Dim result As String = ""
        If num > 99 Then
            digit = num \ 100
            num = num Mod 100
            result = one_to_nineteen(digit) & " hundred"
        End If

        ' If num = 0, we have hundreds only.
        If num = 0 Then Return result.Trim()

        ' See if the rest is less than 20.
        If num < 20 Then
            ' Look up the correct name.
            result &= " " & one_to_nineteen(num)
        Else
            ' Handle the tens digit.
            digit = num \ 10
            num = num Mod 10
            result &= " " & multiples_of_ten(digit - 2)

            ' Handle the final digit.
            If num > 0 Then
                result &= " " & one_to_nineteen(num)
            End If
        End If

        Return result.Trim()
    End Function

    Private Sub Form1_Load(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles MyBase.Load
        RadioButton1.Checked = True
        
        Me.Show()
        TextBox2.Focus()

    End Sub

    ' Return a word representation of a dollar amount.
    Private Function DollarsToString(ByVal money As Double) As String
        ' Get the dollar and cents parts.
        Dim dollars As Double = Int(money)
        Dim cents As Double = CInt(100 * (money - dollars))

        ' Convert the parts into words.
        Dim crore_str As String = NumberToString7plus(dollars)
        Debug.Print("NumberToString7plus(dollars) " + NumberToString7plus(dollars).ToString)
        Dim dollars_str As String = NumberToString7(dollars)
        Debug.Print("NumberToString7(dollars) " + NumberToString7(dollars).ToString)
        'Decimal.Round(cents / 100, 2, MidpointRounding.AwayFromZero)
        Dim cents_str As String = ""
        If GroupToWords(100 * Decimal.Round(cents / 100, 2, MidpointRounding.AwayFromZero)) <> "" Then
            cents_str = " and Paisa " + GroupToWords(100 * Decimal.Round(cents / 100, 2, MidpointRounding.AwayFromZero))
        End If
        'Dim cents_str As String = GroupToWords(cents)
        Dim s As String = ""
        If dollars_str <> "" Then
            s = "Taka " & crore_str & dollars_str
        End If
        s = s + cents_str
        'myString.Replace(",", string.Empty)
        s = s.Replace(",", String.Empty)
        Return s
    End Function
    ' Convert a number between 0 and 999 into words.
    Private Function GroupToWords(ByVal num As Integer) As String
        Static one_to_nineteen() As String = {"zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen", "seventeen", "eightteen", "nineteen"}
        Static multiples_of_ten() As String = {"twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety"}

        ' If the number is 0, return an empty string.
        If num = 0 Then Return ""

        ' Handle the hundreds digit.
        Dim digit As Integer
        Dim result As String = ""
        If num > 99 Then
            digit = num \ 100
            num = num Mod 100
            result = one_to_nineteen(digit) & " hundred"
        End If
        '==================================================================
        ' Learn...............
        'The result is the remainder after number1 is divided by number2. For example, the expression 14 Mod 4 evaluates to 2.
        'If number2 evaluates to zero, the behavior of the Mod operator depends on the data type of the operands. An integral division throws a DivideByZeroException exception. A floating-point division returns NaN.
        'If number1 or number2 evaluates to Nothing, it is treated as zero.

        ' The \ Operator (Visual Basic) returns the integer quotient (vagfol in round) of a division. For example, the expression 14 \ 4 evaluates to 3.

        'The / Operator (Visual Basic) returns the full quotient, including the remainder, as a floating-point number. For example, the expression 14 / 4 evaluates to 3.5.
        'Dim testResult As Double
        'testResult = 10 Mod 5
        'testResult = 10 Mod 3
        'testResult = 12 Mod 4.3
        'testResult = 12.6 Mod 5
        'testResult = 47.9 Mod 9.35

        'firstResult = 2.0 Mod 0.2
        '' Double operation returns 0.2, not 0.
        'secondResult = 2D Mod 0.2D
        '' Decimal operation returns 0.
        ' End learn...............
        '============================================================

        'string.Trim() remove leading and trailing spaces !!!!!!!!!!!!!!!!!!
        ' If num = 0, we have hundreds only.
        If num = 0 Then Return result.Trim()

        ' See if the rest is less than 20.
        If num < 20 Then
            ' Look up the correct name.
            result &= " " & one_to_nineteen(num)
        Else
            ' Handle the tens digit.
            digit = num \ 10
            num = num Mod 10
            result &= " " & multiples_of_ten(digit - 2)

            ' Handle the final digit.
            If num > 0 Then
                result &= " " & one_to_nineteen(num)
            End If
        End If

        Return result.Trim()
    End Function
    ' Return a word representation of the whole number value upto lakh part.
    Private Function NumberToString7(ByVal num As Double) As String

        ' Remove any fractional part.
        num = Int(num)
        Dim ss, s1, s2, s3, s4 As String
        s1 = ""

        ss = num.ToString

        If ss.Trim() <> "" Then
            s1 = Microsoft.VisualBasic.Right(ss, 7)
        End If
        num = Integer.Parse(s1)
        Debug.Print(num)
        Debug.Print("num7 ")
        Debug.Print("upto7 num " + num.ToString)

        Dim dig As Integer = 0
        Dim n As Integer = (num.ToString).Length
        Debug.Print("length " + n.ToString)
        If n <= 3 Then
            dig = 1
        End If
        If n <= 5 And n > 3 Then
            dig = 2
        End If
        If n <= 7 And n > 5 Then
            dig = 3
        End If

        ' If the number is 0, return zero.
        If num = 0 Then Return "" ' "zero"
        'Static groups() As String = {"", "thousand", "million", "billion", "trillion", "quadrillion", "?", "??", "???", "????"}
        Static groups() As String = {"", "thousand", "lac", "crore", "?", "??", "???", "????"}
        Dim result As String = ""

        ' Process the groups, smallest first.
        Dim i, j, k As Integer

        Dim quotient As Double
        Dim remainder As Integer
        Dim group_num As Integer = 0
        'Do While num > 0
        '//If s1.Length <= 9 Then
        For i = 1 To dig
            ' Get the next group of three digits.
            If i = 1 Then
                k = 1000
            Else
                k = 100

            End If
            quotient = Int(num / k)
            remainder = CInt(num - quotient * k)
            num = quotient
            'Console.WriteLine("upto 7 i =" + i.ToString + "quotient =" + quotient.ToString + "remainder =" + remainder.ToString)

            ' Convert the group into words.
            result = " " & GroupToWords(remainder) & " " & groups(group_num) & result & " "

            Debug.Print("upto 7 i =" + i.ToString + "group_num =" + group_num.ToString + "quotient =" + quotient.ToString + "remainder =" + remainder.ToString)
            Debug.Print("upto 7result =" + result.ToString)

            ' Get ready for the next group.
            group_num += 1

        Next i
        'Debug.Print("group_num =0,groups(group_num)=" + groups(0).ToString)
        'Debug.Print("group_num =1,groups(group_num)=" + groups(1).ToString)


        '//End If

        'Loop
        '=======================================
        '        Dim literal As String = "CatDogFence"
        '        literal.Substring(0, 3)=>Cat
        '           literal.Substring(6)=>Fence
        '"me@me.com".Substring(5, 4)=.com at the end (start at position 5 in the string and grab 4 characters). 
        'Email = "me@me.con"
        '"me@me.com".Substring(Email.Length - 4, 4 )=>The starting position for Substring( ) this time is "Email.Length - 4". This is the length of the string variable called Email, minus 4 characters. The other 4 means "grab four characters"
        'FirstName.Substring( i, j )
        '============================================


        ' Remove the trailing ", ".
        'If result.EndsWith(", ") Then
        '    result = result.Substring(0, result.Length - 2)
        'End If

        Return result.Trim()
    End Function

    ' Return a word representation of the whole number value above lakh part.
    Private Function NumberToString7plus(ByVal num As Double) As String

        ' Remove any fractional part.
        num = Int(num)

        Dim ss, s1, s2, s3, s4 As String
        's1 = ""
        'ss = num.ToString
        'If ss.Length > 7 Then
        '    s1 = ss.Substring(0, ss.Length - 7)
        'End If

        'If s1.Trim() <> "" Then
        '    num = Integer.Parse(s1)
        'End If
        Dim quotient As Double
        quotient = Int(num / 10000000)
        num = quotient

        Debug.Print("7plus num " + num.ToString)
        Dim dig As Integer = 0
        Dim n As Integer = (num.ToString).Length
        Debug.Print("length " + n.ToString)
        If n <= 3 Then
            dig = 1
        End If
        If n <= 5 And n > 3 Then
            dig = 2
        End If
        If n <= 7 And n > 5 Then
            dig = 3
        End If
        'Dim text As String = "1234"
        'Dim i As Integer = Convert.ToInt32(text)
        'Dim i As Integer = Integer.Parse(text)

        ' If the number is 0, return zero.
        If num = 0 Then Return "" ' "zero"
        'Static groups() As String = {"", "thousand", "million", "billion", "trillion", "quadrillion", "?", "??", "???", "????"}
        Static groups() As String = {"", "thousand", "lac", "crore", "?", "??", "???", "????"}
        Dim result As String = ""

        ' Process the groups, smallest first.
        Dim i, j, k As Integer


        Dim remainder As Integer
        Dim group_num As Integer = 0
        'Do While num > 0
        '//If s1.Length <= 9 Then
        For i = 1 To dig
            ' Get the next group of three digits.
            If i = 1 Then
                k = 1000
            Else
                k = 100

            End If
            quotient = Int(num / k)
            remainder = CInt(num - quotient * k)
            num = quotient

            ' Convert the group into words.
            result = " " & GroupToWords(remainder) & " " & groups(group_num) & result & " "

            Debug.Print("7plus i =" + i.ToString + "group_num =" + group_num.ToString + "quotient =" + quotient.ToString + "remainder =" + remainder.ToString)
            Debug.Print("7plus result =" + result.ToString)

            ' Get ready for the next group.
            group_num += 1

        Next i
        '//End If

        'Loop
        '=======================================
        '        Dim literal As String = "CatDogFence"
        '        literal.Substring(0, 3)=>Cat
        '           literal.Substring(6)=>Fence
        '"me@me.com".Substring(5, 4)=.com at the end (start at position 5 in the string and grab 4 characters). 
        'Email = "me@me.con"
        '"me@me.com".Substring(Email.Length - 4, 4 )=>The starting position for Substring( ) this time is "Email.Length - 4". This is the length of the string variable called Email, minus 4 characters. The other 4 means "grab four characters"
        'FirstName.Substring( i, j )
        '============================================


        ' Remove the trailing ", ".
        'If result.EndsWith(", ") Then
        '    result = result.Substring(0, result.Length - 2)
        'End If

        Return result.Trim() + " crore "
    End Function


    'END Million.......................................................
    '======================================================================================


















    ' Function for conversion of a Indian Taka into words
    '   Parameter - accept a Currency
    '   Returns the number in words format
    '   You can use this function in Excel, VBA, VB6,.NET
    '====================================================

    '****************************************************
    ' Code Created by Bharat Modha 
    ' Porbandar (Gujarat)-India
    ' Email : bharatmodha@yahoo.com
    '****************************************************

    Function TakaToWord(ByVal MyNumber)
        Dim Temp
        Dim Taka, Paisa As String
        Dim DecimalPlace, iCount
        Dim Hundred, Words As String
        Dim place(9) As String
        place(0) = " Thousand "
        place(2) = " Lakh "
        place(4) = " Crore "
        place(6) = " Arab "
        place(8) = " Kharab "
        On Error Resume Next
        ' Convert MyNumber to a string, trimming extra spaces.
        MyNumber = Trim(Str(MyNumber))

        ' Find decimal place.
        DecimalPlace = InStr(MyNumber, ".")

        ' If we find decimal place...
        If DecimalPlace > 0 Then
            ' Convert Paisa
            Temp = Microsoft.VisualBasic.Left(Mid(MyNumber, DecimalPlace + 1) & "00", 2)
            Paisa = " and " & ConvertTens(Temp) & " Paisa"

            ' Strip off paisa from remainder to convert.
            MyNumber = Trim(Microsoft.VisualBasic.Left(MyNumber, DecimalPlace - 1))
        End If

        '===============================================================
        Dim TM As String  ' If MyNumber between Rs.1 To 99 Only.
        TM = Microsoft.VisualBasic.Right(MyNumber, 2)

        If Len(MyNumber) > 0 And Len(MyNumber) <= 2 Then
            If Len(TM) = 1 Then
                Words = ConvertDigit(TM)
                TakaToWord = "Taka " & Words & Paisa & " Only"

                Exit Function

            Else
                If Len(TM) = 2 Then
                    Words = ConvertTens(TM)
                    TakaToWord = "Taka " & Words & Paisa & " Only"
                    Exit Function

                End If
            End If
        End If
        '===============================================================


        ' Convert last 3 digits of MyNumber to ruppees in word.
        Hundred = ConvertHundred(Microsoft.VisualBasic.Right(MyNumber, 3))
        ' Strip off last three digits
        MyNumber = Microsoft.VisualBasic.Left(MyNumber, Len(MyNumber) - 3)

        iCount = 0
        Do While MyNumber <> ""
            'Strip last two digits
            Temp = Microsoft.VisualBasic.Right(MyNumber, 2)
            If Len(MyNumber) = 1 Then


                If Trim(Words) = "Thousand" Or _
                Trim(Words) = "Lakh  Thousand" Or _
                Trim(Words) = "Lakh" Or _
                Trim(Words) = "Crore" Or _
                Trim(Words) = "Crore  Lakh  Thousand" Or _
                Trim(Words) = "Arab  Crore  Lakh  Thousand" Or _
                Trim(Words) = "Arab" Or _
                Trim(Words) = "Kharab  Arab  Crore  Lakh  Thousand" Or _
                Trim(Words) = "Kharab" Then

                    Words = ConvertDigit(Temp) & place(iCount)
                    MyNumber = Microsoft.VisualBasic.Left(MyNumber, Len(MyNumber) - 1)

                Else

                    Words = ConvertDigit(Temp) & place(iCount) & Words
                    MyNumber = Microsoft.VisualBasic.Left(MyNumber, Len(MyNumber) - 1)

                End If
            Else

                If Trim(Words) = "Thousand" Or _
                   Trim(Words) = "Lakh  Thousand" Or _
                   Trim(Words) = "Lakh" Or _
                   Trim(Words) = "Crore" Or _
                   Trim(Words) = "Crore  Lakh  Thousand" Or _
                   Trim(Words) = "Arab  Crore  Lakh  Thousand" Or _
                   Trim(Words) = "Arab" Then


                    Words = ConvertTens(Temp) & place(iCount)


                    MyNumber = Microsoft.VisualBasic.Left(MyNumber, Len(MyNumber) - 2)
                Else

                    '=================================================================
                    ' if only Lakh, Crore, Arab, Kharab

                    If Trim(ConvertTens(Temp) & place(iCount)) = "Lakh" Or _
                       Trim(ConvertTens(Temp) & place(iCount)) = "Crore" Or _
                       Trim(ConvertTens(Temp) & place(iCount)) = "Arab" Then

                        Words = Words
                        MyNumber = Microsoft.VisualBasic.Left(MyNumber, Len(MyNumber) - 2)
                    Else
                        Words = ConvertTens(Temp) & place(iCount) & Words
                        MyNumber = Microsoft.VisualBasic.Left(MyNumber, Len(MyNumber) - 2)
                    End If

                End If
            End If

            iCount = iCount + 2
        Loop

        TakaToWord = "Taka " & Words & Hundred & Paisa & " Only"
    End Function

    ' Conversion for hundreds
    '*****************************************
    Private Function ConvertHundred(ByVal MyNumber)
        Dim Result As String

        ' Exit if there is nothing to convert.
        If Val(MyNumber) = 0 Then Exit Function

        ' Append leading zeros to number.
        MyNumber = Microsoft.VisualBasic.Right("000" & MyNumber, 3)

        ' Do we have a hundreds place digit to convert?
        If Microsoft.VisualBasic.Left(MyNumber, 1) <> "0" Then
            Result = ConvertDigit(Microsoft.VisualBasic.Left(MyNumber, 1)) & " Hundred "
        End If

        ' Do we have a tens place digit to convert?
        If Mid(MyNumber, 2, 1) <> "0" Then
            Result = Result & ConvertTens(Mid(MyNumber, 2))
        Else
            ' If not, then convert the ones place digit.
            Result = Result & ConvertDigit(Mid(MyNumber, 3))
        End If

        ConvertHundred = Trim(Result)
    End Function

    ' Conversion for tens
    '*****************************************
    Private Function ConvertTens(ByVal MyTens)
        Dim Result As String

        ' Is value between 10 and 19?
        If Val(Microsoft.VisualBasic.Left(MyTens, 1)) = 1 Then
            Select Case Val(MyTens)
                Case 10 : Result = "Ten"
                Case 11 : Result = "Eleven"
                Case 12 : Result = "Twelve"
                Case 13 : Result = "Thirteen"
                Case 14 : Result = "Fourteen"
                Case 15 : Result = "Fifteen"
                Case 16 : Result = "Sixteen"
                Case 17 : Result = "Seventeen"
                Case 18 : Result = "Eighteen"
                Case 19 : Result = "Nineteen"
                Case Else
            End Select
        Else
            ' .. otherwise it's between 20 and 99.
            Select Case Val(Microsoft.VisualBasic.Left(MyTens, 1))
                Case 2 : Result = "Twenty "
                Case 3 : Result = "Thirty "
                Case 4 : Result = "Forty "
                Case 5 : Result = "Fifty "
                Case 6 : Result = "Sixty "
                Case 7 : Result = "Seventy "
                Case 8 : Result = "Eighty "
                Case 9 : Result = "Ninety "
                Case Else
            End Select

            ' Convert ones place digit.
            Result = Result & ConvertDigit(Microsoft.VisualBasic.Right(MyTens, 1))
        End If

        ConvertTens = Result
    End Function

    Private Function ConvertDigit(ByVal MyDigit)
        Select Case Val(MyDigit)
            Case 1 : ConvertDigit = "One"
            Case 2 : ConvertDigit = "Two"
            Case 3 : ConvertDigit = "Three"
            Case 4 : ConvertDigit = "Four"
            Case 5 : ConvertDigit = "Five"
            Case 6 : ConvertDigit = "Six"
            Case 7 : ConvertDigit = "Seven"
            Case 8 : ConvertDigit = "Eight"
            Case 9 : ConvertDigit = "Nine"
            Case Else : ConvertDigit = ""
        End Select
    End Function

    Private Sub Button1_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles Button1.Click
        Dim sr As String = "1234567891"
        Debug.Print(sr.Substring(0, sr.Length - 7))

        If RadioButton1.Checked = True Then

            RichTextBox1.Text = DollarsToString(CDbl(TextBox2.Text))  'TakaToWord(TextBox2.Text)

        End If

        If RadioButton2.Checked = True Then
            ' Convert into words.

            RichTextBox1.Text = MillionToString(CDbl(TextBox2.Text))

        End If

    End Sub


    'Start MICROSOFT 1............................

    ''http://support2.microsoft.com/default.aspx?scid=kb;en-us;213360
    'Option Explicit
    ''Main Function
    'Function SpellNumber(ByVal MyNumber)
    '    Dim Dollars, Cents, Temp
    '    Dim DecimalPlace, Count
    '    Dim Place(9) As String
    '    Place(2) = " Thousand "
    '    Place(3) = " Million "
    '    Place(4) = " Billion "
    '    Place(5) = " Trillion "
    '    ' String representation of amount.
    '    MyNumber = Trim(Str(MyNumber))
    '    ' Position of decimal place 0 if none.
    '    DecimalPlace = InStr(MyNumber, ".")
    '    ' Convert cents and set MyNumber to dollar amount.
    '    If DecimalPlace > 0 Then
    '        Cents = GetTens(Microsoft.VisualBasic.Left(Mid(MyNumber, DecimalPlace + 1) & _
    '                  "00", 2))
    '        MyNumber = Trim(Microsoft.VisualBasic.Left(MyNumber, DecimalPlace - 1))
    '    End If
    '    Count = 1
    '    Do While MyNumber <> ""
    '        Temp = GetHundreds(Microsoft.VisualBasic.Right(MyNumber, 3))
    '        If Temp <> "" Then Dollars = Temp & Place(Count) & Dollars
    '        If Len(MyNumber) > 3 Then
    '            MyNumber = Microsoft.VisualBasic.Left(MyNumber, Len(MyNumber) - 3)
    '        Else
    '            MyNumber = ""
    '        End If
    '        Count = Count + 1
    '    Loop
    '    Select Case Dollars
    '        Case ""
    '            Dollars = "No Dollars"
    '        Case "One"
    '            Dollars = "One Dollar"
    '        Case Else
    '            Dollars = Dollars & " Dollars"
    '    End Select
    '    Select Case Cents
    '        Case ""
    '            Cents = " and No Cents"
    '        Case "One"
    '            Cents = " and One Cent"
    '        Case Else
    '            Cents = " and " & Cents & " Cents"
    '    End Select
    '    SpellNumber = Dollars & Cents
    'End Function

    '' Converts a number from 100-999 into text 
    'Function GetHundreds(ByVal MyNumber)
    '    Dim Result As String
    '    If Val(MyNumber) = 0 Then Exit Function
    '    MyNumber = Microsoft.VisualBasic.Right("000" & MyNumber, 3)
    '    ' Convert the hundreds place.
    '    If Mid(MyNumber, 1, 1) <> "0" Then
    '        Result = GetDigit(Mid(MyNumber, 1, 1)) & " Hundred "
    '    End If
    '    ' Convert the tens and ones place.
    '    If Mid(MyNumber, 2, 1) <> "0" Then
    '        Result = Result & GetTens(Mid(MyNumber, 2))
    '    Else
    '        Result = Result & GetDigit(Mid(MyNumber, 3))
    '    End If
    '    GetHundreds = Result
    'End Function

    '' Converts a number from 10 to 99 into text. 
    'Function GetTens(ByVal TensText)
    '    Dim Result As String
    '    Result = ""           ' Null out the temporary function value.
    '    If Val(Microsoft.VisualBasic.Left(TensText, 1)) = 1 Then   ' If value between 10-19...
    '        Select Case Val(TensText)
    '            Case 10 : Result = "Ten"
    '            Case 11 : Result = "Eleven"
    '            Case 12 : Result = "Twelve"
    '            Case 13 : Result = "Thirteen"
    '            Case 14 : Result = "Fourteen"
    '            Case 15 : Result = "Fifteen"
    '            Case 16 : Result = "Sixteen"
    '            Case 17 : Result = "Seventeen"
    '            Case 18 : Result = "Eighteen"
    '            Case 19 : Result = "Nineteen"
    '            Case Else
    '        End Select
    '    Else                                 ' If value between 20-99...
    '        Select Case Val(Microsoft.VisualBasic.Left(TensText, 1))
    '            Case 2 : Result = "Twenty "
    '            Case 3 : Result = "Thirty "
    '            Case 4 : Result = "Forty "
    '            Case 5 : Result = "Fifty "
    '            Case 6 : Result = "Sixty "
    '            Case 7 : Result = "Seventy "
    '            Case 8 : Result = "Eighty "
    '            Case 9 : Result = "Ninety "
    '            Case Else
    '        End Select
    '        Result = Result & GetDigit _
    '            (Microsoft.VisualBasic.Right(TensText, 1))  ' Retrieve ones place.
    '    End If
    '    GetTens = Result
    'End Function

    '' Converts a number from 1 to 9 into text. 
    'Function GetDigit(ByVal Digit)
    '    Select Case Val(Digit)
    '        Case 1 : GetDigit = "One"
    '        Case 2 : GetDigit = "Two"
    '        Case 3 : GetDigit = "Three"
    '        Case 4 : GetDigit = "Four"
    '        Case 5 : GetDigit = "Five"
    '        Case 6 : GetDigit = "Six"
    '        Case 7 : GetDigit = "Seven"
    '        Case 8 : GetDigit = "Eight"
    '        Case 9 : GetDigit = "Nine"
    '        Case Else : GetDigit = ""
    '    End Select
    'End Function


    'end microsoft.....................................
    ' kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
    'Start MICROSOFT 2............................


    'To create the ConvertCurrencyToEnglish() function, follow these steps:
    '2nd method microsoft ..............................................................................
    'http://support2.microsoft.com/?scid=http://support.microsoft.com:80/support/kb/articles/Q95/6/40.ASP
    'Create a new module and type the following line in the Declarations section if the line is not already there:
    'Option Explicit
    'Type the following four procedures:

    'Function ConvertCurrencyToEnglish(ByVal MyNumber)
    '    Dim Temp
    '    Dim Dollars, Cents
    '    Dim DecimalPlace, Count

    '    Dim Place(9) As String
    '    Place(2) = " Thousand "
    '    Place(3) = " Million "
    '    Place(4) = " Billion "
    '    Place(5) = " Trillion "

    '    ' Convert MyNumber to a string, trimming extra spaces.
    '    MyNumber = Trim(Str(MyNumber))

    '    ' Find decimal place.
    '    DecimalPlace = InStr(MyNumber, ".")

    '    ' If we find decimal place...
    '    If DecimalPlace > 0 Then
    '        ' Convert cents
    '        Temp = Microsoft.VisualBasic.Left(Mid(MyNumber, DecimalPlace + 1) & "00", 2)
    '        Cents = ConvertTens(Temp)

    '        ' Strip off cents from remainder to convert.
    '        MyNumber = Trim(Microsoft.VisualBasic.Left(MyNumber, DecimalPlace - 1))
    '    End If

    '    Count = 1
    '    Do While MyNumber <> ""
    '        ' Convert last 3 digits of MyNumber to English dollars.
    '        Temp = ConvertHundreds(Microsoft.VisualBasic.Right(MyNumber, 3))
    '        If Temp <> "" Then Dollars = Temp & Place(Count) & Dollars
    '        If Len(MyNumber) > 3 Then
    '            ' Remove last 3 converted digits from MyNumber.
    '            MyNumber = Microsoft.VisualBasic.Left(MyNumber, Len(MyNumber) - 3)
    '        Else
    '            MyNumber = ""
    '        End If
    '        Count = Count + 1
    '    Loop

    '    ' Clean up dollars.
    '    Select Case Dollars
    '        Case ""
    '            Dollars = "No Dollars"
    '        Case "One"
    '            Dollars = "One Dollar"
    '        Case Else
    '            Dollars = Dollars & " Dollars"
    '    End Select

    '    ' Clean up cents.
    '    Select Case Cents
    '        Case ""
    '            Cents = " And No Cents"
    '        Case "One"
    '            Cents = " And One Cent"
    '        Case Else
    '            Cents = " And " & Cents & " Cents"
    '    End Select

    '    ConvertCurrencyToEnglish = Dollars & Cents
    'End Function

    'Private Function ConvertHundreds(ByVal MyNumber)
    '    Dim Result As String

    '    ' Exit if there is nothing to convert.
    '    If Val(MyNumber) = 0 Then Exit Function

    '    ' Append leading zeros to number.
    '    MyNumber = Microsoft.VisualBasic.Right("000" & MyNumber, 3)

    '    ' Do we have a hundreds place digit to convert?
    '    If Microsoft.VisualBasic.Left(MyNumber, 1) <> "0" Then
    '        Result = ConvertDigit(Microsoft.VisualBasic.Left(MyNumber, 1)) & " Hundred "
    '    End If

    '    ' Do we have a tens place digit to convert?
    '    If Mid(MyNumber, 2, 1) <> "0" Then
    '        Result = Result & ConvertTens(Mid(MyNumber, 2))
    '    Else
    '        ' If not, then convert the ones place digit.
    '        Result = Result & ConvertDigit(Mid(MyNumber, 3))
    '    End If

    '    ConvertHundreds = Trim(Result)
    'End Function

    'Private Function ConvertTens(ByVal MyTens)
    '    Dim Result As String

    '    ' Is value between 10 and 19?
    '    If Val(Microsoft.VisualBasic.Left(MyTens, 1)) = 1 Then
    '        Select Case Val(MyTens)
    '            Case 10 : Result = "Ten"
    '            Case 11 : Result = "Eleven"
    '            Case 12 : Result = "Twelve"
    '            Case 13 : Result = "Thirteen"
    '            Case 14 : Result = "Fourteen"
    '            Case 15 : Result = "Fifteen"
    '            Case 16 : Result = "Sixteen"
    '            Case 17 : Result = "Seventeen"
    '            Case 18 : Result = "Eighteen"
    '            Case 19 : Result = "Nineteen"
    '            Case Else
    '        End Select
    '    Else
    '        ' .. otherwise it's between 20 and 99.
    '        Select Case Val(Microsoft.VisualBasic.Left(MyTens, 1))
    '            Case 2 : Result = "Twenty "
    '            Case 3 : Result = "Thirty "
    '            Case 4 : Result = "Forty "
    '            Case 5 : Result = "Fifty "
    '            Case 6 : Result = "Sixty "
    '            Case 7 : Result = "Seventy "
    '            Case 8 : Result = "Eighty "
    '            Case 9 : Result = "Ninety "
    '            Case Else
    '        End Select

    '        ' Convert ones place digit.
    '        Result = Result & ConvertDigit(Microsoft.VisualBasic.Right(MyTens, 1))
    '    End If

    '    ConvertTens = Result
    'End Function

    'Private Function ConvertDigit(ByVal MyDigit)
    '    Select Case Val(MyDigit)
    '        Case 1 : ConvertDigit = "One"
    '        Case 2 : ConvertDigit = "Two"
    '        Case 3 : ConvertDigit = "Three"
    '        Case 4 : ConvertDigit = "Four"
    '        Case 5 : ConvertDigit = "Five"
    '        Case 6 : ConvertDigit = "Six"
    '        Case 7 : ConvertDigit = "Seven"
    '        Case 8 : ConvertDigit = "Eight"
    '        Case 9 : ConvertDigit = "Nine"
    '        Case Else : ConvertDigit = ""
    '    End Select
    'End Function

End Class
