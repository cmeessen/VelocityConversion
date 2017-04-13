Attribute VB_Name = "AlphaPT"
Sub ExtracAlpha()
Attribute ExtracAlpha.VB_ProcData.VB_Invoke_Func = " \n14"
'
' Creates table for Alpha(P,T) for defined range in Pressure and Temperature
' Before execution make sure to have a worksheet named "AlphaPT"
' Pressure in GPa
' Temperature in K

'

    Dim i, nP, nT, Col_i As Integer
    Dim dP, Pmin, Pmax, P_i As Double
    Dim dT, Tmin, Tmax, T_i As Double
    Dim progress As Double

    nP = 100
    nT = 100

    Pmin = 0.5
    Pmax = 10

    Tmin = 300
    Tmax = 3000

    dP = (Pmax - Pmin) / nP
    dT = (Tmax - Tmin) / nT

    P_i = Pmin
    T_i = Tmin

    Row_i = 2

    For i = 0 To nP
        For j = 0 To nT
            P_i = Pmin + i * dP
            T_i = Tmin + j * dT
            Worksheets("rocks").Range("B59").Value = P_i
            Worksheets("rocks").Range("B60").Value = T_i
            Worksheets("AlphaPT").Cells(Row_i, 1).Value = P_i
            Worksheets("AlphaPT").Cells(Row_i, 2).Value = T_i
            Call rockProps.HackerandAbers03
            Worksheets("AlphaPT").Range(Cells(Row_i, 3), Cells(Row_i, 54)).Value = _
            WorksheetFunction.Transpose(Worksheets("minerals").Range("H5:H56").Value)
            Row_i = Row_i + 1
        Next j

    Next i
End Sub
