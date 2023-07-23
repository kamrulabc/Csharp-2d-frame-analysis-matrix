using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Diagnostics;
using System.Windows.Forms;
using System.Linq;
using System.Xml.Linq;


namespace _2member__3member_frame_with_spring
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        //====================================================================





        //T Global to Local axex	|	TT Local to Global axex		
        //c	s	0	            |	c	-s	0
        //-s	c	0	            |	s	c	0
        //0	0	1	            |	0	0	1


        //f1 = f1'cosθ − f2'sinθ
        //f2 = f1'sinθ + f2'cosθ
        //f3 = f3'
        //f4 = f4'cosθ − f5'sinθ
        //f5 = f4'sinθ + f5'cosθ
        //f6 = f6'



        //U´= | cosθ sinθ|  |UGX  >>>>T Global to Local axex
        //V´ =|-sinθ cosθ|  |VGY
        //        T

        //UGX= |cosθ-sinθ|  |U´ >>>TT Local to Global axex	  
        //VGY= |sinθ cosθ|  |V´
        //         TT              >>>TT=  Inv(T)=Transpose (T)


        //T Global to Local axex	          |TT Local to Global axex	
        //Local force=T*Global force          |Global f= TT*f local
        //Local d=T*Global d                  |Global d= TT*d local

        //'   kff   | kfs      | Df  |    | Ff |
        //'---------|------- X |-----| =  |----| ALL GLOBAL AXEX
        //'    ksf  |kss       | Ds  |    |Fs  |











        ///////////////////////////////////////////////////////////////////////////////
        //kkkkkkkkkkkkkkkkkkkkkkkkkkkk
        private int totalElement = 50;
        private int totalnode = 100;
        private double[,] Coordinate = new double[100 + 1, 3]; ////The index of the 2nd dimension, 1=X Coordinate, 2=Y Coordinate
        private int[,] Connectivity = new int[50 + 1, 3]; ////The index of the 2nd dimension, 1=Start Node I, 2=End node J of IJ Member
        //kkkkkkkkkkkkkkkkkkkkkkkkkkk
        //kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
        private string MakeDisplayable_2D2D_ij_0_0(double[,] Array2D)
        {
            int row = 0;
            int col = 0;
            string s1 = null;
            string s2 = null;
            string s3 = null;


            s1 = "";


            for (row = 0; row <= Array2D.GetUpperBound(0); row++)
            {
                for (col = 0; col <= Array2D.GetUpperBound(1); col++)
                {
                    //s1 = s1 + row.ToString + col.ToString + "_" + (Array2D(row, col).ToString) + " "
                    s1 = s1 + (Array2D[row, col].ToString()) + " ";
                }
                s1 = s1 + Environment.NewLine;
            }

            return s1;

        }


        ///////////////////////////////////////////////////////////////////////////////
        private int[] Index_6_From_MemID_1D_i_0(int mem_id)
        {
            int[] Index_6 = new int[6];
            Index_6[0] = Connectivity[mem_id, 1] * 3 - 2;
            Index_6[1] = Connectivity[mem_id, 1] * 3 - 1;
            Index_6[2] = Connectivity[mem_id, 1] * 3;

            Index_6[3] = Connectivity[mem_id, 2] * 3 - 2;
            Index_6[4] = Connectivity[mem_id, 2] * 3 - 1;
            Index_6[5] = Connectivity[mem_id, 2] * 3;
            return Index_6;
        }
        private int[] Index_6_From_MemID_1D_i_1(int mem_id)
        {
            int[] Index_6 = new int[7];
            Index_6[0] = 0;
            Index_6[1] = Connectivity[mem_id, 1] * 3 - 2;
            Index_6[2] = Connectivity[mem_id, 1] * 3 - 1;
            Index_6[3] = Connectivity[mem_id, 1] * 3;

            Index_6[4] = Connectivity[mem_id, 2] * 3 - 2;
            Index_6[5] = Connectivity[mem_id, 2] * 3 - 1;
            Index_6[6] = Connectivity[mem_id, 2] * 3;
            return Index_6;
        }

        private double[,] K_G_Element_ij_11i_0j_0(double L, double A, double I, double Ey, double c, double s)
        {
            double[,] matrixD = new double[7, 7];
            double[,] matrixD_0 = new double[6, 6];
            //Dim matrixD As Double(6,6) = New Double(6,6) {}
            Array.Clear(matrixD, 0, matrixD.Length);
            double k_1 = 0;
            double k_2 = 0;
            double k_3 = 0;
            double k_4 = 0;

            k_1 = Ey * A / L;
            k_2 = 12 * Ey * I / Math.Pow(L, 3);
            k_3 = Ey * I / L;
            k_4 = 6 * Ey * I / Math.Pow(L, 2);

            //1st
            matrixD[1, 1] = Math.Pow(c, 2) * k_1 + Math.Pow(s, 2) * k_2; //row, col
            matrixD[1, 2] = s * c * (k_1 - k_2);
            matrixD[1, 3] = -s * k_4;
            matrixD[1, 4] = -(Math.Pow(c, 2) * k_1 + Math.Pow(s, 2) * k_2);
            matrixD[1, 5] = s * c * (-k_1 + k_2);
            matrixD[1, 6] = -s * k_4;
            //2nd
            matrixD[1, 2] = s * c * (k_1 - k_2);
            matrixD[2, 2] = Math.Pow(s, 2) * k_1 + Math.Pow(c, 2) * k_2;
            matrixD[2, 3] = c * k_4;
            matrixD[2, 4] = s * c * (-k_1 + k_2);
            matrixD[2, 5] = -(Math.Pow(s, 2) * k_1 + Math.Pow(c, 2) * k_2);
            matrixD[2, 6] = c * k_4;
            //3rd
            matrixD[3, 1] = -s * k_4;
            matrixD[3, 2] = c * k_4;
            matrixD[3, 3] = 4 * k_3;
            matrixD[3, 4] = s * k_4;
            matrixD[3, 5] = -c * k_4;
            matrixD[3, 6] = 2 * k_3;
            //4th
            matrixD[4, 1] = -(Math.Pow(c, 2) * k_1 + Math.Pow(s, 2) * k_2);
            matrixD[4, 2] = s * c * (-k_1 + k_2);
            matrixD[4, 3] = s * k_4;
            matrixD[4, 4] = Math.Pow(c, 2) * k_1 + Math.Pow(s, 2) * k_2;
            matrixD[4, 5] = s * c * (k_1 - k_2);
            matrixD[4, 6] = s * k_4;
            //5th
            matrixD[5, 1] = s * c * (-k_1 + k_2);
            matrixD[5, 2] = -(Math.Pow(s, 2) * k_1 + Math.Pow(c, 2) * k_2);
            matrixD[5, 3] = -c * k_4;
            matrixD[5, 4] = s * c * (k_1 - k_2);
            matrixD[5, 5] = Math.Pow(s, 2) * k_1 + Math.Pow(c, 2) * k_2;
            matrixD[5, 6] = -c * k_4;
            //6th
            matrixD[6, 1] = -s * k_4;
            matrixD[6, 2] = c * k_4;
            matrixD[6, 3] = 2 * k_3;
            matrixD[6, 4] = s * k_4;
            matrixD[6, 5] = -c * k_4;
            matrixD[6, 6] = 4 * k_3;



            int j = 0;
            int k = 0;
            for (j = 0; j <= 5; j++)
            {
                for (k = 0; k <= 5; k++)
                {
                    matrixD_0[j, k] = matrixD[j + 1, k + 1];
                }
            }
            return matrixD_0;
        }



        //kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
        private double[,] K_G_Element_ij_11(double L, double A, double I, double Ey, double c, double s)
        {
            double[,] matrixD = new double[7, 7];

            Array.Clear(matrixD, 0, matrixD.Length);
            double k_1 = 0;
            double k_2 = 0;
            double k_3 = 0;
            double k_4 = 0;

            k_1 = Ey * A / L;
            k_2 = 12 * Ey * I / Math.Pow(L, 3);
            k_3 = Ey * I / L;
            k_4 = 6 * Ey * I / Math.Pow(L, 2);

            //1st
            matrixD[1, 1] = Math.Pow(c, 2) * k_1 + Math.Pow(s, 2) * k_2; //row, col
            matrixD[1, 2] = s * c * (k_1 - k_2);
            matrixD[1, 3] = -s * k_4;
            matrixD[1, 4] = -(Math.Pow(c, 2) * k_1 + Math.Pow(s, 2) * k_2);
            matrixD[1, 5] = s * c * (-k_1 + k_2);
            matrixD[1, 6] = -s * k_4;
            //2nd
            matrixD[1, 2] = s * c * (k_1 - k_2);
            matrixD[2, 2] = Math.Pow(s, 2) * k_1 + Math.Pow(c, 2) * k_2;
            matrixD[2, 3] = c * k_4;
            matrixD[2, 4] = s * c * (-k_1 + k_2);
            matrixD[2, 5] = -(Math.Pow(s, 2) * k_1 + Math.Pow(c, 2) * k_2);
            matrixD[2, 6] = c * k_4;
            //3rd
            matrixD[3, 1] = -s * k_4;
            matrixD[3, 2] = c * k_4;
            matrixD[3, 3] = 4 * k_3;
            matrixD[3, 4] = s * k_4;
            matrixD[3, 5] = -c * k_4;
            matrixD[3, 6] = 2 * k_3;
            //4th
            matrixD[4, 1] = -(Math.Pow(c, 2) * k_1 + Math.Pow(s, 2) * k_2);
            matrixD[4, 2] = s * c * (-k_1 + k_2);
            matrixD[4, 3] = s * k_4;
            matrixD[4, 4] = Math.Pow(c, 2) * k_1 + Math.Pow(s, 2) * k_2;
            matrixD[4, 5] = s * c * (k_1 - k_2);
            matrixD[4, 6] = s * k_4;
            //5th
            matrixD[5, 1] = s * c * (-k_1 + k_2);
            matrixD[5, 2] = -(Math.Pow(s, 2) * k_1 + Math.Pow(c, 2) * k_2);
            matrixD[5, 3] = -c * k_4;
            matrixD[5, 4] = s * c * (k_1 - k_2);
            matrixD[5, 5] = Math.Pow(s, 2) * k_1 + Math.Pow(c, 2) * k_2;
            matrixD[5, 6] = -c * k_4;
            //6th
            matrixD[6, 1] = -s * k_4;
            matrixD[6, 2] = c * k_4;
            matrixD[6, 3] = 2 * k_3;
            matrixD[6, 4] = s * k_4;
            matrixD[6, 5] = -c * k_4;
            matrixD[6, 6] = 4 * k_3;


            return matrixD;

        }
        private string MakeDisplayable_1D_j_0_INT(int[] Array1D)
        {
            // ----- Prepare a multi-line string that shows the contents
            //       of a matrix, a 1D array of integer.
            string str = new string("xx"[0], 5);
            //Debug.Print(str)
            int row = 0;
            int col = 0;
            string s1 = null;
            string s2 = null;
            string s3 = null;


            s1 = "";
            s2 = "";
            s3 = "Inetger array print as below: ";
            for (col = 0; col <= Array1D.GetUpperBound(0); col++)
            {

                s2 = s2 + col.ToString() + new string(' ', col.ToString().Length);
                s1 = s1 + (Array1D[col].ToString()) + " ";

            }

            s1 = s3 + Environment.NewLine + s2 + Environment.NewLine + s1;
            return s1;

        }
        private int FindRowNUm_i0_1D_INT(int[] ptx, int kk)
        {

            int i = 0;
            int j = 0;
            int k = 0;

            int p = ptx.Length - 1;
            int[] ary = new int[p + 1];
            int vvv = -1;
            string[] pt = new string[p + 1];
            //'copy(Array) as string
            for (j = 0; j <= p; j++)
            {
                pt[j] = ptx[j].ToString();
            }
            //Array.Copy(source, target, target.Length)
            //Array.Copy(ptx, pt, ptx.Length)



            for (j = 0; j <= p; j++)
            {

                if (pt[j] == kk.ToString())
                {
                    vvv = (int)j;
                }

            }

            //For j = 0 To p
            //    ary(j) = CInt(pt(j))
            //Next

            return vvv;

        }

        private int[] sort_Integer_1D(int[] ptx)
        {
            string temp1 = null;
            int i = 0;
            int j = 0;
            int k = 0;

            int p = ptx.Length - 1;
            int[] ary = new int[p + 1];
            string[] pt = new string[p + 1];
            //'copy(Array) as string
            for (j = 0; j <= p; j++)
            {
                pt[j] = ptx[j].ToString();
            }
            //Array.Copy(source, target, target.Length)
            //Array.Copy(ptx, pt, ptx.Length)



            for (j = 0; j <= p; j++)
            {
                for (k = 0; k < p; k++)
                {
                    if (int.Parse(pt[k]) > int.Parse(pt[k + 1]))
                    {
                        temp1 = pt[k];
                        pt[k] = pt[k + 1];
                        pt[k + 1] = temp1;
                        // Debug.Print("hii " & pt(k))
                        //Debug.Print("hii ")
                    }
                }
            }

            for (j = 0; j <= p; j++)
            {
                ary[j] = int.Parse(pt[j]);
            }
            //ReDim Preserve pt(p - 1)
            return ary;


        }
        //kkkkkkkkkkkkkkkkkkkkkkkkkkkkk
        private int FindMAX_Integer_1D(int[] ptx)
        {
            string temp1 = null;
            int i = 0;
            int j = 0;
            int k = 0;

            int p = ptx.Length - 1;
            int[] ary = new int[p + 1];
            string[] pt = new string[p + 1];
            //'copy(Array) as string
            for (j = 0; j <= p; j++)
            {
                pt[j] = ptx[j].ToString();
            }
            //Array.Copy(source, target, target.Length)
            //Array.Copy(ptx, pt, ptx.Length)



            for (j = 0; j <= p; j++)
            {
                for (k = 0; k < p; k++)
                {
                    if (int.Parse(pt[k]) > int.Parse(pt[k + 1]))
                    {
                        temp1 = pt[k];
                        pt[k] = pt[k + 1];
                        pt[k + 1] = temp1;
                        // Debug.Print("hii " & pt(k))
                        //Debug.Print("hii ")
                    }
                }
            }

            for (j = 0; j <= p; j++)
            {
                ary[j] = int.Parse(pt[j]);
            }
            //ReDim Preserve pt(p - 1)
            return ary[p];



        }
        //kkkkkkkkkkkkkkkkkkkkkkkkkkkkk
        private string MakeDisplayable_1D_0_Double(double[] Array1D)
        {
            // ----- Prepare a multi-line string that shows the contents
            //       of a matrix, a 1D array of double.
            string str = new string("xx"[0], 5);
            //Debug.Print(str)
            int row = 0;
            int col = 0;
            string s1 = null;
            string s2 = null;
            string s3 = null;


            s1 = "";
            s2 = "";
            s3 = "Double array print as below: ";
            for (col = 0; col <= Array1D.GetUpperBound(0); col++)
            {

                s2 = s2 + col.ToString() + new string(' ', col.ToString().Length);
                s1 = s1 + (Array1D[col].ToString()) + " ";

            }

            s1 = s3 + Environment.NewLine + s2 + Environment.NewLine + s1;
            return s1;

        }


        //kkkkkkkkkkkkkkkkkkkkkkkkkkkkkk

        private string MakeDisplayable_2D2D_ij_1_1(double[,] Array2D)
        {
            int row = 0;
            int col = 0;
            string s1 = null;
            string s2 = null;
            string s3 = null;

            s1 = "";

            for (row = 1; row <= Array2D.GetUpperBound(0); row++)
            {
                for (col = 1; col <= Array2D.GetUpperBound(1); col++)
                {
                    s1 = s1 + (Array2D[row, col].ToString()) + " ";
                }
                s1 = s1 + Environment.NewLine;
            }

            return s1;

        }

        private string MakeDisplayable_2D3D(int k, double[, ,] Array3D)
        {
            // ----- Prepare a multi-line string that shows the contents
            //       of a matrix, a 3D array.

            int row = 0;
            int col = 0;
            int page = 0;
            string s1 = "";




            //for (page = 1; page <= Array3D.GetUpperBound(0); page++)
            for (page = k; page <= k; page++)
            {
                //s1 = s1 + "   " + page.ToString() + " array print as below: " + Environment.NewLine;
                s1 = s1 + "   " + Environment.NewLine;

                for (row = 1; row <= Array3D.GetUpperBound(1); row++)
                {

                    for (col = 1; col <= Array3D.GetUpperBound(2); col++)
                    {
                        s1 = s1 + (Array3D[page, row, col].ToString()) + " ";

                    }
                    s1 = s1 + Environment.NewLine;
                }

                s1 = s1 + "   " + Environment.NewLine;

            }

            return s1;

        }
        private double[,] mult_2D_2D_Mat(double[,] matrixA, double[,] matrixB)
        {
            double sum = 0;
            double[,] matrixC = new double[matrixB.GetUpperBound(0) + 1, matrixB.GetUpperBound(1) + 1];
            sum = 0;
            int x = 0;
            int y = 0;
            int z = 0;
            sum = 0;

            for (x = 0; x <= matrixB.GetUpperBound(0); x++)
            {

                for (y = 0; y <= matrixB.GetUpperBound(0); y++)
                {
                    for (z = 0; z <= matrixB.GetUpperBound(0); z++)
                    {
                        sum = sum + (matrixA[x, z] * matrixB[z, y]);
                    }
                    matrixC[x, y] = sum;
                    sum = 0;
                }
            }
            return matrixC;
        }
        //vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
        private double[] Multiply_Matrix_2DBY1D_0_0(double[,] matrixA, double[] ZZ)
        {


            //    11   | 12      
            //    21   | 22    |X | 1, 2|
            //    31   | 32       
            //    41   | 42     

            //'Debug.Print(UBound(matrixA, 1).ToString + " =i" + UBound(matrixA, 2).ToString + "=j  " + UBound(ZZ).ToString + "=L  ")
            double sum = 0;

            double[] YZ = new double[matrixA.GetUpperBound(0) + 1];
            int p = matrixA.GetUpperBound(1);
            if (matrixA.GetUpperBound(1) != ZZ.GetUpperBound(0))
            {
                MessageBox.Show("Number of rows of 1st array is not equal to Number of rows of 2nd array while Mutiplying");
            }
            try
            {
                if (matrixA.GetUpperBound(1) == ZZ.GetUpperBound(0))
                {

                    sum = 0;



                    int i = 0;
                    int j = 0;

                    for (i = 0; i <= matrixA.GetUpperBound(0); i++)
                    {
                        sum = 0;
                        for (j = 0; j <= ZZ.GetUpperBound(0); j++) //UBound(ZZ) '
                        {
                            //sum = sum + (matrixA[i, j] * ZZ[j]);

                            YZ[i] = YZ[i] + (matrixA[i, j] * ZZ[j]);


                        }



                    }
                }


            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.ToString());
                Debug.Print("Multiply_Matrix_2DBY1D_0_0 ex " + ex.ToString());

            }

            return YZ;
        }////end Multiply_Matrix_2DBY1D_0_0

        //kkkkkkkkkkkkkkkkkkkkkkkkkkkk
        private double[] Union_Matrix_1D_1D_Double(double[] First_Mat, double[] Second_Mat)
        {


            double[] YZ = null;
            int p = First_Mat.Length - 1;
            int g = Second_Mat.Length - 1;
            //ReDim Preserve YZ(p)
            //4 5 6 7 8 9 0 0 0 0 0
            try
            {


                int i = 0;
                int j = 0;
                Array.Resize(ref YZ, p + g + 2); // i become max val of i+1 ?????????????????????????????????????


                for (i = 0; i <= p; i++)
                {


                    YZ[i] = First_Mat[i];
                    //Debug.Print("SUMlast1.................. i " + i.ToString + "    " + YZ(i).ToString)
                }



                //ReDim Preserve YZ(p + g + 1)


                for (i = p + 1; i <= p + g + 1; i++)
                {
                    YZ[i] = Second_Mat[i - p - 1];
                    //Debug.Print("SUMlast.................. i " + i.ToString + "    " + YZ(i).ToString)
                }

                Array.Resize(ref YZ, p + g + 2); //YZ(p + g + 1)=0 ??????????????????
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.ToString());

                return null;
            }

            //i=6,j=6 ????????????????????????????????????????????????????????????????????
            return YZ;

        }
        //kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
        private int[] Union_Matrix_1D_1D_Intger(int[] First_Mat, int[] Second_Mat)
        {


            int[] YZ = null;
            int p = First_Mat.Length - 1;
            int g = Second_Mat.Length - 1;
            //ReDim Preserve YZ(p)
            //4 5 6 7 8 9 0 0 0 0 0
            try
            {


                int i = 0;
                int j = 0;
                Array.Resize(ref YZ, p + g + 2); // i become max val of i+1 ?????????????????????????????????????


                for (i = 0; i <= p; i++)
                {


                    YZ[i] = First_Mat[i];
                    //Debug.Print("SUMlast1.................. i " + i.ToString + "    " + YZ(i).ToString)
                }




                for (i = p + 1; i <= p + g + 1; i++)
                {
                    YZ[i] = Second_Mat[i - p - 1];
                    //Debug.Print("SUMlast.................. i " + i.ToString + "    " + YZ(i).ToString)
                }

                Array.Resize(ref YZ, p + g + 2);
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.ToString());

                return null;
            }

            //i=6,j=6 ????????????????????????????????????????????????????????????????????
            return YZ;

        }
        //hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh


        ///////////////////////////////////////////////////////////////////////////////
        private string MakeDisplayable_2D(double[,] sourceMatrix)
        {
            // ----- Prepare a multi-line string that shows the contents
            //       of a matrix, a 2D array.
            int rows = 0;
            int cols = 0;
            int eachRow = 0;
            int eachCol = 0;
            System.Text.StringBuilder result = new System.Text.StringBuilder();

            // ----- Process all rows of the matrix, generating one
            //       output line per row.
            rows = sourceMatrix.GetUpperBound(0) + 1;
            cols = sourceMatrix.GetUpperBound(1) + 1;
            for (eachRow = 0; eachRow < rows; eachRow++)
            {
                // ----- Process each column of the matrix on a single
                //       row, separating values by commas.
                if (eachRow > 0)
                {
                    result.AppendLine();
                }
                for (eachCol = 0; eachCol < cols; eachCol++)
                {
                    // ----- Add a single matrix element to the output.
                    if (eachCol > 0) //(",  ")
                    {
                        result.Append("  ");
                    }
                    result.Append(sourceMatrix[eachRow, eachCol].ToString());
                }
            }
            result.AppendLine();
            // ----- Finished.
            return result.ToString();
        }
        private string MakeDisplay_INT(int[,] sourceMatrix)
        {
            // ----- Prepare a multi-line string that shows the contents
            //       of a matrix, a 2D array.
            int rows = 0;
            int cols = 0;
            int eachRow = 0;
            int eachCol = 0;
            System.Text.StringBuilder result = new System.Text.StringBuilder();

            // ----- Process all rows of the matrix, generating one
            //       output line per row.
            rows = sourceMatrix.GetUpperBound(0) + 1;
            cols = sourceMatrix.GetUpperBound(1) + 1;
            for (eachRow = 0; eachRow < rows; eachRow++)
            {
                // ----- Process each column of the matrix on a single
                //       row, separating values by commas.
                if (eachRow > 0)
                {
                    result.AppendLine();
                }
                for (eachCol = 0; eachCol < cols; eachCol++)
                {
                    // ----- Add a single matrix element to the output.
                    if (eachCol > 0) //(",  ")
                    {
                        result.Append("  ");
                    }
                    result.Append(sourceMatrix[eachRow, eachCol].ToString());
                }
            }
            result.AppendLine();
            // ----- Finished.
            return result.ToString();
        }
        private string MakeDisplay_String(string[] sourceMatrix)
        {
            // ----- Prepare a multi-line string that shows the contents
            //       of a matrix, a 2D array.
            int rows = 0;
            int cols = 0;
            int eachRow = 0;
            int eachCol = 0;
            System.Text.StringBuilder result = new System.Text.StringBuilder();

            // ----- Process all rows of the matrix, generating one
            //       output line per row.
            rows = sourceMatrix.GetUpperBound(0) + 1;
            //cols = UBound(sourceMatrix, 2) + 1
            for (eachRow = 0; eachRow < rows; eachRow++)
            {
                // ----- Process each column of the matrix on a single
                //       row, separating values by commas.
                if (eachRow > 0)
                {
                    result.AppendLine();
                }
                //For eachCol = 0 To cols - 1
                // ----- Add a single matrix element to the output.
                //If (eachCol > 0) Then result.Append("  ") '(",  ")
                result.Append(sourceMatrix[eachRow].ToString());
                //Next eachCol
            }
            result.AppendLine();
            // ----- Finished.
            return result.ToString();
        }
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        private double[,] Mat_Transpose2nd_Tranform1st(double c, double s)
        {
            string strArray = "";
            //Create three matrices.     
            double[,] matrixB = new double[6, 6];
            double[,] matrixC = new double[6, 6];
            double[,] matrixD = new double[6, 6];
            double[,] matrixE = new double[6, 6];
            double[,] matrixF = new double[6, 6];

            int intI = 0;
            int intJ = 0;

            double[,] matrixA =
			{
				{c, s, 0, 0, 0, 0},
				{-c, s, 0, 0, 0, 0},
				{0, 0, 1, 0, 0, 0},
				{0, 0, 0, c, s, 0},
				{0, 0, 0, -c, s, 0},
				{0, 0, 0, 0, 0, 1}
			};



            //'''''''''''''''''''''''''''''''''''''''''''''''''
            for (intI = 0; intI <= 5; intI++)
            {
                for (intJ = 0; intJ <= 5; intJ++)
                {
                    //transponse of matrix 

                    matrixB[intI, intJ] = matrixA[intJ, intI];
                }
            }

            return matrixB;

        }
        //kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
        private double[,] Trasform_mat_0_3_3(double c, double s)
        {


            double[,] matrixA =
			{
				{c, s, 0},
				{-c, s, 0},
				{0, 0, 1}
			};



            return matrixA;

        }
        //kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk

        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        private double[,] Trasform_mat_GLOBAL_to_LOCAL_ij_00(double c, double s)
        {




            double[,] matrixA =
			{
				{c, -s, 0, 0, 0, 0},
				{s, c, 0, 0, 0, 0},
				{0, 0, 1, 0, 0, 0},
				{0, 0, 0, c, -s, 0},
				{0, 0, 0, s, c, 0},
				{0, 0, 0, 0, 0, 1}
			};





            return matrixA;

        }

        private double[,] Trasform_mat_G_TO_L_ij_00(double c, double s)
        {




            double[,] matrixA =
			{
				{c, s, 0, 0, 0, 0},
				{-s, c, 0, 0, 0, 0},
				{0, 0, 1, 0, 0, 0},
				{0, 0, 0, c, s, 0},
				{0, 0, 0, -s, c, 0},
				{0, 0, 0, 0, 0, 1}
			};





            return matrixA;

        }

        private double[,] Trasform_mat_TT_L_TO_G_ij_00(double c, double s)
        {


            double[,] TT1 = new double[,]
			{
				{c, -s, 0, 0, 0, 0},
				{s, c, 0, 0, 0, 0},
				{0, 0, 1, 0, 0, 0},
				{0, 0, 0, c, -s, 0},
				{0, 0, 0, s, c, 0},
				{0, 0, 0, 0, 0, 1}
			};





            return TT1;

        }
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        private double[,] GlobalMatrix(double c, double s, double A, double EY, double I, double Ln)
        {
            string strArray = "";
            //Create three matrices.     
            double[,] matrixB = new double[6, 6];
            double[,] matrixC = new double[6, 6];
            double[,] matrixD = new double[6, 6];
            double[,] matrixE = new double[6, 6];
            double[,] matrixF = new double[6, 6];

            int intI = 0;
            int intJ = 0;

            double[,] matrixA =
			{
				{c, -s, 0, 0, 0, 0},
				{s, c, 0, 0, 0, 0},
				{0, 0, 1, 0, 0, 0},
				{0, 0, 0, c, -s, 0},
				{0, 0, 0, s, c, 0},
				{0, 0, 0, 0, 0, 1}
			};






            //'''''''''''''''''''''''''''''''''''''''''''''''''
            //TT
            for (intI = 0; intI <= 5; intI++)
            {
                for (intJ = 0; intJ <= 5; intJ++)
                {
                    //transponse of matrix 

                    matrixB[intI, intJ] = matrixA[intJ, intI];
                }
            }






            //define member stifness matrix.......................



            matrixD[0, 0] = A * EY / Ln; //row, col
            matrixD[0, 1] = 0;
            matrixD[0, 2] = 0;
            matrixD[0, 3] = -A * EY / Ln;
            matrixD[0, 4] = 0;
            matrixD[0, 5] = 0;


            matrixD[1, 0] = 0; //row, col
            matrixD[1, 1] = 12 * EY * I / Math.Pow(Ln, 3);
            matrixD[1, 2] = 6 * EY * I / Math.Pow(Ln, 2);
            matrixD[1, 3] = 0;
            matrixD[1, 4] = -12 * EY * I / Math.Pow(Ln, 3);
            matrixD[1, 5] = 6 * EY * I / Math.Pow(Ln, 2);



            matrixD[2, 0] = 0; //row, col
            matrixD[2, 1] = 6 * EY * I / Math.Pow(Ln, 2);
            matrixD[2, 2] = 4 * EY * I / Ln;
            matrixD[2, 3] = 0;
            matrixD[2, 4] = -6 * EY * I / Math.Pow(Ln, 2);
            matrixD[2, 5] = 2 * EY * I / Ln;




            matrixD[3, 0] = -A * EY / Ln; //row, col
            matrixD[3, 1] = 0;
            matrixD[3, 2] = 0;
            matrixD[3, 3] = A * EY / Ln;
            matrixD[3, 4] = 0;
            matrixD[3, 5] = 0;




            matrixD[4, 0] = 0; //row, col
            matrixD[4, 1] = -12 * EY * I / Math.Pow(Ln, 3);
            matrixD[4, 2] = -6 * EY * I / Math.Pow(Ln, 2);
            matrixD[4, 3] = 0;
            matrixD[4, 4] = 12 * EY * I / Math.Pow(Ln, 3);
            matrixD[4, 5] = -6 * EY * I / Math.Pow(Ln, 2);


            matrixD[5, 0] = 0; //row, col
            matrixD[5, 1] = 6 * EY * I / Math.Pow(Ln, 2);
            matrixD[5, 2] = 2 * EY * I / Ln;
            matrixD[5, 3] = 0;
            matrixD[5, 4] = -6 * EY * I / Math.Pow(Ln, 2);
            matrixD[5, 5] = 4 * EY * I / Ln;




            /////////////////////////////////////////////////////////////////////////
            // product1.2------------------------------
            //T*lock
            double sum = 0;
            for (var x = 0; x <= 5; x++)
            {
                for (var y = 0; y <= 5; y++)
                {
                    for (var z = 0; z <= 5; z++)
                    {
                        sum = sum + (matrixA[x, z] * matrixD[z, y]);
                    }
                    matrixE[x, y] = sum;
                    sum = 0;
                }
            }

            // product(12).3-------------------------------

            //(T*locK)*TT
            sum = 0;
            for (var x = 0; x <= 5; x++)
            {
                for (var y = 0; y <= 5; y++)
                {
                    for (var z = 0; z <= 5; z++)
                    {
                        sum = sum + (matrixE[x, z] * matrixB[z, y]);
                    }
                    matrixF[x, y] = sum;
                    sum = 0;
                }
            }
            /////////////////////////////////////////////////////////


            return matrixF;
        }
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        private double[,] K_Local_Element_Matrix_i_0_j_0(double LL, double AA, double II, double EE)
        {

            double[,] matrixD = new double[6, 6];


            //define member stifness matrix.......................


            matrixD[0, 0] = AA * EE / LL; //row, col
            matrixD[0, 1] = 0;
            matrixD[0, 2] = 0;
            matrixD[0, 3] = -AA * EE / LL;
            matrixD[0, 4] = 0;
            matrixD[0, 5] = 0;


            matrixD[1, 0] = 0; //row, col
            matrixD[1, 1] = 12 * EE * II / Math.Pow(LL, 3);
            matrixD[1, 2] = 6 * EE * II / Math.Pow(LL, 2);
            matrixD[1, 3] = 0;
            matrixD[1, 4] = -12 * EE * II / Math.Pow(LL, 3);
            matrixD[1, 5] = 6 * EE * II / Math.Pow(LL, 2);



            matrixD[2, 0] = 0; //row, col
            matrixD[2, 1] = 6 * EE * II / Math.Pow(LL, 2);
            matrixD[2, 2] = 4 * EE * II / LL;
            matrixD[2, 3] = 0;
            matrixD[2, 4] = -6 * EE * II / Math.Pow(LL, 2);
            matrixD[2, 5] = 2 * EE * II / LL;




            matrixD[3, 0] = -AA * EE / LL; //row, col
            matrixD[3, 1] = 0;
            matrixD[3, 2] = 0;
            matrixD[3, 3] = AA * EE / LL;
            matrixD[3, 4] = 0;
            matrixD[3, 5] = 0;



            matrixD[4, 0] = 0; //row, col
            matrixD[4, 1] = -12 * EE * II / Math.Pow(LL, 3);
            matrixD[4, 2] = -6 * EE * II / Math.Pow(LL, 2);
            matrixD[4, 3] = 0;
            matrixD[4, 4] = 12 * EE * II / Math.Pow(LL, 3);
            matrixD[4, 5] = -6 * EE * II / Math.Pow(LL, 2);


            matrixD[5, 0] = 0; //row, col
            matrixD[5, 1] = 6 * EE * II / Math.Pow(LL, 2);
            matrixD[5, 2] = 2 * EE * II / LL;
            matrixD[5, 3] = 0;
            matrixD[5, 4] = -6 * EE * II / Math.Pow(LL, 2);
            matrixD[5, 5] = 4 * EE * II / LL;

            return matrixD;
        }

        private double[,] Inverse_Mat(double[,] a)
        {
            //Uses Gauss elimination method in order to calculate the inverse matrix [A]-1
            //Method: Puts matrix [A] at the left and the singular matrix [I] at the right:
            //[ a11 a12 a13 | 1 0 0 ]
            //[ a21 a22 a23 | 0 1 0 ]
            //[ a31 a32 a33 | 0 0 1 ]
            //Then using line operations, we try to build the singular matrix [I] at the left.
            //After we have finished, the inverse matrix [A]-1 (bij) will be at the right:
            //[ 1 0 0 | b11 b12 b13 ]
            //[ 0 1 0 | b21 b22 b23 ]
            //[ 0 0 1 | b31 b32 b33 ]

            //On Error GoTo errhandler 'In case the inverse cannot be found (Determinant = 0)
            int dm = a.GetUpperBound(0) + 1;
            double[,] Matrix_A = new double[(dm * 2) + 1, (dm * 2) + 1];
            double[,] Operations_Matrix = new double[(dm * 2) + 1, (dm * 2) + 1];
            double[,] Inverse_Matrix = new double[dm + 1, dm + 1];
            double[,] Inverse_MatrixY = new double[dm, dm];

            double temporary_1 = 0;
            double elem1 = 0;
            double multiplier_1 = 0;

            int Max_Index = 0;
            int m = 0;
            int n = 0;
            int k = 0;
            int line_1 = 0;
            int Double_Size_Max_Index = 0;
            Max_Index = a.GetUpperBound(0) + 1;



            for (n = 0; n < Max_Index; n++)
            {
                for (m = 0; m < Max_Index; m++)
                {
                    Matrix_A[m + 1, n + 1] = a[m, n];
                }
            }
            //Assign values from matrix [A] at the left
            for (n = 1; n <= Max_Index; n++)
            {
                for (m = 1; m <= Max_Index; m++)
                {
                    Operations_Matrix[m, n] = Matrix_A[m, n];
                }
            }

            //Assign values from singular matrix [I] at the right
            for (n = 1; n <= Max_Index; n++)
            {
                for (m = 1; m <= Max_Index; m++)
                {
                    if (n == m)
                    {
                        Operations_Matrix[m, n + Max_Index] = 1;
                    }
                    else
                    {
                        Operations_Matrix[m, n + Max_Index] = 0;
                    }
                }
            }

            //Build the Singular matrix [I] at the left
            for (k = 1; k <= Max_Index; k++)
            {
                //Bring a non-zero element first by changes lines if necessary
                if (Operations_Matrix[k, k] == 0)
                {
                    for (n = k; n <= Max_Index; n++)
                    {
                        if (Operations_Matrix[n, k] != 0) //Finds line_1 with non-zero element
                        {
                            line_1 = n;
                            break;
                        }
                    }
                    //Change line k with line_1

                    for (m = k; m <= Max_Index * 2; m++)
                    {
                        temporary_1 = Operations_Matrix[k, m];
                        Operations_Matrix[k, m] = Operations_Matrix[line_1, m];
                        Operations_Matrix[line_1, m] = temporary_1;
                    }
                }

                elem1 = Operations_Matrix[k, k];

                for (n = k; n <= 2 * Max_Index; n++)
                {
                    Operations_Matrix[k, n] = Operations_Matrix[k, n] / elem1;
                    if (elem1 == 0)
                    {
                        MessageBox.Show("Not invertible matrix, determinant zero !");
                        System.Environment.Exit(1); //Function

                    }

                }

                //For other lines, make a zero element by using:
                //Ai1=Aij-A11*(Aij/A11)
                //and change all the line using the same formula for other elements

                Double_Size_Max_Index = 2 * Max_Index;
                for (n = 1; n <= Max_Index; n++)
                {
                    if (n == k && n == Double_Size_Max_Index) //Finished
                    {
                        break;
                    }
                    if (n == k && n < Double_Size_Max_Index) //Do not change that element (already equals to 1), go for next one
                    {
                        n = n + 1;
                    }
                    if (Operations_Matrix[n, k] != 0) //if it is zero, stays as it is
                    {
                        multiplier_1 = Operations_Matrix[n, k] / Operations_Matrix[k, k];

                        for (m = k; m <= 2 * Max_Index; m++)
                        {
                            Operations_Matrix[n, m] = Operations_Matrix[n, m] - Operations_Matrix[k, m] * multiplier_1;
                        }
                    }
                }
            }

            //Assign the right part to the Inverse_Matrix
            for (n = 1; n <= Max_Index; n++)
            {
                for (k = 1; k <= Max_Index; k++)
                {
                    Inverse_Matrix[n, k] = Operations_Matrix[n, Max_Index + k];
                }
            }

            for (n = 1; n <= Max_Index; n++)
            {
                for (k = 1; k <= Max_Index; k++)
                {
                    Inverse_MatrixY[n - 1, k - 1] = Inverse_Matrix[n, k];
                }
            }
            //Debug.Print("Inverse matrix :  " & vbNewLine)

            //Debug.Print(MakeDisplayable_2D2D_ij_0_0(Inverse_MatrixY))



            return Inverse_MatrixY;


        }
        //////////////////////////////////////////////////////////////////////////////
        //===================================================================================




        private void button1_Click(object sender, EventArgs e)
        {
            Debug.Print("CODE STARTS HERE  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");

            //         ^ 
            //         |
            //         |P=-800kN                                                                                                                                                                
            //     2|--|--2-----|3
            //      |           |
            //      1           3
            //      |           |
            //      |1          |4(Spring)


            //=======================================================
            //------------------------------
            //         P
            //         ^ 
            //         | 
            //I--------|-----------J J>I
            //1 I=2-----a--|(+Py)------b--------------------->J=3
            //L=6
            //a=1m,b=5m, P=-800kN 
            //fixed end reaction force at 2 and 3 =740.740740740741 kN, 59.2592592592592 kN, 
            //fixed end Mz force at 2 and 3 =555.555555555556 kN-m, -111.111111111111 kN-m, 
            //2 I-----UDL +wy----------------------------->J

            //3 I-----a--|(Px)------b--------------------->J

            //3 I-----a--|(M+)------b--------------------->J
            //------------------------------
            //meter, force KN
            int Node_Num_SPRING_Supp_Ydir = 4;
            double SPRING_Constant = 900;//kN/m
            totalnode = 4;
            totalElement = 3;


            int[] DOF_Free = { 3, 4, 5, 6, 7, 8, 9, 11 };
            Debug.Print("DOF_Free  :" + MakeDisplayable_1D_j_0_INT(DOF_Free));

            int[] DOF_Restrained = { 1, 2, 10, 12 }; // DOF X, Y, Z==3 each node

            Debug.Print("DOF_Restrained  :" + MakeDisplayable_1D_j_0_INT(DOF_Restrained));

            int[] ALL_DOF_Free_Restrained = Union_Matrix_1D_1D_Intger(DOF_Free, DOF_Restrained); //= {3, 4, 5, 6, 7, 8, 9, 11,1, 2, 10, 12 }

            Debug.Print("ALL_DOF_Free_Restrained  :" + MakeDisplayable_1D_j_0_INT(ALL_DOF_Free_Restrained));

            //...............


            //Dim Deflection_Support(totalnode, Support * 3) As Integer
            double[,] KG_Structure = new double[(totalnode * 3) + 1, (totalnode * 3) + 1];
            double[,] KElementG = new double[(totalnode * 3) + 1, (totalnode * 3) + 1];


            double[] Ln = new double[totalElement + 1];
            double[] A = new double[totalElement + 1];
            double[] Iz = new double[totalElement + 1];
            double[] Ey = new double[totalElement + 1];
            double[] cosx = new double[totalElement + 1];
            double[] sinx = new double[totalElement + 1];



            double[,] matrixB = null;
            double[] YL = new double[DOF_Free.Length];
            double[] XX = null;
            double[,] FIXED_END_MEMBER_LOCAL_REACTION = new double[totalElement + 1, totalnode * 3 + 1];
            Array.Clear(FIXED_END_MEMBER_LOCAL_REACTION, 0, FIXED_END_MEMBER_LOCAL_REACTION.Length);// set all value to 0
            FIXED_END_MEMBER_LOCAL_REACTION[2, 4] = 0;//4>> local x
            FIXED_END_MEMBER_LOCAL_REACTION[2, 5] = 740.740740740741;//5>local y 
            FIXED_END_MEMBER_LOCAL_REACTION[2, 6] = 555.555555555556;//6>local Mz 

            FIXED_END_MEMBER_LOCAL_REACTION[2, 7] = 0;//4>> local x
            FIXED_END_MEMBER_LOCAL_REACTION[2, 8] = 59.2592592592592;//8>local y 
            FIXED_END_MEMBER_LOCAL_REACTION[2, 9] = -111.111111111111;//9>local Mz 
            //Load Vector along Free_DOF

            //GX dir nodal load - fixed end reaction force from member load (nodal load +Equivalent nodal aplied FIXED END FORCE)
            //GY dir nodal load - fixed end reaction force from member load (nodal load +Equivalent nodal aplied FIXED END FORCE)
            //GZ dir nodal load i.e. XY pane - fixed end reaction moment from member load (nodal load +Equivalent nodal aplied FIXED END FORCE)


            // nodeNo 1 >>degree of freedom(DOF)1,2,3 (X,Y, Mz)>>nodeNo*3-2,nodeNo*3-1,nodeNo*3
            // nodeNo 3 >>degree of freedom(DOF)7,8,9 (X,Y, Mz)>>nodeNo*3-2,nodeNo*3-1,nodeNo*3
            YL[0] = 0;//3>local Mz 
            YL[1] = 0;//4>> local x 
            YL[2] = -740.740740740741;//5>local y 
            YL[3] = -555.555555555556;//6>local Mz 
            YL[4] = 0;//7 >> local x 
            YL[5] = -59.2592592592592;//8>local y 
            YL[6] = 111.111111111111;//9>local Mz 
            YL[7] = 0;//11>> local y 
            Coordinate = new double[totalnode + 1, 3]; //The index of the 2nd dimension, 1=X Coordinate, 2=Y Coordinate

            Coordinate[1, 1] = 0;// X 
            Coordinate[1, 2] = 0;// Y
            Coordinate[2, 1] = 0;
            Coordinate[2, 2] = 6;
            Coordinate[3, 1] = 6;
            Coordinate[3, 2] = 6;
            Coordinate[4, 1] = 6;
            Coordinate[4, 2] = 0;
            //=========================
            Connectivity = new int[totalElement + 1, 3]; //The index of the 2nd dimension, 1=Start Node I, 2=End node J of IJ Member

            Connectivity[1, 1] = 1;
            Connectivity[1, 2] = 2;
            Connectivity[2, 1] = 2;
            Connectivity[2, 2] = 3;
            Connectivity[3, 1] = 3;
            Connectivity[3, 2] = 4;

            int i = 0;
            int j = 0;
            int k = 0;

            for (i = 1; i <= totalElement; i++)
            {

                Ln[i] = Math.Sqrt((Math.Pow(Coordinate[Connectivity[i, 2], 1] - Coordinate[Connectivity[i, 1], 1], 2) + Math.Pow(Coordinate[Connectivity[i, 2], 2] - Coordinate[Connectivity[i, 1], 2], 2)));
                cosx[i] = ((Coordinate[Connectivity[i, 2], 1] - Coordinate[Connectivity[i, 1], 1])) / Ln[i];
                sinx[i] = ((Coordinate[Connectivity[i, 2], 2] - Coordinate[Connectivity[i, 1], 2])) / Ln[i];
                Debug.Print(i + "cosx= " + cosx[i] + "sinx= " + sinx[i] + " Ln(i)= " + Ln[i] + "Coordinate  (X1,Y1), (X2,Y2) = " + Coordinate[Connectivity[i, 1], 1] + "," + Coordinate[Connectivity[i, 1], 2] + "," + Coordinate[Connectivity[i, 2], 1] + "," + Coordinate[Connectivity[i, 2], 2]);
                A[i] = 0.0225;//0.0002;
                Iz[i] = 4.21875E-05;// 0.049;
                Ey[i] = 200000000;// 70000000;


            }

            ////##################################################################################################
            ////##################################################################################################

            double[,] Member_KG = null;
            double[, ,] KGM_ALL2 = new double[totalElement + 1, totalnode * 3 + 1, totalnode * 3 + 1];
            //double[, ,] KGM_ALL2 = new double[100, 100 , 100];

            for (i = 1; i <= totalElement; i++)
            {

                Member_KG = K_G_Element_ij_11(Ln[i], A[i], Iz[i], Ey[i], cosx[i], sinx[i]);
                Debug.Print(i + "= Member#  member global stiffness=" + Environment.NewLine + MakeDisplayable_2D2D_ij_1_1(Member_KG));
                int[] Index_1 = Index_6_From_MemID_1D_i_1(i);
                Debug.Print(i + "= Member# Index_1......................" + MakeDisplayable_1D_j_0_INT(Index_1));

                for (int i1 = 1; i1 <= 6; i1++)
                {

                    for (int j1 = 1; j1 <= 6; j1++)
                    {

                        KGM_ALL2[i, Index_1[i1], Index_1[j1]] = Member_KG[i1, j1];
                    }
                }


                Debug.Print(i + "= Member#  member global stiffness member DOF+all other member DOF where all other member DOF stiffness is zero =" + Environment.NewLine + MakeDisplayable_2D3D(i, KGM_ALL2));
            }


            //adding matrix to get structure stiffness matrix i.e. SSM
            double[,] result = new double[totalnode * 3 + 1, totalnode * 3 + 1];
            for (i = 1; i <= totalElement; i++)
            {

                for (j = 1; j <= totalnode * 3; j++)
                {


                    for (k = 1; k <= totalnode * 3; k++)
                    {
                        result[j, k] = result[j, k] + KGM_ALL2[i, j, k];
                    }
                }
            }
            result[Node_Num_SPRING_Supp_Ydir * 3 - 1, Node_Num_SPRING_Supp_Ydir * 3 - 1] = result[Node_Num_SPRING_Supp_Ydir * 3 - 1, Node_Num_SPRING_Supp_Ydir * 3 - 1] + SPRING_Constant;
            Debug.Print("Global Structure Stiffnes Matrix..SSM = result.................");


            Debug.Print("Global Structure Stiffnes Matrix..SSM = result................." + Environment.NewLine + MakeDisplayable_2D2D_ij_1_1(result));


            //DONE......................................................
            Debug.Print("rearange SSM where index start from 0..........");
            double[,] K_SSM_Rearraned = new double[(totalnode * 3) + 1, (totalnode * 3) + 1];

            for (i = 1; i <= totalnode * 3; i++)
            {

                for (j = 1; j <= totalnode * 3; j++)
                {
                    K_SSM_Rearraned[i, j] = result[ALL_DOF_Free_Restrained[i - 1], ALL_DOF_Free_Restrained[j - 1]];

                }
            }

            Debug.Print("K_SSM_Rearraned Global Structure Stiffnes Matrix:" + Environment.NewLine + MakeDisplayable_2D2D_ij_1_1(K_SSM_Rearraned));
            //End K_SSM_Rearraned................................................................
            //"..rearange ssm DONe.............


            //Kff............................................................................
            Debug.Print("KKK" + DOF_Free.Length.ToString());
            double[,] Kff = new double[DOF_Free.Length, DOF_Free.Length]; //Kff
            for (i = 0; i < DOF_Free.Length; i++)
            {
                for (j = 0; j < DOF_Free.Length; j++)
                {
                    Kff[i, j] = K_SSM_Rearraned[i + 1, j + 1]; //Kff  start from 0, K_SSM_Rearranged start from 1 


                }
            }

            Debug.Print("free dof matrix assembled :" + Environment.NewLine + MakeDisplayable_2D(Kff)); //Kff DOF free

            Debug.Print("START Solution for DOF_Free...................");

            matrixB = Inverse_Mat(Kff);

            Debug.Print("Inverse DOF free matrix assembled: " + Environment.NewLine + MakeDisplayable_2D2D_ij_0_0(matrixB));


            Debug.Print("DOF_Free displacements:");

            XX = Multiply_Matrix_2DBY1D_0_0(matrixB, YL);
            double[] Disp_All_Nodes_Global = new double[(totalnode * 3) + 1];
            double[] Disp_All_Nodes_Global1 = new double[(totalnode * 3) + 1];
            double[] disp_Element = new double[6];
            Array.Clear(Disp_All_Nodes_Global, 0, Disp_All_Nodes_Global.Length);

            for (j = 0; j < XX.Length; j++)
            {
                Debug.Print("XX(" + DOF_Free[j].ToString() + ")" + XX[j]);
                Disp_All_Nodes_Global[DOF_Free[j]] = XX[j];
            }
            for (j = 0; j < XX.Length; j++)
            {
                Debug.Print("DOF wise disp=  " + DOF_Free[j] + "  " + Disp_All_Nodes_Global[DOF_Free[j]]);
            }

            for (j = 0; j <= totalnode * 3; j++)
            {
                Debug.Print("Disp_All_Nodes_Global=  " + j + " , " + +Disp_All_Nodes_Global[j]);
            }


            Debug.Print("START Solution for Ksf LOWER LEFT PART OF K_SSM_Rearraned...................");
            //Debug.Print("KKK" + DOF_Free.Length.ToString)
            double[,] Ksf = new double[DOF_Restrained.Length, DOF_Free.Length]; //DOF_Restrained.Length=6




            for (i = DOF_Free.Length + 1; i <= totalnode * 3; i++)
            {
                for (j = 1; j <= DOF_Free.Length; j++)
                {
                    Ksf[i - DOF_Free.Length - 1, j - 1] = K_SSM_Rearraned[i, j]; ////Ksf start from 0, K_SSM_Rearranged start from 1 
                    //Debug.Print("TEST Ksf LOWER LEFT PART OF K_SSM_Rearraned " + (DOF_Free.Length + 1).ToString + "TO  " + (totalnode * 3).ToString)
                }
            }




            Debug.Print("Ksf Ksf LOWER LEFT PART OF K_SSM_Rearraned Matrix ..........");
            Debug.Print(MakeDisplayable_2D(Ksf));


            //STRAT Reaction Soln.........................
            double[] Reaction;


            Reaction = Multiply_Matrix_2DBY1D_0_0(Ksf, XX);
            for (j = 0; j < Reaction.Length; j++)
            {
                Debug.Print(j.ToString() + "  Reaction Soln.... " + Reaction[j]);
            }
            Debug.Print("Node_Num_SPRING_Supp, DOF_SPRING_Ydir,SPRING_Disp m, SPRING_Force kN  ="
                + Node_Num_SPRING_Supp_Ydir
                + " , " + (Node_Num_SPRING_Supp_Ydir * 3 - 1)
                + ", " + Disp_All_Nodes_Global[Node_Num_SPRING_Supp_Ydir * 3 - 1]
                + ", " + Disp_All_Nodes_Global[Node_Num_SPRING_Supp_Ydir * 3 - 1] * SPRING_Constant);

            //
            //Reaction Soln.... done......................





            //---------------------------------------------
            double[] Disp_j0_global = new double[6];
            double[] Disp_j0_Local;


            //IF MEBER HAS LOAD AND MEBER LOCAL AND GLOBA POSITION NOT SAME............========================

            //step-1.............................................................................................
            //q=[k][d]+q_Fixed
            //[Q joint applied] = [Kff][Dff] + [QFEM as reactio or -QEquivalent nodal aplied FEM] //all all global


            //=>[Kff][D]=[Q joint applied] - [-QFEM as reactio or +QEquivalent nodal aplied FEM] //all all global

            //step2 memberswise..................................................................................

            //[q]1 = [k]1[d]1 + [qF as reaction]1 //all all global

            //[q]1 = [k]1[d]1 - [qF as Equivalent nodal aplied FEM]1 //all all global !!!!!!!!!!!!!!!!!!!!!!!

            ///==================================================================================================================
            double[] Force_Local_Axix = new double[totalnode * 3 + 1];
            Array.Clear(Force_Local_Axix, 0, Force_Local_Axix.Length);// set all value to 0
            for (k = 1; k <= totalElement; k++)
            {
                Disp_j0_global[0] = Disp_All_Nodes_Global[Connectivity[k, 1] * 3 - 2]; //d global FORCE START NODE GX DIR
                Disp_j0_global[1] = Disp_All_Nodes_Global[Connectivity[k, 1] * 3 - 1]; //d global FORCE START NODE GY DIR
                Disp_j0_global[2] = Disp_All_Nodes_Global[Connectivity[k, 1] * 3]; //d global MOMENT START NODE GZ DIR
                Disp_j0_global[3] = Disp_All_Nodes_Global[Connectivity[k, 2] * 3 - 2]; //d global FORCE END NODE GX DIR
                Disp_j0_global[4] = Disp_All_Nodes_Global[Connectivity[k, 2] * 3 - 1]; //d global  FORCE END NODE GY DIR
                Disp_j0_global[5] = Disp_All_Nodes_Global[Connectivity[k, 2] * 3]; //d global  MOMENT END NODE GZ DIR





                // Debug.Print("Disp_j0_global=  " + MakeDisplayable_1D_0_Double(Disp_j0_global));
                //dlocal= T*dglobal
                double[,] matrix1 = null;
                matrix1 = Trasform_mat_G_TO_L_ij_00(cosx[k], sinx[k]);

                Disp_j0_Local = Multiply_Matrix_2DBY1D_0_0(matrix1, Disp_j0_global);
                for (j = 0; j <= 5; j++)
                {
                    // Debug.Print(j + " =j "+   "Disp_j0_Local= "  + Disp_j0_Local[j]);
                }
                double[,] matrix2 = K_Local_Element_Matrix_i_0_j_0(Ln[k], A[k], Iz[k], Ey[k]);

                //matrix22 = K_G_Element_ij_11i_0j_0(Ln[k], A[k], Iz[k], Ey[k], cosx[k], sinx[k]);

                double[] Force_FROM_Local_Displacement = Multiply_Matrix_2DBY1D_0_0(matrix2, Disp_j0_Local);//klocal*dlocal

                Force_Local_Axix[Connectivity[k, 1] * 3 - 2] = Force_FROM_Local_Displacement[0] + FIXED_END_MEMBER_LOCAL_REACTION[k, Connectivity[k, 1] * 3 - 2]; //d global FORCE START NODE GX DIR
                Force_Local_Axix[Connectivity[k, 1] * 3 - 1] = Force_FROM_Local_Displacement[1] + FIXED_END_MEMBER_LOCAL_REACTION[k, Connectivity[k, 1] * 3 - 1]; //d global FORCE START NODE GY DIR
                Force_Local_Axix[Connectivity[k, 1] * 3] = Force_FROM_Local_Displacement[2] + FIXED_END_MEMBER_LOCAL_REACTION[k, Connectivity[k, 1] * 3]; //d global MOMENT START NODE GZ DIR
                Force_Local_Axix[Connectivity[k, 2] * 3 - 2] = Force_FROM_Local_Displacement[3] + FIXED_END_MEMBER_LOCAL_REACTION[k, Connectivity[k, 2] * 3 - 2]; //d global FORCE END NODE GX DIR
                Force_Local_Axix[Connectivity[k, 2] * 3 - 1] = Force_FROM_Local_Displacement[4] + FIXED_END_MEMBER_LOCAL_REACTION[k, Connectivity[k, 2] * 3 - 1]; //d global  FORCE END NODE GY DIR
                Force_Local_Axix[Connectivity[k, 2] * 3] = Force_FROM_Local_Displacement[5] + FIXED_END_MEMBER_LOCAL_REACTION[k, Connectivity[k, 2] * 3]; //d global  MOMENT END NODE GZ DIR


                double[] Force_Local = new double[6];

                Force_Local[0] = Force_Local_Axix[Connectivity[k, 1] * 3 - 2];
                Force_Local[1] = Force_Local_Axix[Connectivity[k, 1] * 3 - 1];
                Force_Local[2] = Force_Local_Axix[Connectivity[k, 1] * 3];
                Force_Local[3] = Force_Local_Axix[Connectivity[k, 2] * 3 - 2];
                Force_Local[4] = Force_Local_Axix[Connectivity[k, 2] * 3 - 1];
                Force_Local[5] = Force_Local_Axix[Connectivity[k, 2] * 3];
                Debug.Print(k + "Force_Local=  " + Environment.NewLine + MakeDisplayable_1D_0_Double(Force_Local));
                ////Debug.Print(k + "Disp_j0_Local=  " + Environment.NewLine + MakeDisplayable_1D_0_Double(Disp_j0_Local));
                double[] Force_Global = null;
                Force_Global = Multiply_Matrix_2DBY1D_0_0(Trasform_mat_TT_L_TO_G_ij_00(cosx[k], sinx[k]), Force_Local);

                Debug.Print("Member Number= " + k + "  Force_Global=  " + Environment.NewLine + MakeDisplayable_1D_0_Double(Force_Global));
                Debug.Print(" SPRING_Force kN  =" + Disp_All_Nodes_Global[Node_Num_SPRING_Supp_Ydir * 3 - 1] * SPRING_Constant);



            }
            //Done..............................
            ////##################################################################################################
            ////##################################################################################################
        }

        private void button2_Click(object sender, EventArgs e)
        {


            Debug.Print("CODE STARTS HERE  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");

            //meter, force KN
            int Node_Num_SPRING_Supp_Ydir = 2;
            double SPRING_Constant = 200;//kN/m
            totalnode = 3;
            totalElement = 2;


            int[] DOF_Free = { 3, 5, 7, 9 };
            Debug.Print("DOF_Free  :" + MakeDisplayable_1D_j_0_INT(DOF_Free));

            int[] DOF_Restrained = { 1, 2, 4, 6, 8 }; // DOF X, Y, Z==3 each node

            Debug.Print("DOF_Restrained  :" + MakeDisplayable_1D_j_0_INT(DOF_Restrained));

            int[] ALL_DOF_Free_Restrained = Union_Matrix_1D_1D_Intger(DOF_Free, DOF_Restrained); //= {4, 5, 6, 7, 8, 9, 1, 2, 3, 10, 11, 12}

            Debug.Print("ALL_DOF_Free_Restrained  :" + MakeDisplayable_1D_j_0_INT(ALL_DOF_Free_Restrained));

            //...............


            //Dim Deflection_Support(totalnode, Support * 3) As Integer
            double[,] KG_Structure = new double[(totalnode * 3) + 1, (totalnode * 3) + 1];
            double[,] KElementG = new double[(totalnode * 3) + 1, (totalnode * 3) + 1];


            double[] Ln = new double[totalElement + 1];
            double[] A = new double[totalElement + 1];
            double[] Iz = new double[totalElement + 1];
            double[] Ey = new double[totalElement + 1];
            double[] cosx = new double[totalElement + 1];
            double[] sinx = new double[totalElement + 1];



            double[,] matrixB = null;
            double[] YL = new double[DOF_Free.Length];
            double[] XX = null;

            //Load Vector along Free_DOF

            //GX dir nodal load - fixed end reaction force from member load (OR +QEquivalent nodal aplied FIXED END FORCE)
            //GY dir nodal load - fixed end reaction force from member load (OR +QEquivalent nodal aplied FIXED END FORCE)
            //GZ dir nodal load i.e. XY pane - fixed end reaction moment from member load (OR +QEquivalent nodal aplied FIXED END FORCE)
            double[,] FIXED_END_MEMBER_GLOBAL_FORCE = new double[totalElement + 1, (totalnode * 3) + 1];
            //Result_Mem_Force_global = new double[totalElement + 1, (totalnode * 3) + 1];
            YL[0] = 0; //3
            YL[1] = -12; //5
            YL[2] = 0; //7
            YL[3] = 0; //9

            Coordinate = new double[totalnode + 1, 3]; //The index of the 2nd dimension, 1=X Coordinate, 2=Y Coordinate

            Coordinate[1, 1] = 0;
            Coordinate[1, 2] = 0;
            Coordinate[2, 1] = 4;
            Coordinate[2, 2] = 0;
            Coordinate[3, 1] = 8;
            Coordinate[3, 2] = 0;
            //=========================
            Connectivity = new int[totalElement + 1, 3]; //The index of the 2nd dimension, 1=Start Node I, 2=End node J of IJ Member

            Connectivity[1, 1] = 1;
            Connectivity[1, 2] = 2;
            Connectivity[2, 1] = 2;
            Connectivity[2, 2] = 3;


            int i = 0;
            int j = 0;
            int k = 0;
            double hz = 0.3;
            for (i = 1; i <= totalElement; i++)
            {

                Ln[i] = Math.Sqrt((Math.Pow(Coordinate[Connectivity[i, 2], 1] - Coordinate[Connectivity[i, 1], 1], 2) + Math.Pow(Coordinate[Connectivity[i, 2], 2] - Coordinate[Connectivity[i, 1], 2], 2)));
                cosx[i] = ((Coordinate[Connectivity[i, 2], 1] - Coordinate[Connectivity[i, 1], 1])) / Ln[i];
                sinx[i] = ((Coordinate[Connectivity[i, 2], 2] - Coordinate[Connectivity[i, 1], 2])) / Ln[i];
                Debug.Print(i + "cosx= " + cosx[i] + "sinx= " + sinx[i] + " Ln(i)= " + Ln[i] + "Coordinate  (X1,Y1), (X2,Y2) = " + Coordinate[Connectivity[i, 1], 1] + "," + Coordinate[Connectivity[i, 1], 2] + "," + Coordinate[Connectivity[i, 2], 1] + "," + Coordinate[Connectivity[i, 2], 2]);
                Iz[i] = 0.0002;
                A[i] = 0.049;
                Ey[i] = 70000000;



            }

            ////##################################################################################################
            ////##################################################################################################

            double[,] Member_KG = null;
            double[, ,] KGM_ALL2 = new double[totalElement + 1, totalnode * 3 + 1, totalnode * 3 + 1];
            //double[, ,] KGM_ALL2 = new double[100, 100 , 100];

            for (i = 1; i <= totalElement; i++)
            {

                Member_KG = K_G_Element_ij_11(Ln[i], A[i], Iz[i], Ey[i], cosx[i], sinx[i]);
                Debug.Print(i + "= Member#  member global stiffness=" + Environment.NewLine + MakeDisplayable_2D2D_ij_1_1(Member_KG));
                int[] Index_1 = Index_6_From_MemID_1D_i_1(i);
                Debug.Print(i + "= Member# Index_1......................" + MakeDisplayable_1D_j_0_INT(Index_1));

                for (int i1 = 1; i1 <= 6; i1++)
                {

                    for (int j1 = 1; j1 <= 6; j1++)
                    {

                        KGM_ALL2[i, Index_1[i1], Index_1[j1]] = Member_KG[i1, j1];
                    }
                }


                Debug.Print(i + "= Member#  member global stiffness member DOF+all other member DOF where all other member DOF stiffness is zero =" + Environment.NewLine + MakeDisplayable_2D3D(i, KGM_ALL2));
            }


            //adding matrix to get structure stiffness matrix i.e. SSM
            double[,] result = new double[totalnode * 3 + 1, totalnode * 3 + 1];
            for (i = 1; i <= totalElement; i++)
            {

                for (j = 1; j <= totalnode * 3; j++)
                {


                    for (k = 1; k <= totalnode * 3; k++)
                    {
                        result[j, k] = result[j, k] + KGM_ALL2[i, j, k];
                    }
                }
            }
            result[Node_Num_SPRING_Supp_Ydir * 3 - 1, Node_Num_SPRING_Supp_Ydir * 3 - 1] = result[Node_Num_SPRING_Supp_Ydir * 3 - 1, Node_Num_SPRING_Supp_Ydir * 3 - 1] + SPRING_Constant;
            Debug.Print("Global Structure Stiffnes Matrix..SSM = result.................");


            Debug.Print("Global Structure Stiffnes Matrix..SSM = result................." + Environment.NewLine + MakeDisplayable_2D2D_ij_1_1(result));


            //DONE......................................................
            Debug.Print("rearange SSM where index start from 0..........");
            double[,] K_SSM_Rearraned = new double[(totalnode * 3) + 1, (totalnode * 3) + 1];

            for (i = 1; i <= totalnode * 3; i++)
            {

                for (j = 1; j <= totalnode * 3; j++)
                {
                    K_SSM_Rearraned[i, j] = result[ALL_DOF_Free_Restrained[i - 1], ALL_DOF_Free_Restrained[j - 1]];

                }
            }

            Debug.Print("K_SSM_Rearraned Global Structure Stiffnes Matrix:" + Environment.NewLine + MakeDisplayable_2D2D_ij_1_1(K_SSM_Rearraned));
            //End K_SSM_Rearraned................................................................
            //"..rearange ssm DONe.............


            //Kff............................................................................
            Debug.Print("KKK" + DOF_Free.Length.ToString());
            double[,] Kff = new double[DOF_Free.Length, DOF_Free.Length]; //Kff
            for (i = 0; i < DOF_Free.Length; i++)
            {
                for (j = 0; j < DOF_Free.Length; j++)
                {
                    Kff[i, j] = K_SSM_Rearraned[i + 1, j + 1]; //Kff  start from 0, K_SSM_Rearranged start from 1 


                }
            }

            Debug.Print("free dof matrix assembled :" + Environment.NewLine + MakeDisplayable_2D(Kff)); //Kff DOF free

            Debug.Print("START Solution for DOF_Free...................");

            matrixB = Inverse_Mat(Kff);

            Debug.Print("Inverse DOF free matrix assembled: " + Environment.NewLine + MakeDisplayable_2D2D_ij_0_0(matrixB));


            Debug.Print("DOF_Free displacements:");

            XX = Multiply_Matrix_2DBY1D_0_0(matrixB, YL);
            double[] Disp_All_Nodes_Global = new double[(totalnode * 3) + 1];
            double[] Disp_All_Nodes_Global1 = new double[(totalnode * 3) + 1];
            double[] disp_Element = new double[6];
            Array.Clear(Disp_All_Nodes_Global, 0, Disp_All_Nodes_Global.Length);

            for (j = 0; j < XX.Length; j++)
            {
                Debug.Print("XX(" + DOF_Free[j].ToString() + ")" + XX[j]);
                Disp_All_Nodes_Global[DOF_Free[j]] = XX[j];
            }
            for (j = 0; j < XX.Length; j++)
            {
                Debug.Print("DOF wise disp=  " + DOF_Free[j] + "  " + Disp_All_Nodes_Global[DOF_Free[j]]);
            }

            for (j = 0; j <= totalnode * 3; j++)
            {
                Debug.Print("Disp_All_Nodes_Global=  " + j + " , " + +Disp_All_Nodes_Global[j]);
            }


            Debug.Print("START Solution for Ksf LOWER LEFT PART OF K_SSM_Rearraned...................");
            //Debug.Print("KKK" + DOF_Free.Length.ToString)
            double[,] Ksf = new double[DOF_Restrained.Length, DOF_Free.Length]; //DOF_Restrained.Length=6




            for (i = DOF_Free.Length + 1; i <= totalnode * 3; i++)
            {
                for (j = 1; j <= DOF_Free.Length; j++)
                {
                    Ksf[i - DOF_Free.Length - 1, j - 1] = K_SSM_Rearraned[i, j]; ////Ksf start from 0, K_SSM_Rearranged start from 1 
                    //Debug.Print("TEST Ksf LOWER LEFT PART OF K_SSM_Rearraned " + (DOF_Free.Length + 1).ToString + "TO  " + (totalnode * 3).ToString)
                }
            }

            //for (i = 0; i <= DOF_Restrained.Length - 1; i++) 
            //{
            //    for (j = 0; j <= DOF_Free.Length - 1; j++) 
            //    {
            //        Ksf[i, j] = K_SSM_Rearraned[DOF_Free.Length - 1 + 1 + i + 1, j + 1]; //Ksf start from 0, K_SSM_Rearranged start from 1 
            //    }
            //}


            Debug.Print("Ksf Ksf LOWER LEFT PART OF K_SSM_Rearraned Matrix ..........");
            Debug.Print(MakeDisplayable_2D(Ksf));


            //STRAT Reaction Soln.........................
            double[] Reaction;


            Reaction = Multiply_Matrix_2DBY1D_0_0(Ksf, XX);
            for (j = 0; j < Reaction.Length; j++)
            {
                Debug.Print(j.ToString() + "  Reaction Soln.... " + Reaction[j]);
            }
            Debug.Print("Node_Num_SPRING_Supp_Ydir, DOF_SPRING_Ydir, SPRING_Force kN  =" + (Node_Num_SPRING_Supp_Ydir * 3 - 1) + ", " + Disp_All_Nodes_Global[Node_Num_SPRING_Supp_Ydir * 3 - 1] + ", " + Disp_All_Nodes_Global[Node_Num_SPRING_Supp_Ydir * 3 - 1] * SPRING_Constant);
            //
            //Reaction Soln.... done......................





            //---------------------------------------------
            double[] Disp_j0_global = new double[6];
            double[] Disp_j0_Local;


            //IF MEBER HAS LOAD AND MEBER LOCAL AND GLOBA POSITION NOT SAME............========================

            //step-1.............................................................................................
            //q=[k][d]+q_Fixed
            //[Q joint applied] = [Kff][Dff] + [QFEM as reactio or -QEquivalent nodal aplied FEM] //all all global


            //=>[Kff][D]=[Q joint applied] - [-QFEM as reactio or +QEquivalent nodal aplied FEM] //all all global

            //step2 memberswise..................................................................................

            //[q]1 = [k]1[d]1 + [qF as reaction]1 //all all global

            //[q]1 = [k]1[d]1 - [qF as Equivalent nodal aplied FEM]1 //all all global !!!!!!!!!!!!!!!!!!!!!!!

            ///==================================================================================================================

            for (k = 1; k <= totalElement; k++)
            {
                Disp_j0_global[0] = Disp_All_Nodes_Global[Connectivity[k, 1] * 3 - 2]; //d global FORCE START NODE GX DIR
                Disp_j0_global[1] = Disp_All_Nodes_Global[Connectivity[k, 1] * 3 - 1]; //d global FORCE START NODE GY DIR
                Disp_j0_global[2] = Disp_All_Nodes_Global[Connectivity[k, 1] * 3]; //d global MOMENT START NODE GZ DIR
                Disp_j0_global[3] = Disp_All_Nodes_Global[Connectivity[k, 2] * 3 - 2]; //d global FORCE END NODE GX DIR
                Disp_j0_global[4] = Disp_All_Nodes_Global[Connectivity[k, 2] * 3 - 1]; //d global  FORCE END NODE GY DIR
                Disp_j0_global[5] = Disp_All_Nodes_Global[Connectivity[k, 2] * 3]; //d global  MOMENT END NODE GZ DIR





                Debug.Print("Disp_j0_global=  " + MakeDisplayable_1D_0_Double(Disp_j0_global));
                //dlocal= T*dgloal
                double[,] matrix1 = null;
                matrix1 = Trasform_mat_GLOBAL_to_LOCAL_ij_00(cosx[k], sinx[k]);

                Disp_j0_Local = Multiply_Matrix_2DBY1D_0_0(matrix1, Disp_j0_global);
                for (j = 0; j <= 5; j++)
                {
                    //Debug.Print(j + " =j "+   "Disp_j0_Local= "  + Disp_j0_Local[j]);
                }

                double[,] matrix2 = null;
                matrix2 = K_G_Element_ij_11i_0j_0(Ln[k], A[k], Iz[k], Ey[k], cosx[k], sinx[k]);
                double[] Force_Global = null;
                //Debug.Print(k + "Disp_j0_Local=  " + Environment.NewLine + MakeDisplayable_1D_0_Double(Disp_j0_Local));
                Force_Global = Multiply_Matrix_2DBY1D_0_0(matrix2, Disp_j0_global);
                Debug.Print("Member Number= " + k + "  Force_Global=  " + Environment.NewLine + MakeDisplayable_1D_0_Double(Force_Global));
                Debug.Print(" SPRING_Force kN  =" + Disp_All_Nodes_Global[Node_Num_SPRING_Supp_Ydir * 3 - 1] * SPRING_Constant);
                double[] Force_Local = Multiply_Matrix_2DBY1D_0_0(Trasform_mat_GLOBAL_to_LOCAL_ij_00(cosx[k], sinx[k]), Force_Global);//ok
                //Debug.Print("Member Number= " + k + "  Force_Local=  " + Environment.NewLine + MakeDisplayable_1D_0_Double(Force_Local));


            }
            //Done..............................
            ////##################################################################################################
            ////##################################################################################################
        }





        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    }
}

