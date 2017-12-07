using System;
using System.Windows.Forms;
using System.IO;
using System.Collections.Generic;

namespace SLAU
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
            SolveButton.Enabled = true;
        }
        public const int a = 1, b = 2;
        public  int  N=5; //Размерности
        public double[] solution,RightPart; //Для решения
        public int[] kol;
        const double E = 1E-8; //Точность
        public double chis;
        List<string> strList = new List<string>();

        public delegate double F_(double x, double y);

        //Сохраняем решение в файл
        public void SaveFile()
        {
            SaveFileDialog saveFile = new SaveFileDialog();
            saveFile.Filter = "Текстовый(*.txt)|*.txt";
            saveFile.DefaultExt = "txt";
            saveFile.Title = "Сохранение решения";

            if (saveFile.ShowDialog() == DialogResult.OK)
            {
                FileStream fs = new FileStream(saveFile.FileName, FileMode.Create);
                StreamWriter st = new StreamWriter(fs);
                foreach (string str in strList)
                {
                        st.WriteLine(str);
                }
                st.Close();
            }
            SolveBox.Items.Add("Решение сохранено в " + saveFile.FileName);
            strList.Clear();
        }       
        public double Y_pr(double x)//yp
        {
            if (rb_13.Checked)
                return 1 + x + 10 * Math.Log(14) * x * x * x * Math.Pow((1 - x), 3);
            return 1 + x + 10 * Math.Log(23) * x * x * x * Math.Pow((1 - x), 3);
        }
        public double Y_pr_derivative(double x)//d1yp
        {
            if (rb_13.Checked)
                return -30 * Math.Pow((1 - x), 2) * x * x * x * Math.Log(14) - 30 * x * x * Math.Pow((1 - x), 3) * Math.Log(14) + 1;
            return -30 * Math.Pow((1 - x), 2) * x * x * x * Math.Log(23) - 30 * x * x * Math.Pow((1 - x), 3) * Math.Log(23) + 1;
        }
        public double Y_pr_derivative_2(double x)//d2yp
        {
            if (rb_13.Checked)
                return -60 * x * (5 * x * x * x - 10 * x * x + 6 * x - 1) * Math.Log(14);
            return -60 * x * (5 * x * x * x - 10 * x * x + 6 * x - 1) * Math.Log(23);
        }
        public double A(double x)
        {
            if (rb_13.Checked)
                return 30 * (x + 1);
            return 50 * (x + 1);
        }
        public double B(double x)
        {
            if (rb_13.Checked)
                return x * x + 1;
            return x * x + 1;
        }
        public double C(double x)
        {
            if (rb_13.Checked)
                return x + 1;
            return 2 * x + 1;
        }
        public double F(double x)//F
        {
            return Y_pr_derivative_2(x) + A(x) * Y_pr_derivative(x) - B(x) * Y_pr(x) + C(x) * Math.Sin(Y_pr(x));
        }
        Func<double, double, double, double> f1 = (x, y, z) => z;
        public double f2(double x, double y, double z)
        {
            return -A(x) * z + B(x) * y - C(x) * Math.Sin(y) + F(x);
        }
        public double RungeKutt(F_ func,double h)
        {
            return 0;
        }
        //Для функции
       
       
        public long Factorial(int n)
        {
            long temp = 1;
            for (int i = 1; i <= n; i++)
                temp *= i;
            return temp;
        }
       
        //Решение уравнения
        private void SolveButton_Click(object sender, EventArgs e)
        {
            SaveFile();
        }

       
    }
}
