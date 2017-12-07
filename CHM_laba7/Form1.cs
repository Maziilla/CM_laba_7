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
                return -30 * Math.Pow((1 - x), 2) * x * x * x * Math.Log(14) - 30 * x * x * Math.Pow((x - 1), 3) * Math.Log(14) + 1;
            return -30 * Math.Pow((1 - x), 2) * x * x * x * Math.Log(23) - 30 * x * x * Math.Pow((x - 1), 3) * Math.Log(23) + 1;
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
        public void RungeKutt()
        {
            string text="";
            double x = 1, x0 = 0, a = 1, b = 2, xn = 1, RungeEps = 1e-5,z0,y0,boundEps = 1e-4;
            Func<double, double> f = alpha => RungeKutt_solve(f1, f2, x0, a, alpha, 1, RungeEps,ref text, 0) - b;

            double start = -100;
            double end = 100;

            if (f(start) * f(end) > 0)
                throw new ArgumentException("f(start) * f(end) > 0");

            double mid, error;

            int i = 1;

            var start_sign = Math.Sign(f(start));
            strList.Add(string.Format("{0, 2};{1, 10};{2, 10};{3, 10};{4, 15};{5, 15};{6, 15};{7, 10}",
                "i", "left", "right", "alpha", "maxKuttaError", "y(1)", "error", "y(1, alpha)-b"));

            do
            {
                string str = "";
                str += string.Format("{0,2};{0, 10:N7};{0, 10:N7}", i++, start, end);
                //пристрелка
                mid = (start + end) / 2;
                var midSign = Math.Sign(f(mid));
                if (start_sign * midSign > 0)
                    start = mid;
                else
                    end = mid;
                error = Math.Abs(f(mid));
                //
                text = "";
                RungeKutt_solve(f1, f2, x0, a, mid, 1, RungeEps, ref text, 2);
                str += string.Format("{0, 10:N7};{0, 15:N7};{0, 15};{0, 10:N7}", mid, text, error, f(mid));
                strList.Add(str);
                SolveBox.Items.Add(str);
            }
            while (error > boundEps);
            strList.Add("Alpha = " + mid);

            strList.Add("Y(1) = "+RungeKutt_solve(f1, f2, x0, a, mid, x, RungeEps, ref text, 1));

        }
        public double RungeKutt_solve(Func<double, double, double, double> f1,
            Func<double, double, double, double> f2,
            double x0,
            double y0,
            double z0,
            double x,
            double eps,
            ref string text,
            int Out)
        {
            var h = 0.1;
            var xj = x0;
            var yj = y0;
            var zj = z0;
            var maxError = double.MinValue;
            if (x == x0)
                return y0;
            if (Out == 1)
                strList.Add(string.Format("{0, 7};{1, 12};{2, 12};{3, 12};{4, 16}", "h", "x", "y(x)", "y*(x)", "error"));
            do
            {
                string str = "";
                if (h > x - xj)
                    h = x - xj;
                var next = For_RungeKut(f1, f2, xj, yj, zj, h);//Рунге для текущего х и h
                var next_stage = For_RungeKut(f1, f2, xj, yj, zj, h / 2.0);//Рунге для текущего х и h/2
                var _next = For_RungeKut(f1, f2, xj + h / 2.0, next_stage.Item1, next_stage.Item2, h / 2.0);//Рунге для х+h/2 и h/2
                var Error = Math.Abs(next.Item1 - _next.Item1);
                if (Error < eps && maxError == double.MinValue)
                    maxError = Error;
                if (Error < eps && Error <= maxError)
                {
                    if(Out == 1)
                        str += string.Format("|{0, 7}|{1, 12}|", h, xj);
                    xj += h;
                    yj = next.Item1;
                    zj = next.Item2;
                    if (Out == 1)
                    {
                        str += string.Format("{0, 12}|{1, 12}|{2, 16}|", yj, _next.Item1, Error);
                        strList.Add(str);
                    }
                }
                else h /= 2;
            } while (xj < x);
            if (Out == 1)
                strList.Add("Max error = " + maxError);
            else
                if (Out == 2)
                    text = string.Format("{0, 15};", maxError);
            return yj;
        }
        //Для функции
        Tuple<double, double> For_RungeKut(
            Func<double, double, double, double> f1,
            Func<double, double, double, double> f2,
            double xi,
            double yi,
            double zi,
            double h)
        {
            var m1 = h * f1(xi, yi, zi);
            var k1 = h * f2(xi, yi, zi);
            var m2 = h * f1(xi + h / 2, yi + m1 / 2, zi + k1 / 2);
            var k2 = h * f2(xi + h / 2, yi + m1 / 2, zi + k1 / 2);
            var m3 = h * f1(xi + h / 2, yi + m2 / 2, zi + k2 / 2);
            var k3 = h * f2(xi + h / 2, yi + m2 / 2, zi + k2 / 2);
            var m4 = h * f1(xi + h, yi + m3, zi + k3);
            var k4 = h * f2(xi + h, yi + m3, zi + k3);

            return Tuple.Create(
                yi + 1.0 / 6.0 * (m1 + 2 * m2 + 2 * m3 + m4), //y
                zi + 1.0 / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4));//z
        }

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
            RungeKutt();
            SaveFile();
        }

       
    }
}
