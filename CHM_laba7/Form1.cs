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
        public const double a = 0.0, b = 1.0;
        public  int  N=5; //Размерности
        public double[] solution; //Для решения
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
            strList.Add(string.Format("{0, 2}|{1, 12}|{2, 11}|{3, 10}|{4, 21}|{5, 18}|{6, 10}|{7, 10}",
                "i", "left", "right", "alpha", "maxKuttaError", "y(1)", "error", "y(1, alpha)-b"));

            do
            {
                string str = "";
                str += string.Format("{0,2}|{1, 12:N7}|{2, 11:N7}|", i++, start, end);
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
                y0=RungeKutt_solve(f1, f2, x0, a, mid, 1, RungeEps, ref text, 2);
                str += string.Format("{0, 10:N7}|{1, 21:N7}|{2, 18}|{3, 10:N7}|{4, 10:N7}", mid, text, y0, error, f(mid));
                strList.Add(str);
                //SolveBox.Items.Add(str);
            }
            while (error > boundEps);
            strList.Add("Alpha = " + mid);

            strList.Add("Y(1) = "+RungeKutt_solve(f1, f2, x0, a, mid, x, RungeEps, ref text, 1));

        }
        public double RungeKutt_solve(Func<double, double, double, double> f1,
            Func<double, double, double, double> f2, double x0, double y0, double z0, double x,
            double eps, ref string text, int Out)
        {
            var h = 0.1;
            var xj = x0;
            var yj = y0;
            var zj = z0;
            var maxError = double.MinValue;
            if (x == x0)
                return y0;
            if (Out == 1)
                strList.Add(string.Format("{0, 7};{1, 12};{2, 16};{3, 16};{4, 16}", "h", "x", "y(x)", "y*(x)", "error"));
            do
            {
                string str = "";
                if (h > x - xj)
                    h = x - xj;
                //Автоматический контроль точности
                var next = For_RungeKut(f1, f2, xj, yj, zj, h);//Рунге для текущего х и h
                var next_stage = For_RungeKut(f1, f2, xj, yj, zj, h / 2.0);//Рунге для текущего х и h/2
                var _next = For_RungeKut(f1, f2, xj + h / 2.0, next_stage.Item1, next_stage.Item2, h / 2.0);//Рунге для х+h/2 и h/2
                //
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
                        str += string.Format("{0, 16}|{1, 16}|{2, 16}|", yj, _next.Item1, Error);
                        strList.Add(str);
                    }
                }
                else h /= 2;//к контролю точности
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
            double xi, double yi, double zi, double h)
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
//--------------------------------------------------------------------------------------------------
//Штуки с теплопроводностью       
        public static double[] progonka(double[,] matrix, double[] rightPart)
        {
            Func<int, double> a = i => -matrix[i, i - 1];
            Func<int, double> c = i => matrix[i, i];
            Func<int, double> b = i => -matrix[i, i + 1];
            Func<int, double> f = i => rightPart[i];
            var n = matrix.GetLength(0);

            var alpha = new double[n];
            var betta = new double[n];
            alpha[1] = b(0) / c(0);
            betta[1] = f(0) / c(0);

            for (var i = 1; i < n - 1; i++)
            {
                alpha[i + 1] = b(i) / (c(i) - a(i) * alpha[i]);
                betta[i + 1] = (f(i) + a(i) * betta[i]) / (c(i) - a(i) * alpha[i]);
            }

            var x = new double[n];
            x[n - 1] = (f(n - 1) + a(n - 1) * betta[n - 1]) / (c(n - 1) - a(n - 1) * alpha[n - 1]);
            for (var i = n - 2; i >= 0; i--)
                x[i] = alpha[i + 1] * x[i + 1] + betta[i + 1];
            return x;
        }
        public double Delta(double[] vector,double t, bool Explic)
        {
            int n = vector.Length - 1;
            double h = 1.0 / n,tau;
            if (Explic)
                tau = h * h / (4 * 0.1);
            else
                tau = h;
            double max = 0;
            double temp;
            for (int j = 0; j < n; j++)
            {
                temp = Math.Abs(u(tau * t, j * h) - vector[j]);
                if (temp > max)
                    max = temp;
            }
            return max;
        }
        public void explicit_()
        {
            double[] vector = null;
            double max_Delta = 0;
            strList.Add("Явный случай");
            for (int n = 8; n <= 32; n *= 2)
            {
                strList.Add("N = " + n);
                max_Delta = 0;
                var end = false;
                for (int i = 0; tn_Ex(n, i) <= 1 || !end; i++)//убрать end
                {
                    vector = explicit_scheme(n, i);
                    var cur_Delta = Delta(vector, i, true);
                    strList.Add(string.Format("{0,10}|{1,15}",tn_Ex(n,i),cur_Delta));
                    if (tn_Ex(n, i) > 1)
                        end = true;
                    if (cur_Delta > max_Delta)
                        max_Delta = cur_Delta;
                }
                strList.Add("В конечный момент времени: ");
                for (int i = 0; i < vector.Length; i++)
                    strList.Add(vector[i].ToString());
                strList.Add("Максимальная делта = " + max_Delta.ToString());
            }
           
        }
        public double[] explicit_scheme(int n, int maxi)
        {
            double[] Un = new double[n + 1];
            double h = 1.0 / n;
            var tau = h * h / (4 * 0.1);
            for (int i = 0; i < Un.Length; i++)
            {
                Un[i] = i * h;//фи
            }
            double[] Un_new =(double[]) Un.Clone();
            for (int i = 0; i < maxi; i++)
            {
                for (int j = 1; j < Un.Length - 1; j++)
                    Un_new[j] = Un[j] + tau * (0.1 * (Un[j + 1] - 2 * Un[j] + Un[j - 1]) / (h * h) + f(tau * i, j * h));
                Un = (double[])Un_new.Clone();
            }
            return Un;
        }
        public void implicit_()
        {
            double[] vector = null;
            double max_Delta = 0;
            strList.Add("Явный случай");
            for (int n = 8; n <= 32; n *= 2)
            {
                strList.Add("N = " + n);
                max_Delta = 0;
                for (int i = 0; tn_Im(n, i) <= 1; i++)//убрать end
                {
                    vector = implicit_scheme(n, i);
                    var cur_Delta = Delta(vector, i, false);
                    strList.Add(string.Format("{0,10}|{1,15}", tn_Im(n, i), cur_Delta));
                    if (cur_Delta > max_Delta)
                        max_Delta = cur_Delta;
                }
                strList.Add("В конечный момент времени: ");
                for (int i = 0; i < vector.Length; i++)
                    strList.Add(vector[i].ToString());
                strList.Add("Максимальная делта = " + max_Delta.ToString());
            }
        }
        public double[] implicit_scheme(int n, int maxi)
        {
            double[] Un = new double[n + 1];
            double h = 1.0 / n;
            var tau = h;
            var d = (tau * 0.1) / (h * h);
            var matrix = new double[n - 1, n - 1];
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                if (i - 1 >= 0)
                    matrix[i, i - 1] = -d;
                matrix[i, i] = 1 + 2 * d;
                if (i + 1 < matrix.GetLength(1))
                    matrix[i, i + 1] = -d;
            }
            for (var j = 0; j < Un.Length; j++)
            {
                Un[j] = j * h;
            }
            Un[0] = a;
            Un[Un.Length - 1] = b;
            var rightPart = new double[Un.Length - 2];
            for (int i = 0; i < maxi; i++)
            {
                for (int j = 0; j < rightPart.Length; j++)
                    rightPart[j] = Un[j + 1] + tau * f(tau * (i + 1), h * (j + 1));
                rightPart[0] += d * a;
                rightPart[rightPart.Length - 1] += d * b;
                rightPart = progonka(matrix, rightPart);
                for (int j = 0; j < rightPart.Length; j++)
                    Un[j + 1] = rightPart[j];
            }
            return Un;
        }
        public double u(double t, double x)
        {
            if (rb_13.Checked)
                return x + 0.1 * t * Math.Sin(Math.PI * x) * 13;
            return x + 0.1 * t * Math.Sin(Math.PI * x) * 22;
        }
        public double f(double t, double x)
        {
            return dt_u(x) - 0.1 * dx2_u(t, x);
        }
        public double dx2_u(double t, double x)
        {
            if (rb_13.Checked)
                return -1.3 * Math.PI * Math.PI * t * Math.Sin(Math.PI * x);
            return -2.2 * Math.PI * Math.PI * t * Math.Sin(Math.PI * x);
        }
        Func<double, double> fi = x => x;
        public double dt_u(double x)
        {
            if (rb_13.Checked)
                return 1.3 * Math.Sin(Math.PI * x);
            return 2.2 * Math.Sin(Math.PI * x);
        }
        public double tn_Ex(int n, double i)
        {
            return 1.0 / (n * n) / (4.0 * 0.1) * i;
        }
        public double tn_Im(int n, double i)
        {
            return i * (1.0 / n);
        }
        //Решение уравнения
        private void SolveButton_Click(object sender, EventArgs e)
        {
            //RungeKutt();
            //explicit_();
            implicit_();
            SaveFile();
        }

       
    }
}
