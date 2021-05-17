using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using MathNet.Numerics.LinearAlgebra;

namespace lab3
{
    class Program
    {
        struct Point{
            public double ksi;
            public double etta;
        }
        static List<Vector<double>> v;
        static List<Vector<double>> x;
        static List<Vector<double>> x_extr;
        static List<Vector<double>> x_correct;
        static Matrix<double> F;
        static Matrix<double> H;
        static Matrix<double> B;
        static Vector<double> u;
        static Matrix<double> P;
        static Matrix<double> P_prev;
        static Matrix<double> P_extr;
        static Matrix<double> R;
        static Matrix<double> Q;
        static Matrix<double> I;
        static Matrix<double> K;
        static int steps, Xmax, Ymax, delta_t = 1;
        static double v_max, omega_max;
        static void Main(string[] args)
        {
            Console.OutputEncoding = Encoding.UTF8;
            x = new List<Vector<double>>();
            x_extr = new List<Vector<double>>();
            x_correct = new List<Vector<double>>();
            v = new List<Vector<double>>();
            Vector<double> c0 = Vector<double>.Build.Dense(3);
            Vector<double> m0 = Vector<double>.Build.Dense(2);
            Point A;
            Console.WriteLine("Введите ограничения области Xmax и Ymax:");
            Console.Write("Xmax = ");
            Xmax = Convert.ToInt32(Console.ReadLine());
            Console.Write("Ymax = ");
            Ymax = Convert.ToInt32(Console.ReadLine());
            Console.WriteLine("Введите ограничениe угловой скорости \u03C9max:");
            Console.Write("\u03C9max = ");
            omega_max = Convert.ToInt32(Console.ReadLine());
            Console.WriteLine("Введите ограничениe линейной скорости vmax:");
            Console.Write("vmax = ");
            v_max = Convert.ToInt32(Console.ReadLine());
            Console.WriteLine("Введите начальные координаты x, y и курс \u03F4");
            Console.Write("x = ");
            c0[0] = Convert.ToDouble(Console.ReadLine());
            Console.Write("y = ");
            c0[1] = Convert.ToDouble(Console.ReadLine());
            Console.Write("\u03F4 = ");
            c0[2] = Convert.ToDouble(Console.ReadLine());
            x.Add(c0);
            Console.Write("Введите начальную линейную скорость:\nv = ");
            m0[0] = Convert.ToDouble(Console.ReadLine());
            Console.Write("Введите начальную угловую скорость:\n\u03C9 = ");
            m0[1] = Convert.ToDouble(Console.ReadLine());
            v.Add(m0);
            Console.WriteLine("Введите начальные координаты \u03BE и \u03B7 точки A");
            Console.Write("\u03BE = ");
            A.ksi = Convert.ToDouble(Console.ReadLine());
            Console.Write("\u03B7 = ");
            A.etta = Convert.ToDouble(Console.ReadLine());
            Console.WriteLine("Введите количество шагов:");
            steps = Convert.ToInt32(Console.ReadLine());
            Initialization(steps);
            for (int i = 1; i < steps; i++){

                //экстраполяция состояния
                u = Vector<double>.Build.DenseOfArray(new[] {v[i][0], v[i][1] });
                F = Calculate_F(v[i], x[i]);
                B = Calculate_B(x[i]);
                Vector<double>  x_ht = (F * x[i - 1]) + (B * u);
                x_extr.Add(x_ht);

                //экстраполяция матрицы ковариации
                P_prev = P;
                P_extr = F * P_prev * F.Transpose();

                //Усиление по Калману
                H = Calculate_H(x[i], A);
                K = P_extr * H.Transpose() * (H * P_extr * H.Transpose() + R).Inverse();

                //Коррекция вектора состояния
                Vector<double> z_k = Calculate_z(A, x[i]);
                Vector<double> x_n = x_ht + K * (z_k - H * x_ht);
                x_correct.Add(x_n);

                //Расчет матрицы ковариации
                P = (I - K * H) * P_extr;
            }
            Save_file(steps, A);
            Print(steps);
            //Console.WriteLine("Запуск окна с графиком.");
            //ProcessStartInfo startInfo = new ProcessStartInfo("python");
            //Process process = new Process();
            //startInfo.Arguments = "main.py";
            //startInfo.UseShellExecute = false;
            //startInfo.CreateNoWindow = true;
            //startInfo.RedirectStandardError = true;
            //startInfo.RedirectStandardOutput = true;
            //process.StartInfo = startInfo;
            //process.Start();
            Console.ReadKey();
        }

        static void Initialization(int cnt){
            Random random = new Random();
            for (int i = 1; i < cnt; i++){
                Vector<double> m_r = Vector<double>.Build.Dense(2);
                Vector<double> a = Vector<double>.Build.Dense(3);
                while (true){
                    m_r[0] = random.NextDouble() * v_max;
                    if (Math.Abs(m_r[0] - v[i - 1][0]) < v[i - 1][0] * 0.1){
                        break;
                    }
                }
                while (true){
                    m_r[1] = -omega_max + random.NextDouble() * 2 * omega_max;
                    if (Math.Abs(m_r[1] - v[i - 1][1]) < v[i - 1][1] * 0.1){
                        v.Add(m_r);
                        break;
                    }
                }
                a[0] = Xmax * random.NextDouble();
                a[1] = Ymax * random.NextDouble();
                a[2] = random.NextDouble();
                x.Add(a);
            }
            P = Calculate_P(Xmax, Ymax);
            R = Calculate_R(Xmax, Ymax);
            Q = R;
            x_correct.Add(x[0]);
            I = Matrix<double>.Build.DenseIdentity(3);
        }

        static Matrix<double> Calculate_F(Vector<double> v, Vector<double> c){
            Matrix<double> calc;
            double a = -v[0] * Math.Sin(c[2]) * delta_t;
            double b = v[0] * Math.Cos(c[2]) * delta_t;
            calc = Matrix<double>.Build.DenseOfArray(new[,] {{1, 0, a},
                                                             {0, 1, b},
                                                             {0, 0, 1}});
            return calc;
        }

        static Matrix<double> Calculate_B(Vector<double> c){
            Matrix<double> calc;
            double a = Math.Sin(c[2]) * delta_t;
            double b = Math.Cos(c[2]) * delta_t;
            calc = Matrix<double>.Build.DenseOfArray(new[,] {{b, 0},
                                                             {a, 0},
                                                             {0, delta_t}});
            return calc;
        }

        static Matrix<double> Calculate_H(Vector<double> C, Point A){
            Matrix<double> calc;
            double a = (-2 * (A.ksi - C[0])) / (Math.Sqrt(Math.Pow(A.ksi - C[0], 2) + Math.Pow(A.etta - C[1], 2)));
            double b = (A.etta - C[1]) / (Math.Pow(A.ksi - C[0], 2) + Math.Pow(A.etta - C[1], 2));
            double c = (-2 * (A.etta - C[1])) / (Math.Sqrt(Math.Pow(A.ksi - C[0], 2) + Math.Pow(A.etta - C[1], 2)));
            double d = -Math.Pow(A.ksi - C[0], 2) / (Math.Pow(A.ksi - C[0], 2) + Math.Pow(A.etta - C[1], 2));
            calc = Matrix<double>.Build.DenseOfArray(new[,] {{a, c, 0},
                                                             {b, d, -1}});
            return calc;
        }

        static Matrix<double> Calculate_P(int xmax, int ymax){
            Matrix<double> calc;
            double sigma2_x = Math.Pow(xmax / 6, 2);
            double sigma2_y = Math.Pow(ymax / 6, 2);
            double sigma2_tetta = Math.Pow(Math.PI / 3, 2);
            calc = Matrix<double>.Build.DenseOfArray(new[,] {{sigma2_x, 0, 0},
                                                             {0, sigma2_y, 0},
                                                             {0, 0, sigma2_tetta}});
            return calc;
        }

        static Matrix<double> Calculate_R(int xmax, int ymax)
        {
            Matrix<double> calc;
            double sigma2_r = (Math.Pow(xmax, 2) + Math.Pow(ymax, 2)) / 36;
            double sigma2_fi = Math.Pow(Math.PI / 3, 2);
            calc = Matrix<double>.Build.DenseOfArray(new[,] {{sigma2_r, 0},
                                                             {0, sigma2_fi}});
            return calc;
        }

        static Vector<double> Calculate_z(Point a, Vector<double> c)
        {
            Vector<double> calc = Vector<double>.Build.Dense(2);
            calc[0] = Math.Sqrt(Math.Pow(a.ksi - c[0], 2) + Math.Pow(a.etta - c[1], 2));
            calc[1] = Math.Atan2(a.etta - c[1], a.ksi - c[0]) - c[2];
            return calc;
        }

        static void Save_file(int N, Point a)
        {
            string Patch = "data.txt";
            string x_coordinate = null;
            string y_coordinate = null;
            StreamWriter sw = null;

            try
            {
                sw = new StreamWriter(Patch, false, System.Text.Encoding.Default);
                Console.WriteLine("Файл создан. Запись данных...");
                x_coordinate = Convert.ToString(a.ksi);
                y_coordinate = Convert.ToString(a.etta);
                sw.WriteLine(x_coordinate);
                sw.WriteLine(y_coordinate);
                x_coordinate = null;
                y_coordinate = null;
                for (int i = 0; i < N - 1; i++)
                {
                    x_coordinate += Convert.ToString(x_extr[i][0]) + " ";
                    y_coordinate += Convert.ToString(x_extr[i][1]) + " ";
                }
                sw.WriteLine(x_coordinate);
                sw.WriteLine(y_coordinate);
                x_coordinate = null;
                y_coordinate = null;
                for (int i = 0; i < N; i++)
                {
                    x_coordinate += Convert.ToString(x[i][0]) + " ";
                    y_coordinate += Convert.ToString(x[i][1]) + " ";
                }
                sw.WriteLine(x_coordinate);
                sw.WriteLine(y_coordinate);
            }
            catch (Exception e)
            {
                Console.WriteLine("Ошибка в файле:" + e.Message);
            }
            finally
            {
                if (sw != null)
                {
                    sw.Close();
                    Console.WriteLine("Данные записаны.");
                }
            }
        }

        static void Print(int N){
            string row;
            Console.WriteLine("x_extr            |  y_extr            |  tetta_extr         |  x_corr            |  y_corr            |  tetta_corr");
            for (int i = 1; i < N - 1; i++){
                row = Convert.ToString(x_extr[i - 1][0]) + "  |  " + Convert.ToString(x_extr[i - 1][1]) + "  |  " +
                      Convert.ToString(x_extr[i - 1][2]) + "  |  " + Convert.ToString(x_correct[i][0]) + "  |  " +
                      Convert.ToString(x_correct[i][1]) + "  |  " + Convert.ToString(x_correct[i][2]);
                Console.WriteLine(row);
            }
        }
    }
}
