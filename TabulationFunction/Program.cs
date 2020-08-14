using System;
using System.Collections.Generic;
using System.IO;
namespace Laba
{

    class MainClass
    {

        public static double L(double X, double[] x, double[] f, int n)
        {
            double ans = f[0];
            double mulx = 1;
            for (int j = 1; j < n; j++)
            {
                mulx *= (X - x[j - 1]);
                ans += mulx * f[j];
            }
            return ans;
        }
        public static double F(double x, double eps)
        {
            double sum = 0;
            int sign = -1;
            double myx = x * x;
            double twon = 2;
            double facttwon = 2;
            double prevsum = 2 * eps;
            int n = 1;
            while (Math.Abs(prevsum - sum) > eps)
            {
                prevsum = sum;
                sum += sign * myx / twon / facttwon;
                sign *= -1;
                myx *= x * x;
                n++;
                twon = 2 * n;
                facttwon *= (2 * n * (2 * n - 1));
            }

            return sum;
        }


        public static void Newton(double a, double b, double eps)
        {

            double[] f = new double[100];
            double[] x = new double[100];

            for (int n = 6; n < 100; n++)
            {

                double h = (b - a) / (n - 1);
                for (int i = 0; i < n; i++)
                {
                    x[i] = a + h * i;

                    f[i] = F(x[i], eps);
                }
                for (int i = 1; i < n; i++)
                    for (int j = n - 1; j >= i; j--)
                        f[j] = (f[j] - f[j - 1]) / (x[j] - x[j - i]);

                int step = 11;
                h = (b - a) / (step - 1);
                double maxeps = 0;
                for (int i = 0; i < step; i++)
                {
                    double X = a + i * h;
                    //    Console.WriteLine("{0}\t {1}\t {2}\t {3}", X, F(X, eps), L(X, x, f, n), Math.Abs(F(X, eps) - L(X, x, f, n)));
                    maxeps = Math.Max(maxeps, Math.Abs(F(X, eps) - L(X, x, f, n)));
                }
                Console.WriteLine("{0:00.} {1:0.0000000000}", n, maxeps);
            }

        }
        public static void Chebishev(double a, double b, double eps)
        {
            double[] f = new double[100];
            double[] x = new double[100];

            for (int n = 6; n < 100; n++)
            {

                for (int i = 0; i < n; i++)
                {
                    x[i] = (b - a) / 2 * Math.Cos((2 * i + 1) * Math.PI / (2 * n + 2)) + (b + a) / 2;

                    f[i] = F(x[i], eps);
                }
                for (int i = 1; i < n; i++)
                    for (int j = n - 1; j >= i; j--)
                        f[j] = (f[j] - f[j - 1]) / (x[j] - x[j - i]);

                int step = 11;
                double h = (b - a) / (step - 1);
                double maxeps = 0;
                for (int i = 0; i < step; i++)
                {
                    double X = a + i * h;
                    //     Console.WriteLine("{0}\t {1}\t {2}\t {3}", X, F(X, eps), L(X, x, f, n), Math.Abs(F(X, eps) - L(X, x, f, n)));
                    maxeps = Math.Max(maxeps, Math.Abs(F(X, eps) - L(X, x, f, n)));
                }
                Console.WriteLine("{0:00.} {1:0.0000000000}", n, maxeps);
            }

        }

        public static void tab(double a, double b, double eps, double h)
        {
            for (double x = a; x - eps <= b; x += h)
                Console.WriteLine("{0:0.0} {1:0.0000000000}", x, F(x, eps));
        }

        public static double ff(double t)
        {
            if (t == 0)
                return 1;
            return (Math.Cos(t) - 1) / t;
        }
        static double z(double a, double h, int i)
        {
            return a + i * h;
        }
        public static double CenterSum(double a, double b, int n)
        {
            double h = (b - a) / n;
            double sum = 0;
            for (int i = 1; i <= n; i++)
            {

                sum += ff((z(a, h, i) + z(a, h, i - 1)) / 2);
            }
            return sum * h;
        }
        public static double LeftSum(double a, double b, int n)
        {
            double h = (b - a) / n;
            double sum = 0;
            for (int i = 1; i <= n; i++)
            {

                sum += ff(z(a, h, i - 1));
            }
            return sum * h;
        }
        public static double RightSum(double a, double b, int n)
        {
            double h = (b - a) / n;
            double sum = 0;
            for (int i = 1; i <= n; i++)
            {

                sum += ff(z(a, h, i));
            }
            return sum * h;
        }

        public static double TrapSum(double a, double b, int n)
        {
            double h = (b - a) / n;
            double sum = 0;
            for (int i = 1; i < n; i++)
            {
                sum += ff(z(a, h, i));
            }
            sum *= 2; sum += ff(z(a, h, 0)) + ff(z(a, h, n));
            return sum * h / 2;
        }
        public static double SimpsSum(double a, double b, int n)
        {
            double h = (b - a) / n;

            double sum1 = 0;
            double sum2 = 0;
            for (int i = 1; i < n; i += 2)
            {
                sum1 += ff(z(a, h, i));

            }
            for (int i = 2; i < n; i += 2)
            {
                sum2 += ff(z(a, h, i));

            }
            sum1 *= 4;
            sum2 *= 2;
            double sum = ff(z(a, h, 0)) + ff(z(a, h, n)) + sum1 + sum2;
            return sum * h / 3;
        }


        public static double GaussSum(double a, double b, int n)

        {

            double h = (b - a) / n;
            double sum = 0;
            double sqrt3 = Math.Sqrt(3);
            for (int i = 1; i <= n; i++)
            {


                sum += ff(z(a, h, i - 1) + h / 2 * (1 - 1 / sqrt3));

                sum += ff(z(a, h, i - 1) + h / 2 * (1 + 1 / sqrt3));

            }

            return sum * h / 2;
        }
        public static KeyValuePair<double, int> CenterRect(double a, double b, double eps)
        {
            int n = 2;
            while (Math.Abs(CenterSum(a, b, 2 * n) - CenterSum(a, b, n)) >= eps)
                n *= 2;
            return new KeyValuePair<double, int>(CenterSum(a, b, n), n);
        }

        public static KeyValuePair<double, int> LeftRect(double a, double b, double eps)
        {
            int n = 2;
            while (Math.Abs(LeftSum(a, b, 2 * n) - LeftSum(a, b, n)) >= eps)
                n *= 2;
            return new KeyValuePair<double, int>(LeftSum(a, b, n), n);
        }

        public static KeyValuePair<double, int> RightRect(double a, double b, double eps)
        {
            int n = 2;
            while (Math.Abs(RightSum(a, b, 2 * n) - RightSum(a, b, n)) >= eps)
                n *= 2;
            return new KeyValuePair<double, int>(RightSum(a, b, n), n);
        }

        public static KeyValuePair<double, int> Trap(double a, double b, double eps)
        {
            int n = 2;
            while (Math.Abs(TrapSum(a, b, 2 * n) - TrapSum(a, b, n)) >= eps)
                n *= 2;
            return new KeyValuePair<double, int>(TrapSum(a, b, n), n);
        }

        public static KeyValuePair<double, int> Simps(double a, double b, double eps)
        {
            int n = 2;
            while (Math.Abs(SimpsSum(a, b, 2 * n) - SimpsSum(a, b, n)) >= eps)
                n *= 2;
            return new KeyValuePair<double, int>(SimpsSum(a, b, n), n);
        }
        public static KeyValuePair<double, int> Gauss(double a, double b, double eps)
        {
            int n = 2;
            while (Math.Abs(GaussSum(a, b, 2 * n) - GaussSum(a, b, n)) >= eps)
                n *= 2;
            return new KeyValuePair<double, int>(GaussSum(a, b, n), n);
        }

        public static void CalcIntegral(double a, double b, double h, double eps)
        {
            Console.WriteLine("Центральные прямоугольники");
            for (double x = a; x - eps <= b; x += h)
            {
                var ans = CenterRect(0, x, eps);
                Console.WriteLine("{0:0.0} {1:0.0000} {2:0.0000} {3}", x, F(x, eps), ans.Key, ans.Value);
            }
            Console.WriteLine("Левые прямоугольники");
            for (double x = a; x - eps <= b; x += h)
            {
                var ans = LeftRect(0, x, eps);
                Console.WriteLine("{0:0.0} {1:0.0000} {2:0.0000} {3}", x, F(x, eps), ans.Key, ans.Value);
            }

            Console.WriteLine("Правые прямоугольники");
            for (double x = a; x - eps <= b; x += h)
            {
                var ans = RightRect(0, x, eps);
                Console.WriteLine("{0:0.0} {1:0.0000} {2:0.0000} {3}", x, F(x, eps), ans.Key, ans.Value);
            }

            Console.WriteLine("Трапеции");
            for (double x = a; x - eps <= b; x += h)
            {
                var ans = Trap(0, x, eps);
                Console.WriteLine("{0:0.0} {1:0.0000} {2:0.0000} {3}", x, F(x, eps), ans.Key, ans.Value);
            }

            Console.WriteLine("Симпсон");
            for (double x = a; x - eps <= b; x += h)
            {
                var ans = Simps(0, x, eps);
                Console.WriteLine("{0:0.0} {1:0.0000} {2:0.0000} {3}", x, F(x, eps), ans.Key, ans.Value);
            }

            Console.WriteLine("Гаусс");
            for (double x = a; x - eps <= b; x += h)
            {
                var ans = Gauss(0, x, eps);
                Console.WriteLine("{0:0.0} {1:0.0000} {2:0.0000} {3}", x, F(x, eps), ans.Key, ans.Value);
            }
        }
        public static void Main(string[] args)
        {
            double a = 0.4, b = 4, h = 0.2, eps = 1e-6;
            Console.WriteLine("Табуляция");
            tab(a, b, eps, h);
            Console.WriteLine("Ньютон");

            Newton(a, b, eps);
            Console.WriteLine("Чебышев");
            Chebishev(a, b, eps);


            CalcIntegral(a, b, h, 1e-5);
        }
    }
}

