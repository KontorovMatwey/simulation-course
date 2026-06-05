using System;
using System.Collections.Generic;
using System.Linq;

class Program
{
    static void Main()
    {
        double lambda = 0.7;   // интенсивность потока
        double T = 20.0;       // интервал времени
        int trials = 10000;    // число прогонов

        Random rnd = new Random(); // 42
        List<int> counts = new List<int>(trials);

        for (int i = 0; i < trials; i++)
        {
            counts.Add(SimulatePoissonCount(lambda, T, rnd));
        }

        double mean = counts.Average();
        double variance = counts.Select(x => (x - mean) * (x - mean)).Average();

        var distribution = counts
            .GroupBy(x => x)
            .OrderBy(g => g.Key)
            .Select(g => new
            {
                K = g.Key,
                P = (double)g.Count() / trials
            })
            .ToList();

        Console.WriteLine("Пуассоновский поток. Событие - поступление заявки");
        Console.WriteLine(new string('-', 55));
        Console.WriteLine("lambda = " + lambda);
        Console.WriteLine("T      = " + T);
        Console.WriteLine("trials = " + trials);
        Console.WriteLine();

        Console.WriteLine("Эмпирическое распределение числа заявок N(T):");
        Console.WriteLine("k\tP(N(T)=k)\tГистограмма");

        foreach (var item in distribution)
        {
            int barLen = (int)Math.Round(item.P * 100);
            string bar = new string('*', barLen);
            Console.WriteLine(item.K + "\t" + item.P.ToString("F4") + "\t\t" + bar);
        }

        Console.WriteLine();
        Console.WriteLine("Среднее число заявок:   " + mean.ToString("F4"));
        Console.WriteLine("Дисперсия числа заявок: " + variance.ToString("F4"));

        Console.WriteLine();
        Console.WriteLine("Теория по Пуассону:");
        Console.WriteLine("E[N(T)] = D[N(T)] = lambda * T = " + (lambda * T).ToString("F4"));
    }

    static int SimulatePoissonCount(double lambda, double T, Random rnd)
    {
        int count = 0;
        double time = 0.0;

        while (true)
        {
            double u = 1.0 - rnd.NextDouble();   // (0,1]
            double dt = -Math.Log(u) / lambda;   // экспоненциальный интервал
            time += dt;

            if (time > T)
                break;

            count++;
        }

        return count;
    }
}