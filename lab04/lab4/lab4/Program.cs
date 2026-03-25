using System;

class Program
{
    // Параметры LCG (как на слайде)
    const ulong M = 1UL << 63;                 // 2^63
    const ulong BETA = (1UL << 32) + 3;        // 2^32 + 3

    static void Main()
    {
        const int N = 100000;

        // --- Сид от времени ---
        ulong seed = (ulong)DateTime.Now.Ticks;

        Console.WriteLine("Seed: " + seed);
        Console.WriteLine();

        // --- Массивы ---
        double[] lcgValues = new double[N];
        double[] builtinValues = new double[N];

        // --- Генерация LCG ---
        ulong x = seed;

        for (int i = 0; i < N; i++)
        {
            x = (BETA * x) % M;
            lcgValues[i] = (double)x / M;
        }

        // --- Встроенный генератор ---
        Random rnd = new Random((int)(seed % int.MaxValue));

        for (int i = 0; i < N; i++)
        {
            builtinValues[i] = rnd.NextDouble();
        }

        // --- Статистика ---
        var (meanLCG, varLCG) = CalculateStats(lcgValues);
        var (meanBuilt, varBuilt) = CalculateStats(builtinValues);

        // --- Теория ---
        double theoreticalMean = 0.5;
        double theoreticalVar = 1.0 / 12.0;

        // --- Вывод ---
        Console.WriteLine("=== LCG ===");
        Console.WriteLine($"Mean: {meanLCG}");
        Console.WriteLine($"Variance: {varLCG}");

        Console.WriteLine("\n=== Built-in Random ===");
        Console.WriteLine($"Mean: {meanBuilt}");
        Console.WriteLine($"Variance: {varBuilt}");

        Console.WriteLine("\n=== Theoretical ===");
        Console.WriteLine($"Mean: {theoreticalMean}");
        Console.WriteLine($"Variance: {theoreticalVar}");

        // --- Первые 10 значений ---
        Console.WriteLine("\nFirst 10 LCG values:");
        for (int i = 0; i < 10; i++)
            Console.WriteLine(lcgValues[i]);

        Console.WriteLine("\nFirst 10 Built-in values:");
        for (int i = 0; i < 10; i++)
            Console.WriteLine(builtinValues[i]);
    }

    static (double mean, double variance) CalculateStats(double[] data)
    {
        int n = data.Length;

        double sum = 0;
        for (int i = 0; i < n; i++)
            sum += data[i];

        double mean = sum / n;

        double varSum = 0;
        for (int i = 0; i < n; i++)
            varSum += Math.Pow(data[i] - mean, 2);

        double variance = varSum / n;

        return (mean, variance);
    }
}