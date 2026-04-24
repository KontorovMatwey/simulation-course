using System;
using System.Drawing;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace lab06
{
    public partial class MainForm : Form
    {
        private readonly Random _rnd = new Random();

        // Задание 1
        private TextBox[] _probBoxes;
        private TextBox _discreteNBox;
        private TextBox _discreteOutput;
        private Label _discreteStatus;
        private Chart _discreteChart;

        // Задание 2
        private TextBox _meanBox;
        private TextBox _varianceBox;
        private TextBox _normalNBox;
        private TextBox _normalOutput;
        private Label _normalStatus;
        private Chart _normalChart;

        public MainForm()
        {
            BuildUi();
            SetDefaults();
        }

        private void BuildUi()
        {
            Text = "Моделирование случайных величин";
            Width = 1700;
            Height = 950;
            StartPosition = FormStartPosition.CenterScreen;
            BackColor = Color.WhiteSmoke;

            Controls.Clear();

            var root = new TableLayoutPanel
            {
                Dock = DockStyle.Fill,
                ColumnCount = 2,
                RowCount = 1,
                BackColor = Color.WhiteSmoke
            };

            root.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 50f));
            root.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 50f));
            root.RowStyles.Add(new RowStyle(SizeType.Percent, 100f));

            root.Controls.Add(BuildDiscretePanel(), 0, 0);
            root.Controls.Add(BuildNormalPanel(), 1, 0);

            Controls.Add(root);
        }

        private Control BuildDiscretePanel()
        {
            var panel = new TableLayoutPanel
            {
                Dock = DockStyle.Fill,
                ColumnCount = 1,
                RowCount = 3,
                Padding = new Padding(8),
                BackColor = Color.WhiteSmoke
            };

            panel.RowStyles.Add(new RowStyle(SizeType.Absolute, 170));
            panel.RowStyles.Add(new RowStyle(SizeType.Percent, 55f));
            panel.RowStyles.Add(new RowStyle(SizeType.Percent, 45f));

            var inputGroup = new GroupBox
            {
                Text = "Задание 1",
                Dock = DockStyle.Fill,
                Padding = new Padding(10)
            };

            var inputGrid = new TableLayoutPanel
            {
                Dock = DockStyle.Fill,
                ColumnCount = 8,
                RowCount = 2
            };

            for (int i = 0; i < 8; i++)
                inputGrid.ColumnStyles.Add(new ColumnStyle(SizeType.AutoSize));

            inputGrid.RowStyles.Add(new RowStyle(SizeType.Absolute, 40));
            inputGrid.RowStyles.Add(new RowStyle(SizeType.Absolute, 40));

            _probBoxes = new TextBox[5];

            for (int i = 0; i < 5; i++)
            {
                inputGrid.Controls.Add(new Label
                {
                    Text = $"P{i + 1}",
                    AutoSize = true,
                    Margin = new Padding(8, 12, 4, 0)
                }, i * 2, 0);

                _probBoxes[i] = new TextBox
                {
                    Width = 70,
                    Margin = new Padding(0, 8, 10, 0)
                };
                inputGrid.Controls.Add(_probBoxes[i], i * 2 + 1, 0);
            }

            inputGrid.Controls.Add(new Label
            {
                Text = "N",
                AutoSize = true,
                Margin = new Padding(8, 12, 4, 0)
            }, 0, 1);

            _discreteNBox = new TextBox
            {
                Width = 90,
                Margin = new Padding(0, 8, 10, 0)
            };
            inputGrid.Controls.Add(_discreteNBox, 1, 1);

            var btnCustom = new Button
            {
                Text = "Run custom",
                Width = 110,
                Height = 28,
                Margin = new Padding(8, 8, 8, 0)
            };
            btnCustom.Click += (_, __) => RunDiscreteCustom();

            var btnBasic = new Button
            {
                Text = "Run basic",
                Width = 110,
                Height = 28,
                Margin = new Padding(8, 8, 8, 0)
            };
            btnBasic.Click += (_, __) => RunDiscreteBasic();

            _discreteStatus = new Label
            {
                Text = "Ready",
                AutoSize = true,
                Margin = new Padding(8, 12, 8, 0)
            };

            inputGrid.Controls.Add(btnCustom, 2, 1);
            inputGrid.Controls.Add(btnBasic, 3, 1);
            inputGrid.Controls.Add(_discreteStatus, 4, 1);

            inputGroup.Controls.Add(inputGrid);

            _discreteChart = CreateChart("Discrete");
            _discreteOutput = CreateOutputBox();

            panel.Controls.Add(inputGroup, 0, 0);
            panel.Controls.Add(_discreteChart, 0, 1);
            panel.Controls.Add(_discreteOutput, 0, 2);

            return panel;
        }

        private Control BuildNormalPanel()
        {
            var panel = new TableLayoutPanel
            {
                Dock = DockStyle.Fill,
                ColumnCount = 1,
                RowCount = 3,
                Padding = new Padding(8),
                BackColor = Color.WhiteSmoke
            };

            panel.RowStyles.Add(new RowStyle(SizeType.Absolute, 140));
            panel.RowStyles.Add(new RowStyle(SizeType.Percent, 55f));
            panel.RowStyles.Add(new RowStyle(SizeType.Percent, 45f));

            var inputGroup = new GroupBox
            {
                Text = "Задание 2",
                Dock = DockStyle.Fill,
                Padding = new Padding(10)
            };

            var inputGrid = new TableLayoutPanel
            {
                Dock = DockStyle.Fill,
                ColumnCount = 6,
                RowCount = 2
            };

            for (int i = 0; i < 6; i++)
                inputGrid.ColumnStyles.Add(new ColumnStyle(SizeType.AutoSize));

            inputGrid.RowStyles.Add(new RowStyle(SizeType.Absolute, 40));
            inputGrid.RowStyles.Add(new RowStyle(SizeType.Absolute, 40));

            inputGrid.Controls.Add(new Label
            {
                Text = "Mean",
                AutoSize = true,
                Margin = new Padding(8, 12, 4, 0)
            }, 0, 0);

            _meanBox = new TextBox
            {
                Width = 80,
                Margin = new Padding(0, 8, 10, 0)
            };
            inputGrid.Controls.Add(_meanBox, 1, 0);

            inputGrid.Controls.Add(new Label
            {
                Text = "Variance",
                AutoSize = true,
                Margin = new Padding(8, 12, 4, 0)
            }, 2, 0);

            _varianceBox = new TextBox
            {
                Width = 80,
                Margin = new Padding(0, 8, 10, 0)
            };
            inputGrid.Controls.Add(_varianceBox, 3, 0);

            inputGrid.Controls.Add(new Label
            {
                Text = "N",
                AutoSize = true,
                Margin = new Padding(8, 12, 4, 0)
            }, 0, 1);

            _normalNBox = new TextBox
            {
                Width = 90,
                Margin = new Padding(0, 8, 10, 0)
            };
            inputGrid.Controls.Add(_normalNBox, 1, 1);

            var btnCustom = new Button
            {
                Text = "Run custom",
                Width = 110,
                Height = 28,
                Margin = new Padding(8, 8, 8, 0)
            };
            btnCustom.Click += (_, __) => RunNormalCustom();

            var btnBasic = new Button
            {
                Text = "Run basic",
                Width = 110,
                Height = 28,
                Margin = new Padding(8, 8, 8, 0)
            };
            btnBasic.Click += (_, __) => RunNormalBasic();

            _normalStatus = new Label
            {
                Text = "Ready",
                AutoSize = true,
                Margin = new Padding(8, 12, 8, 0)
            };

            inputGrid.Controls.Add(btnCustom, 2, 1);
            inputGrid.Controls.Add(btnBasic, 3, 1);
            inputGrid.Controls.Add(_normalStatus, 4, 1);

            inputGroup.Controls.Add(inputGrid);

            _normalChart = CreateChart("Normal");
            _normalOutput = CreateOutputBox();

            panel.Controls.Add(inputGroup, 0, 0);
            panel.Controls.Add(_normalChart, 0, 1);
            panel.Controls.Add(_normalOutput, 0, 2);

            return panel;
        }

        private Chart CreateChart(string title)
        {
            var chart = new Chart
            {
                Dock = DockStyle.Fill,
                BackColor = Color.White,
                MinimumSize = new Size(0, 280)
            };

            var area = new ChartArea("area");
            area.BackColor = Color.White;
            area.AxisX.MajorGrid.LineColor = Color.Gainsboro;
            area.AxisY.MajorGrid.LineColor = Color.Gainsboro;
            area.AxisX.LabelStyle.Font = new Font("Segoe UI", 9f);
            area.AxisY.LabelStyle.Font = new Font("Segoe UI", 9f);
            area.AxisX.IntervalAutoMode = IntervalAutoMode.VariableCount;
            area.AxisY.IntervalAutoMode = IntervalAutoMode.VariableCount;
            chart.ChartAreas.Add(area);

            chart.Legends.Add(new Legend("legend")
            {
                Docking = Docking.Top,
                Font = new Font("Segoe UI", 9f)
            });

            chart.Titles.Add(new Title(title)
            {
                Font = new Font("Segoe UI", 10f, FontStyle.Bold)
            });

            return chart;
        }

        private TextBox CreateOutputBox()
        {
            return new TextBox
            {
                Dock = DockStyle.Fill,
                Multiline = true,
                ScrollBars = ScrollBars.Vertical,
                ReadOnly = true,
                Font = new Font("Consolas", 9.5f),
                BackColor = Color.White
            };
        }

        private void SetDefaults()
        {
            _probBoxes[0].Text = "20";
            _probBoxes[1].Text = "20";
            _probBoxes[2].Text = "20";
            _probBoxes[3].Text = "20";
            _probBoxes[4].Text = "20";

            _discreteNBox.Text = "1000";

            _meanBox.Text = "0";
            _varianceBox.Text = "1";
            _normalNBox.Text = "1000";
        }

        private void RunDiscreteCustom()
        {
            if (!TryReadProbabilities(out var p, out string status))
            {
                _discreteStatus.Text = status;
                return;
            }

            if (!TryReadPositiveInt(_discreteNBox.Text, out int n))
            {
                _discreteStatus.Text = "Bad N";
                return;
            }

            _discreteStatus.Text = status;
            var x = new[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            var result = SimulateDiscrete(n, x, p);

            var sb = new StringBuilder();
            AppendDiscreteReport(sb, n, x, p, result);

            _discreteOutput.Text = sb.ToString();
            DrawDiscreteChart(x, p, result.EmpiricalProbabilities);
        }

        private void RunDiscreteBasic()
        {
            if (!TryReadProbabilities(out var p, out string status))
            {
                _discreteStatus.Text = status;
                return;
            }

            _discreteStatus.Text = status;

            int[] ns = { 10, 100, 1000, 10000 };
            var x = new[] { 1.0, 2.0, 3.0, 4.0, 5.0 };

            var sb = new StringBuilder();
            sb.AppendLine("Discrete basic run");
            sb.AppendLine();

            DiscreteResult last = null;
            int lastN = 0;

            foreach (int n in ns)
            {
                last = SimulateDiscrete(n, x, p);
                lastN = n;
                AppendDiscreteLine(sb, n, x, p, last);
            }

            if (last != null)
            {
                sb.AppendLine();
                sb.AppendLine("Last N result:");
                AppendDiscreteReport(sb, lastN, x, p, last);
                DrawDiscreteChart(x, p, last.EmpiricalProbabilities);
            }

            _discreteOutput.Text = sb.ToString();
        }

        private void RunNormalCustom()
        {
            if (!TryReadDouble(_meanBox.Text, out double mu))
            {
                _normalStatus.Text = "Bad mean";
                return;
            }

            if (!TryReadDouble(_varianceBox.Text, out double variance) || variance <= 0)
            {
                _normalStatus.Text = "Bad variance";
                return;
            }

            if (!TryReadPositiveInt(_normalNBox.Text, out int n))
            {
                _normalStatus.Text = "Bad N";
                return;
            }

            _normalStatus.Text = "Box-Muller";

            double sigma = Math.Sqrt(variance);
            var result = SimulateNormal(n, mu, sigma);

            var sb = new StringBuilder();
            AppendNormalReport(sb, n, mu, variance, result);

            _normalOutput.Text = sb.ToString();
            DrawNormalChart(result.Samples, mu, sigma);
        }

        private void RunNormalBasic()
        {
            if (!TryReadDouble(_meanBox.Text, out double mu))
            {
                _normalStatus.Text = "Bad mean";
                return;
            }

            if (!TryReadDouble(_varianceBox.Text, out double variance) || variance <= 0)
            {
                _normalStatus.Text = "Bad variance";
                return;
            }

            _normalStatus.Text = "Box-Muller";

            int[] ns = { 10, 100, 1000, 10000 };
            double sigma = Math.Sqrt(variance);

            var sb = new StringBuilder();
            sb.AppendLine("Normal basic run");
            sb.AppendLine();

            NormalResult last = null;
            int lastN = 0;

            foreach (int n in ns)
            {
                last = SimulateNormal(n, mu, sigma);
                lastN = n;
                AppendNormalLine(sb, n, mu, variance, last);
            }

            if (last != null)
            {
                sb.AppendLine();
                sb.AppendLine("Last N result:");
                AppendNormalReport(sb, lastN, mu, variance, last);
                DrawNormalChart(last.Samples, mu, sigma);
            }

            _normalOutput.Text = sb.ToString();
        }

        private void AppendDiscreteReport(StringBuilder sb, int n, double[] x, double[] p, DiscreteResult result)
        {
            double theoryMean = CalcDiscreteMean(x, p);
            double theoryVar = CalcDiscreteVariance(x, p, theoryMean);

            sb.AppendLine($"N = {n}");
            sb.AppendLine($"Theory mean = {theoryMean:F6}");
            sb.AppendLine($"Sample mean = {result.SampleMean:F6} (Err = {RelativeError(result.SampleMean, theoryMean)})");
            sb.AppendLine($"Theory var  = {theoryVar:F6}");
            sb.AppendLine($"Sample var  = {result.SampleVariance:F6} (Err = {RelativeError(result.SampleVariance, theoryVar)})");
            sb.AppendLine($"Chi-Squared = {result.ChiSquare:F4} > {result.Critical:F3} => {result.ChiSquare > result.Critical}");
            sb.AppendLine();
        }

        private void AppendDiscreteLine(StringBuilder sb, int n, double[] x, double[] p, DiscreteResult result)
        {
            double theoryMean = CalcDiscreteMean(x, p);
            double theoryVar = CalcDiscreteVariance(x, p, theoryMean);

            sb.AppendLine(
                $"N={n,-6} " +
                $"Avg={result.SampleMean:F6} (Err={RelativeError(result.SampleMean, theoryMean)})  " +
                $"Var={result.SampleVariance:F6} (Err={RelativeError(result.SampleVariance, theoryVar)})  " +
                $"Chi2={result.ChiSquare:F4} > {result.Critical:F3} => {result.ChiSquare > result.Critical}");
        }

        private void AppendNormalReport(StringBuilder sb, int n, double mu, double variance, NormalResult result)
        {
            sb.AppendLine($"N = {n}");
            sb.AppendLine($"Theory mean = {mu:F6}");
            sb.AppendLine($"Sample mean = {result.SampleMean:F6} (Err = {RelativeError(result.SampleMean, mu)})");
            sb.AppendLine($"Theory var  = {variance:F6}");
            sb.AppendLine($"Sample var  = {result.SampleVariance:F6} (Err = {RelativeError(result.SampleVariance, variance)})");
            sb.AppendLine($"Chi-Squared = {result.ChiSquare:F4} > {result.Critical:F3} => {result.ChiSquare > result.Critical}");
            sb.AppendLine();
        }

        private void AppendNormalLine(StringBuilder sb, int n, double mu, double variance, NormalResult result)
        {
            sb.AppendLine(
                $"N={n,-6} " +
                $"Avg={result.SampleMean:F6} (Err={RelativeError(result.SampleMean, mu)})  " +
                $"Var={result.SampleVariance:F6} (Err={RelativeError(result.SampleVariance, variance)})  " +
                $"Chi2={result.ChiSquare:F4} > {result.Critical:F3} => {result.ChiSquare > result.Critical}");
        }

        private bool TryReadProbabilities(out double[] p, out string status)
        {
            p = new double[5];
            double[] raw = new double[5];

            for (int i = 0; i < 5; i++)
            {
                if (!TryReadDouble(_probBoxes[i].Text, out raw[i]))
                {
                    status = $"Bad P{i + 1}";
                    return false;
                }

                if (raw[i] < 0)
                {
                    status = $"P{i + 1} < 0";
                    return false;
                }
            }

            double sum = raw.Sum();
            if (sum <= 0)
            {
                status = "Sum <= 0";
                return false;
            }

            for (int i = 0; i < 5; i++)
                p[i] = raw[i] / sum;

            status = Math.Abs(sum - 1.0) < 1e-9 ? "Ready" : "Normalized";
            return true;
        }

        private DiscreteResult SimulateDiscrete(int n, double[] x, double[] p)
        {
            int k = p.Length;
            int[] counts = new int[k];
            double sum = 0.0;
            double sum2 = 0.0;

            for (int i = 0; i < n; i++)
            {
                double r = _rnd.NextDouble();
                double acc = 0.0;
                int idx = k - 1;

                for (int j = 0; j < k; j++)
                {
                    acc += p[j];
                    if (r < acc)
                    {
                        idx = j;
                        break;
                    }
                }

                counts[idx]++;
                double v = x[idx];
                sum += v;
                sum2 += v * v;
            }

            double mean = sum / n;
            double variance = Math.Max(0.0, sum2 / n - mean * mean);

            double[] empirical = counts.Select(c => c / (double)n).ToArray();

            double chi2 = 0.0;
            for (int i = 0; i < k; i++)
            {
                double expected = n * p[i];
                if (expected > 0)
                {
                    double diff = counts[i] - expected;
                    chi2 += diff * diff / expected;
                }
            }

            return new DiscreteResult
            {
                SampleMean = mean,
                SampleVariance = variance,
                ChiSquare = chi2,
                Critical = ChiSquareCritical05(k - 1),
                EmpiricalProbabilities = empirical
            };
        }

        private NormalResult SimulateNormal(int n, double mu, double sigma)
        {
            int bins = 10;
            int[] counts = new int[bins];
            double[] bounds = new double[bins + 1];
            bounds[0] = double.NegativeInfinity;
            bounds[bins] = double.PositiveInfinity;

            for (int i = 1; i < bins; i++)
                bounds[i] = mu + sigma * NormalQuantile(i / (double)bins);

            double sum = 0.0;
            double sum2 = 0.0;
            double[] samples = new double[n];

            for (int i = 0; i < n; i++)
            {
                double z = BoxMuller();
                double x = mu + sigma * z;
                samples[i] = x;

                sum += x;
                sum2 += x * x;

                int bin = FindBin(x, bounds);
                counts[bin]++;
            }

            double mean = sum / n;
            double variance = Math.Max(0.0, sum2 / n - mean * mean);

            double expected = n / (double)bins;
            double chi2 = 0.0;

            for (int i = 0; i < bins; i++)
            {
                double diff = counts[i] - expected;
                chi2 += diff * diff / expected;
            }

            return new NormalResult
            {
                Samples = samples,
                SampleMean = mean,
                SampleVariance = variance,
                ChiSquare = chi2,
                Critical = ChiSquareCritical05(bins - 1)
            };
        }

        private void DrawDiscreteChart(double[] x, double[] theoretical, double[] empirical)
        {
            _discreteChart.Series.Clear();

            var area = _discreteChart.ChartAreas[0];
            area.AxisX.Minimum = 0.5;
            area.AxisX.Maximum = 5.5;
            area.AxisY.Minimum = 0;
            area.AxisY.Maximum = Math.Max(theoretical.Max(), empirical.Max()) * 1.2;

            var s1 = new Series("Theory")
            {
                ChartType = SeriesChartType.Column,
                Color = Color.DeepSkyBlue,
                BorderColor = Color.RoyalBlue,
                BorderWidth = 1
            };

            var s2 = new Series("Empirical")
            {
                ChartType = SeriesChartType.Column,
                Color = Color.LightBlue,
                BorderColor = Color.SteelBlue,
                BorderWidth = 1
            };

            for (int i = 0; i < x.Length; i++)
            {
                s1.Points.AddXY(x[i], theoretical[i]);
                s2.Points.AddXY(x[i], empirical[i]);
            }

            _discreteChart.Series.Add(s1);
            _discreteChart.Series.Add(s2);
        }

        private void DrawNormalChart(double[] samples, double mu, double sigma)
        {
            _normalChart.Series.Clear();

            double min = samples.Min();
            double max = samples.Max();

            if (Math.Abs(max - min) < 1e-12)
            {
                min -= 1;
                max += 1;
            }

            double pad = (max - min) * 0.15;
            min -= pad;
            max += pad;

            int bins = 20;
            double width = (max - min) / bins;

            int[] counts = new int[bins];
            for (int i = 0; i < samples.Length; i++)
            {
                int idx = (int)((samples[i] - min) / width);
                if (idx < 0) idx = 0;
                if (idx >= bins) idx = bins - 1;
                counts[idx]++;
            }

            var hist = new Series("Histogram")
            {
                ChartType = SeriesChartType.Column,
                Color = Color.LightSkyBlue,
                BorderColor = Color.RoyalBlue,
                BorderWidth = 1
            };

            var pdf = new Series("Theory pdf")
            {
                ChartType = SeriesChartType.Line,
                Color = Color.MidnightBlue,
                BorderWidth = 2
            };

            for (int i = 0; i < bins; i++)
            {
                double center = min + (i + 0.5) * width;
                double density = counts[i] / (samples.Length * width);

                hist.Points.AddXY(center, density);
                pdf.Points.AddXY(center, NormalPdf(center, mu, sigma));
            }

            var area = _normalChart.ChartAreas[0];
            area.AxisX.Minimum = min;
            area.AxisX.Maximum = max;
            area.AxisY.Minimum = 0;
            area.AxisY.Maximum = Math.Max(hist.Points.Max(p => p.YValues[0]), pdf.Points.Max(p => p.YValues[0])) * 1.2;

            _normalChart.Series.Add(hist);
            _normalChart.Series.Add(pdf);
        }

        private int FindBin(double value, double[] bounds)
        {
            for (int i = 0; i < bounds.Length - 1; i++)
            {
                if (value >= bounds[i] && value < bounds[i + 1])
                    return i;
            }
            return bounds.Length - 2;
        }

        private double BoxMuller()
        {
            double u1 = 1.0 - _rnd.NextDouble();
            double u2 = 1.0 - _rnd.NextDouble();
            return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
        }

        private double NormalPdf(double x, double mu, double sigma)
        {
            double z = (x - mu) / sigma;
            return Math.Exp(-0.5 * z * z) / (sigma * Math.Sqrt(2.0 * Math.PI));
        }

        private double NormalQuantile(double p)
        {
            if (p <= 0.0 || p >= 1.0)
                throw new ArgumentOutOfRangeException(nameof(p));

            double[] a =
            {
                -3.969683028665376e+01,
                 2.209460984245205e+02,
                -2.759285104469687e+02,
                 1.383577518672690e+02,
                -3.066479806614716e+01,
                 2.506628277459239e+00
            };

            double[] b =
            {
                -5.447609879822406e+01,
                 1.615858368580409e+02,
                -1.556989798598866e+02,
                 6.680131188771972e+01,
                -1.328068155288572e+01
            };

            double[] c =
            {
                -7.784894002430293e-03,
                -3.223964580411365e-01,
                -2.400758277161838e+00,
                -2.549732539343734e+00,
                 4.374664141464968e+00,
                 2.938163982698783e+00
            };

            double[] d =
            {
                 7.784695709041462e-03,
                 3.224671290700398e-01,
                 2.445134137142996e+00,
                 3.754408661907416e+00
            };

            double plow = 0.02425;
            double phigh = 1.0 - plow;
            double q, r;

            if (p < plow)
            {
                q = Math.Sqrt(-2.0 * Math.Log(p));
                return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
                       ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
            }

            if (p > phigh)
            {
                q = Math.Sqrt(-2.0 * Math.Log(1.0 - p));
                return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
                        ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
            }

            q = p - 0.5;
            r = q * q;

            return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
                   (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
        }

        private double ChiSquareCritical05(int df)
        {
            switch (df)
            {
                case 1: return 3.841;
                case 2: return 5.991;
                case 3: return 7.815;
                case 4: return 9.488;
                case 5: return 11.070;
                case 6: return 12.592;
                case 7: return 14.067;
                case 8: return 15.507;
                case 9: return 16.919;
                case 10: return 18.307;
                default: return 0.0;
            }
        }

        private string RelativeError(double sample, double theory)
        {
            if (Math.Abs(theory) < 1e-12)
                return "n/a";

            double err = Math.Abs(sample - theory) / Math.Abs(theory) * 100.0;
            return $"{err:F2}%";
        }

        private bool TryReadDouble(string text, out double value)
        {
            text = (text ?? "").Trim().Replace(',', '.');
            return double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out value);
        }

        private bool TryReadPositiveInt(string text, out int value)
        {
            text = (text ?? "").Trim();
            return int.TryParse(text, NumberStyles.Integer, CultureInfo.InvariantCulture, out value) && value > 0;
        }

        private double CalcDiscreteMean(double[] x, double[] p)
        {
            double mean = 0.0;
            for (int i = 0; i < x.Length; i++)
                mean += x[i] * p[i];
            return mean;
        }

        private double CalcDiscreteVariance(double[] x, double[] p, double mean)
        {
            double v = 0.0;
            for (int i = 0; i < x.Length; i++)
                v += p[i] * Math.Pow(x[i] - mean, 2);
            return v;
        }

        private class DiscreteResult
        {
            public double SampleMean;
            public double SampleVariance;
            public double ChiSquare;
            public double Critical;
            public double[] EmpiricalProbabilities;
        }

        private class NormalResult
        {
            public double[] Samples;
            public double SampleMean;
            public double SampleVariance;
            public double ChiSquare;
            public double Critical;
        }
    }
}