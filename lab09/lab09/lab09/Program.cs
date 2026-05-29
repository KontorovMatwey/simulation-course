using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Linq;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace PoissonCafeLab
{
    public class CustomerPlan
    {
        public int Id { get; set; }
        public double ArrivalTime { get; set; }
        public double ServiceStart { get; set; }
        public double ServiceEnd { get; set; }
        public double WaitingTime { get; set; }
    }

    static class Program
    {
        [STAThread]
        static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new MainForm());
        }
    }

    public class MainForm : Form
    {
        private NumericUpDown nudLambda;
        private NumericUpDown nudT;
        private NumericUpDown nudTrials;
        private NumericUpDown nudServiceMean;
        private Button btnRun;
        private Button btnAnimate;
        private Label lblStats;
        private Label lblCurrentModel;
        private Chart chartSystem;
        private Chart chartWaiting;
        private CafeSceneControl scene;
        private Timer animTimer;

        private List<CustomerPlan> demoPlans;
        private List<int> systemSamplesAll;
        private List<double> waitingSamplesAll;
        private List<SystemDistributionItem> systemDistribution;
        private List<WaitingBinItem> waitingDistribution;

        private double currentTime;
        private double demoHorizon;
        private bool isAnimating;

        public MainForm()
        {
            Text = "Лабораторная №9 — ПВЗ / очередь / оператор";
            Width = 1600;
            Height = 980;
            StartPosition = FormStartPosition.CenterScreen;
            Font = new Font("Segoe UI", 10);

            demoPlans = new List<CustomerPlan>();
            systemSamplesAll = new List<int>();
            waitingSamplesAll = new List<double>();
            systemDistribution = new List<SystemDistributionItem>();
            waitingDistribution = new List<WaitingBinItem>();

            InitializeControls();
            BuildUi();

            animTimer = new Timer();
            animTimer.Interval = 55;
            animTimer.Tick += AnimTimer_Tick;

            RunSimulation();
        }

        private void InitializeControls()
        {
            nudLambda = new NumericUpDown();
            nudT = new NumericUpDown();
            nudTrials = new NumericUpDown();
            nudServiceMean = new NumericUpDown();
            btnRun = new Button();
            btnAnimate = new Button();
            lblStats = new Label();
            lblCurrentModel = new Label();
            chartSystem = new Chart();
            chartWaiting = new Chart();
            scene = new CafeSceneControl();
        }

        private void BuildUi()
        {
            TableLayoutPanel root = new TableLayoutPanel();
            root.Dock = DockStyle.Fill;
            root.ColumnCount = 2;
            root.RowCount = 1;
            root.Padding = new Padding(10);
            root.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 420f));
            root.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100f));
            root.RowStyles.Add(new RowStyle(SizeType.Percent, 100f));
            Controls.Add(root);

            BuildLeftPanel(root);
            BuildRightPanel(root);
        }

        private void BuildLeftPanel(TableLayoutPanel root)
        {
            Panel leftPanel = new Panel();
            leftPanel.Dock = DockStyle.Fill;
            leftPanel.AutoScroll = true;
            leftPanel.Padding = new Padding(0, 0, 10, 0);
            root.Controls.Add(leftPanel, 0, 0);

            TableLayoutPanel left = new TableLayoutPanel();
            left.Dock = DockStyle.Top;
            left.AutoSize = true;
            left.ColumnCount = 1;
            left.RowCount = 4;
            left.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100f));
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 240f));
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 150f));
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 160f));
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 60f));
            leftPanel.Controls.Add(left);

            GroupBox groupParams = new GroupBox();
            groupParams.Text = "Параметры модели";
            groupParams.Dock = DockStyle.Fill;
            left.Controls.Add(groupParams, 0, 0);

            TableLayoutPanel paramsLayout = new TableLayoutPanel();
            paramsLayout.Dock = DockStyle.Fill;
            paramsLayout.ColumnCount = 2;
            paramsLayout.RowCount = 5;
            paramsLayout.Padding = new Padding(10);
            paramsLayout.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 60f));
            paramsLayout.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 40f));
            groupParams.Controls.Add(paramsLayout);

            AddParam(paramsLayout, "Интенсивность λ:", nudLambda, 0, 0.1m, 10m, 0.70m, 0.05m);
            AddParam(paramsLayout, "Интервал T:", nudT, 1, 0.1m, 1000m, 20.0m, 0.1m);
            AddParam(paramsLayout, "Число прогонов:", nudTrials, 2, 100m, 100000m, 10000m, 100m);
            AddParam(paramsLayout, "Среднее обслуживание:", nudServiceMean, 3, 0.05m, 10m, 0.80m, 0.05m);

            btnRun.Text = "Смоделировать";
            btnRun.Dock = DockStyle.Top;
            btnRun.Height = 38;
            btnRun.Click += delegate { RunSimulation(); };
            paramsLayout.Controls.Add(btnRun, 0, 4);
            paramsLayout.SetColumnSpan(btnRun, 2);

            GroupBox groupStats = new GroupBox();
            groupStats.Text = "Данные расчётов";
            groupStats.Dock = DockStyle.Fill;
            left.Controls.Add(groupStats, 0, 1);

            lblStats.Dock = DockStyle.Fill;
            lblStats.AutoSize = false;
            lblStats.Padding = new Padding(10);
            lblStats.TextAlign = ContentAlignment.TopLeft;
            groupStats.Controls.Add(lblStats);

            GroupBox groupCurrent = new GroupBox();
            groupCurrent.Text = "Текущая модель";
            groupCurrent.Dock = DockStyle.Fill;
            left.Controls.Add(groupCurrent, 0, 2);

            lblCurrentModel.Dock = DockStyle.Fill;
            lblCurrentModel.AutoSize = false;
            lblCurrentModel.Padding = new Padding(10);
            lblCurrentModel.TextAlign = ContentAlignment.TopLeft;
            groupCurrent.Controls.Add(lblCurrentModel);

            btnAnimate.Text = "Анимация заново";
            btnAnimate.Dock = DockStyle.Fill;
            btnAnimate.Height = 36;
            btnAnimate.Click += delegate { StartAnimation(); };
            left.Controls.Add(btnAnimate, 0, 3);
        }

        private void BuildRightPanel(TableLayoutPanel root)
        {
            TableLayoutPanel right = new TableLayoutPanel();
            right.Dock = DockStyle.Fill;
            right.ColumnCount = 1;
            right.RowCount = 3;
            right.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100f));
            right.RowStyles.Add(new RowStyle(SizeType.Absolute, 250f));
            right.RowStyles.Add(new RowStyle(SizeType.Percent, 50f));
            right.RowStyles.Add(new RowStyle(SizeType.Percent, 50f));
            root.Controls.Add(right, 1, 0);

            GroupBox sceneBox = new GroupBox();
            sceneBox.Text = "ПВЗ: вход, очередь и оператор";
            sceneBox.Dock = DockStyle.Fill;
            right.Controls.Add(sceneBox, 0, 0);

            scene.Dock = DockStyle.Fill;
            sceneBox.Controls.Add(scene);

            GroupBox systemBox = new GroupBox();
            systemBox.Text = "Эмпирическое распределение числа клиентов в системе";
            systemBox.Dock = DockStyle.Fill;
            right.Controls.Add(systemBox, 0, 1);

            chartSystem.Dock = DockStyle.Fill;
            chartSystem.ChartAreas.Clear();
            chartSystem.Series.Clear();
            chartSystem.Titles.Clear();
            chartSystem.Legends.Clear();
            chartSystem.ChartAreas.Add(CreateChartArea("System", "N", "P"));
            chartSystem.Legends.Add(new Legend("LegendSystem"));
            chartSystem.Titles.Add("Число клиентов в системе");
            systemBox.Controls.Add(chartSystem);

            GroupBox waitingBox = new GroupBox();
            waitingBox.Text = "Эмпирическое распределение времени ожидания в очереди";
            waitingBox.Dock = DockStyle.Fill;
            right.Controls.Add(waitingBox, 0, 2);

            chartWaiting.Dock = DockStyle.Fill;
            chartWaiting.ChartAreas.Clear();
            chartWaiting.Series.Clear();
            chartWaiting.Titles.Clear();
            chartWaiting.Legends.Clear();
            chartWaiting.ChartAreas.Add(CreateChartArea("Waiting", "Время ожидания", "P"));
            chartWaiting.Legends.Add(new Legend("LegendWaiting"));
            chartWaiting.Titles.Add("Время ожидания в очереди");
            waitingBox.Controls.Add(chartWaiting);
        }

        private ChartArea CreateChartArea(string name, string xTitle, string yTitle)
        {
            ChartArea area = new ChartArea(name);
            area.AxisX.Title = xTitle;
            area.AxisY.Title = yTitle;
            area.AxisX.IsMarginVisible = false;
            area.AxisY.IsMarginVisible = false;
            area.AxisX.MajorGrid.LineColor = Color.LightGray;
            area.AxisY.MajorGrid.LineColor = Color.LightGray;
            area.AxisX.IntervalAutoMode = IntervalAutoMode.VariableCount;
            area.AxisY.IntervalAutoMode = IntervalAutoMode.VariableCount;
            return area;
        }

        private void AddParam(TableLayoutPanel layout, string labelText, NumericUpDown nud, int row,
            decimal min, decimal max, decimal value, decimal increment)
        {
            Label lbl = new Label();
            lbl.Text = labelText;
            lbl.Dock = DockStyle.Fill;
            lbl.TextAlign = ContentAlignment.MiddleLeft;

            nud.Minimum = min;
            nud.Maximum = max;
            nud.Value = value;
            nud.DecimalPlaces = 2;
            nud.Increment = increment;
            nud.Dock = DockStyle.Fill;

            layout.RowStyles.Add(new RowStyle(SizeType.Absolute, 34f));
            layout.Controls.Add(lbl, 0, row);
            layout.Controls.Add(nud, 1, row);
        }

        private void RunSimulation()
        {
            double lambda = (double)nudLambda.Value;
            double T = (double)nudT.Value;
            int trials = (int)nudTrials.Value;
            double meanService = (double)nudServiceMean.Value;

            Random rnd = new Random(); // 42

            systemSamplesAll = new List<int>();
            waitingSamplesAll = new List<double>();

            int samplePointsPerTrial = 120;

            for (int i = 0; i < trials; i++)
            {
                TrialResult result = SimulateTrial(lambda, T, meanService, rnd);

                for (int j = 0; j <= samplePointsPerTrial; j++)
                {
                    double t = T * j / samplePointsPerTrial;
                    systemSamplesAll.Add(CountCustomersInSystem(result.Customers, t));
                }

                for (int j = 0; j < result.Customers.Count; j++)
                {
                    waitingSamplesAll.Add(result.Customers[j].WaitingTime);
                }
            }

            demoPlans = BuildDemoTrial(lambda, T, meanService, new Random(2025));
            demoHorizon = T;
            currentTime = 0.0;

            systemDistribution = BuildSystemDistribution(systemSamplesAll);
            waitingDistribution = BuildWaitingDistribution(waitingSamplesAll);

            DrawSystemChart();
            DrawWaitingChart();
            StartAnimation();

            double meanSystem = systemSamplesAll.Count > 0 ? systemSamplesAll.Average() : 0.0;
            double varSystem = PopulationVariance(systemSamplesAll.Select(x => (double)x).ToList(), meanSystem);

            double meanWaiting = waitingSamplesAll.Count > 0 ? waitingSamplesAll.Average() : 0.0;
            double varWaiting = PopulationVariance(waitingSamplesAll, meanWaiting);

            double rho = lambda * meanService;

            lblStats.Text =
                "λ = " + lambda.ToString("F2") + Environment.NewLine +
                "T = " + T.ToString("F2") + Environment.NewLine +
                "Среднее обслуживание = " + meanService.ToString("F2") + Environment.NewLine +
                "Нагрузка ρ = λ · tср = " + rho.ToString("F3") + Environment.NewLine +
                "Прогонов = " + trials.ToString() + Environment.NewLine + Environment.NewLine +
                "Среднее число клиентов в системе = " + meanSystem.ToString("F4") + Environment.NewLine +
                "Дисперсия числа клиентов в системе = " + varSystem.ToString("F4") + Environment.NewLine +
                "Среднее время ожидания = " + meanWaiting.ToString("F4") + Environment.NewLine +
                "Дисперсия времени ожидания = " + varWaiting.ToString("F4");

            UpdateCurrentModelStats();
        }

        private TrialResult SimulateTrial(double lambda, double T, double meanService, Random rnd)
        {
            List<double> arrivals = GenerateArrivals(lambda, T, rnd);
            List<CustomerPlan> customers = new List<CustomerPlan>();

            double serverFreeAt = 0.0;

            for (int i = 0; i < arrivals.Count; i++)
            {
                double arrival = arrivals[i];
                double serviceStart = Math.Max(arrival, serverFreeAt);
                double serviceTime = SampleExponential(meanService, rnd);
                double serviceEnd = serviceStart + serviceTime;
                double waitingTime = serviceStart - arrival;

                customers.Add(new CustomerPlan
                {
                    Id = i + 1,
                    ArrivalTime = arrival,
                    ServiceStart = serviceStart,
                    ServiceEnd = serviceEnd,
                    WaitingTime = waitingTime
                });

                serverFreeAt = serviceEnd;
            }

            return new TrialResult
            {
                Customers = customers
            };
        }

        private List<CustomerPlan> BuildDemoTrial(double lambda, double T, double meanService, Random rnd)
        {
            return SimulateTrial(lambda, T, meanService, rnd).Customers;
        }

        private List<double> GenerateArrivals(double lambda, double T, Random rnd)
        {
            List<double> arrivals = new List<double>();
            double time = 0.0;

            while (true)
            {
                double u = 1.0 - rnd.NextDouble();
                double dt = -Math.Log(u) / lambda;
                time += dt;

                if (time > T)
                    break;

                arrivals.Add(time);
            }

            return arrivals;
        }

        private double SampleExponential(double mean, Random rnd)
        {
            double u = 1.0 - rnd.NextDouble();
            return -mean * Math.Log(u);
        }

        private int CountCustomersInSystem(List<CustomerPlan> customers, double t)
        {
            int count = 0;
            for (int i = 0; i < customers.Count; i++)
            {
                if (customers[i].ArrivalTime <= t && customers[i].ServiceEnd > t)
                    count++;
            }
            return count;
        }

        private List<SystemDistributionItem> BuildSystemDistribution(List<int> samples)
        {
            List<SystemDistributionItem> result = new List<SystemDistributionItem>();

            if (samples == null || samples.Count == 0)
            {
                result.Add(new SystemDistributionItem { K = 0, Probability = 1.0 });
                return result;
            }

            int total = samples.Count;
            var groups = samples.GroupBy(x => x).OrderBy(g => g.Key);

            foreach (var g in groups)
            {
                result.Add(new SystemDistributionItem
                {
                    K = g.Key,
                    Probability = (double)g.Count() / total
                });
            }

            return result;
        }

        private List<WaitingBinItem> BuildWaitingDistribution(List<double> samples)
        {
            List<WaitingBinItem> result = new List<WaitingBinItem>();

            if (samples == null || samples.Count == 0)
            {
                result.Add(new WaitingBinItem { Left = 0, Right = 1, Center = 0.5, Probability = 1.0 });
                return result;
            }

            double min = samples.Min();
            double max = samples.Max();
            int n = samples.Count;

            int binCount = Math.Max(8, (int)Math.Round(Math.Sqrt(n)));
            binCount = Math.Min(binCount, 30);

            if (Math.Abs(max - min) < 1e-12)
            {
                result.Add(new WaitingBinItem
                {
                    Left = Math.Max(0, min - 0.5),
                    Right = min + 0.5,
                    Center = min,
                    Probability = 1.0
                });
                return result;
            }

            double width = (max - min) / binCount;
            if (width <= 0)
                width = 1e-6;

            double[] freq = new double[binCount];

            for (int i = 0; i < samples.Count; i++)
            {
                int index = (int)((samples[i] - min) / width);
                if (index < 0) index = 0;
                if (index >= binCount) index = binCount - 1;
                freq[index]++;
            }

            for (int i = 0; i < binCount; i++)
            {
                double left = min + i * width;
                double right = left + width;
                result.Add(new WaitingBinItem
                {
                    Left = left,
                    Right = right,
                    Center = (left + right) / 2.0,
                    Probability = freq[i] / n
                });
            }

            return result;
        }

        private void DrawSystemChart()
        {
            chartSystem.Series.Clear();

            ChartArea area = chartSystem.ChartAreas[0];
            area.AxisX.Minimum = systemDistribution.Count > 0 ? systemDistribution.Min(x => x.K) - 1 : 0;
            area.AxisX.Maximum = systemDistribution.Count > 0 ? systemDistribution.Max(x => x.K) + 1 : 1;
            area.AxisY.Minimum = 0;
            area.AxisY.Maximum = systemDistribution.Count > 0
                ? Math.Max(0.1, systemDistribution.Max(x => x.Probability) * 1.25)
                : 1;

            Series bars = new Series("P(N=k)");
            bars.ChartType = SeriesChartType.Column;
            bars.IsValueShownAsLabel = false;
            bars.Color = Color.SteelBlue;

            Series polygon = new Series("Полигон частот");
            polygon.ChartType = SeriesChartType.Line;
            polygon.BorderWidth = 2;
            polygon.MarkerStyle = MarkerStyle.Circle;
            polygon.MarkerSize = 6;
            polygon.Color = Color.DarkOrange;

            for (int i = 0; i < systemDistribution.Count; i++)
            {
                bars.Points.AddXY(systemDistribution[i].K, systemDistribution[i].Probability);
                polygon.Points.AddXY(systemDistribution[i].K, systemDistribution[i].Probability);
            }

            chartSystem.Series.Add(bars);
            chartSystem.Series.Add(polygon);
            chartSystem.Legends[0].Docking = Docking.Top;
        }

        private void DrawWaitingChart()
        {
            chartWaiting.Series.Clear();

            ChartArea area = chartWaiting.ChartAreas[0];

            if (waitingDistribution.Count > 0)
            {
                area.AxisX.Minimum = waitingDistribution.Min(x => x.Left);
                area.AxisX.Maximum = waitingDistribution.Max(x => x.Right);
                area.AxisY.Minimum = 0;
                area.AxisY.Maximum = Math.Max(0.1, waitingDistribution.Max(x => x.Probability) * 1.25);
            }
            else
            {
                area.AxisX.Minimum = 0;
                area.AxisX.Maximum = 1;
                area.AxisY.Minimum = 0;
                area.AxisY.Maximum = 1;
            }

            Series bars = new Series("Гистограмма");
            bars.ChartType = SeriesChartType.Column;
            bars.IsValueShownAsLabel = false;
            bars.Color = Color.MediumSeaGreen;

            Series polygon = new Series("Полигон частот");
            polygon.ChartType = SeriesChartType.Line;
            polygon.BorderWidth = 2;
            polygon.MarkerStyle = MarkerStyle.Circle;
            polygon.MarkerSize = 5;
            polygon.Color = Color.Crimson;

            for (int i = 0; i < waitingDistribution.Count; i++)
            {
                bars.Points.AddXY(waitingDistribution[i].Center, waitingDistribution[i].Probability);
                polygon.Points.AddXY(waitingDistribution[i].Center, waitingDistribution[i].Probability);
            }

            chartWaiting.Series.Add(bars);
            chartWaiting.Series.Add(polygon);
            chartWaiting.Legends[0].Docking = Docking.Top;
        }

        private void StartAnimation()
        {
            currentTime = 0.0;
            isAnimating = true;
            animTimer.Start();
            UpdateScene();
        }

        private void AnimTimer_Tick(object sender, EventArgs e)
        {
            if (!isAnimating)
                return;

            double step = demoHorizon / 420.0;
            if (step <= 0)
                step = 0.03;

            currentTime += step;

            if (currentTime >= demoHorizon)
            {
                currentTime = demoHorizon;
                isAnimating = false;
                animTimer.Stop();
            }

            UpdateScene();
        }

        private void UpdateScene()
        {
            scene.SetData(demoPlans, currentTime, demoHorizon);
            UpdateCurrentModelStats();
        }

        private void UpdateCurrentModelStats()
        {
            if (demoPlans == null || demoPlans.Count == 0)
            {
                lblCurrentModel.Text = "Нет данных.";
                return;
            }

            int arrived = GetArrivedCount();
            int served = GetServedCount();
            int inService = GetInServiceCount();
            int inQueue = Math.Max(0, arrived - served - inService);
            int visible = GetVisibleCount();

            double avgWaitingServed = GetAverageWaitingOfServed();

            lblCurrentModel.Text =
                "Время t = " + currentTime.ToString("F2") + Environment.NewLine +
                "Клиентов всего в демонстрации = " + demoPlans.Count.ToString() + Environment.NewLine +
                "Пришло к текущему моменту = " + arrived.ToString() + Environment.NewLine +
                "Обслужено = " + served.ToString() + Environment.NewLine +
                "В очереди = " + inQueue.ToString() + Environment.NewLine +
                "У оператора = " + inService.ToString() + Environment.NewLine +
                "В системе сейчас = " + visible.ToString() + Environment.NewLine +
                "Среднее ожидание обслуженных = " + avgWaitingServed.ToString("F4");
        }

        private int GetArrivedCount()
        {
            int count = 0;
            for (int i = 0; i < demoPlans.Count; i++)
            {
                if (demoPlans[i].ArrivalTime <= currentTime)
                    count++;
            }
            return count;
        }

        private int GetServedCount()
        {
            int count = 0;
            for (int i = 0; i < demoPlans.Count; i++)
            {
                if (demoPlans[i].ServiceEnd <= currentTime)
                    count++;
            }
            return count;
        }

        private int GetInServiceCount()
        {
            int count = 0;
            for (int i = 0; i < demoPlans.Count; i++)
            {
                if (demoPlans[i].ServiceStart <= currentTime && demoPlans[i].ServiceEnd > currentTime)
                    count++;
            }
            return count;
        }

        private int GetVisibleCount()
        {
            int count = 0;
            for (int i = 0; i < demoPlans.Count; i++)
            {
                if (demoPlans[i].ArrivalTime <= currentTime && demoPlans[i].ServiceEnd > currentTime)
                    count++;
            }
            return count;
        }

        private double GetAverageWaitingOfServed()
        {
            double sum = 0.0;
            int count = 0;

            for (int i = 0; i < demoPlans.Count; i++)
            {
                if (demoPlans[i].ServiceEnd <= currentTime)
                {
                    sum += demoPlans[i].WaitingTime;
                    count++;
                }
            }

            if (count == 0)
                return 0.0;

            return sum / count;
        }

        private double PopulationVariance(List<double> values, double mean)
        {
            if (values == null || values.Count == 0)
                return 0.0;

            double sum = 0.0;
            for (int i = 0; i < values.Count; i++)
            {
                double d = values[i] - mean;
                sum += d * d;
            }
            return sum / values.Count;
        }

        private class TrialResult
        {
            public List<CustomerPlan> Customers { get; set; }
        }

        private class SystemDistributionItem
        {
            public int K { get; set; }
            public double Probability { get; set; }
        }

        private class WaitingBinItem
        {
            public double Left { get; set; }
            public double Right { get; set; }
            public double Center { get; set; }
            public double Probability { get; set; }
        }
    }

    public class CafeSceneControl : Control
    {
        private List<CustomerPlan> plans;
        private double currentTime;
        private double totalTime;

        public CafeSceneControl()
        {
            plans = new List<CustomerPlan>();
            DoubleBuffered = true;
            BackColor = Color.White;
            SetStyle(ControlStyles.AllPaintingInWmPaint | ControlStyles.UserPaint | ControlStyles.OptimizedDoubleBuffer, true);
        }

        public void SetData(List<CustomerPlan> newPlans, double time, double duration)
        {
            plans = newPlans ?? new List<CustomerPlan>();
            currentTime = time;
            totalTime = duration;
            Invalidate();
        }

        protected override void OnPaint(PaintEventArgs e)
        {
            base.OnPaint(e);

            Graphics g = e.Graphics;
            g.SmoothingMode = SmoothingMode.AntiAlias;

            int w = Width;
            int h = Height;

            g.Clear(Color.FromArgb(250, 246, 238));

            using (SolidBrush floor = new SolidBrush(Color.FromArgb(238, 230, 216)))
                g.FillRectangle(floor, 0, (int)(h * 0.68), w, (int)(h * 0.32));

            DrawSign(g, w);
            DrawDoor(g, w, h);
            DrawQueueLine(g, w, h);
            DrawCounter(g, w, h);
            DrawOperator(g, w, h);
            DrawCustomers(g, w, h);

            using (Font f = new Font("Segoe UI", 10, FontStyle.Bold))
            using (Brush b = new SolidBrush(Color.FromArgb(70, 70, 70)))
            {
                g.DrawString("Время: " + currentTime.ToString("F2") + " / " + totalTime.ToString("F2"), f, b, 16, 12);
            }
        }

        private void DrawSign(Graphics g, int w)
        {
            Rectangle sign = new Rectangle(w / 2 - 92, 18, 184, 40);
            using (SolidBrush br = new SolidBrush(Color.FromArgb(95, 60, 31)))
                g.FillRectangle(br, sign);
            using (Pen p = new Pen(Color.FromArgb(60, 35, 20), 2))
                g.DrawRectangle(p, sign);

            using (Font f = new Font("Segoe UI", 14, FontStyle.Bold))
            using (Brush b = new SolidBrush(Color.White))
            {
                StringFormat sf = new StringFormat();
                sf.Alignment = StringAlignment.Center;
                sf.LineAlignment = StringAlignment.Center;
                g.DrawString("ПВЗ", f, b, sign, sf);
            }
        }

        private void DrawDoor(Graphics g, int w, int h)
        {
            RectangleF door = new RectangleF(28, h * 0.42f, 62, 126);
            using (SolidBrush br = new SolidBrush(Color.FromArgb(215, 199, 181)))
                g.FillRectangle(br, door);
            using (Pen p = new Pen(Color.FromArgb(120, 100, 80), 2))
                g.DrawRectangle(p, door.X, door.Y, door.Width, door.Height);

            using (Pen p = new Pen(Color.FromArgb(120, 100, 80), 2))
            {
                g.DrawLine(p, door.X + 8, door.Y, door.X + 8, door.Y + door.Height);
                g.DrawLine(p, door.X + 31, door.Y, door.X + 31, door.Y + door.Height);
            }

            using (Brush b = new SolidBrush(Color.FromArgb(60, 60, 60)))
                g.FillEllipse(b, door.X + 45, door.Y + 56, 6, 6);
        }

        private void DrawQueueLine(Graphics g, int w, int h)
        {
            float y = h * 0.67f;
            float startX = w * 0.12f;
            float endX = w * 0.68f;

            using (Pen p = new Pen(Color.FromArgb(160, 160, 160), 1))
            {
                p.DashStyle = DashStyle.Dot;
                g.DrawLine(p, startX, y, endX, y);

                int slots = 8;
                float spacing = (endX - startX) / slots;
                for (int i = 0; i <= slots; i++)
                {
                    float x = startX + i * spacing;
                    g.DrawEllipse(p, x - 11, y - 11, 22, 22);
                }
            }
        }

        private void DrawCounter(Graphics g, int w, int h)
        {
            float x = w * 0.70f;
            float y = h * 0.40f;
            float cw = w * 0.24f;
            float ch = h * 0.18f;

            RectangleF counter = new RectangleF(x, y, cw, ch);
            using (SolidBrush br = new SolidBrush(Color.FromArgb(175, 119, 68)))
                g.FillRectangle(br, counter);
            using (Pen p = new Pen(Color.FromArgb(120, 76, 38), 2))
                g.DrawRectangle(p, counter.X, counter.Y, counter.Width, counter.Height);

            RectangleF top = new RectangleF(x, y - 12, cw, 14);
            using (SolidBrush br = new SolidBrush(Color.FromArgb(209, 153, 98)))
                g.FillRectangle(br, top);
            using (Pen p = new Pen(Color.FromArgb(120, 76, 38), 2))
                g.DrawRectangle(p, top.X, top.Y, top.Width, top.Height);
        }

        private void DrawOperator(Graphics g, int w, int h)
        {
            float x = w * 0.80f;
            float y = h * 0.29f;

            using (Brush skin = new SolidBrush(Color.FromArgb(244, 207, 182)))
                g.FillEllipse(skin, x, y, 36, 36);

            using (Brush hair = new SolidBrush(Color.FromArgb(90, 60, 30)))
                g.FillEllipse(hair, x - 1, y - 2, 38, 18);

            using (Brush body = new SolidBrush(Color.FromArgb(79, 134, 171)))
                g.FillRectangle(body, x - 4, y + 32, 44, 58);

            using (Pen outline = new Pen(Color.FromArgb(60, 60, 60), 2))
            {
                g.DrawEllipse(outline, x, y, 36, 36);
                g.DrawRectangle(outline, x - 4, y + 32, 44, 58);
            }
        }

        private void DrawCustomers(Graphics g, int w, int h)
        {
            if (plans == null || plans.Count == 0)
                return;

            float entranceX = 34f;
            float entranceY = h * 0.57f;

            float queueStartX = w * 0.18f;
            float queueEndX = w * 0.68f;
            float queueY = h * 0.64f;
            int maxVisible = 8;
            float queueStep = (queueEndX - queueStartX) / (maxVisible - 1);

            float counterX = w * 0.74f;
            float counterY = h * 0.53f;

            float exitX = w - 42f;
            float exitY = h * 0.56f;

            List<CustomerPlan> waiting = plans
                .Where(p => currentTime >= p.ArrivalTime && currentTime < p.ServiceStart)
                .OrderBy(p => p.ArrivalTime)
                .ToList();

            CustomerPlan serving = null;
            for (int i = 0; i < plans.Count; i++)
            {
                if (currentTime >= plans[i].ServiceStart && currentTime < plans[i].ServiceEnd)
                {
                    serving = plans[i];
                    break;
                }
            }

            List<CustomerPlan> leaving = plans
                .Where(p => currentTime >= p.ServiceEnd && currentTime < p.ServiceEnd + 0.5)
                .ToList();

            int hiddenCount = 0;

            for (int i = 0; i < waiting.Count; i++)
            {
                if (i >= maxVisible)
                {
                    hiddenCount = waiting.Count - maxVisible;
                    break;
                }

                CustomerPlan p = waiting[i];
                float targetX = queueEndX - i * queueStep;
                float targetY = queueY;

                float x = targetX;
                float y = targetY - 22f;

                double move = 0.5;
                if (currentTime < p.ArrivalTime + move)
                {
                    float k = (float)((currentTime - p.ArrivalTime) / move);
                    if (k < 0f) k = 0f;
                    if (k > 1f) k = 1f;

                    x = Lerp(entranceX, targetX, k);
                    y = Lerp(entranceY, targetY - 22f, k);
                }

                DrawPerson(g, x, y, p.Id, Color.FromArgb(106, 168, 79));
            }

            if (hiddenCount > 0)
            {
                using (Font f = new Font("Segoe UI", 11, FontStyle.Bold))
                using (Brush b = new SolidBrush(Color.FromArgb(70, 70, 70)))
                {
                    g.DrawString("+" + hiddenCount.ToString(), f, b, queueStartX - 8, queueY - 24);
                }
            }

            if (serving != null)
            {
                float x = counterX;
                float y = counterY;

                if (currentTime < serving.ServiceStart + 0.35)
                {
                    float k = (float)((currentTime - serving.ServiceStart) / 0.35);
                    if (k < 0f) k = 0f;
                    if (k > 1f) k = 1f;

                    float fromX = serving.WaitingTime > 0.0001
                        ? queueEndX
                        : entranceX;

                    float fromY = serving.WaitingTime > 0.0001
                        ? queueY - 22f
                        : entranceY;

                    x = Lerp(fromX, counterX, k);
                    y = Lerp(fromY, counterY, k);
                }

                DrawPerson(g, x, y, serving.Id, Color.FromArgb(244, 162, 97));
            }

            for (int i = 0; i < leaving.Count; i++)
            {
                CustomerPlan p = leaving[i];
                double move = 0.45;
                float k = (float)((currentTime - p.ServiceEnd) / move);
                if (k < 0f) k = 0f;
                if (k > 1f) k = 1f;

                float x = Lerp(counterX + 18f, exitX, k);
                float y = Lerp(counterY, exitY, k);

                DrawPerson(g, x, y, p.Id, Color.FromArgb(153, 102, 204));
            }
        }

        private void DrawPerson(Graphics g, float x, float y, int id, Color bodyColor)
        {
            using (Brush skin = new SolidBrush(Color.FromArgb(245, 220, 193)))
                g.FillEllipse(skin, x, y, 20, 20);

            using (Brush body = new SolidBrush(bodyColor))
                g.FillEllipse(body, x - 2, y + 18, 24, 24);

            using (Pen outline = new Pen(Color.FromArgb(70, 70, 70), 1.5f))
            {
                g.DrawEllipse(outline, x, y, 20, 20);
                g.DrawEllipse(outline, x - 2, y + 18, 24, 24);
            }

            using (Font f = new Font("Segoe UI", 7, FontStyle.Bold))
            using (Brush b = new SolidBrush(Color.White))
                g.DrawString(id.ToString(), f, b, x + 5, y + 24);
        }

        private float Lerp(float a, float b, float t)
        {
            return a + (b - a) * t;
        }
    }
}