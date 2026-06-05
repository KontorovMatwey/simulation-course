using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Text;
using System.Windows.Forms;

namespace BankAgentLab
{
    public enum CustomerState
    {
        Waiting,
        InService,
        Served,
        Left
    }

    public class CustomerAgent
    {
        public int Id { get; set; }
        public double ArrivalTime { get; set; }
        public double ServiceTime { get; set; }
        public double ServiceStart { get; set; }
        public double ServiceEnd { get; set; }
        public double WaitingTime { get; set; }
        public double Waited { get; set; }
        public double Impatience { get; set; }
        public CustomerState State { get; set; }
    }

    public class CashierAgent
    {
        public int Id { get; set; }
        public CustomerAgent Current { get; set; }
        public double FreeAt { get; set; }
        public double BusyTime { get; set; }

        public bool Busy
        {
            get { return Current != null; }
        }

        public double BusyShare(double time)
        {
            if (time <= 0.0) return 0.0;
            return BusyTime / time;
        }
    }

    public class BankModel
    {
        public List<CustomerAgent> Customers { get; private set; }
        public List<CustomerAgent> Waiting { get; private set; }
        public List<CashierAgent> Cashiers { get; private set; }

        public double Time { get; private set; }
        public double ArrivalHorizon { get; private set; }
        public double Lambda { get; private set; }
        public double MeanService { get; private set; }

        public int ServedCount { get; private set; }
        public int LeftCount { get; private set; }

        private double nextArrival;
        private double totalBusyTime;
        private int nextId;
        private Random rnd;

        public BankModel(int cashierCount, double lambda, double meanService, double horizon, int seed)
        {
            Customers = new List<CustomerAgent>();
            Waiting = new List<CustomerAgent>();
            Cashiers = new List<CashierAgent>();

            Lambda = lambda;
            MeanService = meanService;
            ArrivalHorizon = horizon;

            rnd = new Random(seed);
            ResetCashiers(cashierCount);
            ResetSimulation();
        }

        private void ResetCashiers(int count)
        {
            Cashiers.Clear();
            for (int i = 0; i < count; i++)
            {
                Cashiers.Add(new CashierAgent { Id = i + 1 });
            }
        }

        private void ResetSimulation()
        {
            Time = 0.0;
            nextArrival = Exp(Lambda);
            totalBusyTime = 0.0;
            nextId = 1;
            ServedCount = 0;
            LeftCount = 0;

            Customers.Clear();
            Waiting.Clear();

            for (int i = 0; i < Cashiers.Count; i++)
            {
                Cashiers[i].Current = null;
                Cashiers[i].FreeAt = 0.0;
                Cashiers[i].BusyTime = 0.0;
            }
        }

        public void Restart(double lambda, double meanService, double horizon)
        {
            Lambda = lambda;
            MeanService = meanService;
            ArrivalHorizon = horizon;
            ResetSimulation();
        }

        public void Step(double dt)
        {
            Time += dt;

            // 1) Новые приходы
            while (nextArrival <= ArrivalHorizon && nextArrival <= Time)
            {
                AddArrival(nextArrival);
                nextArrival += Exp(Lambda);
            }

            // 2) Освобождение касс
            for (int i = 0; i < Cashiers.Count; i++)
            {
                var cashier = Cashiers[i];
                if (cashier.Busy && Time >= cashier.FreeAt)
                {
                    cashier.Current.State = CustomerState.Served;
                    cashier.Current.ServiceEnd = cashier.FreeAt;
                    ServedCount++;
                    cashier.Current = null;
                }
            }

            // 3) Назначение клиентов свободным кассам
            AssignWaitingCustomers();

            // 4) Нетерпеливость: клиенты в очереди могут уйти
            UpdateImpatience(dt);

            // 5) Учёт времени занятости
            for (int i = 0; i < Cashiers.Count; i++)
            {
                if (Cashiers[i].Busy)
                {
                    Cashiers[i].BusyTime += dt;
                    totalBusyTime += dt;
                }
            }
        }

        private void AddArrival(double arrivalTime)
        {
            var c = new CustomerAgent();
            c.Id = nextId++;
            c.ArrivalTime = arrivalTime;
            c.ServiceTime = Exp(1.0 / MeanService); // среднее время обслуживания
            c.Impatience = 0.01 + rnd.NextDouble() * 0.03;
            c.State = CustomerState.Waiting;

            Customers.Add(c);
            Waiting.Add(c);
        }

        private void AssignWaitingCustomers()
        {
            for (int i = 0; i < Cashiers.Count; i++)
            {
                var cashier = Cashiers[i];
                if (cashier.Busy) continue;
                if (Waiting.Count == 0) break;

                var c = Waiting[0];
                Waiting.RemoveAt(0);

                c.State = CustomerState.InService;
                c.ServiceStart = Time;
                c.WaitingTime = c.ServiceStart - c.ArrivalTime;

                cashier.Current = c;
                cashier.FreeAt = Time + c.ServiceTime;
            }
        }

        private void UpdateImpatience(double dt)
        {
            for (int i = Waiting.Count - 1; i >= 0; i--)
            {
                var c = Waiting[i];
                c.Waited += dt;

                int ahead = i; // 0 = самый близкий к кассе
                double leaveChance = c.Impatience + 0.01 * ahead + 0.002 * c.Waited;

                // Самый первый в очереди почти не уходит
                if (ahead == 0)
                    leaveChance *= 0.02;

                if (leaveChance > 0.95)
                    leaveChance = 0.95;

                if (rnd.NextDouble() < leaveChance * dt)
                {
                    c.State = CustomerState.Left;
                    LeftCount++;
                    Waiting.RemoveAt(i);
                }
            }
        }

        private double Exp(double rate)
        {
            double u = 1.0 - rnd.NextDouble();
            return -Math.Log(u) / rate;
        }

        public int BusyCount
        {
            get
            {
                int count = 0;
                for (int i = 0; i < Cashiers.Count; i++)
                    if (Cashiers[i].Busy) count++;
                return count;
            }
        }

        public int VisibleCount
        {
            get { return Waiting.Count + BusyCount; }
        }

        public double BusyProbability
        {
            get
            {
                if (Time <= 0.0 || Cashiers.Count == 0) return 0.0;
                return totalBusyTime / (Cashiers.Count * Time);
            }
        }

        public double FreeProbability
        {
            get { return 1.0 - BusyProbability; }
        }

        public bool IsFinished
        {
            get
            {
                if (Time < ArrivalHorizon) return false;
                if (Waiting.Count > 0) return false;
                for (int i = 0; i < Cashiers.Count; i++)
                    if (Cashiers[i].Busy) return false;
                return true;
            }
        }

        public string CashiersText()
        {
            var sb = new StringBuilder();
            for (int i = 0; i < Cashiers.Count; i++)
            {
                var c = Cashiers[i];
                sb.Append("Касса ").Append(c.Id).Append(": ");
                sb.Append(c.Busy ? "занята" : "свободна");
                if (c.Busy && c.Current != null)
                    sb.Append(" (клиент #").Append(c.Current.Id).Append(")");
                sb.AppendLine();
            }

            sb.AppendLine();
            sb.AppendLine("Доля занятости касс:");
            for (int i = 0; i < Cashiers.Count; i++)
            {
                sb.Append("К").Append(Cashiers[i].Id).Append(": ");
                sb.Append(Cashiers[i].BusyShare(Time).ToString("F4"));
                sb.AppendLine();
            }

            return sb.ToString();
        }

        public string LiveText()
        {
            return
                "t = " + Time.ToString("F2") + Environment.NewLine +
                "Пришло клиентов = " + Customers.Count.ToString() + Environment.NewLine +
                "Обслужено = " + ServedCount.ToString() + Environment.NewLine +
                "Ушло по нетерпеливости = " + LeftCount.ToString() + Environment.NewLine +
                "В очереди = " + Waiting.Count.ToString() + Environment.NewLine +
                "У касс = " + BusyCount.ToString() + Environment.NewLine +
                "В системе = " + VisibleCount.ToString();
        }
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
        private const int CashierCount = 3;

        private NumericUpDown nudLambda;
        private NumericUpDown nudServiceMean;
        private NumericUpDown nudT;
        private Button btnRun;
        private Button btnReplay;
        private Label lblProb;
        private Label lblLive;
        private Label lblCashiers;
        private BankSceneControl scene;
        private Timer timer;

        private BankModel model;

        public MainForm()
        {
            Text = "Банк — M/M/c с нетерпеливостью";
            Width = 1380;
            Height = 820;
            StartPosition = FormStartPosition.CenterScreen;
            Font = new Font("Segoe UI", 10);

            BuildUi();

            timer = new Timer();
            timer.Interval = 40;
            timer.Tick += Timer_Tick;

            RunModel();
        }

        private void BuildUi()
        {
            var root = new TableLayoutPanel();
            root.Dock = DockStyle.Fill;
            root.ColumnCount = 2;
            root.RowCount = 1;
            root.Padding = new Padding(10);
            root.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 360f));
            root.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100f));
            Controls.Add(root);

            BuildLeft(root);
            BuildRight(root);
        }

        private void BuildLeft(TableLayoutPanel root)
        {
            var left = new TableLayoutPanel();
            left.Dock = DockStyle.Fill;
            left.ColumnCount = 1;
            left.RowCount = 5;
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 210f));
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 120f));
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 140f));
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 120f));
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 40f));
            root.Controls.Add(left, 0, 0);

            var gbParams = new GroupBox();
            gbParams.Text = "Параметры";
            gbParams.Dock = DockStyle.Fill;
            left.Controls.Add(gbParams, 0, 0);

            var p = new TableLayoutPanel();
            p.Dock = DockStyle.Fill;
            p.ColumnCount = 2;
            p.RowCount = 4;
            p.Padding = new Padding(10);
            p.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 60f));
            p.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 40f));
            gbParams.Controls.Add(p);

            nudLambda = MakeNud(0.1m, 10m, 0.70m, 0.05m);
            nudServiceMean = MakeNud(0.05m, 10m, 0.80m, 0.05m);
            nudT = MakeNud(1m, 1000m, 20m, 1m);

            AddRow(p, 0, "Интенсивность λ:", nudLambda);
            AddRow(p, 1, "Среднее обслуживание:", nudServiceMean);
            AddRow(p, 2, "Время прихода T:", nudT);

            btnRun = new Button();
            btnRun.Text = "Пересчитать";
            btnRun.Dock = DockStyle.Fill;
            btnRun.Click += delegate { RunModel(); };
            p.Controls.Add(btnRun, 0, 3);
            p.SetColumnSpan(btnRun, 2);

            var gbProb = new GroupBox();
            gbProb.Text = "Вероятности прибора";
            gbProb.Dock = DockStyle.Fill;
            left.Controls.Add(gbProb, 0, 1);

            lblProb = MakeLabel();
            gbProb.Controls.Add(lblProb);

            var gbLive = new GroupBox();
            gbLive.Text = "Текущая модель";
            gbLive.Dock = DockStyle.Fill;
            left.Controls.Add(gbLive, 0, 2);

            lblLive = MakeLabel();
            gbLive.Controls.Add(lblLive);

            var gbCash = new GroupBox();
            gbCash.Text = "Кассы";
            gbCash.Dock = DockStyle.Fill;
            left.Controls.Add(gbCash, 0, 3);

            lblCashiers = MakeLabel();
            gbCash.Controls.Add(lblCashiers);

            btnReplay = new Button();
            btnReplay.Text = "Анимация заново";
            btnReplay.Dock = DockStyle.Fill;
            btnReplay.Click += delegate { RestartAnimation(); };
            left.Controls.Add(btnReplay, 0, 4);
        }

        private void BuildRight(TableLayoutPanel root)
        {
            var gb = new GroupBox();
            gb.Text = "Банк";
            gb.Dock = DockStyle.Fill;
            root.Controls.Add(gb, 1, 0);

            scene = new BankSceneControl();
            scene.Dock = DockStyle.Fill;
            gb.Controls.Add(scene);
        }

        private NumericUpDown MakeNud(decimal min, decimal max, decimal value, decimal inc)
        {
            var n = new NumericUpDown();
            n.Minimum = min;
            n.Maximum = max;
            n.Value = value;
            n.DecimalPlaces = 2;
            n.Increment = inc;
            n.Dock = DockStyle.Fill;
            return n;
        }

        private void AddRow(TableLayoutPanel p, int row, string text, Control control)
        {
            p.RowStyles.Add(new RowStyle(SizeType.Absolute, 34f));

            var lbl = new Label();
            lbl.Text = text;
            lbl.Dock = DockStyle.Fill;
            lbl.TextAlign = ContentAlignment.MiddleLeft;

            p.Controls.Add(lbl, 0, row);
            p.Controls.Add(control, 1, row);
        }

        private Label MakeLabel()
        {
            var lbl = new Label();
            lbl.Dock = DockStyle.Fill;
            lbl.Padding = new Padding(10);
            lbl.TextAlign = ContentAlignment.TopLeft;
            lbl.AutoSize = false;
            return lbl;
        }

        private void RunModel()
        {
            double lambda = (double)nudLambda.Value;
            double meanService = (double)nudServiceMean.Value;
            double horizon = (double)nudT.Value;

            model = new BankModel(CashierCount, lambda, meanService, horizon, 42);
            timer.Start();
            UpdateUi();
        }

        private void RestartAnimation()
        {
            if (model == null) return;

            double lambda = (double)nudLambda.Value;
            double meanService = (double)nudServiceMean.Value;
            double horizon = (double)nudT.Value;

            model.Restart(lambda, meanService, horizon);
            timer.Start();
            UpdateUi();
        }

        private void Timer_Tick(object sender, EventArgs e)
        {
            if (model == null) return;

            model.Step(0.03);
            scene.SetModel(model);
            UpdateUi();

            if (model.IsFinished)
                timer.Stop();
        }

        private void UpdateUi()
        {
            if (model == null) return;

            lblProb.Text =
                "P(занят) = " + model.BusyProbability.ToString("F4") + Environment.NewLine +
                "P(свободен) = " + model.FreeProbability.ToString("F4") + Environment.NewLine +
                "ρ = λ · tср = " + ((double)nudLambda.Value * (double)nudServiceMean.Value).ToString("F4");

            lblLive.Text = model.LiveText();
            lblCashiers.Text = model.CashiersText();
            scene.SetModel(model);
        }
    }

    public class BankSceneControl : Control
    {
        private BankModel model;

        public BankSceneControl()
        {
            DoubleBuffered = true;
            BackColor = Color.White;
            SetStyle(ControlStyles.AllPaintingInWmPaint | ControlStyles.UserPaint | ControlStyles.OptimizedDoubleBuffer, true);
        }

        public void SetModel(BankModel m)
        {
            model = m;
            Invalidate();
        }

        protected override void OnPaint(PaintEventArgs e)
        {
            base.OnPaint(e);

            Graphics g = e.Graphics;
            g.SmoothingMode = SmoothingMode.AntiAlias;

            int w = Width;
            int h = Height;

            g.Clear(Color.FromArgb(248, 244, 236));

            using (Brush floor = new SolidBrush(Color.FromArgb(236, 228, 216)))
                g.FillRectangle(floor, 0, (int)(h * 0.67), w, (int)(h * 0.33));

            DrawTitle(g, w);
            DrawDoor(g, h);
            DrawCashiers(g, w, h);
            DrawQueue(g, w, h);

            using (Font f = new Font("Segoe UI", 10, FontStyle.Bold))
            using (Brush b = new SolidBrush(Color.FromArgb(60, 60, 60)))
            {
                double t = model == null ? 0.0 : model.Time;
                double T = model == null ? 0.0 : model.ArrivalHorizon;
                g.DrawString("Время: " + t.ToString("F2") + " / " + T.ToString("F2"), f, b, 14, 12);
            }
        }

        private void DrawTitle(Graphics g, int w)
        {
            using (Font f = new Font("Segoe UI", 14, FontStyle.Bold))
            using (Brush b = new SolidBrush(Color.FromArgb(95, 60, 31)))
            {
                g.DrawString("БАНК", f, b, w / 2 - 40, 10);
            }
        }

        private void DrawDoor(Graphics g, int h)
        {
            RectangleF door = new RectangleF(22, h * 0.42f, 60, 120);
            using (Brush br = new SolidBrush(Color.FromArgb(220, 205, 188)))
                g.FillRectangle(br, door);
            using (Pen p = new Pen(Color.FromArgb(120, 100, 80), 2))
                g.DrawRectangle(p, door.X, door.Y, door.Width, door.Height);
        }

        private void DrawCashiers(Graphics g, int w, int h)
        {
            if (model == null) return;

            float[] xs = new float[] { w * 0.72f, w * 0.83f, w * 0.94f };
            float y = h * 0.38f;

            for (int i = 0; i < model.Cashiers.Count; i++)
            {
                float x = xs[i];

                RectangleF desk = new RectangleF(x - 26, y + 30, 52, 58);
                using (Brush br = new SolidBrush(Color.FromArgb(176, 124, 72)))
                    g.FillRectangle(br, desk);
                using (Pen p = new Pen(Color.FromArgb(120, 76, 38), 2))
                    g.DrawRectangle(p, desk.X, desk.Y, desk.Width, desk.Height);

                using (Brush head = new SolidBrush(Color.FromArgb(244, 207, 182)))
                    g.FillEllipse(head, x - 11, y, 22, 22);

                using (Brush hair = new SolidBrush(Color.FromArgb(95, 65, 35)))
                    g.FillEllipse(hair, x - 12, y - 1, 24, 10);

                using (Brush body = new SolidBrush(Color.FromArgb(79, 134, 171)))
                    g.FillRectangle(body, x - 13, y + 20, 26, 34);

                using (Pen p = new Pen(Color.FromArgb(60, 60, 60), 2))
                {
                    g.DrawEllipse(p, x - 11, y, 22, 22);
                    g.DrawRectangle(p, x - 13, y + 20, 26, 34);
                }

                using (Font f = new Font("Segoe UI", 8, FontStyle.Bold))
                using (Brush b = new SolidBrush(Color.FromArgb(70, 70, 70)))
                {
                    g.DrawString("Касса " + model.Cashiers[i].Id.ToString(), f, b, x - 26, y + 82);
                    g.DrawString(model.Cashiers[i].Busy ? "занята" : "свободна", f, b, x - 20, y + 95);
                }

                if (model.Cashiers[i].Busy && model.Cashiers[i].Current != null)
                {
                    DrawCustomer(g, x, y + 56, model.Cashiers[i].Current.Id, Color.FromArgb(244, 162, 97));
                }
            }
        }

        private void DrawQueue(Graphics g, int w, int h)
        {
            if (model == null) return;

            float y = h * 0.66f;
            float right = w * 0.69f;
            float left = w * 0.16f;

            using (Pen p = new Pen(Color.FromArgb(165, 165, 165), 1))
            {
                p.DashStyle = DashStyle.Dot;
                g.DrawLine(p, left, y, right, y);

                int slots = 7;
                float step = (right - left) / (slots - 1);
                for (int i = 0; i < slots; i++)
                    g.DrawEllipse(p, right - i * step - 10, y - 10, 20, 20);
            }

            int waitingCount = model.Waiting.Count;
            int visible = 7;
            int overflow = waitingCount > visible ? waitingCount - visible : 0;
            float step2 = (right - left) / (visible - 1);

            for (int i = 0; i < waitingCount && i < visible; i++)
            {
                var c = model.Waiting[i];
                float targetX = right - i * step2;
                float targetY = y - 20;

                float x = targetX;
                float yy = targetY;

                // короткая анимация входа в очередь
                double delta = model.Time - c.ArrivalTime;
                if (delta < 0.30)
                {
                    float k = (float)(delta / 0.30);
                    if (k < 0f) k = 0f;
                    if (k > 1f) k = 1f;
                    x = Lerp(32f, targetX, k);
                    yy = Lerp(h * 0.56f, targetY, k);
                }

                DrawCustomer(g, x, yy, c.Id, Color.FromArgb(106, 168, 79));
            }

            if (overflow > 0)
            {
                using (Font f = new Font("Segoe UI", 10, FontStyle.Bold))
                using (Brush b = new SolidBrush(Color.FromArgb(70, 70, 70)))
                    g.DrawString("+" + overflow.ToString(), f, b, left - 20, y - 18);
            }

            // рисуем клиентов в обслуживании
            float[] xs = new float[] { w * 0.72f, w * 0.83f, w * 0.94f };
            for (int i = 0; i < model.Cashiers.Count; i++)
            {
                var cashier = model.Cashiers[i];
                if (!cashier.Busy || cashier.Current == null) continue;

                var c = cashier.Current;
                float targetX = xs[i];
                float targetY = h * 0.38f + 56;

                float x = targetX;
                float yy = targetY;

                double sinceStart = model.Time - c.ServiceStart;
                if (sinceStart < 0.25)
                {
                    float k = (float)(sinceStart / 0.25);
                    if (k < 0f) k = 0f;
                    if (k > 1f) k = 1f;

                    float fromX = c.WaitingTime > 0.0001 ? right : 32f;
                    float fromY = c.WaitingTime > 0.0001 ? y - 20f : h * 0.56f;

                    x = Lerp(fromX, targetX, k);
                    yy = Lerp(fromY, targetY, k);
                }

                DrawCustomer(g, x, yy, c.Id, Color.FromArgb(244, 162, 97));
            }
        }

        private void DrawCustomer(Graphics g, float x, float y, int id, Color bodyColor)
        {
            using (Brush skin = new SolidBrush(Color.FromArgb(244, 220, 193)))
                g.FillEllipse(skin, x, y, 20, 20);

            using (Brush body = new SolidBrush(bodyColor))
                g.FillEllipse(body, x - 2, y + 18, 24, 24);

            using (Pen p = new Pen(Color.FromArgb(70, 70, 70), 1.5f))
            {
                g.DrawEllipse(p, x, y, 20, 20);
                g.DrawEllipse(p, x - 2, y + 18, 24, 24);
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