using System;
using System.Collections.Generic;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Windows.Forms;

namespace MM11SimpleLab
{
    public class CustomerEvent
    {
        public int Id { get; set; }
        public double ArrivalTime { get; set; }
        public bool Accepted { get; set; }
        public double ServiceStart { get; set; }
        public double ServiceEnd { get; set; }
        public double RejectExitTime { get; set; }
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
        private NumericUpDown nudMu;
        private NumericUpDown nudT;
        private Button btnRun;
        private Button btnReplay;
        private Label lblProb;
        private Label lblLive;
        private SceneControl scene;
        private Timer timer;

        private List<CustomerEvent> events = new List<CustomerEvent>();
        private double horizon;
        private double currentTime;
        private bool animating;
        private double busyTime;
        private int tickCounter;

        public MainForm()
        {
            Text = "M/M/1";
            Width = 1320;
            Height = 780;
            StartPosition = FormStartPosition.CenterScreen;
            Font = new Font("Segoe UI", 10);

            BuildUi();

            timer = new Timer();
            timer.Interval = 50;
            timer.Tick += Timer_Tick;

            RunSimulation();
        }

        private void BuildUi()
        {
            TableLayoutPanel root = new TableLayoutPanel();
            root.Dock = DockStyle.Fill;
            root.ColumnCount = 2;
            root.RowCount = 1;
            root.Padding = new Padding(10);
            root.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 340f));
            root.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100f));
            Controls.Add(root);

            BuildLeftPanel(root);
            BuildRightPanel(root);
        }

        private void BuildLeftPanel(TableLayoutPanel root)
        {
            TableLayoutPanel left = new TableLayoutPanel();
            left.Dock = DockStyle.Fill;
            left.ColumnCount = 1;
            left.RowCount = 5;
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 190f));
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 90f));
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 250f));
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 90f));
            left.RowStyles.Add(new RowStyle(SizeType.Absolute, 40f));
            root.Controls.Add(left, 0, 0);

            GroupBox gbParams = new GroupBox();
            gbParams.Text = "Параметры";
            gbParams.Dock = DockStyle.Fill;
            left.Controls.Add(gbParams, 0, 0);

            TableLayoutPanel p = new TableLayoutPanel();
            p.Dock = DockStyle.Fill;
            p.ColumnCount = 2;
            p.RowCount = 3;
            p.Padding = new Padding(10);
            p.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 60f));
            p.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 40f));
            gbParams.Controls.Add(p);

            nudLambda = MakeNud(0.05m, 100m, 0.70m, 0.05m);
            nudMu = MakeNud(0.05m, 100m, 0.80m, 0.05m);
            nudT = MakeNud(1m, 1000m, 20m, 1m);

            AddRow(p, 0, "Интенсивность λ:", nudLambda);
            AddRow(p, 1, "Интенсивность обслуживания:", nudMu);
            AddRow(p, 2, "Время моделирования T:", nudT);

            GroupBox gbProb = new GroupBox();
            gbProb.Text = "Эмпирическая вероятность";
            gbProb.Dock = DockStyle.Fill;
            left.Controls.Add(gbProb, 0, 1);

            lblProb = new Label();
            lblProb.Dock = DockStyle.Fill;
            lblProb.Padding = new Padding(10);
            lblProb.TextAlign = ContentAlignment.TopLeft;
            gbProb.Controls.Add(lblProb);

            GroupBox gbLive = new GroupBox();
            gbLive.Text = "Текущая модель";
            gbLive.Dock = DockStyle.Fill;
            left.Controls.Add(gbLive, 0, 2);

            lblLive = new Label();
            lblLive.Dock = DockStyle.Fill;
            lblLive.Padding = new Padding(10);
            lblLive.TextAlign = ContentAlignment.TopLeft;
            gbLive.Controls.Add(lblLive);

            btnRun = new Button();
            btnRun.Text = "Новый расчёт";
            btnRun.Dock = DockStyle.Fill;
            btnRun.Height = 28;
            btnRun.Click += delegate { RunSimulation(); };
            left.Controls.Add(btnRun, 0, 3);

            btnReplay = new Button();
            btnReplay.Text = "Анимация заново";
            btnReplay.Dock = DockStyle.Fill;
            btnReplay.Height = 28;
            btnReplay.Click += delegate { StartAnimation(); };
            left.Controls.Add(btnReplay, 0, 4);
        }

        private void BuildRightPanel(TableLayoutPanel root)
        {
            GroupBox gb = new GroupBox();
            gb.Dock = DockStyle.Fill;
            root.Controls.Add(gb, 1, 0);

            scene = new SceneControl();
            scene.Dock = DockStyle.Fill;
            gb.Controls.Add(scene);
        }

        private NumericUpDown MakeNud(decimal min, decimal max, decimal value, decimal inc)
        {
            NumericUpDown nud = new NumericUpDown();
            nud.Minimum = min;
            nud.Maximum = max;
            nud.Value = value;
            nud.DecimalPlaces = 2;
            nud.Increment = inc;
            nud.Dock = DockStyle.Fill;
            return nud;
        }

        private void AddRow(TableLayoutPanel panel, int row, string text, Control control)
        {
            panel.RowStyles.Add(new RowStyle(SizeType.Absolute, 34f));

            Label lbl = new Label();
            lbl.Text = text;
            lbl.Dock = DockStyle.Fill;
            lbl.TextAlign = ContentAlignment.MiddleLeft;

            panel.Controls.Add(lbl, 0, row);
            panel.Controls.Add(control, 1, row);
        }

        private void RunSimulation()
        {
            double lambda = (double)nudLambda.Value;
            double mu = (double)nudMu.Value;
            horizon = (double)nudT.Value;

            Random rnd = new Random(); //42
            events = Simulate(lambda, mu, horizon, rnd, out busyTime);

            currentTime = 0.0;
            tickCounter = 0;
            animating = true;
            timer.Start();

            UpdateProbabilityLabels(lambda, mu, currentTime);
            UpdateLiveLabels();
            scene.SetData(events, currentTime, horizon);
        }

        private List<CustomerEvent> Simulate(double lambda, double mu, double T, Random rnd, out double busy)
        {
            List<CustomerEvent> list = new List<CustomerEvent>();
            busy = 0.0;

            double time = 0.0;
            double serverFreeAt = 0.0;
            int id = 1;

            while (true)
            {
                time += Exp(lambda, rnd);
                if (time > T) break;

                bool accepted = time >= serverFreeAt;
                if (accepted)
                {
                    double service = Exp(mu, rnd);
                    double end = time + service;

                    list.Add(new CustomerEvent
                    {
                        Id = id++,
                        ArrivalTime = time,
                        Accepted = true,
                        ServiceStart = time,
                        ServiceEnd = end,
                        RejectExitTime = 0.0
                    });

                    busy += Math.Max(0.0, Math.Min(end, T) - time);
                    serverFreeAt = end;
                }
                else
                {
                    list.Add(new CustomerEvent
                    {
                        Id = id++,
                        ArrivalTime = time,
                        Accepted = false,
                        ServiceStart = 0.0,
                        ServiceEnd = 0.0,
                        RejectExitTime = time + 0.7
                    });
                }
            }

            return list;
        }

        private double Exp(double rate, Random rnd)
        {
            double u = 1.0 - rnd.NextDouble();
            return -Math.Log(u) / rate;
        }

        private double GetBusyTimeUpTo(double t)
        {
            if (events == null || events.Count == 0 || t <= 0.0)
                return 0.0;

            double sum = 0.0;

            for (int i = 0; i < events.Count; i++)
            {
                if (!events[i].Accepted)
                    continue;

                double a = events[i].ServiceStart;
                double b = events[i].ServiceEnd;

                if (b <= 0.0 || a >= t)
                    continue;

                double left = a;
                double right = Math.Min(b, t);
                if (right > left)
                    sum += right - left;
            }

            return sum;
        }

        private void UpdateProbabilityLabels(double lambda, double mu, double currentTime)
        {
            double pBusyEmp = 0.0;
            double pFreeEmp = 1.0;

            if (currentTime > 0.0)
            {
                double busyNow = GetBusyTimeUpTo(currentTime);
                pBusyEmp = busyNow / currentTime;
                pFreeEmp = 1.0 - pBusyEmp;
            }

            lblProb.Text =
                "P(занят) = " + pBusyEmp.ToString("F4") + "\r\n" +
                "P(свободен) = " + pFreeEmp.ToString("F4") + "\r\n\r\n";
        }

        private void UpdateLiveLabels()
        {
            int arrived = GetArrivedCount(currentTime);
            int accepted = GetAcceptedCount(currentTime);
            int rejected = GetRejectedCount(currentTime);
            int inService = GetInServiceCount(currentTime);
            int inSystem = GetVisibleCount(currentTime);

            double avgWait = GetAverageWaitingOfAccepted(currentTime);

            lblLive.Text =
                "t = " + currentTime.ToString("F2") + "\r\n" +
                "Всего заявок = " + events.Count.ToString() + "\r\n" +
                "Пришло = " + arrived.ToString() + "\r\n" +
                "Принято = " + accepted.ToString() + "\r\n" +
                "Отказано = " + rejected.ToString() + "\r\n" +
                "В системе = " + inSystem.ToString() + "\r\n";
        }

        private int GetArrivedCount(double t)
        {
            int c = 0;
            for (int i = 0; i < events.Count; i++)
                if (events[i].ArrivalTime <= t) c++;
            return c;
        }

        private int GetAcceptedCount(double t)
        {
            int c = 0;
            for (int i = 0; i < events.Count; i++)
                if (events[i].Accepted && events[i].ArrivalTime <= t) c++;
            return c;
        }

        private int GetRejectedCount(double t)
        {
            int c = 0;
            for (int i = 0; i < events.Count; i++)
                if (!events[i].Accepted && events[i].ArrivalTime <= t) c++;
            return c;
        }

        private int GetInServiceCount(double t)
        {
            int c = 0;
            for (int i = 0; i < events.Count; i++)
                if (events[i].Accepted && events[i].ServiceStart <= t && events[i].ServiceEnd > t) c++;
            return c;
        }

        private int GetVisibleCount(double t)
        {
            int c = 0;
            for (int i = 0; i < events.Count; i++)
            {
                if (events[i].Accepted)
                {
                    if (events[i].ArrivalTime <= t && events[i].ServiceEnd > t)
                        c++;
                }
                else
                {
                    if (events[i].ArrivalTime <= t && events[i].RejectExitTime > t)
                        c++;
                }
            }
            return c;
        }

        private double GetAverageWaitingOfAccepted(double t)
        {
            double sum = 0.0;
            int c = 0;

            for (int i = 0; i < events.Count; i++)
            {
                if (events[i].Accepted && events[i].ServiceEnd <= t)
                {
                    sum += events[i].ServiceStart - events[i].ArrivalTime;
                    c++;
                }
            }

            return c == 0 ? 0.0 : sum / c;
        }

        private void StartAnimation()
        {
            currentTime = 0.0;
            tickCounter = 0;
            animating = true;
            timer.Start();
            scene.SetData(events, currentTime, horizon);
            UpdateProbabilityLabels((double)nudLambda.Value, (double)nudMu.Value, currentTime);
            UpdateLiveLabels();
        }

        private void Timer_Tick(object sender, EventArgs e)
        {
            if (!animating)
                return;

            currentTime += horizon / 500.0;
            if (currentTime >= horizon)
            {
                currentTime = horizon;
                animating = false;
                timer.Stop();
            }

            tickCounter++;

            scene.SetData(events, currentTime, horizon);
            UpdateLiveLabels();

            UpdateProbabilityLabels((double)nudLambda.Value, (double)nudMu.Value, currentTime);
        }
    }

    public class SceneControl : Control
    {
        private List<CustomerEvent> events = new List<CustomerEvent>();
        private double currentTime;
        private double horizon;

        public SceneControl()
        {
            DoubleBuffered = true;
            BackColor = Color.White;
            SetStyle(ControlStyles.AllPaintingInWmPaint | ControlStyles.UserPaint | ControlStyles.OptimizedDoubleBuffer, true);
        }

        public void SetData(List<CustomerEvent> list, double t, double T)
        {
            events = list ?? new List<CustomerEvent>();
            currentTime = t;
            horizon = T;
            Invalidate();
        }

        protected override void OnPaint(PaintEventArgs e)
        {
            base.OnPaint(e);

            Graphics g = e.Graphics;
            g.SmoothingMode = SmoothingMode.AntiAlias;

            int w = Width;
            int h = Height;

            g.Clear(Color.FromArgb(250, 248, 242));

            using (SolidBrush floor = new SolidBrush(Color.FromArgb(238, 230, 220)))
                g.FillRectangle(floor, 0, (int)(h * 0.67), w, (int)(h * 0.33));

            DrawTitle(g, w);
            DrawDoor(g, h);
            DrawCounter(g, w, h);
            DrawOperator(g, w, h);
            DrawCustomers(g, w, h);

            using (Font f = new Font("Segoe UI", 10, FontStyle.Bold))
            using (Brush b = new SolidBrush(Color.FromArgb(70, 70, 70)))
                g.DrawString("Время: " + currentTime.ToString("F2") + " / " + horizon.ToString("F2"), f, b, 14, 12);
        }

        private void DrawTitle(Graphics g, int w)
        {
            using (Font f = new Font("Segoe UI", 14, FontStyle.Bold))
            using (Brush b = new SolidBrush(Color.FromArgb(95, 60, 31)))
                g.DrawString("M/M/1", f, b, w / 2 - 45, 10);
        }

        private void DrawDoor(Graphics g, int h)
        {
            RectangleF door = new RectangleF(24, h * 0.42f, 60, 120);
            using (SolidBrush br = new SolidBrush(Color.FromArgb(215, 200, 185)))
                g.FillRectangle(br, door);
            using (Pen p = new Pen(Color.FromArgb(120, 100, 80), 2))
                g.DrawRectangle(p, door.X, door.Y, door.Width, door.Height);
        }

        private void DrawCounter(Graphics g, int w, int h)
        {
            float x = w * 0.68f;
            float y = h * 0.40f;
            float cw = w * 0.26f;
            float ch = h * 0.18f;

            RectangleF desk = new RectangleF(x, y, cw, ch);
            using (SolidBrush br = new SolidBrush(Color.FromArgb(175, 120, 70)))
                g.FillRectangle(br, desk);
            using (Pen p = new Pen(Color.FromArgb(120, 76, 38), 2))
                g.DrawRectangle(p, desk.X, desk.Y, desk.Width, desk.Height);
        }

        private void DrawOperator(Graphics g, int w, int h)
        {
            float x = w * 0.80f;
            float y = h * 0.29f;

            using (Brush skin = new SolidBrush(Color.FromArgb(244, 207, 182)))
                g.FillEllipse(skin, x, y, 36, 36);
            using (Brush hair = new SolidBrush(Color.FromArgb(85, 55, 25)))
                g.FillEllipse(hair, x - 1, y - 2, 38, 18);
            using (Brush body = new SolidBrush(Color.FromArgb(79, 134, 171)))
                g.FillRectangle(body, x - 4, y + 32, 44, 58);
            using (Pen p = new Pen(Color.FromArgb(60, 60, 60), 2))
            {
                g.DrawEllipse(p, x, y, 36, 36);
                g.DrawRectangle(p, x - 4, y + 32, 44, 58);
            }

            using (Font f = new Font("Segoe UI", 8, FontStyle.Bold))
            using (Brush b = new SolidBrush(Color.FromArgb(70, 70, 70)))
                g.DrawString("оператор", f, b, x - 7, y + 92);
        }

        private void DrawCustomers(Graphics g, int w, int h)
        {
            if (events == null || events.Count == 0) return;

            float entranceX = 32f;
            float entranceY = h * 0.56f;

            float counterX = w * 0.73f;
            float counterY = h * 0.50f;

            float exitX = w - 38f;
            float exitY = h * 0.56f;

            for (int i = 0; i < events.Count; i++)
            {
                CustomerEvent ev = events[i];

                if (ev.Accepted)
                {
                    if (currentTime < ev.ArrivalTime || currentTime >= ev.ServiceEnd)
                        continue;

                    float x = counterX;
                    float y = counterY;

                    if (currentTime < ev.ServiceStart + 0.35)
                    {
                        float k = (float)((currentTime - ev.ServiceStart) / 0.35);
                        if (k < 0f) k = 0f;
                        if (k > 1f) k = 1f;

                        x = Lerp(entranceX, counterX, k);
                        y = Lerp(entranceY, counterY, k);
                    }

                    DrawPerson(g, x, y, ev.Id, Color.FromArgb(244, 162, 97));
                }
                else
                {
                    if (currentTime < ev.ArrivalTime || currentTime >= ev.RejectExitTime)
                        continue;

                    float rejectMidX = w * 0.55f;
                    float rejectMidY = h * 0.45f;

                    float x;
                    float y;

                    if (currentTime < ev.ArrivalTime + 0.25)
                    {
                        float k = (float)((currentTime - ev.ArrivalTime) / 0.25);
                        if (k < 0f) k = 0f;
                        if (k > 1f) k = 1f;

                        x = Lerp(entranceX, rejectMidX, k);
                        y = Lerp(entranceY, rejectMidY, k);
                    }
                    else
                    {
                        double part = ev.RejectExitTime - (ev.ArrivalTime + 0.25);
                        if (part <= 0.0) part = 0.0001;

                        float k = (float)((currentTime - (ev.ArrivalTime + 0.25)) / part);
                        if (k < 0f) k = 0f;
                        if (k > 1f) k = 1f;

                        x = Lerp(rejectMidX, exitX, k);
                        y = Lerp(rejectMidY, exitY, k);
                    }

                    DrawPerson(g, x, y, ev.Id, Color.FromArgb(153, 102, 204));
                }
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