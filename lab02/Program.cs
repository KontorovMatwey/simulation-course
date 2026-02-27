// Program.cs
using System;
using System.Diagnostics;
using System.Drawing;
using System.Globalization;
using System.Reflection;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Heat1D
{
    static class Program
    {
        [STAThread]
        static void Main()
        {
            ApplicationConfiguration.Initialize();
            Application.Run(new MainForm());
        }
    }

    public class MainForm : Form
    {
        // UI элементы
        PictureBox pb;
        Button btnRun, btnStop, btnRunMatrix;
        ComboBox cbMaterial;
        NumericUpDown nudLengthMm, nudInitT, nudLeftT, nudRightT;
        TextBox tb_h_mm, tb_tau_s;
        Label lblCenter, lblStatus;
        DataGridView dgv;
        CancellationTokenSource cts;

        // безопасный предел по числу отрезков (N) и время работы 
        const int NMAX = 100000; const double t = 2.0;

        public MainForm()
        {
            Text = "Моделирование теплопроводности";
            Width = 1120; Height = 840;
            InitUI();
        }

        // ---------------- UI: создаём элементы и компоновку ----------------
        void InitUI()
        {
            // область визуализации
            pb = new PictureBox { Left = 10, Top = 10, Width = 1080, Height = 300, BorderStyle = BorderStyle.FixedSingle };
            Controls.Add(pb);

            int x = 10, y = 320, labelW = 260, ctrlW = 160, dy = 36;

            // выбор материала
            Controls.Add(new Label { Left = x, Top = y, Width = labelW, Text = "Материал (ρ, c, λ)" });
            cbMaterial = new ComboBox { Left = x + labelW, Top = y, Width = ctrlW, DropDownStyle = ComboBoxStyle.DropDownList };
            cbMaterial.Items.AddRange(new string[] { "Glass (2500,750,1.0)", "Wood (600,1700,0.15)", "Aluminium (2700,900,205)", "Steel (7850,460,54)" });
            cbMaterial.SelectedIndex = 0;
            Controls.Add(cbMaterial);

            // длина по умолчанию 1 мм
            y += dy;
            Controls.Add(new Label { Left = x, Top = y, Width = labelW, Text = "Длина L (мм)" });
            nudLengthMm = new NumericUpDown { Left = x + labelW, Top = y, Width = ctrlW, DecimalPlaces = 3, Minimum = 0.1M, Maximum = 1000M, Value = 5.0M, Increment = 0.1M };
            Controls.Add(nudLengthMm);

            // шаг по пространству (текст)
            y += dy;
            Controls.Add(new Label { Left = x, Top = y, Width = labelW, Text = "Шаг по пространству h (мм)" });
            tb_h_mm = new TextBox { Left = x + labelW, Top = y, Width = ctrlW, Text = "0.1" };
            Controls.Add(tb_h_mm);

            // шаг по времени (текст)
            y += dy;
            Controls.Add(new Label { Left = x, Top = y, Width = labelW, Text = "Шаг по времени τ (с)" });
            tb_tau_s = new TextBox { Left = x + labelW, Top = y, Width = ctrlW, Text = "0.001" };
            Controls.Add(tb_tau_s);

            // температуры
            y += dy;
            Controls.Add(new Label { Left = x, Top = y, Width = labelW, Text = "Начальная температура пластины (°C)" });
            nudInitT = new NumericUpDown { Left = x + labelW, Top = y, Width = ctrlW, DecimalPlaces = 2, Minimum = -1000, Maximum = 1000, Value = 0 };
            Controls.Add(nudInitT);

            y += dy;
            Controls.Add(new Label { Left = x, Top = y, Width = labelW, Text = "Температура слева (°C)" });
            nudLeftT = new NumericUpDown { Left = x + labelW, Top = y, Width = ctrlW, DecimalPlaces = 2, Minimum = -1000, Maximum = 1000, Value = 30 };
            Controls.Add(nudLeftT);

            y += dy;
            Controls.Add(new Label { Left = x, Top = y, Width = labelW, Text = "Температура справа (°C)" });
            nudRightT = new NumericUpDown { Left = x + labelW, Top = y, Width = ctrlW, DecimalPlaces = 2, Minimum = -1000, Maximum = 1000, Value = -10 };
            Controls.Add(nudRightT);

            // всякие кнопки
            btnRun = new Button { Left = 760, Top = y - 4, Width = 120, Height = 36, Text = "Запустить" };
            btnRun.Click += BtnRun_Click;
            Controls.Add(btnRun);

            btnStop = new Button { Left = 890, Top = y - 4, Width = 120, Height = 36, Text = "Остановить", Enabled = false };
            btnStop.Click += BtnStop_Click;
            Controls.Add(btnStop);

            btnRunMatrix = new Button { Left = 760, Top = y + 44, Width = 250, Height = 36, Text = "Запустить матрицу" };
            btnRunMatrix.Click += BtnRunMatrix_Click;
            Controls.Add(btnRunMatrix);

            // метка с результатом
            lblCenter = new Label { Left = 10, Top = y + 44, Width = 720, Height = 28, Text = "Центр: —  | Узлы: —" };
            Controls.Add(lblCenter);

            lblStatus = new Label { Left = 10, Top = y + 80, Width = 1080, Height = 28, Text = "Статус: готов." };
            Controls.Add(lblStatus);

            // таблица результатов
            dgv = new DataGridView { Left = 10, Top = y + 120, Width = 1080, Height = 360, AllowUserToAddRows = false, ReadOnly = true };
            Controls.Add(dgv);

            // для уменьшения мерцания включаем двойную буферизацию у DataGridView (reflection)
            typeof(DataGridView).InvokeMember("DoubleBuffered",
                BindingFlags.SetProperty | BindingFlags.Instance | BindingFlags.NonPublic,
                null, dgv, new object[] { true });
        }

        // ---------------- Вспомогательные методы ----------------

        // парсим ввод, поддерживаем 1e-4 и запятую
        bool TryParseInvariant(string s, out double v)
        {
            s = s.Trim().Replace(',', '.');
            return double.TryParse(s, System.Globalization.NumberStyles.Float, CultureInfo.InvariantCulture, out v);
        }

        // задаём параметры материалов
        void GetMaterial(out double rho, out double c, out double lambda)
        {
            switch (cbMaterial.SelectedIndex)
            {
                case 0: rho = 2500; c = 750; lambda = 1.0; break;     // glass
                case 1: rho = 600; c = 1700; lambda = 0.15; break;    // wood
                case 2: rho = 2700; c = 900; lambda = 205.0; break;   // aluminium
                case 3: rho = 7850; c = 460; lambda = 54.0; break;    // steel
                default: rho = 2500; c = 750; lambda = 1.0; break;
            }
        }

        // ---------------- Обработка кнопки "Запустить" ----------------
        async void BtnRun_Click(object sender, EventArgs e)
        {
            // блокируем кнопки, создаём токен
            btnRun.Enabled = false; btnStop.Enabled = true; btnRunMatrix.Enabled = false;
            cts = new CancellationTokenSource();
            lblStatus.Text = "Статус: выполняется...";

            // читаем параметры, проверяем формат
            double L = (double)nudLengthMm.Value / 1000.0;
            if (!TryParseInvariant(tb_h_mm.Text, out double h_mm)) { MessageBox.Show("Неверный формат h"); ResetButtons(); return; }
            if (!TryParseInvariant(tb_tau_s.Text, out double tau)) { MessageBox.Show("Неверный формат τ"); ResetButtons(); return; }
            double h = h_mm / 1000.0;
            double initT = (double)nudInitT.Value;
            double Tleft = (double)nudLeftT.Value;
            double Tright = (double)nudRightT.Value;
            GetMaterial(out double rho, out double c, out double lambda);

            // вычисляем N; если слишком много узлов - отключаем отрисовку для безопасности
            int N = Math.Max(2, (int)Math.Round(L / h));
            bool draw = true;
            if (N > NMAX)
            {
                draw = false;
                lblStatus.Text = $"N={N} превышает {NMAX} - отрисовка отключена, симуляция выполнится.";
            }

            // запускаем симуляцию в фоне
            var sw = Stopwatch.StartNew();
            var res = await Task.Run(() => RunSim(L, N, rho, c, lambda, initT, Tleft, Tright, tau, t, cts.Token, draw));
            sw.Stop();

            // проверяем, не отменили ли
            if (double.IsNaN(res.centerTemp))
            {
                lblStatus.Text = "Остановлено пользователем.";
            }
            else
            {
                // показываем центр и число узлов
                lblCenter.Text = $"Центр: {res.centerTemp:F4} °C  | Узлы: {N + 1}";
                lblStatus.Text = $"Готово. Время (реальное): {sw.ElapsedMilliseconds} ms";
            }

            ResetButtons();
        }

        void ResetButtons()
        {
            btnRun.Enabled = true; btnStop.Enabled = false; btnRunMatrix.Enabled = true;
        }

        void BtnStop_Click(object sender, EventArgs e)
        {
            // запрос отмены
            btnStop.Enabled = false;
            cts?.Cancel();
            lblStatus.Text = "Запрошена остановка...";
        }

        // ---------------- Запуск матрицы 4x4 ----------------
        async void BtnRunMatrix_Click(object sender, EventArgs e)
        {
            // блокируем кнопки и готовим токен
            btnRun.Enabled = false; btnStop.Enabled = true; btnRunMatrix.Enabled = false;
            cts = new CancellationTokenSource();
            lblStatus.Text = "Статус: выполняется матрица...";

            double L = (double)nudLengthMm.Value / 1000.0;
            double initT = (double)nudInitT.Value;
            double Tleft = (double)nudLeftT.Value;
            double Tright = (double)nudRightT.Value;
            GetMaterial(out double rho, out double c, out double lambda);

            // фиксированные наборы
            double[] hs_mm = new double[] { 0.1, 0.01, 0.001, 0.0001 };
            double[] taus = new double[] { 0.1, 0.01, 0.001, 0.0001 };

            // готовим таблицу
            dgv.SuspendLayout();
            dgv.Columns.Clear(); dgv.Rows.Clear();
            dgv.ColumnCount = hs_mm.Length + 1;
            dgv.Columns[0].Name = "τ \\ h (мм)";
            for (int j = 0; j < hs_mm.Length; j++) dgv.Columns[j + 1].Name = hs_mm[j].ToString("G6", CultureInfo.InvariantCulture);
            for (int i = 0; i < taus.Length; i++) dgv.Rows.Add(taus[i].ToString("G6", CultureInfo.InvariantCulture), "", "", "", "");
            dgv.ResumeLayout();

            // Простейшая последовательная прогонка
            for (int i = 0; i < taus.Length; i++)
            {
                double tau = taus[i];
                for (int j = 0; j < hs_mm.Length; j++)
                {
                    // если пользователь запросил отмену между расчётами (вроде уже невозможно)
                    if (cts.IsCancellationRequested)
                    {
                        lblStatus.Text = "Матрица прервана пользователем.";
                        for (int ii = i; ii < taus.Length; ii++)
                            for (int jj = (ii == i ? j : 0); jj < hs_mm.Length; jj++)
                                dgv.Rows[ii].Cells[jj + 1].Value = "stopped";
                        ResetButtons();
                        return;
                    }

                    double h_mm = hs_mm[j];
                    double h = h_mm / 1000.0;
                    int N = Math.Max(2, (int)Math.Round(L / h));
                    bool draw = !(N > NMAX);

                    // прогон без визуализации - записываем только финальную центральную температуру
                    var res = await Task.Run(() => RunSim(L, N, rho, c, lambda, initT, Tleft, Tright, tau, t, cts.Token, false));

                    if (double.IsNaN(res.centerTemp))
                    {
                        // пользователь отменил в середине расчёта
                        dgv.Rows[i].Cells[j + 1].Value = "stopped";
                        lblStatus.Text = "Матрица прервана пользователем.";
                        ResetButtons();
                        return;
                    }
                    else
                    {
                        // записываем температуру, округлённую до 6 знаков
                        dgv.Rows[i].Cells[j + 1].Value = $"{res.centerTemp:F6} °C";
                        lblStatus.Text = $"Матрица: τ={tau}, h={h_mm} мм - рассчитано";
                    }

                    // обновляем UI, но не делаем никаких внутренних анимаций
                    Application.DoEvents();
                }
            }

            lblStatus.Text = "Матрица завершена.";
            ResetButtons();
        }

        // ---------------- Численный решатель ----------------
        // Возвращаем (центр, поле). Если пользователь отменил - возвращаем center = NaN
        (double centerTemp, double[] field) RunSim(double L, int N, double rho, double c, double lambda, double initT, double Tleft, double Tright, double tau, double Tmodel, CancellationToken token, bool draw)
        {
            int nodes = N + 1;
            double h = L / N;

            // инициализация: вся пластина с initT
            double[] Tn = new double[nodes];
            for (int i = 0; i < nodes; i++) Tn[i] = initT;

            int interior = Math.Max(1, nodes - 2);

            // коэффициенты трёхдиагональной системы (неявная схема)
            double a = -lambda / (h * h);
            double cdiag = -lambda / (h * h);
            double b = rho * c / tau + 2.0 * lambda / (h * h);

            double[] cp = new double[interior], dp = new double[interior], rhs = new double[interior], x = new double[interior];

            int steps = Math.Max(1, (int)Math.Ceiling(Tmodel / tau));
            int updateEvery = draw ? Math.Max(1, steps / 60) : int.MaxValue;
            int center = nodes / 2;

            // палитра заранее
            Color[] palette = BuildPalette(512);

            for (int step = 0; step < steps; step++)
            {
                // если юзер нажал stop, то выйти, вернув NaN
                if (token.IsCancellationRequested) return (double.NaN, Tn);

                // формируем правую часть
                for (int i = 1; i <= nodes - 2; i++) rhs[i - 1] = rho * c / tau * Tn[i];

                // учёт граничных условий
                rhs[0] -= a * Tleft;
                rhs[interior - 1] -= cdiag * Tright;

                // Thomas: прямая прогонка
                double denom = b;
                cp[0] = cdiag / denom; dp[0] = rhs[0] / denom;
                for (int i = 1; i < interior; i++)
                {
                    denom = b - a * cp[i - 1];
                    cp[i] = cdiag / denom;
                    dp[i] = (rhs[i] - a * dp[i - 1]) / denom;
                }
                // обратная подстановка
                x[interior - 1] = dp[interior - 1];
                for (int i = interior - 2; i >= 0; i--) x[i] = dp[i] - cp[i] * x[i + 1];

                // собираем T^{n+1}
                double[] Tnp1 = new double[nodes];
                Tnp1[0] = Tleft; Tnp1[nodes - 1] = Tright;
                for (int i = 1; i <= nodes - 2; i++) Tnp1[i] = x[i - 1];
                Tn = Tnp1;

                // обновляем визуализацию время от времени
                if ((step % updateEvery) == 0 && draw) DrawFast(Tn, palette, h);
            }

            // финальная отрисовка
            if (draw) DrawFast(Tn, palette, h);
            return (Tn[center], Tn);
        }

        // палитра blue -> white -> red
        Color[] BuildPalette(int size)
        {
            Color cold = Color.FromArgb(0, 60, 180);
            Color hot = Color.FromArgb(180, 30, 0);
            Color[] p = new Color[size];
            for (int i = 0; i < size; i++)
            {
                double t = (double)i / (size - 1);
                p[i] = (t < 0.5) ? Lerp(cold, Color.White, t / 0.5) : Lerp(Color.White, hot, (t - 0.5) / 0.5);
            }
            return p;
        }
        Color Lerp(Color a, Color b, double t) => Color.FromArgb((int)Math.Round(a.R * (1 - t) + b.R * t), (int)Math.Round(a.G * (1 - t) + b.G * t), (int)Math.Round(a.B * (1 - t) + b.B * t));

        // отрисовка: полосы по узлам, центр отмечен, подписи целыми числами
        void DrawFast(double[] field, Color[] palette, double h)
        {
            if (pb.InvokeRequired) { pb.Invoke(new Action(() => DrawFast(field, palette, h))); return; }

            int W = pb.Width, H = pb.Height, nodes = field.Length;
            double min = double.MaxValue, max = double.MinValue;
            foreach (var v in field) { if (v < min) min = v; if (v > max) max = v; }
            if (Math.Abs(max - min) < 1e-12) max = min + 1.0;

            Bitmap bmp = new Bitmap(W, H);
            using (Graphics g = Graphics.FromImage(bmp))
            using (SolidBrush br = new SolidBrush(Color.Black))
            using (Pen penGrid = new Pen(Color.FromArgb(100, Color.Black), 1f))
            using (Pen penEdge = new Pen(Color.Black, 1.5f))
            using (Pen penCenter = new Pen(Color.Black, 2f))
            {
                g.Clear(Color.White);

                // заполняем полосами соответствующие узлам
                for (int i = 0; i < nodes; i++)
                {
                    double leftD = i * (W - 1) / (double)(nodes - 1);
                    double rightD = (i + 1) * (W - 1) / (double)(nodes - 1);
                    int x0 = (int)Math.Round(leftD), x1 = (int)Math.Round(rightD);
                    if (i == nodes - 1) x1 = W - 1;
                    int rectW = Math.Max(1, x1 - x0 + 1);

                    double norm = (field[i] - min) / (max - min);
                    int idx = Math.Max(0, Math.Min(palette.Length - 1, (int)Math.Round(norm * (palette.Length - 1))));
                    br.Color = palette[idx];
                    g.FillRectangle(br, x0, 0, rectW, H);
                }

                // сетка сверху (но не рисуем слишком много линий)
                int maxLines = 271;
                int step = nodes > maxLines ? Math.Max(1, nodes / maxLines) : 1;
                for (int i = 0; i < nodes; i += step)
                {
                    int xpos = (int)Math.Round(i * (W - 1) / (double)(nodes - 1));
                    g.DrawLine(penGrid, xpos, 0, xpos, 12);
                }

                // границы
                g.DrawLine(penEdge, 0, 0, 0, H);
                g.DrawLine(penEdge, W - 1, 0, W - 1, H);

                // отмечаем центр толстой линией и подпись
                int centerIdx = nodes / 2;
                int cx = (int)Math.Round(centerIdx * (W - 1) / (double)(nodes - 1));
                g.DrawLine(penCenter, cx, 0, cx, H);
                var centerLabel = "Center";
                SizeF csz = g.MeasureString(centerLabel, SystemFonts.DefaultFont);
                g.FillRectangle(Brushes.White, cx - csz.Width / 2 - 2, H / 2 - csz.Height / 2 - 2, csz.Width + 4, csz.Height + 4);
                g.DrawString(centerLabel, SystemFonts.DefaultFont, Brushes.Black, cx - csz.Width / 2, H / 2 - csz.Height / 2);

                // подписи: целые градусы для визуальной читаемости
                string leftS = $"{(int)Math.Round(field[0])}°C";
                string rightS = $"{(int)Math.Round(field[nodes - 1])}°C";
                string minS = $"min {(int)Math.Round(min)}°C";
                string maxS = $"max {(int)Math.Round(max)}°C";

                g.DrawString(leftS, SystemFonts.DefaultFont, Brushes.Black, 6, 6);
                SizeF rs = g.MeasureString(rightS, SystemFonts.DefaultFont);
                g.DrawString(rightS, SystemFonts.DefaultFont, Brushes.Black, W - rs.Width - 6, 6);

                g.DrawString(minS, SystemFonts.DefaultFont, Brushes.Black, 6, H - 18);
                SizeF mxs = g.MeasureString(maxS, SystemFonts.DefaultFont);
                g.DrawString(maxS, SystemFonts.DefaultFont, Brushes.Black, W - mxs.Width - 6, H - 18);
            }

            var old = pb.Image; pb.Image = bmp; old?.Dispose();
        }
    }
}