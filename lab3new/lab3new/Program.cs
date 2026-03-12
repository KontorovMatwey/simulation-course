// Program.cs
// Обновлённая версия: треугольники для деревьев, замедленный рост, старение с риском смерти, уменьшенные ползунки UI
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Windows.Forms;

namespace ForestFire_Clean
{
    // Состояния клетки (явные стадии роста)
    public enum CellState { Empty = 0, Seedling = 1, Young = 2, Old = 3, Burning = 4, Water = 5, Rock = 6 }

    // Цвета и визуальные константы
    static class Visual
    {
        public static readonly Color Background = Color.FromArgb(34, 34, 34);
        public static readonly Color Burnt = Color.FromArgb(90, 82, 79);
        public static readonly Color SeedColor = Color.FromArgb(6, 100, 6);
        public static readonly Color YoungColor = Color.FromArgb(20, 160, 20);
        public static readonly Color OldColor = Color.FromArgb(200, 230, 120);
        public static readonly Color Flame1 = Color.FromArgb(255, 210, 183);
        public static readonly Color Flame2 = Color.FromArgb(255, 154, 74);
        public static readonly Color Flame3 = Color.FromArgb(255, 66, 0);
        public static readonly Color Water = Color.FromArgb(31, 79, 138);
        public static readonly Color Rock = Color.SaddleBrown;
    }

    // Одна клетка поля
    public class Cell
    {
        public CellState State = CellState.Empty;
        public int StageAge = 0;        // тики, проведённые в текущей стадии
        public int BurnTicks = 0;      // оставшиеся тики горения
        public double Fertility = 0.0; // 0..1, после пожара повышается
        public int TimeSinceBurn = 9999; // тики после выгорания (визуал/логи)
        public int Elevation = 0;      // рельеф
    }

    // Модель
    public class ForestModel
    {
        public int Width { get; private set; }
        public int Height { get; private set; }
        public Cell[,] Grid;
        public long Ticks = 0;

        // Среда
        public double AmbientHumidity = 0.30; // базовая (ползунок)
        public double Humidity = 0.30;        // текущая
        public double WindAngle = 0.0;        // градусы
        public double WindStrength = 0.15;    // 0..1

        // Дождь — глобальное событие
        public bool RainActive = false;
        public int RainStrength = 0;       // 0..3
        public int RainTicksRemaining = 0;

        // Параметры модели
        public double BaseIgnitionFromNeighbor = 0.28;
        public double BaseGrowthProb = 0.0012; // понижен рост (медленнее)
        public double BurntFertilityBonus = 0.28;
        public int SeedlingTicks = 20;    // было 6 -> 20
        public int YoungTicks = 200;      // было 40 -> 200
        public int MaxBurnTimeBase = 3;
        public int MaxBurnExtra = 5;

        // Старение и шанс смерти у старых
        public int OldMaxAgeTicks = 800;      // условный "максимальный возраст" в тиках
        public double OldDeathChancePerTick = 0.002; // шанс умереть в тик после достижения OldMaxAgeTicks (0.2%)

        // Фоновые поджоги
        private double spontaneousBase = 4e-7; // уменьшено

        // Влажность / влияния
        private const double fireImpact = 0.28;
        private const double rainHumidityPerStrength = 0.02;
        private const double ambientRelaxRate = 0.006;
        private const double evaporation = 0.00012;

        private Random rnd = new Random();

        public ForestModel(int w, int h)
        {
            Width = Math.Max(10, w);
            Height = Math.Max(10, h);
            Grid = new Cell[Height, Width];
            for (int y = 0; y < Height; y++)
                for (int x = 0; x < Width; x++)
                    Grid[y, x] = new Cell();
            InitializeTerrain();
            Humidity = AmbientHumidity;
        }

        private void InitializeTerrain()
        {
            for (int y = 0; y < Height; y++)
                for (int x = 0; x < Width; x++)
                {
                    var c = Grid[y, x];
                    c.Elevation = (int)(5 * Math.Sin(x / 8.0) + 3 * Math.Cos(y / 10.0) + (rnd.NextDouble() * 4 - 2));
                    c.State = CellState.Empty;
                    c.StageAge = 0; c.BurnTicks = 0; c.Fertility = 0; c.TimeSinceBurn = 9999;
                }

            int lakes = Math.Max(1, (Width * Height) / 12000);
            for (int i = 0; i < lakes; i++)
            {
                int cx = rnd.Next(0, Width), cy = rnd.Next(0, Height);
                int radius = rnd.Next(Math.Max(2, Math.Min(Width, Height) / 30), Math.Max(3, Math.Min(Width, Height) / 14));
                for (int y = Math.Max(0, cy - radius); y <= Math.Min(Height - 1, cy + radius); y++)
                    for (int x = Math.Max(0, cx - radius); x <= Math.Min(Width - 1, cx + radius); x++)
                        if ((x - cx) * (x - cx) + (y - cy) * (y - cy) <= radius * radius)
                            Grid[y, x].State = CellState.Water;
            }

            for (int y = 0; y < Height; y++)
                for (int x = 0; x < Width; x++)
                    if (Grid[y, x].State == CellState.Empty && rnd.NextDouble() < 0.006)
                        Grid[y, x].State = CellState.Rock;

            for (int y = 0; y < Height; y++)
                for (int x = 0; x < Width; x++)
                    if (Grid[y, x].State == CellState.Empty && rnd.NextDouble() < 0.47)
                    {
                        double r = rnd.NextDouble();
                        if (r < 0.20) Grid[y, x].State = CellState.Seedling;
                        else if (r < 0.75) Grid[y, x].State = CellState.Young;
                        else Grid[y, x].State = CellState.Old;
                        Grid[y, x].StageAge = rnd.Next(1, 10);
                    }
        }

        public void Reset(int w, int h)
        {
            Width = Math.Max(10, w); Height = Math.Max(10, h);
            Grid = new Cell[Height, Width];
            for (int y = 0; y < Height; y++)
                for (int x = 0; x < Width; x++)
                    Grid[y, x] = new Cell();
            InitializeTerrain();
            Ticks = 0;
            RainActive = false; RainStrength = 0; RainTicksRemaining = 0;
            Humidity = AmbientHumidity;
        }

        public void Randomize()
        {
            for (int y = 0; y < Height; y++)
                for (int x = 0; x < Width; x++)
                {
                    var c = Grid[y, x];
                    if (c.State == CellState.Water || c.State == CellState.Rock) continue;
                    double r = rnd.NextDouble();
                    if (r < 0.45) { c.State = CellState.Seedling; c.StageAge = 0; }
                    else if (r < 0.80) { c.State = CellState.Young; c.StageAge = rnd.Next(0, 10); }
                    else if (r < 0.92) { c.State = CellState.Old; c.StageAge = rnd.Next(0, 20); }
                    else { c.State = CellState.Empty; c.StageAge = 0; }
                    c.BurnTicks = 0; c.Fertility = 0; c.TimeSinceBurn = 9999;
                }
        }

        public void StartRainManual(int strength)
        {
            if (strength < 1) strength = 1;
            if (strength > 3) strength = 3;
            RainActive = true; RainStrength = strength;
            RainTicksRemaining = 50 + 30 * strength;
        }

        private void TryAutoSpawnRain(int rainMaxStrength)
        {
            if (RainActive || rainMaxStrength <= 0) return;
            double spawnProb = 0.0007 + AmbientHumidity * 0.005;
            if (rnd.NextDouble() < spawnProb) StartRainManual(Math.Max(1, rnd.Next(1, rainMaxStrength + 1)));
        }

        private double StageIgnitionFactor(CellState st)
        {
            if (st == CellState.Seedling) return 0.45;
            if (st == CellState.Young) return 1.0;
            if (st == CellState.Old) return 1.6;
            return 1.0;
        }

        public int BurnTime(Cell c)
        {
            double sf = StageIgnitionFactor(c.State);
            int baseT = MaxBurnTimeBase + (int)(MaxBurnExtra * Math.Min(1.0, (sf - 0.45) / 1.6));
            double factor = 1.0 - 0.45 * Humidity;
            if (RainActive) factor *= (1.0 - 0.45 * (RainStrength / 3.0));
            int t = Math.Max(1, (int)(baseT * factor));
            return t;
        }

        private double NeighborTreeFraction(int x, int y)
        {
            int count = 0, total = 0;
            for (int ny = Math.Max(0, y - 1); ny <= Math.Min(Height - 1, y + 1); ny++)
                for (int nx = Math.Max(0, x - 1); nx <= Math.Min(Width - 1, x + 1); nx++)
                {
                    if (nx == x && ny == y) continue;
                    total++;
                    var st = Grid[ny, nx].State;
                    if (st == CellState.Seedling || st == CellState.Young || st == CellState.Old) count++;
                }
            return total == 0 ? 0.0 : (double)count / total;
        }

        public List<Point> Step(int rainMaxStrength)
        {
            Ticks++;
            var changed = new List<Point>();
            int w = Width, h = Height;
            bool[,] willIgnite = new bool[h, w];
            bool[,] willExtinguish = new bool[h, w];

            TryAutoSpawnRain(rainMaxStrength);
            int burningCount = 0;

            for (int y = 0; y < h; y++)
                for (int x = 0; x < w; x++)
                {
                    var c = Grid[y, x];
                    if (c.TimeSinceBurn < 9999) { c.TimeSinceBurn++; c.Fertility = Math.Max(0.0, c.Fertility - 0.0007); }

                    if (c.State == CellState.Burning)
                    {
                        burningCount++;
                        c.BurnTicks--;
                        if (RainActive && RainStrength > 0)
                        {
                            double extinguish = 0.20 + 0.28 * (RainStrength / 3.0);
                            if (rnd.NextDouble() < extinguish) willExtinguish[y, x] = true;
                            else if (c.BurnTicks <= 0) willExtinguish[y, x] = true;
                        }
                        else if (c.BurnTicks <= 0) willExtinguish[y, x] = true;
                    }
                    else if (c.State == CellState.Seedling || c.State == CellState.Young || c.State == CellState.Old)
                    {
                        var burningNeighbors = new List<Point>();
                        for (int ny = Math.Max(0, y - 1); ny <= Math.Min(h - 1, y + 1); ny++)
                            for (int nx = Math.Max(0, x - 1); nx <= Math.Min(w - 1, x + 1); nx++)
                                if (!(nx == x && ny == y) && Grid[ny, nx].State == CellState.Burning)
                                    burningNeighbors.Add(new Point(nx, ny));

                        if (burningNeighbors.Count > 0)
                        {
                            double pNo = 1.0;
                            foreach (var nb in burningNeighbors)
                            {
                                double contrib = BaseIgnitionFromNeighbor;
                                double stageFactor = StageIgnitionFactor(c.State);
                                contrib *= stageFactor;

                                int dx = x - nb.X, dy = y - nb.Y;
                                double angleToCell = (Math.Atan2(dy, dx) * 180.0 / Math.PI + 360) % 360;
                                double windDiff = Math.Abs(((angleToCell - WindAngle + 180 + 360) % 360) - 180);
                                double cosv = Math.Cos(windDiff * Math.PI / 180.0);
                                double windFactor = cosv > 0 ? 1.0 + WindStrength * 2.6 * cosv : Math.Max(0.01, 1.0 + WindStrength * 1.8 * cosv);
                                contrib *= windFactor;

                                int elevDiff = Grid[y, x].Elevation - Grid[nb.Y, nb.X].Elevation;
                                if (elevDiff > 0) contrib *= (1.0 + Math.Min(0.6, elevDiff * 0.08));
                                else contrib *= (1.0 + elevDiff * 0.02);

                                contrib *= Math.Max(0.0, 1.0 - Math.Pow(Humidity, 1.7));
                                if (RainActive && RainStrength > 0) contrib *= Math.Max(0.01, 1.0 - 0.85 * (RainStrength / 3.0));

                                contrib = Math.Max(0.0, Math.Min(1.0, contrib));
                                pNo *= (1.0 - contrib);
                            }
                            double pIgnite = 1.0 - pNo;
                            pIgnite *= (1.0 + Math.Min(2.0, c.Fertility * 12.0));
                            if (rnd.NextDouble() < pIgnite) willIgnite[y, x] = true;
                        }

                        double ageFactor = StageIgnitionFactor(c.State);
                        double pSpont = spontaneousBase * Math.Pow((1.0 - Humidity), 1.8) * (1.0 + 0.8 * (ageFactor - 1.0));
                        if (rnd.NextDouble() < pSpont) willIgnite[y, x] = true;
                    }
                    else if (c.State == CellState.Empty)
                    {
                        double localTreeFrac = NeighborTreeFraction(x, y);
                        double growthProb = BaseGrowthProb * (1.0 + 4.0 * c.Fertility) * (1.0 + localTreeFrac);
                        growthProb *= (0.4 + 0.6 * Math.Min(1.0, Humidity + 0.1));
                        if (Grid[y, x].State == CellState.Water || Grid[y, x].State == CellState.Rock) growthProb = 0.0;
                        if (rnd.NextDouble() < growthProb)
                        {
                            Grid[y, x].State = CellState.Seedling;
                            Grid[y, x].StageAge = 0;
                            changed.Add(new Point(x, y));
                        }
                    }
                }

            for (int y = 0; y < h; y++)
                for (int x = 0; x < w; x++)
                {
                    if (willIgnite[y, x] && (Grid[y, x].State == CellState.Seedling || Grid[y, x].State == CellState.Young || Grid[y, x].State == CellState.Old))
                    {
                        Grid[y, x].State = CellState.Burning;
                        Grid[y, x].BurnTicks = BurnTime(Grid[y, x]);
                        changed.Add(new Point(x, y));
                    }
                    if (willExtinguish[y, x] && Grid[y, x].State == CellState.Burning)
                    {
                        Grid[y, x].State = CellState.Empty;
                        Grid[y, x].StageAge = 0;
                        Grid[y, x].Fertility = Math.Min(1.0, Grid[y, x].Fertility + BurntFertilityBonus);
                        Grid[y, x].TimeSinceBurn = 0;
                        changed.Add(new Point(x, y));
                    }
                }

            for (int y = 0; y < h; y++)
                for (int x = 0; x < w; x++)
                {
                    var c = Grid[y, x];
                    if (c.State == CellState.Seedling)
                    {
                        c.StageAge++;
                        if (c.StageAge >= SeedlingTicks) { c.State = CellState.Young; c.StageAge = 0; changed.Add(new Point(x, y)); }
                    }
                    else if (c.State == CellState.Young)
                    {
                        c.StageAge++;
                        if (c.StageAge >= YoungTicks) { c.State = CellState.Old; c.StageAge = 0; changed.Add(new Point(x, y)); }
                    }
                    else if (c.State == CellState.Old)
                    {
                        c.StageAge++;
                        // риск смерти при слишком старом возрасте
                        if (c.StageAge >= OldMaxAgeTicks && rnd.NextDouble() < OldDeathChancePerTick)
                        {
                            c.State = CellState.Empty;
                            c.StageAge = 0;
                            c.Fertility = Math.Min(1.0, c.Fertility + 0.12);
                            changed.Add(new Point(x, y));
                        }
                    }
                    else if (c.State == CellState.Burning)
                    {
                        c.Fertility = Math.Min(1.0, c.Fertility + 0.006);
                    }

                    if (c.TimeSinceBurn < 9999 && c.State != CellState.Burning) c.Fertility = Math.Max(0.0, c.Fertility - 0.0008);
                }

            if (RainActive && RainStrength > 0)
            {
                Humidity += rainHumidityPerStrength * RainStrength;
                int ashClear = 3 + 3 * RainStrength;
                for (int y = 0; y < h; y++)
                    for (int x = 0; x < w; x++)
                        if (Grid[y, x].TimeSinceBurn < 9999)
                            Grid[y, x].TimeSinceBurn = Math.Min(9999, Grid[y, x].TimeSinceBurn + ashClear);

                RainTicksRemaining--;
                if (RainTicksRemaining <= 0) { RainActive = false; RainStrength = 0; RainTicksRemaining = 0; }
            }

            double towardAmbient = (AmbientHumidity - Humidity) * ambientRelaxRate;
            double burningFraction = (double)burningCount / (w * h);
            double fireEffect = -burningFraction * fireImpact;
            double evap = -evaporation * Humidity;
            Humidity += towardAmbient + fireEffect + evap;
            Humidity = Math.Max(0.0, Math.Min(1.0, Humidity));

            return changed.Distinct().ToList();
        }
    }

    // GUI
    public class MainForm : Form
    {
        private ForestModel model;
        private PictureBox canvas;
        private Bitmap bmp;
        private Timer timer;
        private int cellSize = 12;
        private int gridW = 80, gridH = 50;

        // UI
        private Panel canvasPanel, controlsPanel;
        private Button startBtn, pauseBtn, stepBtn, resetBtn, randBtn, rainBtn;
        private TrackBar speedBar, windBar, windStrBar, humidityBar, zoomBar, rainStrengthBar;
        private Label speedLabel, windAngleLabel, windStrLabel, humLabel, zoomLabel, rainStrengthLabel, rainStatusLabel;
        private ComboBox brushBox;
        private NumericUpDown widthUpDown, heightUpDown;
        private Label infoLabel;

        public MainForm()
        {
            Text = "Лесной пожар — треугольники, старение";
            ClientSize = new Size(1700, 980);
            InitializeModel();
            InitializeUI();
            DrawFull();
        }

        private void InitializeModel()
        {
            model = new ForestModel(gridW, gridH);
            bmp = new Bitmap(Math.Max(1, gridW * cellSize), Math.Max(1, gridH * cellSize));
            timer = new Timer(); timer.Tick += Timer_Tick; timer.Interval = 140;
        }

        private void InitializeUI()
        {
            canvasPanel = new Panel { Location = new Point(8, 8), Size = new Size(1180, 960), AutoScroll = true, BorderStyle = BorderStyle.FixedSingle };
            Controls.Add(canvasPanel);

            canvas = new PictureBox();
            canvas.Size = new Size(gridW * cellSize, gridH * cellSize);
            canvas.Image = bmp;
            canvas.MouseDown += Canvas_MouseDown; canvas.MouseMove += Canvas_MouseMove;
            canvasPanel.Controls.Add(canvas);

            int ctrlX = canvasPanel.Right + 12;
            controlsPanel = new Panel { Location = new Point(ctrlX, 8), Size = new Size(ClientSize.Width - ctrlX - 12, ClientSize.Height - 16), AutoScroll = true };
            Controls.Add(controlsPanel);

            int cx = 8, cy = 8;
            startBtn = new Button { Text = "Запустить", Location = new Point(cx, cy) }; startBtn.Click += (s, e) => { ApplyUiToModel(); timer.Start(); UpdateControlStates(); };
            pauseBtn = new Button { Text = "Пауза", Location = new Point(cx + 120, cy) }; pauseBtn.Click += (s, e) => { timer.Stop(); UpdateControlStates(); };
            stepBtn = new Button { Text = "Шаг", Location = new Point(cx + 240, cy) }; stepBtn.Click += (s, e) => DoStep();
            resetBtn = new Button { Text = "Сброс (новый рельеф)", Location = new Point(cx, cy + 44) }; resetBtn.Click += (s, e) => { model.Reset(gridW, gridH); ApplyUiToModel(); ResizeCanvasBitmap(); DrawFull(); };
            randBtn = new Button { Text = "Randomize", Location = new Point(cx + 220, cy + 44) }; randBtn.Click += (s, e) => { model.Randomize(); DrawFull(); };
            rainBtn = new Button { Text = "Вызвать дождь", Location = new Point(cx, cy + 92) }; rainBtn.Click += (s, e) => { model.StartRainManual(Math.Max(1, rainStrengthBar.Value)); UpdateInfo(); };

            controlsPanel.Controls.Add(startBtn); controlsPanel.Controls.Add(pauseBtn); controlsPanel.Controls.Add(stepBtn);
            controlsPanel.Controls.Add(resetBtn); controlsPanel.Controls.Add(randBtn); controlsPanel.Controls.Add(rainBtn);

            var sizeLbl = new Label { Text = "Размер поля (cells)", Location = new Point(cx, cy + 140), AutoSize = true };
            widthUpDown = new NumericUpDown { Minimum = 10, Maximum = 300, Value = gridW, Location = new Point(cx, cy + 168) };
            heightUpDown = new NumericUpDown { Minimum = 10, Maximum = 200, Value = gridH, Location = new Point(cx + 140, cy + 168) };
            var applySize = new Button { Text = "Применить", Location = new Point(cx + 260, cy + 166) };
            applySize.Click += (s, e) => { if (timer.Enabled) { MessageBox.Show("Поставьте паузу перед изменением размера"); return; } gridW = (int)widthUpDown.Value; gridH = (int)heightUpDown.Value; model.Reset(gridW, gridH); ResizeCanvasBitmap(); DrawFull(); };
            controlsPanel.Controls.Add(sizeLbl); controlsPanel.Controls.Add(widthUpDown); controlsPanel.Controls.Add(heightUpDown); controlsPanel.Controls.Add(applySize);

            int sy = cy + 220;
            var lblSpeed = new Label { Text = "Скорость (мс/тик)", Location = new Point(cx, sy) }; controlsPanel.Controls.Add(lblSpeed);
            speedBar = new TrackBar { Minimum = 10, Maximum = 1000, Value = timer.Interval, Location = new Point(cx, sy + 18), Size = new Size(480, 45) };
            speedBar.ValueChanged += (s, e) => { timer.Interval = speedBar.Value; speedLabel.Text = $"{speedBar.Value} ms"; };
            speedLabel = new Label { Text = $"{speedBar.Value} ms", Location = new Point(cx + 500, sy + 18), AutoSize = true };
            controlsPanel.Controls.Add(speedBar); controlsPanel.Controls.Add(speedLabel);

            var lblWind = new Label { Text = "Направление ветра (°)", Location = new Point(cx, sy + 78) }; controlsPanel.Controls.Add(lblWind);
            windBar = new TrackBar { Minimum = 0, Maximum = 359, Value = (int)model.WindAngle, Location = new Point(cx, sy + 100), Size = new Size(480, 45) };
            windBar.ValueChanged += (s, e) => { model.WindAngle = windBar.Value; windAngleLabel.Text = $"{windBar.Value}°"; };
            windAngleLabel = new Label { Text = $"{windBar.Value}°", Location = new Point(cx + 500, sy + 100), AutoSize = true };
            controlsPanel.Controls.Add(windBar); controlsPanel.Controls.Add(windAngleLabel);

            var lblWindStr = new Label { Text = "Сила ветра (%)", Location = new Point(cx, sy + 148) }; controlsPanel.Controls.Add(lblWindStr);
            windStrBar = new TrackBar { Minimum = 0, Maximum = 100, Value = (int)(model.WindStrength * 100), Location = new Point(cx, sy + 168), Size = new Size(480, 45) };
            windStrBar.ValueChanged += (s, e) => { model.WindStrength = windStrBar.Value / 100.0; windStrLabel.Text = $"{windStrBar.Value}%"; };
            windStrLabel = new Label { Text = $"{windStrBar.Value}%", Location = new Point(cx + 500, sy + 168), AutoSize = true };
            controlsPanel.Controls.Add(windStrBar); controlsPanel.Controls.Add(windStrLabel);

            var lblHum = new Label { Text = "Базовая влажность (%)", Location = new Point(cx, sy + 218) }; controlsPanel.Controls.Add(lblHum);
            humidityBar = new TrackBar { Minimum = 0, Maximum = 100, Value = (int)(model.AmbientHumidity * 100), Location = new Point(cx, sy + 238), Size = new Size(480, 45) };
            humidityBar.ValueChanged += (s, e) => { model.AmbientHumidity = humidityBar.Value / 100.0; humLabel.Text = $"{humidityBar.Value}%"; };
            humLabel = new Label { Text = $"{humidityBar.Value}%", Location = new Point(cx + 500, sy + 238), AutoSize = true };
            controlsPanel.Controls.Add(humidityBar); controlsPanel.Controls.Add(humLabel);

            var lblBrush = new Label { Text = "Кисть (левая кнопка)", Location = new Point(cx, sy + 298) }; controlsPanel.Controls.Add(lblBrush);
            brushBox = new ComboBox { Location = new Point(cx, sy + 318), Width = 240 };
            brushBox.Items.AddRange(new string[] { "seed/plant", "erase", "water", "rock", "ignite" }); brushBox.SelectedIndex = 0; controlsPanel.Controls.Add(brushBox);

            var lblZoom = new Label { Text = "Zoom (px клетка)", Location = new Point(cx, sy + 358) }; controlsPanel.Controls.Add(lblZoom);
            zoomBar = new TrackBar { Minimum = 6, Maximum = 32, Value = cellSize, Location = new Point(cx, sy + 378), Size = new Size(480, 45) };
            zoomBar.ValueChanged += (s, e) => { zoomLabel.Text = $"{zoomBar.Value}px"; if (!timer.Enabled) ApplyZoom(zoomBar.Value); };
            zoomLabel = new Label { Text = $"{zoomBar.Value}px", Location = new Point(cx + 500, sy + 378), AutoSize = true };
            controlsPanel.Controls.Add(zoomBar); controlsPanel.Controls.Add(zoomLabel);

            var lblRain = new Label { Text = "Макс сила дождя (0 = авт. выкл.)", Location = new Point(cx, sy + 428) }; controlsPanel.Controls.Add(lblRain);
            rainStrengthBar = new TrackBar { Minimum = 0, Maximum = 3, Value = 2, TickFrequency = 1, Location = new Point(cx, sy + 448), Size = new Size(480, 45) };
            rainStrengthBar.ValueChanged += (s, e) => { rainStrengthLabel.Text = $"{rainStrengthBar.Value}"; };
            rainStrengthLabel = new Label { Text = $"{rainStrengthBar.Value}", Location = new Point(cx + 500, sy + 448), AutoSize = true };
            controlsPanel.Controls.Add(rainStrengthBar); controlsPanel.Controls.Add(rainStrengthLabel);

            rainStatusLabel = new Label { Text = "Дождь: нет", Location = new Point(cx, sy + 508), AutoSize = true, ForeColor = Color.Gray };
            controlsPanel.Controls.Add(rainStatusLabel);

            infoLabel = new Label { Text = "информация", Location = new Point(cx, sy + 548), Size = new Size(740, 350) };
            controlsPanel.Controls.Add(infoLabel);

            Panel legend = new Panel { Location = new Point(cx, sy + 908), Size = new Size(740, 200), BorderStyle = BorderStyle.FixedSingle };
            controlsPanel.Controls.Add(legend);
            AddLegendItems(legend);

            ResizeCanvasBitmap();
            canvas.Image = bmp;
            UpdateControlStates();
        }

        private void AddLegendItems(Panel legend)
        {
            void add(string name, Action<Graphics, Rectangle> draw, int idx)
            {
                PictureBox pb = new PictureBox { Location = new Point(6, 6 + idx * 30), Size = new Size(28, 28) };
                Bitmap b = new Bitmap(28, 28);
                using (Graphics g = Graphics.FromImage(b)) { g.Clear(Color.Transparent); draw(g, new Rectangle(0, 0, 28, 28)); }
                pb.Image = b; legend.Controls.Add(pb);
                Label lbl = new Label { Location = new Point(42, 8 + idx * 30), Text = name, AutoSize = true };
                legend.Controls.Add(lbl);
            }
            add("Саженец", (g, r) => { using (Brush br = new SolidBrush(Visual.SeedColor)) { Point[] tri = { new Point(r.Left + r.Width / 2, r.Top + 2), new Point(r.Right - 2, r.Bottom - 2), new Point(r.Left + 2, r.Bottom - 2) }; g.FillPolygon(br, tri); } }, 0);
            add("Молодое дерево", (g, r) => { using (Brush br = new SolidBrush(Visual.YoungColor)) { Point[] tri = { new Point(r.Left + r.Width / 2, r.Top), new Point(r.Right, r.Bottom), new Point(r.Left, r.Bottom) }; g.FillPolygon(br, tri); } }, 1);
            add("Старое (сухое)", (g, r) => { using (Brush br = new SolidBrush(Visual.OldColor)) { Point[] tri = { new Point(r.Left + r.Width / 2, r.Top), new Point(r.Right, r.Bottom), new Point(r.Left, r.Bottom) }; g.FillPolygon(br, tri); } }, 2);
            add("Горит", (g, r) => { using (Brush br = new SolidBrush(Visual.Flame3)) g.FillEllipse(br, r); }, 3);
            add("Вода", (g, r) => { using (Brush br = new SolidBrush(Visual.Water)) g.FillRectangle(br, r); }, 4);
            add("Камень", (g, r) => { using (Brush br = new SolidBrush(Visual.Rock)) g.FillRectangle(br, r); }, 5);
            add("Дождь — глобальное событие", (g, r) => { using (Brush br = new SolidBrush(Color.FromArgb(80, 150, 190))) g.FillEllipse(br, r); }, 6);
        }

        private void ResizeCanvasBitmap()
        {
            bmp = new Bitmap(Math.Max(1, gridW * cellSize), Math.Max(1, gridH * cellSize));
            canvas.Image = bmp;
            canvas.Size = new Size(gridW * cellSize, gridH * cellSize);
        }

        private void ApplyZoom(int newSize)
        {
            cellSize = newSize; ResizeCanvasBitmap(); DrawFull();
        }

        private void Canvas_MouseDown(object sender, MouseEventArgs e) { HandlePaint(e); }
        private void Canvas_MouseMove(object sender, MouseEventArgs e) { if (e.Button == MouseButtons.Left) HandlePaint(e); }

        private void HandlePaint(MouseEventArgs e)
        {
            int x = e.X / cellSize, y = e.Y / cellSize;
            if (x < 0 || y < 0 || x >= gridW || y >= gridH) return;
            string brush = brushBox.SelectedItem?.ToString() ?? "seed/plant";
            var c = model.Grid[y, x];
            if (brush == "seed/plant") { if (c.State != CellState.Water && c.State != CellState.Rock) { c.State = CellState.Seedling; c.StageAge = 0; } }
            else if (brush == "erase") { if (c.State != CellState.Water && c.State != CellState.Rock) { c.State = CellState.Empty; c.StageAge = 0; } }
            else if (brush == "water") { c.State = CellState.Water; }
            else if (brush == "rock") { c.State = CellState.Rock; }
            else if (brush == "ignite") { if (c.State == CellState.Seedling || c.State == CellState.Young || c.State == CellState.Old) { c.State = CellState.Burning; c.BurnTicks = model.BurnTime(c); } }
            DrawCell(x, y);
        }

        private void Timer_Tick(object sender, EventArgs e) { DoStep(); }

        private void DoStep()
        {
            ApplyUiToModel();
            var changed = model.Step(rainStrengthBar.Value);
            using (Graphics g = Graphics.FromImage(bmp))
                foreach (var p in changed) DrawCellInternal(g, p.X, p.Y);
            canvas.Invalidate(); UpdateInfo();
        }

        private void DrawFull()
        {
            using (Graphics g = Graphics.FromImage(bmp))
                for (int y = 0; y < gridH; y++) for (int x = 0; x < gridW; x++) DrawCellInternal(g, x, y);
            canvas.Invalidate();
        }

        private void DrawCell(int x, int y) { using (Graphics g = Graphics.FromImage(bmp)) DrawCellInternal(g, x, y); canvas.Invalidate(new Rectangle(x * cellSize, y * cellSize, cellSize, cellSize)); }

        private void DrawCellInternal(Graphics g, int x, int y)
        {
            var c = model.Grid[y, x];
            int px = x * cellSize, py = y * cellSize;
            Rectangle r = new Rectangle(px, py, cellSize, cellSize);
            using (Brush bg = new SolidBrush(Visual.Background)) g.FillRectangle(bg, r);
            int pad = Math.Max(1, cellSize / 8);
            Rectangle inner = new Rectangle(px + pad, py + pad, Math.Max(1, cellSize - 2 * pad), Math.Max(1, cellSize - 2 * pad));

            switch (c.State)
            {
                case CellState.Empty:
                    Color col = (c.TimeSinceBurn < 200) ? Visual.Burnt : Visual.Background;
                    using (Brush b = new SolidBrush(col)) g.FillRectangle(b, inner);
                    break;
                case CellState.Seedling:
                    Point[] triS = { new Point(inner.Left + inner.Width / 2, inner.Top + inner.Height / 4), new Point(inner.Right - 2, inner.Bottom), new Point(inner.Left + 2, inner.Bottom) };
                    using (Brush b = new SolidBrush(Visual.SeedColor)) g.FillPolygon(b, triS);
                    break;
                case CellState.Young:
                    Point[] triY = { new Point(inner.Left + inner.Width / 2, inner.Top), new Point(inner.Right, inner.Bottom), new Point(inner.Left, inner.Bottom) };
                    using (Brush b = new SolidBrush(Visual.YoungColor)) g.FillPolygon(b, triY);
                    using (Pen p = new Pen(Color.SaddleBrown, Math.Max(1, cellSize / 12))) g.DrawLine(p, inner.Left + inner.Width / 2, inner.Bottom, inner.Left + inner.Width / 2, inner.Bottom - inner.Height / 4);
                    break;
                case CellState.Old:
                    Point[] triO = { new Point(inner.Left + inner.Width / 2, inner.Top), new Point(inner.Right, inner.Bottom), new Point(inner.Left, inner.Bottom) };
                    using (Brush b = new SolidBrush(Visual.OldColor)) g.FillPolygon(b, triO);
                    using (Pen p2 = new Pen(Color.SaddleBrown, Math.Max(1, cellSize / 12))) g.DrawLine(p2, inner.Left + inner.Width / 2, inner.Bottom, inner.Left + inner.Width / 2, inner.Bottom - inner.Height / 4);
                    break;
                case CellState.Burning:
                    int phase = (int)((model.Ticks + x * 7 + y * 13) % 6);
                    Color flame = (phase < 2) ? Visual.Flame1 : (phase < 4) ? Visual.Flame2 : Visual.Flame3;
                    using (Brush b = new SolidBrush(flame)) g.FillEllipse(b, inner);
                    break;
                case CellState.Water:
                    using (Brush b = new SolidBrush(Visual.Water)) g.FillRectangle(b, inner);
                    break;
                case CellState.Rock:
                    using (Brush b = new SolidBrush(Visual.Rock)) g.FillRectangle(b, inner);
                    break;
            }
        }

        private void UpdateInfo()
        {
            string rainText = model.RainActive ? $"Да — сила {model.RainStrength}, тиков осталось {model.RainTicksRemaining}" : "Нет";
            rainStatusLabel.Text = $"Дождь: {rainText}"; rainStatusLabel.ForeColor = model.RainActive ? Color.Blue : Color.Gray;
            infoLabel.Text = $"Тики: {model.Ticks}\nВлажность: {model.Humidity:F3} (базовая {model.AmbientHumidity:F2})\nВетер: {model.WindAngle}°  сила {model.WindStrength:F2}\nДождь активен: {model.RainActive}\nРазмер поля: {gridW}×{gridH}, размер клетки: {cellSize}px";
        }

        private void ApplyUiToModel()
        {
            model.AmbientHumidity = humidityBar.Value / 100.0;
            if (model.Ticks == 0) model.Humidity = model.AmbientHumidity;
            model.WindAngle = windBar.Value;
            model.WindStrength = windStrBar.Value / 100.0;
            timer.Interval = speedBar.Value;
        }

        private void UpdateControlStates()
        {
            bool running = timer.Enabled;
            zoomBar.Enabled = !running; widthUpDown.Enabled = !running; heightUpDown.Enabled = !running;
        }

        [STAThread]
        public static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new MainForm());
        }
    }
}
