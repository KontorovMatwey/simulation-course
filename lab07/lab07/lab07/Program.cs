using System;
using System.Collections.Generic;
using System.Drawing;
using System.Globalization;
using System.IO;
using System.Text;
using System.Windows.Forms;

namespace WeatherMarkovModelWinForms
{
    internal static class Program
    {
        [STAThread]
        private static void Main()
        {
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new MainForm());
        }
    }

    public sealed class MainForm : Form
    {
        private const int StateCount = 5;

        private readonly string[] _stateNames =
        {
            "", "Ясно", "Перем. облачность", "Облачно", "Пасмурно", "Дождь"
        };

        private readonly string[] _stateShortNames =
        {
            "", "Ясн", "П.о.", "Обл", "Пасм", "Дожд"
        };

        private readonly double[,] _defaultQ =
        {
            { 0,    0.30, 0.18, 0.08, 0.04 },
            { 0.22, 0,    0.24, 0.10, 0.05 },
            { 0.10, 0.18, 0,    0.20, 0.10 },
            { 0.06, 0.12, 0.22, 0,    0.18 },
            { 0.03, 0.08, 0.17, 0.30, 0    }
        };

        private readonly double[,] _q = new double[StateCount, StateCount];
        private readonly double[] _empirical = new double[StateCount + 1];
        private double[] _stationary = new double[StateCount];

        private readonly int[] _counts = new int[StateCount + 1];
        private readonly List<int> _dailyStates = new List<int>();
        private readonly List<DateTime> _dailyDates = new List<DateTime>();

        private readonly Random _rng = new Random();
        private readonly Timer _timer = new Timer();

        private readonly DataGridView _qGrid = new DataGridView();
        private readonly DataGridView _statsGrid = new DataGridView();
        private readonly TimelinePanel _timelinePanel = new TimelinePanel();
        private readonly DistributionPanel _distPanel = new DistributionPanel();
        private readonly TextBox _reportBox = new TextBox();
        private readonly ComboBox _monthCombo = new ComboBox();
        private readonly ComboBox _durationCombo = new ComboBox();
        private readonly NumericUpDown _delayUpDown = new NumericUpDown();
        private readonly Button _startButton = new Button();
        private readonly Button _resetButton = new Button();
        private readonly Button _defaultsButton = new Button();
        private readonly Label _currentLabel = new Label();
        private readonly Label _csvLabel = new Label();
        private readonly Label _qHintLabel = new Label();

        private bool _updatingGrid = false;
        private bool _running = false;
        private int _currentState = 1;
        private int _targetHours = 24;
        private int _currentDayIndex = 0;
        private DateTime _startDate = new DateTime(DateTime.Now.Year, DateTime.Now.Month, 1);

        private string _dailyCsvPath = string.Empty;
        private string _summaryCsvPath = string.Empty;
        private string _reportTxtPath = string.Empty;

        public MainForm()
        {
            Text = "Марковская модель погоды";
            StartPosition = FormStartPosition.CenterScreen;
            Width = 1700;
            Height = 1020;
            MinimumSize = new Size(1500, 900);
            Font = new Font("Segoe UI", 9F, FontStyle.Regular, GraphicsUnit.Point);

            BuildUi();
            LoadDefaults();
            ConfigureTimer();
            UpdateCurrentLabel("Готово к запуску");
            UpdateCsvLabel();
        }

        private void BuildUi()
        {
            Panel topPanel = new Panel();
            topPanel.Dock = DockStyle.Top;
            topPanel.Height = 110;
            topPanel.Padding = new Padding(10);
            Controls.Add(topPanel);

            Panel mainPanel = new Panel();
            mainPanel.Dock = DockStyle.Fill;
            mainPanel.Padding = new Padding(10);
            Controls.Add(mainPanel);

            Panel leftPanel = new Panel();
            leftPanel.Dock = DockStyle.Fill;
            leftPanel.Padding = new Padding(0, 0, 8, 0);
            mainPanel.Controls.Add(leftPanel);

            Panel rightPanel = new Panel();
            rightPanel.Dock = DockStyle.Right;
            rightPanel.Width = 520;
            rightPanel.Padding = new Padding(8, 0, 0, 0);
            mainPanel.Controls.Add(rightPanel);

            BuildTopPanel(topPanel);
            BuildLeftPanel(leftPanel);
            BuildRightPanel(rightPanel);

            Shown += MainForm_Shown;
        }

        private void BuildTopPanel(Panel topPanel)
        {
            TableLayoutPanel layout = new TableLayoutPanel();
            layout.Dock = DockStyle.Fill;
            layout.ColumnCount = 8;
            layout.RowCount = 2;
            layout.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 130));
            layout.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 170));
            layout.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 110));
            layout.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 110));
            layout.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 90));
            layout.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 90));
            layout.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 130));
            layout.ColumnStyles.Add(new ColumnStyle(SizeType.Absolute, 130));
            layout.RowStyles.Add(new RowStyle(SizeType.Absolute, 34));
            layout.RowStyles.Add(new RowStyle(SizeType.Absolute, 34));
            topPanel.Controls.Add(layout);

            Label monthLabel = new Label();
            monthLabel.Text = "Стартовый месяц:";
            monthLabel.Dock = DockStyle.Fill;
            monthLabel.TextAlign = ContentAlignment.MiddleLeft;
            layout.Controls.Add(monthLabel, 0, 0);

            _monthCombo.Dock = DockStyle.Fill;
            _monthCombo.DropDownStyle = ComboBoxStyle.DropDownList;
            _monthCombo.Items.AddRange(new object[]
            {
                "Январь", "Февраль", "Март", "Апрель", "Май", "Июнь",
                "Июль", "Август", "Сентябрь", "Октябрь", "Ноябрь", "Декабрь"
            });
            _monthCombo.SelectedIndex = DateTime.Now.Month - 1;
            layout.Controls.Add(_monthCombo, 1, 0);

            Label durationLabel = new Label();
            durationLabel.Text = "Длительность (часы):";
            durationLabel.Dock = DockStyle.Fill;
            durationLabel.TextAlign = ContentAlignment.MiddleLeft;
            layout.Controls.Add(durationLabel, 2, 0);

            _durationCombo.Dock = DockStyle.Fill;
            _durationCombo.DropDownStyle = ComboBoxStyle.DropDownList;
            _durationCombo.Items.AddRange(new object[] { 24, 168, 720, 8760 });
            _durationCombo.SelectedIndex = 0;
            layout.Controls.Add(_durationCombo, 3, 0);

            Label delayLabel = new Label();
            delayLabel.Text = "Пауза, мс:";
            delayLabel.Dock = DockStyle.Fill;
            delayLabel.TextAlign = ContentAlignment.MiddleLeft;
            layout.Controls.Add(delayLabel, 4, 0);

            _delayUpDown.Dock = DockStyle.Fill;
            _delayUpDown.Minimum = 10;
            _delayUpDown.Maximum = 2000;
            _delayUpDown.Value = 50;
            _delayUpDown.Increment = 10;
            layout.Controls.Add(_delayUpDown, 5, 0);

            _startButton.Text = "Старт";
            _startButton.Dock = DockStyle.Fill;
            _startButton.Click += StartButton_Click;
            layout.Controls.Add(_startButton, 6, 0);

            _resetButton.Text = "Сброс";
            _resetButton.Dock = DockStyle.None;
            _resetButton.AutoSize = true;
            _resetButton.Anchor = AnchorStyles.Left;
            _resetButton.Click += ResetButton_Click;
            layout.Controls.Add(_resetButton, 7, 0);

            _defaultsButton.Text = "Q по умолчанию";
            _defaultsButton.Dock = DockStyle.Fill;
            _defaultsButton.Click += DefaultsButton_Click;
            layout.Controls.Add(_defaultsButton, 6, 1);

            Label infoLabel = new Label();
            infoLabel.Text = "1 шаг = 1 час. Вне диагонали вводятся интенсивности переходов qij >= 0. Диагональ считается автоматически.";
            infoLabel.Dock = DockStyle.Fill;
            infoLabel.TextAlign = ContentAlignment.MiddleLeft;
            infoLabel.Margin = new Padding(0, 6, 0, 0);
            layout.Controls.Add(infoLabel, 0, 1);
            layout.SetColumnSpan(infoLabel, 6);
        }

        private void BuildLeftPanel(Panel leftPanel)
        {
            TableLayoutPanel leftLayout = new TableLayoutPanel();
            leftLayout.Dock = DockStyle.Fill;
            leftLayout.ColumnCount = 1;
            leftLayout.RowCount = 2;
            leftLayout.RowStyles.Add(new RowStyle(SizeType.Absolute, 600));
            leftLayout.RowStyles.Add(new RowStyle(SizeType.Percent, 100));
            leftPanel.Controls.Add(leftLayout);

            GroupBox timelineGroup = new GroupBox();
            timelineGroup.Text = "Визуализация по часам";
            timelineGroup.Dock = DockStyle.Fill;
            timelineGroup.Padding = new Padding(8, 28, 8, 8);
            leftLayout.Controls.Add(timelineGroup, 0, 0);

            _timelinePanel.Dock = DockStyle.Fill;
            _timelinePanel.BackColor = Color.White;
            _timelinePanel.AutoScroll = true;
            timelineGroup.Controls.Add(_timelinePanel);

            GroupBox reportGroup = new GroupBox();
            reportGroup.Text = "Итоговый отчёт";
            reportGroup.Dock = DockStyle.Fill;
            leftLayout.Controls.Add(reportGroup, 0, 1);

            _reportBox.Dock = DockStyle.Fill;
            _reportBox.Multiline = true;
            _reportBox.ReadOnly = true;
            _reportBox.ScrollBars = ScrollBars.Vertical;
            _reportBox.Font = new Font("Consolas", 9F, FontStyle.Regular, GraphicsUnit.Point);
            reportGroup.Controls.Add(_reportBox);
        }

        private void BuildRightPanel(Panel rightPanel)
        {
            TableLayoutPanel rightLayout = new TableLayoutPanel();
            rightLayout.Dock = DockStyle.Fill;
            rightLayout.ColumnCount = 1;
            rightLayout.RowCount = 4;
            rightLayout.RowStyles.Add(new RowStyle(SizeType.Absolute, 52));
            rightLayout.RowStyles.Add(new RowStyle(SizeType.Absolute, 300));
            rightLayout.RowStyles.Add(new RowStyle(SizeType.Absolute, 210));
            rightLayout.RowStyles.Add(new RowStyle(SizeType.Percent, 100));
            rightPanel.Controls.Add(rightLayout);

            _currentLabel.Dock = DockStyle.Fill;
            _currentLabel.TextAlign = ContentAlignment.MiddleCenter;
            _currentLabel.Font = new Font("Segoe UI", 17F, FontStyle.Bold, GraphicsUnit.Point);
            _currentLabel.BorderStyle = BorderStyle.FixedSingle;
            rightLayout.Controls.Add(_currentLabel, 0, 0);

            GroupBox qGroup = new GroupBox();
            qGroup.Text = "Матрица интенсивностей Q";
            qGroup.Dock = DockStyle.Fill;
            qGroup.Padding = new Padding(8, 20, 8, 8);
            rightLayout.Controls.Add(qGroup, 0, 1);

            Panel qInner = new Panel();
            qInner.Dock = DockStyle.Fill;
            qInner.Padding = new Padding(0, 4, 0, 0);
            qGroup.Controls.Add(qInner);

            _qHintLabel.Text = "Редактируйте только внедиагональные клетки. Диагональ считается автоматически.";
            _qHintLabel.Dock = DockStyle.Top;
            _qHintLabel.Height = 22;
            _qHintLabel.Padding = new Padding(0, 0, 0, 0);
            qInner.Controls.Add(_qHintLabel);

            _qGrid.Dock = DockStyle.Fill;
            _qGrid.AllowUserToAddRows = false;
            _qGrid.AllowUserToDeleteRows = false;
            _qGrid.RowHeadersWidth = 150;
            _qGrid.AutoSizeColumnsMode = DataGridViewAutoSizeColumnsMode.Fill;
            _qGrid.SelectionMode = DataGridViewSelectionMode.CellSelect;
            _qGrid.MultiSelect = false;
            _qGrid.EditMode = DataGridViewEditMode.EditOnEnter;
            _qGrid.RowTemplate.Height = 34;
            _qGrid.ColumnHeadersHeight = 38;
            _qGrid.Font = new Font("Consolas", 11F, FontStyle.Regular, GraphicsUnit.Point);
            _qGrid.DefaultCellStyle.Alignment = DataGridViewContentAlignment.MiddleCenter;
            _qGrid.CellEndEdit += QGrid_CellEndEdit;
            _qGrid.DataError += QGrid_DataError;
            qInner.Controls.Add(_qGrid);
            qInner.Controls.SetChildIndex(_qGrid, 0);
            qInner.Controls.SetChildIndex(_qHintLabel, 1);

            GroupBox statsGroup = new GroupBox();
            statsGroup.Text = "Статистика";
            statsGroup.Dock = DockStyle.Fill;
            statsGroup.Padding = new Padding(8, 20, 8, 8);
            rightLayout.Controls.Add(statsGroup, 0, 2);

            _statsGrid.Dock = DockStyle.Fill;
            _statsGrid.AllowUserToAddRows = false;
            _statsGrid.AllowUserToDeleteRows = false;
            _statsGrid.ReadOnly = true;
            _statsGrid.RowHeadersVisible = false;
            _statsGrid.AutoSizeColumnsMode = DataGridViewAutoSizeColumnsMode.Fill;
            _statsGrid.SelectionMode = DataGridViewSelectionMode.FullRowSelect;
            _statsGrid.Font = new Font("Consolas", 9.2F, FontStyle.Regular, GraphicsUnit.Point);
            _statsGrid.RowTemplate.Height = 28;
            statsGroup.Controls.Add(_statsGrid);

            GroupBox distGroup = new GroupBox();
            distGroup.Text = "Сравнение распределений";
            distGroup.Dock = DockStyle.Fill;
            distGroup.Padding = new Padding(8, 20, 8, 8);
            rightLayout.Controls.Add(distGroup, 0, 3);

            _distPanel.Dock = DockStyle.Fill;
            _distPanel.BackColor = Color.White;
            _distPanel.BorderStyle = BorderStyle.FixedSingle;
            distGroup.Controls.Add(_distPanel);

            Panel csvPanel = new Panel();
            csvPanel.Dock = DockStyle.Bottom;
            csvPanel.Height = 100;
            csvPanel.Padding = new Padding(8, 4, 8, 4);
            rightPanel.Controls.Add(csvPanel);

            _csvLabel.Dock = DockStyle.Fill;
            _csvLabel.TextAlign = ContentAlignment.TopLeft;
            _csvLabel.Font = new Font("Consolas", 8.3F, FontStyle.Regular, GraphicsUnit.Point);
            _csvLabel.BorderStyle = BorderStyle.FixedSingle;
            csvPanel.Controls.Add(_csvLabel);

            CreateQGridColumns();
            CreateStatsGridColumns();
        }

        private void MainForm_Shown(object sender, EventArgs e)
        {
            _timelinePanel.Invalidate();
            _distPanel.Invalidate();
        }

        private void ConfigureTimer()
        {
            _timer.Interval = (int)_delayUpDown.Value;
            _timer.Tick += Timer_Tick;
        }

        private void CreateQGridColumns()
        {
            _qGrid.Columns.Clear();
            _qGrid.Rows.Clear();

            for (int i = 0; i < StateCount; i++)
            {
                DataGridViewTextBoxColumn col = new DataGridViewTextBoxColumn();
                col.HeaderText = _stateShortNames[i + 1];
                col.SortMode = DataGridViewColumnSortMode.NotSortable;
                _qGrid.Columns.Add(col);
            }

            for (int i = 0; i < StateCount; i++)
            {
                int rowIndex = _qGrid.Rows.Add();
                _qGrid.Rows[rowIndex].HeaderCell.Value = _stateNames[i + 1];
            }
        }

        private void CreateStatsGridColumns()
        {
            _statsGrid.Columns.Clear();
            _statsGrid.Rows.Clear();
            _statsGrid.Columns.Add(new DataGridViewTextBoxColumn { HeaderText = "Состояние" });
            _statsGrid.Columns.Add(new DataGridViewTextBoxColumn { HeaderText = "Count" });
            _statsGrid.Columns.Add(new DataGridViewTextBoxColumn { HeaderText = "Empirical" });
            _statsGrid.Columns.Add(new DataGridViewTextBoxColumn { HeaderText = "Stationary" });
            _statsGrid.Columns.Add(new DataGridViewTextBoxColumn { HeaderText = "|Diff|" });
        }

        private void LoadDefaults()
        {
            _updatingGrid = true;
            try
            {
                for (int i = 0; i < StateCount; i++)
                {
                    for (int j = 0; j < StateCount; j++)
                        _q[i, j] = _defaultQ[i, j];
                }

                WriteQToGrid();
            }
            finally
            {
                _updatingGrid = false;
            }

            _reportBox.Clear();
            _statsGrid.Rows.Clear();
            _timelinePanel.SetData(new List<int>(), new List<DateTime>(), 0, false);
            _distPanel.SetData(_empirical, _stationary, 0);
            UpdateCurrentLabel("Q загружена");
            UpdateCsvLabel();
        }

        private void WriteQToGrid()
        {
            for (int i = 0; i < StateCount; i++)
            {
                for (int j = 0; j < StateCount; j++)
                {
                    if (i == j)
                    {
                        _qGrid.Rows[i].Cells[j].Value = FormatNumber(-RowOffDiagonalSum(i));
                        _qGrid.Rows[i].Cells[j].ReadOnly = true;
                        _qGrid.Rows[i].Cells[j].Style.BackColor = Color.MistyRose;
                        _qGrid.Rows[i].Cells[j].Style.ForeColor = Color.DarkRed;
                    }
                    else
                    {
                        _qGrid.Rows[i].Cells[j].Value = FormatNumber(_q[i, j]);
                        _qGrid.Rows[i].Cells[j].ReadOnly = false;
                        _qGrid.Rows[i].Cells[j].Style.BackColor = Color.White;
                        _qGrid.Rows[i].Cells[j].Style.ForeColor = Color.Black;
                    }
                }
            }
        }

        private void UpdateDiagonalCells()
        {
            for (int i = 0; i < StateCount; i++)
            {
                _qGrid.Rows[i].Cells[i].Value = FormatNumber(-RowOffDiagonalSum(i));
            }
        }

        private void QGrid_CellEndEdit(object sender, DataGridViewCellEventArgs e)
        {
            if (_updatingGrid)
                return;

            string error;
            if (!TryReadQFromGrid(out error))
            {
                MessageBox.Show(error, "Ошибка в Q", MessageBoxButtons.OK, MessageBoxIcon.Error);
                WriteQToGrid();
                return;
            }

            UpdateDiagonalCells();
            UpdateCurrentLabel("Q обновлена");
        }

        private void QGrid_DataError(object sender, DataGridViewDataErrorEventArgs e)
        {
            e.ThrowException = false;
        }

        private bool TryReadQFromGrid(out string error)
        {
            error = string.Empty;

            for (int i = 0; i < StateCount; i++)
            {
                double rowSum = 0;
                for (int j = 0; j < StateCount; j++)
                {
                    if (i == j)
                        continue;

                    string cellText = _qGrid.Rows[i].Cells[j].Value == null ? string.Empty : _qGrid.Rows[i].Cells[j].Value.ToString().Trim();
                    if (string.IsNullOrWhiteSpace(cellText))
                    {
                        error = "Пустая ячейка: строка " + _stateNames[i + 1] + ", столбец " + _stateShortNames[j + 1];
                        return false;
                    }

                    double value;
                    if (!TryParseDouble(cellText, out value) || double.IsNaN(value) || double.IsInfinity(value))
                    {
                        error = "Не число: строка " + _stateNames[i + 1] + ", столбец " + _stateShortNames[j + 1];
                        return false;
                    }

                    if (value < 0)
                    {
                        error = "Интенсивность не может быть отрицательной: " + _stateNames[i + 1] + " -> " + _stateShortNames[j + 1];
                        return false;
                    }

                    _q[i, j] = value;
                    rowSum += value;
                }

                for (int j = 0; j < StateCount; j++)
                {
                    if (i == j)
                        _q[i, j] = -rowSum;
                }
            }

            return true;
        }

        private double RowOffDiagonalSum(int row)
        {
            double sum = 0;
            for (int j = 0; j < StateCount; j++)
            {
                if (j == row)
                    continue;
                sum += _q[row, j];
            }
            return sum;
        }

        private static bool TryParseDouble(string text, out double value)
        {
            return double.TryParse(text, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture, out value)
                   || double.TryParse(text, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.CurrentCulture, out value);
        }

        private static string FormatNumber(double value)
        {
            return value.ToString("0.###", CultureInfo.InvariantCulture);
        }

        private void StartButton_Click(object sender, EventArgs e)
        {
            StartSimulation();
        }

        private void ResetButton_Click(object sender, EventArgs e)
        {
            ResetAll();
        }

        private void DefaultsButton_Click(object sender, EventArgs e)
        {
            LoadDefaults();
        }

        private void StartSimulation()
        {
            if (_running)
                return;

            string validationError;
            if (!TryReadQFromGrid(out validationError))
            {
                MessageBox.Show(validationError, "Нельзя запустить", MessageBoxButtons.OK, MessageBoxIcon.Error);
                WriteQToGrid();
                return;
            }

            _targetHours = (int)_durationCombo.SelectedItem;
            int monthIndex = _monthCombo.SelectedIndex + 1;
            _startDate = new DateTime(DateTime.Now.Year, monthIndex, 1, 0, 0, 0);
            _currentState = 1;
            _currentDayIndex = 0;
            _running = true;

            _dailyStates.Clear();
            _dailyDates.Clear();
            Array.Clear(_counts, 0, _counts.Length);
            Array.Clear(_empirical, 0, _empirical.Length);
            _stationary = new double[StateCount];
            _dailyCsvPath = string.Empty;
            _summaryCsvPath = string.Empty;
            _reportTxtPath = string.Empty;

            _startButton.Enabled = false;
            _defaultsButton.Enabled = false;
            _monthCombo.Enabled = false;
            _durationCombo.Enabled = false;
            _delayUpDown.Enabled = false;
            _qGrid.Enabled = false;

            _timer.Interval = (int)_delayUpDown.Value;
            _timer.Start();

            UpdateCurrentLabel("Старт: " + _stateNames[_currentState] + " | " + _startDate.ToString("MMMM yyyy HH:mm", CultureInfo.CurrentCulture));
            _timelinePanel.SetData(_dailyStates, _dailyDates, _currentDayIndex, true);
            _timelinePanel.Invalidate();
            _distPanel.SetData(_empirical, _stationary, 0);
            _distPanel.Invalidate();
            _statsGrid.Rows.Clear();
            _reportBox.Clear();
        }

        private void Timer_Tick(object sender, EventArgs e)
        {
            if (!_running)
                return;

            if (_currentDayIndex >= _targetHours)
            {
                FinishSimulation();
                return;
            }

            DateTime date = _startDate.AddHours(_currentDayIndex);
            _dailyStates.Add(_currentState);
            _dailyDates.Add(date);
            _counts[_currentState]++;

            _currentDayIndex++;
            UpdateCurrentLabel(string.Format(CultureInfo.CurrentCulture, "Час {0}/{1}: {2} ({3:dd.MM.yyyy HH:00})", _currentDayIndex, _targetHours, _stateNames[_currentState], date));
            _timelinePanel.SetData(_dailyStates, _dailyDates, _currentDayIndex, true);
            _timelinePanel.Invalidate();

            _currentState = Step(_currentState);

            if (_currentDayIndex >= _targetHours)
                FinishSimulation();
        }

        private int Step(int currentState)
        {
            double lambda = 0;
            for (int j = 0; j < StateCount; j++)
            {
                if (j == currentState - 1)
                    continue;
                lambda += _q[currentState - 1, j];
            }

            if (lambda <= 0)
                return currentState;

            double pChange = 1.0 - Math.Exp(-lambda);
            if (_rng.NextDouble() >= pChange)
                return currentState;

            double r = _rng.NextDouble() * lambda;
            double cumulative = 0;
            for (int j = 0; j < StateCount; j++)
            {
                if (j == currentState - 1)
                    continue;

                cumulative += _q[currentState - 1, j];
                if (r <= cumulative)
                    return j + 1;
            }

            return currentState;
        }

        private void FinishSimulation()
        {
            if (!_running)
                return;

            _timer.Stop();
            _running = false;

            int hours = _dailyStates.Count;
            if (hours <= 0)
            {
                RestoreUiAfterRun();
                return;
            }

            for (int i = 1; i <= StateCount; i++)
                _empirical[i] = _counts[i] / (double)hours;

            try
            {
                _stationary = SolveStationaryDistribution(_q);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Не удалось вычислить стационарное распределение: " + ex.Message, "Ошибка", MessageBoxButtons.OK, MessageBoxIcon.Error);
                _stationary = new double[StateCount];
            }

            FillStatsGrid();
            BuildReportText(hours);
            WriteCsvFiles(hours);
            UpdateCsvLabel();

            _timelinePanel.SetData(_dailyStates, _dailyDates, _currentDayIndex, false);
            _timelinePanel.Invalidate();
            _distPanel.SetData(_empirical, _stationary, hours);
            _distPanel.Invalidate();
            UpdateCurrentLabel("Готово. CSV сохранён в output");

            RestoreUiAfterRun();
        }

        private void RestoreUiAfterRun()
        {
            _startButton.Enabled = true;
            _defaultsButton.Enabled = true;
            _monthCombo.Enabled = true;
            _durationCombo.Enabled = true;
            _delayUpDown.Enabled = true;
            _qGrid.Enabled = true;
        }

        private void FillStatsGrid()
        {
            _statsGrid.Rows.Clear();
            for (int i = 1; i <= StateCount; i++)
            {
                _statsGrid.Rows.Add(
                    _stateNames[i],
                    _counts[i],
                    _empirical[i].ToString("F6", CultureInfo.InvariantCulture),
                    _stationary[i - 1].ToString("F6", CultureInfo.InvariantCulture),
                    Math.Abs(_empirical[i] - _stationary[i - 1]).ToString("F6", CultureInfo.InvariantCulture));
            }
        }

        private void BuildReportText(int hours)
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("Марковская модель погоды");
            sb.AppendLine();
            sb.AppendLine("Стартовый месяц: " + _startDate.ToString("MMMM yyyy HH:mm", CultureInfo.CurrentCulture));
            sb.AppendLine("Часов: " + hours);
            sb.AppendLine();
            sb.AppendLine("Состояния:");
            sb.AppendLine("1 — Ясно");
            sb.AppendLine("2 — Переменная облачность");
            sb.AppendLine("3 — Облачно");
            sb.AppendLine("4 — Пасмурно");
            sb.AppendLine("5 — Дождь");
            sb.AppendLine();
            sb.AppendLine("Число попаданий:");
            sb.AppendLine("Ясно:                 " + _counts[1]);
            sb.AppendLine("Перем. облачность:    " + _counts[2]);
            sb.AppendLine("Облачно:              " + _counts[3]);
            sb.AppendLine("Пасмурно:             " + _counts[4]);
            sb.AppendLine("Дождь:                " + _counts[5]);
            sb.AppendLine();
            sb.AppendLine("Эмпирическое распределение:");
            sb.AppendLine("Ясно:                 " + _empirical[1].ToString("F6", CultureInfo.InvariantCulture));
            sb.AppendLine("Перем. облачность:    " + _empirical[2].ToString("F6", CultureInfo.InvariantCulture));
            sb.AppendLine("Облачно:              " + _empirical[3].ToString("F6", CultureInfo.InvariantCulture));
            sb.AppendLine("Пасмурно:             " + _empirical[4].ToString("F6", CultureInfo.InvariantCulture));
            sb.AppendLine("Дождь:                " + _empirical[5].ToString("F6", CultureInfo.InvariantCulture));
            sb.AppendLine();
            sb.AppendLine("Теоретическое стационарное распределение:");
            sb.AppendLine("Ясно:                 " + _stationary[0].ToString("F6", CultureInfo.InvariantCulture));
            sb.AppendLine("Перем. облачность:    " + _stationary[1].ToString("F6", CultureInfo.InvariantCulture));
            sb.AppendLine("Облачно:              " + _stationary[2].ToString("F6", CultureInfo.InvariantCulture));
            sb.AppendLine("Пасмурно:             " + _stationary[3].ToString("F6", CultureInfo.InvariantCulture));
            sb.AppendLine("Дождь:                " + _stationary[4].ToString("F6", CultureInfo.InvariantCulture));

            _reportBox.Text = sb.ToString();
        }

        private void WriteCsvFiles(int hours)
        {
            string outDir = Path.Combine(AppContext.BaseDirectory, "output");
            Directory.CreateDirectory(outDir);

            string stamp = DateTime.Now.ToString("yyyyMMdd_HHmmss");
            _dailyCsvPath = Path.Combine(outDir, "weather_hourly_" + stamp + ".csv");
            _summaryCsvPath = Path.Combine(outDir, "weather_summary_" + stamp + ".csv");
            _reportTxtPath = Path.Combine(outDir, "weather_report_" + stamp + ".txt");

            List<string> dailyLines = new List<string>();
            dailyLines.Add("hour,datetime,state_id,state_name");
            for (int i = 0; i < _dailyStates.Count; i++)
            {
                dailyLines.Add(string.Format(CultureInfo.InvariantCulture, "{0},{1:yyyy-MM-dd HH:mm},{2},{3}", i + 1, _dailyDates[i], _dailyStates[i], _stateNames[_dailyStates[i]]));
            }
            File.WriteAllLines(_dailyCsvPath, dailyLines.ToArray());

            List<string> summaryLines = new List<string>();
            summaryLines.Add("metric,state,values");
            summaryLines.Add(string.Format(CultureInfo.InvariantCulture, "hours,,{0}", hours));
            summaryLines.Add(string.Format(CultureInfo.InvariantCulture, "start_datetime,,{0:yyyy-MM-dd HH:mm}", _startDate));
            for (int i = 1; i <= StateCount; i++)
                summaryLines.Add(string.Format(CultureInfo.InvariantCulture, "count,{0},{1}", _stateNames[i], _counts[i]));
            for (int i = 1; i <= StateCount; i++)
                summaryLines.Add(string.Format(CultureInfo.InvariantCulture, "empirical,{0},{1}", _stateNames[i], _empirical[i].ToString("F6", CultureInfo.InvariantCulture)));
            for (int i = 1; i <= StateCount; i++)
                summaryLines.Add(string.Format(CultureInfo.InvariantCulture, "stationary,{0},{1}", _stateNames[i], _stationary[i - 1].ToString("F6", CultureInfo.InvariantCulture)));
            File.WriteAllLines(_summaryCsvPath, summaryLines.ToArray());
            File.WriteAllText(_reportTxtPath, _reportBox.Text);
        }

        private void UpdateCsvLabel()
        {
            if (string.IsNullOrWhiteSpace(_dailyCsvPath))
            {
                _csvLabel.Text = "CSV: пока нет";
                return;
            }

            _csvLabel.Text =
                "CSV:\r\n" +
                "daily:   " + _dailyCsvPath + "\r\n" +
                "summary: " + _summaryCsvPath + "\r\n" +
                "report:  " + _reportTxtPath;
        }

        private void ResetAll()
        {
            if (_running)
            {
                _timer.Stop();
                _running = false;
            }

            _currentDayIndex = 0;
            _dailyStates.Clear();
            _dailyDates.Clear();
            Array.Clear(_counts, 0, _counts.Length);
            Array.Clear(_empirical, 0, _empirical.Length);
            _stationary = new double[StateCount];
            _reportBox.Clear();
            _statsGrid.Rows.Clear();
            _dailyCsvPath = string.Empty;
            _summaryCsvPath = string.Empty;
            _reportTxtPath = string.Empty;
            UpdateCurrentLabel("Сброшено");
            UpdateCsvLabel();
            _timelinePanel.SetData(_dailyStates, _dailyDates, 0, false);
            _timelinePanel.Invalidate();
            _distPanel.SetData(_empirical, _stationary, 0);
            _distPanel.Invalidate();
            RestoreUiAfterRun();
        }

        private double[] SolveStationaryDistribution(double[,] q)
        {
            int n = StateCount;
            double[,] a = new double[n, n];
            double[] b = new double[n];

            for (int i = 0; i < n - 1; i++)
            {
                for (int j = 0; j < n; j++)
                    a[i, j] = q[j, i];
                b[i] = 0;
            }

            for (int j = 0; j < n; j++)
                a[n - 1, j] = 1;
            b[n - 1] = 1;

            return GaussianSolve(a, b);
        }

        private static double[] GaussianSolve(double[,] a, double[] b)
        {
            int n = b.Length;
            double[,] m = new double[n, n + 1];

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    m[i, j] = a[i, j];
                m[i, n] = b[i];
            }

            for (int col = 0; col < n; col++)
            {
                int pivot = col;
                for (int row = col + 1; row < n; row++)
                {
                    if (Math.Abs(m[row, col]) > Math.Abs(m[pivot, col]))
                        pivot = row;
                }

                if (Math.Abs(m[pivot, col]) < 1e-12)
                    throw new InvalidOperationException("Система вырождена.");

                if (pivot != col)
                    SwapRows(m, pivot, col);

                double div = m[col, col];
                for (int j = col; j <= n; j++)
                    m[col, j] /= div;

                for (int row = 0; row < n; row++)
                {
                    if (row == col)
                        continue;

                    double factor = m[row, col];
                    for (int j = col; j <= n; j++)
                        m[row, j] -= factor * m[col, j];
                }
            }

            double[] x = new double[n];
            for (int i = 0; i < n; i++)
                x[i] = m[i, n];

            double sum = 0;
            for (int i = 0; i < n; i++)
                sum += x[i];

            if (Math.Abs(sum) > 1e-12)
            {
                for (int i = 0; i < n; i++)
                    x[i] /= sum;
            }

            return x;
        }

        private static void SwapRows(double[,] m, int r1, int r2)
        {
            int cols = m.GetLength(1);
            for (int j = 0; j < cols; j++)
            {
                double tmp = m[r1, j];
                m[r1, j] = m[r2, j];
                m[r2, j] = tmp;
            }
        }

        private void UpdateCurrentLabel(string text)
        {
            _currentLabel.Text = text;
        }
    }

    public sealed class TimelinePanel : Panel
    {
        private List<int> _states = new List<int>();
        private List<DateTime> _dates = new List<DateTime>();
        private int _currentIndex;
        private bool _running;

        public TimelinePanel()
        {
            DoubleBuffered = true;
            ResizeRedraw = true;
            AutoScroll = true;
        }

        public void SetData(List<int> states, List<DateTime> dates, int currentIndex, bool running)
        {
            _states = states ?? new List<int>();
            _dates = dates ?? new List<DateTime>();
            _currentIndex = currentIndex;
            _running = running;
            Invalidate();
        }

        protected override void OnPaint(PaintEventArgs e)
        {
            base.OnPaint(e);
            e.Graphics.Clear(Color.White);
            e.Graphics.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;

            const int cellW = 120;
            const int cellH = 80;
            const int gap = 10;
            const int left = 12;
            const int top = 38;
            const int headerH = 40;
            const int cols = 7;

            int total = _states.Count;
            int rows = Math.Max(1, (int)Math.Ceiling(total / (double)cols));
            AutoScrollMinSize = new Size(left * 2 + cols * cellW + (cols - 1) * gap, top * 2 + headerH + rows * cellH + (rows - 1) * gap);

            e.Graphics.TranslateTransform(AutoScrollPosition.X, AutoScrollPosition.Y);

            Font titleFont = new Font(Font.FontFamily, 12F, FontStyle.Bold);
            Font dayFont = new Font(Font.FontFamily, 11F, FontStyle.Bold);
            Font stateFont = new Font(Font.FontFamily, 8.5F, FontStyle.Regular);
            Font dateFont = new Font(Font.FontFamily, 8F, FontStyle.Regular);
            Pen borderPen = new Pen(Color.DimGray, 1f);

            e.Graphics.DrawString("Визуализация по часам", titleFont, Brushes.Black, left, 0);

            if (total == 0)
            {
                e.Graphics.DrawString("Нажмите «Старт», чтобы начать моделирование.", Font, Brushes.Black, left, 34);
                titleFont.Dispose();
                dayFont.Dispose();
                stateFont.Dispose();
                dateFont.Dispose();
                borderPen.Dispose();
                return;
            }

            for (int i = 0; i < total; i++)
            {
                int row = i / cols;
                int col = i % cols;
                int x = left + col * (cellW + gap);
                int y = top + headerH + row * (cellH + gap);
                Rectangle rect = new Rectangle(x, y, cellW, cellH);

                int state = _states[i];
                Brush fill = new SolidBrush(GetColor(state));
                e.Graphics.FillRectangle(fill, rect);
                e.Graphics.DrawRectangle(borderPen, rect);

                if (_running && i == Math.Min(_currentIndex, total - 1))
                {
                    Pen highlightPen = new Pen(Color.OrangeRed, 3f);
                    e.Graphics.DrawRectangle(highlightPen, rect);
                    highlightPen.Dispose();
                }

                e.Graphics.DrawString((i + 1).ToString(), dayFont, Brushes.Black, x + 8, y + 5);
                e.Graphics.DrawString(GetShortName(state), stateFont, Brushes.Black, x + 8, y + 28);

                if (i < _dates.Count)
                    e.Graphics.DrawString(_dates[i].ToString("dd.MM HH:mm", CultureInfo.InvariantCulture), dateFont, Brushes.Black, x + 8, y + 48);

                fill.Dispose();
            }

            titleFont.Dispose();
            dayFont.Dispose();
            stateFont.Dispose();
            dateFont.Dispose();
            borderPen.Dispose();
        }

        private static Color GetColor(int state)
        {
            switch (state)
            {
                case 1: return Color.LightGoldenrodYellow;
                case 2: return Color.LightSkyBlue;
                case 3: return Color.Gainsboro;
                case 4: return Color.Silver;
                case 5: return Color.SteelBlue;
                default: return Color.White;
            }
        }

        private static string GetShortName(int state)
        {
            switch (state)
            {
                case 1: return "Ясно";
                case 2: return "П.о.";
                case 3: return "Обл";
                case 4: return "Пасм";
                case 5: return "Дождь";
                default: return "";
            }
        }
    }

    public sealed class DistributionPanel : Panel
    {
        private double[] _empirical = new double[StateCountHolder.Value];
        private double[] _stationary = new double[StateCountHolder.Value];
        private int _hours;

        public DistributionPanel()
        {
            DoubleBuffered = true;
            ResizeRedraw = true;
        }

        public void SetData(double[] empirical, double[] stationary, int hours)
        {
            _empirical = empirical ?? new double[StateCountHolder.Value];
            _stationary = stationary ?? new double[StateCountHolder.Value];
            _hours = hours;
            Invalidate();
        }

        protected override void OnPaint(PaintEventArgs e)
        {
            base.OnPaint(e);
            e.Graphics.Clear(Color.White);
            e.Graphics.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;

            Font headerFont = new Font(Font.FontFamily, 10F, FontStyle.Bold);
            Font smallFont = new Font(Font.FontFamily, 8F, FontStyle.Regular);
            Pen pen = new Pen(Color.DimGray, 1f);

            const int left = 10;
            const int top = 8;
            const int labelW = 120;
            const int barW = 230;
            const int rowH = 32;

            e.Graphics.DrawString("Распределение", headerFont, Brushes.Black, left, top);
            e.Graphics.DrawString("Синий — эмпирическое, оранжевый — стационарное", smallFont, Brushes.Black, left, top + 18);

            if (_hours <= 0)
            {
                e.Graphics.DrawString("График появится после завершения моделирования.", Font, Brushes.Black, left, top + 44);
                headerFont.Dispose();
                smallFont.Dispose();
                pen.Dispose();
                return;
            }

            string[] names = { "", "Ясно", "П. облач.", "Облачно", "Пасмурно", "Дождь" };
            int startY = top + 44;

            for (int i = 1; i <= 5; i++)
            {
                int y = startY + (i - 1) * rowH;
                e.Graphics.DrawString(names[i], smallFont, Brushes.Black, left, y + 9);

                int empW = (int)(barW * _empirical[i]);
                int statW = (int)(barW * _stationary[i - 1]);

                Rectangle empRect = new Rectangle(left + labelW, y + 5, Math.Max(1, empW), 10);
                Rectangle statRect = new Rectangle(left + labelW, y + 19, Math.Max(1, statW), 10);

                Brush empBrush = new SolidBrush(Color.SteelBlue);
                Brush statBrush = new SolidBrush(Color.DarkOrange);
                e.Graphics.FillRectangle(empBrush, empRect);
                e.Graphics.FillRectangle(statBrush, statRect);
                e.Graphics.DrawRectangle(pen, left + labelW, y + 5, barW, 10);
                e.Graphics.DrawRectangle(pen, left + labelW, y + 19, barW, 10);

                e.Graphics.DrawString(_empirical[i].ToString("F3", CultureInfo.InvariantCulture), smallFont, Brushes.Black, left + labelW + barW + 8, y + 2);
                e.Graphics.DrawString(_stationary[i - 1].ToString("F3", CultureInfo.InvariantCulture), smallFont, Brushes.Black, left + labelW + barW + 8, y + 16);

                empBrush.Dispose();
                statBrush.Dispose();
            }

            headerFont.Dispose();
            smallFont.Dispose();
            pen.Dispose();
        }
    }

    internal static class StateCountHolder
    {
        public const int Value = 5;
    }
}
