 using System;
using System.Drawing;
using System.Windows.Forms;

namespace RandomEventsGui
{
    internal static class Program
    {
        [STAThread]
        static void Main()
        {
            ApplicationConfiguration.Initialize();
            Application.Run(new MainForm());
        }
    }

    public sealed class MainForm : Form
    {
        private readonly LcgRandom _lcg;

        private readonly Label _seedLabel = new();
        private readonly Label _lcgStateLabel = new();
        private readonly Label _yesNoResultLabel = new();
        private readonly Label _eightBallResultLabel = new();
        private readonly Label _lastRangeLabel = new();
        private readonly Label _chanceValueLabel = new();

        private readonly Button _yesNoButton = new();
        private readonly Button _eightBallButton = new();
        private readonly Button _newSeedButton = new();
        private readonly TrackBar _chanceTrackBar = new();

        private struct EightBallItem
        {
            public string Text;
            public double Chance;

            public EightBallItem(string text, double chance)
            {
                Text = text;
                Chance = chance;
            }
        }

        private static readonly EightBallItem[] EightBallAnswers =
        {
            new("Без сомнений — да.", 0.0625),
            new("Да.", 0.0625),
            new("Скорее да.", 0.0625),
            new("Похоже, что да.", 0.625),
            new("Ответ туманный, но шанс есть.", 0.0625),
            new("Спроси позже.", 0.0625),
            new("Сейчас не время.", 0.0625),
            new("Скорее нет.", 0.0625),
            new("Нет.", 0.0625),
            new("Точно нет.", 0.0625),
            new("Очень маловероятно.", 0.0625),
            new("Знаки говорят: да.", 0.0625),
            new("Знаки говорят: нет.", 0.0625),
            new("Да, но с подвохом.", 0.0625),
            new("Нет, и это лучше принять.", 0.0625),
            new("Возможно.", 0.0625),
        };
        
        public MainForm()
        {
            Text = "Моделирование случайных событий";
            StartPosition = FormStartPosition.CenterScreen;
            ClientSize = new Size(980, 520);
            Font = new Font("Segoe UI", 10F, FontStyle.Regular, GraphicsUnit.Point);
            MinimumSize = new Size(900, 480);

            ulong seed = (ulong)DateTime.Now.Ticks;
            _lcg = new LcgRandom(seed);

            BuildUi();
            RefreshSeedInfo();
            UpdateCurrentState();
            UpdateChanceLabel();
        }

        private void BuildUi()
        {
            var title = new Label
            {
                Text = "Моделирование случайных событий",
                AutoSize = true,
                Font = new Font(Font, FontStyle.Bold),
                Location = new Point(18, 16)
            };

            var infoPanel = new Panel
            {
                Location = new Point(18, 84),
                Size = new Size(ClientSize.Width - 36, 62),
                Anchor = AnchorStyles.Left | AnchorStyles.Top | AnchorStyles.Right,
                BorderStyle = BorderStyle.FixedSingle
            };

            var seedTitle = new Label
            {
                Text = "Seed:",
                AutoSize = true,
                Location = new Point(12, 10)
            };

            _seedLabel.AutoSize = true;
            _seedLabel.Location = new Point(62, 10);

            var stateTitle = new Label
            {
                Text = "Текущее состояние LCG:",
                AutoSize = true,
                Location = new Point(12, 33)
            };

            _lcgStateLabel.AutoSize = true;
            _lcgStateLabel.Location = new Point(177, 33);

            _newSeedButton.Text = "Новый seed";
            _newSeedButton.Size = new Size(120, 30);
            _newSeedButton.Anchor = AnchorStyles.Top | AnchorStyles.Right;
            _newSeedButton.Location = new Point(infoPanel.Width - 132, 15);
            _newSeedButton.Click += (_, _) =>
            {
                _lcg.Reseed((ulong)DateTime.Now.Ticks);
                RefreshSeedInfo();
                UpdateCurrentState();
                _yesNoResultLabel.Text = string.Empty;
                _eightBallResultLabel.Text = string.Empty;
                _lastRangeLabel.Text = string.Empty;
            };
            infoPanel.Controls.Add(_newSeedButton);
            infoPanel.Controls.Add(seedTitle);
            infoPanel.Controls.Add(_seedLabel);
            infoPanel.Controls.Add(stateTitle);
            infoPanel.Controls.Add(_lcgStateLabel);

            var left = CreateYesNoPanel();
            var right = CreateEightBallPanel();

            left.Location = new Point(18, 164);
            right.Location = new Point(18 + (ClientSize.Width - 54) / 2, 164);

            left.Anchor = AnchorStyles.Top | AnchorStyles.Bottom | AnchorStyles.Left;
            right.Anchor = AnchorStyles.Top | AnchorStyles.Bottom | AnchorStyles.Left | AnchorStyles.Right;

            Controls.Add(title);
            Controls.Add(infoPanel);
            Controls.Add(left);
            Controls.Add(right);

            Resize += (_, _) =>
            {
                infoPanel.Width = ClientSize.Width - 36;
                _newSeedButton.Left = infoPanel.Width - _newSeedButton.Width - 12;
                left.Width = (ClientSize.Width - 46) / 2;
                right.Width = (ClientSize.Width - 46) / 2;
                right.Left = 28 + left.Width;
                left.Height = ClientSize.Height - 184;
                right.Height = ClientSize.Height - 184;
                ResizeYesNoLayout();
                ResizeEightBallLayout();
            };
        }

        private GroupBox CreateYesNoPanel()
        {
            var box = new GroupBox
            {
                Text = "Часть 1. Скажи да или нет",
                Size = new Size((ClientSize.Width - 46) / 2, ClientSize.Height - 184)
            };

            _yesNoButton.Text = "Сгенерировать ответ";
            _yesNoButton.Size = new Size(190, 40);
            _yesNoButton.Location = new Point(12, 24);
            _yesNoButton.Click += (_, _) =>
            {
                double u = _lcg.NextDouble();
                int yesChance = _chanceTrackBar.Value;
                string answer = u < yesChance / 100.0 ? "Да" : "Нет";
                _yesNoResultLabel.Text = answer;
                _lastRangeLabel.Text = $"LCG.NextDouble() = {u:F6} → порог = {yesChance}%";
                UpdateCurrentState();
            };

            _yesNoResultLabel.AutoSize = false;
            _yesNoResultLabel.TextAlign = ContentAlignment.MiddleCenter;
            _yesNoResultLabel.Font = new Font(Font.FontFamily, 20F, FontStyle.Bold);
            _yesNoResultLabel.Location = new Point(12, 86);
            _yesNoResultLabel.Size = new Size(box.Width - 24, 60);
            _yesNoResultLabel.BorderStyle = BorderStyle.FixedSingle;

            var chanceTitle = new Label
            {
                Text = "Шанс ответа «Да»:",
                AutoSize = true,
                Location = new Point(12, 162)
            };

            _chanceTrackBar.Minimum = 0;
            _chanceTrackBar.Maximum = 100;
            _chanceTrackBar.Value = 50;
            _chanceTrackBar.TickFrequency = 10;
            _chanceTrackBar.SmallChange = 1;
            _chanceTrackBar.LargeChange = 10;
            _chanceTrackBar.Location = new Point(12, 188);
            _chanceTrackBar.Size = new Size(box.Width - 130, 45);
            _chanceTrackBar.Scroll += (_, _) => UpdateChanceLabel();

            _chanceValueLabel.AutoSize = false;
            _chanceValueLabel.TextAlign = ContentAlignment.MiddleLeft;
            _chanceValueLabel.Location = new Point(box.Width - 96, 184);
            _chanceValueLabel.Size = new Size(84, 30);

            _lastRangeLabel.AutoSize = false;
            _lastRangeLabel.Location = new Point(12, 244);
            _lastRangeLabel.Size = new Size(box.Width - 24, 46);
            _lastRangeLabel.ForeColor = Color.DimGray;

            box.Controls.Add(_yesNoButton);
            box.Controls.Add(_yesNoResultLabel);
            box.Controls.Add(chanceTitle);
            box.Controls.Add(_chanceTrackBar);
            box.Controls.Add(_chanceValueLabel);
            box.Controls.Add(_lastRangeLabel);
            return box;
        }

        private GroupBox CreateEightBallPanel()
        {
            var box = new GroupBox
            {
                Text = "Часть 2. Magic 8-Ball",
                Size = new Size((ClientSize.Width - 46) / 2, ClientSize.Height - 184)
            };

            _eightBallButton.Text = "Встряхнуть шар";
            _eightBallButton.Size = new Size(170, 40);
            _eightBallButton.Location = new Point(12, 24);
            _eightBallButton.Click += (_, _) =>
            {
                double total = 0;
                foreach (var item in EightBallAnswers)
                    total += item.Chance;

                double u = _lcg.NextDouble() * total;

                double cumulative = 0.0;

                foreach (var item in EightBallAnswers)
                {
                    cumulative += item.Chance;
                    if (u < cumulative)
                    {
                        _eightBallResultLabel.Text = item.Text;
                        _lastRangeLabel.Text = $"u = {u:F6} / {total:F3}";
                        break;
                    }
                }

                UpdateCurrentState();
            };

            _eightBallResultLabel.AutoSize = false;
            _eightBallResultLabel.TextAlign = ContentAlignment.MiddleCenter;
            _eightBallResultLabel.Font = new Font(Font.FontFamily, 15F, FontStyle.Bold);
            _eightBallResultLabel.Location = new Point(12, 86);
            _eightBallResultLabel.Size = new Size(box.Width - 24, 102);
            _eightBallResultLabel.BorderStyle = BorderStyle.FixedSingle;

            var details = new Label
            {
                Text = "",
                AutoSize = false,
                Size = new Size(box.Width - 24, 30),
                Location = new Point(12, 198)
            };

            box.Controls.Add(_eightBallButton);
            box.Controls.Add(_eightBallResultLabel);
            box.Controls.Add(details);
            return box;
        }

        private void ResizeYesNoLayout()
        {
            if (Controls.Count == 0)
                return;

            if (Controls[2] is GroupBox box)
            {
                _yesNoResultLabel.Width = box.Width - 24;
                _lastRangeLabel.Width = box.Width - 24;
                _chanceTrackBar.Width = box.Width - 130;
                _chanceValueLabel.Left = box.Width - 96;
            }
        }

        private void ResizeEightBallLayout()
        {
            if (Controls.Count == 0)
                return;

            if (Controls[3] is GroupBox box)
            {
                _eightBallResultLabel.Width = box.Width - 24;
            }
        }

        private void UpdateChanceLabel()
        {
            _chanceValueLabel.Text = $"{_chanceTrackBar.Value}%";
        }

        private void RefreshSeedInfo()
        {
            _seedLabel.Text = _lcg.Seed.ToString();
        }

        private void UpdateCurrentState()
        {
            _lcgStateLabel.Text = _lcg.State.ToString();
        }
    }

    public sealed class LcgRandom
    {
        public const ulong M = 1UL << 63;
        private const ulong BETA = (1UL << 32) + 3;

        public ulong Seed { get; private set; }
        public ulong State { get; private set; }

        public LcgRandom(ulong seed)
        {
            Seed = seed;
            State = seed;
        }

        public void Reseed(ulong seed)
        {
            Seed = seed;
            State = seed;
        }

        public ulong NextRaw()
        {
            State = unchecked((BETA * State) % M);
            return State;
        }

        public double NextDouble()
        {
            return (double)NextRaw() / M;
        }
    }
}
