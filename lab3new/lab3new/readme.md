В данно работе представлен клеточный автомат, отражающий 
Структуры данных
Cell (класс — одна клетка)

Поля:

CellState State — состояние клетки.
int StageAge — сколько тиков клетка провела в текущей стадии (используется для переходов Seedling>Young>Old и для старения).
int BurnTicks — оставшиеся тики горения (если клетка в Burning).
double Fertility — плодородие / горючие остатки, диапазон 0..1. Растёт при горении, затем плавно убывает.
int TimeSinceBurn — тики после выгорания; используется для визуальных эффектов и расчётов восстановления.
int Elevation — высота / рельеф (целое).
CellState (enum)
Значения: Empty, Seedling, Young, Old, Burning, Water, Rock.
Основные поля и параметры модели (ForestModel)
int Width, int Height — размеры поля (Grid).
Cell[,] Grid — массив клеток, индексация Grid[y, x].
long Ticks — общий счётчик тиков модели.

Среда
double AmbientHumidity — базовая/целевaя влажность, задаётся ползунком (0..1 в коде, в UI отображается в процентах).
double Humidity — текущая (динамическая) влажность; изменяется от пожаров/дождя/релаксации.
double WindAngle — направление ветра в градусах (0..360). В коде windBar.Value.
double WindStrength — сила ветра (0..1).
Дождь (глобальное событие)
bool RainActive — активен ли дождь.
int RainStrength — 0..3 (0 — автоотключение дождя), сила дождя.
int RainTicksRemaining — сколько тиков дождя ещё идёт.

Параметры правил
double BaseIgnitionFromNeighbor = 0.28 — базовый вклад для вероятности возгорания от одного горящего соседа.
double BaseGrowthProb = 0.0012 — базовая вероятность прорастания пустой клетки (в тик).
double BurntFertilityBonus = 0.28 — сколько плодородия добавляется при выгорании клетки.
int SeedlingTicks = 20 — сколько тиков саженец > young.
int YoungTicks = 200 — сколько тиков young > old.
int MaxBurnTimeBase = 3, int MaxBurnExtra = 5 — параметры для вычисления длительности горения.

Старение
int OldMaxAgeTicks = 800 — возраст (в тиках) от которого начинает действовать шанс смерти у Old.
double OldDeathChancePerTick = 0.002 — шанс в каждом тике умереть у Old после достижения OldMaxAgeTicks.

Спонтанность
double spontaneousBase = 4e-7 — коэффициент фоновых спонтанных поджогов (молния/самовозгорание).
Влажность / вспомогательные константы
const double fireImpact = 0.28 — коэффициент, как сильно коллективный пожар снижает влажность.
const double rainHumidityPerStrength = 0.02 — увеличение Humidity в тик при дожде на 1 силу (в тиках).
const double ambientRelaxRate = 0.006 — скорость, с которой Humidity стремится к AmbientHumidity.
const double evaporation = 0.00012 — базовое испарение (уменьшение влажности).

Генератор случайных чисел
Random rnd — используется повсюду; если нужен воспроизводимый результат — можно добавить возможность задать seed.
GUI — какие элементы управляют какими переменными
speedBar > задаёт timer.Interval (мс/тик). Значение отображается в speedLabel.
windBar > model.WindAngle (градусы).
windStrBar > model.WindStrength = windStrBar.Value / 100.0.
humidityBar > model.AmbientHumidity = humidityBar.Value / 100.0.
rainStrengthBar > максимальная сила дождя (передаётся в Step(rainMaxStrength)).
zoomBar > при паузе вызывает ApplyZoom (меняет cellSize и пересоздаёт bmp).
widthUpDown, heightUpDown > задают gridW, gridH и при применении вызывают model.Reset(gridW, gridH).
brushBox (опции: "seed/plant", "erase", "water", "rock", "ignite") — при клике на поле управляют изменением Cell.State.

Алгоритм шага — подробная логика (ForestModel.Step)
Метод Step(int rainMaxStrength) выполняет один тик симуляции. Порядок операций:
Ticks++.
Подготовка массивов willIgnite[h,w] и willExtinguish[h,w] — флаги, которые будут применены после расчётов.
TryAutoSpawnRain(rainMaxStrength) — попытка автоматически запустить дождь (с вероятностью spawnProb = 0.0007 + AmbientHumidity * 0.005, если дождя ещё нет и rainMaxStrength > 0). Если срабатывает, вызывает StartRainManual(...): RainActive = true, RainStrength = s, RainTicksRemaining = 50 + 30*s.

Перебор всех клеток:
Если клетка Burning:
burningCount++.
c.BurnTicks--.
Если дождь активен: вероятность мгновенного тушения extinguish = 0.20 + 0.28 * (RainStrength / 3.0). Если rnd < extinguish > пометить willExtinguish[y,x] = true. Если c.BurnTicks <= 0, тоже тушить.
Иначе (без дождя) — тушение только когда c.BurnTicks <= 0.
Если клетка — дерево (Seedling / Young / Old):
Собрать список горящих соседей (окрестность Мура, 8 клеток).
Если есть горящие соседи — вычислить вероятность, что клетка загорится от них:
Начать pNo = 1.0.
Для каждого горящего соседа nb вычислить вклад contrib:
contrib = BaseIgnitionFromNeighbor
          * stageFactor
          * windFactor
          * elevationFactor
          * (1 - Humidity^1.7)
          * rainFactor
где:
stageFactor = StageIgnitionFactor(c.State) — 0.45 для Seedling, 1.0 для Young, 1.6 для Old.
windFactor — вычисление по углу:
angleToCell = angle(from neighbor to this cell) в градусах;
windDiff = abs(normalize(angleToCell - WindAngle) - 180); затем cosv = cos(windDiff);
если cosv > 0 (ветер дует в сторону клетки) — windFactor = 1.0 + WindStrength * 2.6 * cosv;
если cosv <= 0 (ветер против) — windFactor = max(0.01, 1.0 + WindStrength * 1.8 * cosv).
Это делает передачу пожара легче по ветру и труднее против.
elevDiff = Grid[y,x].Elevation - Grid[nb.Y, nb.X].Elevation:
если elevDiff > 0 (клетка выше, то есть пожар идёт вверх), elevationFactor *= (1 + min(0.6, elevDiff * 0.08));
иначе elevationFactor *= (1 + elevDiff * 0.02) (понижение уменьшает/не сильно увеличивает).
contrib *= Math.Max(0.0, 1.0 - Math.Pow(Humidity, 1.7)) — влажность подавляет вероятность: экспоненциальное подавление (почему ^1.7 — эмпирический усилитель эффекта).
rainFactor = 1.0 обычно; если дождь активен, contrib *= Math.Max(0.01, 1.0 - 0.85 * (RainStrength / 3.0)) — дождь сильно уменьшает вклад.
Ограничить contrib в [0,1].
pNo *= (1.0 - contrib) (вероятность, что не загорится от данного соседа).

После цикла: pIgnite = 1.0 - pNo.
Умножить pIgnite *= (1.0 + Math.Min(2.0, c.Fertility * 12.0)). То есть плодородие (остатки после предыдущего пожара) увеличивают шанс (логика: там больше горючих остатков).
Если rnd < pIgnite > пометить willIgnite[y,x] = true.

Спонтанное возгорание (молния / фон):
pSpont = spontaneousBase * (1 - Humidity)^1.8 * (1 + 0.8*(ageFactor - 1.0))
где ageFactor = StageIgnitionFactor(c.State). Если rnd < pSpont > пометить willIgnite.

Если клетка Empty:
localTreeFrac = NeighborTreeFraction(x,y) — доля соседей, которые деревья.
growthProb = BaseGrowthProb * (1 + 4.0 * c.Fertility) * (1 + localTreeFrac) * (0.4 + 0.6 * min(1.0, Humidity + 0.1)).
Если клетка Water или Rock — growthProb = 0.
Если rnd < growthProb > Grid[y,x].State = Seedling и Grid[y,x].StageAge = 0 (и changed добавляется).

Применение willIgnite / willExtinguish:
Для всех клеток: если willIgnite[y,x] и клетка деревo > State = Burning, BurnTicks = BurnTime(c), changed.Add.
Если willExtinguish[y,x] и State == Burning > State = Empty, StageAge = 0, Fertility = min(1, Fertility + BurntFertilityBonus), TimeSinceBurn = 0, changed.Add.

Обновление стадий роста / старение / плодородие:
Если Seedling: StageAge++; если StageAge >= SeedlingTicks > State = Young, StageAge = 0.
Если Young: StageAge++; если StageAge >= YoungTicks > State = Old, StageAge = 0.
Если Old: StageAge++; если StageAge >= OldMaxAgeTicks и rnd < OldDeathChancePerTick > дерево умирает (State=Empty, Fertility += 0.12).
Если Burning: Fertility = min(1.0, Fertility + 0.006) (горение повышает плодородие).
Если TimeSinceBurn < 9999 и State != Burning: Fertility = max(0.0, Fertility - 0.0008) (постепенное снижение плодородия).

Эффект дождя (глобально):
Если RainActive:
Humidity += rainHumidityPerStrength * RainStrength (в коде rainHumidityPerStrength = 0.02).
ashClear = 3 + 3 * RainStrength и TimeSinceBurn инкрементируется этим числом для всех клеток (ускорение очистки пепла).
RainTicksRemaining--; если <= 0 — дождь выключается.

Обновление влажности:
towardAmbient = (AmbientHumidity - Humidity) * ambientRelaxRate.
burningFraction = burningCount / (w * h).
fireEffect = -burningFraction * fireImpact.
evap = -evaporation * Humidity.
Humidity += towardAmbient + fireEffect + evap.
Humidity = clamp(Humidity, 0, 1).
Возвращается список changed.Distinct().ToList() — координаты клеток, которые изменились, чтобы GUI перерисовывал только их.

Формулы / пояснения (коротко, но точно)
1. Вероятность горения от соседа (внутри цикла соседей)
contrib (на одного соседа) вычисляется как произведение факторов:
BaseIgnitionFromNeighbor (0.28) ? стадия (0.45/1.0/1.6) ? windFactor ? elevationFactor ? (1 - Humidity^1.7) ? rainFactor.

Итоговая вероятность загореться от всех соседей:
pNo = ?(1 - contrib_i) для всех горящих соседей i.
pIgnite = 1 - pNo.
Затем pIgnite *= (1 + min(2.0, c.Fertility * 12.0)) — плодородие увеличивает шанс.

Пример: при Humidity = 0.3, без ветра и рельефа, стадия Young, один горящий сосед:
contrib ? 0.28 * 1.0 * 1.0 * 1.0 * (1 - 0.3^1.7) ? 0.28 * (1 - 0.3^1.7).
0.3^1.7 ? 0.3^1.7 ? 0.3^1.7 ? 0.3^1.7 ? 0.3^1.7 (примерно 0.3^1.7 ? 0.144) ? 1 - 0.144 = 0.856 ? contrib ? 0.239.

2. Спонтанное возгорание
pSpont = spontaneousBase * (1 - Humidity)^1.8 * (1 + 0.8 * (ageFactor - 1.0)).
При высокой влажности pSpont резко падает.

3. BurnTime (время горения)
baseT = MaxBurnTimeBase + int(MaxBurnExtra * scaledStage) где scaledStage пропорционален стадии.
factor = 1 - 0.45 * Humidity
Если дождь: factor *= (1 - 0.45 * (RainStrength / 3.0)).
BurnTicks = max(1, int(baseT * factor)).

4. Рост деревьев
growthProb = BaseGrowthProb * (1 + 4.0 * c.Fertility) * (1 + localTreeFrac) * (0.4 + 0.6 * min(1.0, Humidity + 0.1)).
Чем больше плодородия и чем больше деревьев вокруг — тем выше шанс.
BaseGrowthProb = 0.0012 — медленный естественный рост.

5. Влажность
Динамика: Humidity += (AmbientHumidity - Humidity) * ambientRelaxRate - burningFraction * fireImpact - evaporation * Humidity.
дождь: Humidity += rainHumidityPerStrength * RainStrength.

6. Дождь
Автоген: spawnProb = 0.0007 + AmbientHumidity * 0.005.
Ручной StartRainManual(strength) даёт RainTicksRemaining = 50 + 30 * strength.

Пока идёт дождь:
Увеличение влажности + 0.02 * RainStrength в тик.
Тушение: вероятность тушения одной горящей клетки в тик ? 0.20 + 0.28 * (RainStrength/3.0).

Параметры по умолчанию (в коде)
(повтор для удобства - взяты значения прямо из ForestModel)
BaseIgnitionFromNeighbor = 0.28
BaseGrowthProb = 0.0012
BurntFertilityBonus = 0.28
SeedlingTicks = 20
YoungTicks = 200
MaxBurnTimeBase = 3
MaxBurnExtra = 5
OldMaxAgeTicks = 800
OldDeathChancePerTick = 0.002
spontaneousBase = 4e-7
fireImpact = 0.28
rainHumidityPerStrength = 0.02
ambientRelaxRate = 0.006
evaporation = 0.00012