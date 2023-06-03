import os
import pathlib

from random import randint

from constellation import Constellation

# Создание объекта типа Constellation
constellation = Constellation()

# Инициализация параметрами группировки Starlink из конфига
try:
	parent_dir_path = pathlib.Path(__file__).parent.parent.resolve()
	example_file_path = os.path.join(parent_dir_path, "constellationsTest.json")
	constellation.load_from_config(example_file_path, "Starlink")
except Exception as ex:
	exit(f"Ошибка загрузки конфига: {ex}")

# Вычисление элементов орбиты для всех КА в начальный момент времени
constellation.initialize_state()

# Определение точек на оси времени, в которые будут производиться расчёты
points_number = 1002
epoch_list = list(range(points_number))

# Расчёт положений всех КА в заданные моменты времени
constellation.propagate_j2(epoch_list)

# Координаты случайного КА
sat_idx = randint(0, constellation.total_sat_count - 1)
epoch_idx = randint(0, points_number - 1)

print(f"Положение КА-{sat_idx} на эпоху {epoch_list[epoch_idx]}:")
print(constellation.state_eci[sat_idx, :, epoch_idx])
