import numpy as np
import os
import pathlib

from constellation import Constellation
from consts import *


# Класс, определяющий параметры, относящиеся к КА
class Satellite:
    def __init__(self, coords, altitude):
        self.coords = coords
        self.altitude = altitude
        self.horde_level = 0.0

    # Обновление уровня (по оси Z) хорды, образованной линией пересечения поверхностей сферы и конуса
    def update_horde_level(self, alpha, radius):
        if ((self.altitude + radius) / radius) * np.sin(np.deg2rad(alpha)) >= 1:
            # Область обзора ограничена линией горизонта
            self.horde_level = radius**2 / (self.altitude + radius)
        else:
            # Область обзора ограничена линией пересечения поверхностей сферы и конуса
            self.horde_level = radius * np.cos(np.arcsin(((self.altitude + radius) / radius) *
                                                         np.sin(np.deg2rad(alpha))) - np.deg2rad(alpha))


# Получение высот КА из созвездия в метрах
def get_sat_altitudes(constellation: Constellation):
    sat_altitudes = [0] * constellation.total_sat_count
    sat_idx = 0

    for group in constellation.groups:
        for _ in range(group.get_total_sat_count()):
            sat_altitudes[sat_idx] = group.altitude * 1000
            sat_idx += 1

    return sat_altitudes


# Получение координат точек на поверхности сферы путём построения икосферы с заданным углом между точками
def get_sphere_points(angle_shift, radius):
    points_coords = []

    slices_count = int(np.ceil(180 / angle_shift)) + 1
    delta_angle = np.pi / (slices_count - 1)

    angle = -0.5 * np.pi

    for slice_idx in range(slices_count):
        circle_radius = np.abs(np.cos(angle))
        coord_z = np.sin(angle) * radius
        circle_slices_count = int(np.ceil(2.0 * np.pi * circle_radius / delta_angle))
        delta_lon_angle = 2.0 * np.pi / circle_slices_count

        if slice_idx % (slices_count - 1) == 0:
            circle_slices_count = 1
            delta_lon_angle = 0.0

        lon_angle = 0.0
        for _ in range(circle_slices_count):
            coord_x = circle_radius * np.cos(lon_angle) * radius
            coord_y = circle_radius * np.sin(lon_angle) * radius
            points_coords.append((coord_x, coord_y, coord_z))
            lon_angle += delta_lon_angle

        angle += delta_angle

    return points_coords


# Представление координат в новой системе координат с осью Z, направленной на KA
def get_transformed_coords(points_coords, sat_coords, radius):
    axis_z = np.array([0, 0, 1.0])
    normal_sat_coords = sat_coords / np.linalg.norm(sat_coords)

    rotation_vector = np.cross(axis_z, normal_sat_coords)
    theta = np.arccos(np.dot(axis_z, normal_sat_coords))

    rotation_matrix = np.array([
        [np.cos(theta) + rotation_vector[0]**2 * (1 - np.cos(theta)),
         rotation_vector[0] * rotation_vector[1] * (1 - np.cos(theta)) - rotation_vector[2] * np.sin(theta),
         rotation_vector[0] * rotation_vector[2] * (1 - np.cos(theta)) + rotation_vector[1] * np.sin(theta)],
        [rotation_vector[1] * rotation_vector[0] * (1 - np.cos(theta)) + rotation_vector[2] * np.sin(theta),
         np.cos(theta) + rotation_vector[1]**2 * (1 - np.cos(theta)),
         rotation_vector[1] * rotation_vector[2] * (1 - np.cos(theta)) - rotation_vector[0] * np.sin(theta)],
        [rotation_vector[2] * rotation_vector[0] * (1 - np.cos(theta)) - rotation_vector[1] * np.sin(theta),
         rotation_vector[2] * rotation_vector[1] * (1 - np.cos(theta)) + rotation_vector[0] * np.sin(theta),
         np.cos(theta) + rotation_vector[2]**2 * (1 - np.cos(theta))]
    ])

    transformed_points_coords = [[]] * len(points_coords)
    for point_idx, coords in enumerate(points_coords):
        transformed_coords = np.dot(rotation_matrix, np.array(coords))
        transformed_coords = transformed_coords / np.linalg.norm(transformed_coords) * radius
        transformed_points_coords[point_idx] = transformed_coords

    return transformed_points_coords


# Получение индексов точек в зоне покрытия в системе ECI
def get_covered_points(transformed_points_coords, horde_level):
    return [point_idx for point_idx, coords in enumerate(transformed_points_coords) if coords[2] >= horde_level]


# Получение процента покрытия сферы по точкам. Возвращаемое значение лежит в диапазоне [0, 1]
def get_coverage_percent(satellites, points_coords, alpha):
    covered_points = set()
    for sat in satellites:
        sat.update_horde_level(alpha, EARTH_RADIUS)
        transformed_points_coords = get_transformed_coords(points_coords, sat.coords, EARTH_RADIUS)

        for point_idx in get_covered_points(transformed_points_coords, sat.horde_level):
            covered_points.add(point_idx)

    coverage_percent = len(covered_points) / len(points_coords)
    return coverage_percent


def main():
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

    # Расчёт положений всех КА в начальный момент времени
    constellation.propagate_j2([0])

    # Установка минимального и максимального углов обзора, точности определения покрытия
    min_alpha = 0
    max_alpha = 60
    accuracy = 3

    alpha_list = np.arange(min_alpha, max_alpha + accuracy, accuracy)
    points_coords = get_sphere_points(accuracy, EARTH_RADIUS)

    sat_altitudes = get_sat_altitudes(constellation)
    sat_coords = constellation.state_eci[:, :, 0]
    sat_count = constellation.total_sat_count
    satellites = [Satellite(sat_coords[sat_idx], sat_altitudes[sat_idx]) for sat_idx in range(sat_count)]

    elements_dict = {}
    print("Выполняется двоичный поиск минимального значения угла обзора для глобального покрытия поверхности Земли...")

    # Выполнение двоичного поиска минимального значения угла обзора для глобального покрытия поверхности Земли
    def binary_search(arr, low, high):
        if high < low:
            print("Глобальное покрытие поверхности Земли невозможно при заданной конфигурации")
            return

        mid = (high + low) // 2

        cur_element = elements_dict.get(mid)
        if cur_element is None:
            cur_element = get_coverage_percent(satellites, points_coords, arr[mid])
            elements_dict[mid] = cur_element
            print(f"Угол обзора,(°): {arr[mid]}. Покрытие,(%): {round(cur_element * 100, 2)}")

        if cur_element != 1:
            binary_search(arr, mid + 1, high)
            return

        prev_element = elements_dict.get(mid - 1)
        if prev_element is None:
            prev_element = get_coverage_percent(satellites, points_coords, arr[mid - 1])
            elements_dict[mid - 1] = prev_element
            print(f"Угол обзора,(°): {arr[mid - 1]}. Покрытие,(%): {round(prev_element * 100, 2)}")

        if prev_element == 1:
            binary_search(arr, low, mid - 1)
            return

        print(f"Глобальное покрытие поверхности Земли возможно при угле обзора,(°): {arr[mid]}")

    binary_search(alpha_list, 0, len(alpha_list)-1)


if __name__ == "__main__":
    main()
