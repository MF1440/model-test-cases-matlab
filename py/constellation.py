import json
import numpy as np

from typing import NamedTuple

from consts import *


class Walker(NamedTuple):
    inclination:    float   # Наклонение орбиты
    sats_per_plane: int     # Число КА в каждой орбитальной плоскости группы
    plane_count:    int     # Число орбитальных плоскостей в группе
    phase_shift:    int     # Фазовый сдвиг по аргументу широты между КА в соседних плоскостях
    altitude:       float   # Высота орбиты
    max_raan:       float   # Максимум прямого восхождения восходящего узла (при распределении орбитальных плоскостей)
    start_raan:     float   # Прямое восхождение восходящего узла для первой плоскости


class WalkerGroup(Walker):

    def get_total_sat_count(self):
        return self.sats_per_plane * self.plane_count

    def calculate_initial_elements(self):
        start_raan = np.deg2rad(self.start_raan)
        max_raan = np.deg2rad(self.max_raan)
        inclination = np.deg2rad(self.inclination)
        altitude = self.altitude * 1000
        sat_count = self.get_total_sat_count()

        sma = EARTH_RADIUS + altitude
        raans = np.linspace(start_raan, start_raan + max_raan, self.plane_count + 1)
        raans = raans[:-1] % (2 * np.pi)

        elements = np.zeros((sat_count, 6))
        element_idx = 0

        for raan_idx, raan in enumerate(raans):
            for sat_idx in range(self.sats_per_plane):
                aol = 2 * np.pi * (sat_idx / self.sats_per_plane + self.phase_shift * raan_idx / sat_count)

                elements[element_idx, :] = [sma, 0, 0, raan, inclination, aol]
                element_idx += 1

        return elements


class Constellation:

    def __init__(self):
        self.total_sat_count = 0
        self.groups = []
        self.elements = []
        self.state_eci = []

    def load_from_config(self, file_path, name_code):
        try:
            with open(file_path) as fp:
                json_data = json.load(fp)
        except FileNotFoundError:
            raise Exception("Файл не найден")

        for constellation_data in json_data:
            current_name = constellation_data.get("name")
            if current_name is None:
                continue

            if current_name.lower() == name_code.lower():
                current_walkers = constellation_data.get("Walkers")
                if not current_walkers:
                    raise Exception(f"Данные о спутниках группировки {name_code} не найдены")

                try:
                    for walker_data in current_walkers:
                        walker_group = WalkerGroup(*walker_data)
                        self.groups.append(walker_group)
                        self.total_sat_count += walker_group.get_total_sat_count()
                except TypeError:
                    raise Exception(f"Формат данных о спутниках группировки {name_code} некорректен")

                print(f"Загружена группировка {name_code}")
                return

        raise Exception("Группировка не найдена в файле")

    def initialize_state(self):
        if self.total_sat_count <= 0:
            return

        self.elements = np.zeros((self.total_sat_count, 6))
        shift = 0

        for single_group in self.groups:
            ending = shift + single_group.get_total_sat_count()
            self.elements[shift:ending, :] = single_group.calculate_initial_elements()
            shift = ending

    def propagate_j2(self, epoch_list):
        if len(self.elements) == 0:
            return

        self.state_eci = np.zeros((self.total_sat_count, 3, len(epoch_list)))

        inclination = self.elements[:, 4]
        sma = self.elements[:, 0]
        raan0 = self.elements[:, 3]
        aol0 = self.elements[:, 5]

        raan_precession_rate = -1.5 * (EARTH_J2 * np.sqrt(EARTH_GM) * EARTH_RADIUS ** 2) / (sma ** (7 / 2)) \
            * np.cos(inclination)

        draconic_omega = np.sqrt(EARTH_GM / sma ** 3) * (1 - 1.5 * EARTH_J2 * (EARTH_RADIUS / sma) ** 2) \
            * (1 - 4 * np.cos(inclination)**2)

        for epoch in epoch_list:
            aol = aol0 + epoch * draconic_omega
            raan_omega = raan0 + epoch * raan_precession_rate

            epoch_state = sma * [
                (np.cos(aol) * np.cos(raan_omega) - np.sin(aol) * np.cos(inclination) * np.sin(raan_omega)),
                (np.cos(aol) * np.sin(raan_omega) + np.sin(aol) * np.cos(inclination) * np.cos(raan_omega)),
                (np.sin(aol) * np.sin(inclination))
            ]

            self.state_eci[:, :, epoch_list.index(epoch)] = np.array(epoch_state).T
