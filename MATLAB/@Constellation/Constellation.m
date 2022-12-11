classdef Constellation < handle
% Класс, описывающий космическую группировку и её состояние в заданные моменты времени
    properties
        totalSatCount = 0;
        groups = {}; % Массив структур для хранения основных параметров группы спутников
        state; % Структура, определяющая состояние космических аппаратов в заданные моменты времени
    end

    methods (Static)
        elements = initGroupElements(group)
    end

    methods

        createConstellationFromJson(constellation, fileName, constellationName)

        initInitialState(constellation)

        function this = Constellation(varargin)
            % Конструктор класса. 
            % varargin{1} - имя конфигурационного .json файла космической группировки
            % varargin{2} - название космической группировки внутри конфигурациионного файла
            if nargin==2 && ischar(varargin{1})
                [~, ~, extension] = fileparts(varargin{1});
                if strcmp(extension, '.json')
                    createConstellationFromJson(this, varargin{1}, varargin{2});
                else
                    error('Неизвеcтное расширение конифгурационного файла')
                end
            elseif isempty(varargin)
                return
            else
                error('Создание экземпляра класса с использованием данного набора аргументов не поддерживается!')
            end
        end

    end
end
