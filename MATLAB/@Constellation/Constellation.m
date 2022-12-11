classdef Constellation < handle

    properties
        totalSatCount = 0;
        groups = {};
        state;
    end

    methods

        createConstellationFromJson(constellation, fileName, constellationName)

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
            end
        end

        function getInitialState(this)
            this.state.elements = zeros(this.totalSatCount, 6);
            shift = 1;

            for group = this.groups
                for i = 1:length(group{1})
                    ending = shift + group{1}(i).totalSatCount - 1;
                    this.state.elements(shift:ending,:) = this.getInitialElements(group{1});
                    shift = ending + 1;
                end
            end
        end

        function elements = getInitialElements(this, group)
            
            earthRadius = Constants.AstroConstants.earthRadius;
            metersPerKm = 1000;

            raans = linspace(group.startRaan, group.startRaan + group.maxRaan, group.planeCount + 1);
            raans = mod(raans(1:end-1), 2 * pi);

            elements = zeros(group.totalSatCount, 6);
            idx = 1;
            raanIDX = 0;
            for raan = raans
                for i = 0:group.satsPerPlane-1
                    sma = earthRadius + group.altitudeKilometers * metersPerKm;
                    aol = 2 * pi / group.satsPerPlane * i + 2 * pi / group.totalSatCount * group.f * raanIDX;

                    elements(idx, :) = [sma, 0, 0, raan, group.inclination, aol];
                    idx = idx + 1;
                end
                raanIDX = raanIDX + 1;
            end
        end        
    end
end
