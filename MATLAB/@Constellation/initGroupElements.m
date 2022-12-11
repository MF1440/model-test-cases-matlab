function elements = initGroupElements(group)
% Инициализация параметров элементов группы космических аппаратов
    metersPerKm = 1000;

    raans = linspace(group.startRaan, group.startRaan + group.maxRaan, group.planeCount + 1);
    raans = mod(raans(1:end-1), 2 * pi);
    
    elements = zeros(group.totalSatCount, 6);
    elementIdx = 1;
    raanIdx = 0;
    for thisRaan = raans
        for satIdx = 1:group.satsPerPlane
            sma = Constants.AstroConstants.earthRadius + group.altitudeKilometers * metersPerKm;
            aol = 2 * pi / group.satsPerPlane * (satIdx - 1) + 2 * pi / group.totalSatCount * group.f * raanIdx;
    
            elements(elementIdx, :) = [sma, 0, 0, thisRaan, group.inclination, aol];
            elementIdx = elementIdx + 1;
        end
        raanIdx = raanIdx + 1;
    end
end