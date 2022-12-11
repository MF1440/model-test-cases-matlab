function initInitialState(constellation)
% Инициализация состояния космической группировки в начальный момент времени
    constellation.state.elements = zeros(constellation.totalSatCount, 6);

    shift = 1;
    for groupIdx = 1 : numel(constellation.groups)
        thisGroup = constellation.groups{groupIdx};
        ending = shift + thisGroup.totalSatCount - 1;
        constellation.state.elements(shift:ending,:) = constellation.initGroupElements(thisGroup);
        shift = ending + 1;
    end
end