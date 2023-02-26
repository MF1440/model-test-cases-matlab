function []  =  plotTotalTraffic(alfaList, satListEcef, coordsEcef, values, constellation)
    %% Строит график зависимости удовлетворённого спроса от угла α.
    % alfaList - массив углов α
    % satListEcef - массив координат спутников в СК, связанной с вращающейся Землей
    % coordsEcef - массив координат точек в осях связанных с вращающейся Землёй
    % values - массив запрашиваемого абонентами трафика
    % constellation - объект типа Constellation
    
    %% Выделим абонентов с ненулевым запросом трафика

    activeUserIdx = find(values);
    activeUserEcef = coordsEcef(activeUserIdx,:);
    activeUserValues = values(activeUserIdx);
    
    %% Найдем для каждого абонента ближайший спутник
    
    distUserSat = zeros(size(activeUserEcef,1),size(satListEcef,1));
    
    for sat = 1:size(satListEcef,1)   
       satPos = repmat(satListEcef(sat,:),size(activeUserEcef,1),1);
       distUserSat(:,sat) = sqrt( sum( (activeUserEcef - satPos).^2, 2 ));
    end
    
    % расстоние до ближайшего спутника и индекс этого спутника для каждого абонента    
    [minDistUserSat,satInd] = min(distUserSat,[],2);

    % удаленность спутников от центра Земли
    satSma = constellation.state.elements(:,1);
    
    % Удаленность абонентов от центра Земли
    userDist = sqrt( sum( (activeUserEcef ).^2, 2 )); % Полагаю, тут можно было использовать радиус Земли из constellation, но вдруг кто-то из абонентов в самолете

    totalValue=zeros(1, length(alfaList));
    
    for alfaIdx = 1:1:length(alfaList)
        for userIdx=1:size(activeUserEcef,1)
            satIdx = satInd(userIdx);
            a = satSma(satIdx);
            b = minDistUserSat(userIdx);
            c = userDist(userIdx);

            angleSatUser = acosd((a^2 + b^2 - c^2)/(2 * a * b));

            if angleSatUser < alfaList(alfaIdx)

                totalValue(alfaIdx) = totalValue(alfaIdx) + activeUserValues(userIdx);

            end % Конец проверки величины угла

        end     % Конец цикла по активным абонентам
        
    end         % Конец цикла по значениям угла
    
    %% Строим график
    
    fig = figure;
    fig.Position = [500,300,900,400];
    plot(alfaList, totalValue, 'linewidth', 2); grid on
    xlabel('\alpha\circ', 'fontsize', 14)
    ylabel('TotalTraffic', 'fontsize', 14)
    title('Зависимость удовлетворённого спроса от угла \alpha\circ', 'fontsize', 14, 'FontWeight', 'normal')

end

