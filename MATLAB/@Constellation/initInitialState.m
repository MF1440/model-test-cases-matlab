function initInitialState(this)
    this.state.elements = zeros(this.totalSatCount, 6);

    shift = 1;
    for groupIdx = 1 : numel(this.groups)
        thisGroup = this.groups{groupIdx};
        ending = shift + thisGroup.totalSatCount - 1;
        this.state.elements(shift:ending,:) = this.initGroupElements(thisGroup);
        shift = ending + 1;
    end
end