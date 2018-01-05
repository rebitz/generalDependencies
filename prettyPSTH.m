% prettyPSTH

figure(gcf);
prop = get(gca);

set(h,'LineWidth',2)

muaStr1 = 'sig\d\d\di';
muaStr2 = 'sig\d\d\dall';

if ~isempty(regexp(sigOI,muaStr1)) || ~isempty(regexp(sigOI,muaStr2))
    ylabel('firing rate (MUA: spikes/s)')
else
    ylabel('firing rate (spikes/s)')
end

if strcmp(evtOI,'tmpEVT')
    xlabel('time from saccade onset (ms)')
elseif strcmp(evtOI,'EVT07')
    xlabel('time from target onset (ms)')
else
    xlabel('time from event (ms)')
end

set(prop.YLabel,'FontSize',16)
set(prop.XLabel,'FontSize',16)
set(prop.Title,'FontSize',16)
set(gca,'FontSize',16)

temp = ylim;
h = line([0 0],[temp(1) temp(2)]);
set(h,'LineWidth',2,'LineStyle','--','Color',[.5 .5 .5])