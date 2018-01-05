figure(gcf);
prop = get(gca);

try,
    set(h,'LineWidth',2)
end

set(prop.YLabel,'FontSize',20)
set(prop.XLabel,'FontSize',20)
set(prop.Title,'FontSize',20)
set(gca,'FontSize',20)

if ~exist('noXTickChange','var') || noXTickChange == 0;
    set(gca,'XTick',[1:length(unique(xpos))])

    tmp={};
    for i = 1:length(unique(xpos))
        tmp = [tmp num2str(xpos(i))];
    end
    set(gca,'XTickLabel',tmp);
end