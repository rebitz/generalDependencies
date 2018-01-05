
idx = [trials.correct] == 1;


fixOn = floor(nanmean([trials(idx).fixOn]));
fixOff = floor(nanmean([trials(idx).fixAcq]+900));

targOn = floor(nanmean([trials(idx).targOn]));
targOff = floor(nanmean([trials(idx).targOff]));

blank = zeros(1,2000);

fixTrace = blank;
fixTrace(fixOn:fixOff) = 1;
targTrace = blank;
targTrace(targOn:targOff) = 1;

figure(); hold on;
plot(fixTrace+10,'-k');
plot(targTrace+8,'-k');

%%
idx = and([trials.theta] == 180,[trials.ecc] == 8);
idx = find(idx,20,'first');
idx = idx(5:10)

starts = [trials(idx).targAcq];

pre = 1500;
post = 500;

figure(); hold on;
these = idx;
for i = 1:length(these)
    tmp = [trials(these(i)).eye];
    eyeId = and(tmp(:,1)>starts(i)-pre,tmp(:,1)<starts(i)+post);
    tstamps = tmp(eyeId,1);
    tmp = tmp(eyeId,2);
%     tmp = (tmp-min(tmp))./(max(tmp)-min(tmp));
    plot(tstamps,gsmooth(tmp,10));
end

% ylim([0 12])