% plots saccades evoked from microstimulation
% first, pull all the files in a directory in together

files = dir;
stimfiles = ~cellfun(@isempty,strfind({files.name},'FIXSTIM'));
stimfiles = and(stimfiles,cellfun(@isempty,strfind({files.name},'site')));
stimfiles = {files(stimfiles).name};

out = [];

for i = 1:length(stimfiles)
    load(stimfiles{i})
    out = [out trials];
end

%%

degToPlot = 10;

selex = and([out.voltage] < 150,[out.duration] > 1);
% selex = and(selex,[out.saccade]==1)
selex = find(selex)

% selex = 4

preWindow = 0.1; %preWindow = preWindow*1000;
postWindow = 0.5; %postWindow = postWindow*1000;

saccadeStart = 0.2*1000;
saccadeEnd = 0.3*1000;

% setup the figure
h = figure(99); clf; hold on;
set(gcf, 'Position', [100 109 509 841]);
subplot(6,1,4:6); hold on;
% h(1) = line([-deg2px(20,env) deg2px(20,env)],[0 0]);
% h(2) = line([0 0],[-deg2px(20,env) deg2px(20,env)]);
h(1) = line([-degToPlot degToPlot],[0 0]);
h(2) = line([0 0],[-degToPlot degToPlot]);

set(h,'Color',[.5 .5 .5],'LineStyle','--');
    
for i = 1:length(selex);

    xPos = px2deg(out(selex(i)).x,env);
    yPos = px2deg(out(selex(i)).y*-1,env);
    tstamps = out(selex(i)).tstamps;
    pulse = out(selex(i)).pulse;
    
    try,
        if out(selex(i)).beforeDuringAfter == 1;
            c = 'k';
        else
            c = [.5 .5 .5];
        end
    catch
        c = 'k';
    end
    
    tFromStim = tstamps - min(tstamps);
    tFromStim = tstamps - out(selex(i)).stimOnTime;
    
    tFromStim = 1:length(tstamps);
    
    % then plotting
    figure(99);
    subplot(6,1,1); hold on; % pulse
    plot(tFromStim,pulse,'Color',c);
    ylabel('stim pulse (cartoon)');
    
    subplot(6,1,2); hold on; % x pos
    plot(tFromStim,xPos-nanmean(xPos(1:preWindow*1000)),'Color',c);
    ylabel('x position (px)');
    
    subplot(6,1,3); hold on; % y pos
    plot(tFromStim,yPos-nanmean(yPos(1:preWindow*1000)),'Color',c);
    ylabel('y position (px)');
    
    subplot(6,1,4:6); hold on;
    plot(xPos(saccadeStart:saccadeEnd)-xPos(saccadeStart),...
        yPos(saccadeStart:saccadeEnd)-yPos(saccadeStart),'Color',c,...
        'LineWidth',1);
    xlabel('xpos (px)'); ylabel('ypos (px)');
    xlim([-degToPlot degToPlot]);
    ylim([-degToPlot degToPlot]);
%     xlim([-env.screenWidth/6 env.screenWidth/6]);
%     ylim([-env.screenHeight/6 env.screenHeight/6]);

    
    % append plot info to trial struct
%    trials(tNum).preEyeTime = preEyeTime;
%    trials(tNum).postEyeTime = postEyeTime;

end