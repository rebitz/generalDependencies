%% GPFA preprocessing: preprocessGPFA

% struct called dat w/ fields
%    trialId (trial number, I guess?)
%    spikes (a neurons X timestamps 0/1 matrix)

targOnEvt = 'EVT05';
% append sacade timestamps
temp = num2cell([data.(targOnEvt)] + (([data.saccadeStart] - [data.targOn])));% ./ 1000));
[data.tmpEVT] = deal(temp{:});

evtOI = targOnEvt;
pre = 300;
post = 100;

tmp = fieldnames(data);
tmp = tmp(and(~cellfun(@isempty,strfind(tmp,'sig')),...
    cellfun(@isempty,strfind(tmp,'all'))));

% initialize output structure
clear dat; started = 0;
dat = struct('trialId',num2cell([1:length(data)]));

for i = 1:length(tmp)
    [~,~,~,~,~,~,A] = getSpikes2(data,evtOI,tmp{i},pre,post);
    mx = full(A)./1000;
    
    % remove low FR neurons
    if sum(full(nanmean(A,2)) > 5)==length(data)
        for j = 1:length(dat)
            if started % if we've already done the first cell
                fu = dat(j).spikes;
                fu = [fu; mx(j,:)];
            else
                started = 1;
                fu = mx(j,:);
            end
            dat(j).spikes = fu;
        end
    end
end

%%

% Results will be saved in mat_results/runXXX/, where XXX is runIdx.
% Use a new runIdx for each dataset.
runIdx = 3;

% Select method to extract neural trajectories:
% 'gpfa' -- Gaussian-process factor analysis
% 'fa'   -- Smooth and factor analysis
% 'ppca' -- Smooth and probabilistic principal components analysis
% 'pca'  -- Smooth and principal components analysis
method = 'gpfa';

% Select number of latent dimensions
xDim = 4;
% NOTE: The optimal dimensionality should be found using 
%       cross-validation (Section 2) below.

% If using a two-stage method ('fa', 'ppca', or 'pca'), select
% standard deviation (in msec) of Gaussian smoothing kernel.
kernSD = 100;
% NOTE: The optimal kernel width should be found using 
%       cross-validation (Section 2) below.

% Extract neural trajectories
result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
                    'kernSDList', kernSD);
% NOTE: This function does most of the heavy lifting.

% Orthonormalize neural trajectories
[estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
% NOTE: The importance of orthnormalization is described on 
%       pp.621-622 of Yu et al., J Neurophysiol, 2009.

% Plot neural trajectories in 3D space
plot3D(seqTrain, 'xorth', 'dimsToPlot', 1:3);
% NOTES:
% - This figure shows the time-evolution of neural population
%   activity on a single-trial basis.  Each trajectory is extracted from
%   the activity of all units on a single trial.
% - This particular example is based on multi-electrode recordings
%   in premotor and motor cortices within a 400 ms period starting 300 ms 
%   before movement onset.  The extracted trajectories appear to
%   follow the same general path, but there are clear trial-to-trial
%   differences that can be related to the physical arm movement. 
% - Analogous to Figure 8 in Yu et al., J Neurophysiol, 2009.
% WARNING:
% - If the optimal dimensionality (as assessed by cross-validation in 
%   Section 2) is greater than 3, then this plot may mask important 
%   features of the neural trajectories in the dimensions not plotted.  
%   This motivates looking at the next plot, which shows all latent 
%   dimensions.

%%
fu = find(and([data.choice]==1,[data.states]~=1));
plotEachDimVsTime(seqTrain, 'xorth', result.binWidth,...
    'redTrials',fu,'nPlotMax',50);