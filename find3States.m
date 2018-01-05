function [states,stateProbs,LL,nParams] = find3States(choices)

% seed the transmission and emissions
T = [(1/4) (1/4) (1/4) (1/4);...
    (1/4) (3/4) (0/4) (0/4);...
    (1/4) (0/4) (3/4) (0/4);...
    (1/4) (0/4) (0/4) (3/4)];

E = [(1/3) (1/3) (1/3);...
        1     0     0;...
        0     1     0;...
        0     0     1];
    
% train the HMM
[Tall,Eall] = hmmtrain(choices, T, E);

% get the most likely sequence of states
states = hmmviterbi(choices, Tall, Eall);
% tmp = num2cell(likelyStates);
% [states] = deal(tmp{:});

% now append the probabilities
[stateProbs,LL] = hmmdecode(choices, Tall, Eall);
% tmp = num2cell(stateProbs);
% [stateProbs] = deal(tmp{1,:});
% only the p(ORE)!!

nParams = sum(sum(E~=0,2)-1) + sum(sum(T~=0,2)-1);