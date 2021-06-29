function [best_state_seq_SSA, best_prob_SSA] = CTMC_decode_most_probable_state_seq_SSA(start_s, end_s, T)

global Q_mat;

%% Set up qi list and Vij_mat (transition prob matrix)
num_state = size(Q_mat, 1);
lambda_list = zeros(num_state, 1);
Vij_mat = zeros(num_state, num_state);
for i = 1:num_state % row
    qi = -Q_mat(i, i);
    lambda_list(i) = qi;
    row = Q_mat(i, :);
    Vij_mat(i, :) = row ./ qi;    
end

%% set up all variables
k = start_s;
l = end_s;

SSAProb.L = lambda_list; % lambda list
SSAProb.T = Vij_mat;  % transition prob: vij

SSAProb.Starts = k;
SSAProb.Time = T;
SSAProb.MaxDom = 0;

SSAProb.HasSpecificEndState = 1;
SSAProb.Ends = l;

% Solve the optimization problem 
%tic;

SSAProb.Q_mat = Q_mat;

SSARes = StateSequenceAnalyze(SSAProb);
SSATime = toc;    

StartStatesOrWeights = k;
EndStatesOrWeights = l;
TimesToDo = T;
MMostProbable = 1;
disp("this is time to do")
disp(TimesToDo)
[MaxSeqsByTime,SeqList] = ExtractMaxSeqs(SSARes,TimesToDo,StartStatesOrWeights,EndStatesOrWeights,MMostProbable);    
best_seq_idx = SeqList(1, 3); % first row, the 3rd component is the best sequence index

best_state_seq_SSA = SSARes.Seqs{k,l}{best_seq_idx}.seq';
best_prob_SSA = SSARes.Seqs{k,l}{best_seq_idx}.p(end);

%tEnd = toc;
%str = sprintf('SSA: %d min, %f sec\n', floor(tEnd/60),rem(tEnd,60))

