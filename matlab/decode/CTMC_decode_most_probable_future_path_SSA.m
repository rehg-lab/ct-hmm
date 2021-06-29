function [best_state_seq_SSA, best_prob_SSA] = CTMC_decode_most_probable_future_path_SSA(start_s, T)

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

SSAProb.Q_mat = Q_mat;

SSAProb.L = lambda_list; % lambda list
SSAProb.T = Vij_mat;  % transition prob: vij

SSAProb.Starts = k;
SSAProb.Time = T;
SSAProb.MaxDom = 0;

SSAProb.HasSpecificEndState = 1;

Pt = expm(Q_mat * T);

SSAProb.Ends = zeros(1, num_state);
num_end = 0;
for i = 1:num_state
    if (Pt(k, i) > 0.0)
        num_end = num_end + 1;
        SSAProb.Ends(num_end) = i;    
    end        
end
SSAProb.Ends = SSAProb.Ends(1:num_end);


num_end

% Solve the optimization problem
%tic;

SSARes = StateSequenceAnalyze(SSAProb);
SSATime = toc;    

StartStatesOrWeights = k;
EndStatesOrWeights = ones(1, num_state);
TimesToDo = T;
MMostProbable = 1;
[MaxSeqsByTime,SeqList] = ExtractMaxSeqs(SSARes,TimesToDo,StartStatesOrWeights,EndStatesOrWeights,MMostProbable);    

best_end_state = SeqList(1, 2); 
best_seq_idx = SeqList(1, 3); % first row, the 3rd component is the best sequence index

best_state_seq_SSA = SSARes.Seqs{k,best_end_state}{best_seq_idx}.seq';
best_prob_SSA = SSARes.Seqs{k,best_end_state}{best_seq_idx}.p(end);

%tEnd = toc;
%str = sprintf('SSA: %d min, %f sec\n', floor(tEnd/60),rem(tEnd,60))

