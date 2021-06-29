%% this tool helps you find the cases when specifying the decoded states
%% please change the following variables:
%% file_decoded_eye_seq_list : change it to the path for your "decoded_eye_seq_list" variable
%% model_top_dir: change it to the top folder of your trained model
%% query_state_list: the target states you want to query the patients with. Please always specify the state range with value from low to high. 
%%    You can specify several states. The found cases shall contain the specified states in their decoded results.

clc;

global eye_seq_list;
global state_list;
addpath('../../learn');

%% load the decoded results
file_decoded_eye_seq_list = 'output_decoding_minvisit5/decoded_eye_seq_list';
temp = load(file_decoded_eye_seq_list); % load the variable "decoded_eye_seq_list"
eye_seq_list = temp.eye_seq_list;

%% load the model variables
model_top_dir = 'output_decoding_minvisit5';
str = sprintf('%s/CV_1/num_iter.txt', model_top_dir);  % read the number of iteration from num_iter.txt file	   
fp = fopen(str, 'rt');
num_iter = fscanf(fp, '%d');
CTHMM_model_dir = sprintf('%s/CV_1/Iter_%d', model_top_dir, num_iter); 
CTHMM_learn_load_para(CTHMM_model_dir);
    

%% here I specifiy two states to query. The found case shall contain these two states in their decoded path
query_state_list = {[96 98; 75 80]; ...  
                    [94 96; 75 80] };

%% this is an example for specifying just one state
%query_state_list = {[96 98; 85 90]};          
                
num_query_state = length(query_state_list);

num_eye = size(eye_seq_list, 1);

for e = 1:num_eye
    
    decoded_state_seq = eye_seq_list{e}.decoded_conti_state_seq;    
    num_decoded_state = length(decoded_state_seq);
            
    for q = 1:num_query_state % for each querying state        
        has_query_state = 0;
        for p = 1:num_decoded_state  % for each state in the decoded path            
            i = decoded_state_seq(p);
            if ((state_list{i}.range(1,1) == query_state_list{q}(1,1)) && (state_list{i}.range(1,2) == query_state_list{q}(1,2)) && ...
                (state_list{i}.range(2,1) == query_state_list{q}(2,1)) && (state_list{i}.range(2,2) == query_state_list{q}(2,2)))   
                has_query_state = 1;
                break;
            end            
        end        
        if (has_query_state == 0)  % this eye doesn't have this querying state
            break;
        end        
    end
    
    if (has_query_state == 1)
        outstr = sprintf('eye index %d, ID: %d, eye: %s has the query states', e, eye_seq_list{e}.ID, eye_seq_list{e}.eye{1});
        disp(outstr);
    end
    
end
