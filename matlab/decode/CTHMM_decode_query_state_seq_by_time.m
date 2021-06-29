%% query the state sequence of this segment (begin_time ~ end_time)
function [seg_state_seq, seg_dwelltime_seq] = CTHMM_decode_query_state_seq_by_time(state_seq, dwelltime_seq, seg_begin_time, seg_end_time)

num_state = length(state_seq);
seg_state_seq = zeros(num_state, 1);
seg_dwelltime_seq = zeros(num_state, 1);

state_begin_time = 0;
seg_state_count = 0;

for s = 1:num_state

    state_end_time = state_begin_time + dwelltime_seq(s);
    
    %% if the state's dwelling time is overlapped in the time segment
    if (state_end_time < seg_begin_time)
         %% continue on checking next state
    elseif (state_begin_time > seg_end_time)
        break;        
    else %% record the overlap        
        if (state_begin_time <= seg_begin_time)
            %% this is the first state for this time segment
            seg_state_seq(1) = state_seq(s);
            
            if (state_end_time <= seg_end_time)
                seg_dwelltime_seq(1) = dwelltime_seq(s) - (seg_begin_time - state_begin_time);
            else
                seg_dwelltime_seq(1) = seg_end_time - seg_begin_time;
            end
            
            seg_state_count = 1;
        else
            %% not the first state
            seg_state_count = seg_state_count + 1;            
            seg_state_seq(seg_state_count) = state_seq(s);
            
            if (state_end_time <= seg_end_time)
                %% entire dwell time is within the segment
                seg_dwelltime_seq(seg_state_count) = dwelltime_seq(s);
            else
                %% part of the state dwell time is within the segment
                seg_dwelltime_seq(seg_state_count) = seg_end_time - state_begin_time;
            end            
        end                
    end
    
    %% update the begin time for the next state
    state_begin_time = state_end_time;
       
end

seg_state_seq = seg_state_seq(1:seg_state_count);
seg_dwelltime_seq = seg_dwelltime_seq(1:seg_state_count);

