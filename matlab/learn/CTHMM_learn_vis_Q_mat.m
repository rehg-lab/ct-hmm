function CTHMM_learn_vis_Q_mat(top_out_folder)

global data_setting;
global is_draw_learn_Q_mat;

if (is_draw_learn_Q_mat == 1)
    if (data_setting.dim == 2)    
        
%         is_vis_nij = 1; % vis nij
%         vis_threshold = 0.1;
%         is_draw_text = 1.0;        
%         CTHMM_vis_2D_Q_mat(top_out_folder, is_vis_nij, vis_threshold, is_draw_text);    
%         
%         is_vis_nij = 0;  % vis qij
%         vis_threshold = 0.01;
%         is_draw_text = 1.0;
%         CTHMM_vis_2D_Q_mat(top_out_folder, is_vis_nij, vis_threshold, is_draw_text);   
      
        vis_mat_type = 1; % nij
        %vis_threshold = 0.000001;
        vis_threshold = 0.0;
        is_draw_text = 1;
        state_is_vis_tau_i = 1; % draw total time at each state as state size
        CTHMM_vis_2D_Q_mat_new(top_out_folder, vis_mat_type, vis_threshold, is_draw_text, state_is_vis_tau_i);


        vis_mat_type = 1; % nij
        %vis_threshold = 0.000001;
        vis_threshold = 0.0;
        is_draw_text = 1;
        state_is_vis_tau_i = 0; % draw n_i = \sum n_ij, as state size
        CTHMM_vis_2D_Q_mat_new(top_out_folder, vis_mat_type, vis_threshold, is_draw_text, state_is_vis_tau_i);
        
        vis_mat_type = 3; % vij
        %vis_threshold = 0.000001;
        vis_threshold = 0.0;
        is_draw_text = 1;
        state_is_vis_tau_i = 0; % draw n_i = \sum n_ij, as state size
        CTHMM_vis_2D_Q_mat_new(top_out_folder, vis_mat_type, vis_threshold, is_draw_text, state_is_vis_tau_i);
        
    elseif (data_setting.dim == 3)    
        is_vis_nij = 1; 
        vis_threshold = 0.1;
        CTHMM_vis_3D_Q_mat(top_out_folder, is_vis_nij, vis_threshold);
        
        is_vis_nij = 0; % vis qij
        vis_threshold = 0.01;
        CTHMM_vis_3D_Q_mat(top_out_folder, is_vis_nij, vis_threshold);        
    end
end


