function [d,CPU_time,J_OPTM] = solve_calibration2(d)
global sim
% Initial Accumulation 
n0 = ones(d.nb_R,d.nb_R).*1*10^(-190);

% Demand matrix 
d.lbw =[];
d.ubw = [];

% initial=  load('TLAdjustmenmt\ULim\initial.mat');
 
alp = ones(d.nb_R*d.nb_R*sum(d.tds_update<=d.t_int),1);

d.TL_llim = 0.70;
d.TL_ulim = 1.30;
d.TL_int = alp;

d.lbw =  alp(:)*d.TL_llim;
d.ubw =  alp(:)*d.TL_ulim;

J_OPTM = manual_OptV2_1(alp(:),n0(:),d)

tic
% Solve NLP
    d.NLP_solution = d.f_NLP('x0', ... % initial guess
        vertcat(alp(:)), ...
        'p', ... % parameters
        vertcat(n0(:)),...
        'lbx', ... % lower bounds for box constraints
        vertcat(d.lbw(:)), ...
        'ubx', ... % upper bounds for box constraints
        vertcat(d.ubw(:)), ...
        'lbg', ... % lower bounds for general constraints
        vertcat(d.lbg(:)), ...
        'ubg', ... % upper bounds for general constraints
        vertcat(d.ubg(:)));
  
 CPU_time = toc    
 
    d.ALPfull = full(d.NLP_solution.x);
    for k = 1:d.t_int  
        if ismember(k,d.tds_update)==1
            tx = find(d.tds_update==k);
            d.ALPIH_t{k,sim} = reshape(d.ALPfull(d.nb_R*d.nb_R*(tx-1)+1:d.nb_R*d.nb_R*tx,1),d.nb_R,d.nb_R);
        else
            d.ALPIH_t{k,sim}= d.ALPIH_t{k-1,sim};
        end
    end
        
  end