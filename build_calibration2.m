function d = build_calibration2(hrs, tmin, d)

addpath('casadi/casadi-windows-matlabR2016a-v3.5.5')
import casadi.*
    d = build_scenario(hrs, tmin, d);

    d = define_dynamics(d);
    
    d = create_prediction_integrator(d);
     
    d = formulate_NLP(d);
    
end

function d = build_scenario(hrs, tmin, d)
addpath('casadi/casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

%todetermine the one hour or two hours 
 d.seed =hrs; 
 
 % final time    
 d.t_f =3600 *d.seed;
 
  %time stepin minutes
 d.t_prd = tmin;
 
 %accumulation aggregation interval 
 d.dt =d.t_prd*60; 
 
 % Number of time intervals 
 d.t_int =d.t_f/d.dt;  
 
 % Number of regions 
 d.nb_R =5;
 
 % Neighborhoods matrix
 d.nbrd = [0 1 1 1 1;
           1 0 1 1 0;
           1 1 0 0 1;
           1 1 0 0 1;
           1 0 1 1 0];
 
 % time step [s]
 d.ts = d.t_prd*60; 
 
 % integration step
 d.is = d.t_prd*60;
 
 %upplerlimit for TL 
 d.TLulim = 5;
 
 % solver type
 d.NLP_solver = 'ipopt';
 %d.NLP_solver = 'sqpmethod';

 % Maximum number of NLP iterations
 d.max_iter_NLP = 150;
 d.max_iter_sqpmethod=50;

 % Choose prediction integrator
 d.pred_int = 'RK';
 %     d.pred_int = 'Euler';

end

function d = define_dynamics(d)
    addpath('casadi/casadi-windows-matlabR2016a-v3.5.5')
    import casadi.*
    
     % Accumulation states n_I(J,H)
    n = SX.sym('n',d.nb_R,d.nb_R);

    % Derivatives of accumulation states
    n_dot = SX.zeros(d.nb_R,d.nb_R);

    %Mtansfer staATE Variable
    M_1 = SX.zeros(d.nb_R,d.nb_R);
    M_2 = SX.zeros(d.nb_R,d.nb_R);
    M_3 = SX.zeros(d.nb_R,d.nb_R);
    M_4 = SX.zeros(d.nb_R,d.nb_R);
    M_5 = SX.zeros(d.nb_R,d.nb_R);
    M_II= SX.zeros(d.nb_R,d.nb_R);
    
    % Flow demands
    q = SX.sym('q',d.nb_R,d.nb_R);
   
    % Route guidance inputs
    Th_1 = SX.sym('th_1',d.nb_R,d.nb_R);
    Th_2 = SX.sym('th_2',d.nb_R,d.nb_R);
    Th_3 = SX.sym('th_3',d.nb_R,d.nb_R);
    Th_4 = SX.sym('th_4',d.nb_R,d.nb_R);
    Th_5 = SX.sym('th_5',d.nb_R,d.nb_R);

    % Average trip length coefficients inputs(only the neighborly related) 
    TL = SX.sym('TL',d.nb_R,d.nb_R);
    
%     % Average trip length variations
     ALP  = SX.sym('ALP',d.nb_R,d.nb_R);


    %Find the Total Region Accumulation
    n_R1 = n(1,1)+n(1,2)+n(1,3)+n(1,4)+n(1,5);
    n_R2 = n(2,1)+n(2,2)+n(2,3)+n(2,4)+n(2,5);
    n_R3 = n(3,1)+n(3,2)+n(3,3)+n(3,4)+n(3,5);
    n_R4 = n(4,1)+n(4,2)+n(4,3)+n(4,4)+n(4,5);
    n_R5 = n(5,1)+n(5,2)+n(5,3)+n(5,4)+n(5,5);
    
    % Model equations
     ProdR1=  d.coff_a(1)*n_R1^3 + d.coff_b(1)*n_R1^2 +d.coff_c(1)*n_R1;
     ProdR2=  d.coff_a(2)*n_R2^3 + d.coff_b(2)*n_R2^2 +d.coff_c(2)*n_R2;
     ProdR3=  d.coff_a(3)*n_R3^3 + d.coff_b(3)*n_R3^2 +d.coff_c(3)*n_R3;
     ProdR4=  d.coff_a(4)*n_R4^3 + d.coff_b(4)*n_R4^2 +d.coff_c(4)*n_R4;
     ProdR5=  d.coff_a(5)*n_R5^3 + d.coff_b(5)*n_R5^2 +d.coff_c(5)*n_R5;
   
for I = 1:d.nb_R

        VI= find(d.nbrd(I,:)==1);
        %CAlculate the region production

    for J = 1:d.nb_R
                
 if I==J %Caculate exit flows  n_dot(I,I) (exit flow) and MIIH (Internal flows with boundary crossing)

             %MII (exit flows)  
                     switch I

                        case 1
                            M_II(1,1) =  n(1,1)/n_R1 * ProdR1 / TL(1,1)*ALP(1,1);  
                        case 2
                            M_II(2,2) =  n(2,2)/n_R2 * ProdR2 / TL(2,2)*ALP(2,2);
                        case 3
                            M_II(3,3) =  n(3,3)/n_R3 * ProdR3 / TL(3,3)*ALP(3,3);
                        case 4
                            M_II(4,4) =  n(4,4)/n_R4 * ProdR4 / TL(4,4)*ALP(4,4);
                        case 5
                            M_II(5,5) =  n(5,5)/n_R5 * ProdR5 / TL(5,5)*ALP(5,5);
                     end
                 
                                  
             else %Caculate transfer flows MIJH 
                    %MIJH 
                    for H = VI                 
                        switch I   
                             case 1
                                 M_1(J,H) =  Th_1(J,H)*n(I,J)/n_R1  * ProdR1 / TL(1,H)*ALP(1,H);
                             case 2
                                 M_2(J,H) =  Th_2(J,H)*n(I,J)/n_R2  * ProdR2 / TL(2,H)*ALP(2,H);
                             case 3
                                 M_3(J,H) =  Th_3(J,H)*n(I,J)/n_R3  * ProdR3 / TL(3,H)*ALP(3,H);
                             case 4
                                 M_4(J,H) =  Th_4(J,H)*n(I,J)/n_R4  * ProdR4 / TL(4,H)*ALP(4,H);
                             case 5
                                 M_5(J,H) =  Th_5(J,H)*n(I,J)/n_R5  * ProdR5 / TL(5,H)*ALP(5,H);
                        end
                    end

             end
     end
end
        
        
         % Calculate n_dot
for I =1:d.nb_R
    for J =1:d.nb_R
        if I==J %Exit flows

           % find incoming EXIT flows    
            switch I
                case 1
                    sum_M_HII= M_2(I,I) + M_3(I,I) + M_4(I,I)+ M_5(I,I);
                case 2
                    sum_M_HII= M_1(I,I) + M_3(I,I) + M_4(I,I);
                case 3
                    sum_M_HII= M_1(I,I) + M_2(I,I) + M_5(I,I);
                case 4
                    sum_M_HII= M_1(I,I) + M_2(I,I) + M_5(I,I);
                case 5 
                    sum_M_HII= M_1(I,I) + M_3(I,I) + M_4(I,I);
            end

            n_dot(I,I) = q(I,I) - M_II(I,I) + sum_M_HII;

        else

            % find incoming transfering flows
            switch I
                case 1
                    sum_M_HJI= M_2(J,I) + M_3(J,I) + M_4(J,I) + M_5(J,I);
                case 2
                    sum_M_HJI= M_1(J,I) + M_3(J,I) + M_4(J,I);
                case 3
                    sum_M_HJI= M_1(J,I) + M_2(J,I) + M_5(J,I);
                case 4
                    sum_M_HJI= M_1(J,I) + M_2(J,I) + M_5(J,I);
                case 5 
                    sum_M_HJI= M_1(J,I) + M_3(J,I) + M_4(J,I);
            end

            % find outgoing transfering flows 
             switch I   
                 case 1
                        sum_M_IJH = M_1(J,2) + M_1(J,3) + M_1(J,4) + M_1(J,5);
                 case 2
                        sum_M_IJH = M_2(J,1) + M_2(J,3) + M_2(J,4);
                 case 3 
                        sum_M_IJH = M_3(J,1) + M_3(J,2) + M_3(J,5);
                 case 4 
                        sum_M_IJH = M_4(J,1) + M_4(J,2) + M_4(J,5);
                 case 5
                        sum_M_IJH = M_5(J,1) + M_5(J,3) + M_5(J,4);
             end
             
             n_dot(I,J) = q(I,J) - sum_M_IJH + sum_M_HJI;
             
        end          
    end
end
        
     d.g = Function('f', ...
        {q,TL,n,ALP,Th_1,Th_2,Th_3,Th_4,Th_5}, ... % inputs
        {n_dot}); % output 
end

function d = create_prediction_integrator(d)
addpath('casadi/casadi-windows-matlabR2016a-v3.5.5')
import casadi.*
    
    % Initial accumulation states
    n0 = MX.sym('n0',d.nb_R,d.nb_R);
    
    % Demand flows
    q_I = MX.sym('q_I',d.nb_R,d.nb_R);
    
%     %Average trip length coefficient
     ALP_I = MX.sym('ALP_I',d.nb_R,d.nb_R);
    
    %Average trip length
     TL_I = MX.sym('TL_I',d.nb_R,d.nb_R);
    
        % Route guidance commands
    Th_1_I = MX.sym('Th_1_I',d.nb_R,d.nb_R);
    Th_2_I = MX.sym('Th_2_I',d.nb_R,d.nb_R);
    Th_3_I = MX.sym('Th_3_I',d.nb_R,d.nb_R);
    Th_4_I = MX.sym('Th_4_I',d.nb_R,d.nb_R);
    Th_5_I = MX.sym('Th_5_I',d.nb_R,d.nb_R);
    
    % Initialize state
    nI = n0;
    tpred = d.ts/d.is;
            
    % For loop over integrator steps
    switch d.pred_int
        
        case 'RK'
            
            for k_I = 1:d.is

                [k1] = d.g(q_I,TL_I,nI,ALP_I,...
                    Th_1_I,Th_2_I,Th_3_I,Th_4_I,Th_5_I);
                
                nk1 = nI +k1.*tpred/2;

                [k2] = d.g(q_I,TL_I,nk1,ALP_I,...
                    Th_1_I,Th_2_I,Th_3_I,Th_4_I,Th_5_I);
                
                nk2 = nk1 +k2.*tpred/2;
                
                [k3] = d.g(q_I,TL_I,nk2,ALP_I,...
                    Th_1_I,Th_2_I,Th_3_I,Th_4_I,Th_5_I);
                
                nk3 = nk2 +k3.*tpred;
                
                [k4] = d.g(q_I,TL_I,nk3,ALP_I,...
                    Th_1_I,Th_2_I,Th_3_I,Th_4_I,Th_5_I);
                
                nI = nI + tpred/6.*...
                        (k1 + 2.*k2 + 2.*k3 + k4);
            end
                    
       case 'Euler'
    
            for k = 1:d.is
                
                [k_I] = d.g(q_I,TL_I,nI,ALP_I,...
                    Th_1_I,Th_2_I,Th_3_I,Th_4_I,Th_5_I);

                nI = nI + tpred .* k_I;
            end
    end

 % Create prediction integrator function F_pred
    d.G_pred = Function('G_pred',{q_I,TL_I,n0,ALP_I,...
        Th_1_I,Th_2_I,Th_3_I,Th_4_I,Th_5_I}, ...
         {nI});
end

function d = formulate_NLP(d)
    addpath('casadi/casadi-windows-matlabR2016a-v3.5.5')
    import casadi.*
    
    % Initialize NLP arguments
    w = {}; % decision variables
    v = {}; % input parameters
    d.lbw = [];
    d.ubw = [];
    J = 0; % objective function
    g = {};
    d.lbg = [];
    d.ubg = [];

    % Initial accumulation states
    n0 = MX.sym('n0_1',d.nb_R,d.nb_R); 
    
    % Initialize states
    nk = n0;

    
    % Append n0k to NLP parameters vector v
    v = {v{:}, n0(:)}; 
    
    % decision variable for the Average trip length
      for k = 1:d.t_int

         % variable for the Average trip length
          TL = d.TLIH_t{k,1};
          
          % decision variable for the Average trip length coefficient
          
           if ismember(k,d.tds_update)==1
            % Decision variable to trip length 
              ALP  = MX.sym(['ALP_' num2str(k)],d.nb_R,d.nb_R);
            % Append Qbhat to NLP variables vector w
              w = {w{:}, ALP(:)};
           end
         
          %icoming flow at current time step
          qk = d.qs{k};
          
          %flow transfer at current time step
          Th_1k = d.Th_1{k};
          Th_2k = d.Th_2{k};
          Th_3k = d.Th_3{k};
          Th_4k = d.Th_4{k};
          Th_5k = d.Th_5{k};
          
          %Acuumulation states  
          ns_k = d.ns{k} ;
            
         [nk]= d.G_pred(qk,TL,nk,ALP,...
                 Th_1k,Th_2k,Th_3k,Th_4k,Th_5k); 
             
         J=J+ sum(((ns_k(:)-nk(:))).^2);

         cons6 = sum(nk(1,:));
         cons7 = sum(nk(2,:));
         cons8 = sum(nk(3,:));
         cons9 = sum(nk(4,:));
         cons10 = sum(nk(5,:));
         
         g = {g{:}, cons6(:), cons7(:),cons8(:),cons9(:),cons10(:)};
          
         lim_lg = [0;0;0;0;0];
         lim_ug = [d.R1_n_jam;d.R2_n_jam;d.R3_n_jam;d.R4_n_jam;d.R5_n_jam];
         
         d.lbg = [d.lbg;lim_lg];
         d.ubg = [d.ubg;lim_ug];
      end
       
   % Create NLP problem structure
    prob = struct(...
        'f', J,... % objective function
        'x', vertcat(w{:}),... % decision variables
        'p', vertcat(v{:}),... % parameters
        'g', vertcat(g{:})); % constraints
    
    % Create function for solving NLP
    switch d.NLP_solver
        
        case 'sqpmethod'
            
            options = struct;
            
            % Choose QP subproblem solver
            %options.qpsol = 'qpoases';
            
            % Choose Hessian update type
            options.hessian_approximation = ...
                'limited-memory'; % BFGS
            
            options.max_iter = d.max_iter_sqpmethod;
            
            d.f_NLP = nlpsol('f_NLP', ...
                'sqpmethod', prob, options);
            
        case 'ipopt'
            
            options = struct;
            
            options.ipopt.max_iter = d.max_iter_NLP;
            
            % Choose Hessian update type
            options.ipopt.hessian_approximation = ...
                'limited-memory'; % BFGS
            
            options.ipopt.tol = 1;
            options.ipopt.dual_inf_tol = 100;
            options.ipopt.compl_inf_tol =1;
            options.ipopt.constr_viol_tol = 1;
            
            options.ipopt.acceptable_iter = 3;
            options.ipopt.acceptable_tol = 1;
            
            %options.ipopt.acceptable_dual_inf_tol = 100;
            options.ipopt.acceptable_compl_inf_tol =100;
            options.ipopt.acceptable_constr_viol_tol = 100;
            
            options.ipopt.limited_memory_max_history = 60;
%             
            options.warn_initial_bounds = true;
            
            d.f_NLP = nlpsol('f_NLP', 'ipopt', prob, options);
            
    end

end

