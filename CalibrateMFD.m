function d = CalibrateMFD(d)
global  sim svdir macro

%Calculate transfer flows given by MFD dynamics

%load coefficient matrix
load (strcat(macro,'\MFD_coff_10.mat'));
 
%load n_jam 
load (strcat(macro,'\n_jam.mat'));

%load simulation based accumulation matrix
n_SIM = d.MCP{sim,1};
%  
%load Split Ratio matrix
Theta_T = d.MCP{sim,2}; 
 
%load q_sim 
q_mod = d.MCP{sim,4}; 

% %load trip lenghth matrix
% AVGTL = d.MCP{sim,3}; 

%load trip lenghth matrix
AVGTLIH = d.MCP{sim,3}; 


 %MFD parameters
 d.coff_a = coff_a; %veh.km/Hr -> veh.m/s
 d.coff_b = coff_b; %veh.km/Hr -> veh.m/s
 d.coff_c = coff_c; %veh.km/Hr -> veh.m/s
 
% Neighborhoods matrix
d.nbrd = [0 1 1 1 1;
          1 0 1 1 0;
          1 1 0 0 1;
          1 1 0 0 1;
          1 0 1 1 0];

 % Maximum acllowed accumulation
 d.R1_n_jam = n_jam_m(1);
 d.R2_n_jam = n_jam_m(2);
 d.R3_n_jam = n_jam_m(3);
 d.R4_n_jam = n_jam_m(4);
 d.R5_n_jam = n_jam_m(5);
 
ThI = cell(d.nb_R,1) ;
Th = cell(d.t_int,1) ;

for t = 1:d.t_int  
    NS = zeros(5,5);
    for i =1:d.nb_R 
          ThI{i} = Theta_T{t}(Theta_T{t}(:,1)==i,3:7);
         NS(i,:)   = sum(n_SIM{d.t_prd*t,i},2)';
    end 
    
    d.ns{t} = NS;
    Th{t} = ThI;
    d.Th_1{t} = Th{t}{1,1};
    d.Th_2{t} = Th{t}{2,1};
    d.Th_3{t} = Th{t}{3,1};
    d.Th_4{t} = Th{t}{4,1};
    d.Th_5{t} = Th{t}{5,1};
    
    %d.qs{t} = q_sim{t}./d.is;
    d.qs{t} = q_mod{t}./d.is;
    d.TLIH_t{t,1} = AVGTLIH{t,1};
    d.TLIH_t{t,1}(d.TLIH_t{t,1}==inf)=0;   
end


% d = build_calibration1(d.seed,d.t_prd, d);
%  [d,CPU_time,J_OPTM] = solve_calibration1(d);

d.tds = 3; %  decision variable is updated at every d.tds minutes 
d.tds_rate = d.tds/d.t_prd;
d.tds_update = (0:d.tds_rate:d.t_int)+1;

% d = build_calibration2_1(d.seed,d.t_prd, d);
% [d,CPU_time,J_OPTM] = solve_calibration2_1(d);

d = build_calibration2(d.seed,d.t_prd, d);
[d,CPU_time,J_OPTM] = solve_calibration2(d);


   %% Solve Optimization 
 Accum_OPT =zeros(d.t_int,d.nb_R );
 Accum_MFD =zeros(d.t_int,d.nb_R );
 Accum_SIM =zeros(d.t_int,d.nb_R );
 n_OPT = cell(d.t_int,1);
 
d.bef_flg=0; % Set to 0 as we usee time dependent Trip lengths

n_MFD = cell(d.t_int,1);
if d.bef_flg ==1 
            J_Before = 0;
              for t = 1:d.t_int

                   if t==1
                        %initial accumulation
                             nt_0 = ones(d.nb_R,d.nb_R).*1*10^(-190);
                  end%

                  q_t = d.qs{t};
                  Th_t = Th{t};

                  dTL = d.TL_t(t,:); 
                  alp = ones(d.nb_R,1);

                  nc= d.F_pred(q_t,dTL,nt_0,alp,...
                               Th_t{1},Th_t{2},Th_t{3},Th_t{4},Th_t{5});

                   n_MFD{t} = full(nc);

                  nt_0 = n_MFD{t}; % accumulation for the next step 

                  for i =1:d.nb_R 
                    Accum_MFD(t,i) = sum(n_MFD{t}(i,:));
                     J_Before = J_Before + sum((n_MFD{t}(i,:)-(sum(n_SIM{d.t_prd*t,i},2))').^2);
                  end  
              end
  
else
         J_Before = 0;
          for t = 1:d.t_int

              if t==1
                    %initial accumulation
                         nt_0 = ones(d.nb_R,d.nb_R).*1*10^(-190);
              end%
               %q_t = q{t};
              q_t = d.qs{t};
              Th_t = Th{t};

              dTL = d.TLIH_t{t,1}; 

              alp = ones(d.nb_R,d.nb_R);

              nc= d.G_pred(q_t,dTL,nt_0,alp,...
                           Th_t{1},Th_t{2},Th_t{3},Th_t{4},Th_t{5});

              n_MFD{t} = full(nc);

              nt_0 = n_MFD{t}; % accumulation for the next step 

              for i =1:d.nb_R 
                Accum_MFD(t,i) = sum(n_MFD{t}(i,:));
                J_Before  = J_Before  + sum((n_MFD{t}(i,:)-(sum(n_SIM{d.t_prd*t,i},2))').^2);

              end 

          end
end
  %%
  J_After = 0;
  for t = 1:d.t_int
      
      if t==1
            %initial accumulation
                 nt_0 = ones(d.nb_R,d.nb_R).*1*10^(-190);
      end%
       %q_t = q{t};
      q_t = d.qs{t};
      Th_t = Th{t};

      dTL = d.TLIH_t{t,1}; 
      
      alp = d.ALPIH_t{t,sim};
         
      nc= d.G_pred(q_t,dTL,nt_0,alp,...
                   Th_t{1},Th_t{2},Th_t{3},Th_t{4},Th_t{5});
               
      n_OPT{t} = full(nc);

      nt_0 = n_OPT{t}; % accumulation for the next step 
          
      for i =1:d.nb_R 
        Accum_OPT(t,i) = sum(n_OPT{t}(i,:));
         Accum_SIM(t,i) = sum(sum(n_SIM{d.t_prd*t,i}));
        J_After = J_After + sum((n_OPT{t}(i,:)-(sum(n_SIM{d.t_prd*t,i},2))').^2);

      end 

  end
  
  d.CLB{sim,1} = Accum_SIM;
  d.CLB{sim,2} = Accum_MFD;
  d.CLB{sim,3} = Accum_OPT;  
  
%%  
  close all

     fig = figure(8); % 
 for i =1: d.nb_R
     dt_sim = d.t_prd*60;
    subplot(2,3,i)
    xOPT = (1:size(Accum_OPT,1))'*d.ts -d.ts/2;
    xMFD = (1:size(Accum_MFD,1))'*d.ts -d.ts/2;
    xSIM = (1:size(Accum_SIM,1))'*dt_sim -dt_sim/2; %Sim Output are set to 60s time intervals 
    plot(xMFD,Accum_MFD(:,i),'b'); hold on
    plot(xSIM,Accum_SIM(:,i),'r'); 
    plot(xOPT,Accum_OPT(:,i),'g');
    plot((0:1)'*d.ts/2,[0; Accum_MFD(1,i)],'--b',[xMFD(end) xMFD(end)+d.ts],[Accum_MFD(end,i);0],'--b');
    plot((0:1)'*dt_sim/2,[0; Accum_SIM(1,i)],'--r',[xSIM(end) xSIM(end)+dt_sim],[Accum_SIM(end,i);0],'--r');
    plot((0:1)'*d.ts/2,[0; Accum_OPT(1,i)],'--g',[xOPT(end) xOPT(end)+d.ts],[Accum_OPT(end,i);0],'--m');
     set(0,'defaultTextInterpreter','latex');
    title(sprintf('%s%i%s','Accumulation from SIM, OPT and MFD for Region ',i,' '));
    xlabel('Time [s]');
    ylabel('Accumulation [veh]');
    legend('No-Calib','SIM','$$L_{IH}(t)$$','location','southeast','interpreter','latex'); 
     xlim([0 3600])
 end
 fig.WindowState = 'maximized';

 figname = sprintf('Calibration_sim_%i.bmp',sim);
 saveas(fig,fullfile(svdir, figname))
end

  