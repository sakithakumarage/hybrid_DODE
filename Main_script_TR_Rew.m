 close all
 clearvars -except Out
%Open Testbed data
%restoredefaultpath
global sim svdir mfddir root testbed macro

% %Office PC 
testbed = 'C:\Users\s4544191\cloudstor\HomePC\testbeddata';
root = 'C:\Users\s4544191\cloudstor\HomePC\regional\';
macro = 'C:\Users\s4544191\cloudstor\HomePC\macro_static';

%Home PC 
% testbed = 'C:\Users\Sakitha\Documents\ODEStimationMatlab\HomePC\testbeddata';
% root = 'C:\Users\Sakitha\Documents\ODEStimationMatlab\HomePC\regional\';
% macro = 'C:\Users\Sakitha\Documents\ODEStimationMatlab\HomePC\macro_static';

svdir = strcat(root,'OD_Estimation_Output\');
mfddir = strcat(root,'MFD\');

load (strcat(testbed,'/Demand.mat'));
load (strcat(testbed,'/True_data.mat'));
load (strcat(testbed,'/Aimsun/scenarioInfo/Origins.txt'));
load (strcat(testbed,'/Aimsun/scenarioInfo/Destinations.txt'));
load (strcat(testbed,'/Aimsun/scenarioInfo/detectors.txt'));

%%
 GG = cell(10,1);
 for k =1:10

    %testbed analytics parameters
     hrs =1; %Simulation time available for OD estimation. 
     d.seed =hrs; %todetermine the one hour or two hours
     d.t_f =3600 *d.seed; % final time       
     d.t_prd = 3; %time stepin minutes
     d.dt =d.t_prd*60; %accumulation aggregation interval 
     d.t_int =d.t_f/d.dt; % Number of time intervals 
     d.nb_R =5; % Number of regions 
     d.ts = d.t_prd*60; % time step [s]
     d.is = d.t_prd*60;% integration step
    
    
    % Original OD and Purturbed OD   
    [d.Prt_OD,d.Org_OD,d.Prt_ROD,d.Org_ROD,d.CentReg] = purturb_OD(ODPattern,Origins,Destinations,4,1);
    
    %Observed accumulation data
    load (strcat(testbed,'\N_Obs.mat'));
    
    %GroundTruth regional OD 
    load(strcat(testbed,'\Q_Pri.mat'));

    for sim =1:10
        %Run the simulation with purturbed OD data and save data 
        d = run_simulator_RG(d,1); %Make second variable =1 of you want to run Aimsun 
    
        %Calibration of MFDs using optimization framework and getting calibrated data
        d = CalibrateMFD(d);
    
        %Run the simulation with purturbed OD data and save data     
        d = run_simulator_CN(d,1);
        
        %run OD estimation
        d = build_ODestimation_combined3(d); 
        d = solve_ODestimation_combined3(d);
    end
    
    %Runn the LSQR on same purturbed OD
    ODPattern_C = cat(2,d.Prt_OD{1, 1}(:),d.Prt_OD{2, 1}(:), d.Prt_OD{3, 1}(:),d.Prt_OD{4, 1}(:));   
    purturb_OD(ODPattern_C,Origins,Destinations,4,0);
    
    for sim = 11
        
        d = run_simulator_RG(d,1);
        
        d = run_simulator_CN(d,1);
        
        d = linkODdist_lsqr_2(d,1);
    end 
     
    GG{k,1} = d;
    filename = sprintf('%s%i%s','OD_Estimation_Output/TR_Rev/OutGGTest_',k,'.mat');
    save(filename,'GG')
 end