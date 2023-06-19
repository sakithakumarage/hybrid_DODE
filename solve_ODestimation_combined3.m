function d = solve_ODestimation_combined3(d)
global mfddir sim
% Initial Accumulation 
n0 = ones(d.nb_R,d.nb_R).*1*10^(-190);

% Demand matrix 
     q_pri =  d.MCP{1,5}(:,5:d.t_int+4);

        q_1 = zeros(d.org*d.dst,1); 
        q_2 = zeros(d.org*d.dst,1);
        q_3 = zeros(d.org*d.dst,1);
        q_4 = zeros(d.org*d.dst,1);
        
pram = [n0(:);q_1(:);q_2(:);q_3(:);q_4(:)];

d.lbw =[];
d.ubw =[];
q_llim = 0.75;
q_ulim = 1.25;

d.lbw =  [d.lbw;q_pri(:)*q_llim];
d.ubw =  [d.ubw;q_pri(:)*q_ulim];

J = manual_ODE_combined3(q_pri(:),pram,d)

tic
% Solve NLP
    d.NLP_solution = d.f_NLP('x0', ... % initial guess
        vertcat(q_pri(:)), ...
        'p', ... % parameters
        vertcat(pram(:)),...
        'lbx', ... % lower bounds for box constraints
        vertcat(d.lbw(:)), ...
        'ubx', ... % upper bounds for box constraints
        vertcat(d.ubw(:)), ...
        'lbg', ... % lower bounds for general constraints
        vertcat(d.lbg(:)), ...
        'ubg', ... % upper bounds for general constraints
        vertcat(d.ubg(:)));
  
 CPU_time = toc  
 
  qsol = full(d.NLP_solution.x);
  for t =1:d.t_int
      qk = qsol((t-1)*d.org*d.dst+1:t*d.org*d.dst,1);
     m =  floor((t-1)/5)+1;   
        switch m
           case 1
               q_1 = q_1 + qk;   
           case 2
              q_2 = q_2 + qk;
           case 3
              q_3  = q_3 + qk;
           case 4
               q_4 = q_4 + qk;
        end
  end
  d.MIN{sim,2} = qsol;
  d.MIN{sim,1}(:,1) = cat(1,q_1,q_2,q_3,q_4);
    x0 = [d.Org_OD{1,1}(:);d.Org_OD{2,1}(:);d.Org_OD{3,1}(:);d.Org_OD{4,1}(:)];
    d.MIN{sim,1}(:,2)=x0;
    
 d= writematrix(d,sim);
 
     xdat= d.MIN{sim,1}(:,2);
     ydat = d.MIN{sim,1}(:,1);
    plotmax= 1:max([xdat(:);ydat(:)]);
    plot(xdat,ydat,'.b'); hold on 
    plot(plotmax,plotmax,'r');%Reference Line
    mdl = fitlm(xdat,ydat,'Intercept',false);
    plot(plotmax,plotmax.*mdl.Coefficients.Estimate,'k','LineWidth',1);
    legend('Data','y=x',sprintf('Fit: y= %g*x',mdl.Coefficients.Estimate),'location','southeast');
    text(100,plotmax(end)-200,sprintf('$$R^2 = %g $$',mdl.Rsquared.Adjusted),'interpreter','latex');
    title('Variation in regional demands for all reagions','interpreter','latex');
    xlabel('Ground truth OD [veh/15min]','interpreter','latex');
    ylabel('Optimized OD-Regional [veh/15min]','interpreter','latex');
 
end


function d= writematrix(d,sim)

global root testbed
load (strcat(testbed,'/Aimsun/scenarioInfo/Origins.txt'));
load (strcat(testbed,'/Aimsun/scenarioInfo/Destinations.txt'));
%Add nearzero unassigned demand valees to matrix 
OD_Pattern =reshape(d.MIN{sim,1}(:,1),[d.org*d.dst,d.t_link]);
OD_Pattern(OD_Pattern==1*10^-64)=0;

%Add zero demand steps for simulator
OD_Pattern(:,5:9) = zeros(size(OD_Pattern,1),5);

for j=1:size(OD_Pattern,2)
     m=reshape(OD_Pattern(:,j),d.org,d.dst);
     d.Prt_OD{j,sim+1}=m;
             %Save the matrix into a .txt file compliant AIMSUN standards
        filename=strcat(strcat(root,'scenarioInfo\matrix\matrix',num2str(j),'.txt'));    
        fid=fopen(filename,'w');
        fprintf(fid,'id\t');
        fprintf(fid,'%i\t',Destinations);
        fprintf(fid,'\n');
        fclose(fid);
        fid=fopen(filename,'a');
        for i=1:length(Origins)
            fprintf(fid,'%i\t',Origins(i));
            fprintf(fid,'%5.2f\t',m(i,:));
            fprintf(fid,'\n');
        end
        fclose(fid);
end  
end