function [P_OD,O_OD,RP_OD,RO_OD,CentReg] = purturb_OD(ODPattern,Origins,Destinations,k,flag) 
global root testbed
load(strcat(testbed,'\GIS\centroids.txt'));
ODPattern(:,5:9) = zeros(size(ODPattern,1),5);

P_OD =cell(size(ODPattern,2),1);
O_OD =cell(size(ODPattern,2),1);
RP_OD = cell(size(ODPattern,2),1);
RO_OD = cell(size(ODPattern,2),1);
    for j=1:size(ODPattern,2)
        % Reshape each vector into an OD matrix
        m=reshape(ODPattern(:,j),length(Origins),length(Destinations)).*1;
        O_OD{j} =m;
        %Purturb the OD
         % Select the Purturbation function 
        if flag ==1 
            for p =1:size(m,1)
                for q =1:size(m,2)
                    switch k
                        case 1
                        %method 1 : D7
                        factorMean = 0.7; 
                        factor = factorMean + 0.3 * rand; 
                        m(p,q)=m(p,q)*factor;
                        case 2
                        % Method 2: D9
                        factorMean = 0.9;
                        factor = factorMean + 0.3 * rand; 
                        m(p,q)=m(p,q)*factor;
                        case 3
                        % Method 2: D11
                        factorMean = 1.1;
                        factor = factorMean + 0.3 * rand; 
                        m(p,q)=m(p,q)*factor;
                        case 4
                        % Method 2: D10
                        factorMean = 1.0;
                        factor = factorMean + 0.3 * rand; 
                        m(p,q)=m(p,q)*factor;
                        case 5
                        %Method 2: D14
                        factorMean = 0.25;
                        factor = factorMean + 2 * rand; 
                        m(p,q)=m(p,q)*factor;
                        case 6
                        %Method 4: trips to region one under estimated by %
                        %25 percent
                        if ismember(centroids(p,2),[2;3;4;5]) && ismember(centroids(q,2),1)
                               factorMean = 0.6;
                               factor = factorMean + 0.3 * rand; 
                               m(p,q)=m(p,q)*factor;
                        elseif ismember(centroids(p,2),[2;3;4;1]) && ismember(centroids(q,2),5)
                               factorMean = 0.9;
                               factor = factorMean + 0.3 * rand; 
                               m(p,q)=m(p,q)*factor;
                        end
                    end
                    
                end
            end
        else
            factor =1;
            m = m.*factor;
        end
        
        P_OD{j} =m;

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
        
        % create regional ODs 
        RP_OD{j} = zeros(5);
        RO_OD{j} = zeros(5);
        for p = 1:length(Origins)
            for q =1: length(Destinations)
                org_R = centroids(centroids(:,1)==Origins(p,1),2);
                des_R = centroids(centroids(:,1)==Destinations(q,1),2);
                %Perturbed regional ODs 
                RP_OD{j}(org_R,des_R) = RP_OD{j}(org_R,des_R)+P_OD{j}(p,q);
                %Original regional ODs
                RO_OD{j}(org_R,des_R) = RO_OD{j}(org_R,des_R)+O_OD{j}(p,q);
            end
        end
        
    end 
    
    maxod =0;
    figure(34)
        subplot(2,1,1)
        for k = 1: size(RO_OD,1)-5
        plot(RO_OD{k}(:),RP_OD{k}(:),'.b'); hold on 
        maxod = max([RO_OD{k}(:);RP_OD{k}(:);maxod]);
        statR(k,1) = mean(abs(RO_OD{k}(:)-RP_OD{k}(:)));
        statR(k,2) = std(abs(RO_OD{k}(:)-RP_OD{k}(:)));
        end
        plot(0:1:maxod,0:1:maxod,'r');
        title('Variation in regional level demands for all reagions');
        xlabel('Ground truth demand [veh/hr]');
        ylabel('Perturbed demand [veh/hr]');

        maxod =0;
        subplot(2,1,2)
        for k = 1: size(O_OD,1)-5
        plot(O_OD{k}(:),P_OD{k}(:),'.b'); hold on 
        ydat=P_OD{k}(:); 
        xdat=O_OD{k}(:);
        subplot(2,1,2)
        maxod= max([xdat(:);ydat(:);maxod]);
        plot(xdat,ydat,'.b'); hold on
        end
        plotmax= 1:maxod;
        plot(plotmax,plotmax,'r');%Reference Line
        mdl = fitlm(xdat,ydat);
        plotmaxY = plotmax.*mdl.Coefficients.Estimate(2)+ones(1,length(plotmax)).*mdl.Coefficients.Estimate(1);
        plot(plotmax,plotmaxY,'k','LineWidth',1);
        legend('Data','y=x',sprintf('Fit: y= %.4f*x  %.4f',mdl.Coefficients.Estimate(2),mdl.Coefficients.Estimate(1)),'location','southeast','EdgeColor',[1 1 1]);
        text(10,plotmax(end)-80,sprintf('$$R^2 = %.4f $$',mdl.Rsquared.Adjusted),'interpreter','latex');
        title('Variation in link level demands for all reagions');
        xlabel('Ground truth demand [veh/hr]');
        ylabel('Perturbed demand [veh/hr]');
      
%regional and centroid connection matrix
  load(strcat(testbed,'\GIS\CentReg.mat'));

end