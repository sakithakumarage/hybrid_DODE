function d = run_simulator_CN(d,flag)
global sim root mfddir
AIMSpath = 'C:\Program Files\Aimsun\Aimsun Next 8.2\aconsole.exe';

 % get initial Information to run the AIMSUN
    fid=fopen(strcat(root,'scenarioInfo\scenario.txt'));
    simData=textscan(fid,'%u %u %s %s %s %s');
    %replID=simData{1,1};
    replID=15058712 %15058724
    dbID=simData{1,2};
%     DBname=cell2mat(simData{1,3});
    DBname='vitoria2.sqlite'
    %angName=cell2mat(simData{1,4});
    angName = 'vitoria2.ang'
    pyPath=cell2mat(simData{1,5});
    assignMatrName=cell2mat(simData{1,6});
    fclose(fid);
    detecPath=strcat('scenarioInfo\detectors.txt'); 
    data=[];
%     dir = 'H:\OD_Estimation_model\testbeddata\Aimsun';
    
    
    %% Run the AIMSUN model
    if flag==1
        delete(DBname);
    %RUN AIMSUN For Windows Users
        commandTerminal= horzcat('"',AIMSpath,'" -script "',cd,'\',pyPath,'" "',cd,'\',angName,'" ',num2str(replID),' ',num2str(dbID),' ',cd,'\',detecPath);
        system(commandTerminal);
        delete('endOfSim.txt');
        delete(strcat(angName,'.old'));
    % After the execution, simulation outputs are taken from the sqlite db
    
    end
   conn = sqlite(DBname);
    
    % Retrieve detectors' IDs and define the SQL Query for outputs'
    sqlQuery='SELECT did, oid, sid, ent, countveh, speed, occupancy, density FROM MEDETECT WHERE did = 4000;';
    out2= fetch(conn,sqlQuery);
    out2= cellfun(@(x) double(x), out2);
    out2= out2(out2(:,3)==1,:);
    out2 = out2(out2(:,4)~=0,:);
    
    close(conn)

    %% Centroid level OD estimation variables 
    %generate detector counts
    [d.CNP{sim,1},d.CNP{sim,3}] = generate_detector_counts(out2);
    %generate detector densities
    [d.CNP{sim,4},d.CNP{sim,6}] = generate_detector_density(out2);
    if sim ==1 %Load only once as it doesnt change over itarations 
    d.CNP{sim,2} = generate_path_assignment(d,1);
    else
    generate_path_assignment(d,1);%Create only the figure 
    end
    %d.CNP{sim,1} = Simulated traffic counts
    %d.CNP{1,2} = Assignment matrix
    %d.CNP{1,3} = Observed traffic counts 
    
end