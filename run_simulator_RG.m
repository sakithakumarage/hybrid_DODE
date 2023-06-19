function d = run_simulator_RG(d,flag)
global sim root mfddir testbed
AIMSpath = 'C:\Program Files\Aimsun\Aimsun Next 8.2\aconsole.exe';

 % get initial Information to run the AIMSUN
    fid=fopen(strcat(root,'scenarioInfo\scenario.txt'));
    simData=textscan(fid,'%u %u %s %s %s %s');
    %replID=simData{1,1};
    replID=15058772 %15058724
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
    % Retrieve detectors' IDs and define the SQL Query for outputs'
    end
   conn = sqlite(DBname);
   
    sqlQuery='SELECT oid, did, sid, ent, input_count, vIn, vOut, vWait, traveltime, travel, totalDistanceTraveledInside, totalTravelTimeInside FROM MESYS ORDER BY oid,ent;';
    out4= fetch(conn,sqlQuery);
    out4= cellfun(@(x) double(x), out4);
    out4= out4(out4(:,3)==1,:);
    
    if out4(end,5)==0
        sprintf('Network is empty at end of simulation')
    else
        sprintf('Network is NOT empty at end of simulation')
        pause
    end
    
   
    sqlQuery='SELECT oid, ent, SectionId, exitTime, travelTime, delayTime FROM MEVEHSECTTRAJECTORY ORDER BY oid,ent;';
    out3= fetch(conn,sqlQuery);
    out3= cellfun(@(x) double(x), out3);
    
    sqlQuery='SELECT oid, did, sid, entranceTime, exitTime, travelledDistance, speed, origin, destination  FROM VEHTRAJECTORY ORDER BY oid,sid;';
    out5= fetch(conn,sqlQuery);
    out5= cellfun(@(x) double(x), out5);
    

    sqlQuery='SELECT oid, did, sid, ent, destination, nbveh, flow, input_count FROM MECENT_O ORDER BY oid,sid;';
    out7= fetch(conn,sqlQuery);
    out7= cellfun(@(x) double(x), out7);
    out7= out7(out7(:,3)==1,:);
    out7 = out7(out7(:,4)~=0,:);
    out7 = out7(out7(:,5)~=0,:);
    
%     sqlQuery='SELECT oid, sid, ent, countveh, speed FROM MEDETECT ORDER BY oid,sid;';
%     out2= fetch(conn,sqlQuery);
%     out2= cellfun(@(x) double(x), out2);
%     out2= out2(out2(:,2)==1,:);
%     out2 = out2(out2(:,3)~=0,:);
    
    close(conn)
    
    TRJ{1,2} =out5;
    TRJ{1,4} =out4;
    TRJ{1,7} =out7;
    
   %Adjust Turn lengths of trajectories with Global trajectories
    [out3,~] = ProcessTrajectory(d,out3,out5);
    
    TRJ{1,1} = out3;
    save(fullfile(mfddir,'TRJ.mat'),'TRJ');
    
    %% Macroscopic variables 

    %generate simulation accumulation profiles
    d.MCP{sim,1} = generate_sim_accum (d,out3,out5);
    
    %generate split ratios 
    d.MCP{sim,2} = generate_split_ratios(d,out3);
    
    %generate averate trip lengths 
    d.MCP{sim,3} = generate_average_triplength(d,out3);
    
    %generate simulation demand profiles
    [d.MCP{sim,4},d.MCP{sim,5}] = generate_sim_demand(d,out3,out7);
    
    
    
    
       
end