% To plot deaths averted under different priority scenarios run this code
% for three different R0 values and save the data as below:
% Out_1p25.m (for R0 =1.25); Out_2.m (for R0 =2) and Out_2p5.m (for R0 = 2.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (First set up the model for Infection preventing vaccine or disease preventing vaccine then run this code) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; load Model_setup_India

% --- Make samples for each parameter -------------------------------------

nsam = 250;
xs = repmat(prm.bounds(1,:),nsam,1) + repmat(diff(prm.bounds,1),nsam,1).*lhsdesign(nsam,size(prm.bounds,2));

opts = odeset('NonNegative',[1:i.nstates], 'Refine', 64, 'AbsTol', 1e-10, 'RelTol', 1e-10);

% Choose as required
R0  = 2.5;     % 1.25 (Out_1p25.m), 2.0 (Out_2.m), 2.5 (Out_2p5.m) 
tf  = 1e3;     % Simulation duration

mk = round(nsam/20);
 
% First find the baseline scenario
    
    for ii = 1:nsam
        
        if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
            
        [p,r] = alloc_parameters(xs(ii,:), xi, p, r, R0, i, s, gps, prm);
                 
        p1 = p; r1 = r;
        p1.vacc1  = [0, 0, 0]; 
        p1.vacc2  = [0, 0, 0]; 
        r1.init = get_init(p1, r1, i, s, gps, prm);

        % --- Perform the simulation
        M1 = make_model2(p1, r1, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M1, i, s, p1, r1, agg, sel, prm);
        [t,soln00] = ode15s(geq, [0:1:tf], r1.init, opts);
        
        % --- Record cumulative deaths (total)
        mor00(:,ii) = sum(soln00(:,i.aux.mort),2); 
                                  
    end
    fprintf('\n');
    
    
 % Stepwise priority groups (PG)in uniform strategy
 stages = [1,2,3];  % Priority groups (PG)
 for jj = 1:length(stages)
     scenario = stages(jj);
      
    pvs = 1; %[1 2]/2; % Vaccination coverages in different groups       
   
for ip = 1:length(pvs)
    
     % CE-> Comorbidity then Elderly; EC-> Elderly then Comorbidity 
    % Uniform vaccination scenarios (in strategy 1)
         PG01_CE = [[0,0.0448,0];[0,0.02,0.4474];[0,0.0448,0.9999]];  % Non risk group
         PG02_CE = [[0,0.0448,0];[0,0.4474,0.4474];[0,0.9999,0.9999]]; % Risk group 
          
     % Priority based Vaccination scenarios (in strategy 1)      
         PG1_CE = [[0,0.0448,0];[0,0.0448,0];[0,0.0448,0.9999]];
         PG2_CE = [[0,0.0448,0];[0,0.9999,0];[0,0.9999,0.9999]];
 
     % Uniform vaccination scenarios (in strategy 2)
        PG01_EC = [[0,0.0448,0];[0,0.0297,0.6628];[0,0.0448,0.9999]];  % Non risk group
        PG02_EC = [[0,0.0448,0];[0,0.6628,0.6628];[0,0.9999,0.9999]];   % Risk group 
  
     % Priority based Vaccination scenarios (in strategy 2) 
        PG1_EC = [[0,0.0448,0];[0,0.0448,0.9999];[0,0.0448,0.9999]];
        PG2_EC = [[0,0.0448,0];[0,0.0448,0.9999];[0,0.9999,0.9999]];

        % CE-> Comorbidity then Elderly
        scen_nonrisk01 = [PG01_CE(1,:)*pvs(ip);PG01_CE(1,:)+(PG01_CE(2,:)-PG01_CE(1,:)).*pvs(ip);PG01_CE(2,:)+(PG01_CE(3,:)-PG01_CE(2,:)).*pvs(ip)];
        scen_risk01    = [PG02_CE(1,:)*pvs(ip);PG02_CE(1,:)+(PG02_CE(2,:)-PG02_CE(1,:)).*pvs(ip);PG02_CE(2,:)+(PG02_CE(3,:)-PG02_CE(2,:)).*pvs(ip)];
                      
        scen_nonrisk1 = [PG1_CE(1,:)*pvs(ip);PG1_CE(1,:)+(PG1_CE(2,:)-PG1_CE(1,:)).*pvs(ip);PG1_CE(2,:)+(PG1_CE(3,:)-PG1_CE(2,:)).*pvs(ip)];
        scen_risk1    = [PG2_CE(1,:)*pvs(ip);PG2_CE(1,:)+(PG2_CE(2,:)-PG2_CE(1,:)).*pvs(ip);PG2_CE(2,:)+(PG2_CE(3,:)-PG2_CE(2,:)).*pvs(ip)];
   
        % EC-> Elderly then Comorbidity 
        scen_nonrisk02 = [PG01_EC(1,:)*pvs(ip);PG01_EC(1,:)+(PG01_EC(2,:)-PG01_EC(1,:)).*pvs(ip);PG01_EC(2,:)+(PG01_EC(3,:)-PG01_EC(2,:)).*pvs(ip)];
        scen_risk02    = [PG02_EC(1,:)*pvs(ip);PG02_EC(1,:)+(PG02_EC(2,:)-PG02_EC(1,:)).*pvs(ip);PG02_EC(2,:)+(PG02_EC(3,:)-PG02_EC(2,:)).*pvs(ip)];
                      
        scen_nonrisk2 = [PG1_EC(1,:)*pvs(ip);PG1_EC(1,:)+(PG1_EC(2,:)-PG1_EC(1,:)).*pvs(ip);PG1_EC(2,:)+(PG1_EC(3,:)-PG1_EC(2,:)).*pvs(ip)];
        scen_risk2    = [PG2_EC(1,:)*pvs(ip);PG2_EC(1,:)+(PG2_EC(2,:)-PG2_EC(1,:)).*pvs(ip);PG2_EC(2,:)+(PG2_EC(3,:)-PG2_EC(2,:)).*pvs(ip)];

        pop = prm.N; %(prm.N*1380004385/1e6); 
        pop_gen = pop.*(1-p.risk);
        pop_risk = pop.*p.risk;
       
        popvac_gen1 = sum(pop_gen.*scen_nonrisk01(scenario,:));
        popvac_risk1 = sum(pop_risk.*scen_risk01(scenario,:));
        pcov_CE(:,ip,jj) = 100*((popvac_gen1+popvac_risk1)./sum(pop)); % percentage coverage 
        
        popvac_gen2 = sum(pop_gen.*scen_nonrisk02(scenario,:));
        popvac_risk2 = sum(pop_risk.*scen_risk02(scenario,:));
        pcov_EC(:,ip,jj) = 100*((popvac_gen2+popvac_risk2)./sum(pop)); % percentage coverage 
        
    
    for ii = 1:nsam
        
        if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
            
        [p,r] = alloc_parameters(xs(ii,:), xi, p, r, R0, i, s, gps, prm);       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % CE-> Comorbidity then Elderly
        p1 = p; r1 = r;
        p1.vacc1 = scen_nonrisk01(scenario,:);
        p1.vacc2 = scen_risk01(scenario,:);  
        r1.init = get_init(p1, r1, i, s, gps, prm);
        % --- Perform the simulation
        M1 = make_model2(p1, r1, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M1, i, s, p1, r1, agg, sel, prm);
        [t,soln0] = ode15s(geq, [0:1:tf], r1.init, opts);
        
        % --- Record cumulative deaths (total)   
        mor01(:,ii,ip) = sum(soln0(:,i.aux.mort),2);  % Uniform strategy
        
        p2 = p; r2 = r;
        p2.vacc1 = scen_nonrisk1(scenario,:);
        p2.vacc2 = scen_risk1(scenario,:);  
        r2.init = get_init(p2, r2, i, s, gps, prm);
        % --- Perform the simulation
        M2 = make_model2(p2, r2, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M2, i, s, p2, r2, agg, sel, prm);
        [t,soln1] = ode15s(geq, [0:1:tf], r2.init, opts);
        
        % --- Record cumulative deaths (total)
        mor1(:,ii,ip) = sum(soln1(:,i.aux.mort),2); % Other vaccination strategy
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % EC-> Elderly then Comorbidity
        p3 = p; r3 = r;
        p3.vacc1 = scen_nonrisk02(scenario,:);
        p3.vacc2 = scen_risk02(scenario,:);  
        r3.init = get_init(p3, r3, i, s, gps, prm);
        % --- Perform the simulation
        M3 = make_model2(p3, r3, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M3, i, s, p3, r3, agg, sel, prm);
        [t,soln2] = ode15s(geq, [0:1:tf], r3.init, opts);
        
        % --- Record cumulative deaths (total)   
        mor02(:,ii,ip) = sum(soln2(:,i.aux.mort),2);  % Uniform strategy
        
        p4 = p; r4 = r;
        p4.vacc1 = scen_nonrisk2(scenario,:);
        p4.vacc2 = scen_risk2(scenario,:);  
        r4.init = get_init(p4, r4, i, s, gps, prm);
        % --- Perform the simulation
        M4 = make_model2(p4, r4, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M4, i, s, p4, r4, agg, sel, prm);
        [t,soln3] = ode15s(geq, [0:1:tf], r4.init, opts);
        
        % --- Record cumulative deaths (total)
        mor2(:,ii,ip) = sum(soln3(:,i.aux.mort),2); % Other vaccination strategy
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    fprintf('\n');
    
end

 for j = 1:length(pvs)
        % CE-> Comorbidity then Elderly
        %pop = (prm.N*1380004385/1e6);                                      % Actual population in three age group
        DA0_CE = (mor00(end,:,1) - mor01(end,:,j));%.*(1380004385/1e6);     % Total DA by uniform strategy
        DA_CE = (mor01(end,:,j) - mor1(end,:,j));%.*(1380004385/1e6);       % Total DA different strategies relative to unifor strategy
        DA1_CE = (mor00(end,:,j) - mor1(end,:,j));%.*(1380004385/1e6);      % Total DA different strategies relative zero vaccine strategy
 
        % Find percentiles for DA (Total)  
        DA0prc_CE(j,:,jj) = prctile(DA0_CE,[2.5,50,97.5],2);%./1e6;  % Uniform strategy
        DAprc_CE(j,:,jj) = prctile(DA_CE,[2.5,50,97.5],2);%./1e6;    % Priority strategy relative to uniform strategy
        DA1prc_CE(j,:,jj) = prctile(DA1_CE,[2.5,50,97.5],2);%./1e6;  % Priority strategy relative to non-vaccinated scenario
        
        % EC-> Elderly then Comorbidity
        DA0_EC = (mor00(end,:,1) - mor02(end,:,j));%.*(1380004385/1e6);     % Total DA by uniform strategy
        DA_EC = (mor02(end,:,j) - mor2(end,:,j));%.*(1380004385/1e6);       % Total DA different strategies relative to unifor strategy
        DA1_EC = (mor00(end,:,j) - mor2(end,:,j));%.*(1380004385/1e6);      % Total DA different strategies relative zero vaccine strategy
 
        % Find percentiles for DA (Total)  
        DA0prc_EC(j,:,jj) = prctile(DA0_EC,[2.5,50,97.5],2);%./1e6;  % Uniform strategy
        DAprc_EC(j,:,jj) = prctile(DA_EC,[2.5,50,97.5],2);%./1e6;    % Priority strategy relative to uniform strategy
        DA1prc_EC(j,:,jj) = prctile(DA1_EC,[2.5,50,97.5],2);%./1e6;  % Priority strategy relative to non-vaccinated scenario 
        
 end
 end
    
DA0prc0_CE =[[0 0 0];DA0prc_CE(:,:,1);DA0prc_CE(:,:,2);DA0prc_CE(:,:,3)];  % [0, PG1, PG2, PG3] % Uniform strategy
DAprc0_CE =[[0 0 0];DAprc_CE(:,:,1);DAprc_CE(:,:,2);DAprc_CE(:,:,3)];      % [0, PG1, PG2, PG3] % Priority strategy relative to uniform strategy
DA1prc0_CE =[[0 0 0];DA1prc_CE(:,:,1);DA1prc_CE(:,:,2);DA1prc_CE(:,:,3)];  % [0, PG1, PG2, PG3] % Priority strategy relative to non-vaccinated scenario

DA0prc0_EC =[[0 0 0];DA0prc_EC(:,:,1);DA0prc_EC(:,:,2);DA0prc_EC(:,:,3)];  % [0, PG1, PG2, PG3] % Uniform strategy
DAprc0_EC =[[0 0 0];DAprc_EC(:,:,1);DAprc_EC(:,:,2);DAprc_EC(:,:,3)];      % [0, PG1, PG2, PG3] % Priority strategy relative to uniform strategy
DA1prc0_EC =[[0 0 0];DA1prc_EC(:,:,1);DA1prc_EC(:,:,2);DA1prc_EC(:,:,3)];  % [0, PG1, PG2, PG3] % Priority strategy relative to non-vaccinated sECnario

pcov0_CE = [0,pcov_CE(:,:,1),pcov_CE(:,:,2),pcov_CE(:,:,3)]';          % [0, PG1, PG2, PG3]
pcov0_EC = [0,pcov_EC(:,:,1),pcov_EC(:,:,2),pcov_EC(:,:,3)]';          % [0, PG1, PG2, PG3]
     
DA0choice_CE  = [pcov0_CE,DA0prc0_CE]; % Uniform strategy
DAchoice_CE   = [pcov0_CE,DAprc0_CE];  % Priority strategy relative to uniform strategy
DA1choice_CE  = [pcov0_CE,DA1prc0_CE]; % Priority strategy relative to non-vaccinated scenario
  
DA0choice_EC  = [pcov0_EC,DA0prc0_EC]; % Uniform strategy
DAchoice_EC   = [pcov0_EC,DAprc0_EC];  % Priority strategy relative to uniform strategy
DA1choice_EC  = [pcov0_EC,DA1prc0_EC]; % Priority strategy relative to non-vaccinated scenario

save('Out_1p25.m') 




