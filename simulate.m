% Propagate uncertainty to incidence projections
% First run "Setup_model.m" by choosing appropriate values for p.c1 and p.c3 run this code and save the data as described in line 73-75.  

clear all; load Model_setup_India

% --- Make samples for each parameter -------------------------------------

nsam = 250;
xs = repmat(prm.bounds(1,:),nsam,1) + repmat(diff(prm.bounds,1),nsam,1).*lhsdesign(nsam,size(prm.bounds,2));

p.seropos = 0.071*ones(1,3); % To adjust current sero-positivity in India

opts = odeset('NonNegative',[1:i.nstates], 'Refine', 64, 'AbsTol', 1e-10, 'RelTol', 1e-10);

R0  = 2;
tf  = 2e3;          % Simulation duration
pvs = [0 1];        % Vaccination coverages, whole popn        

mk = round(nsam/20);

for ip = 1:length(pvs)
    
    for ii = 1:nsam
        
        if mod(ii,mk)==0; fprintf('%0.5g ', ii/mk); end
            
        [p,r] = alloc_parameters(xs(ii,:), xi, p, r, R0, i, s, gps, prm);
              
        % Vaccination scenarios in different age groups for figure 3 
        scen_nonrisk = [[0,0.045,0];[0,0.045,0];[0,0.045,0.9999]];   % proporion coverage among [HCW+FW, HCW+FW+Co-morbid, HCW+FW+Co-morbid+Elderly]
        scen_risk    = [[0,0.045,0];[0,0.9999,0];[0,0.9999,0.9999]]; % proporion coverage among [HCW+FW, HCW+FW+Co-morbid, HCW+FW+Co-morbid+Elderly]
 
        scenario = 1; % Choose the scenario 1: HCW+FW; 2:HCW+FW+Co-morbid; 3:HCW+FW+Co-morbid+Elderly 
        
        p1 = p; r1 = r;
        p1.vacc1 = pvs(ip)*scen_nonrisk(scenario,:);  
        p1.vacc2 = pvs(ip)*scen_risk(scenario,:);  
        r1.init = get_init(p1, r1, i, s, gps, prm);

        % --- Perform the simulation
        M1 = make_model2(p1, r1, i, s, gps, prm);
        geq = @(t,in) goveqs_basis3(t, in, M1, i, s, p1, r1, agg, sel, prm);
        [t,soln1] = ode15s(geq, [0:1:tf], r1.init, opts);
        
 
        % --- Record symptomatic incidence (total)
        inc(:,ii,ip) = sum(diff(soln1(:,i.aux.inc),1),2); % Daily incidence
        cinc(:,ii,ip) = sum(soln1(:,i.aux.inc),2);        % Cumulative incidence 
        % --- Record cumulative deaths (total)
        mor(:,ii,ip) = sum(soln1(:,i.aux.mort),2); 
                             
         % --- Record cumulative incidence (group-wise)
        inc_gr1(:,ii,ip) = sum(soln1(:,i.aux.inc(1)),2);
        inc_gr2(:,ii,ip) = sum(soln1(:,i.aux.inc(2)),2);
        inc_gr3(:,ii,ip) = sum(soln1(:,i.aux.inc(3)),2);
        
        % --- Record cumulative deaths (group-wise)
        mor_gr1(:,ii,ip) = sum(soln1(:,i.aux.mort(1)),2); 
        mor_gr2(:,ii,ip) = sum(soln1(:,i.aux.mort(2)),2); 
        mor_gr3(:,ii,ip) = sum(soln1(:,i.aux.mort(3)),2); 
              
        % --- Effective reproduction number
        Reff(:,ii,ip) = R0*sum(soln1(:,1:12),2)./sum(soln1(:,1:72),2);    
        
        
    end
    fprintf('\n');
    
end

% Save the simulation outcome as below for each of the three scenarios (Assigne proper value for line 32).
% Then keep all saved files in the the Fig3 folder to get figure 3
% and the estimates for Table 1 (BMJ Open).

% Susceptibility reducing vaccine  (p.c1 = 0.6 or 0.9; p.c3 = 0)
save simulate1_sus;    % scenario 1: HCW+FW; 
% save simulate2_sus;  % scenario 2: HCW+FW+Co-morbid 
% save simulate3_sus;  % scenario 3: HCW+FW+Co-morbid+Elderly 

% Disease severity reducing vaccine (p.c1 = 0; p.c3 = 0.6 or 0.9)
% save simulate1_sev;  % scenario 1: HCW+FW 
% save simulate2_sev;  % scenario 2: HCW+FW+Co-morbid 
% save simulate3_sev;  % scenario 3: HCW+FW+Co-morbid+Elderly 





