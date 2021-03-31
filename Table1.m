%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load simulate1_sus; 
%load simulate1_sev;

%load simulate2_sus; 
%load simulate2_sev;

load simulate3_sus; 
%load simulate3_sev;

% Vaccinated total
% vacc_tot = [28769343, 96313092, 341798410];   % [Scenario 1, Scenario 2, Scenario 3]
 vacc_tot = [27235304, 110612567, 247233001];   % [Scenario 1, Scenario 2, Scenario 3]

 % choose approriate 
target_pop =  vacc_tot(3); % 1: simulate1; 2: simulate 2; simulate3: 3
 
% Find percentiles of reduction of peak incidence
for ii = 1:1:250
inc_b = inc(:,ii,1);
inc_i = inc(:,ii,2);
Red_peakinc(ii) = 100*((max(inc_b)- max(inc_i))/ max(inc_b)); % percentage reduction
%Red_peakinc(ii) = (max(inc_b)- max(inc_i)); % Absulate reduction
end
perRed_peakinc=prctile(Red_peakinc,[2.5,50,97.5],2) 

% Find percentiles of reduction of cumulative mortality
for ii = 1:1:250
mor_b = mor(end,ii,1);
mor_i = mor(end,ii,2);
Red_cummor(ii) = 100*((mor_b- mor_i)/ mor_b);
%Red_cummor(ii) = ((mor_b- mor_i));
end
perRed_cummor = prctile(Red_cummor,[2.5,50,97.5],2) 

% Deaths averted per individual vaccinated
 % pop = (prm.N*1380004385/1e6);                             % Actual population in three age group
 DA = (mor(end,:,1) - mor(end,:,2)).*(1380004385/1e6);     % Total DA 

 %  % Find percentiles for DA (Total)  
%  DA1 = prctile(DA/target_pop,[2.5,50,97.5],2)

 % Number needed to vaccinate to avert one death
 DA2 = prctile(target_pop./DA,[2.5,50,97.5],2)
 
 
 
 
 
 
 
 
 
 
%  %%%%%%%%%%%%%%%%%%%%%%%%
% % load DA1_R01p25.mat
% %load DA1_R02.mat
% load DA1_R02p5.mat
% cov = round(DAchoice(:,1),1); % Coverage
% mat = DAchoice(:,2:end);      % Deaths averted
% y = permute(mat,[2,1,3]); 
% scenario1 = y(2,end);
% 
% load DA2_R02p5.mat
% cov = round(DAchoice(:,1),1); % Coverage
% mat = DAchoice(:,2:end);      % Deaths averted
% y = permute(mat,[2,1,3]); 
% scenario2 = y(2,end);
%  
% load DA3_R02p5.mat
% cov = round(DAchoice(:,1),1); % Coverage
% mat = DAchoice(:,2:end);      % Deaths averted
% y = permute(mat,[2,1,3]); 
% scenario3 = y(2,end);
%  
% load DA4_R02p5.mat
% cov = round(DAchoice(:,1),1); % Coverage
% mat = DAchoice(:,2:end);      % Deaths averted
% y = permute(mat,[2,1,3]); 
% scenario4 = y(2,end);
% 
% load DA5_R02p5.mat
% cov = round(DAchoice(:,1),1); % Coverage
% mat = DAchoice(:,2:end);      % Deaths averted
% y = permute(mat,[2,1,3]); 
% scenario5 = y(2,end);
% 
% load DA6_R02p5.mat
% cov = round(DAchoice(:,1),1); % Coverage
% mat = DAchoice(:,2:end);      % Deaths averted
% y = permute(mat,[2,1,3]); 
% scenario6 = y(2,end);
%  
% [scenario1,scenario2,scenario3,scenario4,scenario5,scenario6]
