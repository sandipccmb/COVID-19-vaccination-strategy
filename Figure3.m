clear all; 

% Use function 'linespecer':It basically provides colour selections that are suitable for publication.
cols = linspecer(9); 

figure; fs = 14; lw = 1.5;transparency1=0.1;transparency2=0.1; 
xl0= 350;
yl0 = 11000;
yl1 = 12000;
subplot(2,2,1)

load simulate1_sus; %(Key workers)
% Find percentiles for incidence
inc_prc = prctile(inc,[2.5,50,97.5],2);                
y = permute(inc_prc,[2,1,3]);
h1=plot(y(2,:,1),'Color',cols(2,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:,1), y(1,:,1),cols(2,:), 'None', 1, 0.3); hold on;

ylabel({'Daily symptomatic cases','(numbers/million)'}); %(Key workers)
h2=plot(y(2,:,2),'Color',cols(1,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:,2), y(1,:,2),cols(1,:), 'None', 1, 0.3); hold on;

load simulate2_sus; % (Co-morbidity) 
% Find percentiles for incidence
inc_prc = prctile(inc,[2.5,50,97.5],2);                % Incidence per day 
y = permute(inc_prc,[2,1,3]);
h3=plot(y(2,:,2),'Color',cols(4,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:,2), y(1,:,2),cols(4,:), 'None', 1, 0.3); hold on;

load simulate3_sus;  % (Elderly)
% Find percentiles for incidence
inc_prc = prctile(inc,[2.5,50,97.5],2);                % Incidence per day 
y = permute(inc_prc,[2,1,3]);
h4 =plot(y(2,:,2),'Color',cols(3,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:,2), y(1,:,2),cols(3,:), 'None', 1, 0.3); 
%legend([h1 h2 h3 h4],'No vaccination','Scenario 1','Scenario 2','Scenario 3')
%legend([h1 h4],'No vaccine','With vaccine')
xlim([0 xl0]);
ylim([0 yl0]);




subplot(2,2,2)
load simulate1_sus;
% Find percentiles for deaths
mor_prc = prctile(mor(1:end-1,:,:),[2.5,50,97.5],2); % to adjust the length of the matrix I remove last row
y = permute(mor_prc,[2,1,3]);
h1 = plot(y(2,:,1),'Color',cols(2,:),'linewidth', lw); %cols{ii}
jbfill(1:tf, y(3,:,1), y(1,:,1),cols(2,:), 'None', 1, 0.3);hold on; 
ylabel({'Cumulative deaths','(numbers/million)'});

h2 = plot(y(2,:,2),'Color',cols(1,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:,2), y(1,:,2),cols(1,:), 'None', 1, 0.3); hold on;

load simulate2_sus;
mor_prc = prctile(mor(1:end-1,:,:),[2.5,50,97.5],2); % to adjust the length of the matrix I remove last row
y = permute(mor_prc,[2,1,3]);
h3= plot(y(2,:,2),'Color',cols(4,:),'linewidth', lw);  hold on;%cols{ii}
jbfill(1:tf, y(3,:,2), y(1,:,2),cols(4,:), 'None', 1, 0.3);  hold on;

load simulate3_sus;
mor_prc = prctile(mor(1:end-1,:,:),[2.5,50,97.5],2); % to adjust the length of the matrix I remove last row
y = permute(mor_prc,[2,1,3]);
h4 =plot(y(2,:,2),'Color',cols(3,:),'linewidth', lw);  hold on;
jbfill(1:tf, y(3,:,2), y(1,:,2),cols(3,:), 'None', 1, 0.3); 
legend([h1 h2 h3 h4],'No vaccination','KeyW','KeyW + Co-M','KeyW + Co-M + >60')
%legend([h1 h4],'No vaccine','With vaccine')
xlim([0 xl0]);
ylim([0 yl1]);

%suptitle('Infection preventing vaccine')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3)

load simulate1_sev; %(Key workers)
% Find percentiles for incidence
inc_prc = prctile(inc,[2.5,50,97.5],2);                
y = permute(inc_prc,[2,1,3]);
h1=plot(y(2,:,1),'Color',cols(2,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:,1), y(1,:,1),cols(2,:), 'None', 1, 0.3); hold on;

ylabel({'Daily symptomatic cases','(numbers/million)'}); %(Key workers)
h2=plot(y(2,:,2),'Color',cols(1,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:,2), y(1,:,2),cols(1,:), 'None', 1, 0.3); hold on;

load simulate2_sev; % (Co-morbidity) 
% Find percentiles for incidence
inc_prc = prctile(inc,[2.5,50,97.5],2);                % Incidence per day 
y = permute(inc_prc,[2,1,3]);
h3=plot(y(2,:,2),'Color',cols(4,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:,2), y(1,:,2),cols(4,:), 'None', 1, 0.3); hold on;

load simulate3_sev;  % (Elderly)
% Find percentiles for incidence
inc_prc = prctile(inc,[2.5,50,97.5],2);                % Incidence per day 
y = permute(inc_prc,[2,1,3]);
h4 =plot(y(2,:,2),'Color',cols(3,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:,2), y(1,:,2),cols(3,:), 'None', 1, 0.3); 
%legend([h1 h2 h3 h4],'No vaccination','Scenario 1','Scenario 2','Scenario 3')
%legend([h1 h4],'No vaccine','With vaccine')
xlim([0 xl0]);
ylim([0 yl0]);




subplot(2,2,4)
load simulate1_sev;
% Find percentiles for deaths
mor_prc = prctile(mor(1:end-1,:,:),[2.5,50,97.5],2); % to adjust the length of the matrix I remove last row
y = permute(mor_prc,[2,1,3]);
h1 = plot(y(2,:,1),'Color',cols(2,:),'linewidth', lw); %cols{ii}
jbfill(1:tf, y(3,:,1), y(1,:,1),cols(2,:), 'None', 1, 0.3);hold on; 
ylabel({'Cumulative deaths','(numbers/million)'});

h2 = plot(y(2,:,2),'Color',cols(1,:),'linewidth', lw); hold on;
jbfill(1:tf, y(3,:,2), y(1,:,2),cols(1,:), 'None', 1, 0.3); hold on;

load simulate2_sev;
mor_prc = prctile(mor(1:end-1,:,:),[2.5,50,97.5],2); % to adjust the length of the matrix I remove last row
y = permute(mor_prc,[2,1,3]);
h3= plot(y(2,:,2),'Color',cols(4,:),'linewidth', lw);  hold on;%cols{ii}
jbfill(1:tf, y(3,:,2), y(1,:,2),cols(4,:), 'None', 1, 0.3);  hold on;

load simulate3_sev;
mor_prc = prctile(mor(1:end-1,:,:),[2.5,50,97.5],2); % to adjust the length of the matrix I remove last row
y = permute(mor_prc,[2,1,3]);
h4 =plot(y(2,:,2),'Color',cols(3,:),'linewidth', lw);  hold on;
jbfill(1:tf, y(3,:,2), y(1,:,2),cols(3,:), 'None', 1, 0.3); hold on;
%legend([h1 h2 h3 h4],'No vaccination','Scenario 1','Scenario 2','Scenario 3')
%legend([h1 h4],'No vaccine','With vaccine')
xlim([0 xl0]);
ylim([0 yl1]);

%suptitle('Disease preventing vaccine')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % IMPORTANT: See Table1.mat for peak incidence reduction, cumulative mortality and
% % deaths averted per vaccinated 






