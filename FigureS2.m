% mat = DA0choice(:,2:end); % DA: Uniform strategy
% mat = DA1choice(:,2:end); % DA: Priority strategy relative to non-vaccinated scenario
% mat = DAchoice(:,2:end);  % DA: Priority strategy relative to uniform strategy

cols = linspecer(3); 
load Out_1p25.mat
cov1 = round(DA1choice_CE(:,1),1); % Coverage (comorbidity then, >60y) 
mat1 = DA1choice_CE(:,2:end); % DA: Priority strategy relative to non-vaccinated scenario
mat01 = DA0choice_CE(:,2:end); % DA: Uniform strategy

cov2 = round(DA1choice_EC(:,1),1); % Coverage (>60y then, comorbidity then) 
mat2 = DA1choice_EC(:,2:end); % DA: Priority strategy relative to non-vaccinated scenario
mat02 = DA0choice_EC(:,2:end); % DA: Uniform strategy

load Out_2.mat
mat3 = DA1choice_CE(:,2:end); % DA: Priority strategy relative to non-vaccinated scenario
mat03 = DA0choice_CE(:,2:end); % DA: Uniform strategy

mat4 = DA1choice_EC(:,2:end); % DA: Priority strategy relative to non-vaccinated scenario
mat04 = DA0choice_EC(:,2:end); % DA: Uniform strategy

load Out_2p5.mat
mat5 = DA1choice_CE(:,2:end); % DA: Priority strategy relative to non-vaccinated scenario
mat05 = DA0choice_CE(:,2:end); % DA: Uniform strategy

mat6 = DA1choice_EC(:,2:end); % DA: Priority strategy relative to non-vaccinated scenario
mat06 = DA0choice_EC(:,2:end); % DA: Uniform strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; fs = 10; lw = 1.5;
yl0 = 3200;
xl0 = 18;

subplot(3,2,1)
y0 = permute(mat01,[2,1,3]); 
% Median plot
plot(cov1(1:2),y0(2,1:2),':k','linewidth', lw); hold on %'Color',[0.5,0.5,0.5]
plot(cov1(2:3),y0(2,2:3),':k','linewidth', lw); hold on
plot(cov1(3:4),y0(2,3:4),':k','linewidth', lw); hold on
% Uncertainty intervals
% jbfill(cov1(1:2)', y0(3,1:2), y0(1,1:2),cols(1,:), 'None', 1, 0.3); hold on
% jbfill(cov1(2:3)', y0(3,2:3), y0(1,2:3),[0.5,0.5,0.5], 'None', 1, 0.3); hold on
% jbfill(cov1(3:4)', y0(3,3:4), y0(1,3:4),[0.5,0.5,0.5], 'None', 1, 0.3); hold on
y = permute(mat1,[2,1,3]); 
% Median plot
h1=plot(cov1(1:2),y(2,1:2),'Color',cols(1,:),'linewidth', lw); hold on
h2=plot(cov1(2:3),y(2,2:3),'Color',cols(2,:),'linewidth', lw); hold on
h3=plot(cov1(3:4),y(2,3:4),'Color',cols(3,:),'linewidth', lw); hold on
% Uncertainty intervals
jbfill(cov1(1:2)', y(3,1:2), y(1,1:2),cols(1,:), 'None', 1, 0.3); hold on
jbfill(cov1(2:3)', y(3,2:3), y(1,2:3),cols(2,:), 'None', 1, 0.3); hold on
jbfill(cov1(3:4)', y(3,3:4), y(1,3:4),cols(3,:), 'None', 1, 0.3); hold on
%xlabel('Coverage (percent of total population)');
%legend([h1 h2 h3],'Key workers','Co-morbidity (<60 year)','>60 year')
ylabel({'Deaths averted','(per million population)'});
set(gca,'fontsize',fs);
set(gca,'XTickLabelRotation',45)
title('R0 = 1.25')
xlim([0 xl0]);
ylim([0 yl0]);    

subplot(3,2,2)
y0 = permute(mat02,[2,1,3]); 
% Median plot
plot(cov2(1:2),y0(2,1:2),':k','linewidth', lw); hold on
plot(cov2(2:3),y0(2,2:3),':k','linewidth', lw); hold on
plot(cov2(3:4),y0(2,3:4),':k','linewidth', lw); hold on

y = permute(mat2,[2,1,3]); 
% Median plot
h1=plot(cov2(1:2),y(2,1:2),'Color',cols(1,:),'linewidth', lw); hold on
h2=plot(cov2(2:3),y(2,2:3),'Color',cols(3,:),'linewidth', lw); hold on
h3=plot(cov2(3:4),y(2,3:4),'Color',cols(2,:),'linewidth', lw); hold on
% Uncertainty intervals
jbfill(cov2(1:2)', y(3,1:2), y(1,1:2),cols(1,:), 'None', 1, 0.3); hold on
jbfill(cov2(2:3)', y(3,2:3), y(1,2:3),cols(3,:), 'None', 1, 0.3); hold on
jbfill(cov2(3:4)', y(3,3:4), y(1,3:4),cols(2,:), 'None', 1, 0.3); hold on
% xlabel('Coverage (percent of total population)');
% ylabel({'Deaths averted','(million)'});
legend([h1 h2 h3],'Key workers','>60 year','Co-morbidity (<60 year)')
set(gca,'fontsize',fs);
%set(gca,'xtick',[0:9],'xticklabel',cov)
set(gca,'XTickLabelRotation',45)
% line([3 3],[0 3.5], 'linestyle', ':','Color',0.2*[1 1 1],'linewidth', lw)
% line([6 6],[0 3.5], 'linestyle', ':','Color',0.2*[1 1 1],'linewidth', lw)
title('R0 = 1.25')
xlim([0 xl0]);
ylim([0 yl0]);    
%ylim([-0.15  0]);    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,2,3)
y0 = permute(mat03,[2,1,3]); 
% Median plot
plot(cov1(1:2),y0(2,1:2),':k','linewidth', lw); hold on
plot(cov1(2:3),y0(2,2:3),':k','linewidth', lw); hold on
plot(cov1(3:4),y0(2,3:4),':k','linewidth', lw); hold on
% Uncertainty intervals
% jbfill(cov1(1:2)', y0(3,1:2), y0(1,1:2),cols(1,:), 'None', 1, 0.3); hold on
% jbfill(cov1(2:3)', y0(3,2:3), y0(1,2:3),[0.5,0.5,0.5], 'None', 1, 0.3); hold on
% jbfill(cov1(3:4)', y0(3,3:4), y0(1,3:4),[0.5,0.5,0.5], 'None', 1, 0.3); hold on
y = permute(mat3,[2,1,3]); 
% Median plot
h1=plot(cov1(1:2),y(2,1:2),'Color',cols(1,:),'linewidth', lw); hold on
h2=plot(cov1(2:3),y(2,2:3),'Color',cols(2,:),'linewidth', lw); hold on % 'Color',cols(2,:)
h3=plot(cov1(3:4),y(2,3:4),'Color',cols(3,:),'linewidth', lw); hold on % 'Color',cols(3,:)
% Uncertainty intervals
jbfill(cov1(1:2)', y(3,1:2), y(1,1:2),cols(1,:), 'None', 1, 0.3); hold on
jbfill(cov1(2:3)', y(3,2:3), y(1,2:3),cols(2,:), 'None', 1, 0.3); hold on
jbfill(cov1(3:4)', y(3,3:4), y(1,3:4),cols(3,:), 'None', 1, 0.3); hold on
%xlabel('Coverage (percent of total population)');
%legend([h1 h2 h3],'Key workers','Co-morbidity (<60 year)','>60 year')
ylabel({'Deaths averted','(per million population)'});
set(gca,'fontsize',fs);
set(gca,'XTickLabelRotation',45)
title('R0 = 2')
xlim([0 xl0]);
ylim([0 yl0]);    

subplot(3,2,4)
y0 = permute(mat04,[2,1,3]); 
% Median plot
plot(cov2(1:2),y0(2,1:2),':k','linewidth', lw); hold on
plot(cov2(2:3),y0(2,2:3),':k','linewidth', lw); hold on
plot(cov2(3:4),y0(2,3:4),':k','linewidth', lw); hold on

y = permute(mat4,[2,1,3]); 
% Median plot
h1=plot(cov2(1:2),y(2,1:2),'Color',cols(1,:),'linewidth', lw); hold on
h2=plot(cov2(2:3),y(2,2:3),'Color',cols(3,:),'linewidth', lw); hold on
h3=plot(cov2(3:4),y(2,3:4),'Color',cols(2,:),'linewidth', lw); hold on
% Uncertainty intervals
jbfill(cov2(1:2)', y(3,1:2), y(1,1:2),cols(1,:), 'None', 1, 0.3); hold on
jbfill(cov2(2:3)', y(3,2:3), y(1,2:3),cols(3,:), 'None', 1, 0.3); hold on
jbfill(cov2(3:4)', y(3,3:4), y(1,3:4),cols(2,:), 'None', 1, 0.3); hold on
% xlabel('Coverage (percent of total population)');
% ylabel({'Deaths averted','(million)'});
%legend([h1 h2 h3],'Key workers','>60 year','Co-morbidity (<60 year)')
set(gca,'fontsize',fs);
%set(gca,'xtick',[0:9],'xticklabel',cov)
set(gca,'XTickLabelRotation',45)
% line([3 3],[0 3.5], 'linestyle', ':','Color',0.2*[1 1 1],'linewidth', lw)
% line([6 6],[0 3.5], 'linestyle', ':','Color',0.2*[1 1 1],'linewidth', lw)
title('R0 = 2')
xlim([0 xl0]);
ylim([0 yl0]);    
%ylim([-0.15  0]);    
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(3,2,5)
y0 = permute(mat05,[2,1,3]); 
% Median plot
plot(cov1(1:2),y0(2,1:2),':k','linewidth', lw); hold on
plot(cov1(2:3),y0(2,2:3),':k','linewidth', lw); hold on
plot(cov1(3:4),y0(2,3:4),':k','linewidth', lw); hold on
% Uncertainty intervals
% jbfill(cov1(1:2)', y0(3,1:2), y0(1,1:2),cols(1,:), 'None', 1, 0.3); hold on
% jbfill(cov1(2:3)', y0(3,2:3), y0(1,2:3),[0.5,0.5,0.5], 'None', 1, 0.3); hold on
% jbfill(cov1(3:4)', y0(3,3:4), y0(1,3:4),[0.5,0.5,0.5], 'None', 1, 0.3); hold on
y = permute(mat5,[2,1,3]); 
% Median plot
h1=plot(cov1(1:2),y(2,1:2),'Color',cols(1,:),'linewidth', lw); hold on
h2=plot(cov1(2:3),y(2,2:3),'Color',cols(2,:),'linewidth', lw); hold on % 'Color',cols(2,:)
h3=plot(cov1(3:4),y(2,3:4),'Color',cols(3,:),'linewidth', lw); hold on % 'Color',cols(3,:)
% Uncertainty intervals
jbfill(cov1(1:2)', y(3,1:2), y(1,1:2),cols(1,:), 'None', 1, 0.3); hold on
jbfill(cov1(2:3)', y(3,2:3), y(1,2:3),cols(2,:), 'None', 1, 0.3); hold on
jbfill(cov1(3:4)', y(3,3:4), y(1,3:4),cols(3,:), 'None', 1, 0.3); hold on
xlabel('Coverage (percent of total population)');
%legend([h1 h2 h3],'Key workers','Co-morbidity (<60 year)','>60 year')
ylabel({'Deaths averted','(per million population)'});
set(gca,'fontsize',fs);
set(gca,'XTickLabelRotation',45)
title('R0 = 2.5')
xlim([0 xl0]);
ylim([0 yl0]);    

subplot(3,2,6)
y0 = permute(mat06,[2,1,3]); 
% Median plot
plot(cov2(1:2),y0(2,1:2),':k','linewidth', lw); hold on
plot(cov2(2:3),y0(2,2:3),':k','linewidth', lw); hold on
plot(cov2(3:4),y0(2,3:4),':k','linewidth', lw); hold on

y = permute(mat6,[2,1,3]); 
% Median plot
h1=plot(cov2(1:2),y(2,1:2),'Color',cols(1,:),'linewidth', lw); hold on
h2=plot(cov2(2:3),y(2,2:3),'Color',cols(3,:),'linewidth', lw); hold on
h3=plot(cov2(3:4),y(2,3:4),'Color',cols(2,:),'linewidth', lw); hold on
% Uncertainty intervals
jbfill(cov2(1:2)', y(3,1:2), y(1,1:2),cols(1,:), 'None', 1, 0.3); hold on
jbfill(cov2(2:3)', y(3,2:3), y(1,2:3),cols(3,:), 'None', 1, 0.3); hold on
jbfill(cov2(3:4)', y(3,3:4), y(1,3:4),cols(2,:), 'None', 1, 0.3); hold on
xlabel('Coverage (percent of total population)');
% ylabel({'Deaths averted','(million)'});
%legend([h1 h2 h3],'Key workers','>60 year','Co-morbidity (<60 year)')
set(gca,'fontsize',fs);
%set(gca,'xtick',[0:9],'xticklabel',cov)
set(gca,'XTickLabelRotation',45)
% line([3 3],[0 3.5], 'linestyle', ':','Color',0.2*[1 1 1],'linewidth', lw)
% line([6 6],[0 3.5], 'linestyle', ':','Color',0.2*[1 1 1],'linewidth', lw)
title('R0 = 2.5')
xlim([0 xl0]);
ylim([0 yl0]);    
%ylim([-0.15  0]);    
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 


sgtitle({'Priority based strategy relative to no-vaccination','(Symptomatic disease-preventing vaccine of efficacy 60%)'}) 
