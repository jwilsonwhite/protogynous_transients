function Fig2_open_pop_trajectories

% Plot population trajectories, as a function of F and phi, for open
% population model

% read data
load transient_open_may2019.mat
figure(1)
clf
set(gcf,'units','cent','position',[10 10 17 18])

Years = 0:1:20 ;
Fs = [2 6 21]; % F = 0.05, 0.5, 1
F(Fs)
Phis = [2, 18];

Cols1 = flipud([0.8 0.2 0.2; 0.2 0.2 0.8; 0 0 0]);
Cols2 = winter(3);
Cols2(2,:) = [0 0 0]; 
Cols2(:,1) = [1 0 0];

% subplot 1-2: Sex ratio
SPs = [1 3];

for ph = 1:2
subplot(4,2,SPs(ph))
hold on

for f = 1:length(Fs)

  plot(Years,ones(length(Years),1)*0.5,'k-','linewidth',1.5)
  pF(f) = plot(Years,Transient.SC1.F(Fs(f)).Lf(2).PHI(Phis(ph)).NumSRatio(51:71));
  pFF(f) = plot(Years,Transient.SC2.F(Fs(f)).Lf(2).PHI(Phis(ph)).NumSRatio(51:71));

set(pF(f),'color',Cols1(f,:),'linewidth',1.5,'linestyle','--') ; % Red
set(pFF(f),'color',Cols1(f,:),'linewidth',1.5,'linestyle',':') ; % Red
end % end loop over Fs

if ph == 2
xlabel(gca,'Time (years after MPA implementation)','fontsize',14)
ylabel(gca,'Sex ratio','fontsize',14)
end

set(gca,'xlim',[0 15])
set(gca,'ylim',[0 0.55],'ytick',0:0.1:1)
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])
end



%subplot 3-4: Egg production
SPs = [5 7];

for ph = 1:2
subplot(4,2,SPs(ph))
hold on

for f = 1:length(Fs)

  pG(f) = plot(Years,Transient.GON.F(Fs(f)).Lf(2).PHI(Phis(ph)).FertEggs(51:71)./Transient.GON.F(Fs(f)).Lf(2).PHI(phi).FertEggs(end));

  pF(f) = plot(Years,Transient.SC1.F(Fs(f)).Lf(2).PHI(Phis(ph)).FertEggs(51:71)./Transient.SC1.F(Fs(f)).Lf(2).PHI(phi).FertEggs(end));
  pFF(f) = plot(Years,Transient.SC2.F(Fs(f)).Lf(2).PHI(Phis(ph)).FertEggs(51:71)./Transient.SC2.F(Fs(f)).Lf(2).PHI(phi).FertEggs(end));

set(pG(f),'color',Cols1(f,:),'linewidth',1.5,'linestyle','-') ; % Red
set(pF(f),'color',Cols1(f,:),'linewidth',1.5,'linestyle','--') ; % Red
set(pFF(f),'color',Cols1(f,:),'linewidth',1.5,'linestyle',':') ; % Red
end % end loop over Fs

set(gca,'ylim',[0 1.05],'ytick',0.0:0.25:1)
if ph == 2
xlabel(gca,'Time (years after MPA implementation)','fontsize',14)
ylabel(gca,'Egg production','fontsize',14)
end
set(gca,'xlim',[0 15])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])

end


%--------------------------------------------
% Plots on right: Duration vs. F

subplot(4,2,[2 4])
hold on

phi = 10;

for l = 1:3
for f = 1:length(Transient.SC1.F)
  F(f) = Transient.SC1.F(f).F;
  TC0(f) = Transient.GON.F(f).Lf(l).PHI(phi).TimeConv;
  TC1(f) = Transient.SC1.F(f).Lf(l).PHI(phi).TimeConv;
  TC2(f) = Transient.SC2.F(f).Lf(l).PHI(phi).TimeConv;
  
end

plot(F,TC0,'color',Cols2(l,:),'linewidth',1.5,'linestyle','-')
plot(F,TC1,'color',Cols2(l,:),'linewidth',1.5,'linestyle','--')
plot(F,TC2,'color',Cols2(l,:),'linewidth',1.5,'linestyle',':')
end

xlabel(gca,'Fishing rate (y-1)','fontsize',14)
ylabel(gca,'Transient duration (y)','fontsize',14)
ylim([0 20])
xlim([0 1])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])



subplot(4,2,[6 8])
hold on

ff = 8;
clear TC0 TC1 TC2
for l = 1:3
for f = 1:length(Transient.SC1.F(1).Lf(1).PHI)
  Phi(f) = Transient.SC1.F(ff).Lf(l).PHI(f).PHI;
  TC0(f) = Transient.GON.F(ff).Lf(l).PHI(f).TimeConv;
  TC1(f) = Transient.SC1.F(ff).Lf(l).PHI(f).TimeConv;
  TC2(f) = Transient.SC2.F(ff).Lf(l).PHI(f).TimeConv;
  
end

plot(Phi,TC0,'color',Cols2(l,:),'linewidth',1.5,'linestyle','-')
plot(Phi,TC1,'color',Cols2(l,:),'linewidth',1.5,'linestyle','--')
plot(Phi,TC2,'color',Cols2(l,:),'linewidth',1.5,'linestyle',':')
end % end l

xlabel(gca,'Male importance (phi)','fontsize',14)
ylabel(gca,'Transient duration (y)','fontsize',14)
ylim([0 20])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])


