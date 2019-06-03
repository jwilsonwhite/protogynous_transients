function Fig3_closed_pop_trajectories

% Plot population trajectories and other statistics, as a function of F and
% phi, for closed population model

% read data
load('transient_feb2019.mat')


% First figure: plot trajectories of N, sex ratio, egg production
figure(1)
clf
set(gcf,'units','cent','position',[10 10 18 21])

T = 21;
Years = 0:1:20 ;
Fs = [1 2 11]; % F = 0.05, 0.5, 1
Phis = [2, 18];

Cols1 = flipud([0.8 0.2 0.2; 0.2 0.2 0.8; 0 0 0]);
Cols2 = winter(3);
Cols2(2,:) = [0 0 0]; 

% subplot 1: Abundance
SPs = [1 2];

for ph = 1:2
subplot(3,2,SPs(ph))
hold on

for f = 1:length(Fs)

  pG(f) = plot(Years,Transient.GON.F(Fs(f)).Lf(2).PHI(Phis(ph)).Ntotal(51:(50+T))/Transient.GON.F(Fs(f)).Lf(2).PHI(Phis(ph)).Ntotal(51));

  pF(f) = plot(Years,Transient.SC1.F(Fs(f)).Lf(2).PHI(Phis(ph)).Ntotal(51:(50+T))/Transient.SC1.F(Fs(f)).Lf(2).PHI(Phis(ph)).Ntotal(51));
  pFF(f) = plot(Years,Transient.SC2.F(Fs(f)).Lf(2).PHI(Phis(ph)).Ntotal(51:(50+T))/Transient.SC2.F(Fs(f)).Lf(2).PHI(Phis(ph)).Ntotal(51));

set(pG(f),'color',Cols1(f,:),'linewidth',1.5,'linestyle','-') ; % Red
set(pF(f),'color',Cols1(f,:),'linewidth',1.5,'linestyle','--') ; % Red
set(pFF(f),'color',Cols1(f,:),'linewidth',1.5,'linestyle',':') ; % Red
end % end loop over Fs

if ph == 1
xlabel(gca,'Time (years after MPA implementation)','fontsize',14)
ylabel(gca,'Abundance ratio','fontsize',14)
end

set(gca,'xlim',[0 T-1])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])
end

% subplot 1: Sex ratio
SPs = [3 4];

for ph = 1:2
subplot(3,2,SPs(ph))
hold on

for f = 1:length(Fs)

  pF(f) = plot(Years,Transient.SC1.F(Fs(f)).Lf(1).PHI(Phis(ph)).NumSRatio(51:(50+T)));
  pFF(f) = plot(Years,Transient.SC2.F(Fs(f)).Lf(1).PHI(Phis(ph)).NumSRatio(51:(50+T)));

set(pF(f),'color',Cols1(f,:),'linewidth',1.5,'linestyle','--') ; % Red
set(pFF(f),'color',Cols1(f,:),'linewidth',1.5,'linestyle',':') ; % Red
end % end loop over Fs

plot([0 T],[0.5 0.5],'k-') % gonochores

if ph == 1
xlabel(gca,'Time (years after MPA implementation)','fontsize',14)
ylabel(gca,'Sex ratio','fontsize',14)
end

set(gca,'xlim',[0 T-1])
set(gca,'ylim',[0 0.6])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])
end



%subplot 3-4: Egg production
SPs = [5 6];

for ph = 1:2
subplot(3,2,SPs(ph))
hold on

for f = 1:length(Fs)

  pG(f) = plot(Years,Transient.GON.F(Fs(f)).Lf(2).PHI(Phis(ph)).FertEggs(51:(50+T))./Transient.GON.F(Fs(f)).Lf(2).PHI(Phis(ph)).FertEggs(51));

  pF(f) = plot(Years,Transient.SC1.F(Fs(f)).Lf(2).PHI(Phis(ph)).FertEggs(51:(50+T))./Transient.SC1.F(Fs(f)).Lf(2).PHI(Phis(ph)).FertEggs(51));
  pFF(f) = plot(Years,Transient.SC2.F(Fs(f)).Lf(2).PHI(Phis(ph)).FertEggs(51:(50+T))./Transient.SC2.F(Fs(f)).Lf(2).PHI(Phis(ph)).FertEggs(51));


set(pG(f),'color',Cols1(f,:),'linewidth',1.5,'linestyle','-') ; % Red
set(pF(f),'color',Cols1(f,:),'linewidth',1.5,'linestyle','--') ; % Red
set(pFF(f),'color',Cols1(f,:),'linewidth',1.5,'linestyle',':') ; % Red
end % end loop over Fs

set(gca,'ylim',[1 2.5],'ytick',0:0.5:3)
if ph == 1
xlabel(gca,'Time (years after MPA implementation)','fontsize',14)
ylabel(gca,'Egg production','fontsize',14)
end
set(gca,'xlim',[0 T-1])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])

end


%--------------------------------------------
% Plots on right: Duration vs. F
figure(2)
clf
set(gcf,'units','cent','position',[20 10 9 18])


subplot(2,1,1)
hold on

phi = 10;

for l = 1:3
for f = 1:length(Transient.SC1.F)
  F(f) = Transient.SC1.F(f).F;
  TC0(f) = Transient.GON.F(f).Lf(l).PHI(phi).TimeConv;
  TC1(f) = Transient.SC1.F(f).Lf(l).PHI(phi).TimeConv;
  TC2(f) = Transient.SC2.F(f).Lf(l).PHI(phi).TimeConv;
  
end
%keyboard
plot(F,TC0,'color',Cols2(l,:),'linewidth',1.5,'linestyle','-')
plot(F,TC1,'color',Cols2(l,:),'linewidth',1.5,'linestyle','--')
plot(F,TC2,'color',Cols2(l,:),'linewidth',1.5,'linestyle',':')
end


xlabel(gca,'Fishing rate (y-1)','fontsize',14)
ylabel(gca,'Transient duration (y)','fontsize',14)
set(gca,'xlim',[0 1])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])



subplot(2,1,2)
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
ylim([7 50])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])

% --------------------------------------------
% Fig 3, initial trajectory
figure(3)
clf
set(gcf,'units','cent','position',[40 10 9 18])


subplot(2,1,1)
hold on

phi = 10;

for l = 1:3
for f = 1:length(Transient.SC1.F)
  F(f) = Transient.SC1.F(f).F;
  TC0(f) = Transient.GON.F(f).Lf(l).PHI(phi).LambInit;
  TC1(f) = Transient.SC1.F(f).Lf(l).PHI(phi).LambInit;
  TC2(f) = Transient.SC2.F(f).Lf(l).PHI(phi).LambInit;
  
end
%keyboard
plot(F,TC0,'color',Cols2(l,:),'linewidth',1.5,'linestyle','-')
plot(F,TC1,'color',Cols2(l,:),'linewidth',1.5,'linestyle','--')
plot(F,TC2,'color',Cols2(l,:),'linewidth',1.5,'linestyle',':')
end

xlabel(gca,'Fishing rate (y-1)','fontsize',14)
ylabel(gca,'Initial trajectory','fontsize',14)
set(gca,'xlim',[0 1])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])



subplot(2,1,2)
hold on

ff = 8;
clear TC0 TC1 TC2
for l = 1:3
for f = 1:length(Transient.SC1.F(1).Lf(1).PHI)
  Phi(f) = Transient.SC1.F(ff).Lf(l).PHI(f).PHI;
  TC0(f) = Transient.GON.F(ff).Lf(l).PHI(f).LambInit;
  TC1(f) = Transient.SC1.F(ff).Lf(l).PHI(f).LambInit;
  TC2(f) = Transient.SC2.F(ff).Lf(l).PHI(f).LambInit;
  
end

plot(Phi,TC0,'color',Cols2(l,:),'linewidth',1.5,'linestyle','-')
plot(Phi,TC1,'color',Cols2(l,:),'linewidth',1.5,'linestyle','--')
plot(Phi,TC2,'color',Cols2(l,:),'linewidth',1.5,'linestyle',':')
end % end l

xlabel(gca,'Male importance (phi)','fontsize',14)
ylabel(gca,'Initial trajectory','fontsize',14)
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])

% --------------------------------------------
% Fig 4, Theta
figure(4)
clf
set(gcf,'units','cent','position',[50 10 9 18])


subplot(2,1,1)
hold on

phi = 10;

for l = 1:3
for f = 1:length(Transient.SC1.F)
  F(f) = Transient.SC1.F(f).F;
  TC0(f) = Transient.GON.F(f).Lf(l).PHI(phi).Theta2;
  TC1(f) = Transient.SC1.F(f).Lf(l).PHI(phi).Theta2;
  TC2(f) = Transient.SC2.F(f).Lf(l).PHI(phi).Theta2;
  
end

plot(F,TC0,'color',Cols2(l,:),'linewidth',1.5,'linestyle','-')
plot(F,TC1,'color',Cols2(l,:),'linewidth',1.5,'linestyle','--')
plot(F,TC2,'color',Cols2(l,:),'linewidth',1.5,'linestyle',':')
end


xlabel(gca,'Fishing rate (y-1)','fontsize',14)
ylabel(gca,'Theta','fontsize',14)
set(gca,'xlim',[0 1])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])



subplot(2,1,2)
hold on

ff = 8;
clear TC0 TC1 TC2
for l = 1:3
for f = 1:length(Transient.SC1.F(1).Lf(1).PHI)
  Phi(f) = Transient.SC1.F(ff).Lf(l).PHI(f).PHI;
  TC0(f) = Transient.GON.F(ff).Lf(l).PHI(f).Theta2;
  TC1(f) = Transient.SC1.F(ff).Lf(l).PHI(f).Theta2;
  TC2(f) = Transient.SC2.F(ff).Lf(l).PHI(f).Theta2;
  
end

plot(Phi,TC0,'color',Cols2(l,:),'linewidth',1.5,'linestyle','-')
plot(Phi,TC1,'color',Cols2(l,:),'linewidth',1.5,'linestyle','--')
plot(Phi,TC2,'color',Cols2(l,:),'linewidth',1.5,'linestyle',':')
end % end l

xlabel(gca,'Male importance (phi)','fontsize',14)
ylabel(gca,'Theta','fontsize',14)
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.02 0.02])

