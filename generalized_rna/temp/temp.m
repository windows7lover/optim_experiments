% 

clear all

offline = dlmread('offline.log');
online = dlmread('online.log');
% online2 = dlmread('online2.log');

vanilla = offline(:,1);
offline = offline(:,2);
online2 = online(:,2);
online = online(:,1);


figure
hold on
plot(vanilla,'linewidth',1.5)
plot(offline,'linewidth',1.5)
plot(online,'linewidth',1.5)
plot(online2,'linewidth',1.5)
legend({'SGD+momentum','RNA (offline)','RNA (online, before acc)','RNA (online, after acc)'})

display([ min(vanilla), min(offline), min(online), min(online2)])

axis([-inf inf 0 25]);