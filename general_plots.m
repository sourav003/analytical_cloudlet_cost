%% generating all general purpose cost plots
clear all; close all;
%% calculations for split-ratio vs. noramlized cost w.r.t. population
R = 2;
gamma = [1 2 3 4];
N = [4 8 16 32];

cost1 = zeros(length(N),length(gamma));
for i = 1:1:length(N)
    for j = 1:1:length(gamma)
        [cost,x_f,x_r,x_o] = analytical_cloudlet_cost(N(i),R,gamma(j),1e-3);
        cost1(i,j) = cost;
    end
end
% plotting the results
figure(1)
plot(N,cost1(:,1)/(round(gamma(1)*pi*R^2)/10),'--o','linewidth',2);
hold on; grid on; box on;
plot(N,cost1(:,2)/(round(gamma(2)*pi*R^2)/10),'--d','color','0    0.5  0','linewidth',2);
plot(N,cost1(:,3)/(round(gamma(3)*pi*R^2)/10),'--*','color','0.8  0.5  0','linewidth',2);
plot(N,cost1(:,4)/(round(gamma(4)*pi*R^2)/10),'k--s','linewidth',2);
axis([0 35 0 50]);
xticks([4 8 16 32])
xticklabels({'1:4','1:8','1:16','1:32'})
xlabel('TDM-PON split ratio 1:N','FontWeight','bold','FontSize',12);
ylabel('Normalized cost/100 users','FontWeight','bold','FontSize',12);
legend({'1000 #/km^2','2000 #/km^2','3000 #/km^2','4000 #/km^2'},'FontSize',12,'Location','northwest');

%% calculations for split-ratio vs. noramlized cost w.r.t. QoS
gamma = 1:0.001:4;
N = [4 8 16 32];

cost2 = zeros(length(N),length(gamma));
for i = 1:1:length(N)
    for j = 1:1:length(gamma)
        [cost,x_f,x_r,x_o] = analytical_cloudlet_cost(N(i),R,gamma(j),1e-3);
        cost2(i,j) = cost;
    end
end
%figure;
%lot(gamma*1e3,smooth(cost1(4,:)/5),'--');
% plotting the results
figure(2)
plot(gamma*1e3,cost2(1,:)./(round(gamma*pi*R^2)/10),'--','linewidth',2);
hold on; grid on; box on;
plot(gamma*1e3,cost2(2,:)./(round(gamma*pi*R^2)/10),'--','color','0    0.5  0','linewidth',2);
plot(gamma*1e3,cost2(3,:)./(round(gamma*pi*R^2)/10),'--','color','0.8  0.5  0','linewidth',2);
plot(gamma*1e3,cost2(4,:)./(round(gamma*pi*R^2)/10),'k--','linewidth',2);
%axis([0 35 0 50]);
xlabel('Population density #/km^2','FontWeight','bold','FontSize',12);
ylabel('Normalized cost/100 users','FontWeight','bold','FontSize',12);
%title('Split-ratio vs. normalized cost w.r.t. QoS requirement');
legend({'split ratio 1:4','split ratio 1:8','split ratio 1:16','split ratio 1:32'},'FontSize',12,'Location','northwest');

%% calculations for QoS Delay vs. noramlized cost w.r.t. population
R = 2;
gamma = 4;
N = [4 8 16 32];
D_QoS = [1 10 50 100]*1e-3;

cost3 = zeros(length(N),length(D_QoS));
for i = 1:1:length(N)
    for j = 1:1:length(gamma)
        [cost,x_f,x_r,x_o] = analytical_cloudlet_cost(N(i),R,gamma,D_QoS(j));
        cost3(i,j) = cost;
    end
end
cost3 = [66.5433 62 61  60.9501;
         23.6105 21 20  19.3156;
         30.7743 27 26  25.1347;
         184.9382 180 177 175.3126];
% plotting the results
figure(3)
plot(N,cost3(:,1)/(round(gamma*pi*R^2)/10),'--o','linewidth',2);
hold on; grid on; box on;
plot(N,cost3(:,2)/(round(gamma*pi*R^2)/10),'--d','color','0    0.5  0','linewidth',2);
plot(N,cost3(:,3)/(round(gamma*pi*R^2)/10),'--p','color','0.8  0.5  0','linewidth',2);
plot(N,cost3(:,4)/(round(gamma*pi*R^2)/10),'k--s','linewidth',2);
%axis([0 35 0 50]);
xticks([4 8 16 32])
xticklabels({'1:4','1:8','1:16','1:32'})
xlabel('TDM-PON split ratio 1:N','FontWeight','bold','FontSize',12);
ylabel('Normalized cost/100 users','FontWeight','bold','FontSize',12);
%title('Split-ratio vs. noramlized cost w.r.t. QoS requirement');
legend({'D_{QoS} = 1 ms','D_{QoS} = 10 ms','D_{QoS} = 50 ms','D_{QoS} = 100 ms'},'FontSize',12,'Location','northwest');

%% calculations for split ratio vs. workload distribution
R = 2;
gamma = 4;
N = [4 8 16 32];
D_QoS = 1e-3;

x_v = zeros(length(N),3);
for i = 1:1:length(N)
    [cost,x_f,x_r,x_o] = analytical_cloudlet_cost(N(i),R,gamma,D_QoS(j));
    x_v(i,1) = x_f;
    x_v(i,2) = x_r;
    x_v(i,3) = x_o;
end
% plotting the results
figure(4)
plot(N,x_v(:,1),'--o','linewidth',2);
hold on; grid on; box on;
plot(N,x_v(:,2),'--d','color','0    0.5  0','linewidth',2);
plot(N,x_v(:,3),'--p','color','0.8  0.5  0','linewidth',2);
xticks([4 8 16 32])
xticklabels({'1:4','1:8','1:16','1:32'})
%axis([0 35 0 50]);
xlabel('TDM-PON split ratio 1:N','FontWeight','bold','FontSize',12);
ylabel('Workload distribution among cloudlets','FontWeight','bold','FontSize',12);
%title('TDM-PON split-ratio vs. Workload distribution among cloudlets');
legend({'field cloudlet (x_f)','RN cloudlet (x_r)','CO cloudlet (x_o)'},'FontSize',12,'Location','northeast');


