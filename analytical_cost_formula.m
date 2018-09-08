%% cloudlet deployment cost calculation using analytical formulae
% Author: Sourav Mondal, EEE, The University of Melbourne
% =======================================================
close all;clear all;clc;
% general network parameters
R = 2;                              % radius of circular area (km)
N = 8;                              % split ratio of TDM-PON
pon = ceil(50/N);                   % total number of PON branches
theta = 2*pi/pon;                   % inscribed angle of the circular wedge
phi = theta/2;
varphi = 2500;                      % no. of VMs per cloudlet rack
centroid = 2*R*sin(phi)/(3*phi);    % distance of centroid
cent_avg = (16*R/27)*(sin(phi)/phi)^3*(1-(2*sin(phi)/(3*phi)));     % avg. distance from centroid
c = 3*10^8;                         % velocity of light (m/s)
D_r = cent_avg*10^3/c;
D_f = cent_avg*10^3/c;
D_o = (cent_avg+centroid)*10^3/c;
BW_f = 1*10^9;                      % maximum bandwidth in the up/down-link field
BW_r = 10*10^9;                     % maximum bandwidth in the up/down-link RN
BW_o = 10*10^9;             		% maximum bandwidth in the up/down-link CO
bl_ul = 5*10^9;                     % background load in uplink of TDM-PON
bl_dl = 7*10^9;                     % background load in downlink of TDM-PON
n_lambda = 1;                       % number of wavelengths to RN cloudelet
zeta_f = 4;                         % normalized cost of cloudlet installation
zeta_r = 3;
zeta_o = 2;
eta = 50;                           % cost of fiber installation per km
sigma_ul = 100*10^3;            	    % total computation load from ONU in uplink
sigma_dl = 100;                	% total computation load from ONU in downlink
alpha = 1;                          % cost of a single rack in a cloudlet
D_QoS = 1*10^-3;                    % QoS latency constraint
p = 1000;                           % no. of users served by each ONU

% computing the latency constraint parameters
A = D_r+(sigma_ul*N/(n_lambda*BW_r))+(sigma_dl*N/(n_lambda*BW_r));
B = D_o+(sigma_ul*N/(BW_o-bl_ul))+(sigma_dl*N/(BW_o-bl_dl));

%% cloudlet deployment cost calculation - STAGE : 1
x_r = (1/(B-A))*((B-D_QoS) - sqrt((2*alpha*(B-A))/(varphi*(zeta_r-zeta_o))))
x_o = 1-x_r
mu = (p*N + 2*sqrt(varphi*(zeta_r-zeta_o)/(2*alpha*(B-A))))
if (x_r<=0)
    x_r = 0
    x_o = 1-x_r
    mu = (p*N + (1/(D_QoS-B)))
elseif (x_r >= 1)
    x_r = 1
    x_o = 1-x_r
    mu = (p*N + (1/(D_QoS-A)))
end
cost_1 = (alpha*mu/varphi)+zeta_r*x_r+zeta_o*x_o
total_cost_1 = cost_1*pon
rc_1 = (alpha*mu/varphi)*pon
ic_1 = (zeta_r*x_r+zeta_o*x_o)*pon

%% cloudlet deployment cost calculation - STAGE : 2
if (x_r > 0)
    % combined average cloudlet installation cost at RN and CO
    zeta_ro = zeta_r*x_r + zeta_o*x_o;
    % computing the latency constraint parameters
    P = D_f+(sigma_ul/BW_f)+(sigma_dl/BW_f);
    Q = (D_r+(sigma_ul*N/(n_lambda*BW_r))+(sigma_dl*N/(n_lambda*BW_r))) ...
        +(D_o+(sigma_ul*N/(BW_o-bl_ul))+(sigma_dl*N/(BW_o-bl_dl)));

    % computing the analytical values
    x_f = (1/(P-Q))*((Q-D_QoS) - sqrt((3*alpha*(Q-P))/(varphi*(zeta_f-zeta_ro))))
    x_r_new = (1-x_f)*x_r
    x_o_new = (1-x_f)*(1-x_r)
    if(x_f<=0)
        x_f = 0
        x_r_new = (1-x_f)*x_r
        x_o_new = (1-x_f)*(1-x_r)
    elseif(x_f>=1)
        x_f = 1
        x_r_new = 0
        x_o_new = 0
        mu = p*N+(1/(D_QoS-P))
    else
        mu = (p*N - 3*sqrt(varphi*(zeta_f-zeta_ro)/(3*alpha*(Q-P))))
    end
        cost_3 = (alpha*mu/varphi)+eta*N*D_f*10^3+zeta_f*x_f+zeta_ro*(1-x_f)
        total_cost_2 = cost_3*pon
        rc_2 = (alpha*mu/varphi)*pon
        fc_2 = (eta*N*D_f*10^3)*pon
        ic_2 = (zeta_f*x_f+zeta_ro*(1-x_f))*pon
end
