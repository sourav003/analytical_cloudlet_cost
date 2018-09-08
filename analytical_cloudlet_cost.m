%% cloudlet deployment cost calculation using analytical formulae
% Author: Sourav Mondal, EEE, The University of Melbourne
% =======================================================
function [cost,x_f,x_r,x_o] = analytical_cloudlet_cost(N,R,gamma,D_QoS)
    % ---------------------------------------------------------
    % N = split ratio of TDM-PON
    % R = radius of circular area (km)
    % gamma = population density (ONUs/sq.km)
    % D_QoS = QoS latency constraint (sec)
    % ---------------------------------------------------------
    pon = ceil(round(gamma*pi*R^2)/N);                   % total number of PON branches
    theta = 2*pi/pon;                             % inscribed angle of the circular wedge
    phi = theta/2;
    varphi = 2500;                                % no. of VMs per cloudlet rack
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
    zeta_f = 8;                         % normalized cost of cloudlet installation
    zeta_r = 6;
    zeta_o = 5;
    eta = 20;                           % cost of fiber installation per km
    sigma_ul = 1*10^6;            	    % total computation load from ONU in uplink
    sigma_dl = 50*10^3;                	% total computation load from ONU in downlink
    alpha = 1;                          % cost of a single rack in a cloudlet
    p = 1000;                           % no. of users served by each ONU
    
    % computing the latency10 constraint parameters for Stage - 1
    A = D_r+(sigma_ul*N/(n_lambda*BW_r))+(sigma_dl*N/(n_lambda*BW_r));
    B = D_o+(sigma_ul*N/(BW_o-bl_ul))+(sigma_dl*N/(BW_o-bl_dl));
    
    if ((B < (D_QoS+sqrt(2*alpha*(B-A)/(varphi*(zeta_r-zeta_o))))) && (A <= (D_QoS+sqrt(2*alpha*(B-A)/(varphi*(zeta_r-zeta_o))))))
        cost_pon = (alpha/varphi)*(p*N+(1/(D_QoS-B)))+zeta_o
        x_f = 0; x_r = 0; x_o = 1-x_r;
        disp('stage-1 new if');
        disp('printing total cloudlet installation cost:');
        cost = pon*(cost_pon-zeta_o)+zeta_o;
        
    elseif (B >= (D_QoS+sqrt(2*alpha*(B-A)/(varphi*(zeta_r-zeta_o))))) && (A <= (D_QoS+sqrt(2*alpha*(B-A)/(varphi*(zeta_r-zeta_o)))))
        cost_pon = (alpha/varphi)*(p*N-2*sqrt(varphi*(zeta_r-zeta_o)/(2*alpha*(B-A))))+ ...
            ((zeta_r-zeta_o)/(B-A))*(B-D_QoS-sqrt(2*alpha*(B-A)/(varphi*(zeta_r-zeta_o))))+zeta_o;
        x_f = 0; x_r = (1/(B-A))*(B-D_QoS-sqrt(2*alpha*(B-A)/(varphi*(zeta_r-zeta_o)))); x_o = 1-x_r;
        disp('stage-1 new elseif');
        disp('printing total cloudlet installation cost:');
        cost = pon*(cost_pon-zeta_o)+zeta_o;
    end
    % entry point - Stage - 2
    if ((A > (D_QoS+sqrt(2*alpha*(B-A)/(varphi*(zeta_r-zeta_o))))) && (B >= (D_QoS+sqrt(2*alpha*(B-A)/(varphi*(zeta_r-zeta_o))))))
        eps = A/B;
        % computing the latency constraint parameters for Stage - 2
        zeta_ro = zeta_r*eps + zeta_o*(1-eps);
        P = D_f+(sigma_ul/(BW_f))+(sigma_dl/(BW_f));
        Q = eps*(D_r+(sigma_ul*N/(n_lambda*BW_r))+(sigma_dl*N/(n_lambda*BW_r)))...
            +(1-eps)*(D_o+(sigma_ul*N/(BW_o-bl_ul))+(sigma_dl*N/(BW_o-bl_dl)));
    
        x_f = (1/(Q-P))*(Q-D_QoS-sqrt(3*alpha*(Q-P)/(varphi*(zeta_f-zeta_ro))));
        x_r = (1-x_f)*eps;
        x_o = (1-x_f)*(1-eps);
        cost_pon = (alpha/varphi)*(p*N-3*sqrt(varphi*(zeta_f-zeta_ro)/(3*alpha*(Q-P))))...
                    +x_f*eta*cent_avg*N + ((zeta_f-zeta_ro)/(Q-P))*(Q-D_QoS-sqrt(3*alpha*(Q-P)/(varphi*(zeta_f-zeta_ro))))+zeta_ro;
        disp('stage-2 second if');
        disp('printing total cloudlet installation cost:');
        cost = pon*(cost_pon-zeta_o)+zeta_o;
    end
end