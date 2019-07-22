T = 1;
r = 0.1;
sigma = 0.5;
K = 5;
S0 = 5;
N_sample = 100000;

% 1.1 Antithetic variates
W_T = normrnd(0,sqrt(T),N_sample,1);
W_T_neg = -W_T;
S_T_1 = S0*exp((r-0.5*sigma^2)*T+sigma.*W_T);
S_T_2 = S0*exp((r-0.5*sigma^2)*T+sigma.*W_T_neg);
g_S_T_1 = zeros(N_sample,1);
g_S_T_2 = zeros(N_sample,1);
g_S_T = zeros(N_sample,1);
for i = 1:N_sample
    g_S_T_1(i,1) = exp(-r*T)*(subplus(S_T_1(i,1)-K)+subplus(K-S_T_1(i,1)));
    g_S_T_2(i,1) = exp(-r*T)*(subplus(S_T_2(i,1)-K)+subplus(K-S_T_2(i,1)));
    g_S_T(i,1) = (g_S_T_1(i,1)+g_S_T_2(i,1))/2;
end
[mean_an,sigma_an,muci_an,sigmaci_an] = normfit(g_S_T);
std_error_an = sigma_an/sqrt(N_sample);

% 1.2 Control variates
W_T = normrnd(0,sqrt(T),N_sample,1);
S_T = S0*exp((r-0.5*sigma^2)*T+sigma.*W_T);
g_S_T = zeros(N_sample,1);
for i = 1:N_sample
    g_S_T(i,1) = exp(-r*T)*(subplus(S_T(i,1)-K)+subplus(K-S_T(i,1)));
end
[mean_s,std_s] = normfit(S_T);
cov_S_g_matrix = cov(S_T,g_S_T);
cov_S_g = cov_S_g_matrix(1,2);
c = -cov_S_g/(std_s)^2;
x_c = zeros(N_sample,1);
for i = 1:N_sample
    x_c(i,1) = g_S_T(i,1) + c*(S_T(i,1)-mean_s);
end
[mean_x_c,sigma_x_c,muci_x_c,sigmaci_x_c] = normfit(x_c);
std_error_x_c = sigma_x_c/sqrt(N_sample);

% 1.3 Halton
p = haltonset(1);
U = p(2:N_sample+1,:);
U_1 = U(1:N_sample/2);
U_2 = U((1+N_sample/2):N_sample);
X_i = zeros(N_sample/2,1);
Y_i = zeros(N_sample/2,1);
for i = 1:(N_sample/2)
    X_i(i,1) = sqrt(-2*log(U_1(i,1)))*cos(2*pi*U_2(i,1));
    Y_i(i,1) = sqrt(-2*log(U_1(i,1)))*sin(2*pi*U_2(i,1));
end
W_T_h = sqrt(T)*[X_i;Y_i];
S_T_h = S0*exp((r-0.5*sigma^2)*T+sigma.*W_T_h);
g_S_T_h = zeros(N_sample,1);
for i = 1:N_sample
    g_S_T_h(i,1) = exp(-r*T)*(subplus(S_T_h(i,1)-K)+subplus(K-S_T_h(i,1)));
end
[mean_h,sigma_h,muci_h,sigmaci_h] = normfit(g_S_T_h);
std_error_h = sigma_h/sqrt(N_sample);
    



