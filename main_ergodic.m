%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Ma15cap] Capacity analysis of lattice reduction aided equalizers for MIMO systems
% this program verifies p1 and p2 in [Ma15cap]
% p1: C^{ml} = C_LR^{ml} >= C_LR^{zf}
% p2: ergodic capacity
%     E_H(C^{ml}) >= E_H(C_LR^{zf}) >= E_H(C^{zf})
%
% Written by: Yiming Kong
% Date: 3/1/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
addpath('./lr', './lib');
SNR = 0 : 2 : 20;

sim_n = 1e6;
Nr = 8; Nt = 8;

% to compute ergodic capacity
E_ml = zeros(1, length(SNR));
E_zf = zeros(1, length(SNR));
E_clll_zf = zeros(1, length(SNR));
E_Dclll_zf = zeros(1, length(SNR));
E_SA_zf = zeros(1, length(SNR));
E_elr_zf = zeros(1, length(SNR));

% rand seed
randn('state', sum(100 .* clock));
rand('state', sum(100 .* clock));

tic
for runs = 1 : sim_n
    H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) / sqrt(2);
    C = inv(H' * H);
    R_eta_v = diag(C);
    invR_eta = diag(1 ./ R_eta_v);
    
    % clll
    [H_t, T] = clll_c(H);
    C_hat = inv(H_t' * H_t);
    R_eps_v = diag(C_hat);
    invR_eps = diag(1 ./ R_eps_v);
    invTh = inv(T');
    invT = inv(T);

    % dclll
    [H_t2, T2] = dclll_c(H);
    C_hat2 = inv(H_t2' * H_t2);
    R_eps_v2 = diag(C_hat2);
    invR_eps2 = diag(1 ./ R_eps_v2);
    invTh2 = inv(T2');
    invT2 = inv(T2);

    % SA
    [H_t3, T3] = sa_c(H);
    C_hat3 = inv(H_t3' * H_t3);
    R_eps_v3 = diag(C_hat3);
    invR_eps3 = diag(1 ./ R_eps_v3);
    invTh3 = inv(T3');
    invT3 = inv(T3);

    % ELR dual
    [H_t4, T4] = elr_dual_c(H);
    C_hat4 = inv(H_t4' * H_t4);
    R_eps_v4 = diag(C_hat4);
    invR_eps4 = diag(1 ./ R_eps_v4);
    invTh4 = inv(T4');
    invT4 = inv(T4);
        
    for SNR_ind = 1 : length(SNR)
        sigma = sqrt(1 / 10 .^ (SNR(SNR_ind) ./ 10));
        
        a = log2(det(eye(Nt) + 1 / sigma^2 .* H' * H));
        b = log2(det(eye(Nt) + 1 / sigma^2 .* invR_eta));
        
        c1 = log2(det(eye(Nt) + 1 / sigma^2 .* invTh * invR_eps * invT));
        c2 = log2(det(eye(Nt) + 1 / sigma^2 .* invTh2 * invR_eps2 * invT2));
        c3 = log2(det(eye(Nt) + 1 / sigma^2 .* invTh3 * invR_eps3 * invT3));
        c4 = log2(det(eye(Nt) + 1 / sigma^2 .* invTh4 * invR_eps4 * invT4));
        
        E_ml(1, SNR_ind) = E_ml(1, SNR_ind) + a;
        E_zf(1, SNR_ind) = E_zf(1, SNR_ind) + b;
        
        E_clll_zf(1, SNR_ind) = E_clll_zf(1, SNR_ind) + c1;
        E_Dclll_zf(1, SNR_ind) = E_Dclll_zf(1, SNR_ind) + c2;
        E_SA_zf(1, SNR_ind) = E_SA_zf(1, SNR_ind) + c3;
        E_elr_zf(1, SNR_ind) = E_elr_zf(1, SNR_ind) + c4;
        
    end
    if(mod(runs, 5000) == 0)
        fprintf('left: %0.2f\n', toc / runs * (sim_n - runs));
    end
end

figure(1)
plot(SNR, real(E_ml) / sim_n, 'k>-', 'LineWidth', 1.5); hold on
plot(SNR, real(E_zf) / sim_n, 'r<-', 'LineWidth', 1.5); hold on
plot(SNR, real(E_clll_zf) / sim_n, 'b*-', 'LineWidth', 1.5); hold on
plot(SNR, real(E_Dclll_zf) / sim_n, 'bo-', 'LineWidth', 1.5); hold on
plot(SNR, real(E_SA_zf) / sim_n, 'y^-', 'LineWidth', 1.5); hold on
plot(SNR, real(E_elr_zf) / sim_n, 'ms-', 'LineWidth', 1.5); hold on
legend('MLE', 'ZF', 'CLLL-ZF','DCLLL-ZF','SA-ZF','ELR-ZF');
name = sprintf('Ecap_%dx%d_sim%0.1e', Nr, Nt, sim_n); 
xlabel('SNR(dB)'); 
ylabel('Ergodic Capacity'); 
title(name);
tmp_fn = sprintf('_%s', datestr(now, 'yymmddTHHMMSS'));
saveas(gcf, strcat(name, tmp_fn, '.fig'));