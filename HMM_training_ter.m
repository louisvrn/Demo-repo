clear
close all
clc

%% Load data to workspace

datasets.parameters = dir('model_parameters_nstd0.mat'); 
for i = 1:length(datasets.parameters)
    data(i) = load(datasets.parameters(i).name, '-mat');
    cycles{i,:} = 1:length(data(i).nSTD);
    nSTD{i}(:,cycles{i,:}) = data(i).nSTD;
end

%% Gaussian fitting
%Now let us introduce the initial Gaussian parameters for Q = 2 hidden states for participant 1 to
%base first iteration of HMM off
%The mean and standard deviation of the Gaussians are used to define
%emission probabilities of the HMM 
Q = 2; 

mu0 = [data(i).m1 data(i).m2];
Sigma0 = [data(1).s1 data(1).s2];

mu0 = reshape(mu0, [1 Q 1]);
Sigma0 = reshape(Sigma0, [1 1 Q 1]);

%Define prior and transition probabilities, the mixmatrix is null
%in this case since we are not using a mixture of Gaussian inputs to 
%represent each state but only one

mixmat0 = [];
prior0 = [1; 0];
transmat0 = data(1).transmat;

%% Iteratively train model
%Finally, let us iteratively improve these parameter estimates using EM
%and the computed STD per contraction

[LL_HMM, prior_HMM, transmat_HMM, mu_HMM, Sigma_HMM, mixmat_HMM] = mhmm_em(nSTD, prior0, transmat0, mu0, (Sigma0).^2, mixmat0, 'max_iter', 20);
 
save('HMM_parameters', 'prior_HMM', 'transmat_HMM', 'mu_HMM', 'Sigma_HMM', 'mixmat_HMM');