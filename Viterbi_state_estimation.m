clear 
close all
clc
%% Estimating state at every observation for a given data sequence
%load observation data for one participant to test the model
observation = load('model_parameters_nstd0.mat');
data = observation.nSTD;
LT2 = observation.LT_per_contraction; %LT found previously


%%
%load model parameters
model = load('HMM_parameters.mat');
mu = model.mu_HMM;
Sigma = model.Sigma_HMM;
mixmat = model.mixmat_HMM;
transmat = model.transmat_HMM;

%% Use Viterbi algorithm to predict hidden state for every contraction in the input data
B = mixgauss_prob(data, mu, Sigma);
[path] = viterbi_path([1;0], transmat, B);

%% Find contraction at which hidden state transition is predicted to have been made (if at all) 
curr_state = path(end);
transition = find(path>1, 1, 'first');

