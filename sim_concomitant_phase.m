% Simulate Concomitant FIeld Corrections at 0.55T
% clear, close all, clc

load('Experiment_gradient_waveforms.mat')
% whos

%%

% EXPERIMENT 1 
%   No First-echo comp. | No Inter-echo comp. | 60 degree step between spiral readouts
EXP1_none_none_step60       = sim_tse_maxwell_phase(EXP1_none_none_step60);

% EXPERIMENT 2 
%   Spiral First-echo comp. | No Inter-echo comp. | 0 degree step between spiral readouts
EXP2_spiral_none_step0      = sim_tse_maxwell_phase(EXP2_spiral_none_step0);

% EXPERIMENT 3
%   Bipolar First-echo comp. | No Inter-echo comp. | 60 degree step between spiral readouts
EXP3_bipolar_none_step60    = sim_tse_maxwell_phase(EXP3_bipolar_none_step60);

% EXPERIMENT 4
%   Bipolar First-echo comp. | Inter-echo comp. | 60 degree step between spiral readouts
EXP4_bipolar_full_step60    = sim_tse_maxwell_phase(EXP4_bipolar_full_step60);

% EXPERIMENT 5
%   Spiral First-echo comp. | No Inter-echo comp. | 60 degree step between spiral readouts
EXP5_spiral_none_step60     = sim_tse_maxwell_phase(EXP5_spiral_none_step60);

% EXPERIMENT 6
%   Spiral First-echo comp. | Inter-echo comp. | 60 degree step between spiral readouts
EXP6_spiral_full_step60     = sim_tse_maxwell_phase(EXP6_spiral_full_step60);

% EXPERIMENT 7
%   Spiral First-echo comp. | Inter-echo comp. | 15 degree step between spiral readouts
EXP7_spiral_full_step15     = sim_tse_maxwell_phase(EXP7_spiral_full_step15);

  