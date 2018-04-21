function [A] = createOperatorFFT(num_measurements, dimension)

% *************************************************************************
% Parameters

% Percentage of measurements taken around 0 freq. The rest is uniform samp.
PERC_NONUNIFORM = 0.1;    
% *************************************************************************


% *************************************************************************
% Sample frequencies randomly

% Number of freq. near 0 (non uniform) and elsewhere (unif)
nonUnifMeas = round(PERC_NONUNIFORM*num_measurements);
UnifMeas    = round((1-PERC_NONUNIFORM)*num_measurements);

lim_inf = ceil(0.5*nonUnifMeas);
lim_sup = dimension - floor(0.5*nonUnifMeas);

% Near the origin: no sampling
coeff_near_origin = [1: lim_inf, lim_sup + 1: dimension];

% Far from origin: uniform sampling
possible_choices_frequencies = lim_inf + 1 : lim_sup;

len_possible_choices_frequencies = ...
    length(possible_choices_frequencies);

perm_vec = randperm(len_possible_choices_frequencies);

coeff_far_origin = possible_choices_frequencies(perm_vec(1:UnifMeas));

ind_f = sort([coeff_near_origin, coeff_far_origin])';

Phi = opFoG(opRestriction(dimension, ind_f), opFFT(dimension));
% *************************************************************************

A = Phi;
