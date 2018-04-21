%PhaseTransitionLookup.m
%
%Garrett Warnell
%July 2013
%
%DESCRIPTION:
%
%INPUTS:
%    s: signal sparsity
%
%    N: signal dimension
%
%    delta_samples: the values of delta used to generate pt_diagram,i.e.,
%    pt_diagram(i,:) corresponds to the grid points generated using
%    delta=delta_samples(i)
%
%    rho_samples: the values of rho used to generate pt_diagram, i.e.,
%    pt_diagram(:,j) corresponds to the grid points generated using the
%    rho=rho_samples(j)
%
%    pt_diagram: the phase transition diagram, i.e., a matrix where the
%    (i,j)th entry corresponds to the probability of success for the
%    problem instance defined by delta_samples(i) and rho_samples(j)
%
%OUTPUTS:
%
%    M: the number of compressive samples needed to ensure recovery for
%    this (s,N) pair based on the phase transition analysis
%
%NOTES:

function M = PhaseTransitionLookup(s,N,delta_pts,rho_pts,pt_diagram)

%define success threshold
success_thresh = 0.9;

%translate deltas to M's using the specific value of N and find possible
%values to consider based on "s"
M_pts = round(N*delta_pts);
M_possible_idx = find(M_pts>=s);

%loop through M_possible_idx until one is found where "s" was successful
M = N;
for i=1:length(M_possible_idx)
    %translate the rho values to sparsities
    s_samples = round(M_pts(M_possible_idx(i))*rho_pts);
    %determine if any of these satisfy the recovery probability threshold
    if(any( (s_samples>=s)&(pt_diagram(M_possible_idx(i),:)>=success_thresh) ))
        M = M_pts(M_possible_idx(i));
        break;
    end
end