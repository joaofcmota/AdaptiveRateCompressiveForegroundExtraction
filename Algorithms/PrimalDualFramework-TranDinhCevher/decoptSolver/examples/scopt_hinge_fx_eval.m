%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION: [fx, gradfx, qx_out, px_out] = scopt_logistic_fx_eval(x, ...
%%%                              Ax, What, bhat, N, yhat, NWThat, isFxEval) 
%%% PURPOSE:  Evaluate the objective value/gradient of the function.
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f_val,fx,gx,accuracy] = scopt_hinge_fx_eval(Xtest, W, bias,lambda, y, l1)



fx      = lambda*sum(max(1 - y.*(Xtest*W + bias),0));
if l1
        gx =        norm(W,1);
else
        gx =        0.5 * norm(W,2)^2;
end

f_val       =   lambda*fx +gx;

accuracy    = 1 - sum(sign(Xtest*W + bias) ~= y) / numel(y);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF THE IMPLEMENTATION.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%