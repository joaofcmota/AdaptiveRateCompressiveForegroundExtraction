%SPG_Interface.m
%
%LAST MODIFIED: Garrett Warnell, 1/2011
%
%DESCRIPTION:
%   switches between passing an input to two function handles based on a
%   second input.  useful for the (x,mode) interface necessary for the SPGL
%   library
%
%INPUTS:
%   x: the input to be passed to either function
%
%   mode: interface function will evaluate to f1(x) if mode==1 and f2(x) if
%   mode==2
%
%   f1: function handle to the first function
%
%   f2: function handle to the second function
%
%OUTPUTS:
%   y: (output format dictated by output of selected function)
%
%NOTES:
%
%REFERENCES:

function y = SPG_Interface(x,mode,f1,f2)

if(mode==1)
    y = f1(x);
elseif(mode==2)
    y = f2(x);
else
    error('SPG_Interface.m: mode should be 1 or 2');
end