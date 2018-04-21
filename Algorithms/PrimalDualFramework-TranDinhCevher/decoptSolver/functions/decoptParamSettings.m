% FUNCTION: param = decoptParamSettings(varargin)
%
% PURPOSE:  Generate optional parameters for algorithms.
% 
% CALL: 
%    help    decoptParamSettings: see optional parameters.   
%    param = decoptParamSettings() to see the default values.
%    param = decoptParamSettings('Name', Value, ...) to set new param.
%    param = decoptParamSettings(param, 'Name', Value, ...) to add new param.
%
% NOTES:
%    MaxIters        : The maximum number of iterations.
%    RelTolX         : The relative tolerance for the search direction.
%    RelTolFun       : The relative tolerance for the objective value.
%    TerminationType : Type of termination conditions.', char(13), ...
%                      (0 = default, 1 = norm of the search
%                      direction, 2 - objective value change, 3 = 1 + 2).
%    isFxEval        : Evaluating the objective value if required.
%    isFGapCheck     : Check the increase of the feasibility gap
%    numFxCheck      : Number of checking iterations in the objective value.
%    Verbosity       : Output printing level.
%    isStoreBest     : Store the best value so far or not.
%    saveHistMode    : Save the history information (0 = non and 1 = details, 
%                      2 = include the iterative vector).
%    Gamma           : The initial value for smoothness parameter.
%                      Default value is 1e-6.
%    GammaRatio      : The ratio for increasing 'Gamma', default: 1.05.
%    MaxGamma        : The upper bound on 'Gamma', default: 1e30.
%    BetaRatio       : The ratio for decreasing 'Beta', default: 0.5.
%    PrintStep       : Number of iterations for printing output (50).
%    PwMaxIters      : Maximum number of iterations for Power method.
%
% INFORMATION:
%    By Quoc Tran Dinh, Laboratory for Informations and Inference Systems,
%       (LIONS), EPFL, Lausanne, Switzerland.
%    Joint work with Volkan Cevher.
%    Date: 20.12.2013.
%    Last modified: 20.12.2013.
%    Contact: quoc.trandinh@epfl.ch
%
function param = decoptParamSettings(varargin)

% Define the default param.
defaultopt = struct( ...
                     'Algorithm',       1,          ...
                     'MaxIters',        5000,       ...
                     'RelTolX',         1e-5,       ...
                     'RelTolFun',       1e-7,       ...
                     'RelTolFeas',      1e-7,       ...
                     'isFxEval',        true,       ...
                     'isFGapCheck',     false,      ...
                     'isRestart',       true,       ...
                     'numFxCheck',      5,          ...
                     'TerminationType', 1,          ...
                     'Verbosity',       2,          ...
                     'saveHistMode',    0,          ...
                     'Gamma',           1e-6,       ...
                     'GamRatio',        1.05,       ...
                     'MaxGamma',        1e100,      ...
                     'isStoreBest',     false,      ...
                     'BetaRatio',       0.5,        ...
                     'PrintStep',       50,         ...
                     'PwMaxIters',      30,         ...
                     'PwRelTol',        1e-4,       ...
                     'defDualScale',    1,          ...
                     'isOper2Mat',      true,       ...
                     'maxMatSize',      2048,       ...
                     'adaptStepSize',   true,       ...
                     'isRescaleData',   false,      ...
                     'MinMaxNormA',     [0.75, 100],...
                     'InnerMaxIters',   100,        ...
                     'InnerRelTol',     1e-7,       ...
                     'isRestartFISTA',  false       ...
                   );
             
% In case no input or empty input.
if nargin < 1,        
    param = defaultopt; 
    printOptions();
    return; 
end
if isempty(varargin{1})
    param = defaultopt; 
    return; 
end

% In case varargin{1} is a structure.
if isstruct(varargin{1})
    param = varargin{1};
    if nargin > 1
        inputs = varargin(2:end);
    else
        inputs = [];
    end
else
   param = defaultopt; 
   inputs = varargin;
end
if mod(length(inputs), 2) ~= 0
    error('Optional parameters must be in pair (Name, Value)!');
end

% Now, add the new param.
param = addNewOptions(param, inputs);

% Check if inputs are correct.
param = checkfield(param, defaultopt, 'Algorithm' );
param = checkfield(param, defaultopt, 'MaxIters' );
param = checkfield(param, defaultopt, 'RelTolX' );
param = checkfield(param, defaultopt, 'RelTolFun' );
param = checkfield(param, defaultopt, 'RelTolFeas' );
param = checkfield(param, defaultopt, 'isFxEval' );
param = checkfield(param, defaultopt, 'numFxCheck' );
param = checkfield(param, defaultopt, 'Verbosity' );
param = checkfield(param, defaultopt, 'TerminationType');
param = checkfield(param, defaultopt, 'saveHistMode');
param = checkfield(param, defaultopt, 'Gamma' );
param = checkfield(param, defaultopt, 'GamRatio' );
param = checkfield(param, defaultopt, 'MaxGamma' );
param = checkfield(param, defaultopt, 'BetaRatio' );
param = checkfield(param, defaultopt, 'PrintStep' );
param = checkfield(param, defaultopt, 'PwMaxIters' );
param = checkfield(param, defaultopt, 'isFGapCheck' );
param = checkfield(param, defaultopt, 'isStoreBest' );
param = checkfield(param, defaultopt, 'PwRelTol' );
param = checkfield(param, defaultopt, 'defDualScale' );
param = checkfield(param, defaultopt, 'isOper2Mat' );
param = checkfield(param, defaultopt, 'maxMatSize' );
param = checkfield(param, defaultopt, 'isRescaleData' );
param = checkfield(param, defaultopt, 'MinMaxNormA' );
param = checkfield(param, defaultopt, 'isRestart' );
param = checkfield(param, defaultopt, 'adaptStepSize' );
param = checkfield(param, defaultopt, 'InnerMaxIters' );
param = checkfield(param, defaultopt, 'InnerRelTol' );
param = checkfield(param, defaultopt, 'isRestartFISTA' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION: outopt = checkfield(outopt, default, fieldname) 
%%% PURPOSE:  This function checks the field of optional parameters.
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outopt = checkfield(outopt, default, fieldname)
    
    fieldvalue = getfield( default, fieldname );
    if isfield(outopt, fieldname )
        optfieldval = getfield(outopt, fieldname );
        if isempty( optfieldval )
            outopt = setfield( outopt, fieldname, fieldvalue );
        end;
    else
        outopt = setfield( outopt, fieldname, fieldvalue );
    end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION: param = addNewOptions(param, inputs)
%%% PURPOSE:  Add new param to the structure.
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function param = addNewOptions(defin, inputs)

param = defin;
if isempty(inputs), return; end

for k=1:length(inputs)/2
    name  = inputs{2*k-1};
    value = inputs{2*k};
    if ~ischar(name)
        error('Name of an option must be a string!');
    end
    if ~isfield(param, name)
        error('This optional parameter does not exist! Add a new one.');
    end
    param = setfield(param, name, value);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FUNCTION: printOptions()
%%% PURPOSE:  Print the param.
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printOptions()

fprintf('  MaxIters      : The maximum number of iterations.\n');
fprintf(['  RelTolX       : The relative tolerance for the search ', ...
         'direction.\n']);
fprintf(['  RelTolFun     : The relative tolerance for the objective ', ...
         'values.\n']);
fprintf('  isFxEval      : Evaluating the objective value if required.\n');
fprintf('  isFGapCheck   : Check the increase of the feasibility gap.\n');
fprintf('  isStoreBest   : Store the best value so far or not?\n');
fprintf(['  numFxCheck    : Number of checking iterations in the ', ...
         'objective value (using in the stopping criterion).\n']);
fprintf( '  Verbosity     : Option for printing output.\n');
fprintf(['  TerminateType : Type of termination conditions.', char(13), ...
         '                  (0 = default, 1 = norm of the search', ...
         ' direction, 2 - objective value change, 3 = 1 + 2).\n']);
fprintf(['  saveHistMode  : Save the history information mode (0 = non', ...
         ', 1 = details, 2 - iterative vector).\n']);
fprintf(['  Gamma         : The initial value for smoothness parameter.', ...
         ' Default value is 1e-6.\n']);
fprintf('  GammaRatio    : The ratio for increasing Gamma, default: 1.05.\n');
fprintf('  MaxGamma      : The upper bound on Gamma, default: 1e30.\n');
fprintf('  BetaRatio     : The ratio for decreasing Beta, default: 0.5.\n');
fprintf('  PrintStep     : Number of iterations for printing output (50).\n');
fprintf('  PwMaxIters    : Maximum number of iterations for Power method.\n');
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% END OF THE IMPLEMENTATION.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%