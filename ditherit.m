
function x = ditherit(x,db,type)
% DITHERIT Dithers a signal by db dBs.
%   x = ditherit(x,db,type)
%
%   Dithers a signal by db dBs
%
%   x:    input signal
%   db:   dithering amount, e.g. -96
%   type: type of dithering
%   x:    dithered signal

% ------- ditherit.m ---------------------------------------
% Marios Athineos, marios@ee.columbia.edu
% http://www.ee.columbia.edu/~marios/
% Copyright (c) 2002 by Columbia University.
% All rights reserved.
% ----------------------------------------------------------

% Check arguments
if (nargin < 2); db =   -96;  end
if (nargin < 3); type = 'db'; end

switch lower(type)
    case {'db'}
%        randn('state',0);
        x = x + (10^(db/20))*randn(size(x));
    case {'bit'}
        x = x + round(2*db*rand(size(x))-db);
end