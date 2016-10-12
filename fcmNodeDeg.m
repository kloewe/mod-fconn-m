function k = fcmNodeDeg(data,thr,varargin)
%FCMNODEDEG Compute the degree of each node in a func. connectivity graph
%
%   K = FCMNODEDEG(F,T) is a vector consisting of the degrees of each node
%   in the functional connectivity graph (implicitly) constructed based on
%   the functional data F and the binarization threshold T.
%   F is an M-by-N array, where rows and columns correspond to observations
%   and variables (nodes), respectively, i.e., there are M observations for
%   each of the N variables (nodes).
%
%   K = FCMNODEDEG(F,T,M) specifies the measure of association M, which is
%   used to estimate pairwise functional connectivity (during graph
%   construction).
%
%     Value                  Description
%     -----                  -----------
%     'pcc'  (default)       Pearson sample correlation coefficient
%     'tetracc'              tetrachoric correlation coefficient
%
%   K = FCMNODEDEG(F,T,P) specifies the number of threads P to be used for
%   parallelization through multi-threading.
%
%     Value                  Description
%     -----                  -----------
%     -1     (default)       auto-determine
%      0                     single-threaded version
%      1-n                   multi-threaded version with P threads
%
%   K = FCMNODEDEG(F,T,M,P) specifies both M and P.
%
%   K = FCMNODEDEG(F,T,M,P,C) also specifies the tile/block size C.
%   If C is not specified, it is auto-determined.
%
%     Value                  Description
%     -----                  -----------
%     -1     (default)       auto-determine
%      0                     on-demand computation (no cache)
%      0 < C < N             cache-based computation with tile size C
%      N                     precompute all (half-stored)
%
%   Author: Kristian Loewe

% inits and defaults
m       = 'pcc';              % FC measure
c       = -1;                 % tile/block size for caching
p       = -1;                 % number of threads

% input
if nargin == 3                % fcmNodedeg(F,T,M) or fcmNodeDeg(F,T,P)
  if ischar(varargin{1})
    m = varargin{1};          % fcmNodedeg(F,T,M)
  else
    p = varargin{1};          % fcmNodeDeg(F,T,P)
  end
elseif nargin > 3             % fcmNodeDeg(F,T,M,P) or fcmNodeDeg(F,T,M,P,C)
  m = varargin{1};
  p = varargin{2};
end
if nargin > 4                 % fcmNodeDeg(F,T,M,P,C)
  c = varargin{3};
end

% threshold
assert(isfloat(thr) && thr >= -1.0 && thr <= 1.0);

% FC measure
assert(ischar(m));
switch m
  case 'pcc'
    m = int32(1); % FCM_PCC
  case {'tetracc', 'tcc'}
    m = int32(2); % FCM_TCC
  otherwise
    error('Unknown measure of association.');
end

% number of threads
assert(isscalar(p) && isnumeric(p) && isreal(p) && p >= -1 && p <= 1024);
p = min(p, size(data,2));
p = int32(p);

% cache parameter (tile/block size)
assert(isscalar(c) && isnumeric(c) && isreal(c) && c >= -1 && ...
  c <= size(data,2));
c = min(c, size(data,2));
c = int32(c);

% compute degrees
switch class(data)
  case 'single'
    k = mxNodeDegFlt(data, single(thr), m, p, c); % 2 -1 -1
  case 'double'
    k = mxNodeDegDbl(data, double(thr), m, p, c);
  otherwise
    error('Unexpected data type.');
end
k = cast(k, class(data));

end
