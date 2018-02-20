function varargout = refwrap(func,varargin)
%REFWRAP Wraps reference functions used for tests
%
%  Author: Kristian Loewe

switch func

  case 'mdiff'
    s = mean(varargin{1}) - mean(varargin{2});
  
  case {'tstat', 'pairedt'}
    [~,p,~,stats] = ttest(varargin{:});
    s = stats.tstat;
  
  case 'tstat2'
    [~,p,~,stats] = ttest2(varargin{:});
    s = stats.tstat;
    
  case 'corrv'
    [s,p] = corrcoef(varargin{:});
    s = s(2);
    p = p(2);
  otherwise
    error('Unexpected input function handle.');
end

varargout{1} = s;

if nargout == 2
  varargout{2} = p;
end

end
