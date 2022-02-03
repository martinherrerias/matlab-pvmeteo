function [Z,x,y,F] = series2daymatrix(t,z,varargin)
% [Z,x,y,F] = SERIES2DAYMATRIX(t,z) - Arrange time series (t,z) as a (24 hour / dt)·M array Z,
%   where dt is the auto-detected time-step resolution of t; the M columns of Z correspond to
%   data-points of the same day; and the (24 hour / dt) rows  of Z correspond to data-points at
%   the same time of day.
%
% [..] = SERIES2DAYMATRIX(..,'name',value) - pass additional options to PARSETIME, such as a 
%   known time-step (e.g. 'step',minutes(5),..), interval summarization (e.g. 'interval','c'),
%   time-zone, date-string format, etc.
%
% INPUT:
%   t - N-vector of timestamps (DATENUM, DATETIME, DATEVEC,...) passed to PARSETIME
%   z - N-vector of values
%
% OUTPUT: 
%   Z - (24·60/dt)·M array with the rearranged values of z.
%   x - horizontal M-vector of DATENUM day timestamps (integer days), for the columns of Z
%   y - vertical 24·60/dt-vector of day-times (i.e. day-fractions [0,1]) for the rows of Z
%   F - (24·60/dt)·M x N sparse matrix, mapping z into Z, i.e. Z(:) = F*z
%
% See also: SERIESHEATMAP, DATENUM

    narginchk(2,Inf);
    assert(isequal(numel(t),numel(z)),'Incompatible t and z sizes');
    
    [t,dt,idx,ia] = parsetime(t,'-regular','latticeunit',days(1),varargin{:});
    N = numel(ia);
        
    t0 = dateshift(t(1),'start','day');
    padding = fliplr(t(1)-dt:-dt:t0)';
    t = [padding;t];
    idx = idx + numel(padding);
    
    TZ = t.TimeZone;
    t.TimeZone = 'UTC'; % must use UTC, otherwise it might skip steps due to DST

    % offset = t(1) - t0 - dt/2;
    t_end = t(1) + days(ceil(days(t(end)-t(1)+dt))) - dt/2;
    
    t = [t;(t(end)+dt:dt:t_end)'];
    
    t.TimeZone = TZ;   
    M = numel(t);
    
    n = accumarray(ia,1);
    F = sparse(idx(ia),1:N,1./n(ia),M,N) + sparse(setdiff(1:M,idx),1,NaN,M,N);
    z = full(F*double(z(:)));

    % t = datenum(t);
    y = (datenum(t(1)):days(dt):datenum(t(1)+1 - dt/2)) - datenum(t0);
    x = floor(datenum(t(1:1/days(dt):end)));
    Z = reshape(z,[],numel(x));
end
