function [k,offset] = checksummarization(K)
% [K,OFFSET] = CHECKSUMMARIZATION(K) - Check that K is a known summarization type {Center/Middle, 
%   End, Instant, or Beginning/Start). Parse into a single char {c,e,i,b}, and return OFFSET 
%   (distance from the start of the interval to the time label, in units of the time step dT).
%
% INPUT: K - char key ('C','b','e','i',..), word ('center','start','EOI','Instant',...), or cellstr.
% See also: PARSETIME, CHECKTIMEZONE

    KEYS = 'bice';
    OFF = [0,0,0.5,1];

    if iscell(K), [k,offset] = cellfun(@checksummarization,K); return; end
    
    if isempty(K) || ~ischar(K), K = '?'; end  % crash below
    if numel(K) == 2, [k,offset] = arrayfun(@checksummarization,K); return; end
    
    K = lower(K);
    if ~isscalar(K)
        switch K
        case {'beginning','start','soi','boi','begin'}, K = 'b';
        case {'instant','instantaneous','sample','point'}, K = 'i'; 
        case {'end','eoi'},K = 'e';
        case {'center','coi','middle','moi','mid'}, K = 'c';
            otherwise, K = '?';
        end
    end
    
    match = (K == KEYS);
    assert(any(match),'Unknown interval summarization key');
    
    k = KEYS(match);
    offset = OFF(match);
end