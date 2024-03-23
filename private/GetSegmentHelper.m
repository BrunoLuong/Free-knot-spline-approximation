function segidx = GetSegmentHelper(j,l,breaks,jmin)
% segidx = GetSegmentHelper(j,l,breaks,jmin)
% Return segment index of the query points located inside the intervals
% indices by (j:j+l)
% j is sub-interval relative to full knot t
% l is spline order
% breaks is array as returned by GetBreaksArray
% jmin is started index of knot where breaks is computed
%
% See also: GetBreaksArray

jM = j-jmin+1; % jM is sub-interval relative to restricted knot t(jmin:...)
segidx = breaks(jM):breaks(jM+l)-1; 

end % GetSegmentHelper
