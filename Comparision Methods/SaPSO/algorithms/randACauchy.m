function r= randACauchy(varargin)
% 

medianVal = 0;
upperQuartileVal = 1;
n = rand();
% Check the arguments
if(nargin >= 1)
    medianVal=	varargin{1};
    if(nargin >= 2)
        upperQuartileVal=			varargin{2};
        upperQuartileVal(upperQuartileVal <= 0)=	NaN;	% Make NaN of out of range values.
        if(nargin >= 3),	n=	[varargin{3:end}];		
        end
    end
end

r = upperQuartileVal^2/(pi*(n^2+upperQuartileVal^2));