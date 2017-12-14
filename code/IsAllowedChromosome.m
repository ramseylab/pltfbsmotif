function allowed=IsAllowedChromosome(chrName)
% IsAllowedChromosome - test for allowed chromosome names, for ChIP-seq peaks
%
%   allowed = IsAllowedChromosome(chrName)
%
% This function returns 0 if the argument chrName matches "chrM" or 
% "chrUn_random".  Otherwise, it returns 1.  This function is used 
% to centralize the filtering for allowed chromsome names, for
% computing genome-wide fragment count statistics and for
% identifying ChIP-seq peaks.  We filter out mitochondrial DNA
% because it tends to have very high fragment counts, which can
% potentially bias our genome-wide fragment count statistics. 
%
% Inputs:  
%  chrName - a string containing the name of the chromosome
% Outputs:
%  allowed - a Boolean (1 or 0) indicating whether the chromosome
%     named in the string variable "chrName", is allowed
%
% Subfunctions:            none
% MAT-files required:      none
% See also:                none
%

% =================================================
% IsAllowedChromosome.m
% 
% Author: Stephen A. Ramsey
%         Institute for Systems Biology
%         1441 N 34th St
%         Seattle, WA 98103 USA
%
% Copyright (C) 2009 by Institute for Systems Biology.
% All rights reserved.
% =================================================
disallowedChrs = {'chrM'; ...
		  'chrUn_random'};

allowed = ~length(find(strcmp(chrName, disallowedChrs)));


