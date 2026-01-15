function [Gt,temporal_model] = calc_gt_yagi(n_as,n_dd,n_tstp,n_diff)

% function [Gt,temporal_model] = calc_gt_yagi(n_as,n_dd,n_tstp,n_diff)
%
% calc_gt_yagi.m - a program that calculates the smoothing matrix Gt,
%                  the prior information needed to calculate ABIC
%                   
%                  Gt is the matrix that represents smoothing in time
%                  - this version of the routine uses the numbering
%                  scheme used by Yuji Yagi for his inversions, which
%                  is also the format he supplies me the data in...
%
%                  unlike smoothing in space (2D), smoothing in time
%                  is one dimensional - between the basis functions before 
%                  and after and at the present time, for the same patch. 
%
%                  you have the choice here of 1st or 2nd difference
%                  smoothing - i prefer the latter, Yuki the former.
%
%                  take your pick
%
% gjf, 01-nov-2004
% geniusinaction.com
%
% change history
%
%   28-may-2007   gjf   functionalised for inclusion in the main program



% set up the temporal model (relates model parameters to each other in time)

cntm=0; 			% model parameter counter

for jt=1:n_tstp

   cntc=0;			% column counter

   for n=1:n_dd

      for m=1:n_as

         cntm=cntm+1;
         cntc=cntc+1;

         temporal_model(jt,cntc)=cntm;

      end

   end

end




% and do it!

if (n_diff==1)

	Gt=calc_smoothing_1stdiff_time(temporal_model) ;	
else

	Gt=calc_smoothing_2nddiff_time(temporal_model) ;

end

