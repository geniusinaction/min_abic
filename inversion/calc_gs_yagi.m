function [Gs,spatial_model]=calc_gs_yagi(ingeomfile,n_tstp,n_diff)

% function [Gs,spatial_model]=calc_gs_yagi(ingeomfile,n_tstp,n_diff)
%
% calc_gs_yagi.m - a program that calculates the smoothing matrix Gs,
%                  the prior information needed to calculate ABIC
%                   
%                  Gs is the matrix that represents smoothing in space
%                  - this version of the routine uses the numbering
%                  scheme used by Yuji Yagi for his inversions, which
%                  is also the format he supplies me the data in...
%
%                  you have the choice here of 1st or 2nd difference
%                  smoothing - i prefer the latter, Yuki the former
%
%                  take your pick
%
% gjf, 31-oct-2004 (woooooo)
% geniusinaction.com
%
% modifications:
%
% 	15-dec-2004	gjf	sparse matrices now handled, yagi-style smoothing too
% 	29-may-2007	gjf	now in functional form


variable_patch_size = 0 ; 	% 1=yes, 0=no

% read the input fault geometry and set up the spatial models

[disloc_model, model1, model2, model3, model4, seg_coords, oksar_out] =...
  	proc_fdata_cntr_g_yagi(ingeomfile);

% set up the inputs for the calculation of G...

if (variable_patch_size==1) 

	inmodel2=model2;
	inmodel3=model3;

else 
	inmodel2=ones(size(model2));
	inmodel3=ones(size(model3));

end

% and do it!

if (n_diff==1)

	Gs=calc_smoothing_1stdiff_space(model1, inmodel2, inmodel3, model4, n_tstp) ;	% haven't modified this yet

else

	Gs=calc_smoothing_2nddiff_space(model1, inmodel2, inmodel3, model4, n_tstp) ;

end


spatial_model=model1;
