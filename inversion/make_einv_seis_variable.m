function [einv, log_det_e] = make_einv_seis_variable_20130114(stasamp,stasigw)

% function [einv, log_det_e] = make_einv_seis_variable(stasamp,stasigw)
%
% calculates the inverse covariance einv and the relevant log-determinant
% to calculate ABIC harmoniously
%
% stasamp - number of time samples per station component
% stasigw - uncertainty in station records
%
% gjf, 29-may-2007
% genius in action (tm)

[n_sta,woo]=size(stasamp);
n_data=sum(stasamp);

errs=zeros(n_data,1);
cnt=0;

for i=1:n_sta

   errs(cnt+1:cnt+stasamp(i))=ones(stasamp(i),1)/stasigw(i)/stasigw(i);
   cnt=cnt+stasamp(i);

end

einv=spdiags(errs,0,n_data,n_data);
log_det_e=2*sum(double(stasamp).*log(stasigw));
