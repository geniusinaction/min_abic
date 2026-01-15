
% calc_e.m  calculate the 'E' matrix required by ABIC
%
% calculates a covariance matrix for multiple independent displacement 
% datasets using a given exponential covariance relation:
%
% e = maxvar^2*exp(-alpha*r)
%
% where 
%
% e		= covariance
%
% r		= separation between points
% maxvar	= maximum covariance (autocorrelation function)
%		  in radians/km
% alpha		= exponential coefficient calculated by cvdcalc
%  
% this is a bastardised version of 'pcmc.m' by tim wright
%
% he won't mind
%
%
% gjf, 11-oct-2004
% geniusinaction.com
%

% -------------------
% Set up the inputs and outputs
% -------------------

in_locs = './locs.dat'	% 3 column matrix containing (x,y,z) locations 
			% of observations (in metres)

n_data=load('ndata.dat');
[jnk1,jnk2]=size(n_data);
n_datasets=jnk1*jnk2;


maxvar = [4.82e-5 4.69e-5 4.81e-5 4.81e-5];		 	% variances in m^2
alpha =  [0.134   0.199   0.085   0.085];		% exponent coefficients in km^-1

%maxvar = (0.028333/(2*pi))^2*maxvar ; 		% convert variance from radians^2/km 
                   				% to metres^2/km

% -------------------
% read in the data and locations
% -------------------

fin=fopen(in_locs);
locs=fscanf(fin,'%lf',[3,inf]);
fclose(fin);

locs = locs';
x = locs(:,1)/1000;
y = locs(:,2)/1000;


% -------------------
% set up the output; start the loop
% -------------------

total_data=sum(n_data);
datacounter=1;

e_ins=zeros(sum(n_data));
einv_ins=zeros(sum(n_data));

log_det_e=0;

for i=1:n_datasets

	x_sub=x(datacounter:n_data(i)+datacounter-1);
	y_sub=y(datacounter:n_data(i)+datacounter-1);

% -------------------
% calculate distance between points in a matrix
% -------------------

	[xx1,xx2] = meshgrid(x_sub,x_sub); % matrices containing coordinates of points
	[yy1,yy2] = meshgrid(y_sub,y_sub);
	rgrid = sqrt((xx1-xx2).^2+(yy1-yy2).^2);
	clear xx1 xx2 yy1 yy2;

% -------------------
% Use expcos function to calculate vcm
% -------------------

	vcm = maxvar(i)*exp(-alpha(i)*rgrid);
	invvcm=inv(vcm);
        [u,s,v]=svd(vcm,0);
        log_det_e=log_det_e+sum(log10(abs(diag(s))));
	clear rgrid;

%disp('willy');
%pause;
	
% -------------------
% Copy this vcm and its inverse to the E and E_inv matrices
% -------------------

	e_ins(datacounter:n_data(i)+datacounter-1,datacounter:n_data(i)+datacounter-1)=vcm;
	einv_ins(datacounter:n_data(i)+datacounter-1,datacounter:n_data(i)+datacounter-1)=invvcm;

	clear vcm;

	datacounter=datacounter+n_data(i);

end

%E_inv=inv(E);

log_det_e

save e_ins.mat e_ins;
save einv_ins.mat einv_ins
save log_det_e.dat log_det_e -ascii
