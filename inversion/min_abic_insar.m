% min_abic_insar.m  - a program for the computation of
% akaike's bayesian information criterion, (abic) a statistical
% wotsit that allows you to estimate objectively the degree
% of smoothing necessary in a distributed slip inversion.
%
% inputs are as follows:
%
% h_infile     - ascii file containing the h (kernel) matrix
% d_infile     - ascii file containing the d (data) vector
% sta_infile   - ascii file containing station information
% einv_infile  - ascii file containing the einv (inverse covariance) matrix
%
% this is the state of the art; robert muir wood would be agog
%
% based on min_ABIC2.f by Yukitoshi Fukahata and someone called 'Honsho'
% (any new errors are purely my fault)
%
% gjf, 13-oct-2004
% geniusinaction.com
%
% change history:
%
%   28-may-2007   gjf   adapted for variable trace length and errors

clear;
!\rm abic.out
 
% ----------------------
% set up inputs
% ----------------------

h_insr_infile='h_insr.dat';
d_insr_infile='d_insr.dat';


geom_infile='manyi_planar_6km.dat';

n_as=30;
n_dd=3;
n_tstp=1;
n_diff=2;

load locs.dat;
load ndata.dat;

load einv_insr.mat;
load e_insr.mat;
log_det_e_insr=-3.5497e4;

% ----------------------
% nuisance parameters here:
% ----------------------

solve_offset=[1 2 3 4]; % solve for a LOS offset for each dataset? (0=no)
solve_x_grad=[1 2 3 3]; % solve for a LOS x gradient for each dataset? (0=no)
solve_y_grad=[1 2 3 3]; % solve for a LOS y gradient for each dataset? (0=no)
offset=[-0.3 -0.3 -0.3 -0.4];      % static offset to add to the data (for each dataset)
x_grad=[-1 -1 0 0];      % x-gradient to add to the data (ppm)
y_grad=[-1 -1 -3 -3];      % y-gradient to add to the data (ppm)


% ----------------------
% preliminaries - set up some variables, load inputs
% ----------------------

bignum=1e20;            % a suitably big number 

h_insr=load(h_insr_infile);
d_insr=load(d_insr_infile);


[N_insr,M_insr]=size(h_insr);


% ----------------------
% assemble the kernel and smoothings next
% ----------------------

h_top=[];

for i=1:n_tstp

   h_top=[h_top h_insr];

end

[h_top_nuis,d_top_nuis]=add_offramps(h_top,d_insr,4,ndata,solve_offset,solve_x_grad,solve_y_grad,offset,x_grad,y_grad,locs,497964,3894725);

[woo,M_insr_nuis]=size(h_top_nuis);

H=h_top_nuis;
d=d_top_nuis;

M_nuis=M_insr_nuis-M_insr;

Gs=zeros(M_insr_nuis);

[Gs(1:M_insr,1:M_insr),spatial_model]=calc_gs_yagi(geom_infile,n_tstp,n_diff);
Gt=zeros(size(Gs));

% ----------------------
% bodge it here to see why it isn't working...
% ----------------------

%einv=speye(N);
%g2=speye(M);
%G=g'*g;

% ----------------------
% set up three arrays with various values of alphasq, betasq and gammasq
% ----------------------

alphasq=[1e-27 1e-21 1e-15 1e-12 1e-9 1e-6 1e-3 1e-2 3e-2 5e-2 7e-2 1e-1 1.5e-1 2e-1 2.5e-1 3e-1 3.5e-1 4e-1 4.5e-1 5e-1 7e-1 1e0 2e0 3e0 5e0 7e0 1e1 3e1 5e1 7e1 1e2 3e2 5e2 7e2 1e3 3e3 5e3 7e3 1e4 3e4 5e4 7e4 1e5 1e6 1e7 1e8 1e9];
%betasq=[10 20 30 40 50];
%gammasq=[0.01 0.02 0.03 0.04 0.05];


%alphasq=0.7;
betasq=0;
gammasq=1;

% octave and matlab have a different convention for vectors - we
% want a column here, so let's make it so

[rws,cls]=size(alphasq);

if (rws<cls)

   alphasq=alphasq';
   betasq=betasq';
   gammasq=gammasq';

end

[n_loops,woo]=size(gammasq)


% ----------------------
% loop through gammasq values; assemble weighted inverse-covariance matrix
% ----------------------

for i=1:n_loops

   disp(sprintf('Assembling inverse covariance/weight matrix, gamma^2 = %f',gammasq(i)));

   Einv(1:N_insr,1:N_insr)=einv_insr;

%   log_det_e_gammasq=2*N_seis*log10(gammasq(i));
   log_det_e_gammasq=0;
  log_det_e=log_det_e_insr+log_det_e_gammasq;

% ----------------------
% and call calc_abic to calculate the values of abic
% ----------------------

[abic,astar,s_astar,covar]=calc_abic_2prior(H,Gs,Gt,d,log_det_e,Einv,alphasq,betasq);

        move_output=sprintf('mv abic.out abic%f_insar.dat',gammasq(i));
        system(move_output);

end 

if ((rws==1)&(cls==1))

   slip=zeros(90,1);
   sum_covar=zeros(90);

   for i=1:1

      slip=slip+astar(((i-1)*90)+1:i*90);
      sum_covar=sum_covar+covar(((i-1)*90)+1:i*90,((i-1)*90)+1:i*90);

   end

imagesc(reshape(slip,30,3)'); colorbar; axis image; axis xy;


end
