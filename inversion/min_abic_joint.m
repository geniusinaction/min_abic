% min_abic_seis.m  - a program for the computation of
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


h_seis_infile='A.matrix';
d_seis_infile='B.vector';
sta_infile='stations.txt';
geom_infile='manyi_planar_6km.dat';

n_as=30;
n_dd=3;
n_tstp=6;
n_diff=2;

load locs.dat;
load ndata.dat;

load einv_insr.mat;
load e_insr.mat;
log_det_e_insr=-3.5497e4;   # only needs to be calculated once

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

[staname,stasigw,stasamp,statstp,stalocn]=read_stations(sta_infile);
[nsta,jnk]=size(stasigw);
[h_seis,d_seis]=read_yagi_data(h_seis_infile,d_seis_infile);

[N_insr,M_insr]=size(h_insr);
[N_seis,M_seis]=size(h_seis);

[einv_seis, log_det_e_seis] = make_einv_seis_variable(stasamp,stasigw);


% ----------------------
% assemble the kernel and smoothings next
% ----------------------

h_top=[];

for i=1:n_tstp

   h_top=[h_top h_insr];

end

[h_top_nuis,d_top_nuis]=add_offramps(h_top,d_insr,4,ndata,solve_offset,solve_x_grad,solve_y_grad,offset,x_grad,y_grad,locs,497964,3894725);

[woo,M_insr_nuis]=size(h_top_nuis);

H=zeros(N_insr+N_seis,M_insr_nuis);
H(1:N_insr,1:M_insr_nuis)=h_top_nuis;
H(N_insr+1:N_insr+N_seis,1:M_seis)=h_seis;

d=[d_top_nuis; d_seis];

M_nuis=M_insr_nuis-M_seis;

Gs=zeros(M_insr_nuis);
Gt=zeros(M_insr_nuis);

[Gs(1:M_seis,1:M_seis),spatial_model]=calc_gs_yagi(geom_infile,n_tstp,n_diff);
[Gt(1:M_seis,1:M_seis),temporal_model]=calc_gt_yagi(n_as,n_dd,n_tstp,n_diff);


% ----------------------
% bodge it here to see why it isn't working...
% ----------------------

%einv=speye(N);
%g2=speye(M);
%G=g'*g;

% ----------------------
% set up three arrays with various values of alphasq, betasq and gammasq
% ----------------------

%alphasq=[1e-5 1e-4 1e-3 1e-2 1e-1 0.5 0.9 1e0 1e1 1e2 1e3 1e4 1e5];
%betasq=[1e-5 1e-4 1e-3 1e-2 1e-1 1 1e1 35 1e2 1e3 1e4 1e5];
%gammasq=[0.01 0.02 0.025 0.03 0.035 0.04 0.05];
%gammasq=[0.001 0.005 0.01 0.025 0.05 0.1 1 5];
%gammasq=[0.5];

alphasq=0.9;
betasq=35;
gammasq=0.025;

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

   Einv=eye(N_insr+N_seis);
   Einv(1:N_insr,1:N_insr)=einv_insr;
   Einv(N_insr+1:N_insr+N_seis,N_insr+1:N_insr+N_seis)=einv_seis/gammasq(i);

   log_det_e_gammasq=2*N_seis*log(gammasq(i));
   log_det_e=log_det_e_insr+log_det_e_seis+log_det_e_gammasq;

% ----------------------
% and call calc_abic to calculate the values of abic
% ----------------------

[abic,astar,s_astar,covar]=calc_abic_2prior(H,Gs,Gt,d,log_det_e,Einv,alphasq,betasq);

        move_output=sprintf('mv abic.out abic%f_joint.dat',gammasq(i));
        system(move_output);

end 

if ((rws==1)&(cls==1))

   slip=zeros(90,1);
   sum_covar=zeros(90);

   for i=1:6

      slip=slip+astar(((i-1)*90)+1:i*90);
      sum_covar=sum_covar+covar(((i-1)*90)+1:i*90,((i-1)*90)+1:i*90);

   end

imagesc(reshape(slip,30,3)'); colorbar; axis image; axis xy;

WRSS_seis = (d_seis-h_seis*astar(1:M_seis))'*(d_seis-h_seis*astar(1:M_seis))
NM_seis = (d_seis-h_seis*astar(1:M_seis))'*(d_seis-h_seis*astar(1:M_seis))/(d_seis'*d_seis)

   synthetic_waveforms = h_seis*astar(1:M_seis);

   cnt=1;

% output the waveforms

   for i=1:nsta

      cntend=cnt+stasamp(i)-1;

      synthwaveform=synthetic_waveforms(cnt:cntend);
      outrename=sprintf('mv output_wave.dat %s_joint.dat',staname{i,1});
   
      timesamples=(1:stasamp(i))'-1;
      timestamp=double(timesamples)*statstp(i);
      
      output_wave=[timestamp synthwaveform];
      save output_wave.dat output_wave -ascii;
      system(outrename);

      datawaveform=d_seis(cnt:cntend);
      datrename=sprintf('mv output_wave.dat %s_data.dat',staname{i,1});

      output_wave=[timestamp datawaveform];
      save output_wave.dat output_wave -ascii;
      system(datrename);
      
      cnt=cnt+stasamp(i);

   end


end
