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



% ----------------------
% preliminaries - set up some variables, load inputs
% ----------------------

bignum=1e20;            % a suitably big number 

[staname,stasigw,stasamp,statstp,stalocn]=read_stations(sta_infile);
[nsta,jnk]=size(stasigw);
[h_seis,d_seis]=read_yagi_data(h_seis_infile,d_seis_infile);

[N_seis,M_seis]=size(h_seis);

[einv_seis, log_det_e_seis] = make_einv_seis_variable(stasamp,stasigw);


% ----------------------
% assemble the kernel and smoothings next
% ----------------------

H=h_seis;
d=d_seis;

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

%alphasq=[1e-5 1e-4 1e-3 3e-3 5e-3 7e-3 1e-2 2e-2 3e-2 4e-2 5e-2 7e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5];
%betasq=[1e-5 1e-4 1e-3 1e-2 1e-1 1e0 5e0 6e0 6.5e0 7e0 7.5 8e0 9e0 1e1 2e1 3e1 5e1 7e1 1e2 1e3 1e4 1e5];
gammasq=[1];


alphasq=0.03;
betasq=6.5;
%gammasq=0.03;

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

   Einv=einv_seis;

   log_det_e_gammasq=2*N_seis*log(gammasq(i));
   log_det_e=log_det_e_seis+log_det_e_gammasq;

% ----------------------
% and call calc_abic to calculate the values of abic
% ----------------------

[abic,astar,s_astar,covar]=calc_abic_2prior(H,Gs,Gt,d,log_det_e,Einv,alphasq,betasq);



end 

if ((rws==1)&(cls==1))

   slip=zeros(90,1);
   sum_covar=zeros(90);

   for i=1:6

      slip=slip+astar(((i-1)*90)+1:i*90);
      sum_covar=sum_covar+covar(((i-1)*90)+1:i*90,((i-1)*90)+1:i*90);

   end

   imagesc(reshape(slip,30,3)'); colorbar; axis image; axis xy;


   WRSS_seis = (d-H*astar)'*Einv*(d-H*astar)
   NM_seis = (d-H*astar)'*(d-H*astar)/(d'*d)

   synthetic_waveforms = H*astar;

   cnt=1;

% output the waveforms

   for i=1:nsta

      cntend=cnt+stasamp(i)-1;
      waveform=synthetic_waveforms(cnt:cntend);

      outrename=sprintf('mv output_wave.dat %s_seis.dat',staname{i,1});
   
      timesamples=(1:stasamp(i))'-1;

      timestamp=double(timesamples)*statstp(i);
      
      output_wave=[timestamp waveform];

      save output_wave.dat output_wave -ascii;

      system(outrename);

      cnt=cnt+stasamp(i);

   end

else

        move_output=sprintf('mv abic.out abic%f_seis.dat',gammasq(i));
        system(move_output);

end

