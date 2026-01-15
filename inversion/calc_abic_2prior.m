function [abic,astar,s_astar,covar]=calc_abic_2prior(h,gs,gt,d,loge,einv,alphasq,betasq)

% function [abic,astar,s_astar,covar]=calc_abic_2prior(h,gs,gt,d,e,einv,alphasq,betasq)
%
% calculates ABIC and related important values [a* and s(a*)]
% given a kernel matrix, h, data vector, d, prior information matrices, gs
% (spatial smoothing) and gt (temporal smoothing),inverse covariance matrix, 
% einv, and trial hyperparameters, alphasq and betasq
%
% this here is the 'proper' formulation of ABIC, as expounded in Fukahata et al.
% (2003, 2004) - this gives better estimates of slip distributions that are less
% likely to be over-simplified than the 'improper' formulation used by Yoshida
% and others
%
% this is loosely based on code written by Yukitoshi Fukahata
% (see also Fukahata et al., 2003, Geophys. Res. Lett, 30, no. 6, 1305,
% and Fukahata et al., 2004, Geophys. J. Int., 156, 140-153)
%
% any mistakes are bound to be mine, in other words
%
% gjf, 01-nov-2004
% geniusinaction.com
%
% modifications:
%
%	15-dec-2004	gjf	optimised with sparse functions

% ----------------------
% preliminaries - calculate some necessary things
% ----------------------

% count number of values in alphasq

  [n_loopsa,junk]=size(alphasq);
  [n_loopsb,junk]=size(betasq);


% calculate N and M - numbers of datapoints and model parameters

  [N,M]=size(h);

%            T -1
% calculate H E  H

  hth=h'*einv*h;

%           T -1
% and also H E  d

  htd=h'*einv*d;



% ----------------------
% loop through values of alphasq and betasq...
% ----------------------

   for i=1:n_loopsa

      for j=1:n_loopsb


% ----------------------
% calculate the sum of the smoothing matrices, scaled by the hyperparameters, and its determinant
% ----------------------

         ags_bgt=(alphasq(i)*gs)+(betasq(j)*gt);
         p=rank(ags_bgt);

 %	 [p,ps]=size(ags_bgt);

         [u,s,v]=svd(ags_bgt,0);

%         logdeter1=sum(log(abs(diag(s))));
         logdeter1=sum(log(abs(diag(s(1:p,1:p)))));  % determinant of a square matrix is 
					      % equal to the product
                                              % of all of its non-zero singular values
                                              % -- this term is the log of that product

% ----------------------          T -1         2         2    -1  T -1
% ...and hence calculate a* ( = [H E  H + alpha G  + beta G  ]   H E  d
%                                                s         t 
% ----------------------

         hth_ags_bgt=hth+ags_bgt;
         pp=rank(hth_ags_bgt);

         [U,S,V]=svd(hth_ags_bgt,0);
         covar=inv(hth_ags_bgt);

         astar=fnnls(hth_ags_bgt,htd);
%          astar=covarhtd;

% ----------------------
% and then s(a*)
% ----------------------
            
% calculate d - Ha*

         d_ha=(d-(h*astar));

% and plug that into this equation for s(a*)

         s_astar(i,j)=(d_ha'*einv*d_ha)+(astar'*ags_bgt*astar);


% ----------------------
% and finally, abic...
% ----------------------

% now for the crunch...

% bung it all into one bundle to calculate ABIC

         logdeter2=sum(log(abs(diag(S))));  % determinant of a square matrix is equal to the product
                                              % of all of its non-zero singular values
                                              % -- this term is the log of that product


	 nlogsastar=N*log(s_astar(i,j));

         abic(i,j)=nlogsastar-logdeter1+logdeter2+loge;

         disp(sprintf('%0.5g %0.5g %d %d %f %f %f %f %f %f %f %1f %f',alphasq(i), betasq(j),p, pp, min(astar), max(astar), s_astar(i,j), log(s_astar(i,j)), nlogsastar, logdeter1, logdeter2, loge, abic(i,j)));                 

         outrep=sprintf('echo %0.5g %0.5g %d %d %f %f %f %f %f %f %f %f %f >> abic.out\n',alphasq(i), betasq(j),p, pp, min(astar), max(astar), s_astar(i,j), log(s_astar(i,j)), nlogsastar, logdeter1, logdeter2, loge, abic(i,j));
          
          system(outrep);


      end

   end

