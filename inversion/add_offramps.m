function [kernel,obs]=add_offramps(raw_kernel, raw_obs, ndsets, ndata, soff, sxg, syg, offs, xgs, ygs, locs, xref, yref)

% function [kernel,obs]=add_offramps(raw_kernel, raw_obs, ndsets, ndata, soff, sxg, syg, offs, xgs, ygs, locs, xref, yref)
%
% function that adds optional extras to the kernel, such as the ability to
% solve for los offsets or los gradients (ramps) in the x and y directions
% 
% means you only need to calculate the kernel once
%
% gjf, 16-sep-2004
% geniusinaction.com
%
% added gradient solving capability, 27-oct-2004
% fixed addition of a priori offsets+gradients 09-mar-2005

% --------
% Work out how big the output will be, and define it
% --------

[nrows,ncols]=size(raw_kernel);
nsoff=sum(soff==(1:ndsets));
nsxg=sum(sxg==(1:ndsets));
nsyg=sum(syg==(1:ndsets));
extracols=nsoff+nsxg+nsyg;

obs=zeros(nrows,1);
kernel=zeros(nrows,ncols+extracols);
kernel(1:nrows,1:ncols)=raw_kernel;

% ---------
% Calculate the x and y distances to the reference observation
% ---------

x_dist=(locs(:,1)-(xref))*1e-6;
y_dist=(locs(:,2)-(yref))*1e-6;

% ---------
% Add the a priori offsets and gradients, if desired
% ---------

% add any offsets desired to the data

linecount=1;

for i=1:ndsets

%    disp(['dataset ' num2str(i) ' offset ' num2str(offs(i)) ' start ' num2str(linecount) ' end ' num2str(linecount+ndata(i)-1) ]);
    obs(linecount:linecount+ndata(i)-1)=raw_obs(linecount:linecount+ndata(i)-1)+offs(i);

%    for j=linecount:linecount+9
%       disp([raw_obs(j) obs(j) offs(i)])
%    end

    linecount=linecount+ndata(i);

end

% add any x gradients desired to the data

linecount=1;

for i=1:ndsets
   
    obs(linecount:linecount+ndata(i)-1)=obs(linecount:linecount+ndata(i)-1)+...
          x_dist(linecount:linecount+ndata(i)-1)*xgs(i);
    linecount=linecount+ndata(i);

end

% add any y gradients desired to the data

linecount=1;

for i=1:ndsets

    obs(linecount:linecount+ndata(i)-1)=obs(linecount:linecount+ndata(i)-1)+...
          y_dist(linecount:linecount+ndata(i)-1)*ygs(i);
    linecount=linecount+ndata(i);

end

% ---------
% Add offset and gradient-solving terms to the kernel, if desired
% ---------

% if required, add the bit that solves for offsets


colcount=1;
linecount=1;
offcount=0;


    for j=1:ndsets

        if(soff(j)~=0)

            if(soff(j)==j)

                kernel(linecount:linecount+ndata(j)-1,ncols+colcount)=-1;               
                linecount=linecount+ndata(j);
                colcount=colcount+1;
                offcount=offcount+1;

             else
   
                kernel(linecount:linecount+ndata(j)-1,ncols+soff(j))=-1;               
                linecount=linecount+ndata(j);

             end

        end

    end


% and now for the x gradients...

linecount=1;
xgcount=0;

    for j=1:ndsets

        if(sxg(j)~=0)
 
            if(sxg(j)==j)

               xgcount=xgcount+1;
               kernel(linecount:linecount+ndata(j)-1,ncols+colcount)=x_dist(linecount:linecount+ndata(j)-1)*-1; 
               colcount=colcount+1;
               linecount=linecount+ndata(j);

            else

               kernel(linecount:linecount+ndata(j)-1,ncols+offcount+sxg(j))=...
                  x_dist(linecount:linecount+ndata(j)-1)*-1; 
               linecount=linecount+ndata(j);

            end

        end     

    end


% and the y...

linecount=1;

    for j=1:ndsets

        if(syg(j)~=0)

            if(syg(j)==j)

               kernel(linecount:linecount+ndata(j)-1,ncols+colcount)=y_dist(linecount:linecount+ndata(j)-1)*-1;  
               colcount=colcount+1;
               linecount=linecount+ndata(j);

            else

               kernel(linecount:linecount+ndata(j)-1,ncols+offcount+xgcount+syg(j))=...
                  y_dist(linecount:linecount+ndata(j)-1)*-1;  
               linecount=linecount+ndata(j);

            end

        end 

    end
