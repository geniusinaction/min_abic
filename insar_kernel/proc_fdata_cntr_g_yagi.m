function [disloc_model, spatial_model1, spatial_model2, spatial_model3, spatial_model4,...
     seg_coords, oksar_out] =  proc_fdata_cntr_g_yagi(input_filename)

% function [disloc_model, spatial_model1, spatial_model2, spatial_model3,...
%    spatial_model4,fault_coords, faultseg_coords, oksar_out] =...
%    proc_fdata_cntr_g_yagi(input_filename)
%
% This is a function which formats fault model data, from the input
% file specified, into a form that the Okada function 'disloc3d'
% can read. It also generates a spatial model which relates the 
% various fault patches spatially, and a dataset that will be used 
% to generate 'oksar'-friendly output later on.
%
% --
%
% this particular version uses the numbering scheme of Fukahata, Yagi and other
% luminaries, which is different to the one i traditionally used
%
% modified at least five times from the original classic 'process_faultdata'
% written by me in 2001
%
% gjf, 31-oct-2004 (wooooo)
% geniusinaction.com


% input data file
% disp(input_filename)
  fin = fopen(input_filename,'r');
  fault_input = fscanf(fin,'%lf',[11,inf]);
  fault_input = fault_input';
  fclose(fin);

% set up variables

strike=fault_input(:,1) ;
dip=fault_input(:,2) ;
rake=fault_input(:,3) ;
xc=fault_input(:,4)*1000 ;
yc=fault_input(:,5)*1000 ;
length=fault_input(:,6) ;
top=fault_input(:,7) ;
bottom=fault_input(:,8) ;
alongstrike=fault_input(:,9) ;
downdip=fault_input(:,10) ;
affil=fault_input(:,11);

%seg_coords = fault_input(:,1:4);

deg2rad=pi/180 ;

% calculate number of model parameters

[m,n] = size(xc);  %number of fault 'segments'

% work out dimensions of the spatial models, and set 'em up

rws=max(downdip);
clms=sum(alongstrike);
spatial_model1=zeros(rws,clms);
spatial_model2=zeros(rws,clms);
spatial_model3=zeros(rws,clms);
spatial_model4=ones(rws,clms)*-1;


% set count variables to zero

l=0 ;
p=0 ;

% set reversal flag to zero

reversed=0 ;

% set up loops in which to calculate the desired parameters
% (e.g. strike, dip, length, width...) and put them into the
% rows of a matrix 'disloc_model' which is in the correct
% format for the Okada routine 'disloc3d'


% loop through entries in the input data

for i=1:m ;

% check the dip of the fault plane... if it is greater than 90
% degrees, take appropriate action - switching end point
% coordinates, adjusting the dip, and reversing the direction of
% strike-slip (this is to allow reversals in dip along strike)

    if (dip(i)>90)
 
% reversal flag        
        
        reversed=1 ;
        
% reverse strike

        strike(i)=strike(i)+180 ;
        if (strike(i)>=360)
		strike(i)=strike(i)-360;
	end
        
% dip adjustment:

       dip(i)=180-dip(i) ;   
       
% rake reversal:        
        
       rake(i)=-1*rake(i) ;
    
    elseif (dip(i)<=90)

        
% reversal flag        
        
        reversed=0 ;
        
    end
    
% calculate fault patch lengths

    kmlength(i)=length(i)/alongstrike(i) ;
    length(i)=kmlength(i)*1000 ;
    
% calculate fault patch widths

    width(i)=(faultwidth(dip(i), top(i), bottom(i)))/downdip(i); 
    vertwidth(i) = (bottom(i)-top(i))*1000/downdip(i);

% calculate strike-slip and dip-slip components
    
    rrake(i)=(rake(i)+90)*deg2rad ;
    dipslip(i)=cos(rrake(i)) ;
    strslip(i)=sin(rrake(i)) ;
    
% calculate the x and y increments

    [xinctmp,yinctmp]=strike2inc(strike(i));
    xinc(i)=xinctmp*length(i);
    yinc(i)=yinctmp*length(i);

% calculate the fault segment end coordinates
    
    seg_coords(i,1)=xc(i)-(xinctmp*length(i)*alongstrike(i)/2);
    seg_coords(i,2)=yc(i)-(yinctmp*length(i)*alongstrike(i)/2);
    seg_coords(i,3)=xc(i)+(xinctmp*length(i)*alongstrike(i)/2);   
    seg_coords(i,4)=yc(i)+(yinctmp*length(i)*alongstrike(i)/2);


% set down-dip counter to zero

    ddc=0;

% loop through down-dip subdivisions        
        
        for k=downdip(i):-1:1
           

           ddc=ddc+1;

% loop through along-strike subdivisions   


% if reversed, loop through coordinates backwards

    if (reversed==1)
     startj = alongstrike(i)-(alongstrike(i)+1)/2;
     incj   = -1;
     endj   = 1-(alongstrike(i)+1)/2;
    elseif (reversed==0)
     startj = 1-(alongstrike(i)+1)/2;
     incj   = 1;
     endj   = alongstrike(i)-(alongstrike(i)+1)/2;
    end


% reset along-strike counter

     p=0;

    for j=startj:incj:endj 

% increment along-strike division counter

        p=p+1 ;
        
% increment fault element counter

        l=l+1;

% calculate end coordinates for surface trace of the subdivided
% fault patches:       
        
%        tempx1(i,j)=((x2(i)-x1(i))*(j-1)/alongstrike(i))+x1(i) ;
%        tempy1(i,j)=((y2(i)-y1(i))*(j-1)/alongstrike(i))+y1(i) ;
%        tempx2(i,j)=((x2(i)-x1(i))*j/alongstrike(i))+x1(i) ;
%        tempy2(i,j)=((y2(i)-y1(i))*j/alongstrike(i))+y1(i) ;
        
%       tempxcentre(i,j)=tempx1(i,j)+(tempx2(i,j)-tempx1(i,j))/2 ;
%        tempycentre(i,j)=tempy1(i,j)+(tempy2(i,j)-tempy1(i,j))/2 ;

        tempxcentre=xc(i)+j*xinc(i);
        tempycentre=yc(i)+j*yinc(i) ;
              
            
% disloc_model output:            
 

% x location (centre of projection of fault plane)
            disloc_model(1,l)=tempxcentre;

% y location (centre of projection of fault plane)
            disloc_model(2,l)=tempycentre;

% strike
            disloc_model(3,l)=strike(i);
           
% dip 
            disloc_model(4,l)=dip(i);

% rake
            disloc_model(5,l)=rake(i);

% slip
            disloc_model(6,l)=1;

% length
            disloc_model(7,l)=length(i);

% hmin      
            disloc_model(8,l)=top(i)+(k-1)*vertwidth(i);
% hmax
            disloc_model(9,l)=disloc_model(8,l)+vertwidth(i);

% spatial_model output: % NB octave only deals with 2D arrays not 3D           
            
% set fault element number (layer 1)

            spatial_model1(ddc,p)=l ;
            
% set fault element length (layer 2)            
            
            spatial_model2(ddc,p)=length(i) ;
            
% set fault element width (layer 3)     
            
            spatial_model3(ddc,p)=width(i);

% set smoothing 'affiliation' (layer 4)

            spatial_model4(ddc,p)=affil(i);            
            
% oksar_out output:

% set x coordinate (row 1)

            oksar_out(1,l)=tempxcentre ;

% set y coordinate (row 2)            
            
            oksar_out(2,l)=tempycentre;
         
% set strike (row 3)            
            
            oksar_out(3,l)=strike(i) ;
            
% set dip (row 4)
            
            oksar_out(4,l)=dip(i) ;
            
% set rake (row 5) 
            
            oksar_out(5,l)=rake(i) ;
            
% set length (row 6)

            oksar_out(6,l)=kmlength(i) ;
            
% set top (row 7)
            
            oksar_out(7,l)=top(i)+((k-1)*...
                (bottom(i)-top(i))/downdip(i)) ;
 
% set bottom (row 8)
            
            oksar_out(8,l)=bottom(i)-((downdip(i)-k)*...
                (bottom(i)-top(i))/downdip(i)) ;
            
        end
    end
end

    
