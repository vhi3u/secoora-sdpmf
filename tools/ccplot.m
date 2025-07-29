function ccHndl = ccplot(x,y,z,c_axis,symbol,marker_size)
% ccplot        Creates a color-coded plot.
%
%       CCPLOT(X,Y,Z,C_AXIS,SYMBOL,MARKER_SIZE) creates a color-coded plot of
%               vector Z with respect to vectors X and Y using C_AXIS to map
%               values of Z to the colormap of the current figure. Data is
%               plotted using SYMBOL's of MARKER_SIZE.
%
%               If C_AXIS is an empty vector CCPLOT calculates C_AXIS to be
%               [min(z) max(z)].
%
%       example...
%                       x = [0:.1:100];               
%                       y = sin(x);                   
%                       z = rand(1,length(y)) + 10.^y;
%                       subplot(2,1,1);plot3(x,y,z)                  
%                       subplot(2,1,2);ccplot(x,y,z,[],'.',20);      
%
  
  % Author:       Trevor Cooper, Marine Life Research Group/SIO
  %               tcooper@ucsd.edu
  %               December 8, 1995.
  
  % Modified by Sara Haines to return handles to objects created by this
  %  current call to this function, so that attributes can be set or changed
  % Modified hold behaviour to return it to the state it was before ccplot
  % called.
  % Modified to change marker edge and face attributes for color
  
  ccHndl = [];
  prev_hold = ishold;
  
  if(length(c_axis) == 0)
    isfin = find(finite(z));
    c_axis = [min(z(isfin)) max(z(isfin))];
  end
  
   cmap = get(gcf,'colormap');
%  colormap('default');
%  cmap=colormap;
  
  index = floor( (z-c_axis(1)) / ((c_axis(2)-c_axis(1)) / size(cmap,1)));
  
  too_small = find(index<1);
  if(length(too_small) > 0)
    index(too_small)=ones(length(too_small),1);
  end
  
  too_big = find(index>size(cmap,1));
  if(length(too_big) > 0)
    index(too_big)=ones(length(too_big),1)*size(cmap,1);
  end
  
  for j=1:size(cmap,1)
    matched=find(index == j);
    if(length(matched) > 0)
      ccHndl = [ccHndl plot(x(matched), y(matched), symbol, ...
			    'markeredgecolor', cmap(j,:), ...
			    'markerfacecolor', cmap(j,:), ...			    
			    'markersize', marker_size, 'tag', 'ccplot')];
      hold on;
    end
  end

  % return hold state to previous
  if prev_hold
    hold on;
  else
    hold off;
  end
    