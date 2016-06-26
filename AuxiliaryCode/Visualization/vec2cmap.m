function [ color2 ] = vec2cmap( I, cmap, varargin )
%VEC2CMAP Returns a colormap 'color' from a vector 'I' for various cmaps  
%   Optionallly specify minimum and maximum value of I for which color saturates   
%   
%   Example: 
%   I = (1:10)';    
%   color = vec2cmap(I,'jet')
%   
%   Example:
%   color = vec2cmap(I,'jet', 3,9)
%
%   (c) David Kelbe, Rochester Institute of Technology 
%   December 2013

%load('D:\Users\djk2312\Documents\Harvard\master\031\17\mat\data_ieq.mat');
%I = data_ieq;'
if size(I,1)==1;
    I = I';
end
n_el = numel(I);
isvalid = ~isinf(I);
I = I(isvalid);

if nargin == 4;
    Imin = varargin{1};
    Imax = varargin{2};
    I(I>Imax) = Imax;
    I(I<Imin) = Imin;
    Iz = (I-Imin)/(Imax-Imin);
%Iz(Iz>1) = 1;
%Iz(Iz<0) = 0;
elseif nargin ==2;
    Imin = min(I(:));
    Imax = max(I(:)); 
Iz = (I-Imin)/(Imax-Imin);
end


    
switch cmap
    case 'jet'
        c = jet(256); 
    case 'HSV'
        c = colormap(hsv(256));
    case 'hsv'
        c = colormap(hsv(256));
    case 'hot'
        c = colormap(hot(256));
    case 'cool'
        c = colormap(cool(256));
    case 'spring'
        c = colormap(spring(256));
    case 'summer'
        c = colormap(summer(256));
    case 'autumn'
        c = colormap(autumn(256));
    case 'winter'
        c = colormap(winter(256));
    case 'gray'
        c = colormap(gray(256));
    case 'bone'
        c = colormap(bone(256));
    case 'copper'
        c = colormap(copper(256));
    case 'pink'
        c = colormap(pink(256));
    case 'BrBG4'
        c = colormap(othercolor('BrBG4',256));
    otherwise
        error('Must supply valid colormap string');
end

%delete(gcf)
%a = 1;
%c  = (c + a)./(1+a);
R = 0*I;
G = 0*I;
B = 0*I;

c = uint8(round(255*c));

J = round(255*Iz);
for i=0:255;
    ix = find(J==i);
    R(ix)=c(i+1,1);
    G(ix)=c(i+1,2);
    B(ix)=c(i+1,3);
end

color = uint8([R,G,B]);
color2 = zeros(n_el,3);
color2(isvalid,1) = R;
color2(isvalid,2) = G;
color2(isvalid,3) = B;
color2 = uint8(color2);


%color2(isvalid,:) = color

