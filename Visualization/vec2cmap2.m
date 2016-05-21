function [ color ] = vec2cmap2( I, cmap, varargin )
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
%I = data_ieq;
if nargin == 4;
    Imin = varargin{1};
    Imax = varargin{2};
elseif nargin ==2;
    Imin = min(I(:));
    Imax = max(I(:));
end

if size(I,1)==1;
    I = I';
end
    
switch cmap
    case 'jet'
        c = colormap(jet(256)); 
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

delete(gcf)
%a = 1;
%c  = (c + a)./(1+a); %Changed in v2
R = 0*I;
G = 0*I;
B = 0*I;

c = uint8(round(255*c));

Iz = (I-Imin)/(Imax-Imin);
%Iz(Iz>1) = 1;
%Iz(Iz<0) = 0;
J = round(255*Iz) + 1;
for i=1:256;
    ix = find(J==i);
    R(ix)=c(i,1);
    G(ix)=c(i,2);
    B(ix)=c(i,3);
end

color = uint8([R,G,B]);

