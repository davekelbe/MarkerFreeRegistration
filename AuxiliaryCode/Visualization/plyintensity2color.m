function [ color ] = plyintensity2color( I, cmap )
%PLYINTENSITY2COLOR Returns a colormap 'color' from a vector 'I' for various cmaps  
%   
%   Example: 
%   I = (1:10)';    
%   color = plyintensity2color(I,'jet')
%   figure;
%   scatter(I,I,100,color,'filled');
%
%   (c) David Kelbe, Rochester Institute of Technology 
%   December 2013

if size(I,1)==1;
    I = I';
end

switch cmap
    case 'jet'
        c = jet(256);
    case 'HSV'
        c = hsv(256);
    case 'hsv'
        c = hsv(256);
    case 'hot'
        c = hot(256);
    case 'cool'
        c = cool(256);
    case 'spring'
        c = spring(256);
    case 'summer'
        c = summer(256);
    case 'autumn'
        c = autumn(256);
    case 'winter'
        c = winter(256);
    case 'gray'
        c = gray(256);
    case 'bone'
        c = bone(256);
    case 'copper'
        c = copper(256);
    case 'pink'
        c = pink(256);
    case 'BuDRd_12'
        c = colormap(othercolor('BuDRd_12',256));
    otherwise
        error('Must supply valid colormap string');
end

R = 0*I;
G = 0*I;
B = 0*I;

c = uint8(round(255*c));

J = I-min(I(:));
J = round(255*I/max(I));
for i=1:255;
    ix = find(J==i);
    R(ix)=c(i,1);
    G(ix)=c(i,2);
    B(ix)=c(i,3);
end

color = [R,G,B];

