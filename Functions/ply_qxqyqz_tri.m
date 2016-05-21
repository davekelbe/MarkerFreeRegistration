function [ vertices, indices ] = ply_qxqyqz_tri( qx,qy,qz )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

nr = size(qx,1);
nc = size(qx,2);

vertices = zeros((nc-1)*(nr-1),3);
indices = zeros((nc-1)*(nr-1),3);
vertices2 = zeros((nc-1)*(nr-1),3);
indices2 = zeros((nc-1)*(nr-1),3);

i=1;
for r = 1:nr-1;
    for c = 1:nc-1;
        indices(i,:) = [nc*(r-1)+c nc*(r-1)+c+1 nc*(r-1)+c+nc ]; %nc*(r-1)+c+nc];
        i = i+1;
    end
end
 
i=1;
for r = 1:nr;
    for c = 1:nc;
        vertices(i,:) = [qx(r,c) qy(r,c) qz(r,c)];
        i = i+1;
    end
end

i=1;
for r = 1:nr-1;
    for c = 1:nc-1;
        indices2(i,:) = [ nc*(r-1)+c+1 nc*(r-1)+c+nc+1 nc*(r-1)+c+nc ];
        i = i+1;
    end
end
 
i=1;
for r = 1:nr;
    for c = 1:nc;
        vertices2(i,:) = [qx(r,c) qy(r,c) qz(r,c)];
        i = i+1;
    end
end

indices = [indices; indices2];
vertices = [vertices; vertices2];

end
