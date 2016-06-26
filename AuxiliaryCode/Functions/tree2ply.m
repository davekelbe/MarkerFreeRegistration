function tree2ply( file, tree, tess )
%CYL2PLY exports a tree model to a ply
%
% DESCRIPTION:
%   Converts an array of tree objects (locations, radii, and colors) to a
%   ply file
%
% INPUTS:
%   file (string) the path containing the ply file to write to.
%   tree( i ).loc ( 3 x m_i array )
%   tree( i ).r ( 1 x m_i array )
%   tree( i ).color ( 3 x m_i array )
%   clr [optional] ( 3 x 1 or 3 x n array ) the color (in range [0, 255])
%       for the cylinder. default is [ 127; 127; 127 ].
%   tess [optional] the number of tessellations to use for radius and 
%       height. default is 6.
%       tess : the tessellations for radius (must be >= 3)
%       tess( 2 ) : the tessellations for height (must be >= 1)
%
% OUTPUTS:
%   <NONE>
%
% HISTORY:
%   2013-05-27: written by Paul Romanczyk (par4249 at rit dot edu)

% check the tessellations
if nargin < 3
    tess = 6;
end

if numel( tess ) ~= 1
    error( 'bole2ply:InvalidTess', ...
        'Invalid tess: must have 1 element' );
end
if tess < 3
    error( 'bole2ply:InvalidTessR', ...
        'There must be at least 3 tessellations of radius' );
end

fid = fopen( file, 'w' );
if fid <= 0
    % we cannot write to standard input (fid=0), we can write to std output
    % (fid=1), std error (fid=2), or a file (fid>2)
    error( 'bole2ply:InvalidFile', [ '"' file '" is not a valid file.' ] );
end

n = numel( tree );
if n==1;
    n= 0;
end
% pre-flight check
cylsPerTree = zeros( size( tree ) );
for i = 1:n
    cylsPerTree( i ) = numel( tree( i ).r ) - 1;
    if cylsPerTree( i ) < 1
        error( 'bole2ply:InvalidTree', ...
            [ 'tree ' num2str( i ) ' has too few cylinders' ] );
    end
    
    % check colors
    sc = size( tree( i ).color );
    if numel( sc ) == 0
        tree( i ).color =  [ 127; 127; 127 ] * ...
            ones( 1, cylsPerTree( i ) + 1 );
    elseif numel( sc ) ~= 2
        error( 'bole2ply:InvalidColor', ...
            [ 'tree ' num2str( i ) ' has an invalid color' ] );
    end
    if sc( 1 ) ~= 3
        error( 'bole2ply:InvalidColor', ...
            [ 'tree ' num2str( i ) ' has an invalid color' ] );
    end
    if sc( 2 ) == 1
        tree( i ).color = tree( i ).color * ...
            ones( 1, cylsPerTree( i ) + 1 );
    elseif sc( 2 ) ~= cylsPerTree( i ) + 1
        error( 'bole2ply:InvalidColor', ...
            [ 'tree ' num2str( i ) ' has an invalid color' ] );
    end
    % make sure the tree is a uint8 (so meshlab won't crash on Dave)
    tree( i ).color = uint8( tree( i ).color );
end

numVerts = sum( cylsPerTree + 1 ) * tess;
numFaces = sum( cylsPerTree ) * tess;

% write the header
fprintf( fid, 'ply\n' );
fprintf( fid, 'format ascii 1.0\n' );
fprintf( fid, 'comment author: Paul Romanczyk\n' );
fprintf( fid, 'element vertex %d\n', numVerts );
fprintf( fid, 'property float x\n' );
fprintf( fid, 'property float y\n' );
fprintf( fid, 'property float z\n' );
fprintf( fid, 'property uchar red\n' );
fprintf( fid, 'property uchar green\n' );
fprintf( fid, 'property uchar blue\n' );
fprintf( fid, 'element face %d\n', numFaces *2 );
fprintf( fid, 'property list uchar int vertex_index\n' );
fprintf( fid, 'end_header\n' );

% construct the base cylinder

theta = linspace( 0, 2 * pi, tess + 1 );
theta = theta( 1:end - 1 ); % use only the minimum number
xyz = [ cos( theta ); sin( theta ); ones( size( theta ) ) ];

% write the vertices
for i = 1:n
    % write the bottom circle
    xyzp = xyz .* tree( i ).r( 1 );
    
    % compute the rotations
    delta = tree( i ).loc( :, 2 ) - tree( i ).loc( :, 1 );
    H = sqrt( sum( delta .^2 ) );
        if H == 0
            error( 'bole2ply:DegenerateCase', ...
                [ 'Dave: You have a degererate case in tree ' ...
                num2str( i ) ' between circle 1 and 2.' ] );
        end
    
    rz = atan2( delta( 2 ), delta( 1 ) ) + pi / 2;
    rx = acos( delta( 3 ) / H );

    % compute the rotation matrix
    RZ = [ cos( rz ) -sin( rz ) 0; sin( rz ) cos( rz ) 0; 0 0 1 ];
    RX = [ 1 0 0; 0 cos( rx ) -sin( rx ); 0 sin( rx ) cos( rx ) ];

    % rotate the points
    xyzp = RZ * RX * xyzp + tree( i ).loc( :, 1 ) * ones( 1, tess );
    
    idx = find( xyzp( 1, : ) == max( xyzp( 1, : ) ), 1, 'first' );
    idx = undistortPts( idx, tess );
    xyzp = xyzp( :, idx );
    
    % write the points
    for k = 1:tess
        fprintf( fid, '%f %f %f %d %d %d\n', ...
            xyzp( :, k ), ...
            tree( i ).color( :, end ) );
      %  fprintf( '%f \n', ...
      %      xyzp( 3, k ) );
    end
    
    % write the middle circle
    for j = 2:cylsPerTree( i )   

        xyzp = xyz .* tree( i ).r( j );

        % compute the rotations
        delta = tree( i ).loc( :, j + 1 ) - tree( i ).loc( :, j - 1 );
        H = sqrt( sum( delta .^2 ) );
        if H == 0
            error( 'bole2ply:DegenerateCase', ...
                [ 'Dave: You have a degererate case in tree ' ...
                num2str( i ) ' between circle ' num2str( j - 1 ) ...
                ' and ' num2str( j + 1 ) '.'] );
        end
        rz = atan2( delta( 2 ), delta( 1 ) ) + pi / 2;
        rx = acos( delta( 3 ) / H );

        % compute the rotation matrix
        RZ = [ cos( rz ) -sin( rz ) 0; sin( rz ) cos( rz ) 0; 0 0 1 ];
        RX = [ 1 0 0; 0 cos( rx ) -sin( rx ); 0 sin( rx ) cos( rx ) ];

        % rotate the points
        xyzp = RZ * RX * xyzp + tree( i ).loc( :, j ) * ones( 1, tess );

        idx = find( xyzp( 1, : ) == max( xyzp( 1, : ) ), 1, 'first' );
        idx = undistortPts( idx, tess );
        xyzp = xyzp( :, idx );

        % write the points
        for k = 1:tess
            fprintf( fid, '%f %f %f %d %d %d\n', ...
                xyzp( :, k ), ...
                tree( i ).color( :, end ) );
        end
   
    end
    
    % write the top circle
    xyzp = xyz .* tree( i ).r( end );
    % need to re-order the top points since they are different than the
    % rest
    xyzp( 1, : ) = -1 .* xyzp( 1, : );
    
    % compute the rotations
    delta = tree( i ).loc( :, end - 1 ) - tree( i ).loc( :, end );
    H = sqrt( sum( delta .^2 ) );
        if H == 0
            error( 'bole2ply:DegenerateCase', ...
                [ 'Dave: You have a degererate case in tree ' ...
                num2str( i ) ' between circle ' ...
                num2str( cylsPerTree( i ) ) ' and ' ...
                num2str( cylsPerTree( i ) + 1 ) '.' ] );
        end
    
    rz = atan2( delta( 2 ), delta( 1 ) ) + pi / 2;
    rx = acos( delta( 3 ) / H );

    % compute the rotation matrix
    RZ = [ cos( rz ) -sin( rz ) 0; sin( rz ) cos( rz ) 0; 0 0 1 ];
    RX = [ 1 0 0; 0 cos( rx ) -sin( rx ); 0 sin( rx ) cos( rx ) ];

    % rotate the points
    xyzp = RZ * RX * xyzp + tree( i ).loc( :, end ) * ones( 1, tess );
    
    idx = find( xyzp( 1, : ) == max( xyzp( 1, : ) ), 1, 'first' );
    idx = undistortPts( idx, tess );
    xyzp = xyzp( :, idx );
    
    % write the points
    for k = 1:tess
        fprintf( fid, '%f %f %f %d %d %d\n', ...
            xyzp( :, k ), ...
            tree( i ).color( :, end ) );
    end
end


% these values start at 0!
lr = 0:tess - 1;
ll = [ lr( 2:end ) lr( 1 ) ];
ul = tess + [ lr( 2:end ) lr( 1 ) ];
ur = lr + tess;

% write the faces
% compute the base vertex indices

base = 0;
for i = 1:n
    for j = 1:cylsPerTree( i )
        for k = 1:tess
            fprintf( fid, '3 %d %d %d\n', ...
                ( base + ( j - 1 ) * tess ) + ...
                [ lr( k ) ll( k ) ur( k ) ] );
            fprintf( fid, '3 %d %d %d\n', ...
                ( base + ( j - 1 ) * tess ) + ...
                [ ll( k ) ul( k ) ur( k ) ] );
        end
    end
    base = base + ( cylsPerTree( i ) + 1 ) * tess;
end


fclose( fid );

end

% This is used to "un rotate" the circls caused by the rz rotation. The
% points are the same, just in a different (better) order for display.
% Leaving this out causes "twists" in the tree model.
function out = undistortPts( idx, tess )
    l = 1:tess;
    out = [ l( idx:end ), l( 1:idx-1 ) ];
end

