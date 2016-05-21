function isComplete = write2ply( outputFileName, pointCloud1, color1, pointCloud2, color2 )
% Name: write2ply
% Description: write a .ply file (with optional colors) from point cloud
% data.
% Inputs:
%   outputFileName - a string containing the name of the output file. Can
%       optionally include a path. Must contain a .ply on the end.
%   pointCloud1 - an n x 3 matrix of points.  The first column contains the
%       x locations, the second the y locations, and the third the z
%       locations.
%   color1 - (optional) an n x 3 or 1 x 3 matrix of colors.  The values are
%       scaled in the range [0,255].  If the size is n x 3, each point in 
%       pointCloud1 has its own color.  If the size is 1 x 3, all points in
%       pointCloud1 have the same color.
%   pointCloud2 - (optional) an m x 3 matrix of points.  The first column 
%       contains the x locations, the second the y locations, and the third 
%       the z locations. If this is set, color2 must also be set.
%   color2 - (optional) an m x 3 or 1 x 3 matrix of colors.  The values are
%       scaled in the range [0,255].  If the size is m x 3, each point in 
%       pointCloud1 has its own color.  If the size is 1 x 3, all points in
%       pointCloud1 have the same color.
% Outputs:
%   isComplete - true if complete, false otherwise.
%       The file specified by outputFileName is created if there is
%       properly formatted input data
% History:
%   ????-??-?? Dave Nilosek (drn2369@rit.edu) Wrote initial version for
%       writing a single point cloud with no color information.
%   2012-03-30 Paul Romanczyk (par4249@rit.edu) Added ability to write out
%       with color, and to have 2 point clouds in the same file. 

% setup for the writing of the file. Also check that the data is properly
% formatted
isComplete = false;

if nargin == 0
    % no options set
    fprintf( 'write2ply\n' );
    fprintf( [ 'Usage: write2ply( outputFileName, pointCloud1, ' ...
        'color1, pointCloud2, color2 ).\n' ] );
    writeFile = false;
elseif nargin == 1
    % no point cloud set
    fprintf( 'write2ply: No point cloud set.\n' );
    fprintf( [ 'Usage: write2ply( outputFileName, pointCloud1, ' ...
        'color1, pointCloud2, color2 ).\n' ] );
    writeFile = false;
elseif nargin == 2
    % only 1 point cloud
    % no colors set
    writeFile = true;
    hasColors = false;
    s = size( pointCloud1 );
    if size( s ) ~= 2
        writeFile = false;
        fprintf( 'write2ply: invalid pointCloud1 size.\n' );
    elseif s( 2 ) ~= 3
        writeFile = false;
        fprintf( 'write2ply: invalid pointCloud1 size.\n' );
    else
        numPoints = s( 1 );
        pointCloud = pointCloud1;
    end
elseif nargin == 3
    % 1 point cloud with color
    writeFile = true;
    hasColors = true;
    numPoints = length( pointCloud1 );
    pointCloud = pointCloud1;
    s = size( pointCloud1 );
    if size( s ) ~= 2
        writeFile = false;
        fprintf( 'write2ply: invalid pointCloud1 size.\n' );
    elseif s( 2 ) ~= 3
        writeFile = false;
        fprintf( 'write2ply: invalid pointCloud1 size.\n' );
    else
        numPoints = s( 1 );
        pointCloud = pointCloud1;
    end
    s = size( color1 );
    if ( ( size( s ) ~= 2 ) & ( writeFile ) )
        writeFile = false;
        fprintf( 'write2ply: invalid color size.\n' );
    elseif ( ( s( 2 ) ~= 3 ) & ( writeFile ) )
        writeFile = false;
        fprintf( 'write2ply: invalid color size.\n' );
    else
        if ( ( s( 1 ) == 1 ) & ( writeFile ) )
            colors = ones( numPoints, 1 ) * color1;
        elseif ( ( s( 1 ) == numPoints ) & ( writeFile ) );
            colors = color1;
        else
            writeFile = false;
            fprintf( 'write2ply: invalid color size.\n' );
        end
    end
elseif nargin == 4
    writeFile = false;
    fprintf( 'write2ply: Invalid Useage\n' );
    fprintf( [ 'Usage: write2ply( outputFileName, pointCloud1, ' ...
        'color1, pointCloud2, color2 ).\n' ] );
elseif nargin == 5
    % 1 point cloud with color
    writeFile = true;
    hasColors = true;
    s1 = size( pointCloud1 );
    s2 = size( pointCloud2 );
    if ( ( size( s1 ) ~= 2 ) | ( size( s2 ) ~= 2 ) )
        writeFile = false;
        fprintf( 'write2ply: invalid pointCloud size.\n' );
    elseif ( ( s1( 2 ) ~= 3 ) | ( s2( 2 ) ~= 3 ) )
        writeFile = false;
        fprintf( 'write2ply: invalid pointCloud size.\n' );
    else
        numPoints = s1( 1 ) + s2( 1 );
        pointCloud = [ pointCloud1; pointCloud2 ];
    end
    s1 = size( color1 );
    s2 = size( color2 );
    %% add in color1
    if ( ( size( s1 ) ~= 2 ) & ( writeFile ) )
        writeFile = false;
        fprintf( 'write2ply: invalid color1 size.\n' );
    else
        if ( ( s1( 1 ) == 1 ) & ( writeFile ) )
            colors = ones( size( pointCloud1, 1 ), 1 ) * color1;
        elseif ( ( s1( 1 ) == size( pointCloud1, 1 ) ) & ( writeFile ) );
            colors = color1;
        else
            writeFile = false;
            fprintf( 'write2ply: invalid color1 size.\n' );
        end
    end
    
    if ( ( size( s2 ) ~= 2 ) & ( writeFile ) )
        writeFile = false;
        fprintf( 'write2ply: invalid color2 size.\n' );
    else
        if ( ( s2( 1 ) == 1 ) & ( writeFile ) )
            colors = [ colors; ones( length( pointCloud2 ), 1 ) * color2 ];
        elseif ( ( s2( 1 ) == size( pointCloud2, 1 ) ) & ( writeFile ) );
            colors = [ colors; color2 ];
        else
            writeFile = false;
            fprintf( 'write2ply: invalid color2 size.\n' );
        end
    end
else
    % This should not happen
    disp( 'write2ply: Some error.' );
    writeFile = false;
end

% write the file
if ( writeFile )
    ID = fopen( outputFileName, 'w+' );
    
    fprintf( ID, 'ply\n');
    fprintf( ID, 'format ascii 1.0\n' );
    fprintf( ID, 'element vertex %d\n', numPoints );
    fprintf( ID, 'property float x\n' );
    fprintf( ID, 'property float y\n' );
    fprintf( ID, 'property float z\n' );
    if hasColors
        fprintf( ID, 'property uchar red\n' );
        fprintf( ID, 'property uchar green\n' );
        fprintf( ID, 'property uchar blue\n' );
    end
    fprintf( ID, 'end_header\n' );
    if hasColors
        % point cloud with colors
        for i = 1:numPoints
            fprintf( ID, '%f %f %f %d %d %d\n', pointCloud( i, 1:3 ),uint8(colors(i,1:3)));
        end
    else
        % points only
        for i = 1:numPoints
            fprintf( ID, '%f %f %f\n', pointCloud( i, 1:3 ));
        end
    end
    fprintf( ID, '0 0 0 0\n' );
    fclose( ID );
    
    isComplete = true;
end