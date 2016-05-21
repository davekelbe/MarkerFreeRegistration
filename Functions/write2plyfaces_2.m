function isComplete = write2plyfaces_2( outputFileName, pointCloud1, color1, faces1 )
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
%   faces1 - (optional) a p x m matrix of vertex indices. There are p faces
%   and m vertex indices per face
% Outputs:
%   isComplete - true if complete, false otherwise.
%       The file specified by outputFileName is created if there is
%       properly formatted input data
% History:
%   ????-??-?? Dave Nilosek (drn2369@rit.edu) Wrote initial version for
%       writing a single point cloud with no color information.
%   2012-03-30 Paul Romanczyk (par4249@rit.edu) Added ability to write out
%       with color, and to have 2 point clouds in the same file. 
%   2012-03-30 Dave Kelbe (djk2312@rit.edu) Added ability to write out
%       with faces and removed 2 point cloud addition

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
    hasFaces = false;
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
    hasFaces = false;
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
    % 1 point cloud with color and faces
    writeFile = true;
    hasColors = true;
    hasFaces = true;
    numPoints = length( pointCloud1 );
    %pointCloud = pointCloud1;
    s = size( pointCloud1 );
    %faces = faces1;
    f = size(faces1);
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
    if size( f ) ~= 2
        writeFile = false;
        fprintf( 'write2ply: invalid faces1 size.\n' );
    elseif s( 2 ) > s( 1 )
        writeFile = false;
        fprintf( 'write2ply: invalid faces1 size.\n' );
    else
        numFaces = f( 1 );
        numIndices = f( 2);
        faces = faces1;
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
end
% write the file
colors = colors/max(max(colors));
colors = colors*255;
colors = uint8(colors);
colors = num2str(colors);
%pointCloud = num2str(pointCloud);

if ( writeFile )
    ID = fopen( outputFileName, 'w+' );
    
    fprintf( ID, 'ply\n');
    fprintf( ID, 'format ascii 1.0\n' );
    fprintf( ID, 'element vertex %d\n', numPoints );
    fprintf( ID, 'property float x\n' );
    fprintf( ID, 'property float y\n' );
    fprintf( ID, 'property float z\n' );
    if hasColors
        fprintf( ID, 'property uint red\n' );
        fprintf( ID, 'property uint green\n' );
        fprintf( ID, 'property uint blue\n' );
    end
    if hasFaces
        fprintf( ID, 'element face %d\n', numFaces );
        fprintf( ID, 'property list uchar uint vertex_indices \n' );
    end
    fprintf( ID, 'end_header\n' );
    if (hasColors && hasFaces)
        % point cloud with colors
        for i = 1:numPoints
            fprintf( ID, ' %f', pointCloud( i, 1 ) );
            fprintf( ID, ' %f', pointCloud( i, 2 ) );
            fprintf( ID, ' %f', pointCloud( i, 3 ) );
            fprintf( ID, ' %u', color1( i, 1 ) );
            fprintf( ID, ' %u', color1( i, 2 ) );
            fprintf( ID, ' %u\n', color1( i, 3 ) );
        end
        for i = 1:numFaces
            fprintf( ID, ' %d', numIndices );
            for j = 1:numIndices 
            fprintf( ID, ' %u', faces( i, j )-1 );
            if j == numIndices
                fprintf( ID, '\n');
            end
            end
        end
    elseif (hasColors && ~hasFaces)
        % point cloud with colors
        for i = 1:numPoints
            fprintf( ID, ' %f', pointCloud( i, 1 ) );
            fprintf( ID, ' %f', pointCloud( i, 2 ) );
            fprintf( ID, ' %f', pointCloud( i, 3 ) );
            fprintf( ID, ' %u', color1( i, 1 ) );
            fprintf( ID, ' %u', color1( i, 2 ) );
            fprintf( ID, ' %u\n', color1( i, 3 ) );
        end
    elseif (~hasColors && hasFaces)
        % points only
        for i = 1:numPoints
            fprintf( ID, ' %f', pointCloud( i, 1 ) );
            fprintf( ID, ' %f', pointCloud( i, 2 ) );
            fprintf( ID, ' %f\n', pointCloud( i, 3 ) );
        end
        for i = 1:numFaces
            fprintf( ID, ' %d', numIndices );
            for j = 1:numIndices 
            fprintf( ID, ' %d', faces( i, j ) );
            if j == numIndices
                fprintf( ID, '\n');
            end
            end
        end
    elseif (~hasColors && ~hasFaces)
        % points only
        for i = 1:numPoints
            fprintf( ID, ' %f', pointCloud( i, 1 ) );
            fprintf( ID, ' %f', pointCloud( i, 2 ) );
            fprintf( ID, ' %f\n', pointCloud( i, 3 ) );
        end
    end
    fprintf( ID, '0 0 0 0\n' ); % not sure ? 
    fclose( ID );
    
    isComplete = true;
end