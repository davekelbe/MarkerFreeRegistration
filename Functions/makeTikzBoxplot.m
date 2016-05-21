function makeTikzBoxplot( data, varargin )

fid = 1;
labels = {};
labelPrefix = '';
title = '';
name = '';
boxsep = 2;
boxwidth = 1.75;
fill = '';
figwidth = '0.6\textwidth';
figheight = '';
legendEntries = '';
legendPrefix = '';
makelegend = true;
legendAnchor = 'south';
legendPos = [ 0.5, 1.02 ];
outlierMark = 'o';
meanMark = 'asterisk';
axisLoc = ''; %at=(discrete.right of north east),%
axisAnchor = ''; %anchor=left of north west,%
ylabel = '';
xlabel = '';
        

i = 1;
while i < numel( varargin )
    switch lower( varargin{ i } )
        case 'fid'
            fid = varargin{ i + 1 };
            i = i + 2;
        case 'labels'
            labels = varargin{ i + 1 };
            i = i + 2;
        case 'labelprefix'
            labelPrefix = varargin{ i + 1 };
            i = i + 2;
        case 'title'
            title = varargin{ i + 1 };
            i = i + 2;
        case 'name'
            name = varargin{ i + 1 };
            i = i + 2;
        case 'boxsep'
            boxsep = varargin{ i + 1 };
            i = i + 2;
        case 'boxwidth'
            boxwidth = varargin{ i + 1 };
            i = i + 2;
        case 'fill'
            fill = varargin{ i + 1 };
            i = i + 2;
        case 'figwidth'
            figwidth = varargin{ i + 1 };
            i = i + 2;
        case 'figheight'
            figheight = varargin{ i + 1 };
            i = i + 2;
        case 'legendlabels'
            legendEntries = varargin{ i + 1 };
            i = i + 2;
        case 'legend'
            switch varargin{ i + 1 }
                case 'off'
                    makelegend = false;
                otherwise
                    makelegend = true;
            end
            i = i + 2;
        case 'nolegend'
            makelegend = false;
            i = i + 1;
        case 'legendprefix'
            legendPrefix = varargin{ i + 1 };
            i = i + 2;
        case 'legendanchor'
            legendAnchor = varargin{ i + 1 };
            i = i + 2;
        case 'legendpos'
            legendPos = varargin{ i + 1 };
            i = i + 2;
        case 'xlegendpos'
            legendPos( 1 ) = varargin{ i + 1 };
            i = i + 2;
        case 'ylegendpos'
            legendPos( 2 ) = varargin{ i + 1 };
            i = i + 2;
        case 'meanmark'
            meanMark = varargin{ i + 1 };
            i = i + 2;
        case 'outliermark'
            outlierMark = varargin{ i + 1 };
            i = i + 2;
        case 'axisloc'
            axisLoc = varargin{ i + 1 };
            i = i + 2;
        case 'axisanchor'
            axisAnchor = varargin{ i + 1 };
            i = i + 2;
        case 'ylabel'
            ylabel = varargin{ i + 1 };
            i = i + 2;
        case 'xlabel'
            xlabel = varargin{ i + 1 };
            i = i + 2;
        otherwise
            
            error( [ 'Unknown input ''', varargin{ i }, ...
                ''' to makeTikzBoxPlot' ] );
    end
end

q25 = quantile( data, 0.25 );
q50 = quantile( data, 0.50 );
q75 = quantile( data, 0.75 );
qr = ( q75 - q25 ) .* 1.5;
m = mean( data );
lt = q25 - qr;
ut = q75 + qr;
n = size( data, 2 );

fprintf( fid, '\\begin{tikzpicture}\n' );
fprintf( fid, '\t\\begin{axis}\n' );
fprintf( fid, '\t\t[%%\n' );
if ~isempty( ylabel )
    fprintf( fid, '\t\tylabel=%s,%%\n', ylabel );
    fprintf( fid, '\t\ty label style={%%\n' );
	fprintf( fid, '\t\t\tat={(axis description cs:0.025,.5)},%%\n' );
	fprintf( fid, '\t\t},%%\n' );
end
if ~isempty( xlabel )
    fprintf( fid, '\t\txlabel=%s,%%\n', xlabel );
end
if ~isempty( axisLoc )
    fprintf( fid, '\t\tat=(%s),%%\n', axisLoc );
end
if ~isempty( axisAnchor )
    fprintf( fid, '\t\tanchor=%s,%%\n', axisAnchor );
end
if ~isempty( figwidth )
    fprintf( fid, '\t\twidth=%s,%%\n', figwidth );
end
if ~isempty( figheight )
    fprintf( fid, '\t\theight=%s,%%\n', figheight );
end
fprintf( fid, '\t\txmin=0,%%\n' );
fprintf( fid, '\t\txmax=%f,%%\n', n * boxsep );
fprintf( fid, '\t\tymin=0,%%\n' );
fprintf( fid, '\t\tymax=1,%%\n' );
if ~isempty( title );
    fprintf( fid, '\t\ttitle={%s},%%\n', title );
end
if ~isempty( name );
    fprintf( fid, '\t\tname={%s},%%\n', name );
end
fprintf( fid, '\t\tytick pos=left,%%\n' );
fprintf( fid, '\t\ttick label style={%%\n' );
fprintf( fid, '\t\t\tmajor tick length=0pt%%\n' );
fprintf( fid, '\t\t},%%\n' );

% fprintf( fid, '\t\tboxplot/every boxplot/.style={%%\n' );
% fprintf( fid, '\t\t\tdraw=black,%%\n' );
% fprintf( fid, '\t\t},%%\n' );

if ~isempty( fill )
    fprintf( fid, '\t\tboxplot/every box/.style={%%\n' );
    fprintf( fid, '\t\t\tfill=%s,%%\n', fill );
    fprintf( fid, '\t\t},%%\n' );
end
fprintf( fid, '\t\tboxplot/every average/.style={%%\n' );
fprintf( fid, '\t\t\t/tikz/mark=%s,%%\n', meanMark );
% fprintf( fid, '\t\t\t/tikz/solid,%%\n' );
fprintf( fid, '\t\t},%%\n' );
fprintf( fid, '\t\tboxplot/draw direction=y,%%\n' );
fprintf( fid, '\t\txtick={1' );
for i = 2:n
    fprintf( fid, ',%d', boxsep * ( i - 1/2 ) );
end
fprintf( fid, '},%%\n' );
if numel( labels ) > 0
    fprintf( fid, '\t\txticklabels={%s%s', labelPrefix, labels{ 1 } );
    for i = 2:n
        fprintf( fid, ',%s%s', labelPrefix, labels{ i } );
    end
    fprintf( fid, '},%%\n' );
end
if makelegend
    fprintf( fid, '\t\tlegend style={%%\n' );
    % fprintf( fid, '\t\t\tdraw=none,%%\n' );
    fprintf( fid, [ '\t\t\t/tikz/every even column/.append style={', ...
        'column sep=2ex},%%\n' ] );
    % fprintf( fid, '\t\t\ttext width=2em,%%\n' );
    % fprintf( fid, '\t\t\ttext depth=.5ex,%%\n' );
    fprintf( fid, '\t\t\tat={(%f,%f)},%%\n', legendPos );
    fprintf( fid, '\t\t\tanchor=%s,%%\n', legendAnchor );
    fprintf( fid, '\t\t},%%\n' );
    fprintf( fid, '\t\tlegend entries={%smean,%soutlier', ... 
        legendPrefix, legendPrefix );
    if ~isempty( legendEntries )
        for i = 1:numel( legendEntries )
            fprintf( fid, ',%s%s', legendPrefix, legendEntries{ i } );
        end
    end
    fprintf( fid, '},%%\n' );
    fprintf( fid, '\t\tlegend columns=-1,%%\n' );
end

fprintf( fid, '\t\t]\n' );
if makelegend
    fprintf( fid, '\t\t\\addlegendimage{%%\n' );
    fprintf( fid, '\t\t\tonly marks,%%\n' );
    fprintf( fid, '\t\t\t/tikz/mark=%s,%%\n', meanMark );
    fprintf( fid, '\t\t}\n' );
    fprintf( fid, '\t\t\\addlegendimage{%%\n' );
    fprintf( fid, '\t\t\tonly marks,%%\n' );
    fprintf( fid, '\t\t\tmark=%s\n', outlierMark );
    fprintf( fid, '\t\t}\n' );
    if ~isempty( legendEntries )
        fprintf( fid, '\t\t\\addlegendimage{%%\n' );
        fprintf( fid, '\t\t\tblack,%%\n' );
        if ~isempty( fill ) 
            fprintf( fid, '\t\t\tfill=%s,%%\n', fill );
        end
        fprintf( fid, '\t\t\tarea legend\n' );
        fprintf( fid, '\t\t}\n' );
    end
end

for i = 1:n
    % compute the whiskers 
    lw = min( data( data( :, i ) > lt( i ), i ) );
    uw = max( data( data( :, i ) < ut( i ), i ) );
    
    
    idxl = find( data( :, i ) < lt( i ) );
    idxu = find( data( :, i ) > ut( i ) );
    ol = data( [ idxl, idxu ], i );

    if numel( ol ) == 0
        ols = '';
    else
        ols = '(0.0,%f)';
        for j = 2:numel( ol )
            ols = strcat( ols, ' (0.0,%f)' );
        end
    end
    
    
    
    fprintf( fid, '\t\t\\addplot[%%\n' );
    fprintf( fid, '\t\t\tmark=%s,%%\n', outlierMark );
    fprintf( fid, '\t\t\tboxplot prepared={%%\n' );
    fprintf( fid, '\t\t\t\tdraw position=%d,%%\n', boxsep * ( i - 0.5 ) );
    fprintf( fid, '\t\t\t\tmedian=%f,%%\n', q50( i ) );
    fprintf( fid, '\t\t\t\taverage=%f,%%\n', m( i ) );
    fprintf( fid, '\t\t\t\tupper quartile=%f,%%\n', q75( i ) );
    fprintf( fid, '\t\t\t\tlower quartile=%f,%%\n', q25( i ) );
    fprintf( fid, '\t\t\t\tupper whisker=%f,%%\n', uw );
    fprintf( fid, '\t\t\t\tlower whisker=%f,%%\n', lw );
    fprintf( fid, '\t\t\t\tbox extend=%f,%%\n', boxwidth );
    fprintf( fid, '\t\t\t},%%\n' );
    fprintf( fid, [ '\t\t] coordinates {', ols '};\n' ], ol );
end
fprintf( fid, '\t\\end{axis}\n' );
fprintf( fid, '\\end{tikzpicture}\n' );
end