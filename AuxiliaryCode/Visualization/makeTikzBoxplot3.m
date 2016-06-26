function makeTikzBoxplot3( data_1, data_2, data_3, varargin )

% 2014-06-06 Dave Kelbe received from Paul Romanczyk
% 2014-06-06 Dave Kelbe converted to use cells, where data_1 and data_2 need
% not have identical number of elements


fid = 1;
labels = {};
labelPrefix = '';
title = '';
name = '';
boxsep = 2;
boxinnersep = 0;
boxwidth = 0.6;
fill1 = 'gray!25';
fill2 = '';
fill3 = 'gray!60';
figwidth = '0.9\textwidth';
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
ylim = [-1 1];


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
        case 'boxinnersep'
            boxinnersep = varargin{ i + 1 };
            i = i + 2;
        case 'boxwidth'
            boxwidth = varargin{ i + 1 };
            i = i + 2;
        case 'fill'
            tmp = varargin{ i + 1 };
            fill1 = tmp{ 1 };
            fill2 = tmp{ 2 };
            fill3 = tmp{ 3 };
            i = i + 2;
        case 'fill1'
            fill1 = varargin{ i + 1 };
            i = i + 2;
        case 'fill2'
            fill2 = varargin{ i + 1 };
            i = i + 2;
        case 'fill3'
            fill3 = varargin{ i + 1 };
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
        case 'ylim'
            ylim = varargin{ i + 1 };
            i = i + 2;
        otherwise
            error( 'Unknown input to makeTikzBoxPlot3' ); 
    end
end

if iscell(data_1);
    n_p = numel(data_1);
    
    q25_1 = zeros(n_p,1);
    q50_1 = zeros(n_p,1);
    q75_1 = zeros(n_p,1); 
    m_1 = zeros(n_p,1);
    n_1 = zeros(n_p,1);
    
    q25_2 = zeros(n_p,1);
    q50_2 = zeros(n_p,1);
    q75_2 = zeros(n_p,1);
    m_2 = zeros(n_p,1);
    
    q25_3 = zeros(n_p,1);
    q50_3 = zeros(n_p,1);
    q75_3 = zeros(n_p,1);
    m_3 = zeros(n_p,1);
    
    for p = 1:n_p;
    q25_1(p) = quantile( data_1{p}, 0.25 );
    q50_1(p) = quantile( data_1{p}, 0.50 );
    q75_1(p) = quantile( data_1{p}, 0.75 );
    m_1(p) = mean( data_1{p} );
    n_1(p) = size( data_1{p}, 1 );
    
    q25_2(p) = quantile( data_2{p}, 0.25 );
    q50_2(p) = quantile( data_2{p}, 0.50 );
    q75_2(p) = quantile( data_2{p}, 0.75 );
    m_2(p) = mean( data_2{p} );

    % n_2 = size( data_2, 2 );
    
    q25_3(p) = quantile( data_3{p}, 0.25 );
    q50_3(p) = quantile( data_3{p}, 0.50 );
    q75_3(p) = quantile( data_3{p}, 0.75 );
    m_3(p) = mean( data_3{p} );

    % n_3 = size( data_3, 2 );
    end
    qr_1 = ( q75_1 - q25_1 ) .* 1.5;
    qr_2 = ( q75_2 - q25_2 ) .* 1.5;
    qr_3 = ( q75_3 - q25_3 ) .* 1.5;
    lt_1 = q25_1 - qr_1; 
    ut_1 = q75_1 + qr_1;
    lt_2 = q25_2 - qr_2;
    ut_2 = q75_2 + qr_2;
    lt_3 = q25_3 - qr_3;
    ut_3 = q75_3 + qr_3;
else 
    q25_1 = quantile( data_1, 0.25 );
    q50_1 = quantile( data_1, 0.50 );
    q75_1 = quantile( data_1, 0.75 );
    qr_1 = ( q75_1 - q25_1 ) .* 1.5;
    m_1 = mean( data_1 );
    lt_1 = q25_1 - qr_1;
    ut_1 = q75_1 + qr_1;
    n_1 = size( data_1, 2 );
     
    q25_2 = quantile( data_2, 0.25 );
    q50_2 = quantile( data_2, 0.50 );
    q75_2 = quantile( data_2, 0.75 );
    qr_2 = ( q75_2 - q25_2 ) .* 1.5;
    m_2 = mean( data_2 );
    lt_2 = q25_2 - qr_2;
    ut_2 = q75_2 + qr_2;
    % n_2 = size( data_2, 2 );
    
    q25_3 = quantile( data_3, 0.25 );
    q50_3 = quantile( data_3, 0.50 );
    q75_3 = quantile( data_3, 0.75 );
    qr_3 = ( q75_3 - q25_3 ) .* 1.5;
    m_3 = mean( data_3 );
    lt_3 = q25_3 - qr_3;
    ut_3 = q75_3 + qr_3;
    % n_3 = size( data_3, 2 );
end



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
fprintf( fid, '\t\txmax=%f,%%\n', n_p * boxsep );
fprintf( fid, '\t\tymin=%3.2f,%%\n',ylim(1) );
fprintf( fid, '\t\tymax=%3.2f,%%\n', ylim(2) );
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
fprintf( fid, '\t\tboxplot/every average/.style={%%\n' );
fprintf( fid, '\t\t\t/tikz/mark=%s,%%\n', meanMark );
fprintf( fid, '\t\t},%%\n' );
fprintf( fid, '\t\tboxplot/draw direction=y,%%\n' );
fprintf( fid, '\t\txtick={1' );
for i = 2:numel(n_1)
    fprintf( fid, ',%d', boxsep * ( i - 1/2 ) );
end
fprintf( fid, '},%%\n' );
if numel( labels ) > 0
    fprintf( fid, '\t\txticklabels={%s%s', labelPrefix, labels{ 1 } );
    for i = 2:numel(n_1);
        fprintf( fid, ',%s', labels{ i } );
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
        if ~isempty( fill1 )
            fprintf( fid, '\t\t\tfill=%s,%%\n', fill1 );
        end
        fprintf( fid, '\t\t\tarea legend\n' );
        fprintf( fid, '\t\t}\n' );
        fprintf( fid, '\t\t\\addlegendimage{%%\n' );
        fprintf( fid, '\t\t\tblack,%%\n' );
        if ~isempty( fill2 )
            fprintf( fid, '\t\t\tfill=%s,%%\n', fill2 );
        end
        fprintf( fid, '\t\t\tarea legend\n' );
        fprintf( fid, '\t\t}\n' );
        fprintf( fid, '\t\t\\addlegendimage{%%\n' );
        fprintf( fid, '\t\t\tblack,%%\n' );
        if ~isempty( fill3 )
            fprintf( fid, '\t\t\tfill=%s,%%\n', fill3 );
        end
        fprintf( fid, '\t\t\tarea legend\n' );
        fprintf( fid, '\t\t}\n' );
    end
end 

for i = 1:n_p
    fprintf('\ni = %g\n',i)
    if i==3;
        foo = 1;
    end
    
    if isempty(data_1{i});
        continue
    end
    
    % compute the whiskers
    lw_1 = min( data_1{i}( data_1{i} > lt_1( i )) );
    uw_1 = max( data_1{i}( data_1{i} < ut_1( i ) ) );
    
    
    idxl = find( data_1{i} < lt_1( i ) );
    idxu = find( data_1{i} > ut_1( i ) );
    if isempty(idxl) && isempty(idxu);
        ol_1 = [];
    else
        ol_1 = data_1{i}( [ idxl; idxu ] );
    end
    
    if numel( ol_1 ) == 0
        ols_1 = '';
    else
        ols_1 = '(0.0,%f)';
        for j = 2:numel( ol_1 )
            ols_1 = strcat( ols_1, ' (0.0,%f)' );
        end
    end
    
    lw_2 = min( data_2{i}( data_2{i} > lt_2( i ) ) );
    uw_2 = max( data_2{i}( data_2{i} < ut_2( i ) ) );
    
    
    idxl = find( data_2{i} < lt_2( i ) );
    idxu = find( data_2{i} > ut_2( i ) );
    if size(idxl) ~= size(idxu);
        foo = 1;
    end
    ol_2 = data_2{i}( [ idxl; idxu ]);
    
    if numel( ol_2 ) == 0
        ols_2 = '';
    else
        ols_2 = '(0.0,%f)';
        for j = 2:numel( ol_2 )
            ols_2 = strcat( ols_2, ' (0.0,%f)' );
        end
    end
    
    lw_3 = min( data_3{i}( data_3{i} > lt_3( i ) ) );
    uw_3 = max( data_3{i}( data_3{i} < ut_3( i ) ) );
    
    
    idxl = find( data_3{i} < lt_3( i ) );
    idxu = find( data_3{i} > ut_3( i ) );
    ol_3 = data_3{i}( [ idxl; idxu ]);
    
    if numel( ol_3 ) == 0
        ols_3 = '';
    else
        ols_3 = '(0.0,%f)';
        for j = 2:numel( ol_3 )
            ols_3 = strcat( ols_3, ' (0.0,%f)' );
        end
    end
    
    fprintf( fid, '\t\t\\addplot[%%\n' );
    fprintf( fid, '\t\t\tcolor=black,%%\n' );
    if ~isempty( fill1 )
        fprintf( fid, '\t\t\tfill=%s,%%\n', fill1 );
    end
    fprintf( fid, '\t\t\tmark=%s,%%\n', outlierMark );
    fprintf( fid, '\t\t\tboxplot prepared={%%\n' );
    fprintf( fid, '\t\t\t\tdraw position=%d,%%\n', ...
        boxsep * ( i - 0.5 ) - boxinnersep - boxwidth );
    fprintf( fid, '\t\t\t\tmedian=%f,%%\n', q50_1( i ) );
    fprintf( fid, '\t\t\t\taverage=%f,%%\n', m_1( i ) );
    fprintf( fid, '\t\t\t\tupper quartile=%f,%%\n', q75_1( i ) );
    fprintf( fid, '\t\t\t\tlower quartile=%f,%%\n', q25_1( i ) );
    fprintf( fid, '\t\t\t\tupper whisker=%f,%%\n', uw_1 );
    fprintf( fid, '\t\t\t\tlower whisker=%f,%%\n', lw_1 );
    fprintf( fid, '\t\t\t\tbox extend=%f,%%\n', boxwidth );
    fprintf( fid, '\t\t\t},%%\n' );
    fprintf( fid, [ '\t\t] coordinates {', ols_1 '};\n' ], ol_1 );
    
    fprintf( fid, '\t\t\\addplot[%%\n' );
    fprintf( fid, '\t\t\tcolor=black,%%\n' );
    if ~isempty( fill2 )
        fprintf( fid, '\t\t\tfill=%s,%%\n', fill2 );
    end
    fprintf( fid, '\t\t\tmark=%s,%%\n', outlierMark );
    fprintf( fid, '\t\t\tboxplot prepared={%%\n' );
    fprintf( fid, '\t\t\t\tdraw position=%d,%%\n', boxsep * ( i - 0.5 ) );
    fprintf( fid, '\t\t\t\tmedian=%f,%%\n', q50_2( i ) );
    fprintf( fid, '\t\t\t\taverage=%f,%%\n', m_2( i ) );
    fprintf( fid, '\t\t\t\tupper quartile=%f,%%\n', q75_2( i ) );
    fprintf( fid, '\t\t\t\tlower quartile=%f,%%\n', q25_2( i ) );
    fprintf( fid, '\t\t\t\tupper whisker=%f,%%\n', uw_2 );
    fprintf( fid, '\t\t\t\tlower whisker=%f,%%\n', lw_2 );
    fprintf( fid, '\t\t\t\tbox extend=%f,%%\n', boxwidth );
    fprintf( fid, '\t\t\t},%%\n' );
    fprintf( fid, [ '\t\t] coordinates {', ols_2 '};\n' ], ol_2 );
    
    fprintf( fid, '\t\t\\addplot[%%\n' );
    fprintf( fid, '\t\t\tcolor=black,%%\n' );
    if ~isempty( fill3 )
        fprintf( fid, '\t\t\tfill=%s,%%\n', fill3 );
    end
    fprintf( fid, '\t\t\tmark=%s,%%\n', outlierMark );
    fprintf( fid, '\t\t\tboxplot prepared={%%\n' );
    fprintf( fid, '\t\t\t\tdraw position=%d,%%\n', ...
        boxsep * ( i - 0.5 ) + boxinnersep + boxwidth );
    fprintf( fid, '\t\t\t\tmedian=%f,%%\n', q50_3( i ) );
    fprintf( fid, '\t\t\t\taverage=%f,%%\n', m_3( i ) );
    fprintf( fid, '\t\t\t\tupper quartile=%f,%%\n', q75_3( i ) );
    fprintf( fid, '\t\t\t\tlower quartile=%f,%%\n', q25_3( i ) );
    fprintf( fid, '\t\t\t\tupper whisker=%f,%%\n', uw_3 );
    fprintf( fid, '\t\t\t\tlower whisker=%f,%%\n', lw_3 );
    fprintf( fid, '\t\t\t\tbox extend=%f,%%\n', boxwidth );
    fprintf( fid, '\t\t\t},%%\n' );
    fprintf( fid, [ '\t\t] coordinates {', ols_3 '};\n' ], ol_3 );
    
end
fprintf( fid, '\t\\end{axis}\n' );
fprintf( fid, '\\end{tikzpicture}\n' );