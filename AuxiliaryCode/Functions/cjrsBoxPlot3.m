function cjrsBoxPlot3( data1, data2, data3, labels, isColor )

if ( nargin < 5 )
    isColor = false;
end

w = 0.9/6;
s = 0.6;
w2 = w * 2;


q25 = quantile( data1, 0.25 );
q50 = quantile( data1, 0.50 );
q75 = quantile( data1, 0.75 );
qr = ( q75 - q25 ) .* 1.5;
m = mean( data1 );
lt = q25 - qr;
ut = q75 + qr;

h1 = [];
h2 = [];
h3 = [];
if ( ~isColor )
    c1 = [ 0.5 0.5 0.5 ];
    c2 = [ 0.95 0.95 0.95 ];
    c3 = [ 0.8 0.8 0.8 ];
    C1 = 'k';
    C2 = 'k';
    C3 = 'k';
else
    c1 = [ 0.7 0.2 0.2 ];
    c2 = [ 0.2 0.5 0.2 ];
    c3 = [ 0.2 0.2 0.7 ];
    C1 = c1;
    C2 = c2;
    C3 = c3;
end

hold on;
for i = 1:size( data1, 2 );


    % the box
    h1 = [ h1 fill( [ i - w, i - w, i + w, i + w] - w2, ...
        [ q25( i ), q75( i ), q75( i ), q25( i ) ], ...
        c1 ) ];

    % the median
    plot( [ i - w, i + w ] - w2, ...
        [ q50( i ), q50( i ) ], ...
        'k' );
%             'k', 'LineWidth', 2 );

    % compute the whiskers 
    lw = min( data1( data1( :, i ) > lt( i ), i ) );
    uw = max( data1( data1( :, i ) < ut( i ), i ) );

   
    % the lower whisker
    if ( numel( lw ) > 0 )
        plot( [ i, i ] - w2, ...
            [ q25( i ), lw ], ...
            'k' );

        plot( [ i - w * s, i + w * s ] - w2, ...
            [ lw, lw ], ...
            'k' );
    end
    
    % the upper whisker
    if ( numel( uw ) > 0 )
        plot( [ i, i ] - w2, ...
            [ q75( i ), uw ], ...
            'k' );

        plot( [ i - w * s, i + w * s ] - w2, ...
            [ uw, uw ], ...
            'k' );
    end

    % the mean
    plot( i - w2, m( i ), 'k*' );

    % compute the outliers
    idx = find( data1( :, i ) < lt( i ) );
    for j = 1:numel( idx )
       plot( i - w2, data1( idx( j ), i ), 'o', 'Color', C1 );
    end
    idx = find( data1( :, i ) > ut( i ) );
    for j = 1:numel( idx )
       plot( i - w2, data1( idx( j ), i ), 'o', 'Color', C1 );
    end

%         % plot the data
%         for j = 1:size( data, 1 )
%             plot( i  - w2, data( j, i ), '.r' );
%         end
end

%hatchfill( h1, 'fill' );

q25 = quantile( data2, 0.25 );
q50 = quantile( data2, 0.50 );
q75 = quantile( data2, 0.75 );
qr = ( q75 - q25 ) .* 1.5;
m = mean( data2 );
lt = q25 - qr;
ut = q75 + qr;

for i = 1:size( data2, 2 );


    % the box
    h2 = [ h2 fill( [ i - w, i - w, i + w, i + w ], ...
        [ q25( i ), q75( i ), q75( i ), q25( i ) ], ...
        c2 ) ];

    % the median
    plot( [ i - w, i + w ], ...
        [ q50( i ), q50( i ) ], ...
        'k' );
%             'k', 'LineWidth', 2 );

    % compute the whiskers 
    lw = min( data2( data2( :, i ) > lt( i ), i ) );
    uw = max( data2( data2( :, i ) < ut( i ), i ) );

    % the lower whisker
    if ( numel( lw ) > 0 )
        plot( [ i, i ], ...
            [ q25( i ), lw ], ...
            'k' );

        plot( [ i - w * s, i + w * s ], ...
            [ lw, lw ], ...
            'k' );
    end
    
    % the upper whisker
    if ( numel( uw ) > 0 )
        plot( [ i, i ], ...
            [ q75( i ), uw ], ...
            'k' );

        plot( [ i - w * s, i + w * s ], ...
            [ uw, uw ], ...
            'k' );
    end

    % the mean
    plot( i, m( i ), 'k*' );

    % compute the outliers
    idx = find( data2( :, i ) < lt( i ) );
    for j = 1:numel( idx )
       plot( i, data2( idx( j ), i ), 'o', 'Color', C2 );
    end
    idx = find( data2( :, i ) > ut( i ) );
    for j = 1:numel( idx )
       plot( i, data2( idx( j ), i ), 'o', 'Color', C2 );
    end

%         % plot the data
%         for j = 1:size( data, 1 )
%             plot( i, data( j, i ), '.r' );
%         end
end



q25 = quantile( data3, 0.25 );
q50 = quantile( data3, 0.50 );
q75 = quantile( data3, 0.75 );
qr = ( q75 - q25 ) .* 1.5;
m = mean( data3 );
lt = q25 - qr;
ut = q75 + qr;
for i = 1:size( data3, 2 );

    % the box
    h3 = [ h3 fill( [ i - w, i - w, i + w, i + w ] + w2, ...
        [ q25( i ), q75( i ), q75( i ), q25( i ) ], ...
        c3 ) ] ;

    % the median
    plot( [ i - w, i + w ] + w2, ...
        [ q50( i ), q50( i ) ], ...
        'k' );
%             'k', 'LineWidth', 2 );


    % compute the whiskers 
    lw = min( data3( data3( :, i ) > lt( i ), i ) );
    uw = max( data3( data3( :, i ) < ut( i ), i ) );

    % the lower whisker
    if ( numel( lw ) > 0 )
        plot( [ i, i ] + w2, ...
            [ q25( i ), lw ], ...
            'k' );

        plot( [ i - w * s, i + w * s ] + w2, ...
            [ lw, lw ], ...
            'k' );
    end

    % the upper whisker
    if ( numel( uw ) > 0 )
        plot( [ i, i ] + w2, ...
            [ q75( i ), uw ], ...
            'k' );

        plot( [ i - w * s, i + w * s ] + w2, ...
            [ uw, uw ], ...
            'k' );
    end

    % the mean
    plot( i + w2, m( i ), 'k*' );

    % compute the outliers
    idx = find( data3( :, i ) < lt( i ) );
    for j = 1:numel( idx )
       plot( i + w2, data3( idx( j ), i ), 'o', 'Color', C3 );
    end
    idx = find( data3( :, i ) > ut( i ) );
    for j = 1:numel( idx )
       plot( i + w2, data3( idx( j ), i ), 'o', 'Color', C3 );
    end

%         % plot the data
%         for j = 1:size( data, 1 )
%             plot( i + w2, data( j, i ), '.r' );
%         end
end
%h3hatchfill( h3, 'single' );
hold off;

bb = axis();
hold on
h4 = plot( -1, -1, 'ko' );
h5 = plot( -1, -1, 'k*' );
lgnd{ 1 } = labels{ 1 };
lgnd{ 2 } = labels{ 2 };
lgnd{ 3 } = labels{ 3 };
lgnd{ 4 } = 'Outlier';
lgnd{ 5 } = 'Mean';
hold off
axis( bb );
l = legend( [ h1( 1 ), h2( 1 ), h3( 1 ), h4, h5 ], lgnd, ...
    'Location', 'NorthOutside', 'Orientation', 'Horizontal' );
axis( [ 0.5, size( data1, 2 ) + 0.5, 0, 1 ] );
end

