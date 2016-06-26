function cjrsBoxPlot( data )

q25 = quantile( data, 0.25 );
q50 = quantile( data, 0.50 );
q75 = quantile( data, 0.75 );
qr = ( q75 - q25 ) .* 1.5;
m = mean( data );
lt = q25 - qr;
ut = q75 + qr;


w2 = 0.75/2;
s = 0.6;
hold on;
    for i = 1:size( data, 2 );

        
        % the box
        fill( [ i - w2, i - w2, i + w2, i + w2 ], ...
            [ q25( i ), q75( i ), q75( i ), q25( i ) ], ...
            [ 0.95 0.95 0.95 ] );
        
        % the median
        plot( [ i - w2, i + w2 ], ...
            [ q50( i ), q50( i ) ], ...
            'k' );
%             'k', 'LineWidth', 2 );
        
        % compute the whiskers 
        lw = min( data( data( :, i ) > lt( i ), i ) );
        uw = max( data( data( :, i ) < ut( i ), i ) );
        
        % the lower whisker
        if ( numel( lw ) > 0 )
            plot( [ i, i ], ...
                [ q25( i ), lw ], ...
                'k' );

            plot( [ i - w2 * s, i + w2 * s ], ...
                [ lw, lw ], ...
                'k' );
        end
        
        % the upper whisker
        if ( numel( uw ) > 0 )
            plot( [ i, i ], ...
                [ q75( i ), uw ], ...
                'k' );

            plot( [ i - w2 * s, i + w2 * s ], ...
                [ uw, uw ], ...
                'k' );
        end
        
        % the mean
        plot( i, m( i ), 'k*' );
        
        % compute the outliers
        idx = find( data( :, i ) < lt( i ) );
        for j = 1:numel( idx )
           plot( i, data( idx( j ), i ), 'ko' );
        end
        idx = find( data( :, i ) > ut( i ) );
        for j = 1:numel( idx )
           plot( i, data( idx( j ), i ), 'ko' );
        end
        
%         % plot the data
%         for j = 1:size( data, 1 )
%             plot( i + 0.1, data( j, i ), '.r' );
%         end
    end
hold off;
axis( [ 0.5, size( data, 2 ) + 0.5, 0, 1 ] );
end

