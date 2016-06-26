function testcode()
h.myfig = figure;

h.myax(1) = axes('Parent', h.myfig, 'Units', 'Normalized', 'Position', [0.1 0.1 0.35 0.35]);
h.myax(2) = axes('Parent', h.myfig, 'Units', 'Normalized', 'Position', [0.55 0.1 0.35 0.35]);
h.myax(3) = axes('Parent', h.myfig, 'Units', 'Normalized', 'Position', [0.55 0.55 0.35 0.35]);

plot(h.myax(1), 0:10);
plot(h.myax(2), 10:20);

hold(h.myax(3), 'on');
plot(h.myax(3), 20:30);
plot(h.myax(3), 30:-1:20);
hold(h.myax(3), 'off');

h.mybrush = brush;
set(h.mybrush, 'Enable', 'on', 'ActionPostCallback', @displayBrushData); 

