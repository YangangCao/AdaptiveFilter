%plot out the different weights found through simulations
figure
plot_mvdr('run1'); hold on
plot_mvdr('run2');
plot_mvdr('run3');

hold off
title('figure 15.12')
print -dpsc 15_11