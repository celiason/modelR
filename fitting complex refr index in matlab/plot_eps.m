% plot comparison between fitting and tabular data
function plot_eps(x,eps_re_data,eps_im_data,eps_re_fit,eps_im_fit,arg,plot_title)

[x_nrow,x_ncol]=size(x);

x1=x(1);
x2=x(x_nrow);

xlabelsize=10;
ylabelsize=10;
axsz=8;
tsz=10;

subplot(3,1,1);
plot(x,eps_re_data,x,eps_re_fit);
title(plot_title,'FontSize',tsz);
legend('Data','Fitting');
xlabel(arg,'FontSize',xlabelsize);
ylabel('Re(eps)','FontSize',ylabelsize);
fig=get(0,'CurrentFigure');
ax=get(fig,'CurrentAxes');
set(ax, 'FontSize', axsz);
xlim([x1,x2]);

subplot(3,1,2);
plot(x,eps_im_data,x,eps_im_fit);
legend('Data','Fitting');
xlabel(arg,'FontSize',xlabelsize);
ylabel('Im(eps)','FontSize',ylabelsize);
fig=get(0,'CurrentFigure');
ax=get(fig,'CurrentAxes');
set(ax, 'FontSize', axsz);
xlim([x1,x2]);

subplot(3,1,3);
plot(x,abs(eps_im_data),x,abs(eps_im_fit));
legend('Data','Fitting');
xlabel(arg,'FontSize',xlabelsize);
ylabel('Im(eps)','FontSize',ylabelsize);
fig=get(0,'CurrentFigure');
ax=get(fig,'CurrentAxes');
set(ax, 'FontSize', axsz);
xlim([x1,x2]);
set(ax, 'yScale', 'log');

fig=get(0,'CurrentFigure');
set(fig,'Position',[300 100 400 800]);
