function [] = plotOptimizationPath(x,functionID)

lb=-5;
up=5;
xi=lb:0.1:up;
f=zeros(length(xi),length(xi));
for i=1:length(xi)
    for j=1:length(xi)
        f(j,i)=getObjFVal([xi(i);xi(j)],functionID);
    end
end
figure('name','Optimization path')
hold on
axis equal
xlabel('x_1')
ylabel('x_2')
axis([lb up lb up])
title('Optimization path')
switch functionID
    case 1
        [C,h]=contour(xi,xi,f,[0:1:2 2:6:20 20:20:80 80:40:200 200:100:2000]);
        clabel(C,h);
    case 2
        [C,h]=contour(xi,xi,f,([-10:2:2 2:4:100]));
        clabel(C,h);
end
hold on

ind=1;
for i=1:size(x,2)-1
    plot(x(1,i),x(2,i),'.c','markersize',30)
    plot([x(1,i) x(1,i+1)],[x(2,i) x(2,i+1)],'c','linewidth',2)
    text(x(1,i),x(2,i),num2str(ind-1),'horizontalalignment','center','verticalalignment','middle')
    ind=ind+1;
end
plot(x(1,end),x(2,end),'.c','markersize',30)
text(x(1,end),x(2,end),num2str(ind-1),'horizontalalignment','center','verticalalignment','middle')

end