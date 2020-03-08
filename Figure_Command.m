CASE1FM((-30+75)/resolution+1:(30+75)/resolution+1,:)=0;
CASE1FM(:,(-30+75)/resolution+1:(30+75)/resolution+1)=0;

CASE2FM((-30+75)/resolution+1:(30+75)/resolution+1,:)=0;
CASE2FM(:,(-30+75)/resolution+1:(30+75)/resolution+1)=0;

%%Colormap transformation%%

colormap(parula)

%%TA Figure%%

mesh(XX,YY,CASE2TA)
title('Twist Angle','Fontname', 'Times New Roman','FontSize',24);
x1=xlabel('Alpha (degree)','Fontname', 'Times New Roman','FontSize',24);
x2=ylabel('Beta (degree)','Fontname', 'Times New Roman','FontSize',24);
x3=zlabel('Twist angle (degree)','Fontname', 'Times New Roman','FontSize',24);
set(x1,'Rotation',15);
set(x2,'Rotation',-25);

saveas(gcf,'Twist Angle.png')

%%TA plane figure%%

contour(XX,YY,CASE2TA,10,'ShowText','on')
title('Twist Angle (degree)','Fontname', 'Times New Roman','FontSize',24);
x1=xlabel('Alpha (degree)','Fontname', 'Times New Roman','FontSize',24);
x2=ylabel('Beta (degree)','Fontname', 'Times New Roman','FontSize',24);

saveas(gcf,'Twist Angle (degree).jpg')

%%TA Figure without Failure%%

mesh(XX,YY,TAFM)
title('Twist Angle without Failure','Fontname', 'Times New Roman','FontSize',24);
x1=xlabel('Alpha (degree)','Fontname', 'Times New Roman','FontSize',24);
x2=ylabel('Beta (degree)','Fontname', 'Times New Roman','FontSize',24);
x3=zlabel('Twist angle (degree)','Fontname', 'Times New Roman','FontSize',24);
set(x1,'Rotation',15);
set(x2,'Rotation',-25);

saveas(gcf,'Twist Angle without Failure.jpg')

%%Colormap transformation%%

colormap(summer)

%%Devided failure analysis figure%%

gca = pcolor(XX,YY,CASE1FM);
set(gca, 'LineStyle','none');
title('Failure Analysis of Case One','Fontname', 'Times New Roman','FontSize',24);
x1=xlabel('Alpha (degree)','Fontname', 'Times New Roman','FontSize',24);
x2=ylabel('Beta (degree)','Fontname', 'Times New Roman','FontSize',24);

saveas(gcf,'Failure Analysis of Case One.jpg')

gca = pcolor(XX,YY,CASE2FM);
set(gca, 'LineStyle','none');
title('Failure Analysis of Case Two','Fontname', 'Times New Roman','FontSize',24);
x1=xlabel('Alpha (degree)','Fontname', 'Times New Roman','FontSize',24);
x2=ylabel('Beta (degree)','Fontname', 'Times New Roman','FontSize',24);

saveas(gcf,'Failure Analysis of Case Two.jpg')

%%Assembled failure analysis figure%%

gca = pcolor(XX,YY,AssemblyFM);
set(gca, 'LineStyle','none');
title('Failure Analysis','Fontname', 'Times New Roman','FontSize',24);
x1=xlabel('Alpha (degree)','Fontname', 'Times New Roman','FontSize',24);
x2=ylabel('Beta (degree)','Fontname', 'Times New Roman','FontSize',24);

saveas(gcf,'Failure Analysis.jpg')





colormap(parula)
colormap(jet)
colormap(gray)
colormap(summer)


