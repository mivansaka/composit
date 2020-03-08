%%Colormap transformation%%

colormap(parula)

%%TA Figure%%

mesh(XX,YY,CASE2TA)
title('Twist Angle','Fontname', 'Times New Roman','FontSize',18);
x1=xlabel('Alpha (degree)','Fontname', 'Times New Roman','FontSize',18);
x2=ylabel('Beta (degree)','Fontname', 'Times New Roman','FontSize',18);
x3=zlabel('Twist angle (degree)','Fontname', 'Times New Roman','FontSize',18);
set(x1,'Rotation',17);
set(x2,'Rotation',-22);

saveas(gcf,'Twist Angle.jpg')

%%TA plane figure%%

contour(XX,YY,CASE2TA,10,'ShowText','on')
title('Twist Angle (degree)','Fontname', 'Times New Roman','FontSize',18);
x1=xlabel('Alpha (degree)','Fontname', 'Times New Roman','FontSize',18);
x2=ylabel('Beta (degree)','Fontname', 'Times New Roman','FontSize',18);

saveas(gcf,'Twist Angle (degree).jpg')

%%TA Figure without Failure%%

mesh(XX,YY,TAFM)
title('Twist Angle without Failure','Fontname', 'Times New Roman','FontSize',18);
x1=xlabel('Alpha (degree)','Fontname', 'Times New Roman','FontSize',18);
x2=ylabel('Beta (degree)','Fontname', 'Times New Roman','FontSize',18);
x3=zlabel('Twist angle (degree)','Fontname', 'Times New Roman','FontSize',18);
set(x1,'Rotation',17);
set(x2,'Rotation',-22);

saveas(gcf,'Twist Angle without Failure.jpg')

%%Colormap transformation%%

colormap(summer)

%%Devided failure analysis figure%%

gca=pcolor(XX,YY,CASE1FM)
set(gca, 'LineStyle','none');
title('Failure Analysis of Case One','Fontname', 'Times New Roman','FontSize',18);
x1=xlabel('Alpha (degree)','Fontname', 'Times New Roman','FontSize',18);
x2=ylabel('Beta (degree)','Fontname', 'Times New Roman','FontSize',18);

saveas(gcf,'Failure Analysis of Case One.jpg')

gca=pcolor(XX,YY,CASE2FM)
set(gca, 'LineStyle','none');
title('Failure Analysis of Case Two','Fontname', 'Times New Roman','FontSize',18);
x1=xlabel('Alpha (degree)','Fontname', 'Times New Roman','FontSize',18);
x2=ylabel('Beta (degree)','Fontname', 'Times New Roman','FontSize',18);

saveas(gcf,'Failure Analysis of Case Two.jpg')

%%Assembled failure analysis figure%%

pcolor(XX,YY,AssemblyFM)
title('Failure Analysis','Fontname', 'Times New Roman','FontSize',18);
x1=xlabel('Alpha (degree)','Fontname', 'Times New Roman','FontSize',18);
x2=ylabel('Beta (degree)','Fontname', 'Times New Roman','FontSize',18);

saveas(gcf,'Failure Analysis.jpg')