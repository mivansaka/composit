format long
%% Angle of Twist loading case two%%
Delete=[0];
q=3e6;
P=-2.5e4;
R=25e-3; %%Modify later%%
thick=2.5e-4;
L=0.3;
E1=236e9;
E2=5e9;
G=2.6e9;
v=0.25;
Xt=3800e6;
Xc=689e6;
Yt=41e6;
Yc=107e6;
S=69e6;
Z=[0;0;0]
TAT=zeros(151);
for i=-75:1:75
	Alpha=i
	AlphaF=i+76
	for i=-75:1:75
		Beta=i;
		BetaF=i+76;
		Theta=[Alpha;Beta;Alpha;Beta];
		t=[thick;thick;thick;thick];
		ctraQ=1-((v^2)*(E2/E1));
		Q=[E1/ctraQ v*E2/ctraQ 0; v*E2/ctraQ E2/ctraQ 0; 0 0 G];
		TT=[0;0;0];
		for i = 1:4
		straTAlpha=sind(Theta(i));
		ctraTAlpha=cosd(Theta(i));
		TAlpha=[ctraTAlpha^2 straTAlpha^2 -2*ctraTAlpha*straTAlpha;
		straTAlpha^2 ctraTAlpha^2 2*ctraTAlpha*straTAlpha;
		ctraTAlpha*straTAlpha -ctraTAlpha*straTAlpha ctraTAlpha^2-straTAlpha^2];
		TT=[TT TAlpha];
		eval(['TT',num2str(i),'=TAlpha;'])
				end
		TT=TT(:, 2:13);
		A=zeros(3);
		for i=1:4
		m=(i-1)*3+1;
		n=i*3;
		TTa=TT(:,m:n);
		A=A+thick*TTa*Q*TTa';
		Delete=Delete';
		end
		a=A^(-1);
		Q0=Q;
		N=[0.5*q*R; q*R; 0];
		Epsilonx=a*N;
		Phi=Epsilonx(3)*L/R*180/pi;
		TAT(AlphaF,BetaF)=Phi;

		for i=1:4
		m=(i-1)*3+1;
		n=i*3;
		TTa=TT(:,m:n);
		Epsilon= TTa'*Epsilonx;
		Delete=Delete';
		Sigma=Q0*Epsilon;
			if Sigma(1)>=0
			w1=Sigma(1)/Xt;
			else 
			w1=-Sigma(1)/Xc;
			end
			if Sigma(2)>=0
			w2=Sigma(2)/Yt;
			else
			w2=-Sigma(2)/Yc;
			end
			w3=abs(Sigma(3)/S);
			w=[w1; w2; w3];
		eval(['wcone',num2str(i),'=w;'])
		end
		wwassembly=[wcone1 wcone2 wcone3 wcone4];
		Fail_case_one=all(wwassembly(:)<=1);
		TATT(AlphaF,BetaF)=Fail_case_one;
		TATT=TATT*1;
		end
	end
	XX=-75:75;
	YY=-75:75;
	ZZ=TAT;
	mesh(XX,YY,ZZ)
	contour(XX,YY,ZZ,20,'ShowText','on')
	contour(XX,YY,TATT,1,'ShowText','on')