%%Term CASE1FM loading case 1 Failure matrix%%
%%Term CASE2FM loading case 2 Failure matrix%%
%%Term AssemblyFM Assembled Falure matrix [-30,30] has been removed%%
%%Term CASE2TA loading case 2 twist angle%%
%%Term BP Best point%%

format short %%Modify according to the requirement%%

%% Angle of Twist loading case one%%

resolution=0.1; %%Modify according to the requirement%%
Delete=[0]; %%leave it alone%%
q=3e6;
P=-2.5e4;
R=25e-3; 
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
TAT=zeros(150/resolution+1); %%Initial matrix of Twist Angle%%

for i=-75/resolution:1:75/resolution %%2 for loots control the Alpha and Beta%%
	Alpha=i*resolution
	AlphaF=i+75/resolution+1;
	for i=-75/resolution:1:75/resolution
		Beta=i*resolution;
		BetaF=i+75/resolution+1;
		Theta=[Alpha;Beta;Alpha;Beta];
		t=[thick;thick;thick;thick];
		ctraQ=1-((v^2)*(E2/E1));
		Q=[E1/ctraQ v*E2/ctraQ 0; v*E2/ctraQ E2/ctraQ 0; 0 0 G];
		TT=[0;0;0];
			for i = 1:4 %%Find related equition in sheet%%
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
CASE1FM=TATT;

%% Angle of Twist loading case two%%

TAT=zeros(150/resolution+1);

for i=-75/resolution:1:75/resolution
	Alpha=i*resolution
	AlphaF=i+75/resolution+1;
	for i=-75/resolution:1:75/resolution
		Beta=i*resolution;
		BetaF=i+75/resolution+1;
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
		N=[P/(2*pi*R); 0; 0];
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
XX=-75/resolution:75/resolution;
YY=-75/resolution:75/resolution;
CASE2FM=TATT;
CASE2TA=TAT;

%%Assemble the Failure matrix%%

AssemblyFM=zeros(150/resolution+1);
	for i=-75/resolution:1:75/resolution
		i=i+75/resolution+1;
		for j=-75/resolution:1:75/resolution
			j=j+75/resolution+1;
			m=CASE2FM(i,j);
			n=CASE1FM(i,j);
			if m==1 && n==1
				AssemblyFM(i,j)=1;
			else 
				AssemblyFM(i,j)=0;
			end
		end
	end
AssemblyFM((-30+75+1)/resolution:(30+75+1)/resolution,(-30+75+1)/resolution:(30+75+1)/resolution)=0;

%%Find the best point%%

BPX=0;
BPY=0;
BPZ=0;
for i=-75/resolution:1:75/resolution
	i=i+75/resolution+1;
		for j=-75/resolution:1:75/resolution
			j=j+75/resolution+1;
			k=AssemblyFM(i,j);
			BPZS=CASE2TA(i,j);
				if k==1 && BPZS>BPZ
					BPX=(i-75/resolution-1)*resolution;
					BPY=(j-75/resolution-1)*resolution;
					BPZ=BPZS;
				else 
					BPZ=BPZ;
				end
		end
end
BP=[BPX,BPY,BPZ]

%%Figure commands for observation and later using%%
%%mesh(XX,YY,CASE2TA)%%
%%contour(XX,YY,CASE1FM,2,'ShowText','on')%%
%%contour(XX,YY,CASE2FM,2,'ShowText','on')%%
%%contour(XX,YY,AssemblyFM,2)%%
%%mesh(XX,YY,AssemblyFM)%%

