tic;
clear all;
clc;
close all;
disp('Passive Mechanism');
global n_d
n_q=1; n_d=7; n_j=7; g=9.81; 
% sinus
% a1=27.06; b1=9.932; c1=-3.288; a2=31.08; b2=16.79;
% c2=-1.364; a3=23.92; b3=23.86; c3=0.3736; a4=8.795; 
% b4=38.86; c4=-3.286; a5=15.36; b5=31.21; c5=-4.452;
% a6=126.1; b6=0.1404; c6=3.046; a7=2.454; b7=54.49;
% c7=-7.564; a8=4.931; b8=46.66; c8=-8.572;
% Gaussian
a1=222.2; b1=0.5822; c1=0.07462; a2=4.864; b2=0.4984;
c2=0.00144; a3=121.5; b3=0.6325; c3=0.05247; c5=0.1196;
a4=-1885; b4=0.5458; c4=0.1249; a5=1772; b5=0.5371; 
param={'L45','L5C','L3C','L23','L38','L18','L56','L67','L7C','L6','x4','y4','x2','y2'};
param_eval={0.096,0.239,0.115,0.115,0.115,0.038,0.239,0.239,0.239,0.8365,0.038,0.41,-0.089,0};
PARAM=[param;param_eval];
G={'x2+L23*cos(th2)-L3C*cos(th2+th3)-x4-L45*cos(th4)-L5C*cos(th4+th5)';
    'y2+L23*sin(th2)-L3C*sin(th2+th3)-y4-L45*sin(th4)-L5C*sin(th4+th5)';
    'L56*cos(th4)+L67*cos(th4+th6)+L7C*cos(th4+th6+th7)-L5C*cos(th4+th5)';
    'L56*sin(th4)+L67*sin(th4+th6)+L7C*sin(th4+th6+th7)-L5C*sin(th4+th5)';
    'x2+L23*cos(th2)+L38*cos(th2+th3)-L18*cos(th1)';
    'y2+L23*sin(th2)+L38*sin(th2+th3)-L18*sin(th1)';
%     'y4+(L45+L56)*sin(th(4))+L6*sin(th(4)+th(6))';
    'q-th1'};
G=sym(G);
for ii=1:length(G)
    eval(sprintf('G(ii)=subs(G(ii),''%s'',%f); ',PARAM{:}));
end
step=0.025; eps1=eye(n_d)*1e-6; eps2=eye(n_q)*1e-6;  eps=1e-6; 
Q(1,1)=0*pi/180; Qdot(1,1)=4.7124; q_eval=0*pi/180; qdot_eval=4.7124; Ang=[];
eval(sprintf('%s=%f;',PARAM{:}));
% options=optimoptions('fsolve','TolFun',1e-13,'TolX',1e-10,'MaxFunEvals',500*n_d,'MaxIter',1000);
options=optimoptions('fmincon','Algorithm','sqp','TolCon',1e-6,'MaxFunEvals',50*n_d,'MaxIter',50);
cc=0; init=[0,0.0011,-0.0023,0.003,-0.0046,-0.0029,0.0028]+1e-1;
%% gradiant projection:
for tt=0:step:5
    fprintf('time=%d\n',tt);
    for ll=1:4
        flag=-2;
        while flag~=1 && flag~=2
            [Sol,~,flag]=fmincon(@(x)0,init,[],[],[],[],[],[],@(x)NLC(x,q_eval,PARAM),options);
%             [Sol,~,flag]=fsolve(@(x)NLC(x,q_eval,PARAM),init);
            init=randn(7,1); 
        end
        init=Sol;
        [~,fval]=NLC(Sol,q_eval,PARAM);
        ss=sprintf('%d,',fval);
        fprintf('fval=[%s] && exitflag=%i\n',ss,flag);
        if (flag~=1 && flag~=2)
            error('flag~=1');
        end
%         ss=[num2cell(1:n_j);angle];
        ss=[num2cell(1:n_j);num2cell(1:n_j)];
%         eval(sprintf('angle_eval(%d)=double(Sol.%s);',ss{:}));
        eval(sprintf('angle_eval(%d)=double(Sol(%d));',ss{:}));
        while any(angle_eval>pi+1e-1)
            angle_eval(angle_eval>pi+1e-1)=angle_eval(angle_eval>pi+1e-1)-2*pi;
        end
        while any(angle_eval<-pi-1e-1)
            angle_eval(angle_eval<-pi-1e-1)=angle_eval(angle_eval<-pi-1e-1)+2*pi;
        end
        if ll==1
            Ang=[Ang;angle_eval];
        end
        eval(sprintf('th(%d)=angle_eval(%d);',[1:n_d;1:n_d]));
        for ii=1:n_d
            for jj=1:n_d
                var1=Runge_Kutta(G(ii),[angle_eval+eps1(jj,:),q_eval]);
                var2=Runge_Kutta(G(ii),[angle_eval,q_eval]);
                ghat(ii,jj)=(var1-var2)/eps;
                dghattoq=zeros(7,7);
%                 for kk=1:n_q
%                     var1=Runge_Kutta(G(ii),[angle_eval+eps1(jj,:),q_eval+eps2(kk,:)]);
%                     var2=Runge_Kutta(G(ii),[angle_eval+eps1(jj,:),q_eval]);
%                     var3=Runge_Kutta(G(ii),[angle_eval,q_eval+eps2(kk,:)]);
%                     var4=Runge_Kutta(G(ii),[angle_eval,q_eval]);
%                     dghattoq((ii-1)*n_q+kk,jj)=(var1-var2-var3+var4)/eps^2;
%                 end
            end
        end
        for ii=1:n_d
            for jj=1:n_q
                var1=Runge_Kutta(G(ii),[angle_eval,q_eval+eps2(jj,:)]);
                var2=Runge_Kutta(G(ii),[angle_eval,q_eval]);
                dGtoq((ii-1)*n_q+jj,1)=(var1-var2)/eps;
                for kk=1:n_q
                    var1=Runge_Kutta(G(ii),[angle_eval,q_eval+eps2(jj,:)+eps2(kk,:)]);
                    var2=Runge_Kutta(G(ii),[angle_eval,q_eval+eps2(jj,:)]);
                    var3=Runge_Kutta(G(ii),[angle_eval,q_eval+eps2(kk,:)]);
                    var4=Runge_Kutta(G(ii),[angle_eval,q_eval]);
                    ddGtoq((ii-1)*n_q^2+(jj-1)*n_q+kk,1)=(var1-var2-var3+var4)/eps^2;
                end
            end
        end
        gradf=Prediction(ghat,dGtoq,n_q);
        var=zeros(n_d*n_q^2,1);
        for ii=1:n_d  
            var=var+Diadic(eye(n_d),FU(n_q,n_q))*Diadic(dghattoq(:,ii),gradf(:,ii));
        end
        hessian=Prediction(ghat,ddGtoq+var,n_q^2);
        %% Transformations: 
        Len={1,2,[2,3],4,[4,5],[4,6],[4,6,7]};
        [A0toj,Ajtoi,Ahatj]=Transform(th,PARAM);
        %% DBtoq:
        for ii=1:n_j
            for jj=length(Len{ii})
                dbtoq=zeros(4*n_q,4);
                for kk=1:jj
                    W(:,:,kk,jj,ii)=A0toj(:,:,kk,ii)*Ahatj*Ajtoi(:,:,kk,jj,ii);
                    dbtoq=dbtoq+Diadic(W(:,:,kk,jj,ii)',gradf(:,Len{ii}(kk)));
                end
            dBtoq(:,:,jj,ii)=dbtoq;
            end
        end
        %% DDBtoq:
        for ii=1:n_d
            for kk=length(Len{ii})
%                 var=zeros(4*n_q,4);
                ddbtoq=zeros(4*n_q^2,4);
                for jj=1:kk
                    var=zeros(4*n_q,4);
                    ddbtoq=ddbtoq+Diadic((W(:,:,jj,kk,ii))',hessian(:,Len{ii}(jj)));
                    for ss=1:kk  
                        if (ss<=jj)
                            What=W(:,:,ss,jj,ii)*Ahatj*Ajtoi(:,:,jj,kk,ii);
                            var=var+Diadic(What',gradf(:,Len{ii}(jj)));
                        elseif (ss>jj)
                            What=W(:,:,jj,ss,ii)*Ahatj*Ajtoi(:,:,ss,kk,ii);
                            var=var+Diadic(What',gradf(:,Len{ii}(ss)));
                        end
%                         ddbtoq=ddbtoq+Diadic(What',Diadic(gradf(:,Len{ii}(jj)),gradf(:,Len{ii}(ss))));
                    end
                    ddbtoq=ddbtoq+Diadic(eye(4),FU(n_q,n_q))*Diadic(var,gradf(:,Len{ii}(jj)));
                end
                ddBtoq(:,:,kk,ii)=ddbtoq;
            end
        end
        %% Component of dynamic equation
        m_c=[1.5,3,5,4.2,3,10.5,3];
        Ip(:,:,1)=[0,0,0;0,0,0;0,0,0.722e-3]; Ip(:,:,2)=[0,0,0;0,0,0;0,0,13.271e-3];
        Ip(:,:,3)=[0,0,0;0,0,0;0,0,2.21e-2]; Ip(:,:,4)=[0,0,0;0,0,0;0,0,156.74e-3];
        Ip(:,:,5)=[0,0,0;0,0,0;0,0,57.12e-3]; Ip(:,:,6)=[0,0,0;0,0,0;0,0,2.4491];
        Ip(:,:,7)=[0,0,0;0,0,0;0,0,57.12e-3]; 
        for ii=1:n_j
            Ipp(:,:,ii)=-(Ip(:,:,ii)-0.5*trace(Ip(:,:,ii))*eye(3));
        end
        C1=Diadic(Vec((eye(4))),(eye(n_q)));
        C2=Diadic((eye(4)),FU(n_q,4))*C1;
        C3=Diadic(Diadic(eye(4),qdot_eval),FU(n_q,4))*C1;
        rc_c=[[0.019;0;0],[0.0576;0;0],[0;0;0],[0.1673;0;0],[0.1195;0;0],[0.4183;0;0],[0.1195;0;0]];
        for ii=1:n_j
            J(:,:,ii)=[Ipp(:,:,ii),m_c(ii)*rc_c(:,ii);m_c(ii)*rc_c(:,ii)',m_c(ii)];
            for jj=length(Len{ii})
                M(:,:,jj,ii)=C1'*Diadic(J(:,:,ii),FU(n_q,4)')*Diadic(dBtoq(:,:,jj,ii)*dBtoq(:,:,jj,ii)',eye(4))*C2;
                N(:,:,jj,ii)=C1'*Diadic(J(:,:,ii),FU(n_q,4)')*Diadic(dBtoq(:,:,jj,ii)*ddBtoq(:,:,jj,ii)',eye(4))*C3;
                GG(:,jj,ii)=m_c(ii)*Diadic([rc_c(:,ii)',1],eye(n_q))*dBtoq(:,:,jj,ii)*[0;-g;0;0];
            end
        end
        N_total=zeros(n_q,n_q);
        M_total=zeros(n_q,n_q);
        GG_total=zeros(n_q,1);
        for ii=1:n_j
            finish=length(Len{ii});
%             for jj=1:finish
%                 M_total=M_total+M(:,:,jj,ii);
                M_total=M_total+M(:,:,finish,ii);
%                 N_total=N_total+N(:,:,jj,ii);
                N_total=N_total+N(:,:,finish,ii);
%                 GG_total=GG_total+GG(:,jj,ii);
                GG_total=GG_total+GG(:,finish,ii);
%             end
        end
        %% Matrix of arm and moment arm of external force
        for ii=1:n_j
            finish=length(Len{ii});
            gama(1:n_q,1:3,1:finish,ii)=zeros(n_q,3,finish);
            beta(1:n_q,1:3,1:finish,ii)=zeros(n_q,3,finish);
            for jj=1:finish
                for kk=1:jj
                    rltoj=Ajtoi(1:3,4,kk,jj,ii);
                    Vteta=Ahatj(1:3,1:3);
                    HH=[Vteta(3,2);Vteta(1,3);Vteta(2,1)]; 
                        gama(:,:,jj,ii)=gama(:,:,jj,ii)+Diadic(rltoj'*Vteta'*Ajtoi(1:3,1:3,kk,jj,ii),gradf(:,Len{ii}(kk)));
                    if kk==1
                        beta(:,:,jj,ii)=beta(:,:,jj,ii)+gradf(:,Len{ii}(kk))*(A0toj(1:3,1:3,jj,ii)'*HH)';
                    else
                        beta(:,:,jj,ii)=beta(:,:,jj,ii)+gradf(:,Len{ii}(kk))*(Ajtoi(1:3,1:3,kk-1,jj,ii)'*HH)';
                    end
                end 
            end
        end
        M_ext=0;
%         str=sprintf('a%d*sin(b%d*tt+c%d)*heaviside(tt-.4)*heaviside(.8-tt)+',repmat(1:7,3,1));
%         str=sprintf('a%d*exp(-((tt-b%d)/c%d)^2)+',repmat(1:5,3,1));
%         M_ext=eval(str(1:end-1));
        T_ext=beta(:,:,1,1)*[0;0;1]*M_ext;
        %% Forward Dynamic  
        if ll==1
            cc=cc+1;
            fprintf('cc=%d\n',cc);
%             h=PlotStructure(A0toj,PARAM,cc);
%             title(sprintf('Time=%f',tt));
%             xlabel('X');
%             ylabel('Y');
%             figure(cc);
%             set(gcf, 'PaperPosition', [0 0 4 2]);
%             subplot(1,2,2);
%             plot(0:step:tt,Q,'k');
%             xlabel('Time(s)');
%             ylabel('Generalized coordinate(rad)');
%             xlim([0,5]);
%             ylim([-1 3.5]);
%             axis square;
%             cd('D:\Scrapping\Figures2');
%             print(h,'-dtiff',sprintf('Fig%d',cc));
%             pause(0.001);
%             subplot(1,2,1);
%             cla;
            K1_q=qdot_eval*step; K1_qd=M_total\(T_ext+GG_total-N_total*qdot_eval)*step;
%             q_eval=Q(fix(tt/step)+1,:)'+K1_q/2; qdot_eval=Qdot(fix(tt/step)+1,:)'+K1_qd/2;
            q_eval=Q(end,:)'+K1_q/2; qdot_eval=Qdot(end,:)'+K1_qd/2;
        elseif ll==2
            K2_q=qdot_eval*step; K2_qd=M_total\(T_ext+GG_total-N_total*qdot_eval)*step;
            q_eval=Q(end,:)'+K2_q/2; qdot_eval=Qdot(end,:)'+K2_qd/2;
        elseif ll==3
            K3_q=qdot_eval*step; K3_qd=M_total\(T_ext+GG_total-N_total*qdot_eval)*step;
            q_eval=Q(end,:)'+K3_q; qdot_eval=Qdot(end,:)'+K3_qd;
        else
            K4_q=qdot_eval*step; K4_qd=M_total\(T_ext+GG_total-N_total*qdot_eval)*step;
        end
    end
    q_eval=Q(end,:)'+(K1_q+2*K2_q+2*K3_q+K4_q)/6;
    qdot_eval=Qdot(end,:)'+(K1_qd+2*K2_qd+2*K3_qd+K4_qd)/6;
    fprintf('Q =');
    fprintf('%d',Q(end,:));
    fprintf('Qdot =');
    fprintf('%d',Qdot(end,:));
%     if q_eval>pi
%         q_eval=q_eval-2*pi;
%     elseif q_eval<-pi
%         q_eval=q_eval+2*pi;
%     end
    Q(end+1,:)=q_eval; Qdot(end+1,:)=qdot_eval;
end
toc;