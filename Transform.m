function [ A0toj,Ajtoi,Ahatj ] = Transform( th,PARAM )
    global n_d
    eval(sprintf('%s=%f;',PARAM{:}));
    eval(sprintf('th%d=th(%d);',[1:n_d;1:n_d]));
    A0toj(:,:,1,1)=[cos(th1),-sin(th1),0,0;sin(th1),cos(th1),0,0;0,0,1,0;0,0,0,1];
    Ajtoi(:,:,1,1,1)=eye(4);
%     Ahatj(:,:,1)=[0,-sign(th1),0,0;sign(th1),0,0,0;0,0,0,0;0,0,0,0];
    A0toj(:,:,1,2)=[cos(th2),-sin(th2),0,x2;sin(th2),cos(th2),...
        0,y2;0,0,1,0;0,0,0,1];
    Ajtoi(:,:,1,1,2)=eye(4);
%     Ahatj(:,:,2)=[0,-sign(th2),0,0;sign(th2),0,0,0;0,0,0,0;0,0,0,0];
    A0toj(:,:,1,3)=A0toj(:,:,1,2);
    Ajtoi(:,:,1,2,3)=[cos(th3),-sin(th3),0,L23;sin(th3),cos(th3),...
        0,0;0,0,1,0;0,0,0,1];
    A0toj(:,:,2,3)=A0toj(:,:,1,3)*Ajtoi(:,:,1,2,3);
    Ajtoi(:,:,1,1,3)=eye(4);
    Ajtoi(:,:,2,2,3)=eye(4); 
%     Ahatj(:,:,3)=[0,-sign(th3),0,0;sign(th3),0,0,0;0,0,0,0;0,0,0,0];
    A0toj(:,:,1,4)=[cos(th4),-sin(th4),0,x4;sin(th4),cos(th4),...
        0,y4;0,0,1,0;0,0,0,1];
    Ajtoi(:,:,1,1,4)=eye(4);
%     Ahatj(:,:,4)=[0,-sign(th4),0,0;sign(th4),0,0,0;0,0,0,0;0,0,0,0];
    A0toj(:,:,1,5)=A0toj(:,:,1,4);
    Ajtoi(:,:,1,2,5)=[cos(th5),-sin(th5),0,L45;sin(th5),cos(th5),...
        0,0;0,0,1,0;0,0,0,1];
    A0toj(:,:,2,5)=A0toj(:,:,1,5)*Ajtoi(:,:,1,2,5);
    Ajtoi(:,:,1,1,5)=eye(4);
    Ajtoi(:,:,2,2,5)=eye(4);
%     Ahatj(:,:,5)=[0,-sign(th5),0,0;sign(th5),0,0,0;0,0,0,0;0,0,0,0];
    A0toj(:,:,1,6)=A0toj(:,:,1,4);
    Ajtoi(:,:,1,2,6)=[cos(th6),-sin(th6),0,L56+L45;sin(th6),cos(th6),...
        0,0;0,0,1,0;0,0,0,1];
    A0toj(:,:,2,6)=A0toj(:,:,1,6)*Ajtoi(:,:,1,2,6);
    Ajtoi(:,:,1,1,6)=eye(4);
    Ajtoi(:,:,2,2,6)=eye(4);
%     Ahatj(:,:,6)=[0,-sign(th6),0,0;sign(th6),0,0,0;0,0,0,0;0,0,0,0];
    A0toj(:,:,1,7)=A0toj(:,:,1,6);
    A0toj(:,:,2,7)=A0toj(:,:,2,6);
    Ajtoi(:,:,2,3,7)=[cos(th7),-sin(th7),0,L67;sin(th7),cos(th7),...
        0,0;0,0,1,0;0,0,0,1];
    A0toj(:,:,3,7)=A0toj(:,:,2,7)*Ajtoi(:,:,2,3,7);
    Ajtoi(:,:,1,1,7)=eye(4);
    Ajtoi(:,:,2,2,7)=eye(4);
    Ajtoi(:,:,3,3,7)=eye(4);
%     Ahatj(:,:,7)=[0,-sign(th7),0,0;sign(th7),0,0,0;0,0,0,0;0,0,0,0];
    Ahatj=[0,-1,0,0;1,0,0,0;0,0,0,0;0,0,0,0];
end

