function  [C,Ceq]  = NLC( th,q,PARAM )
    global n_d
    eval(sprintf('%s=%d;',PARAM{:}));
    G={'x2+L23*cos(th(2))-L3C*cos(th(2)+th(3))-x4-L45*cos(th(4))-L5C*cos(th(4)+th(5))';
    'y2+L23*sin(th(2))-L3C*sin(th(2)+th(3))-y4-L45*sin(th(4))-L5C*sin(th(4)+th(5))';
    'L56*cos(th(4))+L67*cos(th(4)+th(6))+L7C*cos(th(4)+th(6)+th(7))-L5C*cos(th(4)+th(5))';
    'L56*sin(th(4))+L67*sin(th(4)+th(6))+L7C*sin(th(4)+th(6)+th(7))-L5C*sin(th(4)+th(5))';
    'x2+L23*cos(th(2))+L38*cos(th(2)+th(3))-L18*cos(th(1))';
    'y2+L23*sin(th(2))+L38*sin(th(2)+th(3))-L18*sin(th(1))';
%     'y4+(L45+L56)*sin(th(4))+L6*sin(th(4)+th(6))';
    'q-th(1)'};
    for ii=1:length(G)
        F(ii,1)=eval(G{ii});
    end
    C=[];
    Ceq=F;
%     eval(sprintf('C=[C;th(%d)-2*pi-1e-1];',1:n_d));
    eval(sprintf('C=[C;pi+1e-1-th(%d)];',3:n_d));
end