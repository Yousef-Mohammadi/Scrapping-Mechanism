function h = PlotStructure(A0toj,PARAM,cc)
    eval(sprintf('%s=%f;',PARAM{:}));
    h=figure(cc);
    set(gcf, 'PaperPosition', [0 0 4 2]);
    subplot(1,2,1);
%     axes('Position',[0.05 0.05 0.7 0.7]);
    hold on;
    plot([A0toj(1,4,1,1),A0toj(1,:,2,3)*[L38;0;0;1]],[A0toj(2,4,1,1),A0toj(2,:,2,3)*[L38;0;0;1]],...
        'k-o','linewidth',2,'markersize',5);
    plot([A0toj(1,4,1,3),A0toj(1,4,2,3)],[A0toj(2,4,1,3),A0toj(2,4,2,3)],'k-o','linewidth',2,...
        'markersize',5);
    plot([A0toj(1,:,2,3)*[-L3C;0;0;1],A0toj(1,:,2,3)*[L38;0;0;1]],...
        [A0toj(2,:,2,3)*[-L3C;0;0;1],A0toj(2,:,2,3)*[L38;0;0;1]],'k-o','linewidth',2,...
        'markersize',5);
    plot([A0toj(1,4,1,6),A0toj(1,4,2,6)],[A0toj(2,4,1,6),A0toj(2,4,2,6)],'k-o','linewidth',2,...
        'markersize',5);
    plot([A0toj(1,4,2,5),A0toj(1,:,2,3)*[-L3C;0;0;1]],[A0toj(2,4,2,5),A0toj(2,:,2,3)*[-L3C;0;0;1]],...
        'k-o','linewidth',2,'markersize',5);
    plot([A0toj(1,4,2,6),A0toj(1,:,2,6)*[L6;0;0;1]],[A0toj(2,4,2,6),A0toj(2,:,2,6)*[L6;0;0;1]],...
        'k-o','linewidth',2,'markersize',5);
    plot([A0toj(1,4,3,7),A0toj(1,:,2,3)*[-L3C;0;0;1]],[A0toj(2,4,3,7),A0toj(2,:,2,3)*[-L3C;0;0;1]],...
        'k-o','linewidth',2,'markersize',5); 
    xlim([-.8 0.35]);
    ylim([-0.5 0.5]);
    axis square
%     cla;
end

