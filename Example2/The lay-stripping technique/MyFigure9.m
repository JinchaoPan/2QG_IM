%% Draw a Error variation graph with respect to the iteration time (Figure 9)
% if the error value is infinite, replace it with e^{450}.
for Tki=1:1:L
    QQ1=ERROR_ABSQ1(Tki,:);
    QQ2=ERROR_ABSQ2(Tki,:);

    pp1=ERROR_ABSP1(Tki,:);
    pp2=ERROR_ABSP2(Tki,:);

    kn=size(QQ1,2);
    Kqk=1:kn;

    QQPP=[QQ1;QQ2;pp1;pp2];
    akak=isfinite(QQPP);
    if any(~akak(:))
        [~,kn]=find(akak==0);kn=kn(1);
        Kqk=1:kn;
        QQ1=QQ1(1:kn); QQ1(kn)=exp(450);
        QQ2=QQ2(1:kn); QQ2(kn)=exp(450);
        pp1=pp1(1:kn); pp1(kn)=exp(450);
        pp2=pp2(1:kn); pp2(kn)=exp(450);
    end

    LnQ1=log(QQ1);
    LnQ2=log(QQ2);
    Lnp1=log(pp1);
    Lnp2=log(pp2);

    figure(11+Tki-1)
    plot(Kqk,LnQ1,'r-^','LineWidth',2,'MarkerSize',8)
    hold on
    plot(Kqk,LnQ2,'r--^','LineWidth',2,'MarkerSize',8)
    hold on
    plot(Kqk,Lnp1,'b-o','LineWidth',2,'MarkerSize',8)
    hold on
    plot(Kqk,Lnp2,'b--o','LineWidth',2,'MarkerSize',8)
    hold off

    % legend('ErQ1(n)','ErQ2(n)','ErP1(n)','ErP2(n)','FontSize',10)
    legend('ln(ErQ1(n))','ln(ErQ2(n))','ln(ErP1(n))','ln(ErP2(n))','FontSize',10, 'Location', 'northeast')
    xlabel('Iteration times n','FontSize',15)
    % ylabel('E','FontSize',15)
    % xticks([1 3 5 7 9 12]);
    xlim([1,kn])
    xticks(Kqk);
    % ylim([-1,6]);
    % yticks([-50 0 50 100 150 200 250 300 350 400 450])
    % yticklabels({'-50','0','50','100','150','200','250','300','350','400','Inf'});
    if any(~akak(:))
        yticks([-50 0 50  150  250  350  450])
        yticklabels({'-50','0','50','150','250','350','Inf'});
        ylim([-50 450]);
    end
    % QP=[LnQ1 LnQ2 Lnp1 Lnp2];
    % cbI=[min(QP) max(QP)];
    % save 2CBI8 cbI
    % CBI=[-11.562663965250556,2.876044109730560];
    % CBI=[-12.562663965250556,2.876044109730560];
    % CBI=[-10,5.8];
    % ylim(CBI);
end