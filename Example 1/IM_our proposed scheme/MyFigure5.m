%% Draw a Error variation graph with respect to the iteration time(Figure 5)
% if the error value is infinite, replace it with e^{450}.
QQ1=ERROR_absQ1;
QQ2=ERROR_absQ2;
pp1=ERROR_absp1;
pp2=ERROR_absp2;

Kqk=size(QQ1,2);
Kqk=1:Kqk;

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

figure(11)
plot(Kqk,LnQ1,'r-^','LineWidth',2,'MarkerSize',8)
hold on
plot(Kqk,LnQ2,'r--^','LineWidth',2,'MarkerSize',8)
hold on
plot(Kqk,Lnp1,'b-o','LineWidth',2,'MarkerSize',8)
hold on
plot(Kqk,Lnp2,'b--o','LineWidth',2,'MarkerSize',8)
hold off

legend('ln(ErQ1(n))','ln(ErQ2(n))','ln(ErP1(n))','ln(ErP2(n))','FontSize',10, 'Location', 'northeast')
xlabel('Iteration times n','FontSize',15)
xticks(Kqk);
if any(~akak(:))
    yticks([-50 0 50  150  250  350  450])
    yticklabels({'-50','0','50','150','250','350','Inf'});
    ylim([-50 450]);
end

