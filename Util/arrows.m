function [] = arrows(t1, t2,tend, ax)

Xlim=[t1 t2];
Ylim=[-3 3];

% arrow right 
P1 = [tend(2)/2+tend(2)*0.1 Ylim(1)+(Ylim(2)-Ylim(1))*0.09]; % from point
P2 = [tend(2) Ylim(1)+(Ylim(2)-Ylim(1))*0.09]; % to point

Pos = ax.Position;
X_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P1(1)-Xlim(1));
X_conv(2)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P2(1)-Xlim(1));
Y_conv(1)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P1(2)-Ylim(1));
Y_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P2(2)-Ylim(1));
an = annotation('textarrow', X_conv, Y_conv,'String',' act. ', 'HeadStyle','vback3','Interpreter', 'latex');
an.HeadLength = 5;
hold on

% arrow left 
P1 = [tend(2)/2-tend(2)*0.1 Ylim(1)+(Ylim(2)-Ylim(1))*0.09];
P2 = [0 Ylim(1)+(Ylim(2)-Ylim(1))*0.09];
X_conv(1)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P1(1)-Xlim(1));
X_conv(2)=Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P2(1)-Xlim(1));
Y_conv(1)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P1(2)-Ylim(1));
Y_conv(2)=Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P2(2)-Ylim(1));
an = annotation('textarrow', X_conv, Y_conv,'String',' ', 'HeadStyle','vback3');
an.HeadLength = 4;
an.HeadWidth = 7;

% vertical line at beginning and end of muscle activation
plot([0 0],[-5000 5000],'--k')
plot([tend(2) tend(2)],[-5000 5000],'--k')
end