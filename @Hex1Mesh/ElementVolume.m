function V = ElementVolume(T,iel)
% Efficient Computation of Volume of Hexhedral Cells by J.Grandy
node = T.Connectivity(iel,:);
X = T.X;
v11 = X(node(6),:)-X(node(2),:);
v12 = X(node(1),:)-X(node(2),:);
v13 = X(node(4),:)-X(node(5),:);

v21 = v11;
v22 = X(node(8),:)-X(node(2),:);
v23 = X(node(5),:)-X(node(7),:);

v31 = v11;
v32 = X(node(3),:)-X(node(2),:);
v33 = X(node(7),:)-X(node(4),:);

V = 1/6*(det([v11',v12',v13'])+det([v21',v22',v23'])+det([v31',v32',v33']));
end