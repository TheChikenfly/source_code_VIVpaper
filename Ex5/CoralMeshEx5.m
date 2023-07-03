%create Coral Mesh
numElements = 50;
l = 0.15; % 15 cm % height of the coral 
Nlevel = 1;
if (Nlevel == 1);numElemCentral = numElements; else; numElemCentral = floor(numElements/(Nlevel+1)); end
numElemBranch = numElements/2;
%NElemLevel = numElemCentral+2*numElemBranch;
NElemLevel = numElemCentral+2*numElemBranch;
numElements = NElemLevel*Nlevel + numElemCentral % + last trunk
% Angle of branches
angle = 60; % degrees
alpha = pi*angle/180; % radians
% Length of the branches
lb = 0.08; % 4 cm for real corals
% Space between the branches
L1 = l/(Nlevel+1); % The levels compress to keep l as the height
%L2 = l/sqrt(2); % Length of the branches
%L2 = l/(sqrt(2)*numElemBranch*(Nlevel+1)); % Length of the branches: they scale with L1!
L2 = lb/(numElemBranch); % Length of the branches: they scale with L1!
if Nlevel == 0;
    mesh.nodesCoords = [  zeros(numElements+1,1) zeros(numElements+1,1)  (0:(numElements))'*l/numElements] ;
    %mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of nodes that compose the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
    mesh.conecCell = { } ;
    %md then the first welded node is defined with material (M) zero since nodes don't have material, the first element (E) type (the first entry of the `elements` struct), and (B) is the first entry of the the `boundaryConds` struct. For (I) no non-homogeneous initial condition is considered (then zero is used) and finally the node is assigned:
    mesh.conecCell{ 1, 1 } = [ 0 1 1  1 ] ;
    %md Next the frame elements MEBI parameters are set. The frame material is the first material of `materials` struct, then $1$ is assigned. The second entry of the `elements` struct correspond to the frame element employed, so $2$ is set. Finally no BC and no IC is required for this element, then $0$ is used.  Consecutive nodes build the element so then the `mesh.conecCell` is:
    for i=1:numElements,
      mesh.conecCell{ i+1,1 } = [ 1 2 0  i i+1 ] ;
    end
else
    for level = 1:Nlevel
    % level
    % Build level 
    aux1 = [  zeros(numElemCentral+1,1)  L1*(level-1)+(0:(numElemCentral))'*L1/numElemCentral zeros(numElemCentral+1,1) ] ; % trunk
    aux2 = [  zeros(numElemBranch +1,1)  (0:(numElemBranch) )'*L2*cos(alpha)  (0:(numElemBranch))'*L2*sin(alpha) ]; % right branch
    aux3 = [  zeros(numElemBranch +1,1)  (0:(numElemBranch) )'*L2*cos(alpha) -(0:(numElemBranch))'*L2*sin(alpha) ] ; % left branch
    %[ aux1 ; aux2(2:end, :)+aux1(end, :); aux3(2:end, :)+aux1(end, :)]
    if level == 1
        mesh.nodesCoords = [ aux1 ; aux2(2:end, :)+aux1(end, :); aux3(2:end, :)+aux1(end, :)] ;
        mesh.conecCell = { } ;
    else
    %    mesh.nodesCoords((level-1)*NElemLevel+2:level*NElemLevel+1, :) = [ aux1(2:end, :) ; aux2(2:end, :)+aux1(end, :); aux3(2:end, :)+aux1(end, :)] ; % Remove the first trunk node
        mesh.nodesCoords((level-1)*NElemLevel+2:level*NElemLevel+1, :) = [ aux1(2:end, :) ; aux2(2:end, :)+aux1(end, :); aux3(2:end, :)+aux1(end, :)] ; % Remove the first trunk node
    end
    % mesh.nodesCoords
    %md then the first welded node is defined with material (M) zero since nodes don't have material, the first element (E) type (the first entry of the `elements` struct), and (B) is the first entry of the the `boundaryConds` struct. For (I) no non-homogeneous initial condition is considered (then zero is used) and finally the node is assigned:
    mesh.conecCell{ 1, 1 } = [ 0 1 1  1 ] ;
    %
    for i=1:numElemCentral+numElemBranch
      j = NElemLevel*(level-1)+i;
      if i == 1 && level>1
          mesh.conecCell{ j+1,1 } = [ 1 2 0  NElemLevel*(level-2)+numElemCentral+1 j+1 ] ; % no accroche l'étage à celui d'en dessous
      else
          mesh.conecCell{ j+1,1 } = [ 1 2 0  j j+1 ] ;
      end
      %myCell2Mat( mesh.conecCell ) 
    end
    i = i +1;
    mesh.conecCell{ NElemLevel*(level-1)+i+1,1 } = [ 1 2 0  NElemLevel*(level-1)+numElemCentral+1 NElemLevel*(level-1)+i+1 ] ;
    %myCell2Mat( mesh.conecCell ) 
    for k=i+1:NElemLevel
      j = NElemLevel*(level-1)+k;
      mesh.conecCell{ j+1,1 } = [ 1 2 0  j j+1 ] ;
      %myCell2Mat( mesh.conecCell ) 
    end
    % myCell2Mat( mesh.conecCell ) 
    end
    % ADD LAST TRUNK
    aux1 = [  zeros(numElemCentral+1,1)  L1*(level+1-1)+(0:(numElemCentral))'*L1/numElemCentral zeros(numElemCentral+1,1) ];
    % mesh.nodesCoords
    mesh.nodesCoords(Nlevel*NElemLevel+2:Nlevel*NElemLevel+numElemCentral+1, :) = aux1(2:end, :);
    for i=1:numElemCentral
      j = NElemLevel*(level+1-1)+i;
      if i == 1 % link fork and fisrt node
        mesh.conecCell{ j+1,1 } = [ 1 2 0  NElemLevel*(level+1-2)+numElemCentral+1 j+1 ] ; % no accroche l'étage à celui d'en dessous
      else
        mesh.conecCell{ j+1,1 } = [ 1 2 0  j j+1 ] ; % no accroche l'étage à celui d'en dessous    
      end

    end
end
% Inverts Y and Z
qpa = mesh.nodesCoords(:,2);
mesh.nodesCoords(:,2) = mesh.nodesCoords(:,3);
mesh.nodesCoords(:,3) = qpa;
%
nodei = 76;
figure(2)
plot(mesh.nodesCoords(:,2), mesh.nodesCoords(:,3), 'o')
hold on 
plot(mesh.nodesCoords(nodei,2), mesh.nodesCoords(nodei,3), 'r*') % 76 and 101 and 150
axis([ -0.1 0.1 0 0.2])
%   Show
% myCell2Mat( mesh.conecCell )  
% mesh.nodesCoords
%save('N=1_Nelem=150mesh.mat', 'mesh')
