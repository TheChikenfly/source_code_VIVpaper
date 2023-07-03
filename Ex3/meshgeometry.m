% Geometry
%-----------------------Cantilever Beam---------------------------------------
%mdThe coordinates of the mesh nodes are given by the matrix:
mesh.nodesCoords = [  zeros(numElements+1,1) (0:(numElements))'*l/numElements zeros(numElements+1,1) ] ;
%mdThe connectivity is introduced using the _conecCell_. Each entry of the cell contains a vector with the four indexes of the MEBI parameters, followed by the indexes of nodes that compose the element (node connectivity). For didactical purposes each element entry is commented. First the cell is initialized:
mesh.conecCell = { } ;
%md then the first welded node is defined with material (M) zero since nodes don't have material, the first element (E) type (the first entry of the `elements` struct), and (B) is the first entry of the the `boundaryConds` struct. For (I) no non-homogeneous initial condition is considered (then zero is used) and finally the node is assigned:
mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
%md Next the frame elements MEBI parameters are set. The frame material is the first material of `materials` struct, then $1$ is assigned. The second entry of the `elements` struct correspond to the frame element employed, so $2$ is set. Finally no BC and no IC is required for this element, then $0$ is used.  Consecutive nodes build the element so then the `mesh.conecCell` is:
for i=1:numElements,
  mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
end
% %-----------------------Fork---------------------------------------
% numElements = 4;
% numElemCentral = 6 ;
% numElemBranch = 6;
% L1 = l;
% L2 = l/sqrt(2);
% 
% aux1 = [  zeros(numElemCentral+1,1)  (0:(numElemCentral))'*L1/numElemCentral zeros(numElemCentral+1,1) ];
% aux2 = [  zeros(numElemBranch +1,1)  (0:(numElemBranch) )'*L2/numElemBranch  (0:(numElemBranch))'*L2/numElemBranch ];
% aux3 = [  zeros(numElemBranch +1,1)  (0:(numElemBranch) )'*L2/numElemBranch  -(0:(numElemBranch))'*L2/numElemBranch ] ;
% [ aux1 ; aux2(2:end, :)+aux1(end, :); aux3(2:end, :)+aux1(end, :)]
% mesh.nodesCoords = [ aux1 ; aux2(2:end, :)+aux1(end, :); aux3(2:end, :)+aux1(end, :)] ;
% mesh.conecCell = { } ;
% %md then the first welded node is defined with material (M) zero since nodes don't have material, the first element (E) type (the first entry of the `elements` struct), and (B) is the first entry of the the `boundaryConds` struct. For (I) no non-homogeneous initial condition is considered (then zero is used) and finally the node is assigned:
% mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1 ] ;
% %
% for i=1:numElemCentral+numElemBranch,
%   mesh.conecCell{ i+1,1 } = [ 1 2 0 0  i i+1 ] ;
% end
% i = i +1;
% mesh.conecCell{ i+1,1 } = [ 1 2 0 0  numElemCentral+1 i+1 ] ;
% for j=i+1:numElemCentral+2*numElemBranch,
%   mesh.conecCell{ j+1,1 } = [ 1 2 0 0  j j+1 ] ;
% end
% myCell2Mat( mesh.conecCell )  ;
% %---------------------Coral--------------------------------------------------------
% numElements = 8;
% %mesh.nodesCoords =[ 0 0 0; 0.5 0 0; 1 0 0; 2 0 0; 3 0 0; 2 0 1; 2 0 -1; 3 0 1; 3 0 -1]; % coral along x
% mesh.nodesCoords =[ 0 0 0; 0 0.5*l 0; 0 l 0; 0 2*l 0; 0 3*l 0; 0 2*l l; 0 2*l -l; 0 3*l l; 0 3*l -l]; % coral along x
% %mesh.nodesCoords =[ 0 0 0; 0 0 0.5; 0 0 1; 0 0 2; 0 0 3; 0 1 2; 0 -1 2; 0 1 3 ; 0 -1 3]; % coral along z
% mesh.conecCell = { } ;
% mesh.conecCell{ 1, 1 } = [ 0 1 1 0  1   ] ;
% mesh.conecCell{ 2,1 } = [ 1 2 0 0  1 2 ] ;
% mesh.conecCell{ 3,1 } = [ 1 2 0 0  2 3 ] ;
% mesh.conecCell{ 4,1 } = [ 1 2 0 0  3 4 ] ;
% mesh.conecCell{ 5,1 } = [ 1 2 0 0  4 5 ] ;
% mesh.conecCell{ 6,1 } = [ 1 2 0 0  3 6 ] ;
% mesh.conecCell{ 7,1 } = [ 1 2 0 0  3 7 ] ;
% mesh.conecCell{ 8,1 } = [ 1 2 0 0  4 8 ] ;
% mesh.conecCell{ 9,1 } = [ 1 2 0 0  4 9 ] ;
% %----------------------------------------------------------------------------------