function [] = BuildFVMElementLoop(EWPMesh,Bounded)
% Define size options

E = length(EWPMesh.Elements);
N = length(EWPMesh.NodePos);



Size1x1 = coder.typeof(ones(1));
Size9x6 = coder.typeof(ones(9,6));

if Bounded == 1
    SizeEx1 = coder.typeof(ones(E,1));
    SizeNx1 = coder.typeof(ones(N,1));
    SizeNx9 = coder.typeof(ones(N,9));
    SizeEx6 = coder.typeof(ones(E,6));
    SizeEx9 = coder.typeof(ones(E,9));
    Size6xE = coder.typeof(ones(6,E));

else


    SizeEx1 = coder.typeof(ones(E,1),[inf,1],1);
    SizeNx1 = coder.typeof(ones(N,1),[inf,1],1);
    SizeNx9 = coder.typeof(ones(N,9),[inf,9],1);
    SizeEx6 = coder.typeof(ones(E,6),[inf,6],1);
    SizeEx9 = coder.typeof(ones(E,9),[inf,9],1);
    Size6xE = coder.typeof(ones(6,E),[6,inf],1);

end


codegen FVMTPETransPoreElementLoop.m -args {Size1x1, Size1x1, SizeEx6, SizeEx9, SizeEx9, SizeEx9, SizeEx1, Size9x6, Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE Size6xE, SizeNx1, Size1x1, Size1x1, SizeEx9, SizeEx9, SizeEx9, SizeEx9, SizeEx9, SizeEx9, SizeEx9, SizeEx9, SizeEx9, SizeEx9, SizeEx9, SizeEx9, SizeEx9}

end

