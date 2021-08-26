%Gauss-Seidel
function [x, err] = prog2_Gauss_FB(A,b,x0,niter)
  dimension = size(A);  %Necesario para saber nxn
  D = diag(diag(A));    %Diagonal
  U = triu(A);          %Triangular Superior
  L = tril(A);          %Triangular Inferior
  
  %limpieza de la diagonal de U y L
  for i = 1 : dimension(1)
    U(i,i) = 0;
    L(i,i) = 0;
  end
  
  %calculo radio espectral
  matrAux = -(inv(D+L))*U;
  radioEspectral = max(abs(eig(matrAux)));
  if (radioEspectral < 1)
    display("Gauss: El radio espectral es menor a 1");
  else
    display("Gauss: El radio espectral es mayor o igual a 1");
  end
  
  %verificar si es diagonalmente dominante
  diagDom = true;
  
  for i = 1 : dimension(1)
   sumatoria(i) = 0;
   for j = 1 : dimension(2)
      if(i ~= j)
         sumatoria(i) = sumatoria(i) + abs(A(i,j));
      end
   end
   if(A(i,i)<=sumatoria(i))
    diagDom = false;
   end
  end
  
  if(diagDom)
    display("Gauss: A es diagonalmente dominante");
  else
    display("Gauss: A no es diagonalmente dominante");
  end
  
  %metodo iterativo principal
  for i = 1 : niter
    x0 = -(inv(D+L))*U*x0 + inv(D+L)*b;
  end
  
  %x
  x=x0;
  
  %error
  err = norm(A*x - b);
end