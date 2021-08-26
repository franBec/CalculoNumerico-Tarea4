%Jacobi
function [x, err] = prog1_FB(A,b,x0,niter)
  dimension = size(A);      %Necesario para saber nxn
  D1 = inv(diag(diag(A)));  %D^-1
  I = eye(dimension(1));    %identidad
  
  %calculo de radio espectral
  radioEspectral = max(abs(eig((I-D1*A))));
  if (radioEspectral < 1)
    display("Jacob: El Radio Espectral es menor a 1");
  else
    display("Jacob: El Radio Espectral es mayor o igual a 1");
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
    break;
   end
  end

  if(diagDom)
    display("Jacob: A es diagonalmente dominante");
  else
    display("Jacob: A no es diagonalmente dominante");
  end
  
  %metodo iterativo principal
  for i = 1 : niter
    x0 = (I-D1*A)*x0 + D1*b;
  end
  
  %x
  x=x0;
  
  %error
  err = norm(A*x - b);
end