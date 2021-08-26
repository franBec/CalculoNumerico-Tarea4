%M Frobenius
function [x, err] = prog2_mF_FB(A,b,x0,niter)
  dimension = size(A);        %Necesario para saber nxn
  I = eye(dimension(1));      %identidad
  N = vecnorm(transpose(A));  %N(i) = norma de la fila i de A
  
  %creacion de la matriz m
  for i = 1 : dimension(1)
    m(i,i) = A(i,i)/(N(i))^2;
  end

  %calculo de radio espectral
  radioEspectral = max(abs(eig((I-m*A))));
  if (radioEspectral < 1)
    display("mFrobenius: El Radio Espectral es menor a 1");
  else
    display("mFrobenius: El Radio Espectral es mayor o igual a 1");
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
    display("mFrobenius: A es diagonalmente dominante");
  else
    display("mFrobenius: A no es diagonalmente dominante");
  end

  %metodo iterativo principal
  for i = 1 : niter
    x0 = (I-m*A)*x0 + m*b;
  end
  
  %x
  x = x0;
  
  %error
  err = norm(A*x - b);
end