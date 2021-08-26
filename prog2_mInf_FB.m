%M infinito
function [x, err] = prog2_mInf_FB(A,b,x0,niter)
  dimension = size(A);   %Necesario para saber nxn
  I = eye(dimension(1)); %identidad
  
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
    display("mInf: A es diagonalmente dominante");
  else
    display("mInf: A no es diagonalmente dominante");
  end

  %calculos de los gi
  for i = 1 : dimension(1)
    gi(i) = A(i,i) - sumatoria(i);
  end
  
  %sg = gi menor
  sg = min(gi);
  
  %A Norma Infinita
  A_Norma_Inf = norm(A,'inf');
  
  %alfa 0
  alfa0 = 2/(A_Norma_Inf + sg);
  
  %calculo de radio espectral
  radioEspectral = max(abs(eig((I-alfa0*A))));
  if (radioEspectral < 1)
    display("mInf: El Radio Espectral es menor a 1");
  else
    display("mInf: El Radio Espectral es mayor o igual a 1");
  end

  %proceso iterativo principal
  for i = 1 : niter
    x0 = (I-alfa0*A)*x0 + alfa0*b;
  end
  
  %x
  x=x0;
  
  %error
  err = norm(A*x - b);
end