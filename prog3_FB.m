function [respuesta, errores] = prog3_FB(A,b,x0,niter)
  
  %Jacobi
  [x, err] = prog1_FB(A,b,x0,niter);
  for i=1 : size(x);
    respuesta(1,i) = x(i);
  end
  errores(1,1) = err;
  
  %m_Frobenius
  [x, err] = prog2_mF_FB(A,b,x0,niter);
  for i=1 : size(x);
    respuesta(2,i) = x(i);
  end
  errores(2,1) = err;
  
  %m_Infinito
  [x, err] = prog2_mInf_FB(A,b,x0,niter);
  for i=1 : size(x);
    respuesta(3,i) = x(i);
  end
  errores(3,1) = err;
  
end