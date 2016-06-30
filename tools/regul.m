function Nu = regul(u, nodes, boundary, entity)
 % This function computes a regularization term on a boundary
 % N is such as u'Nu messes of the gradient of u
 
 choice = 1;
 if choice == 1
     Nu = Mgrad( size(u,1), nodes, boundary, entity )*u;
 elseif choice == 2 % norm H 1/2
     Nu = H1demi(size(u,1), nodes, boundary, entity)*u;
 end

end