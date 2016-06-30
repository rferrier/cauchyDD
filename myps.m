function ps = myps( u1, u2, K, belem, M, nodes)
 % This function computes the inner product
 index = 5;
 
 if index == 1 % The most simple IP
     ps = u1'*u2;
 elseif index == 2 % K IP
     ps = u1'*K*u2;
 elseif index == 3 % u1'*u2 on the boundary 3
     u11 = keepField( u1, 3, belem );
     u22 = keepField( u2, 3, belem );
     ps = u11'*u22;
 elseif index == 4 % int2d(u1*u2)
     ps = u1'*M*u2;
 elseif index == 5 % int1d(3)(u1*u2)
     ps = u1'*norm_bound(u2, nodes, belem, 3);
 elseif index == 6; % K IP on the bound 3
     ps = u1'*regul(u2, nodes, belem, 3);
 end
end

