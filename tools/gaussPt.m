function [ Xg, Wg ] = gaussPt( Ng )
% This function returns the gauss points coordinates, and weight that match
% the number Ng asked
% Values from the code of Pierre-Eric Allier, taken from :
% http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF

if Ng == 1
    Xg = [0.33333333333333 0.33333333333333];
    Wg = 1.;
elseif Ng == 2
    Xg = [0.16666666666667 0.16666666666667
          0.16666666666667 0.66666666666667
          0.66666666666667 0.16666666666667];
    Wg = [ 0.33333333333333, 0.33333333333333, 0.33333333333333 ];
elseif Ng == 3
    Xg = [0.33333333333333 0.33333333333333
          0.20000000000000 0.20000000000000
          0.20000000000000 0.60000000000000
          0.60000000000000 0.20000000000000];
    Wg = [ -0.56250000000000 ; 0.52083333333333 ;...
          0.52083333333333 ; 0.52083333333333 ];
elseif Ng == 4
    Xg = [0.44594849091597 0.44594849091597
          0.44594849091597 0.10810301816807
          0.10810301816807 0.44594849091597
          0.09157621350977 0.09157621350977
          0.09157621350977 0.81684757298046
          0.81684757298046 0.09157621350977];
    Wg = [ 0.22338158967801 ; 0.22338158967801 ; 0.22338158967801 ;...
           0.10995174365532 ; 0.10995174365532 ; 0.10995174365532 ];
elseif Ng == 5
    Xg = [0.33333333333333 0.33333333333333
          0.47014206410511 0.47014206410511
          0.47014206410511 0.05971587178977
          0.05971587178977 0.47014206410511
          0.10128650732346 0.10128650732346
          0.10128650732346 0.79742698535309
          0.79742698535309 0.10128650732346];
    Wg = [ 0.22500000000000 ; 0.13239415278851 ; 0.13239415278851 ;...
           0.13239415278851 ; 0.12593918054483 ; 0.12593918054483 ;...
           0.12593918054483 ];
elseif Ng == 6
    Xg = [0.24928674517091 0.24928674517091
          0.24928674517091 0.50142650965818
          0.50142650965818 0.2492867451709
          0.06308901449150 0.06308901449150
          0.06308901449150 0.87382197101700
          0.87382197101700 0.06308901449150
          0.31035245103378 0.63650249912140
          0.63650249912140 0.05314504984482
          0.05314504984482 0.31035245103378
          0.63650249912140 0.31035245103378
          0.31035245103378 0.05314504984482
          0.05314504984482 0.63650249912140];
    Wg = [ 0.11678627572638 ; 0.11678627572638 ; 0.11678627572638 ;...
           0.05084490637021 ; 0.05084490637021 ; 0.05084490637021 ;...
           0.08285107561837 ; 0.08285107561837 ; 0.08285107561837 ;...
           0.08285107561837 ; 0.08285107561837 ; 0.08285107561837];
elseif Ng == 7
    Xg = [0.33333333333333 0.33333333333333
          0.26034596607904 0.26034596607904
          0.26034596607904 0.47930806784192
          0.47930806784192 0.26034596607904
          0.06513010290222 0.06513010290222
          0.06513010290222 0.86973979419557
          0.86973979419557 0.06513010290222
          0.31286549600487 0.63844418856981
          0.63844418856981 0.04869031542532
          0.04869031542532 0.31286549600487
          0.63844418856981 0.31286549600487
          0.31286549600487 0.04869031542532
          0.04869031542532 0.63844418856981];
    Wg = [ -0.14957004446768 ; 0.17561525743321 ; 0.17561525743321 ;...
           0.17561525743321 ; 0.05334723560884 ; 0.05334723560884 ;...
           0.05334723560884 ; 0.07711376089026 ; 0.07711376089026 ;...
           0.07711376089026 ; 0.07711376089026 ; 0.07711376089026 ;...
           0.07711376089026];
else %Ng == 8
    Xg = [0.33333333333333 0.33333333333333
          0.45929258829272 0.45929258829272
          0.45929258829272 0.08141482341455
          0.08141482341455 0.45929258829272
          0.17056930775176 0.17056930775176
          0.17056930775176 0.65886138449648
          0.65886138449648 0.17056930775176
          0.05054722831703 0.05054722831703
          0.05054722831703 0.89890554336594
          0.89890554336594 0.05054722831703
          0.26311282963464 0.72849239295540
          0.72849239295540 0.00839477740996
          0.00839477740996 0.26311282963464
          0.72849239295540 0.26311282963464
          0.26311282963464 0.00839477740996
          0.00839477740996 0.72849239295540];
    Wg = [ 0.14431560767779 ;...
           0.09509163426728 ; 0.09509163426728 ; 0.09509163426728 ;...
           0.10321737053472 ; 0.10321737053472 ; 0.10321737053472 ;...
           0.03245849762320 ; 0.03245849762320 ; 0.03245849762320 ;...
           0.02723031417443 ; 0.02723031417443 ; 0.02723031417443 ;...
           0.02723031417443 ; 0.02723031417443 ; 0.02723031417443];
end

if Ng > 8
    warning('The required Gauss number points is too high : 8 is the max')
end


end

