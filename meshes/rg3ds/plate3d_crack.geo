
//Include "plate3d_crack0.geo";

// Yay, it works !
Plugin(Crack).Dimension = 2 ; 
Plugin(Crack).PhysicalGroup = (7) ; 
//Plugin(Crack).OpenBoundaryPhysicalGroup = (107) ; 
Plugin(Crack).Run ;
