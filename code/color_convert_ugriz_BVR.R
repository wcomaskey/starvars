

#### this code transforms the SDSS ugriz magnitudes into Johnson B, V, R

# transformations from http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php
## Lupton (2005)
## These equations that Robert Lupton derived by matching DR4 photometry to Peter Stetson's published photometry for stars.

## Stars

##    B = u - 0.8116*(u - g) + 0.1313;  sigma = 0.0095
##    B = g + 0.3130*(g - r) + 0.2271;  sigma = 0.0107

##    V = g - 0.2906*(u - g) + 0.0885;  sigma = 0.0129
##    V = g - 0.5784*(g - r) - 0.0038;  sigma = 0.0054

##    R = r - 0.1837*(g - r) - 0.0971;  sigma = 0.0106
##    R = r - 0.2936*(r - i) - 0.1439;  sigma = 0.0072


color.transform = function(u,g,r,i,z){
  # transform SDSS ugriz to BVR
  
  B = u - 0.8116*(u - g) + 0.1313
  V = g - 0.5784*(g - r) - 0.0038
  R = r - 0.2936*(r - i) - 0.1439

  return(list(B=B,V=V,R=R))
}
