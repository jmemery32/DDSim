John,

The RodriguesParams class is just a type of Vec3D.

-wash

//
//  RodriguesToMtx - Given the Rodrigues coordinates this routine
//      formulates a 3x3 rotation matrix.
//
//      void RodriguesToMtx(double x,double y,double z,double mtx[3][3])
//           x   - (in) first Rodrigues coordinate
//           y   - (in) second Rodrigues coordinate
//           z   - (in) third Rodrigues coordinate
//           mtx - (out) rotation matrix
//

static void RodriguesToMtx(
     const RodriguesParams& r,
     double mtx[3][3])
{
     double mag = r.Magnitude() ;
     double ang = 2.0 * atan(mag) ;
     Vec3D  n = r.Normalize() ;

     if (ang == 0.0) {
         for (int i=0 ; i<3 ; ++i) {
             for (int j=0 ; j<3 ; ++j) mtx[i][j] = 0.0 ;
             mtx[i][i] = 1.0 ;
         }
     } else {
         mtx[0][0] = (1-n[0]*n[0])*cos(ang) + n[0]*n[0] ;
         mtx[0][1] = n[0]*n[1]*(1-cos(ang)) + n[2]*sin(ang) ;
         mtx[0][2] = n[0]*n[2]*(1-cos(ang)) - n[1]*sin(ang) ;

         mtx[1][0] = n[1]*n[0]*(1-cos(ang)) - n[2]*sin(ang) ;
         mtx[1][1] = (1-n[1]*n[1])*cos(ang) + n[1]*n[1] ;
         mtx[1][2] = n[1]*n[2]*(1-cos(ang)) + n[0]*sin(ang) ;

         mtx[2][0] = n[2]*n[0]*(1-cos(ang)) + n[1]*sin(ang) ;
         mtx[2][1] = n[2]*n[1]*(1-cos(ang)) - n[0]*sin(ang) ;
         mtx[2][2] = (1-n[2]*n[2])*cos(ang) + n[2]*n[2] ;
     }
}

