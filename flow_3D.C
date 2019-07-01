#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <math.h>
enum direction{
  ZERO_ZERO_ZERO,
  ZERO_ZERO_PLUS,
  ZERO_ZERO_MINUS,
  PLUS_ZERO_ZERO,
  PLUS_ZERO_PLUS,
  PLUS_ZERO_MINUS,
  ZERO_PLUS_ZERO,
  ZERO_PLUS_PLUS,
  ZERO_PLUS_MINUS,
  MINUS_ZERO_ZERO,
  MINUS_ZERO_PLUS,
  MINUS_ZERO_MINUS,
  ZERO_MINUS_ZERO,
  ZERO_MINUS_PLUS,
  ZERO_MINUS_MINUS,
  PLUS_PLUS_ZERO,
  PLUS_PLUS_PLUS,
  PLUS_PLUS_MINUS,
  MINUS_PLUS_ZERO,
  MINUS_PLUS_PLUS,
  MINUS_PLUS_MINUS,
  MINUS_MINUS_ZERO,
  MINUS_MINUS_PLUS,
  MINUS_MINUS_MINUS,
  PLUS_MINUS_ZERO,
  PLUS_MINUS_PLUS,
  PLUS_MINUS_MINUS,
};

#define N_DV 27
#define SIZE_X 128
#define SIZE_Y 128
#define SIZE_Z 2
const int r=10;

//define structure for the D2Q9 model
template <typename dataType>
struct ModelD2Q9
{
  dataType ci_x[N_DV];
  dataType ci_y[N_DV];
  dataType ci_z[N_DV];

  dataType dVopposite[N_DV];
  dataType cc[N_DV];
  dataType wt[N_DV];
  dataType theta0,oneBytheta0,sqrtTheta0;
  dataType c2[N_DV];
  dataType rho;
  dataType u1,u2,u3,theta;
  dataType inlet;
  dataType Ma;
  dataType fEQ[N_DV];
  dataType denominator;
};

//define grid class
class Grid
{
      public:
//       double data[ (SIZE_X+2)*(SIZE_Y+2)*N_DV ];

	double *data;
	Grid();
      
      double  operator()(const int i,const int j, const int k,const int dv) const  { return data[  k*(SIZE_X+2)*(SIZE_Y+2)*N_DV  +  j*(SIZE_X+2)*N_DV +  i*N_DV + dv ];}
      double& operator()(const int i,const int j,const int k, const int dv)        { return data[  k*(SIZE_X+2)*(SIZE_Y+2)*N_DV  +  j*(SIZE_X+2)*N_DV +  i*N_DV + dv];}
	
      void initializeGrid()
      {      
	      for (int i = 0; i < (SIZE_X + 2) * (SIZE_Y + 2) * (SIZE_Z + 2)*(N_DV); i++)
		      data[i] = i;
    
      }
  
};
//constructor
Grid::Grid()
    {
        data = new double [(SIZE_X+2)*(SIZE_Y+2)*(SIZE_Z+2)*N_DV ];
    }

void advection(Grid &gridLB) {
    for (int k = 1; k <= SIZE_Z; k++) {
        for (int j = 1; j <= SIZE_Y; j++) {
            for (int i = 1; i <= SIZE_X; i++) {//ZERO_PLUS
                gridLB(i, j, k, MINUS_ZERO_ZERO) = gridLB(i + 1, j, k, MINUS_ZERO_ZERO);
                gridLB(i, j, k, MINUS_MINUS_ZERO) = gridLB(i + 1, j + 1, k, MINUS_MINUS_ZERO);
                gridLB(i, j, k, PLUS_MINUS_ZERO) = gridLB(i - 1, j + 1, k, PLUS_MINUS_ZERO);
                gridLB(i, j, k, ZERO_MINUS_ZERO) = gridLB(i, j + 1, k, ZERO_MINUS_ZERO);

                gridLB(i, j, k, PLUS_MINUS_MINUS) = gridLB(i - 1, j + 1, k + 1, PLUS_MINUS_MINUS);
                gridLB(i, j, k, MINUS_MINUS_MINUS) = gridLB(i + 1, j + 1, k + 1, MINUS_MINUS_MINUS);
                gridLB(i, j, k, ZERO_MINUS_MINUS) = gridLB(i, j + 1, k + 1, ZERO_MINUS_MINUS);
                gridLB(i, j, k, PLUS_PLUS_MINUS) = gridLB(i - 1, j - 1, k + 1, PLUS_PLUS_MINUS);
                gridLB(i, j, k, MINUS_PLUS_MINUS) = gridLB(i + 1, j - 1, k + 1, MINUS_PLUS_MINUS);
                gridLB(i, j, k, ZERO_PLUS_MINUS) = gridLB(i, j - 1, k + 1, ZERO_PLUS_MINUS);
                gridLB(i, j, k, PLUS_ZERO_MINUS) = gridLB(i - 1, j, k + 1, PLUS_ZERO_MINUS);
                gridLB(i, j, k, MINUS_ZERO_MINUS) = gridLB(i + 1, j, k + 1, MINUS_ZERO_MINUS);
                gridLB(i, j, k, ZERO_ZERO_MINUS) = gridLB(i , j, k + 1, ZERO_ZERO_MINUS);
            }
        }
    }
   for (int k = SIZE_Z; k >= 1; k--) {
       for (int j = SIZE_Y; j >= 1; j--) {
            for (int i = SIZE_X; i >= 1; i--) {

                gridLB(i, j, k, PLUS_ZERO_ZERO)  = gridLB(i - 1, j, k, PLUS_ZERO_ZERO);
                gridLB(i, j, k, PLUS_PLUS_ZERO)  = gridLB(i - 1, j - 1, k, PLUS_PLUS_ZERO);
                gridLB(i, j, k, ZERO_PLUS_ZERO)  = gridLB(i, j - 1, k, ZERO_PLUS_ZERO);
                gridLB(i, j, k, MINUS_PLUS_ZERO) = gridLB(i + 1, j - 1, k, MINUS_PLUS_ZERO);
//
                gridLB(i, j, k, MINUS_PLUS_PLUS) = gridLB(i + 1, j - 1, k - 1, MINUS_PLUS_PLUS);
                gridLB(i, j, k, PLUS_PLUS_PLUS)  = gridLB(i - 1, j - 1, k - 1, PLUS_PLUS_PLUS);
                gridLB(i, j, k, ZERO_PLUS_PLUS)  = gridLB(i, j - 1, k - 1, ZERO_PLUS_PLUS);
                gridLB(i, j, k, MINUS_MINUS_PLUS)= gridLB(i + 1, j + 1, k - 1, MINUS_MINUS_PLUS);
                gridLB(i, j, k, PLUS_MINUS_PLUS) = gridLB(i - 1, j + 1, k - 1, PLUS_MINUS_PLUS);
                gridLB(i, j, k, ZERO_MINUS_PLUS) = gridLB(i, j + 1, k - 1, ZERO_MINUS_PLUS);
                gridLB(i, j, k, MINUS_ZERO_PLUS) = gridLB(i + 1, j, k - 1, MINUS_ZERO_PLUS);
                gridLB(i, j, k, PLUS_ZERO_PLUS)  = gridLB(i - 1, j, k - 1, PLUS_ZERO_PLUS);
                gridLB(i, j, k, ZERO_ZERO_PLUS)  = gridLB(i , j, k - 1, ZERO_ZERO_PLUS);
            }
        }
    }
}

void getfEQ(ModelD2Q9<double> &myModel)
{
  double uDotC, uSq, tempX, tempY, term1, term2, term3, term4, term5, oneByTerm3, oneByTerm4 ;	
	
  for(int i=0;i<N_DV;i++)
  {
	 uDotC = myModel.u1*myModel.ci_x[i] + myModel.u2*myModel.ci_y[i]+myModel.u3*myModel.ci_z[i];
	 uSq   = myModel.u2*myModel.u2+myModel.u1*myModel.u1 + myModel.u3*myModel.u3;
	 
//	 myModel.fEQ[i]= myModel.wt[i]*myModel.rho*(1.0 + uDotC*myModel.oneBytheta0 - 1.5*uSq + 4.5*uDotC*uDotC );

// 	 myModel.fEQ[i]=myModel.wt[i]*myModel.rho*(1.0 + uDotC*myModel.oneBytheta0 -1.5*(myModel.u2*myModel.u2+myModel.u1*myModel.u1) +4.5*uDotC*uDotC  );
 	 myModel.fEQ[i]= myModel.wt[i]*myModel.rho*(1.0 + uDotC*myModel.oneBytheta0 - 1.5*uSq + 4.5*uDotC*uDotC  + 4.5*uDotC*uDotC*uDotC -4.5*uSq*uDotC );
  }
  
//   tempX = sqrt(1+3.0*myModel.u1*myModel.u1) ;
//   tempY = sqrt(1+3.0*myModel.u2*myModel.u2) ;
//   
//   term1 = 2.0 - tempX;
//   term2 = 2.0 - tempY;  
//   term5 = term1*term2;
//   
//   term3 = (2.0*myModel.u1 + tempX)/(1.0-myModel.u1);
//   term4 = (2.0*myModel.u2 + tempY)/(1.0-myModel.u2);
//   oneByTerm3 = 1.0/term3;
//   oneByTerm4 = 1.0/term4;
//                                                                            ;
//   myModel.fEQ[ZERO_ZERO  ] = myModel.rho*myModel.wt[ZERO_ZERO  ]*term5             ;
//   myModel.fEQ[PLUS_ZERO  ] = myModel.rho*myModel.wt[PLUS_ZERO  ]*term5*term3       ;
//   myModel.fEQ[ZERO_PLUS  ] = myModel.rho*myModel.wt[ZERO_PLUS  ]*term5*term4       ;
//   myModel.fEQ[MINUS_ZERO ] = myModel.rho*myModel.wt[MINUS_ZERO ]*term5*oneByTerm3       ;
//   myModel.fEQ[ZERO_MINUS ] = myModel.rho*myModel.wt[ZERO_MINUS ]*term5*oneByTerm4       ;
//   myModel.fEQ[PLUS_PLUS  ] = myModel.rho*myModel.wt[PLUS_PLUS  ]*term5*term3*term4 ;
//   myModel.fEQ[MINUS_PLUS ] = myModel.rho*myModel.wt[MINUS_PLUS ]*term5*oneByTerm3*term4 ;
//   myModel.fEQ[MINUS_MINUS] = myModel.rho*myModel.wt[MINUS_MINUS]*term5*oneByTerm3*oneByTerm4 ;
//   myModel.fEQ[PLUS_MINUS ] = myModel.rho*myModel.wt[PLUS_MINUS ]*term5*term3*oneByTerm4 ;
  
  // Exact in temperature

//   double dot, K, ux2(myModel.u1*myModel.u1), uy2(myModel.u2*myModel.u2), uSq(ux2 + uy2), theta0(1.0/3.0), oneByTheta0(3.0);
//   double eta(myModel.theta/theta0-1.0);
//   double theta2(myModel.theta*myModel.theta);
//   double oneByTheta(1.0/myModel.theta);
//   double oneByTheta2(oneByTheta*oneByTheta) ;
//  
//   myModel.fEQ[0] = myModel.rho*(1.0-myModel.theta)*(1.0-myModel.theta);
//   for (int dv = 1; dv < 5; dv++)
//     myModel.fEQ[dv] = myModel.rho*(1.0-myModel.theta)*0.5*myModel.theta;
//   for (int dv = 5; dv < 9; dv++)
//     myModel.fEQ[dv] = myModel.rho*0.25*myModel.theta*myModel.theta;
// 
//   for (int dv = 0; dv < 9; dv++)
//   {
//     dot = (myModel.u1*myModel.ci_x[dv] + myModel.u2*myModel.ci_y[dv])*oneByTheta;
//     K = 2.0*myModel.theta*myModel.theta/(1.0-myModel.theta) + 0.5*myModel.cc[dv]*(1.0-3.0*myModel.theta)/(1.0-myModel.theta);
//     myModel.fEQ[dv] = myModel.wt[dv]*(1.0 + dot + 0.5*dot*dot - 0.5*uSq*oneByTheta2*K ); 
//   }
//   
  
  
}

//void periodicAll(ModelD2Q9<double> &myModel,Grid &gridLB)
//{
//  for(int dv=0;dv<N_DV;dv++)
//  {
//	  for(int j=0;j<=SIZE_Y+1;j++)
//	   {
//	    //copy from right wall to left ghost
//	    gridLB(0, j, dv)   = gridLB(SIZE_X, j, dv);
//	    //copy from left wall to right ghost
//	    gridLB(SIZE_X+1, j, dv)    = gridLB(1, j, dv);
//	   }
//	   for(int i=0;i<=SIZE_X+1;i++)
//	   {
//	      //copy from bottom wall to top ghost
//	      gridLB(i, SIZE_Y+1, dv)=gridLB(i, 1, dv);
//	      //copy from bottom wall to top ghost
//	      gridLB(i, 0, dv)=gridLB(i, SIZE_Y, dv);
//	   }
//  }
//}

void prepareWallTopBottom(ModelD2Q9<double> &myModel,Grid &gridLB)
{
    for(int dv=0;dv<N_DV;dv++) {
        for (int k = 0; k <=SIZE_Z+1; k++) {
            for (int i = 0; i <= SIZE_X +1; i++) {  //copy from bottom wall to bottom ghost
                gridLB(i, 0, k, dv) = gridLB(i, 1, k, dv);
                //copy from top wall at to top ghost
                gridLB(i, SIZE_Y + 1, k, dv) = gridLB(i, SIZE_Y, k, dv);
            }
        }
    }
}

void prepareWallInletOutlet(ModelD2Q9<double> &myModel,Grid &gridLB)
{
    for(int dv=0;dv<N_DV;dv++)
	  {
        for (int k=0; k<= SIZE_Z+1; k++) {
            for (int j = 0; j <= SIZE_Y+1; j++) { //copy from right wall to right ghost
                gridLB(SIZE_X+1,j,k,dv) = gridLB(SIZE_X,j,k,dv);
                //copy from left wall to left ghost
                gridLB(0, j, k, dv) = gridLB(1, j, k, dv);
            }
        }
	  }
}

void prepareWallFrontBack(ModelD2Q9<double> &myModel,Grid &gridLB)
{
    for(int dv=0;dv<N_DV;dv++)
    {
        for (int i=0; i<= SIZE_X+1; i++) {
            for (int j = 0; j <= SIZE_Y+1; j++) { //copy from right wall to right ghost
                gridLB(i,j,0,dv) = gridLB(i,j,1,dv);
                //copy from left wall to left ghost
                gridLB(i, j, SIZE_Z+1, dv) = gridLB(i, j, SIZE_Z, dv);
            }
        }
    }
}

//void prepareWallObject(ModelD2Q9<double> &myModel,Grid &gridLB, double fLeft[][9], double fRight[][9], double fTop[][9], double fBottom[][9], int objectOriginX,int objectOriginY,int objectLengthX,int objectLengthY)
//{
//
//   int objectEndX = objectOriginX + objectLengthX - 1;
//   int objectEndY = objectOriginY + objectLengthY - 1;
//
//   for(int i2 = objectOriginY-1; i2 <= objectEndY+1; i2++)
//   for(int dv=0;dv<N_DV;dv++)
//   {
//      fLeft[i2-objectOriginY+1][dv] = gridLB(objectOriginX-1,i2,dv);
//     fRight[i2-objectOriginY+1][dv] = gridLB(objectEndX+1   ,i2,dv);
//   }
//
//   for(int i1 = objectOriginX-1; i1 <= objectEndX+1; i1++)
//   for(int dv=0;dv<N_DV;dv++)
//   {
//     fBottom[i1-objectOriginX+1][dv] = gridLB(i1,objectOriginY-1,dv);
//        fTop[i1-objectOriginX+1][dv] = gridLB(i1,objectEndY+1   ,dv);
//   }
//
//}
//
//void prepareWallCylinder(ModelD2Q9<double> &myModel,Grid &gridLB,double fwall[][2*r+3][N_DV], int x,int y,int r)
//{
//
//
//    for (int i=0;i<=2*r+2;i++) {
//        for (int j = 0; j <= 2 * r + 2; j++) {
//            for (int dv = 0; dv < N_DV; dv++) {
//                fwall[i][j][dv] = 0.0;
//            }
//        }
//    }
//    int xCentre, yCentre;
//    xCentre = r+1;
//    yCentre = r+1;
//
//    for (int i=0; i<=2*r+2; i++) {
//        for (int j=0; j <= 2*r + 2; j++) {
//            if( (i-xCentre)*(i-xCentre) + (j-yCentre)*(j-yCentre) <= r*r)
//	    {
//                for (int dv=0; dv < N_DV; dv++) {
//
//                if ( (i-xCentre-myModel.ci_x[dv])*(i-xCentre-myModel.ci_x[dv]) + (j-yCentre-myModel.ci_y[dv])*(j-yCentre-myModel.ci_y[dv]) > r*r) {
//                    fwall[i-(int)myModel.ci_x[dv]][j-(int)myModel.ci_y[dv]][(int)myModel.dVopposite[dv]] = gridLB( x+i-myModel.ci_x[dv], y+j-myModel.ci_y[dv], dv);
//                }
//
//                }
//
//
//            }
//
//	}}
//
//
//}

//
//
//void applyWallBC_old(ModelD2Q9<double> &myModel,Grid &gridLB)
//{
//      myModel.rho=1.0;
//      myModel.u1=myModel.inlet;
//      myModel.u2=0.0;
//
//      getfEQ(myModel);
//      myModel.denominator=1.0/(myModel.fEQ[ZERO_MINUS]+myModel.fEQ[PLUS_MINUS]+myModel.fEQ[MINUS_MINUS]);
//
//      double rhoTemp;
//      double factor,denom;
//
//      for(int j=2;j <= SIZE_Y-1;j++)
//	   {
//	     //bounceback from right wall
//	     factor =  gridLB(SIZE_X+1, j, PLUS_ZERO) + gridLB(SIZE_X+1, j, PLUS_PLUS) + gridLB(SIZE_X+1, j, PLUS_MINUS)   ;
//	    denom  =  1.0/(myModel.wt[PLUS_ZERO] + myModel.wt[PLUS_PLUS] + myModel.wt[PLUS_MINUS] );
//	    gridLB(SIZE_X, j, MINUS_ZERO ) = factor*denom*myModel.wt[MINUS_ZERO ];
//	    gridLB(SIZE_X, j, MINUS_MINUS) = factor*denom*myModel.wt[MINUS_MINUS];
//	    gridLB(SIZE_X, j, MINUS_PLUS ) = factor*denom*myModel.wt[MINUS_PLUS ];
//
//	    //bounceback from left wall
//            factor =  gridLB(0, j, MINUS_ZERO) + gridLB(0, j, MINUS_PLUS) + gridLB(0, j, MINUS_MINUS)   ;
//	    denom  =  1.0/(myModel.wt[MINUS_ZERO] + myModel.wt[MINUS_PLUS] + myModel.wt[MINUS_MINUS] );
//	    gridLB(1, j, PLUS_ZERO ) = factor*denom*myModel.wt[PLUS_ZERO ];
//	    gridLB(1, j, PLUS_MINUS) = factor*denom*myModel.wt[PLUS_MINUS];
//	    gridLB(1, j, PLUS_PLUS ) = factor*denom*myModel.wt[PLUS_PLUS ];
//	   }
//
//      for(int i=2;i<=SIZE_X-1;i++)
//	   {
//	      //classical diffuse from bottom wall
//	      factor =  gridLB(i, 0, ZERO_MINUS) + gridLB(i, 0, PLUS_MINUS) + gridLB(i, 0, MINUS_MINUS)   ;
//	      denom  =  1.0/(myModel.wt[ZERO_MINUS] + myModel.wt[PLUS_MINUS] + myModel.wt[MINUS_MINUS] );
//	      gridLB(i, 1, ZERO_PLUS ) = factor*denom*myModel.wt[ZERO_PLUS ] ;
//	      gridLB(i, 1, MINUS_PLUS) = factor*denom*myModel.wt[MINUS_PLUS] ;
//	      gridLB(i, 1, PLUS_PLUS ) = factor*denom*myModel.wt[PLUS_PLUS ] ;
//
//	      //classical diffuse top wall
//	      factor =  gridLB(i, SIZE_Y+1, ZERO_PLUS) + gridLB(i, SIZE_Y+1, PLUS_PLUS) + gridLB(i, SIZE_Y+1, MINUS_PLUS)   ;
//	      gridLB(i, SIZE_Y, ZERO_MINUS) = factor*myModel.denominator*myModel.fEQ[ZERO_MINUS ];
//	      gridLB(i, SIZE_Y, PLUS_MINUS) = factor*myModel.denominator*myModel.fEQ[PLUS_MINUS ];
//	      gridLB(i, SIZE_Y, MINUS_MINUS)= factor*myModel.denominator*myModel.fEQ[MINUS_MINUS];
//	   }
//
//      //BOTTOM LEFT CORNER
//	      gridLB(1 , 1, ZERO_PLUS )=  gridLB(1 , 0, ZERO_MINUS) ;
//	      gridLB(1 , 1, PLUS_PLUS )=  gridLB(1 , 0, MINUS_MINUS);
//	      gridLB(1 , 1, PLUS_ZERO )=  gridLB(1 , 0, MINUS_ZERO) ;
//	      gridLB(1 , 1, PLUS_MINUS)=  gridLB(1 , 0, MINUS_PLUS) ;
//	      gridLB(1 , 1, MINUS_PLUS)=  gridLB(1 , 0, PLUS_MINUS) ;
//      //BOTTOM RIGHT CORNER
//	      gridLB(SIZE_X, 1 , ZERO_PLUS  ) = gridLB(SIZE_X, 0 , ZERO_MINUS) ;
//	      gridLB(SIZE_X, 1 , MINUS_PLUS ) = gridLB(SIZE_X, 0 , PLUS_MINUS) ;
//	      gridLB(SIZE_X, 1 , MINUS_ZERO ) = gridLB(SIZE_X, 0 , PLUS_ZERO) ;
//	      gridLB(SIZE_X, 1 , PLUS_PLUS  ) = gridLB(SIZE_X, 0 , MINUS_MINUS);
//	      gridLB(SIZE_X, 1 , MINUS_MINUS )= gridLB(SIZE_X, 0 , PLUS_PLUS) ;
//      //TOP LEFT CORNER
//	      gridLB(1 , SIZE_Y, PLUS_PLUS  ) =  gridLB(1 , SIZE_Y+1, MINUS_MINUS) ;
//	      gridLB(1 , SIZE_Y, MINUS_MINUS) =  gridLB(1 , SIZE_Y+1, PLUS_PLUS);
//	      gridLB(1 , SIZE_Y, PLUS_ZERO  ) =  gridLB(1 , SIZE_Y+1, MINUS_ZERO )  ;
//	      gridLB(1 , SIZE_Y, PLUS_MINUS ) =  gridLB(1 , SIZE_Y+1, MINUS_PLUS)  ;
//	      gridLB(1 , SIZE_Y, ZERO_MINUS ) =  gridLB(1 , SIZE_Y+1, ZERO_PLUS)  ;
//
//	      rhoTemp = 0.0;
//	      for(int dv=0;dv<N_DV;dv++)
//		      rhoTemp += gridLB(1, SIZE_Y, dv);
//	      for(int dv=0;dv<N_DV;dv++)
//		      gridLB(1, SIZE_Y, dv) = rhoTemp*myModel.fEQ[dv];
//
//      //TOP RIGHT CORNER
//	      gridLB(SIZE_X , SIZE_Y, PLUS_MINUS ) =  gridLB(SIZE_X , SIZE_Y+1, MINUS_PLUS) ;
//	      gridLB(SIZE_X , SIZE_Y, MINUS_PLUS ) =  gridLB(SIZE_X , SIZE_Y+1, PLUS_MINUS);
//	      gridLB(SIZE_X , SIZE_Y, MINUS_ZERO ) =  gridLB(SIZE_X , SIZE_Y+1, PLUS_ZERO ) ;
//	      gridLB(SIZE_X , SIZE_Y, MINUS_MINUS) =  gridLB(SIZE_X , SIZE_Y+1, PLUS_PLUS) ;
//	      gridLB(SIZE_X , SIZE_Y, ZERO_MINUS ) =  gridLB(SIZE_X , SIZE_Y+1, ZERO_PLUS) ;
//
//	      rhoTemp = 0.0;
//	      for(int dv=0;dv<N_DV;dv++)
//		      rhoTemp += gridLB(SIZE_X, SIZE_Y, dv);
//	      for(int dv=0;dv<N_DV;dv++)
//		      gridLB(SIZE_X, SIZE_Y, dv) = rhoTemp*myModel.fEQ[dv];
//
//
//}

void getMoments(ModelD2Q9<double> &myModel,Grid &gridLB,int x_coord,int y_coord, int z_coord)
{
    myModel.rho = 0.0;
    myModel.u1  = 0.0;
    myModel.u2  = 0.0;
    myModel.u3  = 0.0;
    myModel.theta= 0.0;
    for(int k=0;k<N_DV;k++)
    {
        myModel.rho   += gridLB(x_coord, y_coord, z_coord,k);
        myModel.u1    += gridLB(x_coord, y_coord, z_coord,k) * myModel.ci_x[k];
        myModel.u2    += gridLB(x_coord, y_coord, z_coord,k) * myModel.ci_y[k];
        myModel.u3    += gridLB(x_coord, y_coord, z_coord,k) * myModel.ci_z[k];
        myModel.theta += gridLB(x_coord, y_coord, z_coord,k) * myModel.cc[k];
    }

    myModel.u1 /= myModel.rho;
    myModel.u2 /= myModel.rho;
    myModel.u3 /= myModel.rho;
//    if(x_coord==10 && y_coord==10&&z_coord==10)
//    {
//        std::cout<<"  "<<myModel.rho<<std::endl;
//    }
//     if(x_coord>=x && x_coord<x+10&&y_coord>=y && y_coord<y+10){
//         myModel.u1 = 0.0;
//         myModel.u2 = 0.0;
//     }
// 	myModel.theta /= myModel.rho;
// 	myModel.theta -= (myModel.u1*myModel.u1 + myModel.u2*myModel.u2);
// 	myModel.theta *= 0.5;

    myModel.theta = myModel.theta0;

}

void applyWallTopBottom(ModelD2Q9<double> &myModel,Grid &gridLB) {
    for (int k = 1; k <= SIZE_Z; k++) {
        for (int i = 1; i <= SIZE_X ; i++) {
//            //BOTTOM WALL-----------SPECULAR CONDITION
//            gridLB(i, 1,k,ZERO_PLUS_ZERO) = gridLB(i, 0, k, ZERO_MINUS_ZERO);
//            gridLB(i, 1, k,MINUS_PLUS_ZERO) = gridLB(i, 0, k, MINUS_MINUS_ZERO);
//            gridLB(i, 1, k,PLUS_PLUS_ZERO) = gridLB(i, 0, k, PLUS_MINUS_ZERO);
//
//            gridLB(i, 1,k,ZERO_PLUS_PLUS) = gridLB(i, 0, k, ZERO_MINUS_PLUS);
//            gridLB(i, 1, k,MINUS_PLUS_PLUS) = gridLB(i, 0, k, MINUS_MINUS_PLUS);
//            gridLB(i, 1, k,PLUS_PLUS_PLUS) = gridLB(i, 0, k, PLUS_MINUS_PLUS);
//
//            gridLB(i, 1,k,ZERO_PLUS_MINUS) = gridLB(i, 0, k, ZERO_MINUS_MINUS);
//            gridLB(i, 1, k,MINUS_PLUS_MINUS) = gridLB(i, 0, k, MINUS_MINUS_MINUS);
//            gridLB(i, 1, k,PLUS_PLUS_MINUS) = gridLB(i, 0, k, PLUS_MINUS_MINUS);
//BOTTOM WALL------BOUNCE BACK
            gridLB(i, 1, k, ZERO_PLUS_ZERO) =   gridLB(i, 0, k, ZERO_MINUS_ZERO);
            gridLB(i, 1, k, MINUS_PLUS_ZERO) =  gridLB(i, 0, k, PLUS_MINUS_ZERO);
            gridLB(i, 1, k, PLUS_PLUS_ZERO) =   gridLB(i, 0, k, MINUS_MINUS_ZERO);
            gridLB(i, 1, k, ZERO_PLUS_PLUS) =   gridLB(i, 0, k, ZERO_MINUS_MINUS);
            gridLB(i, 1, k, MINUS_PLUS_PLUS) =  gridLB(i, 0, k, PLUS_MINUS_MINUS);
            gridLB(i, 1, k, PLUS_PLUS_PLUS) =   gridLB(i, 0, k, MINUS_MINUS_MINUS);
            gridLB(i, 1, k, ZERO_PLUS_MINUS) =  gridLB(i, 0, k, ZERO_MINUS_PLUS);
            gridLB(i, 1, k, MINUS_PLUS_MINUS) = gridLB(i, 0, k, PLUS_MINUS_PLUS);
            gridLB(i, 1, k, PLUS_PLUS_MINUS) =  gridLB(i, 0, k, MINUS_MINUS_PLUS);

            //TOP WALL----------SPECULAR CONDITIONS
//            gridLB(i, SIZE_Y, k, ZERO_MINUS_ZERO) = gridLB(i, SIZE_Y + 1,  k, ZERO_PLUS_ZERO);
//            gridLB(i, SIZE_Y,  k, PLUS_MINUS_ZERO) = gridLB(i, SIZE_Y + 1, k,  PLUS_PLUS_ZERO);
//            gridLB(i, SIZE_Y,  k, MINUS_MINUS_ZERO) = gridLB(i, SIZE_Y + 1, k,  MINUS_PLUS_ZERO);
//
//            gridLB(i, SIZE_Y, k, ZERO_MINUS_PLUS) = gridLB(i, SIZE_Y + 1,  k, ZERO_PLUS_PLUS);
//            gridLB(i, SIZE_Y,  k, PLUS_MINUS_PLUS) = gridLB(i, SIZE_Y + 1, k,  PLUS_PLUS_PLUS);
//            gridLB(i, SIZE_Y,  k, MINUS_MINUS_PLUS) = gridLB(i, SIZE_Y + 1, k,  MINUS_PLUS_PLUS);
//
//            gridLB(i, SIZE_Y, k, ZERO_MINUS_MINUS) = gridLB(i, SIZE_Y + 1,  k, ZERO_PLUS_MINUS);
//            gridLB(i, SIZE_Y,  k, PLUS_MINUS_MINUS) = gridLB(i, SIZE_Y + 1, k,  PLUS_PLUS_MINUS);
//            gridLB(i, SIZE_Y,  k, MINUS_MINUS_MINUS) = gridLB(i, SIZE_Y + 1, k,  MINUS_PLUS_MINUS);

//TOP WALL___________LID DRIVEN CAVITY
            getMoments(myModel, gridLB, i, SIZE_Y+1, k);
            myModel.u1 = myModel.inlet;
            getfEQ(myModel);
            for (int dv = 0; dv < N_DV; dv++) {
                gridLB(i, SIZE_Y, k, dv) = myModel.fEQ[dv];
            }
        }
    }
}


void applyWallFrontBack(ModelD2Q9<double> &myModel,Grid &gridLB)
{
    for (int j=1;j<=SIZE_Y;j++) {
        for (int i = 2; i <= SIZE_X - 1; i++) {
            //BOTTOM WALL
            gridLB(i, j, 1,ZERO_PLUS_ZERO) = gridLB(i, j, 0,  ZERO_PLUS_ZERO);
            gridLB(i, j, 1,MINUS_PLUS_ZERO) = gridLB(i, j, 0, MINUS_PLUS_ZERO);
            gridLB(i, j, 1,PLUS_PLUS_ZERO) = gridLB(i, j, 0,  PLUS_PLUS_ZERO);
            gridLB(i, j, 1,ZERO_PLUS_PLUS) = gridLB(i, j, 0,  ZERO_PLUS_MINUS);
            gridLB(i, j, 1,MINUS_PLUS_PLUS) = gridLB(i, j, 0, MINUS_PLUS_MINUS);
            gridLB(i, j, 1,PLUS_PLUS_PLUS) = gridLB(i, j, 0,  PLUS_PLUS_MINUS);
            gridLB(i, j, 1,ZERO_PLUS_MINUS) = gridLB(i, j, 0, ZERO_PLUS_PLUS);
            gridLB(i, j, 1,MINUS_PLUS_MINUS) = gridLB(i, j, 0,MINUS_PLUS_PLUS);
            gridLB(i, j, 1,PLUS_PLUS_MINUS) = gridLB(i, j, 0, PLUS_PLUS_PLUS);

            //TOP WALL
            gridLB(i, j,SIZE_Z, ZERO_MINUS_ZERO)  = gridLB(i, j, SIZE_Z+1,   ZERO_MINUS_ZERO);
            gridLB(i, j,SIZE_Z, PLUS_MINUS_ZERO)  = gridLB(i, j,SIZE_Z+1,    PLUS_MINUS_ZERO);
            gridLB(i, j,SIZE_Z, MINUS_MINUS_ZERO) = gridLB(i, j, SIZE_Z+1,   MINUS_MINUS_ZERO);
            gridLB(i, j,SIZE_Z, ZERO_MINUS_PLUS)  =    gridLB(i, j, SIZE_Z+1,ZERO_MINUS_MINUS);
            gridLB(i, j,SIZE_Z, PLUS_MINUS_PLUS)  =    gridLB(i, j,SIZE_Z+1, PLUS_MINUS_MINUS);
            gridLB(i, j,SIZE_Z, MINUS_MINUS_PLUS) =   gridLB(i, j, SIZE_Z+1, MINUS_MINUS_MINUS);
            gridLB(i, j,SIZE_Z, ZERO_MINUS_MINUS) =   gridLB(i, j, SIZE_Z+1, ZERO_MINUS_PLUS);
            gridLB(i, j,SIZE_Z, PLUS_MINUS_MINUS) =   gridLB(i, j,SIZE_Z+1,  PLUS_MINUS_PLUS);
            gridLB(i, j,SIZE_Z, MINUS_MINUS_MINUS)=  gridLB(i, j, SIZE_Z+1,  MINUS_MINUS_PLUS);
        }
    }
}


void applyWallInletOutlet(ModelD2Q9<double> &myModel,Grid &gridLB) {

    double rhoTemp;
    double factor, denom;
    for (int k = 1; k <= SIZE_Z; k++) {
        for (int j = 1; j <= SIZE_Y; j++) {
            //OUTLET-------FLOW
//            getMoments(myModel, gridLB, SIZE_X + 1, j, k);
//            getfEQ(myModel);
////              gridLB(SIZE_X, j, MINUS_ZERO) = myModel.fEQ[MINUS_ZERO];
////              gridLB(SIZE_X, j, MINUS_MINUS) = myModel.fEQ[MINUS_MINUS];
////              gridLB(SIZE_X, j, MINUS_PLUS) = myModel.fEQ[MINUS_PLUS];
//            for (int dv = 0; dv < N_DV; dv++) {
//                gridLB(SIZE_X, j, k, dv) = myModel.fEQ[dv];
//            }
//OUTLET-----RIGHT WALL----BOUNCE BACK
            gridLB(SIZE_X, j, k, MINUS_ZERO_ZERO) =  gridLB(SIZE_X+1, j, k, PLUS_ZERO_ZERO);
            gridLB(SIZE_X, j, k, MINUS_MINUS_ZERO) = gridLB(SIZE_X+1, j, k, PLUS_PLUS_ZERO);
            gridLB(SIZE_X, j, k, MINUS_PLUS_ZERO) =  gridLB(SIZE_X+1, j, k, PLUS_MINUS_ZERO);
            gridLB(SIZE_X, j, k, MINUS_ZERO_PLUS) =  gridLB(SIZE_X+1, j, k, PLUS_ZERO_MINUS);
            gridLB(SIZE_X, j, k, MINUS_MINUS_PLUS) = gridLB(SIZE_X+1, j, k, PLUS_PLUS_MINUS);
            gridLB(SIZE_X, j, k, MINUS_PLUS_PLUS) =  gridLB(SIZE_X+1, j, k, PLUS_MINUS_MINUS);
            gridLB(SIZE_X, j, k, MINUS_ZERO_MINUS) = gridLB(SIZE_X+1, j, k, PLUS_ZERO_PLUS);
            gridLB(SIZE_X, j, k, MINUS_MINUS_MINUS) =gridLB(SIZE_X+1, j, k, PLUS_PLUS_PLUS);
            gridLB(SIZE_X, j, k, MINUS_PLUS_MINUS) = gridLB(SIZE_X+1, j, k, PLUS_MINUS_PLUS);
//INLET --------LEFT WjLL-----BOUNCE BACK
            gridLB(1, j, k, PLUS_ZERO_ZERO) =   gridLB(0, j, k,   MINUS_ZERO_ZERO);
            gridLB(1, j, k, PLUS_MINUS_ZERO) =  gridLB(0, j, k,  MINUS_PLUS_ZERO);
            gridLB(1, j, k, PLUS_PLUS_ZERO) =   gridLB(0, j, k,   MINUS_MINUS_ZERO);
            gridLB(1, j, k, PLUS_ZERO_PLUS) =   gridLB(0, j, k,   MINUS_ZERO_MINUS);
            gridLB(1, j, k, PLUS_MINUS_PLUS) =  gridLB(0, j, k,  MINUS_PLUS_MINUS);
            gridLB(1, j, k, PLUS_PLUS_PLUS) =   gridLB(0, j, k,   MINUS_MINUS_MINUS);
            gridLB(1, j, k, PLUS_ZERO_MINUS) =  gridLB(0, j, k,  MINUS_ZERO_PLUS);
            gridLB(1, j, k, PLUS_MINUS_MINUS) = gridLB(0, j, k, MINUS_PLUS_PLUS);
            gridLB(1, j, k, PLUS_PLUS_MINUS) =  gridLB(0, j, k,  MINUS_MINUS_PLUS);
///////////////////////////////////////////////////////////////////////////////////
//INLET----------------FLOW
//            myModel.rho = 1.0;
//            myModel.u1 = myModel.inlet;
//            myModel.u2 = 0.0;
//            myModel.u3 = 0.0;
//
//            getfEQ(myModel);
//
//            for (int dv = 0; dv < N_DV; dv++) {
//                gridLB(1, j, k, dv) = myModel.fEQ[dv];
//            }
//

//           myModel.denominator=1.0/(myModel.fEQ[ZERO_MINUS]+myModel.fEQ[PLUS_MINUS]+myModel.fEQ[MINUS_MINUS]);
//            factor = gridLB(0,j,ZERO_MINUS)+ gridLB(0,j,PLUS_MINUS)+ gridLB(0,j,MINUS_MINUS);
//            gridLB(1,j,PLUS_ZERO)=factor*myModel.denominator*myModel.fEQ[PLUS_ZERO];
//            gridLB(1,j,PLUS_MINUS)=factor*myModel.denominator*myModel.fEQ[PLUS_MINUS];
//            gridLB(1,j,PLUS_PLUS)=factor*myModel.denominator*myModel.fEQ[PLUS_PLUS];
        }
    }
}

//void applyWallObject(ModelD2Q9<double> &myModel,Grid &gridLB, double fLeft[][9], double fRight[][9], double fTop[][9], double fBottom[][9], int objectOriginX,int objectOriginY,int objectLengthX,int objectLengthY)
//{
//
//
//   int objectEndX = objectOriginX + objectLengthX - 1;
//   int objectEndY = objectOriginY + objectLengthY - 1;
//
//   for(int i2 = objectOriginY; i2 <= objectEndY; i2++)
//   {
//     gridLB(objectOriginX-1,i2,MINUS_PLUS ) =  fLeft[i2-objectOriginY+1][PLUS_MINUS ] ;
//     gridLB(objectOriginX-1,i2,MINUS_ZERO ) =  fLeft[i2-objectOriginY+1][PLUS_ZERO  ] ;
//     gridLB(objectOriginX-1,i2,MINUS_MINUS) =  fLeft[i2-objectOriginY+1][PLUS_PLUS  ] ;
//
//     gridLB(objectEndX+1   ,i2,PLUS_ZERO  ) = fRight[i2-objectOriginY+1][MINUS_ZERO ] ;
//     gridLB(objectEndX+1   ,i2,PLUS_MINUS ) = fRight[i2-objectOriginY+1][MINUS_PLUS ] ;
//     gridLB(objectEndX+1   ,i2,PLUS_PLUS  ) = fRight[i2-objectOriginY+1][MINUS_MINUS] ;
//   }
//
//   for(int i1 = objectOriginX; i1 <= objectEndX; i1++)
//   {
//     gridLB(i1,objectOriginY-1,ZERO_MINUS ) = fBottom[i1-objectOriginX+1][ZERO_PLUS  ] ;
//     gridLB(i1,objectOriginY-1,MINUS_MINUS) = fBottom[i1-objectOriginX+1][PLUS_PLUS  ] ;
//     gridLB(i1,objectOriginY-1,PLUS_MINUS ) = fBottom[i1-objectOriginX+1][MINUS_PLUS ] ;
//
//     gridLB(i1,objectEndY+1   ,PLUS_PLUS  ) =    fTop[i1-objectOriginX+1][MINUS_MINUS] ;
//     gridLB(i1,objectEndY+1   ,ZERO_PLUS  ) =    fTop[i1-objectOriginX+1][ZERO_MINUS ] ;
//     gridLB(i1,objectEndY+1   ,MINUS_PLUS ) =    fTop[i1-objectOriginX+1][PLUS_MINUS ] ;
//   }
//
//
//     gridLB(objectOriginX-1,objectOriginY-1,MINUS_MINUS) =  fLeft[objectOriginY-1-objectOriginY+1][PLUS_PLUS  ] ;
//
//     gridLB(objectEndX+1   ,objectEndY+1,PLUS_PLUS  ) = fRight[objectEndY+1-objectOriginY+1][MINUS_MINUS] ;
//
//     gridLB(objectEndX+1,objectOriginY-1,PLUS_MINUS ) = fBottom[objectEndX+1-objectOriginX+1][MINUS_PLUS ] ;
//
//     gridLB(objectOriginX-1,objectEndY+1   ,MINUS_PLUS ) =    fTop[objectOriginX-1-objectOriginX+1][PLUS_MINUS ] ;
//
//}

//void applyWallCylinder(ModelD2Q9<double> &myModel,Grid &gridLB,double fwall[][2*r+3][N_DV],int x,int y,int r)
//{
//
//
////The square enclosing the circle.
//    for (int i =0; i<=2*r+2; i++){
//        for (int j=0; j<=2*r+2; j++){
//            for (int dv=0; dv<N_DV; dv++){
//
//                if(fwall[i][j][dv] != 0)
//                    gridLB(x+i,y+j,dv) = fwall[i][j][dv];
//
//    }}}
//
//
//
//}

//collision
//f=f+2*beta*(fEQ-f)
void Collide(ModelD2Q9<double> &myModel,Grid &gridLB,double twoBeta, int step)
{
	double f_i[27];
	double x_i[27];
	double ximax(0.0);
  
        for(int i=1;i<=SIZE_X;i++) {
            for (int j = 1; j <= SIZE_Y; j++) {
                for (int k = 1; k <= SIZE_Z; k++) {
                    getMoments(myModel, gridLB, i, j, k);
                    getfEQ(myModel);
                    for (int l = 0; l < N_DV; l++)
                        gridLB(i, j, k, l) = gridLB(i, j, k, l) + twoBeta * (myModel.fEQ[l] - gridLB(i, j, k, l));
                }
            }
        }
}

void copyToGrid(ModelD2Q9<double> &myModel,Grid &gridLB,int x_coord,int y_coord ,int z_coord)
{
    for(int k=0;k<N_DV;k++)
	gridLB(x_coord, y_coord, z_coord, k)= myModel.fEQ[k];
}

void InitialConditions(ModelD2Q9<double> &myModel,Grid &gridLB)
{
    double x;
    for (int k=1;k<=SIZE_Z;k++) {
        for (int i = 1; i <= SIZE_X; i++)
            for (int j = 1; j <= SIZE_Y; j++) {
//                double x = ((double)i-0.5)/((double)SIZE_X );
//                double y = ((double)j-0.5)/((double)SIZE_Y );
//
//                if(y<0.5)
//                    myModel.u1  = 0.02886*tanh( (4.0*y - 1.0)/0.05 );
//                else
//                    myModel.u1  = 0.02886*tanh( (3.0 - 4.0*y)/0.05 );
//
//                myModel.u2  = 0.02886*0.05*sin( 2.0*M_PI*(x + 0.25 ) );
               // x = ((double) i) / ((double) SIZE_X);
                myModel.rho = 1.0;
                //myModel.u1 = 0.01*sin(2*M_PI*(double)i/(double)SIZE_X)*cos(2*M_PI*(double)j/(double)SIZE_Y);
                //myModel.u2 = 0.01*sin(2*M_PI*(double)j/(double)SIZE_Y)*cos(2*M_PI*(double)i/(double)SIZE_X);
                myModel.u3 = 0.0;
                myModel.u1 = 0.01*sin(2*M_PI*(double)i/(double)SIZE_X);
                myModel.u2 = 0.02*sin(2*M_PI*(double)j/(double)SIZE_Y);
//                myModel.u3 = 0.03*sin(2*M_PI*(double)k/(double)SIZE_Z);
                myModel.theta = myModel.theta0;
                getfEQ(myModel);
                copyToGrid(myModel, gridLB, i, j,k);
                //std::cout<<i<<"\t"<<j<<std::endl;
            }
    }
}

void massConservationCheck(Grid &gridLB, int time) {
    double mass = 0.0;
    for (int k = 1; k <= SIZE_Z; k++) {
        for (int i = 1; i <= SIZE_X; i++) {
            for (int j = 1; j <= SIZE_Y; j++) {
                double count = 0.0;
                for (int l = 0; l < N_DV; l++) {
                    mass += gridLB(i, j, k, l);
                    //count +=gridLB(i,j,k,l);
                }
//                if(count<0.99999) {
//                    std::cout << i<< "  "<<j<< "  "<< k<<std::endl;
//                }
            }
        }
    }
    //mass /= SIZE_X * SIZE_Y*SIZE_Z;
    std::cout << "at time step=" << time << " global mass=" << mass << std::endl;
}
//    for (int k =0;k<=SIZE_Z+1;k++){
//        std::cout<<gridLB(SIZE_X+1,0,k,MINUS_PLUS_ZERO)<<std::endl;
//    }
//std::cout<<"Now the opposite end : ";
//    for (int k =0;k<=SIZE_Z+1;k++){
//        std::cout<<gridLB(1,SIZE_Y,k,MINUS_PLUS_ZERO)<<std::endl;
//    }

void setModelParameters(ModelD2Q9<double> &myModel)
{   	
    myModel.ci_x[ZERO_ZERO_ZERO]=0.0;
	myModel.ci_x[PLUS_ZERO_ZERO]=1.0;
	myModel.ci_x[ZERO_PLUS_ZERO]=0.0;
	myModel.ci_x[MINUS_ZERO_ZERO]=-1.0;
	myModel.ci_x[ZERO_MINUS_ZERO]=0.0;
	myModel.ci_x[PLUS_PLUS_ZERO]=1.0;
	myModel.ci_x[MINUS_PLUS_ZERO]=-1.0;
	myModel.ci_x[MINUS_MINUS_ZERO]=-1.0;
	myModel.ci_x[PLUS_MINUS_ZERO]=1.0;


    myModel.ci_x[ZERO_ZERO_PLUS]=0.0;
    myModel.ci_x[PLUS_ZERO_PLUS]=1.0;
    myModel.ci_x[ZERO_PLUS_PLUS]=0.0;
    myModel.ci_x[MINUS_ZERO_PLUS]=-1.0;
    myModel.ci_x[ZERO_MINUS_PLUS]=0.0;
    myModel.ci_x[PLUS_PLUS_PLUS]=1.0;
    myModel.ci_x[MINUS_PLUS_PLUS]=-1.0;
    myModel.ci_x[MINUS_MINUS_PLUS]=-1.0;
    myModel.ci_x[PLUS_MINUS_PLUS]=1.0;

    myModel.ci_x  [ZERO_ZERO_MINUS]=0.0;
    myModel.ci_x  [PLUS_ZERO_MINUS]=1.0;
    myModel.ci_x  [ZERO_PLUS_MINUS]=0.0;
    myModel.ci_x [MINUS_ZERO_MINUS]=-1.0;
    myModel.ci_x [ZERO_MINUS_MINUS]=0.0;
    myModel.ci_x  [PLUS_PLUS_MINUS]=1.0;
    myModel.ci_x [MINUS_PLUS_MINUS]=-1.0;
    myModel.ci_x[MINUS_MINUS_MINUS]=-1.0;
    myModel.ci_x [PLUS_MINUS_MINUS]=1.0;
    /////////////////////////////////

	myModel.ci_y  [ZERO_ZERO_ZERO]=0.0;
	myModel.ci_y  [PLUS_ZERO_ZERO]=0.0;
	myModel.ci_y  [ZERO_PLUS_ZERO]=1.0;
	myModel.ci_y [MINUS_ZERO_ZERO]=0.0;
	myModel.ci_y [ZERO_MINUS_ZERO]=-1.0;
	myModel.ci_y  [PLUS_PLUS_ZERO]=1.0;
	myModel.ci_y [MINUS_PLUS_ZERO]=1.0;
	myModel.ci_y[MINUS_MINUS_ZERO]=-1.0;
	myModel.ci_y [PLUS_MINUS_ZERO]=-1.0;

    myModel.ci_y  [ZERO_ZERO_MINUS]=0.0;
    myModel.ci_y  [PLUS_ZERO_MINUS]=0.0;
    myModel.ci_y  [ZERO_PLUS_MINUS]=1.0;
    myModel.ci_y [MINUS_ZERO_MINUS]=0.0;
    myModel.ci_y [ZERO_MINUS_MINUS]=-1.0;
    myModel.ci_y  [PLUS_PLUS_MINUS]=1.0;
    myModel.ci_y [MINUS_PLUS_MINUS]=1.0;
    myModel.ci_y[MINUS_MINUS_MINUS]=-1.0;
    myModel.ci_y [PLUS_MINUS_MINUS]=-1.0;

    myModel.ci_y  [ZERO_ZERO_PLUS]=0.0;
    myModel.ci_y  [PLUS_ZERO_PLUS]=0.0;
    myModel.ci_y  [ZERO_PLUS_PLUS]=1.0;
    myModel.ci_y [MINUS_ZERO_PLUS]=0.0;
    myModel.ci_y [ZERO_MINUS_PLUS]=-1.0;
    myModel.ci_y  [PLUS_PLUS_PLUS]=1.0;
    myModel.ci_y [MINUS_PLUS_PLUS]=1.0;
    myModel.ci_y[MINUS_MINUS_PLUS]=-1.0;
    myModel.ci_y [PLUS_MINUS_PLUS]=-1.0;

    /////////////////////////////////////////
    myModel.ci_z  [ZERO_ZERO_ZERO]= 0.0;
    myModel.ci_z  [PLUS_ZERO_ZERO]= 0.0;
    myModel.ci_z  [ZERO_PLUS_ZERO]= 0.0;
    myModel.ci_z [MINUS_ZERO_ZERO]= 0.0;
    myModel.ci_z [ZERO_MINUS_ZERO]= 0.0;
    myModel.ci_z  [PLUS_PLUS_ZERO]= 0.0;
    myModel.ci_z [MINUS_PLUS_ZERO]= 0.0;
    myModel.ci_z[MINUS_MINUS_ZERO]= 0.0;
    myModel.ci_z [PLUS_MINUS_ZERO]= 0.0;
    myModel.ci_z  [ZERO_ZERO_MINUS]=-1.0;
    myModel.ci_z  [PLUS_ZERO_MINUS]=-1.0;
    myModel.ci_z  [ZERO_PLUS_MINUS]=-1.0;
    myModel.ci_z [MINUS_ZERO_MINUS]=-1.0;
    myModel.ci_z [ZERO_MINUS_MINUS]=-1.0;
    myModel.ci_z  [PLUS_PLUS_MINUS]=-1.0;
    myModel.ci_z [MINUS_PLUS_MINUS]=-1.0;
    myModel.ci_z[MINUS_MINUS_MINUS]=-1.0;
    myModel.ci_z [PLUS_MINUS_MINUS]=-1.0;
    myModel.ci_z  [ZERO_ZERO_PLUS]= 1.0;
    myModel.ci_z  [PLUS_ZERO_PLUS]= 1.0;
    myModel.ci_z  [ZERO_PLUS_PLUS]= 1.0;
    myModel.ci_z [MINUS_ZERO_PLUS]= 1.0;
    myModel.ci_z [ZERO_MINUS_PLUS]= 1.0;
    myModel.ci_z  [PLUS_PLUS_PLUS]= 1.0;
    myModel.ci_z [MINUS_PLUS_PLUS]= 1.0;
    myModel.ci_z[MINUS_MINUS_PLUS]= 1.0;
    myModel.ci_z [PLUS_MINUS_PLUS]= 1.0;

//////////////////////////////////////////////////////

    myModel.wt[ZERO_ZERO_ZERO]=16.0/36.0*4.0/6.0;
	myModel.wt[PLUS_ZERO_ZERO]=4.0/36.0*4.0/6.0;
        myModel.wt[ZERO_PLUS_ZERO]=4.0/36.0*4.0/6.0;
	myModel.wt[MINUS_ZERO_ZERO]=4.0/36.0*4.0/6.0;
	myModel.wt[ZERO_MINUS_ZERO]=4.0/36.0*4.0/6.0;
	myModel.wt[PLUS_PLUS_ZERO]=1.0/36.0*4.0/6.0;
	myModel.wt[MINUS_PLUS_ZERO]=1.0/36.0*4.0/6.0;
	myModel.wt[MINUS_MINUS_ZERO]=1.0/36.0*4.0/6.0;
	myModel.wt[PLUS_MINUS_ZERO]=1.0/36.0*4.0/6.0;

    myModel.wt[ZERO_ZERO_MINUS]=16.0/36.0*1.0/6.0;
    myModel.wt[PLUS_ZERO_MINUS]=4.0/36.0*1.0/6.0;
    myModel.wt[ZERO_PLUS_MINUS]=4.0/36.0*1.0/6.0;
    myModel.wt[MINUS_ZERO_MINUS]=4.0/36.0*1.0/6.0;
    myModel.wt[ZERO_MINUS_MINUS]=4.0/36.0*1.0/6.0;
    myModel.wt[PLUS_PLUS_MINUS]=1.0/36.0*1.0/6.0;
    myModel.wt[MINUS_PLUS_MINUS]=1.0/36.0*1.0/6.0;
    myModel.wt[MINUS_MINUS_MINUS]=1.0/36.0*1.0/6.0;
    myModel.wt[PLUS_MINUS_MINUS]=1.0/36.0*1.0/6.0;

    myModel.wt[ZERO_ZERO_PLUS]=16.0/36.0*1.0/6.0;
    myModel.wt[PLUS_ZERO_PLUS]=4.0/36.0*1.0/6.0;
    myModel.wt[ZERO_PLUS_PLUS]=4.0/36.0*1.0/6.0;
    myModel.wt[MINUS_ZERO_PLUS]=4.0/36.0*1.0/6.0;
    myModel.wt[ZERO_MINUS_PLUS]=4.0/36.0*1.0/6.0;
    myModel.wt[PLUS_PLUS_PLUS]=1.0/36.0*1.0/6.0;
    myModel.wt[MINUS_PLUS_PLUS]=1.0/36.0*1.0/6.0;
    myModel.wt[MINUS_MINUS_PLUS]=1.0/36.0*1.0/6.0;
    myModel.wt[PLUS_MINUS_PLUS]=1.0/36.0*1.0/6.0;
   ///////////////////////////////////////////////

	myModel.dVopposite[ZERO_ZERO_ZERO  ] = ZERO_ZERO_ZERO   ;
	myModel.dVopposite[PLUS_ZERO_ZERO  ] = MINUS_ZERO_ZERO  ;
	myModel.dVopposite[ZERO_PLUS_ZERO  ] = ZERO_MINUS_ZERO  ;
	myModel.dVopposite[MINUS_ZERO_ZERO ] = PLUS_ZERO_ZERO   ;
	myModel.dVopposite[ZERO_MINUS_ZERO ] = ZERO_PLUS_ZERO   ;
	myModel.dVopposite[PLUS_PLUS_ZERO ]  = MINUS_MINUS_ZERO ;
	myModel.dVopposite[MINUS_PLUS_ZERO ] = PLUS_MINUS_ZERO  ;
	myModel.dVopposite[MINUS_MINUS_ZERO] = PLUS_PLUS_ZERO   ;
	myModel.dVopposite[PLUS_MINUS_ZERO]  = MINUS_PLUS_ZERO  ;

    myModel.dVopposite[ZERO_ZERO_PLUS  ] = ZERO_ZERO_MINUS ;
    myModel.dVopposite[PLUS_ZERO_PLUS  ] = MINUS_ZERO_MINUS ;
    myModel.dVopposite[ZERO_PLUS_PLUS  ] = ZERO_MINUS_MINUS ;
    myModel.dVopposite[MINUS_ZERO_PLUS ] = PLUS_ZERO_MINUS ;
    myModel.dVopposite[ZERO_MINUS_PLUS ] = ZERO_PLUS_MINUS ;
    myModel.dVopposite[PLUS_PLUS_PLUS]  = MINUS_MINUS_MINUS;
    myModel.dVopposite[MINUS_PLUS_PLUS ] = PLUS_MINUS_MINUS ;
    myModel.dVopposite[MINUS_MINUS_PLUS] = PLUS_PLUS_MINUS ;
    myModel.dVopposite[PLUS_MINUS_PLUS]  = MINUS_PLUS_MINUS ;

    myModel.dVopposite[ZERO_ZERO_MINUS ] = ZERO_ZERO_PLUS   ;
    myModel.dVopposite[PLUS_ZERO_MINUS  ] = MINUS_ZERO_PLUS  ;
    myModel.dVopposite[ZERO_PLUS_MINUS ] = ZERO_MINUS_PLUS  ;
    myModel.dVopposite[MINUS_ZERO_MINUS ] = PLUS_ZERO_PLUS   ;
    myModel.dVopposite[ZERO_MINUS_MINUS ] = ZERO_PLUS_PLUS  ;
    myModel.dVopposite[PLUS_PLUS_MINUS ]  = MINUS_MINUS_PLUS ;
    myModel.dVopposite[MINUS_PLUS_MINUS ] = PLUS_MINUS_PLUS  ;
    myModel.dVopposite[MINUS_MINUS_MINUS] = PLUS_PLUS_PLUS   ;
    myModel.dVopposite[PLUS_MINUS_MINUS]  = MINUS_PLUS_PLUS  ;


	////////////////////////////////////////

	myModel.theta0=1.0/3.0;
	myModel.oneBytheta0=3.0;
	myModel.sqrtTheta0=sqrt(myModel.theta0);

	for(int dv=0;dv<N_DV;dv++)
	  myModel.cc[dv] = myModel.ci_x[dv]*myModel.ci_x[dv]+ myModel.ci_y[dv]*myModel.ci_y[dv]+myModel.ci_z[dv]*myModel.ci_z[dv];

}

void printVtk(ModelD2Q9<double> &myModel,Grid &gridLB,int step, double beta, double dt)
{
    double alpha(2.0), ximax(0.0), ximin(0.0), alphaMax(0.0),temp(0.0),tempPsi(0.0);
    double f_i[9];double x_i[9];
    std::ofstream file;
    char fileName[250];
    sprintf(fileName,"./results/x-velocity_%d.vtk",step);
    file.open(fileName);

    //   vtk file header
    file<<"# vtk DataFile Version 3.0"<<std::endl<<"Velocity"<<std::endl<<"ASCII"<<std::endl<<"DATASET STRUCTURED_GRID"<<std::endl;
    file<<"DIMENSIONS "<<(SIZE_X)<<" "<<(SIZE_Y)<<" "<<"1"<<std::endl;
    file<<"POINTS "<<(SIZE_Y)*(SIZE_X)*(SIZE_Z)<<" double"<<std::endl;
    for (int i3=1 ; i3<=SIZE_Z; i3++)
    for (int i2=1 ; i2<=SIZE_Y;  i2++)
        for (int i1=1 ; i1<=SIZE_X;  i1++)
            file<<i1<<" "<<i2<<" "<<"1"<<std::endl;

    file<<"POINT_DATA "<<(SIZE_X)*(SIZE_Y)*(SIZE_Z)<<std::endl;
    file<<"VECTORS"<<" "<<"velocity"<<" "<<"double"<<std::endl;
    for (int i3=1 ; i3<=SIZE_Z; i3++)
        for (int i2=1 ; i2<=SIZE_Y;  i2++)
            for (int i1=1 ; i1<=SIZE_X;  i1++)
            {
                if(i1==1){
                    getMoments(myModel,gridLB, i1, i2,i3);
                    file << myModel.u1 << " " << myModel.u2 << " " <<"0"<< std::endl;
                    tempPsi = 0.0;
                }
                else{
                    getMoments(myModel,gridLB, i1-1, i2,i3);
                    temp = myModel.u2;
                    getMoments(myModel,gridLB, i1, i2,i3);
                    file << myModel.u1 << " " << myModel.u2 << " " <<"0"<< std::endl;
                    tempPsi = tempPsi-0.5*(temp+myModel.u2);
                }
            }
    file.close();
}

//void printVtk(ModelD2Q9<double> &myModel,Grid &gridLB,int step, double beta, double dt) {
//    double alpha(2.0), ximax(0.0), ximin(0.0), alphaMax(0.0), temp(0.0), tempPsi(0.0);
//    double f_i[9];
//    double x_i[9];
//    std::ofstream file;
//    char fileName[250];
//    sprintf(fileName, "./results/x-velocity_%d.vtk", step);
//    file.open(fileName);
//
//    //   vtk file header
//    file << "# vtk DataFile Version 3.0" << std::endl << "Velocity" << std::endl << "ASCII" << std::endl
//         << "DATASET STRUCTURED_GRID" << std::endl;
//    file << "DIMENSIONS " << (SIZE_X) << " " << (SIZE_Y) << " " << (SIZE_Z) << std::endl;
//    file << "POINTS " << (SIZE_Y) * (SIZE_X) * (SIZE_Z) << " double" << std::endl;
//    for (int i2 = 1; i2 <= SIZE_Y; i2++)
//        for (int i1 = 1; i1 <= SIZE_X; i1++)
//            file << i1 << " " << i2 << " " << "1" << std::endl;
//
//    file << "POINT_DATA " << (SIZE_X) * (SIZE_Y) * (SIZE_Z) << std::endl;
//    file << "VECTORS" << " " << "velocity" << " " << "double" << std::endl;
//    for (int i3 = 1; i3 <= SIZE_Z; i3++) {
//        for (int i2 = 1; i2 <= SIZE_Y; i2++) {
//            for (int i1 = 1; i1 <= SIZE_X; i1++) {
//                if (i1 == 1) {
//                    getMoments(myModel, gridLB, i1, i2);
//                    file << myModel.u1 << "     " << myModel.u2 << "     " << "0.0" << std::endl;
//                    tempPsi = 0.0;
//                } else {
//                    getMoments(myModel, gridLB, i1 - 1, i2);
//                    temp = myModel.u2;
//                    getMoments(myModel, gridLB, i1, i2);
//                    file << myModel.u1 << "     " << myModel.u2 << "     " << "0.0" << std::endl;
//                    tempPsi = tempPsi - 0.5 * dt * (temp + myModel.u2);
//                }
////                if(i2 == (0.5*SIZE_Y) && i1 == (0.5*SIZE_X)) {
////                    std::ofstream outdata;
////                    outdata.open("velocity.txt", std::ios_base::app);
////                    outdata << sqrt(myModel.u1 * myModel.u1 + myModel.u2 * myModel.u2) <<std::endl;
////                    outdata.close();
////            }
//
//            file.close();
//
//            }
//        }
//    }
//}
void periodicAll(ModelD2Q9<double> &myModel,Grid &gridLB)
{
    for(int dv=0;dv<N_DV;dv++)
    {
        for (int k=0;k<=SIZE_Z+1;k++) {
            for (int j = 0; j <= SIZE_Y+1; j++) {
                //copy from right wall to left ghost
                gridLB(0, j, k,dv) = gridLB(SIZE_X, j, k,dv);
                //copy from left wall to right ghost
                gridLB(SIZE_X + 1, j,k, dv) = gridLB(1, j, k,dv);
            }
        }
        for (int k=0;k<=SIZE_Z+1;k++) {
            for (int i = 0; i <= SIZE_X+1 ; i++) {
                //copy from bottom wall to top ghost
                gridLB(i, SIZE_Y + 1, k,dv) = gridLB(i, 1,k, dv);
                //copy from bottom wall to top ghost
                gridLB(i, 0,k, dv) = gridLB(i, SIZE_Y, k,dv);
            }
        }
        for (int j = 0; j <= SIZE_Y + 1; j++) {
            for (int i = 0; i <= SIZE_X + 1; i++) {
                //copy from bottom wall to top ghost
                gridLB(i, j, 0,dv) = gridLB(i, j,SIZE_Z, dv);
                //copy from bottom wall to top ghost
                gridLB(i, j,SIZE_Z+1, dv) = gridLB(i, j, 1,dv);
            }
        }

    }
}

void printstrouhal(ModelD2Q9<double> &myModel,Grid &gridLB,int step, double beta, double dt) {
    for (int i3 = 1; i3 <= SIZE_Z; i3++) {
        for (int i2 = 1; i2 <= SIZE_Y; i2++) {
            for (int i1 = 1; i1 <= SIZE_X; i1++) {
                if (i2 == (0.5 * SIZE_Y) && i1 == (0.5 * SIZE_X)) {
                    std::ofstream outdata;
                    outdata.open("velocity.txt", std::ios_base::app);
                    outdata << sqrt(myModel.u1 * myModel.u1 + myModel.u2 * myModel.u2) << std::endl;
                    outdata.close();
                }

            }

            }
        }
    }
void check(ModelD2Q9<double> &myModel, Grid &gridLB, double time){
    std::cout<<" It did come here!"<<std::endl;
    for (int k=1;k<SIZE_Z+1;k++) {
        for (int i = 0; i < SIZE_X + 1; i++) {
            for (int j = 0; j < SIZE_Y + 1; j++) {
                double dens = 0.0;
                for (int dv = 0; dv < N_DV; dv++) {
                    dens = dens + gridLB(i, j, k, dv);
                }
                if (dens > 10.0) {
                    std::cout << "Error at " << i <<" "<<j <<"  "<<k<<std::endl;
                }
            }
        }
    }
}
void print(ModelD2Q9<double> &myModel, Grid &gridLB){
    for (int j=0;j<SIZE_Y;j++){
        int i = SIZE_X/2;
        int k = SIZE_Z/2;
        getMoments(myModel,gridLB,i,j,k);
       std::cout<<myModel.u1<<std::endl;
    }
}
int main()
{
  ModelD2Q9 <double> lbModel;
  Grid lbGrid;
  setModelParameters(lbModel);
  InitialConditions(lbModel,lbGrid);
  //massConservationCheck(lbGrid, 0);
  //simulation parameters

  lbModel.Ma         = 0.05;
  lbModel.inlet = lbModel.Ma*lbModel.sqrtTheta0;

  double Re      = 300.0;
  double Kn      = lbModel.Ma/Re;
  double refLen  = SIZE_X;
  double kinVisc = lbModel.inlet*refLen/Re;
  double tau     = kinVisc*lbModel.oneBytheta0;
  double dt      = 1.0;
  double tauNdim = tau/dt;
  double beta    = 1.0/(2.0*tauNdim+1.0);

  int simulationTime = 1*((int)(SIZE_X/lbModel.inlet));

  //const int objectOriginX = 100;
  //  const int objectOriginY = 120;
  //  const int objectLengthX = 10;
  //  const int objectLengthY = 10;
  
//  double   fWallLeft[objectLengthY+2][9];
//  double  fWallRight[objectLengthY+2][9];
//  double    fWallTop[objectLengthX+2][9];
//  double fWallBottom[objectLengthX+2][9];
  //double fwall[2*r+3][2*r+3][N_DV] = {0.0};
  
  printVtk(lbModel,lbGrid,0,beta,dt);
  for(int time=1;time<=simulationTime;time++)
  {
      Collide(lbModel,lbGrid,2.0*beta,time);
    // entropicCollide(lbModel,lbGrid,beta);
   //printstrouhal(lbModel,lbGrid,time,beta,dt);
      if(time%1000==0) {
          //check(lbModel,lbGrid,time);
          printVtk(lbModel, lbGrid, time, beta, dt);
          massConservationCheck(lbGrid, time);
          //print(lbModel, lbGrid);
      }
      periodicAll(lbModel,lbGrid);
//     periodicAll(lbModel,lbGrid);
      prepareWallTopBottom(lbModel,lbGrid);
      prepareWallInletOutlet(lbModel,lbGrid);
    //prepareWallFrontBack(lbModel,lbGrid);
//     prepareWallObject(lbModel,lbGrid,fWallLeft,fWallRight,fWallTop,fWallBottom,objectOriginX,objectOriginY,objectLengthX,objectLengthX);
   // prepareWallCylinder(lbModel, lbGrid, fwall, objectOriginX, objectOriginY, r);

    //massConservationCheck(lbGrid, time);
      advection(lbGrid);
      //printVtk(lbModel, lbGrid, time, beta, dt);
      applyWallTopBottom(lbModel,lbGrid);
    //massConservationCheck(lbGrid, time);
      applyWallInletOutlet(lbModel,lbGrid);
  //  applyWallFrontBack(lbModel,lbGrid);
      //printVtk(lbModel, lbGrid, time, beta, dt);
//     applyWallObject(lbModel,lbGrid,fWallLeft,fWallRight,fWallTop,fWallBottom,objectOriginX,objectOriginY,objectLengthX,objectLengthX);
    //applyWallCylinder(lbModel, lbGrid, fwall, objectOriginX, objectOriginY, r);

  }
    printVtk(lbModel, lbGrid, 100000000, beta, dt);

    // print(lbModel,lbGrid);
  return 0;
}
/////////////////////////////////////
//
//template<typename dataType>
//inline void  prepareOutletFgrads(grid2D<dataType> &lbNode, dataType Pout, dataType dt, dataType tau){
//
//    dataType fEq[NUM_DV],fTemp[NUM_DV];
//    dataType rho, ux, uy, T, Pxx, Pyy, Pxy;
//
//    dataType dtByTwoTau = 0.5*dt/tau;
//
//    for (int i2=0 ;i2<= lbNode.numX2()+1;  i2++) {
//
//        for(int dv=0;dv<NUM_DV;dv++)
//            fTemp[dv]=lbNode(lbNode.numX1()-1,i2,dv);
//
//        getMoment(rho, ux, uy, T,  fTemp);
//        getStress( Pxx, Pyy, Pxy, fTemp);
//
//        rho = Pout/T;
//        getGradDist( fEq, rho, ux, uy, 0.5*rho*(ux*ux+uy*uy) + Pout , Pxx-Pyy, Pxy)    ;
////           getGradDist( fEq, rho, ux, uy, Pout, Pxx-Pyy, Pxy)    ;
//
//        for(int dv=0;dv<NUM_DV;dv++)
//        {
//            lbNode(lbNode.numX1()+1,i2,dv)= fEq[dv];
//        }
//
//    }
//
//
//}

//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//template<typename dataType>
//inline void getGradDist(dataType (&fEq)[N_DV], dataType rho, dataType  ux, dataType  uy, dataType P, dataType pXXminuspYY, dataType pXY){
//
//    dataType wZ = 4.0/9.0, wS = 1.0/9.0, wB = 1.0/36.0, tmp1 = 4.5;
////   dataType tmp = rho*2.5 - 1.5*(pXX + pYY)  ;
////   dataType jx = 3.0*rho*ux;
////   dataType jy = 3.0*rho*uy;
////       fEq[dv_ZERO_ZERO] =  wZ * tmp;
////       tmp -= 1.5*rho ;
////       fEq[dv_P1_ZERO]   = wS * (tmp + jx + tmp1*pXX);
////       fEq[dv_ZERO_P1]   = wS * (tmp + jy + tmp1*pYY);
////       fEq[dv_M1_ZERO]   = wS * (tmp - jx + tmp1*pXX);
////       fEq[dv_ZERO_M1]   = wS * (tmp - jy + tmp1*pYY);
////
////       tmp -= 1.5*rho ;
////       double dot = jx + jy;
////       double dot1= pXX + pYY + 2.0*pXY ;
////       fEq[dv_P1_P1]     = wB *(tmp + dot + tmp1*dot1);
////       fEq[dv_M1_M1]     = wB *(tmp - dot + tmp1*dot1);
////
////       dot = jx - jy;
////       dot1= pXX + pYY - 2.0*pXY ;
////       fEq[dv_P1_M1]    = wB *(tmp + dot + tmp1*dot1);
////       fEq[dv_M1_P1]    = wB *(tmp - dot + tmp1*dot1);
//
//
//    dataType weight[NUM_DV], cx[NUM_DV], cy[NUM_DV], cc, cx2, cy2 ;
//    weight[dv_ZERO_ZERO] = wZ ;  cx[dv_ZERO_ZERO] =  0.0 ;   cy[dv_ZERO_ZERO] =  0.0 ;
//    weight[dv_P1_ZERO  ] = wS ;  cx[dv_P1_ZERO  ] =  1.0 ;   cy[dv_P1_ZERO  ] =  0.0 ;
//    weight[dv_ZERO_P1  ] = wS ;  cx[dv_ZERO_P1  ] =  0.0 ;   cy[dv_ZERO_P1  ] =  1.0 ;
//    weight[dv_M1_ZERO  ] = wS ;  cx[dv_M1_ZERO  ] = -1.0 ;   cy[dv_M1_ZERO  ] =  0.0 ;
//    weight[dv_ZERO_M1  ] = wS ;  cx[dv_ZERO_M1  ] =  0.0 ;   cy[dv_ZERO_M1  ] = -1.0 ;
//    weight[dv_P1_P1    ] = wB ;  cx[dv_P1_P1    ] =  1.0 ;   cy[dv_P1_P1    ] =  1.0 ;
//    weight[dv_M1_P1    ] = wB ;  cx[dv_M1_P1    ] = -1.0 ;   cy[dv_M1_P1    ] =  1.0 ;
//    weight[dv_M1_M1    ] = wB ;  cx[dv_M1_M1    ] = -1.0 ;   cy[dv_M1_M1    ] = -1.0 ;
//    weight[dv_P1_M1    ] = wB ;  cx[dv_P1_M1    ] =  1.0 ;   cy[dv_P1_M1    ] = -1.0 ;
//
//    dataType dot(0.0), T(1.0/3.0), Tinv(3.0);
//
//    for(int dv=0;dv<NUM_DV;dv++)
//    {
//        dot = (ux*cx[dv] + uy*cy[dv])*Tinv*rho;
//        cx2 = cx[dv]*cx[dv];
//        cy2 = cy[dv]*cy[dv];
//        cc  = cx2+cy2;
//        fEq[dv] = weight[dv]*(rho + dot + 0.5*Tinv*Tinv*((P-rho*T)*(cc-2.0*T)+0.5*pXXminuspYY*(cx2-cy2)+2.0*pXY*cx[dv]*cy[dv]));
//    }
//}
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//template<typename dataType>
//inline void getStress( dataType  &Pxx,dataType  &Pyy,dataType  &Pxy, dataType (&f)[NUM_DV]){
//
//
//    Pxx = f[dv_P1_ZERO]+f[dv_M1_ZERO];
//    Pyy=  f[dv_ZERO_P1]+f[dv_ZERO_M1];
//
//    dataType  minus = f[dv_P1_P1]-f[dv_M1_M1];
//    Pxy = f[dv_P1_P1]+f[dv_M1_M1];
//    dataType    den =Pxy+f[dv_P1_M1]+f[dv_M1_P1];
//    Pxx += den;
//    Pyy += den;
//
//    minus = -f[dv_P1_M1]-f[dv_M1_P1];
//    Pxy +=minus;
//
//}
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//template<typename double>
//inline void getMoment(ModelD2Q9<double> &myModel,Grid &gridLB,int x_coord,int y_coord, int z_coord){
//
//    dataType den, energy;
//    rho  = gridLB[ZERO_ZERO];
//    energy = gridLB[PLUS_ZERO]+gridLB[MINUS_ZERO]+gridLB[ZERO_PLUS]+gridLB[ZERO_MINUS];
//    rho +=energy;
//    ux  = gridLB[PLUS_ZERO] - gridLB[MINUS_ZERO];
//    uy  = gridLB[ZERO_PLUS] - gridLB[ZERO_MINUS];
//    dataType  minus = gridLB[PLUS_PLUS] - gridLB[MINUS_MINUS];
//    den = gridLB[PLUS_PLUS] + gridLB[MINUS_MINUS] + gridLB[PLUS_MINUS] + gridLB[MINUS_PLUS];
//    rho +=den;
//    energy += 2.0*den;
//    ux += minus;
//    uy += minus;
//    minus = gridLB[PLUS_MINUS]-gridLB[MINUS_PLUS];
//
//    ux += minus;
//    uy -= minus;
//    dataType rhoInv = 1.0/rho;
//    dataType j2 = (ux*ux+uy*uy)*rhoInv;
//    T = (energy -j2)*rhoInv*0.5 ;
//    ux *=rhoInv;
//    uy *=rhoInv;
//
//}
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//template<typename dataType>
//inline void  applyOutlet(  grid2D<dataType> &lbNode ){
//
//    for (int i2=1 ; i2<= lbNode.numX2();  i2++) for(int dv=0;dv<NUM_DV;dv++){
////           lbNode(lbNode.numX1(),i2,dv) = lbNode(lbNode.numX1()+1,i2,dv) ;
//            lbNode(lbNode.numX1(),i2,dv_M1_ZERO) = lbNode(lbNode.numX1()+1,i2,dv_M1_ZERO) ;
//            lbNode(lbNode.numX1(),i2,dv_M1_P1  ) = lbNode(lbNode.numX1()+1,i2,dv_M1_P1  ) ;
//            lbNode(lbNode.numX1(),i2,dv_M1_M1  ) = lbNode(lbNode.numX1()+1,i2,dv_M1_M1  ) ;
//
////           lbNode(lbNode.numX1(),i2,dv_ZERO_ZERO) = lbNode(lbNode.numX1()+1,i2,dv_ZERO_ZERO) ;
////           lbNode(lbNode.numX1(),i2,dv_ZERO_P1) = lbNode(lbNode.numX1()+1,i2,dv_ZERO_P1) ;
////               lbNode(lbNode.numX1(),i2,dv_ZERO_M1) = lbNode(lbNode.numX1()+1,i2,dv_ZERO_M1) ;
//
//            //           lbNode(lbNode.numX1(),i2,dv_P1_ZERO) = lbNode(lbNode.numX1()+1,i2,dv_P1_ZERO) ;
////               lbNode(lbNode.numX1(),i2,dv_P1_P1  ) = lbNode(lbNode.numX1()+1,i2,dv_P1_P1  ) ;
////               lbNode(lbNode.numX1(),i2,dv_P1_M1  ) = lbNode(lbNode.numX1()+1,i2,dv_P1_M1  ) ;
//
//        }
//}