/*
  --------------------------------------------------
  Weighting functions for ESPIC
  AGRT 2019
  
  Basic linear particle weighting routines 
  --------------------------------------------------
*/

// Function for weighting particles to grid
// particle positions parxVec of size Npar
// charge density grid of size Nx rhoVec
//-----------------------------------------
int WeightRho(double parxVec[], double rhoVec[])
{   
  const double idx = 1.0/dx;
  unsigned long i,j;

  // Initialize
  for (j=0;j<Nx;++j)
    {
      rhoVec[j]=0.0; 
    }

  double weight;

   for (i=0;i<Npar;++i)
    {
      weight = parxVec[i];      
      // This is because I have chosen to start my boundary at -dx/2
      if (weight<0.0) weight += L_sys;  //////WHAT IF MY BOUNDARY IS AT 0???
      
      weight *= idx;
      j= (int) (weight); // cast float -> int
      
      weight -= (double) j; // now weight is difference between xi and xj
      
      // Increment rho
      rhoVec[j] += (1.0 - weight)*Wpar;    
      ++j;
      if (j >= Nx) 
        {
          rhoVec[0] += weight*Wpar;
        }
      else 
        {
          rhoVec[j] += weight*Wpar;
	}
    }
   return 0;
}
//-----------------------------------------


// Function for weighting electric field from 
// Evec to Npar size vector of forces Fvec  
// given particle positions parxVec
//-----------------------------------------
int WeightF(double EVec[],double parxVec[],double FVec[])
{   
  const double idx = 1.0/dx;
  unsigned long i,j;

  double weight;

  for (i=0;i<Npar;++i)
    {
      weight = parxVec[i];      
      // This is because I have chosen to start my boundary at -dx/2
      if (weight<0.0) weight += L_sys; 


      weight *= idx;
      j= (int) (weight); // cast float -> int

      weight -= (double) j; // now weight is difference between xi and xj

      if (j+1 >= Nx) 
        {
          FVec[i] = (1.0 - weight)*Qpar*EVec[j]+weight*Qpar*EVec[0];
        }
      else 
        {
	  FVec[i] = (1.0 - weight)*Qpar*EVec[j]+weight*Qpar*EVec[j+1];
        }
    }
  return 0;
}
//-----------------------------------------

