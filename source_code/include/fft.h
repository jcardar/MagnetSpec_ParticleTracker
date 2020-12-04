/*
****************************************************************************
*                                                                          *
*   FFT Header file, based on code from www.codeproject.com, originally    *
*   based on code from Numerical Recipes with refinement in speed.         *
*   Added functionality to load in complex and real arrays          	   *
*   separately. This requires a number of values which is a power of 2.    *
*   This version is also normalized so that the 1/N is taken into account. *
*   AGRT 2010                                                              *
*                                                                          *
****************************************************************************
 
    data_re -> float array that represent the real array of complex samples
    data_im -> float array that represent the imag array of complex samples
    NVALS -> length of real or imaginary arrays (N^2 order number) 
    isign -> 1 to calculate FFT and -1 to calculate Reverse FFT

    The function returns an integer, 1 if FFT ran, 0 otherwise. It will be 
    zero if NVALS is not a power of 2
*/

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
//tempr is a variable from our FFT function

namespace fftconsts // defined so as not to intefere with other definitions
{
  const double PI=acos(-1.0);
};

int FFT (double data_re[], double data_im[], const unsigned long NVALS, int isign)
{
  // TEST THAT NVALS IS A POWER OF 2
  const int test=!(NVALS&(NVALS-1));
  if (test)
    {
      //variables for trigonometric recurrences
      unsigned long mmax,m,j,istep,i;
      double wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
      /*
	the complex array is real+complex so the array 
	as a size n = 2* number of samples
	real part is data[index] and the complex part is data[index+1]
      */
      const unsigned long n=NVALS << 1; // bitwise multiply by 2
      
      /*
	Pack components into double array of size n - the ordering 
	is data[0]=real[0], data[1]=imag[0] .......
      */
      double *data;
      data = new double[n];
      for (i=0,j=0;j<n;++i,j+=2)
	{
	  data[j]=data_re[i];
	  data[j+1]=data_im[i];
	}
      /*
	binary inversion (note that 
	the indexes start from 1 which means that the
	real part of the complex is on the odd-indexes
	and the complex part is on the even-indexes
      */
      
      j=0;
      for (i=0;i<n/2;i+=2) 
	{
	  if (j > i) 
	    {
	      //swap the real part
	      SWAP(data[j],data[i]);
	      
	      //swap the complex part
	      SWAP(data[j+1],data[i+1]);
	      
	      // checks if the changes occurs in the first half
	      // and use the mirrored effect on the second half
	      if((j/2)<(n/4))
		{
		  //swap the real part
		  SWAP(data[(n-(i+2))],data[(n-(j+2))]);
		  
		  //swap the complex part
		  SWAP(data[(n-(i+2))+1],data[(n-(j+2))+1]);
		}
	    }
	  
	  m=NVALS;
	  while (m >= 2 && j >= m) 
	    {
	      j -= m;
	      m >>= 1; 
	    }
	  j += m;
	}
      
      //Danielson-Lanzcos routine 
      mmax=2;
      while (n > mmax)
	{
	  istep = mmax << 1;
	  theta=isign*(2.0*fftconsts::PI/mmax);
	  wtemp= sin(0.5*theta);
	  wpr = -2.0*wtemp*wtemp;
	  wpi=sin(theta);
	  wr=1.0;
	  wi=0.0;
	  //internal loops
	  
	  for (m=1;m<mmax;m+=2) 
	    {
	      for (i=m;i<=n;i+=istep) 
		{
		  j=i+mmax;
		  tempr=wr*data[j-1]-wi*data[j];
		  tempi=wr*data[j]+wi*data[j-1];
		  data[j-1]=data[i-1]-tempr;
		  data[j]=data[i]-tempi;
		  data[i-1] += tempr;
		  data[i] += tempi;
		}
	      wr=(wtemp=wr)*wpr-wi*wpi+wr;
	      wi=wi*wpr+wtemp*wpi+wi;
	    }
	  mmax=istep;
	}  
      /*
	Return elements to real and complex components
      */
      if (isign>0)
	{
	  for (i=0,j=0;j<n;++i,j+=2)
	    {
	      data_re[i]=data[j];
	      data_im[i]=data[j+1];
	    }
	}
      else 
	{
	  double invNVALs=1.0/NVALS;
	  for (i=0,j=0;j<n;++i,j+=2)
	    {
	      data_re[i]=data[j]*invNVALs;
	      data_im[i]=data[j+1]*invNVALs;
	    }
	}
    }
  return test;
}
