require,"waly.i";

func check_waly(void)
{
  image = fits_read("lena.fits");
  
  daub  = waly_dwt2d(image,1,family="daubechie");
  daubr = waly_dwt2d(daub,-1,family="daubechie");
  write,format="Daubechie reconstruction error : %f\n",
    sqrt(sum(abs(image-daubr)^2)/numberof(image));

  bspl  = waly_dwt2d(image,1,family="bspline");
  bsplr = waly_dwt2d(bspl,-1,family="bspline");
  write,format="BSpline reconstruction error : %f\n",
    sqrt(sum(abs(image-bsplr)^2)/numberof(image));

  /*
  haar  = waly_dwt2d(image,1,family="haar");
  haarr = waly_dwt2d(haar,-1,family="haar");
  write,format="Haar reconstruction error : %f\n",
    sqrt(sum(abs(image-haarr)^2)/numberof(image));
  */

  // heu je sais pas trop comment ca marche ce denoising ... a voir ...
  imnoise  = image+random_n(dimsof(image))*max(image)/10.;
  denoise1 = waly_denoise_dwt2d(imnoise,percentile=0.3);
  denoise2 = waly_denoise_dwt2d(imnoise,threshold=0.3);
  denoise3 = waly_denoise_dwt2d(imnoise,ratio=0.3);
  denoise4 = waly_denoise_dwt2d(imnoise,hardCutOff=0.3);
}


/*
object=png_read("lena.png");
object=object(avg,,);
object=object(,::-1);
fits_write,"lena.fits",object;
*/
