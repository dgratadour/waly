plug_in,"waly";

extern dwt2d_forward_bspline
/*PROTOTYPE
  void dwt2d_forward_bspline(double array data,int order,int n,int ns)
*/
extern dwt2d_inverse_bspline
/*PROTOTYPE
  void dwt2d_inverse_bspline(double array data,int order,int n,int ns)
*/
extern dwt2d_forward_haar
/*PROTOTYPE
  void dwt2d_forward_haar(double array data,int n,int ns)
*/
extern dwt2d_inverse_haar
/*PROTOTYPE
  void dwt2d_inverse_haar(double array data,int n,int ns)
*/
extern dwt2d_forward_daub
/*PROTOTYPE
  void dwt2d_forward_daub(double array data,int order,int n,int ns)
*/
extern dwt2d_inverse_daub
/*PROTOTYPE
  void dwt2d_inverse_daub(double array data,int order,int n,int ns)
*/
extern dwt2d_filt_daub
/*PROTOTYPE
  void dwt2d_filt_daub(double array data,long array mask,int order,int n,int ns)
*/

func waly_dwt2d(data,direction,family=,order=,ns=)
{
  /* Direction = 1 -> Forward
     Direction = -1 -> Inverse
     Default family: Daubechie and default order: 20
  */

  if (ns == []) ns = 0;
  if (direction == []) direction = 1;
  if (abs(direction) != 1)
    error,"direction must be 1 (forward DWT) or -1 (inverse DWT)";

  if (family == []) family = "daubechie";
  if (allof(["daubechie","bspline","haar"] != family))
    error,"family nust be daublechie, haar or bspline";
  
  if (order == []) {
    if (family == "haar") order = 2;
    if (family == "daubechie") order = 20;
    if (family == "bspline") order = 103;
  }

  size = dimsof(data);
  if (size(1) != 2) error,"data must be a 2D array";
  n1 = size(2);
  n2 = size(3);

  if (n1 != n2) error,"data must be a square array";
  n = n1;
  
  if ((2^(long(log(n)/log(2))) != n) || (2^(long(log(n2)/log(2))) != n2))
    error,"The Dimensions of Data must be a power of 2";

  data = double(data);
  
  if (family == "daubechie") {
    if ((order < 4) || (order > 20) || (order % 2 != 0))
      error,"For Daubechie DWT, order must be > 4 and < 20 and even";
    if (direction == 1) res = dwt2d_forward_daub(data,int(order),n,int(ns));
    else res = dwt2d_inverse_daub(data,int(order),n,int(ns));
  }
  if (family == "bspline") {
    tabOrder = [103,105,202,204,206,208,301,303,305,307,309];
    if (allof(tabOrder != order))
      error,"For bspline DWT, order must be 103, 105, 202, 204, 206, 208, 301, 303, 305 307 or 309";
    if (direction == 1) res = dwt2d_forward_bspline(data,int(order),n,int(ns));
    else res = dwt2d_inverse_bspline(data,int(order),n,int(ns));
  }
  if (family == "haar") {
    if (order != 2)
      error,"For Haar DWT, order must be 2";
    if (direction == 1) res = dwt2d_forward_haar(data,n,int(ns));
    else res = dwt2d_inverse_haar(data,n,int(ns));
  }
  
  return res;
}

func waly_denoise_dwt2d(data,family=,order=,ns=,percentile=,threshold=,ratio=,hardCutOff=)
{

  if (percentile == []) percentile = 0;
  if (threshold == []) threshold = 0;
  if (ratio == []) ratio = 0;
  if (hardCutOff == []) hardCutOff = 0;
  
  waveCoeff = waly_dwt2d(data,1,family=family,order=order,ns=ns);

  if (percentile != 0) {
    index = (sort(waveCoeff(*))(::-1))(1:long(numberof(waveCoeff)*percentile));
    mask = waveCoeff*0.;
    mask(index) = 1.;
  } else {
    if (threshold != 0) mask = abs(waveCoeff) > threshold;
    else {
      if (ratio != 0) mask = abs(waveCoeff) > max(abs(waveCoeff))/ratio;
      else {
        if (hardCutOff != 0) {
          size = dimsof(waveCoeff)(2);
          // fixme : may need recentering using a roll
          xx = ((indgen(size)/size)-0.5)(,-:1:size);
          mask = sqrt(xx^2 + transpose(xx)^2) > hardCutOff;
          //mask = roll(dist(size))/(size/2.) > hardCutOff;
        } else mask = 1.;
      }
    }
  }

  waveCoeff *= mask;

  dataFilt = waly_dwt2d(waveCoeff,-1,family=family,order=order,ns=ns);

  return dataFilt;
}

func waly_atrous(image,filter=,nScales=)
{
  if (filter == []) filter = 1./[16, 4, 8/3., 4, 16];

  fmat = filter(,-)(,+)*transpose(filter(,-))(+,);
  sz = dimsof(image);

  nScales = long(floor(log(min([sz(2),sz(3)])/numberof(filter))/log(2.)));
  decomp = array(float,[3,sz(2),sz(3),nScales+1]);

  im = image;
  for (k=0;k<=nScales-1;k++) {
    smooth = yeti_convolve(im, kernel=fmat); // fixme : require yeti !
    decomp(,,nScales+1-k) = im-smooth;
    im = smooth;
    newfilter = array(float,2*numberof(filter)-1); 
    newfilter(2*indgen(numberof(filter))-1) = filter;
    fmat = newfilter(,-)(,+)*transpose(newfilter(,-))(+,);
    filter = newfilter;
  }

  decomp(,,1) = smooth;

  return decomp;
}

