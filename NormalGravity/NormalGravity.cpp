// Example of using the GeographicLib::GravityModel class
// This requires that the egm96 gravity model be installed; see
// https://geographiclib.sourceforge.io/html/gravity.html#gravityinst

#include <iostream>
#include <exception>
#include <GeographicLib/GravityModel.hpp>

using namespace std;
using namespace GeographicLib;

int main() {
  try {
    GravityModel grav("egm2008");
    double lat = 27.99, lon = 86.93, h = 8820; // Mt Everest
    double gx, gy, gz;
    // grav.Gravity(lat,lon, h, gx, gy, gz);
    // cout << gx << " " << gy << " " << gz << "\n";
    while (cin>>lon>>lat>>h)
    {
      grav.Gravity(lat,lon, h, gx, gy, gz);
      cout<<lon<<" "<<lat<<" "<<gz<<endl;
    }
    
    
  }
  catch (const exception& e) {
    cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
}