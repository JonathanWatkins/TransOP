/*	
 * 	This program takes a space delimited text file of
 * 	x and y positions and returns a file containing
 *  the translational order parameter - TOPresult.txt.
 * 
 *  The calculations are based on real space lattice
 *  spacings of a0=1e-7m
 * 
 *  The range/resolution in reciprocal space is 
 *  -3pi to 3pi in qx and qy in steps of pi/20.0
 * 
 * 
 *  Modifications need to be made to this program to allow
 *  for different densities and parameters to be specified.
 * 
 *  
 *  
 * 	This was written by Jon Watkins - Jan 2014
 * 
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <list>
#include <sstream>
#include <iomanip>
#include <string>
//#include <stdio>
#include <cilk/cilk.h>


using namespace std;

double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;
double a0=1;
int t=0;
class CVortex {
    double x, y, z;
    
    
  public:
	CVortex() {
		x=0;
		y=0;
		z=0;
		
	}
	
	~CVortex() {}
	
	
	
	void set_pos (double a,double b,double c)	{
			x = a;
			y = b;
			z = c;
	};
	
	   
    double get_x ()	{
			return x;
	};
    
    double get_y ()	{
			return y;
	};
   
    double get_z ()	{
			return z;
	};
	  
};

struct recipPoint
{
	double qx;
	double qy;
	double sf;
	
};


void initialiseVortices(list<CVortex>& vorticesList, bool& file, std::string pos_file_name_){
	
	ostringstream oss;
	char renderChar[100];
	string renderStr;
	oss.str("");
	oss << pos_file_name_;
	std::cout << "oss: " << oss.str() << endl;
	renderStr = oss.str();
	cout << renderStr << endl;
	strcpy(renderChar,renderStr.c_str());
			
  ifstream myfile (renderChar);
  cout << "initialise Vortices" << endl;
  
	if (myfile.is_open()) {
		cout << "Initial Vortex Positions From File" << endl;
		file = true;
		double xval;
		double yval;
    
    while ( myfile.good() )
    {
      myfile >> xval;
      myfile >> yval;
     
      
      CVortex newVortex;
      //if (xval > 10e-7 && xval < 60e-7 && yval > 2e-7 && yval < 50e-7 ) {
			newVortex.set_pos(xval/a0,yval/a0,0);
			
			vorticesList.push_back(newVortex);
			
    }
    myfile.close();

	}
	
   
  cout << "initialised" << endl;
}

bool readSingleDataStep(std::list<CVortex> &posdata_, std::string &pos_file_name_)
{
	posdata_.clear();
	static std::ifstream readStepDataFile;
	
	if (!readStepDataFile.is_open())
	{
		//std::ostringstream jobstr;
		//jobstr.str("");
		//jobstr << "guidata.txt"; 
	
		readStepDataFile.open(pos_file_name_.c_str());
		std::cout << "Opening guidata file...";	
	}
	
	std::string temp;
	std::string input;
	std::string in;
	std::string in2;
	std::string tline;
	std::string vortline;
	int timestep;
	int numVortices;
	std::stringstream sinput;
	
	double xval;
	double yval;
	int coord_num;
	double a0;
				
							
	if (readStepDataFile.is_open() && readStepDataFile.good()) 
	{
		// read timestep header
	
		std::getline(readStepDataFile, input);
		
		sinput.str("");
		sinput << input;
		sinput >> in >> timestep >> in2 >> numVortices;
		t=timestep;
		//std::cout << in << " " << timestep << " " << in2 << " " << numVortices << std::endl;
		
		
		for (int i =1;i<=numVortices;i++)
		{
			
			std::getline(readStepDataFile, input);
		
				
			
			std::stringstream tinput;
			tinput << input;
			tinput >> xval >> yval >> coord_num >> a0;
				
		
			// add new vortex
			CVortex newVortex;
			newVortex.set_pos(xval,yval,0);
			posdata_.push_back(newVortex);
		
			
		}
		
	
	}
	else return false;
	
	return true;
}



int main(int argc, char *argv[])
{
	std::string pos_file_name;
	if(argc!=2)
	{
		std::cout << "The program should be run with the command line \n"
							<< "   TransOP <pos_file_name>\n"
							<< " where <pos_file_name> should be replaced by the name of the position data.\n";
		return 1;
	}
	else if (argc==2)
	{ 
			pos_file_name = argv[1];
	}
	
	std::list<CVortex> posdata;
	
  double Gx= 4.08406;  // These should be calculated from the structure factor plots. 
	double Gy= 6.04137;
			
	
	
	std::cout << "G: (" << Gx << ", " << Gy << ")" << std::endl;
	
	std::ofstream TOPfile("TOPresults.txt");
	
	TOPfile << "Translational Order parameter calculated at each timestep\n"
					<< "t TOP\n";
	
	while ( readSingleDataStep(posdata, pos_file_name) )
	{
		int numParticles=posdata.size();
		
		if (numParticles==0) continue;
		
		std::cout << t << std::endl;
		
		double r0x=posdata.begin()->get_x();
		double r0y=posdata.begin()->get_y();
		
		double real_part=0;
		double imag_part=0;
		
		for (std::list<CVortex>::iterator p = posdata.begin(); 
				p != posdata.end(); ++p)
		{
			//std::cout << std::distance(posdata.begin(),p)+1 << std::endl;
			
			real_part+=cos(Gx*(p->get_x()-r0x)+Gy*(p->get_y()-r0y));
			imag_part+=sin(Gx*(p->get_x()-r0x)+Gy*(p->get_y()-r0y));
			
			//std::cout << real_part << ", " << imag_part << std::endl;
		
		}
		real_part/=numParticles;
		imag_part/=numParticles;
		
		double sqrd=real_part*real_part+imag_part*imag_part;
		
		std::cout << sqrd << std::endl;
		TOPfile << t << " " << sqrd << std::endl;
	
	
		
	}
  TOPfile.close();
	
	return 0;
	
}
