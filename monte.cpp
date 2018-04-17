#include "monte.h"
#include "parameter.h"
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <stdlib.h>
void atom::setx(double x1){
	x=x1;
}
void atom::sety(double x2){
	y=x2;
}
void atom::setr(double ra){
	radius=ra;
}
void atom::setstress_tensor(std::vector<double> a){
        stresstensor=a;
}
void atom::setm(double m){
   mass=m;
}
std::vector<double> atom::getstress(){
	return stresstensor;
}
void atom::setv(std::vector<double> v){
        speed=v;
}
void atom::setf(std::vector<double> f){
   force=f;
}
double atom::getx(){
	return x;
}
double atom::gety(){
	return y;
}
double atom::getr(){
	return radius;
}
void atom::printneighbor(){
		std::cout<<"the neighbors are:"<<std::endl;
		for(std::list<int>::iterator a=neighbor.begin();a!=neighbor.end();a++){
			std::cout<<*a<<" "<<std::endl;
		}
	}
void atom::printstress(){
	std::cout<<stresstensor[0]<<" "<<stresstensor[1]<<" "<<stresstensor[2]<<" "<<stresstensor[3]<<std::endl;
}
void atom::printinfo(){
	std::cout<<x<<" "<<y<<" "<<speed[0]<<" "<<speed[1]<<" "<<force[0]<<" "<<force[1]<<std::endl;
}
void atom::printforce(){
    std::cout<<"F(x):"<<force[0]<<std::endl;
    std::cout<<"F(y):"<<force[1]<<std::endl;
}
double distance(atom& one,atom& two){
	double r=(one.x-two.x)*(one.x-two.x)+(one.y-two.y)*(one.y-two.y);
	return sqrt(r);
}
double dxx(atom& one,atom& two){
	return one.x-two.x;
}
double dyy(atom& one,atom& two){
	return one.y-two.y;
}
double potential(atom& one,atom& two){
	double r=(one.x-two.x)*(one.x-two.x)+(one.y-two.y)*(one.y-two.y);
	r=sqrt(r);
	if(r<r0){
		return 4*eps*(pow(sigma/r,12)-pow(sigma/r,6));
	}
	else if(r>r_cut){
		return 0;
	}
	else{
		return A*pow(r-r_cut,3)+B*pow(r-r_cut,2);
	}
}
double allpotential(std::vector<atom>& allatom){
    int size=allatom.size();
    double sum=0.0;
    double temp=0.0;
    for(size_t i=0;i<size;i++){
       for(std::list<int>::iterator a=allatom[i].neighbor.begin();a!=allatom[i].neighbor.end();a++){
        temp=potential(allatom[i],allatom[*a]);
        sum=sum+temp;
       }
    }
    return sum/2;
}
atom& atom::operator =(atom& one){
	this->mass=one.mass;
	this->x=one.x;
	this->y=one.y;
	this->speed=one.speed;
  this->force=one.force;
  this->neighbor=one.neighbor;
	return *this;
}
double allener(std::vector<atom>& allatom){
    double poten=allpotential(allatom);
    double kinet=0.0;
    for(size_t i=0;i<allatom.size();i++){
        kinet=kinet+0.5*(allatom[i].mass)*(allatom[i].speed[0]*allatom[i].speed[0]+allatom[i].speed[1]*allatom[i].speed[1]);
    }
    return kinet+poten;
}
//the force exerted on one by two
std::vector<double> str_tensor(atom& one,atom& two){
    std::vector<double> a(4,0);
    double r=distance(one,two);
    double deri=0;
    if(r<r0){
        deri=24*eps/r*(-2*pow(sigma/r,12)+pow(sigma/r,6));
    }
    else if(r<r_cut){
        deri=3*A*(r-r_cut)*(r-r_cut)+2*B*(r-r_cut);
    }
    else{
        deri=0.0;
    }
    a[0]=deri*(one.x-two.x)*(one.x-two.x)/r;
    a[1]=deri*(one.x-two.x)*(one.y-two.y)/r;
    a[2]=deri*(one.x-two.x)*(one.y-two.y)/r;
    a[3]=deri*(one.y-two.y)*(one.y-two.y)/r;
    return a;
}
std::vector<double> cal_force(atom& one,atom& two){
   std::vector<double> a(2,0);
   double r=distance(one,two);
   double deri=0;
   if(r<r0){
      deri=24*eps/r*(-2*pow(sigma/r,12)+pow(sigma/r,6));
   }
   else if(r<r_cut){
      deri=3*A*(r-r_cut)*(r-r_cut)+2*B*(r-r_cut);
   }
   else{
      deri=0.0;
   }
   a[0]=-1*deri*(one.x-two.x)/r;
   a[1]=-1*deri*(one.y-two.y)/r;
   return a;
}
std::vector<double>& operator +=(std::vector<double>& one,std::vector<double>& two){
   for(size_t i=0;i<one.size();i++){
      one[i]=one[i]+two[i];
   }
   return one;
}
double hydropressure(std::vector<atom> atomall){
    int size=atomall.size();
    double xx=0.0;
    double yy=0.0;
    for(size_t i=0;i<size;i++){
        xx=xx+atomall[i].stresstensor[0];
        yy=yy+atomall[i].stresstensor[3];
    }
    return xx+yy;
}
void updatetensor(std::vector<atom>& atomall){
    std::vector<double> temp(4,0);
    std::vector<double> all(4,0);
    int size=atomall.size();
    for(size_t i=0;i<size;i++){
        for(size_t k=0;k<4;k++){
            all[k]=0.0;
        }
        for(std::list<int>::iterator a=atomall[i].neighbor.begin();a!=atomall[i].neighbor.end();a++){
            temp=str_tensor(atomall[i],atomall[*a]);
            for(size_t k=0;k<4;k++){
                all[k]=all[k]+temp[k];
            }
        }
    atomall[i].setstress_tensor(all);
    }
}
std::ostream& operator<<(std::ostream& os,atom& output){
		os<<"for atom position ("<<output.x<<","<<output.y<<")"<<std::endl;
		os<<"stress tensor(xx,xy,yx,yy): "<<output.stresstensor[0]<<" "<<output.stresstensor[1]<<" "<<output.stresstensor[2]<<" "<<output.stresstensor[3]<<std::endl;
		return os;
}
std::fstream& operator<<(std::fstream& os,atom& output){
		os<<output.x<<" "<<output.y<<" "<<output.stresstensor[0]<<" "<<output.stresstensor[1]<<" "<<output.stresstensor[2]<<" "<<output.stresstensor[3];
		return os;
}
int count(std::vector<atom> all,atom& input,double r){
    int c=0;
    int size=all.size();
    for(int i=0;i<size;i++){
        if(distance(all[i],input)<r){
            c++;
        }
    }
    return c-1;
}
void print_radial_dis(double r_start,double r_stop,std::vector<atom>& atomall,std::string name){
	 // double r_start=0.0000001;
  //  double r_stop=15;
    int N=10000;
    int size=atomall.size();
    std::vector<double> ra_dis(N,0.0);
    std::vector<double> ra_dis_all(N,0.0);
    std::vector<double> rinter(N,0.0);
    double r_delta=(r_stop-r_start)/N;
    double r_inter=0.0;
    int count_old=0;
    int count_new=0;
    int count_delta=0;
    std::fstream radis;
    radis.open(name,std::fstream::out);
    for(size_t j=0;j<size;j++){
    for(size_t i=0;i<N;i++){
        r_inter=i*r_delta+r_start;
        rinter[i]=r_inter;//record the radius.
        count_new=count(atomall,atomall[j],r_inter);
        count_delta=count_new-count_old;
        ra_dis[i]=count_delta/2/Pi/r_inter/r_delta;
        count_old=count_new;
       // radis<<r_inter<<" "<<ra_dis[i]<<std::endl;
    }
    std::cout<<j<<std::endl;
    ra_dis_all+=ra_dis;
    }
    for(size_t j=0;j<size;j++){
        ra_dis_all[j]=ra_dis_all[j]/size;
    }
    for(size_t i=0;i<N;i++){
        radis<<rinter[i]<<" "<<ra_dis_all[i]<<std::endl;
    }
    radis.close();
}
bool accept(double energydiff,double temp){
    double possi;
    double rand_num;
    if(energydiff>0){
        possi=exp(-1*energydiff/temp/Kb);
				std::cout<<possi<<std::endl;
				rand_num=((double)rand()/(RAND_MAX));
        if(possi>rand_num){
            return true;
        }
        else{
            return false;
        }
    }
    else{
        return true;
    }
}

