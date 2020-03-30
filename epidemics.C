#include "epidemics.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TMath.h"

TRandom3 rGen(0);

//2D maxwell distribution
TF1 maxwell("maxwell","x/[0]*TMath::Exp(-0.5*x*x/[0])",0,100);
TF1 logNormal("lognormal","TMath::LogNormal(x,[0],0.,[1])",5,25);

//Model a person as a particle in a 2D gas at temperature T
Person::Person(float pInfect, float initialP, float dieP, float recoT, float temperature)
{
  _position[0] = rGen.Rndm();
  _position[1] = rGen.Rndm();
  float theta = rGen.Rndm()*2*TMath::Pi();
  maxwell.SetParameter(0,temperature);
  float velocity = maxwell.GetRandom();
  _vel[0] = velocity*TMath::Cos(theta);
  _vel[1] = velocity*TMath::Sin(theta);
  _P = initialP;
  _ill = false;
  setAsIll();
  _P = pInfect;
  //logNormal.SetParameter(1,recoT);
  //_recoveryTime = logNormal.GetRandom();
  _recoveryTime = rGen.Gaus(recoT,recoT*0.5); //randomise recovery
  if (_recoveryTime<0)
    _recoveryTime=0;
  if (_ill)
    _illtime=rGen.Rndm()*_recoveryTime; //simulate uniform initial illtime [0,recoT]
  _dieP = dieP;
  _nMeets = 0;
  _nInfected = 0;
}

void Person::moveToTown(Town* town)
{
  _town = town;
  int s = town->_size - 1;
  _position[0] = _position[0] * s;
  _position[1] = _position[1] * s;
  _town->putPerson(this);
}
  
void Person::setAsIll()
{
  if (!_ill && rGen.Rndm() < _P)
    {
      _ill = 1;
      _illtime = 0;
      if (rGen.Rndm() < _dieP)
	{
	  //rough numbers 
	  logNormal.SetParameter(0,0.2);
	  logNormal.SetParameter(1,11);
	  _dieTime = logNormal.GetRandom();
	}
      else
	{
	  _dieTime = 0;
	  _dieP =0;
	}
    }
}

void Person::setAsHealthy()
{
  //once recovered it cannot be infected anymore                                                                                                          
  _ill = 0;
  _illtime = 0;  
  _P = 0;
}

void Person::tryToRecover(float dt)
{
  //check if it recovers                                                                                                                                  
  if (_ill)
    {
      _illtime+=dt;
      if (_dieTime==0 && _illtime>_recoveryTime) //recovery, if not dying
	setAsHealthy();
    }
}

void Person::bounce()
{
  float s = _town->_size-1;
  for (int id=0; id<2;++id)
    {
      //check if it is at the border                                                                                                                          
      if (_position[id] >= s)
	{
	  _position[id] = s;
	  _vel[id] = -_vel[id];
	}
      if (_position[id] < 0)
	{
	  _position[id] = 0;
	  _vel[id] = -_vel[id];
	}
    }
}

void Person::meet(Person* p,float dt)
{
  if (p != NULL)
    {
      float d=0;
      for (int id=0;id<2;++id)
	  d+=TMath::Power(_position[id]-p->_position[id],2);
      d=TMath::Sqrt(d);
      
      float vrel=0;
      for (int id=0;id<2;++id)
	vrel+=TMath::Power(_vel[id]-p->_vel[id],2);
      vrel=TMath::Sqrt(vrel);

      float vcm[2];
      float v1[2];
      float v2[2];
      
      for (int id=0;id<2;++id)
	{
	  vcm[id]= 0.5*(_vel[id]+p->_vel[id]);
	  v1[id]=-(_vel[id]-vcm[id]);
	  v2[id]=-(p->_vel[id]-vcm[id]);
	  v1[id]+=vcm[id];
	  v2[id]+=vcm[id];
	  _vel[id] = v1[id];
	  p->_vel[id] = v2[id];
	}

      _nMeets++;
      //check if it  can transmit infection
      if (_ill)
	{
	  bool pStatus=p->_ill;
	  p->setAsIll();
	  if (!pStatus && p->_ill)
	    _nInfected++;
	}
    }
}

void Person::step(float dt)
{
  // remove it from the current position                                                                                                                   
  _town->removePerson(this);
  
  // compute next position
  for (int id=0;id<2;++id)
    _position[id] += dt*_vel[id];
  bounce();

  // move it to the current position                                                                                                                       
  _town->putPerson(this);
  
  // check if there is anyone in the vicinity                                                                                                              
  int x = int(_position[0]);
  int y = int(_position[1]);
  int s = _town->_size - 1;
  for (int ix=max(0, x-1);ix<min(s,x+1);++ix)
    for (int iy=max(0, y-1);iy<min(s,y+1);++iy)
      if (ix != x && iy != y)
	{
	  Person* p = (*_town)._map[ix*_town->_size+iy];
	  meet(p,dt);
	}

  // perform a simulation step
  tryToRecover(dt);  
  //check if dead
  if (_ill && _dieTime>0)
    if (_illtime >_dieTime)
      _town->deletePerson(this);
}


Town::Town(int size,std::string name)
{
  _name = name;
  _size=size;
  _map.resize(_size*_size); //double precision
}

Town::Town(std::string name, int N, float initInfected, float temperature, int size, float recoT, float P, float dieP) :
Town(size,name)
{
  for (int i=0;i<N;++i)
    addPerson(new Person(P,initInfected,dieP,recoT,temperature));
}

void Town::addPerson(Person* person)
{
  _population.insert(person);
  person->moveToTown(this);
}

void Town::deletePerson(Person* person)
{
  _population.erase(person);
  removePerson(person);
  //  delete person;
}

void Town::putPerson(Person* person)
{
  int x = int(person->_position[0]);
  int y = int(person->_position[1]);
  _map[x*_size+y] = person;
}

void Town::removePerson(Person* person)
{
  int x = int(person->_position[0]);
  int y = int(person->_position[1]);
  _map[x*_size+y] = NULL;
}

void Town::step(float dt)
{
  for (auto& p: _population)
    p->step(dt);
}

int Town::nIll()
{
  int n = 0;
  for (auto& p: _population)
    if (p->_ill)
      n++;
  return n;
}

int Town::nRecovered()
{
  int n = 0;
  for (auto& p: _population)
    if (p->_P==0)
      n++;
  return n;
}

int Town::nSusceptible()
{
  int n = 0;
  for (auto& p: _population)
    if (p->_ill==0 && p->_P>0)
      n++;
  return n;
}

long Town::nMeet()
{
  long n = 0;
  for (auto& p: _population)
      n += p->_nMeets;
  return n;
}

float Town::nTransmissionPerInfected()
{
  float nTot = 0;
  int n = 0;
  for (auto& p: _population)
    if (p->_ill && p->_illtime>0. && p->_recoveryTime>0)
      {
	if ((p->_illtime/p->_recoveryTime)>0.1)
	  {
	    nTot += p->_nInfected/(p->_illtime/p->_recoveryTime);
	    n++;
	  }
      }
  
  if (n>0)
    return float(nTot)/float(n);
  else
    return 0.;
}
