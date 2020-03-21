#include <set>
#include <string>

class Town;

class Person
{
 public:
  Person(float P, float initialP = 0.03, float dieP = 0.03, float recoveryTime = 100, float temperature = 1);
  
  void moveToTown(Town* town);
  void setAsIll();
  void setAsHealthy();
  void bounce();
  void meet(Person* p,float dt);
  void step(float dt);
  void tryToRecover(float dt);
  
  float _position[2];
  float _vel[2];
  float _temperature;
  Town* _town;
  bool _ill;
  float _illtime;
  float _recoveryTime;
  float _P;
  float _dieP;
  float _dieTime;
  long _nMeets;
  int _nInfected;
};

class Town
{
 public:
  Town(int size, std::string name);
  Town(std::string name="COVIDVille", int N=500, float initInfected=0.05, float temperature=1, int size=100, float recoT=100, float P=0.3, float dieP=0.03);
  
  void addPerson(Person* person);
  void deletePerson(Person* person);
  
  void putPerson(Person* person);
  void removePerson(Person* person);

  void step(float dt);
  int nIll();
  int nRecovered();
  int nSusceptible();
  long nMeet();
  float nTransmissionPerInfected();
  
  std::set<Person*> _population;
  std::string _name;
  //linearised position map
  std::vector<Person*> _map;
  int _size;
};

