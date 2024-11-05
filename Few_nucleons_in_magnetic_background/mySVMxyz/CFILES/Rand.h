#ifndef RAND_H
#define RAND_H
class Rand
{
private:
	unsigned long long int u, v, w;
public:
	Rand(unsigned long long int j);
	~Rand();
	unsigned long long int int64();
	double doub();
	unsigned int int32();
};
#endif 
