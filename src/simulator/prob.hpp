// Class representing empirical distribution as histogram array
class EmpiricalDistribution {  
public:
  EmpiricalDistribution(double a, double b, size_t N);

  // given a sample x in [a,b] returns the index where x must be counted
  size_t indexForSample(double x) const;
  // returns the relative frequenxy (i.e., count[i] / total) for the i-th index
  double valueAtIndex(size_t i) const;
  // returns the number of elements at the i-th index
  double countAtIndex(size_t i) const { return f[i]; }

  // returns the number of intervals 'n'
  double getIntervalCount() const { return n; }
  
  // add a new sample x to the distribution
  void addSample(double x);

  // returns the emprical cumulative distribution derived from
  // currently stored samples.
  void getCDF(std::vector<double>& cdf)  const;
    
  
private:
  std::vector<double> f;
  size_t n;
  double xa;
  double xb;
  double step; // defined as (b-a)/n
  double total;
};

struct ScoreSumFreq
{
  double sSum;
  int sFreq;
};

void initProbabilities();
void clearProbabilities();
double randomReadsOverlapProbNoErr(const std::string& s1, const std::string& s2, size_t s);
double overlappingStringsSum(const std::string & s1, const std::string& s2);
double overlappingStringsSumWithErr(const std::string& s1, const std::string& s2);
size_t generateInterReadDistance();

size_t percentileIndex(const std::vector<double>& cdf, double perc);

double approximatedScore(size_t s, double* num_den);
double approximatedScore(size_t s);

double score(const std::string& r1, const std::string& r2, size_t s);
double scoreExt(const std::string& r1, const std::string& r2, size_t s, double* num_den);
