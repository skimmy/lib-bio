#ifndef FILTER_BASE_H
#define FILTER_BASE_H

class FilterBase {
public:
  virtual bool apply(Read& read);
};

#endif
