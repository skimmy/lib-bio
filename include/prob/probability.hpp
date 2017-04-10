// probabilitiesx.hpp

// Copyright 2017 Michele Schimd

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LBIO_PROBABILITY_HPP
#define LBIO_PROBABILITY_HPP

#include <lbio.h>

#include <include/util/foreach.hpp>

#include <algorithm>
#include <unordered_map>
#include <random>

DEFAULT_NAMESPACE_BEGIN

/**
   \brief Representation of a discrete probability space.

   A discrete probability space is a set of events (\c event_type),
   with a probability measure P. The space is \a immutable in the
   sense that the assigned values can not be changed after the
   construction (for this reason no default empty constructor is
   provided).

   \note This class does not check that the values passed actually
   define a valid probability function (or even a measure). This is
   done for performance reasons, the user should guarantee that all
   the invariants hold. In any case the class acting as a container
   will continue to work even when the content is not semantically
   correct.

   \tparam _Event  The type of the events
   \tparam _Prob   The numeric type used for probability (\c double
                   by default)
 */
template <typename _Event, typename _Prob = double>
class DiscreteProbability
{
public:
  using event_type = _Event;
  using probability_type = _Prob;

  using container = std::unordered_map<event_type, probability_type>;

  template <typename _Iter>
  DiscreteProbability(_Iter _beg, _Iter _end)
    : _events{_beg,_end} { }

  /**
     \brief Gives the probability of the elementary event \c e (0
     if such event does not exist)

     \note Rather than failing for events \c e not present as 
     elementary events, the function returns 0.
   */
  probability_type inline
  probability_of(event_type e) const {
    return _events.count(e) ? _events.at(e) : probability_type{0};
  }


  /**
     \brief Returns the sum of probabilities for every element in 
     the range defined by the iterators.

     \tparam _Iter  The iterator type used to define tha range

     \note There is no check about the validity of value in the range
     for example duplicated contribute to the final sum as many times
     as they appear in the range. This unchecked behavior is kept in
     the spirit of efficiency, users must take care of correct semantic.
   */
  template <typename _Iter>
  probability_type inline
  probability_of(_Iter _beg, _Iter _end) const {
    probability_type p {0};
    std::for_each(_beg, _end,
    		  [&p,this](event_type e) { p += probability_of(e); });
    return p;
  }

  /**
     \brief Returns the sum of probabilities for events that satisfy 
     the given predicate.

     \tparam _Pred  The predicate used as test

     \note This member function can be used to implement various
     concept like marginal probability and conditional probability.
     This makes use of \c lbio::for_each_such_that template function.
   */
  template <typename _Pred>
  probability_type inline
  probability_if(_Pred predicate) const {
    probability_type p {0};
    lbio::for_each_such_that(_events.cbegin(), _events.cend(),
			     [=](entry_type e) {return predicate(e.first);},
     			     [&p,this](entry_type e) {
			       p += probability_of(e.first);
			     }); 
    return p;
  }

  /**
     \brief Returns a copy of the underlying container which is of
     type of \c container.
   */
  const container
  get_container() const {
    return _events;
  }

private:
  using entry_type = std::pair<event_type, probability_type>;

  container _events;
};

DEFAULT_NAMESPACE_END

#endif
