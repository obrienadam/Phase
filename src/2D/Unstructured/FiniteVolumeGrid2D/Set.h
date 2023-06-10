#ifndef PHASE_SET_H
#define PHASE_SET_H

#include <string>
#include <unordered_set>
#include <vector>

#include "Types/Types.h"

template <class T> class Set {
public:
  typedef typename std::vector<Ref<const T>>::iterator iterator;

  typedef typename std::vector<Ref<const T>>::const_iterator const_iterator;

  struct EqualTo {
    bool operator()(const T &lhs, const T &rhs) const {
      return lhs.id() == rhs.id();
    }
  };

  struct Hash {
    std::size_t operator()(const T &item) const {
      return std::hash<Label>()(item.id());
    }
  };

  Set(const std::string &name = "") : name_(name) {}

  Set &operator=(const std::vector<Ref<const T>> &rhs);

  //- name
  const std::string &name() const { return name_; }

  //- allocation
  void reserve(Size size);

  //- Size info
  size_t size() const { return items_.size(); }

  bool empty() const { return items_.empty(); }

  //- Access
  const T &operator[](Label i) const { return items_[i]; }

  const std::vector<Ref<const T>> &items() const { return items_; }

  bool isInSet(const T &item) const {
    return itemSet_.find(std::cref(item)) != itemSet_.end();
  }

  //- Add
  virtual bool add(const T &item);

  virtual void add(const_iterator begin, const_iterator end);

  virtual void add(typename std::vector<T>::const_iterator begin,
                   typename std::vector<T>::const_iterator end);

  virtual void add(const Set<T> &set);

  virtual bool remove(const T &item);

  virtual void remove(const_iterator begin, const_iterator end);

  virtual void remove(const Set<T> &set);

  //- Clear
  virtual void clear();

  //- Operators
  Set &operator+=(const Set &rhs);

  Set &operator-=(const Set &rhs);

  //- Iterators
  iterator begin() { return items_.begin(); }

  iterator end() { return items_.end(); }

  const_iterator begin() const { return items_.begin(); }

  const_iterator end() const { return items_.end(); }

protected:
  std::string name_;

  std::vector<Ref<const T>> items_; // Used for faster iteration over all cells

  std::unordered_set<Ref<const T>, Hash, EqualTo>
      itemSet_; // Allows cell lookup via an id
};

#include "Set.tpp"

#endif
