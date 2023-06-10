#ifndef PHASE_SET_H
#define PHASE_SET_H

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

  Set(const std::string &name = "") : _name(name) {}

  const std::string &name() const { return _name; }

  virtual bool add(const T &item);

  virtual bool remove(const T &item);

  //- Iterators

  const_iterator begin() const { return _items.begin(); }

  const_iterator end() const { return _items.end(); }

  //- Access
  const T &operator[](int i) const { return _items[i]; }

  const T &operator[](Label i) const { return _items[i]; }

protected:
  std::string _name;

  std::unordered_set<Ref<const T>, Hash, EqualTo> _itemSet;

  std::vector<Ref<const T>> _items;
};

#include "Set.tpp"

#endif
