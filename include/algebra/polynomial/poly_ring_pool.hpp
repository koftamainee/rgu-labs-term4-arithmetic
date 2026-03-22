#pragma once
#include <map>

#include "poly_ring.hpp"

struct NoLock final {
  void lock() noexcept {}
  void unlock() noexcept {}
};


template <typename T, typename Mutex = NoLock>
class PolyRingPool final {
public:
  using size_type = std::size_t;
  using ring_type = PolyRing<T>;
  using ring_ptr = std::shared_ptr<const ring_type>;
  using weak_ptr = std::weak_ptr<const ring_type>;
  using key_type = std::vector<std::string>;

  ring_ptr get(std::vector<std::string> var_names) {
    std::lock_guard<Mutex> guard(m_mutex);

    auto it = m_pool.find(var_names);
    if (it != m_pool.end()) {
      if (auto sp = it->second.lock()) {
        return sp;
      }
    }

    auto ring = std::make_shared<const ring_type>(var_names);
    m_pool[std::move(var_names)] = ring;
    return ring;
  }

  std::size_t size() const {
    std::lock_guard<Mutex> guard(m_mutex);
    return m_pool.size();
  }

  void purge_expired() {
    std::lock_guard<Mutex> guard(m_mutex);
    for (auto it = m_pool.begin(); it != m_pool.end();) {
      it = it->second.expired() ? m_pool.erase(it) : std::next(it);
    }
  }

private:
  mutable Mutex m_mutex;
  std::map<key_type, weak_ptr> m_pool;
};
