#ifndef RGU_LABS_TERM4_ARITHMETIC_POLY_RING_HPP
#define RGU_LABS_TERM4_ARITHMETIC_POLY_RING_HPP

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

template <typename T>
class TPolyRing {
public:
  using coeff_type = T;
  using size_type  = std::size_t;

  explicit TPolyRing(std::vector<std::string> var_names)
      : m_var_names(std::move(var_names)) {
    if (m_var_names.empty()) {
      throw std::invalid_argument("PolyRing: var_names must not be empty");
    }
    for (const auto& name : m_var_names) {
      if (name.empty()) {
        throw std::invalid_argument("PolyRing: variable name must not be empty");
      }
    }
  }

  size_type n_vars() const noexcept { return m_var_names.size(); }

  const std::string& var_name(size_type i) const { return m_var_names.at(i); }

  const std::vector<std::string>& var_names() const noexcept { return m_var_names; }

  bool operator==(const TPolyRing& other) const noexcept {
    return m_var_names == other.m_var_names;
  }

  bool operator!=(const TPolyRing& other) const noexcept {
    return !(*this == other);
  }

private:
  std::vector<std::string> m_var_names;
};

template <typename T>
std::shared_ptr<const TPolyRing<T>> make_ring(std::vector<std::string> var_names) {
  return std::make_shared<const TPolyRing<T>>(std::move(var_names));
}

#endif //RGU_LABS_TERM4_ARITHMETIC_POLY_RING_HPP
