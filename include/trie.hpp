#pragma once

#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

template <typename Key, typename V>
class Trie final {
public:
  using key_type = Key;
  using symbol_type = typename Key::value_type;
  using value_type = V;
  using size_type = std::size_t;

  struct Node {
    std::map<symbol_type, std::unique_ptr<Node>> children;
    std::optional<value_type> value;
  };

private:
  template <bool IsConst>
  class basic_iterator {
    friend class Trie;

    using node_ptr = std::conditional_t<IsConst, const Node*, Node*>;
    using value_ref = std::conditional_t<IsConst, const value_type&, value_type&>;
    using map_iterator = std::conditional_t<
      IsConst,
      typename std::map<symbol_type, std::unique_ptr<Node>>::const_iterator,
      typename std::map<symbol_type, std::unique_ptr<Node>>::iterator
    >;

    struct Frame {
      node_ptr node;
      map_iterator it;
    };

    std::vector<Frame> m_stack;
    key_type m_current_key;

    void advance() {
      while (!m_stack.empty()) {
        auto& [node, it] = m_stack.back();
        if (it != node->children.end()) {
          node_ptr child = it->second.get();
          m_current_key.push_back(it->first);
          ++it;
          m_stack.push_back({child, child->children.begin()});
          if (child->value.has_value()) {
            return;
          }
        }
        else {
          m_stack.pop_back();
          if (!m_current_key.empty()) {
            m_current_key.pop_back();
          }
        }
      }
    }

    explicit basic_iterator(node_ptr root) {
      if (root != nullptr) {
        m_stack.reserve(8);
        m_stack.push_back({root, root->children.begin()});
        if (!root->value.has_value()) {
          advance();
        }
      }
    }

    basic_iterator(node_ptr root, node_ptr target, key_type key)
      : m_current_key(std::move(key)) {
      if (root != nullptr && target != nullptr) {
        m_stack.reserve(m_current_key.size() + 1);
        node_ptr current = root;
        m_stack.push_back({root, root->children.begin()});
        for (const symbol_type& symbol : m_current_key) {
          auto it = current->children.find(symbol);
          m_stack.back().it = std::next(it);
          current = it->second.get();
          m_stack.push_back({current, current->children.begin()});
        }
      }
    }

    basic_iterator() = default;

  public:
    struct entry {
      const key_type& key;
      value_ref value;
    };

    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;

    const key_type& key() const {
      return m_current_key;
    }

    value_ref value() const {
      return m_stack.back().node->value.value();
    }

    entry operator*() const {
      return {m_current_key, m_stack.back().node->value.value()};
    }

    basic_iterator& operator++() {
      advance();
      return *this;
    }

    basic_iterator operator++(int) {
      basic_iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    bool operator==(const basic_iterator& other) const {
      if (m_stack.empty() && other.m_stack.empty()) {
        return true;
      }
      if (m_stack.empty() || other.m_stack.empty()) {
        return false;
      }
      return m_stack.back().node == other.m_stack.back().node;
    }

    bool operator!=(const basic_iterator& other) const {
      return !(*this == other);
    }

    explicit operator basic_iterator<true>() const {
      basic_iterator<true> result;
      result.m_current_key = m_current_key;
      result.m_stack.reserve(m_stack.size());
      for (const auto& frame : m_stack) {
        result.m_stack.push_back({frame.node, frame.it});
      }
      return result;
    }
  };

public:
  using iterator = basic_iterator<false>;
  using const_iterator = basic_iterator<true>;

  Trie() = default;
  Trie(const Trie& other);
  Trie(Trie&&) noexcept = default;
  ~Trie() = default;

  Trie& operator=(const Trie& other);
  Trie& operator=(Trie&&) noexcept = default;

  [[nodiscard]] bool empty() const noexcept;
  size_type size() const noexcept;

  iterator begin() noexcept;
  iterator end() noexcept;
  const_iterator begin() const noexcept;
  const_iterator end() const noexcept;
  const_iterator cbegin() const noexcept;
  const_iterator cend() const noexcept;

  iterator insert(const key_type& key, const value_type& value);
  iterator insert(const key_type& key, value_type&& value);
  iterator erase(const key_type& key);
  void clear() noexcept;

  value_type& at(const key_type& key);
  const value_type& at(const key_type& key) const;
  value_type& operator[](const key_type& key);

  iterator find(const key_type& key);
  const_iterator find(const key_type& key) const;
  bool contains(const key_type& key) const;

  const Node* root() const noexcept;

  bool operator==(const Trie&) const;
  bool operator!=(const Trie&) const;

private:
  template <typename W>
  iterator insert_impl(const key_type& key, W&& value, bool overwrite);

  Node* find_node_ptr(const key_type& key);
  const Node* find_node_ptr(const key_type& key) const;
  bool erase_impl(Node* node, const key_type& key, size_type depth);
  static std::unique_ptr<Node> clone_impl(const Node* source);

  std::unique_ptr<Node> m_root = nullptr;
  size_type m_size = 0;
};

template <typename Key, typename V>
Trie<Key, V>::Trie(const Trie& other) : m_size(other.m_size) {
  m_root = clone_impl(other.m_root.get());
}

template <typename Key, typename V>
Trie<Key, V>& Trie<Key, V>::operator=(const Trie& other) {
  if (this != &other) {
    m_root = clone_impl(other.m_root.get());
    m_size = other.m_size;
  }
  return *this;
}

template <typename Key, typename V>
bool Trie<Key, V>::empty() const noexcept {
  return m_size == 0;
}

template <typename Key, typename V>
typename Trie<Key, V>::size_type Trie<Key, V>::size() const noexcept {
  return m_size;
}

template <typename Key, typename V>
typename Trie<Key, V>::iterator Trie<Key, V>::begin() noexcept {
  return iterator(m_root.get());
}

template <typename Key, typename V>
typename Trie<Key, V>::iterator Trie<Key, V>::end() noexcept {
  return iterator();
}

template <typename Key, typename V>
typename Trie<Key, V>::const_iterator Trie<Key, V>::begin() const noexcept {
  return const_iterator(m_root.get());
}

template <typename Key, typename V>
typename Trie<Key, V>::const_iterator Trie<Key, V>::end() const noexcept {
  return const_iterator();
}

template <typename Key, typename V>
typename Trie<Key, V>::const_iterator Trie<Key, V>::cbegin() const noexcept {
  return const_iterator(m_root.get());
}

template <typename Key, typename V>
typename Trie<Key, V>::const_iterator Trie<Key, V>::cend() const noexcept {
  return const_iterator();
}

template <typename Key, typename V>
typename Trie<Key, V>::iterator Trie<Key, V>::insert(const key_type& key, const value_type& value) {
  return insert_impl(key, value, true);
}

template <typename Key, typename V>
typename Trie<Key, V>::iterator Trie<Key, V>::insert(const key_type& key, value_type&& value) {
  return insert_impl(key, std::move(value), true);
}

template <typename Key, typename V>
typename Trie<Key, V>::iterator Trie<Key, V>::erase(const key_type& key) {
  iterator next = find(key);
  if (next == end()) {
    return end();
  }
  ++next;
  erase_impl(m_root.get(), key, 0);
  return next;
}

template <typename Key, typename V>
void Trie<Key, V>::clear() noexcept {
  m_root.reset();
  m_size = 0;
}

template <typename Key, typename V>
typename Trie<Key, V>::value_type& Trie<Key, V>::at(const key_type& key) {
  Node* node = find_node_ptr(key);
  if (node == nullptr) {
    throw std::out_of_range("Trie::at: key not found");
  }
  return node->value.value();
}

template <typename Key, typename V>
const typename Trie<Key, V>::value_type& Trie<Key, V>::at(const key_type& key) const {
  const Node* node = find_node_ptr(key);
  if (node == nullptr) {
    throw std::out_of_range("Trie::at: key not found");
  }
  return node->value.value();
}

template <typename Key, typename V>
typename Trie<Key, V>::value_type& Trie<Key, V>::operator[](const key_type& key) {
  return insert_impl(key, value_type{}, false).value();
}

template <typename Key, typename V>
typename Trie<Key, V>::iterator Trie<Key, V>::find(const key_type& key) {
  Node* node = find_node_ptr(key);
  if (node == nullptr) {
    return end();
  }
  return iterator(m_root.get(), node, key);
}

template <typename Key, typename V>
typename Trie<Key, V>::const_iterator Trie<Key, V>::find(const key_type& key) const {
  const Node* node = find_node_ptr(key);
  if (node == nullptr) {
    return cend();
  }
  return const_iterator(m_root.get(), node, key);
}

template <typename Key, typename V>
bool Trie<Key, V>::contains(const key_type& key) const {
  return find_node_ptr(key) != nullptr;
}

template <typename Key, typename V>
const typename Trie<Key, V>::Node* Trie<Key, V>::root() const noexcept {
  return m_root.get();
}

template <typename Key, typename V>
bool Trie<Key, V>::operator==(const Trie& other) const {
  if (m_size != other.m_size) return false;
  auto it2 = other.cbegin();
  for (const auto& elem : *this) {
    if (elem.key != it2.key()) return false;
    if (elem.value != it2.value()) return false;
    ++it2;
  }
  return true;
}

template <typename Key, typename V>
bool Trie<Key, V>::operator!=(const Trie& other) const {
  return !(*this == other);
}

template <typename Key, typename V>
template <typename W>
typename Trie<Key, V>::iterator Trie<Key, V>::insert_impl(const key_type& key, W&& value, bool overwrite) {
  if (m_root == nullptr) {
    m_root = std::make_unique<Node>();
  }
  Node* current = m_root.get();
  for (const symbol_type& symbol : key) {
    auto it = current->children.find(symbol);
    if (it == current->children.end()) {
      auto [inserted, _] = current->children.emplace(symbol, std::make_unique<Node>());
      current = inserted->second.get();
    }
    else {
      current = it->second.get();
    }
  }
  if (!current->value.has_value()) {
    current->value = std::forward<W>(value);
    ++m_size;
  }
  else if (overwrite) {
    current->value = std::forward<W>(value);
  }
  return iterator(m_root.get(), current, key);
}

template <typename Key, typename V>
typename Trie<Key, V>::Node* Trie<Key, V>::find_node_ptr(const key_type& key) {
  if (m_root == nullptr) {
    return nullptr;
  }
  Node* current = m_root.get();
  for (const symbol_type& symbol : key) {
    auto it = current->children.find(symbol);
    if (it == current->children.end()) {
      return nullptr;
    }
    current = it->second.get();
  }
  if (!current->value.has_value()) {
    return nullptr;
  }
  return current;
}

template <typename Key, typename V>
const typename Trie<Key, V>::Node* Trie<Key, V>::find_node_ptr(const key_type& key) const {
  if (m_root == nullptr) {
    return nullptr;
  }
  const Node* current = m_root.get();
  for (const symbol_type& symbol : key) {
    auto it = current->children.find(symbol);
    if (it == current->children.end()) {
      return nullptr;
    }
    current = it->second.get();
  }
  if (!current->value.has_value()) {
    return nullptr;
  }
  return current;
}

template <typename Key, typename V>
bool Trie<Key, V>::erase_impl(Node* node, const key_type& key, size_type depth) {
  if (node == nullptr) {
    return false;
  }
  if (depth == key.size()) {
    if (!node->value.has_value()) {
      return false;
    }
    node->value.reset();
    --m_size;
    return node->children.empty();
  }
  auto it = node->children.find(key[depth]);
  if (it == node->children.end()) {
    return false;
  }
  bool should_delete = erase_impl(it->second.get(), key, depth + 1);
  if (should_delete) {
    node->children.erase(it);
  }
  return node->children.empty() && !node->value.has_value();
}

template <typename Key, typename V>
std::unique_ptr<typename Trie<Key, V>::Node> Trie<Key, V>::clone_impl(const Node* source) {
  if (source == nullptr) {
    return nullptr;
  }

  struct Frame {
    const Node* source;
    Node* dest;
    typename std::map<symbol_type, std::unique_ptr<Node>>::const_iterator it;
  };

  auto new_root = std::make_unique<Node>();
  new_root->value = source->value;

  std::vector<Frame> stack;
  stack.reserve(8);
  stack.push_back({source, new_root.get(), source->children.begin()});

  while (!stack.empty()) {
    auto& [src, dst, it] = stack.back();
    if (it != src->children.end()) {
      auto child = std::make_unique<Node>();
      child->value = it->second->value;
      Node* child_ptr = child.get();
      dst->children.emplace(it->first, std::move(child));
      const Node* src_child = it->second.get();
      ++it;
      stack.push_back({src_child, child_ptr, src_child->children.begin()});
    }
    else {
      stack.pop_back();
    }
  }

  return new_root;
}
