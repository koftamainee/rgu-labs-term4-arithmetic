#pragma once

#include <string>

class ILatexSerializable {
public:
  virtual std::string to_latex() const = 0;
  virtual ~ILatexSerializable() = default;
};
