#pragma once

#include <string>

class IReportObject {
public:
  virtual std::string to_html() const = 0;
  virtual void prepare(const std::string& assets_dir) { (void)assets_dir; }
  virtual ~IReportObject() = default;
};
