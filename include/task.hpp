#pragma once
#include <memory>
#include <string>

#include "publisher.hpp"
#include "reporter.hpp"

class Task final {
public:
  explicit Task(const std::string& task_name)
    : m_reporter(std::make_unique<Reporter>(task_name)),
      m_publisher(std::make_unique<Publisher>("localhost", 8080)) {}

  void build_and_publish(const std::function<void(Reporter&)>& callback) const {
    try {
      callback(*m_reporter);
      std::string artifact_dir = m_reporter->build();
      std::cout << "Report built at: " << artifact_dir << "\n";

      Publisher pub("localhost", 8080);
      pub.publish(artifact_dir);
      std::cout << "Published to http://localhost:8080\n";
    }
    catch (const std::exception& e) {
      std::cerr << "error: " << e.what() << "\n";
    }
  }

private:
  std::unique_ptr<Reporter> m_reporter;
  std::unique_ptr<Publisher> m_publisher;
};
