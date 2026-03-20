#pragma once

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "latex_serializable.hpp"
#include "report_object.hpp"

class Reporter {
public:
  explicit Reporter(std::string title)
    : m_title(std::move(title)) {}

  Reporter& section(const std::string& text) {
    push(std::make_unique<HeadingBlock>(2, text));
    return *this;
  }

  Reporter& subsection(const std::string& text) {
    push(std::make_unique<HeadingBlock>(3, text));
    return *this;
  }

  Reporter& subsubsection(const std::string& text) {
    push(std::make_unique<HeadingBlock>(4, text));
    return *this;
  }

  Reporter& text(const std::string& paragraph) {
    push(std::make_unique<TextBlock>(paragraph));
    return *this;
  }

  template <typename... Args>
  Reporter& text(const std::string& fmt, const Args&... args) {
    return text(format(fmt, args...));
  }

  Reporter& raw_latex(const std::string& latex) {
    push(std::make_unique<MathBlock>(latex));
    return *this;
  }

  Reporter& math(const ILatexSerializable& obj) {
    push(std::make_unique<MathBlock>(obj.to_latex()));
    return *this;
  }

  Reporter& add(const IReportObject& obj) {
    push(std::make_unique<RawBlock>(obj.to_html()));
    return *this;
  }

  std::string build() const {
    std::string dir = make_artifact_dir();
    if (dir == "/tmp/report_") {
      throw std::invalid_argument("Reporter::build: title produces empty directory name");
    }
    std::string assets_dir = dir + "/assets";
    std::filesystem::create_directories(dir);
    std::filesystem::create_directories(assets_dir);

    for (auto& block : m_blocks) {
      block->prepare(assets_dir);
    }

    std::string path = dir + "/index.html";
    std::ofstream out(path);
    if (!out.is_open()) {
      throw std::runtime_error("Reporter::build: cannot open file: " + path);
    }
    out << render_shell();
    return dir;
  }

private:
  struct HeadingBlock : IReportObject {
    int level;
    std::string text;

    HeadingBlock(int level, std::string text)
      : level(level), text(std::move(text)) {}

    std::string to_html() const override {
      std::string tag = "h" + std::to_string(level);
      return "<" + tag + ">" + escape_html(text) + "</" + tag + ">\n";
    }
  };

  struct TextBlock : IReportObject {
    std::string content;

    explicit TextBlock(std::string content)
      : content(std::move(content)) {}

    std::string to_html() const override {
      return "<p>" + content + "</p>\n";
    }
  };

  struct MathBlock : IReportObject {
    std::string latex;

    explicit MathBlock(std::string latex)
      : latex(std::move(latex)) {}

    std::string to_html() const override {
      return R"(<p class="math">\[)" + latex + "\\]</p>\n";
    }
  };

  struct RawBlock : IReportObject {
    std::string html;

    explicit RawBlock(std::string html)
      : html(std::move(html)) {}

    std::string to_html() const override {
      return html + "\n";
    }
  };

  static std::string escape_html(const std::string& s) {
    std::string result;
    result.reserve(s.size());
    for (const char c : s) {
      switch (c) {
      case '&': result += "&amp;";
        break;
      case '<': result += "&lt;";
        break;
      case '>': result += "&gt;";
        break;
      case '"': result += "&quot;";
        break;
      case '\'': result += "&#39;";
        break;
      default: result += c;
        break;
      }
    }
    return result;
  }

  static std::string format(const std::string& fmt) {
    if (fmt.find("{}") != std::string::npos) {
      throw std::invalid_argument("Reporter::format: too few arguments for format string");
    }
    return fmt;
  }

  template <typename T, typename... Rest>
  static std::string format(const std::string& fmt, const T& first, const Rest&... rest) {
    auto pos = fmt.find("{}");
    if (pos == std::string::npos) {
      throw std::invalid_argument("Reporter::format: too many arguments for format string");
    }
    std::string result = fmt.substr(0, pos);
    if constexpr (std::is_base_of_v<ILatexSerializable, T>) {
      result += "\\(" + first.to_latex() + "\\)";
    }
    else {
      std::ostringstream oss;
      oss << first;
      result += oss.str();
    }
    result += format(fmt.substr(pos + 2), rest...);
    return result;
  }

  void push(std::unique_ptr<IReportObject> block) {
    m_blocks.push_back(std::move(block));
  }

  static std::string sanitize_title(const std::string& title) {
    std::string result;
    for (const char c : title) {
      if (std::isalnum(static_cast<unsigned char>(c))) {
        result += c;
      }
      else if (std::isspace(static_cast<unsigned char>(c))) {
        result += '_';
      }
    }
    return result;
  }

  std::string make_artifact_dir() const {
    return "/tmp/report_" + sanitize_title(m_title);
  }

  static std::string current_timestamp() {
    auto now = std::chrono::system_clock::now();
    auto t = std::chrono::system_clock::to_time_t(now);
    std::ostringstream oss;
    oss << std::put_time(std::gmtime(&t), "%Y-%m-%d %H:%M:%S UTC");
    return oss.str();
  }

  std::string render_shell() const {
    std::ostringstream html;
    html << R"(<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>)" << escape_html(m_title) << R"(</title>
</head>
<body>
<header>
  <h1>)" << escape_html(m_title) << R"(</h1>
  <div class="timestamp">)" << current_timestamp() << R"(</div>
</header>
<main>
)";
    for (const auto& block : m_blocks) {
      html << block->to_html();
    }
    html << R"(</main>
</body>
</html>
)";
    return html.str();
  }

  std::string m_title;
  std::vector<std::unique_ptr<IReportObject>> m_blocks;
};
