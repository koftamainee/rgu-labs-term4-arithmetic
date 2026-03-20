#pragma once

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "report_object.hpp"
#include "latex_serializable.hpp"

class LatexList : public IReportObject {
public:
    enum class Kind { Itemize, Enumerate };

    static LatexList itemize() { return LatexList(Kind::Itemize); }
    static LatexList enumerate() { return LatexList(Kind::Enumerate); }

    LatexList& item(const std::string& text) {
        m_items.push_back(text);
        return *this;
    }

    template <typename... Args>
    LatexList& item(const std::string& fmt, const Args&... args) {
        m_items.push_back(format(fmt, args...));
        return *this;
    }

    std::string to_html() const override {
        std::string tag = (m_kind == Kind::Itemize) ? "ul" : "ol";
        std::string result = "<" + tag + ">\n";
        for (const auto& item : m_items) {
            result += "<li>" + item + "</li>\n";
        }
        result += "</" + tag + ">\n";
        return result;
    }

private:
    explicit LatexList(Kind kind) : m_kind(kind) {}

    static std::string format(const std::string& fmt) {
        if (fmt.find("{}") != std::string::npos) {
            throw std::invalid_argument("LatexList::format: too few arguments for format string");
        }
        return fmt;
    }

    template <typename T, typename... Rest>
    static std::string format(const std::string& fmt, const T& first, const Rest&... rest) {
        auto pos = fmt.find("{}");
        if (pos == std::string::npos) {
            throw std::invalid_argument("LatexList::format: too many arguments for format string");
        }
        std::string result = fmt.substr(0, pos);
        if constexpr (std::is_base_of_v<ILatexSerializable, T>) {
            result += "\\(" + first.to_latex() + "\\)";
        } else {
            std::ostringstream oss;
            oss << first;
            result += oss.str();
        }
        result += format(fmt.substr(pos + 2), rest...);
        return result;
    }

    Kind m_kind;
    std::vector<std::string> m_items;
};