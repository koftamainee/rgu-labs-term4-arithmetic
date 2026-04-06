#include "algebra/polynomial/ideals.hpp"

#include <imgui.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>
#include <GLFW/glfw3.h>

#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <cmath>

struct DivStep {
  int m, n;
  int gen_idx;
  bool is_remainder;
};

struct DivAnim {
  std::vector<DivStep> steps;
  int current = -1;
  double timer = 0.0;
  double step_duration = 0.5;
  bool playing = false;
  std::vector<std::pair<int, int>> trail;
};

static DivAnim div_anim;

static std::vector<DivStep> compute_div_steps(
  const std::vector<Gen>& gens, int m0, int n0) {
  std::vector<DivStep> steps;
  int m = m0, n = n0;
  for (int iter = 0; iter < 64; ++iter) {
    bool divided = false;
    for (int i = 0; i < static_cast<int>(gens.size()); ++i) {
      if (m >= gens[i].a && n >= gens[i].b) {
        steps.push_back({m, n, i, false});
        m -= gens[i].a;
        n -= gens[i].b;
        divided = true;
        break;
      }
    }
    if (!divided) {
      bool zero = (m == 0 && n == 0);
      steps.push_back({m, n, -1, !zero});
      break;
    }
  }
  return steps;
}

static void draw_grid(
  ImDrawList* dl,
  ImVec2 origin,
  float cell,
  int N,
  const std::vector<Gen>& gens,
  bool show_remainder,
  int hover_m, int hover_n,
  const DivAnim& anim) {
  const ImU32 col_ideal = IM_COL32(100, 160, 240, 200);
  const ImU32 col_remainder = IM_COL32(80, 200, 140, 200);
  const ImU32 col_neutral = IM_COL32(235, 232, 225, 255);
  const ImU32 col_grid = IM_COL32(180, 180, 180, 80);
  const ImU32 col_hover = IM_COL32(255, 220, 80, 180);
  const ImU32 col_gen = IM_COL32(220, 60, 60, 255);
  const ImU32 col_gen_ring = IM_COL32(255, 255, 255, 200);
  const ImU32 col_label = IM_COL32(100, 100, 100, 200);
  const ImU32 col_axis = IM_COL32(80, 80, 80, 180);
  const ImU32 col_trail = IM_COL32(255, 200, 60, 80);

  for (int m = 0; m <= N; ++m) {
    for (int n = 0; n <= N; ++n) {
      float x0 = origin.x + static_cast<float>(m) * cell;
      float y0 = origin.y - static_cast<float>(n + 1) * cell;
      ImVec2 p0(x0, y0), p1(x0 + cell, y0 + cell);

      bool ideal = in_ideal(gens, m, n);
      ImU32 fill;
      if (m == hover_m && n == hover_n)
        fill = col_hover;
      else if (show_remainder)
        fill = ideal ? col_neutral : col_remainder;
      else
        fill = ideal ? col_ideal : col_neutral;

      dl->AddRectFilled(p0, p1, fill);
      dl->AddRect(p0, p1, col_grid, 0.0f, 0, 0.5f);
    }
  }

  for (auto& [tm, tn] : anim.trail) {
    if (tm > N || tn > N) continue;
    float x0 = origin.x + static_cast<float>(tm) * cell;
    float y0 = origin.y - static_cast<float>(tn + 1) * cell;
    dl->AddRectFilled({x0 + 1, y0 + 1}, {x0 + cell - 1, y0 + cell - 1}, col_trail);
  }

  if (anim.current >= 0 && anim.current < static_cast<int>(anim.steps.size())) {
    const auto& step = anim.steps[anim.current];
    if (step.m <= N && step.n <= N) {
      float x0 = origin.x + static_cast<float>(step.m) * cell;
      float y0 = origin.y - static_cast<float>(step.n + 1) * cell;

      ImU32 ring_col = step.is_remainder
                         ? IM_COL32(80, 220, 140, 255)
                         : IM_COL32(255, 220, 60, 255);
      dl->AddRect({x0, y0}, {x0 + cell, y0 + cell}, ring_col, 0.f, 0, 2.5f);

      float cx = x0 + cell * 0.5f;
      float cy = y0 + cell * 0.5f;
      auto t = static_cast<float>(fmod(ImGui::GetTime() * 4.0, 1.0));
      float r = cell * 0.3f + cell * 0.15f * t;
      ImU32 pulse_col = step.is_remainder
                          ? IM_COL32(80, 220, 140, static_cast<int>(180 * (1.f - t)))
                          : IM_COL32(255, 220, 60, static_cast<int>(180 * (1.f - t)));
      dl->AddCircle({cx, cy}, r, pulse_col, 0, 2.0f);

      char lbl[32];
      snprintf(lbl, sizeof(lbl), "x^%d y^%d", step.m, step.n);
      dl->AddText({x0, y0 - 14.f}, IM_COL32(255, 255, 255, 220), lbl);

      if (anim.current + 1 < static_cast<int>(anim.steps.size())) {
        const auto& next = anim.steps[anim.current + 1];
        if (next.m <= N && next.n <= N) {
          float nx = origin.x + static_cast<float>(next.m) * cell + cell * 0.5f;
          float ny = origin.y - static_cast<float>(next.n + 1) * cell + cell * 0.5f;
          ImU32 arrow_col = IM_COL32(255, 200, 60, 160);
          dl->AddLine({cx, cy}, {nx, ny}, arrow_col, 1.5f);
          ImVec2 d = {nx - cx, ny - cy};
          float len = sqrtf(d.x * d.x + d.y * d.y);
          if (len > 0) {
            d.x /= len;
            d.y /= len;
            ImVec2 p = {nx - d.x * 8.f, ny - d.y * 8.f};
            ImVec2 perp = {-d.y * 4.f, d.x * 4.f};
            dl->AddTriangleFilled(
              {nx, ny},
              {p.x + perp.x, p.y + perp.y},
              {p.x - perp.x, p.y - perp.y},
              arrow_col);
          }
        }
      }
    }
  }

  for (int m = 0; m <= N; m += 2) {
    float x = origin.x + static_cast<float>(m) * cell + cell * 0.5f;
    float y = origin.y + 4.0f;
    char buf[8];
    snprintf(buf, sizeof(buf), "%d", m);
    dl->AddText(ImVec2(x - 4, y), col_label, buf);
  }
  for (int n = 0; n <= N; n += 2) {
    float x = origin.x - 20.0f;
    float y = origin.y - static_cast<float>(n) * cell - cell * 0.5f - 7.0f;
    char buf[8];
    snprintf(buf, sizeof(buf), "%d", n);
    dl->AddText(ImVec2(x, y), col_label, buf);
  }

  ImVec2 ax0(origin.x, origin.y);
  ImVec2 ax1(origin.x + static_cast<float>(N + 1) * cell + 10, origin.y);
  dl->AddLine(ax0, ax1, col_axis, 1.5f);
  dl->AddTriangleFilled(
    ImVec2(ax1.x + 8, ax1.y),
    ImVec2(ax1.x, ax1.y - 4),
    ImVec2(ax1.x, ax1.y + 4), col_axis);
  dl->AddText(ImVec2(ax1.x + 10, ax1.y - 7), col_axis, "m");

  ImVec2 ay0(origin.x, origin.y);
  ImVec2 ay1(origin.x, origin.y - static_cast<float>(N + 1) * cell - 10);
  dl->AddLine(ay0, ay1, col_axis, 1.5f);
  dl->AddTriangleFilled(
    ImVec2(ay1.x, ay1.y - 8),
    ImVec2(ay1.x - 4, ay1.y),
    ImVec2(ay1.x + 4, ay1.y), col_axis);
  dl->AddText(ImVec2(ay1.x + 4, ay1.y - 4), col_axis, "n");

  for (const auto& g : gens) {
    if (g.a > N || g.b > N) continue;
    float cx = origin.x + static_cast<float>(g.a) * cell + cell * 0.5f;
    float cy = origin.y - static_cast<float>(g.b) * cell - cell * 0.5f;
    dl->AddCircleFilled(ImVec2(cx, cy), 7.0f, col_gen);
    dl->AddCircle(ImVec2(cx, cy), 9.0f, col_gen_ring, 0, 1.5f);
  }
}

static void draw_division_panel(
  const std::vector<Gen>& gens,
  const RingPtr& ring,
  int N,
  DivAnim& anim) {
  ImGui::SeparatorText("Division demo");
  ImGui::TextDisabled("Divide x^m * y^n by the generators using Lex order");

  static int dm = 4, dn = 4;
  ImGui::SetNextItemWidth(80);
  ImGui::InputInt("m##div", &dm);
  dm = std::clamp(dm, 0, N);
  ImGui::SameLine();
  ImGui::SetNextItemWidth(80);
  ImGui::InputInt("n##div", &dn);
  dn = std::clamp(dn, 0, N);

  if (gens.empty()) {
    ImGui::TextDisabled("Add generators first.");
    return;
  }

  ImGui::Spacing();

  if (ImGui::Button("Animate")) {
    anim.steps = compute_div_steps(gens, dm, dn);
    anim.current = 0;
    anim.trail.clear();
    anim.playing = true;
    anim.timer = 0.0;
    if (!anim.steps.empty())
      anim.trail.emplace_back(anim.steps[0].m, anim.steps[0].n);
  }
  ImGui::SameLine();
  if (ImGui::Button(anim.playing ? "Pause" : "Play"))
    anim.playing = !anim.playing;
  ImGui::SameLine();
  if (ImGui::Button("Reset")) {
    anim.current = 0;
    anim.trail.clear();
    anim.playing = false;
  }

  ImGui::SetNextItemWidth(120);
  auto spd = static_cast<float>(anim.step_duration);
  if (ImGui::SliderFloat("speed##anim", &spd, 0.1f, 2.0f, "%.1fs"))
    anim.step_duration = spd;

  if (anim.current >= 0 && !anim.steps.empty()) {
    int s = anim.current;
    ImGui::SetNextItemWidth(150);
    if (ImGui::SliderInt("step##anim", &s, 0, static_cast<int>(anim.steps.size()) - 1)) {
      anim.current = s;
      anim.trail.clear();
      for (int i = 0; i <= s; ++i)
        anim.trail.emplace_back(anim.steps[i].m, anim.steps[i].n);
    }
  }

  ImGui::Spacing();

  auto res = compute_division(gens, ring, dm, dn);

  ImGui::Text("f = x^%d * y^%d", dm, dn);
  ImGui::Spacing();

  for (size_t i = 0; i < gens.size(); ++i) {
    ImGui::Text("q[%zu] (div by x^%d y^%d) = %s",
                i, gens[i].a, gens[i].b,
                res.quotient_strs[i].c_str());
  }

  ImGui::Spacing();
  if (!res.remainder_zero) {
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.3f, 0.8f, 0.5f, 1.0f));
    ImGui::Text("remainder = %s", res.remainder_str.c_str());
    ImGui::PopStyleColor();
  }
  else {
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.6f, 0.6f, 0.6f, 1.0f));
    ImGui::Text("remainder = 0  (monomial is in I)");
    ImGui::PopStyleColor();
  }

  if (anim.current >= 0 && !anim.steps.empty()) {
    ImGui::Spacing();
    ImGui::SeparatorText("Steps");
    for (int i = 0; i < static_cast<int>(anim.steps.size()); ++i) {
      const auto& st = anim.steps[i];
      bool is_current = (i == anim.current);
      if (is_current)
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.85f, 0.25f, 1.0f));
      if (st.gen_idx >= 0)
        ImGui::Text("%s[%d] x^%d y^%d  / x^%d y^%d",
                    is_current ? "> " : "  ",
                    i, st.m, st.n,
                    gens[st.gen_idx].a, gens[st.gen_idx].b);
      else if (st.is_remainder)
        ImGui::Text("%s[%d] remainder x^%d y^%d",
                    is_current ? "> " : "  ",
                    i, st.m, st.n);
      else
        ImGui::Text("%s[%d] done (= 0)",
                    is_current ? "> " : "  ", i);
      if (is_current)
        ImGui::PopStyleColor();
    }
  }
}

int main() {
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  GLFWwindow* window = glfwCreateWindow(
    1100, 750,
    "Ideal Visualizer -- Cox-Little-O'Shea ss4 Ex.3",
    nullptr, nullptr);
  glfwMakeContextCurrent(window);
  glfwSwapInterval(1);

  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGui::StyleColorsDark();
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init("#version 330");

  ImGuiStyle& style = ImGui::GetStyle();
  style.WindowRounding = 6.0f;
  style.FrameRounding = 4.0f;
  style.ItemSpacing = ImVec2(8, 6);

  std::vector<Gen> gens = {{6, 0}, {2, 3}, {1, 7}};
  int N = 14;
  bool show_remainder = false;
  static int new_a = 0, new_b = 0;

  auto ring = make_ring<double>({"x", "y"});

  float splitter_x = 310.0f;
  bool dragging_split = false;

  float pan_x = 0.0f, pan_y = 0.0f;
  float zoom = 1.0f;
  bool panning = false;
  ImVec2 pan_start = {0, 0};
  ImVec2 pan_origin = {0, 0};

  while (!glfwWindowShouldClose(window)) {
    const float splitter_w = 5.0f;
    glfwPollEvents();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    double dt = ImGui::GetIO().DeltaTime;
    if (div_anim.playing && div_anim.current >= 0) {
      div_anim.timer += dt;
      if (div_anim.timer >= div_anim.step_duration) {
        div_anim.timer = 0.0;
        if (div_anim.current + 1 < static_cast<int>(div_anim.steps.size())) {
          div_anim.current++;
          const auto& s = div_anim.steps[div_anim.current];
          div_anim.trail.emplace_back(s.m, s.n);
        }
        else {
          div_anim.playing = false;
        }
      }
    }

    int ww, wh;
    glfwGetFramebufferSize(window, &ww, &wh);
    auto fww = static_cast<float>(ww);
    auto fwh = static_cast<float>(wh);

    splitter_x = std::clamp(splitter_x, 150.0f, fww - 200.0f);

    ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(splitter_x, fwh), ImGuiCond_Always);
    ImGui::Begin("Controls", nullptr,
                 ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                 ImGuiWindowFlags_NoCollapse);

    ImGui::SeparatorText("Ideal I = < generators >");

    int remove_idx = -1;
    for (int i = 0; i < static_cast<int>(gens.size()); ++i) {
      ImGui::PushID(i);
      ImGui::BulletText("x^%d * y^%d", gens[i].a, gens[i].b);
      ImGui::SameLine();
      if (ImGui::SmallButton("x")) remove_idx = i;
      ImGui::PopID();
    }
    if (remove_idx >= 0)
      gens.erase(gens.begin() + remove_idx);

    ImGui::Spacing();
    ImGui::SetNextItemWidth(70);
    ImGui::InputInt("a##new", &new_a);
    new_a = std::max(0, new_a);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(70);
    ImGui::InputInt("b##new", &new_b);
    new_b = std::max(0, new_b);
    ImGui::SameLine();
    if (ImGui::Button("Add")) {
      bool dup = false;
      for (const auto& g : gens)
        if (g.a == new_a && g.b == new_b) {
          dup = true;
          break;
        }
      if (!dup) gens.push_back({new_a, new_b});
    }

    ImGui::Spacing();
    ImGui::SeparatorText("Grid size");
    ImGui::SetNextItemWidth(120);
    ImGui::SliderInt("N##grid", &N, 5, 20);

    ImGui::Spacing();
    ImGui::SeparatorText("View mode");
    if (ImGui::RadioButton("(a)  Monomials in I", !show_remainder))
      show_remainder = false;
    if (ImGui::RadioButton("(b)  Remainder monomials", show_remainder))
      show_remainder = true;

    ImGui::Spacing();
    ImGui::SeparatorText("Legend");
    if (!show_remainder) {
      ImGui::ColorButton("##c1", ImVec4(0.39f, 0.63f, 0.94f, 0.78f),
                         ImGuiColorEditFlags_NoTooltip, ImVec2(14, 14));
      ImGui::SameLine();
      ImGui::Text("in I  (x^m y^n divisible by some gen)");
      ImGui::ColorButton("##c2", ImVec4(0.92f, 0.91f, 0.88f, 1.0f),
                         ImGuiColorEditFlags_NoTooltip, ImVec2(14, 14));
      ImGui::SameLine();
      ImGui::Text("not in I");
    }
    else {
      ImGui::ColorButton("##c3", ImVec4(0.31f, 0.78f, 0.55f, 0.78f),
                         ImGuiColorEditFlags_NoTooltip, ImVec2(14, 14));
      ImGui::SameLine();
      ImGui::Text("remainder  (not in I)");
      ImGui::ColorButton("##c4", ImVec4(0.92f, 0.91f, 0.88f, 1.0f),
                         ImGuiColorEditFlags_NoTooltip, ImVec2(14, 14));
      ImGui::SameLine();
      ImGui::Text("in I");
    }
    ImGui::ColorButton("##c5", ImVec4(0.86f, 0.24f, 0.24f, 1.0f),
                       ImGuiColorEditFlags_NoTooltip, ImVec2(14, 14));
    ImGui::SameLine();
    ImGui::Text("generator");

    ImGui::Spacing();
    draw_division_panel(gens, ring, N, div_anim);

    ImGui::End();

    ImVec2 mouse = ImGui::GetMousePos();
    bool mouse_down = ImGui::IsMouseDown(ImGuiMouseButton_Left);

    float sx0 = splitter_x;
    float sx1 = splitter_x + splitter_w;
    bool over_split = (mouse.x >= sx0 && mouse.x <= sx1 &&
      mouse.y >= 0 && mouse.y <= fwh);

    if (over_split || dragging_split)
      ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeEW);

    if (over_split && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
      dragging_split = true;
    if (!mouse_down)
      dragging_split = false;
    if (dragging_split)
      splitter_x = mouse.x;

    ImDrawList* bg = ImGui::GetBackgroundDrawList();
    bg->AddRectFilled(
      ImVec2(sx0, 0), ImVec2(sx1, fwh),
      over_split || dragging_split
        ? IM_COL32(180, 180, 200, 180)
        : IM_COL32(100, 100, 110, 120));

    float canvas_x = splitter_x + splitter_w;
    float canvas_w = fww - canvas_x;

    ImGui::SetNextWindowPos(ImVec2(canvas_x, 0), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(canvas_w, fwh), ImGuiCond_Always);
    ImGui::Begin("Grid", nullptr,
                 ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                 ImGuiWindowFlags_NoCollapse |
                 ImGuiWindowFlags_NoScrollbar |
                 ImGuiWindowFlags_NoScrollWithMouse);

    if (ImGui::IsWindowHovered()) {
      float wheel = ImGui::GetIO().MouseWheel;
      if (wheel != 0.0f) {
        float old_zoom = zoom;
        zoom *= (wheel > 0 ? 1.12f : 1.0f / 1.12f);
        zoom = std::clamp(zoom, 0.2f, 8.0f);
        float mx = mouse.x - canvas_x;
        float my = mouse.y;
        pan_x = mx - (mx - pan_x) * (zoom / old_zoom);
        pan_y = my - (my - pan_y) * (zoom / old_zoom);
      }

      if (ImGui::IsMouseClicked(ImGuiMouseButton_Middle) ||
        (ImGui::IsMouseClicked(ImGuiMouseButton_Left) &&
          ImGui::GetIO().KeyAlt)) {
        panning = true;
        pan_start = mouse;
        pan_origin = {pan_x, pan_y};
      }
    }
    if (!ImGui::IsMouseDown(ImGuiMouseButton_Middle) &&
      !(ImGui::IsMouseDown(ImGuiMouseButton_Left) && ImGui::GetIO().KeyAlt))
      panning = false;

    if (panning) {
      pan_x = pan_origin.x + (mouse.x - pan_start.x);
      pan_y = pan_origin.y + (mouse.y - pan_start.y);
      ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);
    }

    float pad = 50.0f;
    float base_cell = std::min(
      (canvas_w - pad * 2) / static_cast<float>(N + 1),
      (fwh - pad * 2 - 30) / static_cast<float>(N + 1));
    base_cell = std::max(base_cell, 8.0f);
    float cell = base_cell * zoom;

    ImVec2 win_pos = ImGui::GetWindowPos();
    ImVec2 origin(
      win_pos.x + pad + pan_x,
      win_pos.y + 30 + pad + static_cast<float>(N + 1) * cell + pan_y);

    ImDrawList* dl = ImGui::GetWindowDrawList();

    int hover_m = -1, hover_n = -1;
    float gx0 = origin.x;
    float gy1 = origin.y;
    float gx1 = gx0 + static_cast<float>(N + 1) * cell;
    float gy0 = gy1 - static_cast<float>(N + 1) * cell;
    if (!panning &&
      mouse.x >= gx0 && mouse.x < gx1 &&
      mouse.y >= gy0 && mouse.y < gy1) {
      hover_m = static_cast<int>((mouse.x - gx0) / cell);
      hover_n = N - static_cast<int>((mouse.y - gy0) / cell);
      hover_n = std::clamp(hover_n, 0, N);
    }

    draw_grid(dl, origin, cell, N, gens, show_remainder, hover_m, hover_n, div_anim);

    if (hover_m >= 0) {
      bool ideal = in_ideal(gens, hover_m, hover_n);
      ImGui::SetTooltip("x^%d y^%d  ->  %s",
                        hover_m, hover_n,
                        ideal ? "in I" : "not in I (remainder)");
    }

    const char* title = show_remainder
                          ? "(b) Monomials that can appear in the remainder"
                          : "(a) Monomials x^m y^n that belong to I";
    ImGui::SetCursorPos(ImVec2(pad, 8));
    ImGui::TextDisabled("%s", title);

    ImGui::SetCursorPos(ImVec2(pad, 22));
    ImGui::TextDisabled("Scroll: zoom  |  Alt+drag or MMB: pan");

    ImGui::End();

    ImGui::Render();
    glViewport(0, 0, ww, wh);
    glClearColor(0.12f, 0.12f, 0.13f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    glfwSwapBuffers(window);
  }

  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();
  glfwDestroyWindow(window);
  glfwTerminate();
  return 0;
}
