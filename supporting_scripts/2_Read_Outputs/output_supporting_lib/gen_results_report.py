# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of the FLUXOS model.

# This program, FLUXOS, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
HTML results-report generator for FLUXOS flood simulations.

Dark-theme, self-contained HTML file with Plotly charts (CDN-loaded) and the
same visual aesthetic as the 1_Model_Config report:
- Inter + JetBrains Mono fonts
- #0066cc / #00a86b / #ff6b35 palette
- Per-line copy buttons always visible
- Sidebar navigation + theme toggle

Sections (complementary to the WebGL dynamic view):
  1. Summary KPIs
  2. Time-series (volume / flooded area / max depth / mean wet depth / max velocity)
  3. Maximum inundation map (per-cell max depth, always the key output)
  4. Flood-hazard classification map (ARR-2019 H×V criterion)
  5. Depth histogram
  6. Time-of-first-inundation map (flood-front propagation)
  7. Chemical concentration panel (only if ADE_TRANSPORT output is present)
"""

from __future__ import annotations

import html
import json
import os
from datetime import datetime


# ---------------------------------------------------------------------------
#  CSS (matches 1_Model_Config gen_report.py aesthetic)
# ---------------------------------------------------------------------------

_CSS = r"""
:root,[data-theme="light"]{
  --primary:#0066cc;--primary-dark:#004499;
  --secondary:#00a86b;--accent:#ff6b35;
  --dark:#1a1a2e;--light:#f8f9fa;
  --bg:#f0f2f5;--surface:#fff;--text:#2d3436;--text2:#636e72;--text3:#a0aec0;
  --border:#e2e8f0;--border2:#d0d4da;
  --shadow:0 4px 20px rgba(0,0,0,.06);--shadow-lg:0 10px 40px rgba(0,0,0,.08);
  --code-bg:#1e1f2e;--code-text:#e2e8f0;
  --plot-bg:#fff;--plot-grid:#f0f0f0;--plot-font:#333;
  --glass:rgba(255,255,255,.85);
}
[data-theme="dark"]{
  --primary:#4d9ee8;--primary-dark:#3a8bd4;
  --secondary:#34d399;--accent:#fb923c;
  --dark:#1a1a2e;--light:#1e1f2e;
  --bg:#0f1117;--surface:#1a1b2e;--text:#e2e8f0;--text2:#a0aec0;--text3:#636e72;
  --border:#2a2b3d;--border2:#3a3b4d;
  --shadow:0 4px 20px rgba(0,0,0,.3);--shadow-lg:0 10px 40px rgba(0,0,0,.4);
  --code-bg:#0c0d15;--code-text:#e2e8f0;
  --plot-bg:#1a1b2e;--plot-grid:#2d2e42;--plot-font:#ccc;
  --glass:rgba(15,17,23,.92);
}
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:'Inter',-apple-system,BlinkMacSystemFont,sans-serif;background:var(--bg);
  color:var(--text);line-height:1.65;transition:background .3s,color .3s}
a{color:var(--primary);text-decoration:none}

.layout{display:flex;min-height:100vh}
.sidebar{width:220px;position:sticky;top:0;height:100vh;background:var(--glass);
  backdrop-filter:blur(12px);border-right:1px solid var(--border);padding:1.2rem .8rem;
  display:flex;flex-direction:column;z-index:10;overflow-y:auto}
.sidebar .logo{font-weight:700;font-size:1.1rem;color:var(--text);margin-bottom:1.2rem;
  padding-left:.7rem;letter-spacing:-.3px}
.sidebar .logo span{color:var(--primary)}
.sidebar nav{display:flex;flex-direction:column;gap:2px}
.sidebar nav a{display:block;padding:.4rem .7rem;border-radius:6px;font-size:.82rem;
  color:var(--text2);transition:all .15s}
.sidebar nav a:hover{background:rgba(77,158,232,.08);color:var(--primary)}
.sidebar nav a.active{background:rgba(77,158,232,.15);color:var(--primary);font-weight:600}
.theme-toggle{margin-top:auto;padding-top:1rem;border-top:1px solid var(--border)}
.theme-btn{display:flex;align-items:center;gap:.5rem;width:100%;padding:.45rem .7rem;
  border:1px solid var(--border);border-radius:8px;background:var(--surface);
  color:var(--text2);font-size:.78rem;cursor:pointer;transition:all .15s;font-family:inherit}
.theme-btn:hover{border-color:var(--primary);color:var(--primary)}

.main{flex:1;min-width:0;overflow-x:hidden}
.header{background:linear-gradient(135deg,#0f3460 0%,#16213e 50%,#1a1a2e 100%);
  color:#fff;padding:3rem 2.5rem 2rem;position:relative;overflow:hidden}
.header::before{content:'';position:absolute;inset:0;
  background:radial-gradient(circle at 80% 20%,rgba(255,255,255,.08) 0%,transparent 60%)}
.header h1{font-size:2rem;font-weight:700;letter-spacing:-.5px;position:relative}
.header h1 span{color:var(--primary)}
.header .subtitle{opacity:.85;font-size:.92rem;margin-top:.4rem;position:relative;max-width:800px}
.header .meta{display:flex;gap:.6rem;margin-top:1rem;font-size:.78rem;opacity:.8;
  flex-wrap:wrap;position:relative}
.header .meta span{background:rgba(255,255,255,.1);padding:.25rem .7rem;border-radius:20px}

.container{max-width:1180px;margin:0 auto;padding:2rem 2rem 3rem}
.section{margin-bottom:2.5rem;scroll-margin-top:1rem}
.section h2{font-size:1.2rem;font-weight:700;color:var(--text);margin-bottom:1rem;
  display:flex;align-items:center;gap:.5rem;letter-spacing:-.3px}
.section h2::before{content:'';width:5px;height:1.2em;border-radius:3px;
  background:linear-gradient(180deg,var(--primary),var(--primary-dark));flex-shrink:0}
.section h3{font-size:1rem;font-weight:600;margin:1rem 0 .6rem;color:var(--text)}

.card{background:var(--surface);border:1px solid var(--border);border-radius:16px;
  padding:1.5rem;box-shadow:var(--shadow);transition:transform .3s,box-shadow .3s}

.kpi-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:1rem;
  margin-bottom:1.5rem}
.kpi{background:var(--surface);border:1px solid var(--border);border-radius:16px;
  padding:1.2rem 1rem;text-align:center;box-shadow:var(--shadow);
  position:relative;overflow:hidden}
.kpi::before{content:'';position:absolute;top:0;left:0;right:0;height:3px;
  background:linear-gradient(90deg,var(--primary),var(--secondary))}
.kpi .icon{font-size:1.6rem;margin-bottom:.4rem}
.kpi .value{font-size:1.4rem;font-weight:700;color:var(--primary);
  font-family:'JetBrains Mono',monospace;overflow-wrap:break-word;word-break:break-word}
.kpi .label{font-size:.7rem;color:var(--text3);text-transform:uppercase;
  letter-spacing:.8px;margin-top:.3rem}
.kpi .sub{font-size:.7rem;color:var(--text3);margin-top:.2rem}

.table-wrap{overflow-x:auto;border-radius:12px;border:1px solid var(--border)}
table{width:100%;border-collapse:collapse;font-size:.88rem;background:var(--surface)}
th{background:var(--dark);color:#fff;padding:.65rem .9rem;text-align:left;font-weight:600;
  font-size:.72rem;text-transform:uppercase;letter-spacing:.5px}
[data-theme="dark"] th{background:#252640}
td{padding:.55rem .9rem;border-bottom:1px solid var(--border)}
tr:last-child td{border-bottom:none}
td.num{text-align:right;font-family:'JetBrains Mono',monospace;font-size:.82rem}

.plot{width:100%;min-height:380px;margin:1rem 0 0;
  border-radius:12px;overflow:hidden}
.plot-wide{min-height:480px}

/* Per-line copy buttons in code blocks (for the shell snippet sections) */
.code-block{position:relative;margin:.6rem 0}
.code-block pre{background:var(--code-bg);color:var(--code-text);padding:1rem;
  border-radius:10px;overflow-x:auto;font-family:'JetBrains Mono',monospace;
  font-size:.82rem;line-height:1.6}
.code-block pre .code-cmd{display:block;position:relative;border-radius:3px;padding-right:2.2rem}
.code-block pre .code-cmd:hover{background:rgba(255,255,255,.06)}
.code-block pre .code-blank,.code-block pre .code-comment-line{display:block}
.code-block pre .line-copy-btn{position:absolute;top:0;right:.25rem;padding:0 .45rem;
  font-size:.68rem;background:rgba(255,255,255,.1);border:1px solid rgba(255,255,255,.18);
  border-radius:4px;color:#e2e8f0;cursor:pointer;opacity:1;
  transition:background .15s,border-color .15s;font-family:inherit;line-height:1.55}
.code-block pre .line-copy-btn:hover{background:rgba(255,255,255,.22);border-color:rgba(255,255,255,.35)}
.code-block pre .line-copy-btn.copied{background:var(--secondary);border-color:var(--secondary);color:#fff}

.highlight-box{padding:1rem 1.3rem;border-radius:0 12px 12px 0;margin:1rem 0}
.highlight-box.info{background:rgba(77,158,232,.1);border-left:5px solid var(--primary)}
.highlight-box.success{background:rgba(52,211,153,.1);border-left:5px solid var(--secondary)}
.highlight-box.warning{background:rgba(251,146,60,.1);border-left:5px solid var(--accent)}
.highlight-box.error{background:rgba(220,38,38,.12);border-left:5px solid #dc2626}

.footer{text-align:center;padding:2rem;color:var(--text3);font-size:.78rem;
  border-top:1px solid var(--border);margin-top:2rem}
.footer .logo{font-weight:700;font-size:.9rem;color:var(--text);margin-bottom:.3rem}
.footer .logo span{color:var(--primary)}

@media(max-width:900px){
  .sidebar{display:none}
  .container{padding:1.2rem}
  .kpi-grid{grid-template-columns:repeat(2,1fr)}
}
"""


# ---------------------------------------------------------------------------
#  JS (theme toggle, per-line copy, Plotly theme)
# ---------------------------------------------------------------------------

_JS = r"""
// Theme — dark by default, overridden by localStorage if user changed it
(function(){
  const saved = localStorage.getItem('fluxos_theme');
  if (saved) document.documentElement.setAttribute('data-theme', saved);
  const btn = document.getElementById('theme-btn');
  if(btn){
    const update = () => {
      const cur = document.documentElement.getAttribute('data-theme');
      btn.textContent = cur === 'dark' ? '☀ Light mode' : '🌙 Dark mode';
    };
    update();
    btn.addEventListener('click', () => {
      const cur = document.documentElement.getAttribute('data-theme');
      const next = cur === 'dark' ? 'light' : 'dark';
      document.documentElement.setAttribute('data-theme', next);
      localStorage.setItem('fluxos_theme', next);
      update();
      // Re-theme all Plotly charts
      if (window.Plotly) {
        document.querySelectorAll('.plot').forEach(p => {
          if (p._fluxos_fig) {
            Plotly.relayout(p, plotlyLayoutPatch());
          }
        });
      }
    });
  }
})();

// Plotly theme helper — reads CSS variables so chart style follows the toggle
function plotlyLayoutPatch(){
  const cs = getComputedStyle(document.documentElement);
  return {
    'paper_bgcolor': cs.getPropertyValue('--plot-bg').trim(),
    'plot_bgcolor':  cs.getPropertyValue('--plot-bg').trim(),
    'font.color':    cs.getPropertyValue('--plot-font').trim(),
    'xaxis.gridcolor': cs.getPropertyValue('--plot-grid').trim(),
    'yaxis.gridcolor': cs.getPropertyValue('--plot-grid').trim(),
  };
}

// Per-line copy buttons (always visible)
(function(){
  const scratch = document.createElement('span');
  const plainOf = (h) => { scratch.innerHTML = h; return scratch.textContent; };
  document.querySelectorAll('.code-block pre').forEach(pre => {
    if (pre.dataset.perLineProcessed === '1') return;
    pre.dataset.perLineProcessed = '1';
    const rawLines = pre.innerHTML.split('\n');
    const parts = [];
    let i = 0;
    while (i < rawLines.length) {
      const lineHTML = rawLines[i];
      const trimmed = plainOf(lineHTML).trim();
      if (trimmed === '') { parts.push('<span class="code-blank">\u00A0</span>'); i++; continue; }
      if (trimmed.startsWith('#')) { parts.push('<span class="code-comment-line">'+lineHTML+'</span>'); i++; continue; }
      const groupH = [lineHTML];
      const groupP = [plainOf(lineHTML)];
      while (groupP[groupP.length-1].trimEnd().endsWith('\\') && i+1 < rawLines.length) {
        i++; groupH.push(rawLines[i]); groupP.push(plainOf(rawLines[i]));
      }
      parts.push('<span class="code-cmd">'+groupH.join('\n')+
        '<button class="line-copy-btn" type="button" title="Copy this command">\u{1F4CB}</button></span>');
      i++;
    }
    pre.innerHTML = parts.join('');
  });
  document.addEventListener('click', (e) => {
    const btn = e.target.closest('.line-copy-btn');
    if (!btn) return;
    e.stopPropagation(); e.preventDefault();
    const cmd = btn.closest('.code-cmd');
    if (!cmd) return;
    const clone = cmd.cloneNode(true);
    clone.querySelectorAll('.line-copy-btn').forEach(b => b.remove());
    navigator.clipboard.writeText(clone.textContent).then(() => {
      const orig = btn.innerHTML;
      btn.textContent = '\u2713';
      btn.classList.add('copied');
      setTimeout(() => { btn.innerHTML = orig; btn.classList.remove('copied'); }, 1500);
    });
  });
})();

// Plotly chart rendering — each .plot element has data-fig JSON on it
(function(){
  function waitForPlotly(cb, tries = 80){
    if (window.Plotly) return cb();
    if (tries <= 0) { console.warn('Plotly failed to load'); return; }
    setTimeout(() => waitForPlotly(cb, tries-1), 100);
  }
  waitForPlotly(() => {
    document.querySelectorAll('.plot').forEach(el => {
      const raw = el.getAttribute('data-fig');
      if (!raw) return;
      const fig = JSON.parse(raw);
      Plotly.newPlot(el, fig.data || [], Object.assign({}, fig.layout || {}, plotlyLayoutPatch()),
                     {displaylogo:false, responsive:true});
      el._fluxos_fig = fig;
    });
  });
})();
"""


# ---------------------------------------------------------------------------
#  Plot builders (return Plotly figure dicts)
# ---------------------------------------------------------------------------

def _fig_time_series(stats: dict) -> dict:
    t = [x / 3600.0 for x in stats["time_s"]]  # seconds → hours
    vol = stats["volume_m3"]
    area = stats["area_flooded_m2"]
    max_h = stats["max_h_m"]
    max_v = stats["max_v_ms"]
    mean_h = stats["mean_h_wet_m"]

    # Four stacked panels. Bottom axis is shared ("x"); each panel has its own y.
    def _tr(y, name, yaxis, color, dash=None):
        tr = dict(x=t, y=y, name=name, mode="lines",
                  line=dict(color=color, width=2), xaxis="x", yaxis=yaxis)
        if dash:
            tr["line"]["dash"] = dash
        return tr

    data = [
        _tr(vol,    "Flood volume",            "y",  "#4d9ee8"),
        _tr(area,   "Flooded area",            "y2", "#34d399"),
        _tr(max_h,  "Max depth",               "y3", "#fb923c"),
        _tr(mean_h, "Mean wet-cell depth",     "y3", "#d4d4d8", dash="dot"),
        _tr(max_v,  "Max velocity",            "y4", "#f472b6"),
    ]
    # Stack four panels vertically with explicit domains; xaxis shared below.
    layout = {
        "xaxis":  {"domain": [0, 1], "anchor": "y4",
                   "title": "Time (hours since start)"},
        "yaxis":  {"domain": [0.78, 1.00], "title": "Volume (m³)"},
        "yaxis2": {"domain": [0.52, 0.74], "title": "Area (m²)"},
        "yaxis3": {"domain": [0.26, 0.48], "title": "Depth (m)"},
        "yaxis4": {"domain": [0.00, 0.22], "title": "Velocity (m/s)"},
        "height": 780,
        "margin": {"l": 80, "r": 30, "t": 30, "b": 60},
        "showlegend": True,
        "legend": {"orientation": "h", "y": 1.05},
    }
    return {"data": data, "layout": layout}


def _fig_max_inundation(stats: dict) -> dict | None:
    import numpy as np
    x = stats.get("map_cell_x")
    y = stats.get("map_cell_y")
    h = stats.get("map_max_h")
    if x is None or h is None or len(h) == 0:
        return None
    h_arr = np.asarray(h)
    # 95th percentile to clip single-cell numerical spikes out of the colorbar
    cmax = float(np.percentile(h_arr[h_arr > 0], 95)) if (h_arr > 0).any() else 1.0
    cmax = max(cmax, 0.1)
    data = [{
        "type": "scattergl",
        "mode": "markers",
        "x": list(map(float, x)),
        "y": list(map(float, y)),
        "marker": {
            "color": list(map(float, h)),
            "colorscale": "Blues",
            "cmin": 0,
            "cmax": cmax,
            "size": 4,
            "colorbar": {"title": "Max depth (m)", "thickness": 14},
            "line": {"width": 0},
        },
        "hovertemplate": "x: %{x:.1f}<br>y: %{y:.1f}<br>max h: %{marker.color:.3f} m<extra></extra>",
        "name": "Max water depth",
    }]
    layout = {
        "xaxis": {"title": "X (m, projected CRS)", "scaleanchor": "y", "scaleratio": 1},
        "yaxis": {"title": "Y (m, projected CRS)"},
        "height": 520,
        "margin": {"l": 80, "r": 30, "t": 30, "b": 50},
        "annotations": [{
            "text": f"Colorbar capped at 95th percentile ({cmax:.2f} m); "
                    f"absolute max = {float(h_arr.max()):.2f} m.",
            "xref": "paper", "yref": "paper", "x": 0, "y": -0.12,
            "showarrow": False, "font": {"size": 10, "color": "#888"},
            "xanchor": "left",
        }] if float(h_arr.max()) > cmax * 1.1 else [],
    }
    return {"data": data, "layout": layout}


def _fig_hazard(stats: dict) -> dict | None:
    x = stats.get("map_cell_x")
    y = stats.get("map_cell_y")
    cls = stats.get("hazard_class")
    if x is None or cls is None or len(cls) == 0:
        return None
    cls = list(map(int, cls))
    labels = stats["hazard_labels"]
    palette = ["#60a5fa", "#fbbf24", "#fb923c", "#ef4444"]
    data = []
    for i, lab in enumerate(labels):
        mask = [c == i for c in cls]
        if not any(mask):
            continue
        data.append({
            "type": "scattergl", "mode": "markers", "name": lab,
            "x": [float(xi) for xi, m in zip(x, mask) if m],
            "y": [float(yi) for yi, m in zip(y, mask) if m],
            "marker": {"color": palette[i], "size": 4, "line": {"width": 0}},
        })
    layout = {
        "xaxis": {"title": "X (m)", "scaleanchor": "y", "scaleratio": 1},
        "yaxis": {"title": "Y (m)"},
        "height": 520,
        "legend": {"orientation": "h", "y": -0.15},
        "margin": {"l": 80, "r": 30, "t": 30, "b": 80},
    }
    return {"data": data, "layout": layout}


def _fig_depth_histogram(stats: dict) -> dict | None:
    edges = stats.get("depth_hist_edges") or []
    counts = stats.get("depth_hist_counts") or []
    if not edges or not counts:
        return None
    centres = [(edges[i] + edges[i + 1]) / 2 for i in range(len(counts))]
    widths = [edges[i + 1] - edges[i] for i in range(len(counts))]
    data = [{
        "type": "bar", "x": centres, "y": counts, "width": widths,
        "marker": {"color": "#4d9ee8", "line": {"color": "#0066cc", "width": 1}},
        "hovertemplate": "h: %{x:.2f} m<br>cells: %{y}<extra></extra>",
        "name": "Max depth",
    }]
    layout = {
        "xaxis": {"title": "Maximum water depth (m)"},
        "yaxis": {"title": "Cells"},
        "height": 360,
        "margin": {"l": 70, "r": 30, "t": 20, "b": 50},
        "bargap": 0.05,
    }
    return {"data": data, "layout": layout}


def _fig_first_wet(stats: dict) -> dict | None:
    x = stats.get("map_cell_x")
    y = stats.get("map_cell_y")
    fw = stats.get("map_first_wet_time_s")
    if x is None or fw is None or len(fw) == 0:
        return None
    fw_h = [float(t) / 3600.0 if t == t else None for t in fw]  # NaN-safe
    # Separate never-wet from wet
    x_wet, y_wet, t_wet = [], [], []
    for xi, yi, ti in zip(x, y, fw_h):
        if ti is not None:
            x_wet.append(float(xi)); y_wet.append(float(yi)); t_wet.append(ti)
    if not t_wet:
        return None
    data = [{
        "type": "scattergl", "mode": "markers",
        "x": x_wet, "y": y_wet,
        "marker": {
            "color": t_wet, "colorscale": "Viridis",
            "cmin": min(t_wet), "cmax": max(t_wet),
            "size": 4, "line": {"width": 0},
            "colorbar": {"title": "First wet (h)", "thickness": 14},
        },
        "hovertemplate": "x: %{x:.1f}<br>y: %{y:.1f}<br>first wet: %{marker.color:.2f} h<extra></extra>",
        "name": "First inundation",
    }]
    layout = {
        "xaxis": {"title": "X (m)", "scaleanchor": "y", "scaleratio": 1},
        "yaxis": {"title": "Y (m)"},
        "height": 520,
        "margin": {"l": 80, "r": 30, "t": 30, "b": 50},
    }
    return {"data": data, "layout": layout}


def _fig_chemical(stats: dict) -> dict | None:
    if not stats.get("has_conc"):
        return None
    t = [x / 3600.0 for x in stats["time_s"]]
    max_c = stats.get("max_conc_SW_t") or []
    mean_c = stats.get("mean_conc_SW_t") or []
    data = [
        dict(x=t, y=max_c,  name="Max conc_SW",  mode="lines",
             line={"color": "#fb923c", "width": 2}),
        dict(x=t, y=mean_c, name="Mean conc_SW (wet)", mode="lines",
             line={"color": "#d4d4d8", "width": 2, "dash": "dot"}),
    ]
    layout = {
        "xaxis": {"title": "Time (hours)"},
        "yaxis": {"title": "Concentration (mg/L)"},
        "height": 320,
        "legend": {"orientation": "h", "y": -0.2},
        "margin": {"l": 70, "r": 30, "t": 20, "b": 60},
    }
    return {"data": data, "layout": layout}


# ---------------------------------------------------------------------------
#  HTML building blocks
# ---------------------------------------------------------------------------

def _esc(s) -> str:
    return html.escape(str(s), quote=True)


def _fmt_int(x):
    try:
        return f"{int(round(float(x))):,}"
    except Exception:
        return str(x)


def _fmt_num(x, fmt=".3f"):
    try:
        return format(float(x), fmt)
    except Exception:
        return str(x)


def _fmt_duration(seconds: int) -> str:
    s = int(seconds)
    h, rem = divmod(s, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h} h {m:02d} min"
    return f"{m} min {s:02d} s"


def _kpi(icon: str, value, label: str, sub: str = "") -> str:
    sub_html = f'<div class="sub">{_esc(sub)}</div>' if sub else ""
    return (f'<div class="kpi"><div class="icon">{icon}</div>'
            f'<div class="value">{_esc(value)}</div>'
            f'<div class="label">{_esc(label)}</div>{sub_html}</div>')


def _plot_div(plot_id: str, fig: dict, wide: bool = False) -> str:
    klass = "plot plot-wide" if wide else "plot"
    encoded = _esc(json.dumps(fig, default=float))
    return f'<div id="{plot_id}" class="{klass}" data-fig="{encoded}"></div>'


# ---------------------------------------------------------------------------
#  Main report
# ---------------------------------------------------------------------------

def generate_results_report(report_data: dict, output_path: str) -> str:
    """
    Build a self-contained HTML results report.

    ``report_data`` is expected to have:
      - ``config``        : the user's _config dict (for project name / metadata)
      - ``modset``        : parsed modset.json
      - ``mesh_type``     : "regular" or "triangular"
      - ``results_dir``   : directory that was analysed
      - ``stats``         : dict from ``flood_statistics.analyse``
      - ``kpis``          : dict from ``flood_statistics.summary_kpis``
      - ``generated_at``  : ISO timestamp (optional)
      - ``errors``        : list (optional)
    """
    config = report_data.get("config", {})
    kpis = report_data["kpis"]
    stats = report_data["stats"]
    modset = report_data.get("modset", {})
    generated_at = report_data.get("generated_at") \
        or datetime.now().isoformat(timespec="seconds")

    authors = ", ".join(config.get("authors", [])) or "—"

    # --- KPIs grid --------------------------------------------------------
    kpi_html = ''.join([
        _kpi("🌊", _fmt_int(kpis["peak_volume_m3"]), "Peak flood volume", "m³"),
        _kpi("🗺️", _fmt_int(kpis["peak_area_m2"]), "Peak flooded area", "m²"),
        _kpi("📏", _fmt_num(kpis["peak_depth_m"], ".2f"), "Max water depth", "metres"),
        _kpi("🌀", _fmt_num(kpis["peak_velocity_ms"], ".2f"), "Max velocity", "m/s"),
        _kpi("⏱️", _fmt_duration(kpis["peak_time_s"]), "Time to peak", "since start"),
        _kpi("🕰️", _fmt_duration(kpis["sim_duration_s"]), "Sim duration", ""),
        _kpi("🔢", _fmt_int(kpis["wet_cells"]), "Wet cells", "any time during run"),
        _kpi("📊", _fmt_int(stats["n_timesteps"]), "Timesteps analysed", ""),
    ])

    # --- Hazard distribution table ---------------------------------------
    haz_rows = ''.join(
        f'<tr><td>{_esc(lab)}</td><td class="num">{_fmt_int(c)}</td></tr>'
        for lab, c in zip(stats.get("hazard_labels", []),
                          stats.get("hazard_counts", []))
    )

    # --- Plots ------------------------------------------------------------
    fig_ts = _fig_time_series(stats)
    fig_max = _fig_max_inundation(stats)
    fig_haz = _fig_hazard(stats)
    fig_hist = _fig_depth_histogram(stats)
    fig_first = _fig_first_wet(stats)
    fig_chem = _fig_chemical(stats)

    # --- Navigation -------------------------------------------------------
    nav_items = [("summary", "Summary"), ("timeseries", "Time series"),
                 ("maxmap", "Max inundation"), ("hazard", "Hazard"),
                 ("histogram", "Depth histogram"),
                 ("firstwet", "First inundation")]
    if fig_chem:
        nav_items.append(("chemical", "Chemical transport"))
    nav_items.append(("meta", "Run info"))
    nav_html = ''.join(f'<a href="#{sid}">{_esc(label)}</a>'
                       for sid, label in nav_items)

    # --- Sections ---------------------------------------------------------
    summary = f"""
    <section class="section" id="summary">
      <h2>Summary</h2>
      <div class="kpi-grid">{kpi_html}</div>
      <div class="card">
        <h3>Hazard classification (ARR-2019 H·V criterion)</h3>
        <div class="table-wrap"><table>
          <thead><tr><th>Class</th><th>Cells</th></tr></thead>
          <tbody>{haz_rows}</tbody>
        </table></div>
      </div>
    </section>"""

    timeseries = f"""
    <section class="section" id="timeseries">
      <h2>Time series</h2>
      <div class="card">
        <p style="color:var(--text2);margin-bottom:.3rem">
          Per-timestep aggregates over the full domain. Volume and area quantify
          the overall flood extent; max/mean depth and max velocity characterise
          the severity.</p>
        {_plot_div("ts-plot", fig_ts, wide=True)}
      </div>
    </section>"""

    maxmap = f"""
    <section class="section" id="maxmap">
      <h2>Maximum inundation map</h2>
      <div class="card">
        <p style="color:var(--text2);margin-bottom:.3rem">
          Per-cell maximum water depth observed at any timestep — the single
          most important output for floodplain delimitation and post-event
          analysis.</p>
        {_plot_div("max-plot", fig_max) if fig_max else '<em>No spatial data.</em>'}
      </div>
    </section>"""

    hazard = f"""
    <section class="section" id="hazard">
      <h2>Flood-hazard classification</h2>
      <div class="card">
        <p style="color:var(--text2);margin-bottom:.3rem">
          Cells classified by the hazard factor <code>D = h·(v + 0.5)</code>
          (Australian Rainfall &amp; Runoff 2019). H1 (Low) to H4 (Extreme).</p>
        {_plot_div("haz-plot", fig_haz) if fig_haz else '<em>No velocity data.</em>'}
      </div>
    </section>"""

    histogram = f"""
    <section class="section" id="histogram">
      <h2>Max-depth distribution</h2>
      <div class="card">
        <p style="color:var(--text2);margin-bottom:.3rem">
          Histogram of max depth across cells that were ever wet — shows
          whether the flood is predominantly shallow spreading or deep channelised.</p>
        {_plot_div("hist-plot", fig_hist) if fig_hist else '<em>No data.</em>'}
      </div>
    </section>"""

    firstwet = f"""
    <section class="section" id="firstwet">
      <h2>Time-of-first-inundation</h2>
      <div class="card">
        <p style="color:var(--text2);margin-bottom:.3rem">
          First timestep at which each cell exceeded the wet-depth threshold —
          reveals the flood-front propagation that the WebGL animation shows
          dynamically.</p>
        {_plot_div("fw-plot", fig_first) if fig_first else '<em>No data.</em>'}
      </div>
    </section>"""

    chemical_section = ""
    if fig_chem:
        chemical_section = f"""
        <section class="section" id="chemical">
          <h2>Chemical transport (conc_SW)</h2>
          <div class="card">
            <p style="color:var(--text2);margin-bottom:.3rem">
              Domain-wide statistics of the transported solute concentration.
              Enabled because <code>ADE_TRANSPORT.STATUS = true</code> in the modset.</p>
            {_plot_div("chem-plot", fig_chem)}
          </div>
        </section>"""

    meta_rows = [
        ("Project",        _esc(config.get("project_name") or "—")),
        ("Authors",        _esc(authors)),
        ("Mesh type",      _esc(report_data.get("mesh_type"))),
        ("Results folder", f'<code>{_esc(report_data.get("results_dir"))}</code>'),
        ("Modset",         f'<code>{_esc(report_data.get("modset_path") or "—")}</code>'),
        ("Timesteps",      _fmt_int(stats["n_timesteps"])),
        ("Depth threshold (wet)",
                           f'{_fmt_num(config.get("h_threshold_flooded_m") or 0.01, ".3f")} m'),
        ("Generated",      _esc(generated_at)),
    ]
    meta_html = ''.join(f'<tr><td>{k}</td><td>{v}</td></tr>' for k, v in meta_rows)

    meta = f"""
    <section class="section" id="meta">
      <h2>Run info</h2>
      <div class="card"><div class="table-wrap"><table><tbody>{meta_html}</tbody></table></div></div>
    </section>"""

    subtitle = (f"Flood statistics report — "
                f"{stats['n_timesteps']} timesteps, {_fmt_duration(kpis['sim_duration_s'])}")

    # --- Assemble HTML ----------------------------------------------------
    doc = f"""<!DOCTYPE html>
<html lang="en" data-theme="dark">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>FLUXOS results — {_esc(config.get('project_name') or 'Simulation')}</title>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
  <script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
  <style>{_CSS}</style>
</head>
<body>
  <div class="layout">
    <aside class="sidebar">
      <div class="logo">FLUX<span>OS</span></div>
      <nav>{nav_html}</nav>
      <div class="theme-toggle">
        <button id="theme-btn" class="theme-btn" type="button">☀ Light mode</button>
      </div>
    </aside>
    <main class="main">
      <header class="header">
        <h1>FLUX<span>OS</span> results — {_esc(config.get('project_name') or 'Simulation')}</h1>
        <p class="subtitle">{_esc(subtitle)}</p>
        <div class="meta">
          <span>📅 {_esc(generated_at[:10])}</span>
          <span>👤 {_esc(authors)}</span>
          <span>🔺 {_esc(report_data.get('mesh_type'))} mesh</span>
        </div>
      </header>
      <div class="container">
        {summary}
        {timeseries}
        {maxmap}
        {hazard}
        {histogram}
        {firstwet}
        {chemical_section}
        {meta}
      </div>
      <footer class="footer">
        <div class="logo">FLUX<span>OS</span></div>
        Generated {_esc(generated_at)} —
        <a href="https://github.com/ue-hydro/FLUXOS_cpp">github.com/ue-hydro/FLUXOS_cpp</a>
      </footer>
    </main>
  </div>
  <script>{_JS}</script>
</body>
</html>
"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(doc)
    return output_path
