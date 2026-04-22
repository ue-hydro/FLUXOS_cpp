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
FLUXOS model-config HTML report generator.

Produces a single self-contained HTML file summarising what the driver built
(DEM, mesh, modset.json) and giving the user copy-paste Docker commands to
build the image, run the simulation, check outputs, and visualise results.

Style: Inter + JetBrains Mono, #0066cc primary / #00a86b secondary /
#ff6b35 accent, dark/light theme toggle, KPI grid, sidebar navigation.
"""

from __future__ import annotations

import html
import json
import os
import platform
from datetime import datetime


# ---------------------------------------------------------------------------
#  CSS
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
.sidebar .logo{font-weight:700;font-size:1.1rem;color:var(--dark);margin-bottom:1.2rem;
  padding-left:.7rem;letter-spacing:-.3px}
[data-theme="dark"] .sidebar .logo{color:var(--text)}
.sidebar .logo span{color:var(--primary)}
.sidebar nav{display:flex;flex-direction:column;gap:2px}
.sidebar nav a{display:block;padding:.4rem .7rem;border-radius:6px;font-size:.82rem;
  color:var(--text2);transition:all .15s}
.sidebar nav a:hover{background:rgba(0,102,204,.06);color:var(--primary)}
.sidebar nav a.active{background:rgba(0,102,204,.1);color:var(--primary);font-weight:600}
.theme-toggle{margin-top:auto;padding-top:1rem;border-top:1px solid var(--border)}
.theme-btn{display:flex;align-items:center;gap:.5rem;width:100%;padding:.45rem .7rem;
  border:1px solid var(--border);border-radius:8px;background:var(--surface);
  color:var(--text2);font-size:.78rem;cursor:pointer;transition:all .15s;font-family:inherit}
.theme-btn:hover{border-color:var(--primary);color:var(--primary)}

.main{flex:1;min-width:0;overflow-x:hidden}
.header{background:linear-gradient(135deg,var(--dark) 0%,#16213e 50%,#0f3460 100%);
  color:#fff;padding:3rem 2.5rem 2rem;position:relative;overflow:hidden}
.header::before{content:'';position:absolute;inset:0;
  background:radial-gradient(circle at 80% 20%,rgba(255,255,255,.08) 0%,transparent 60%)}
.header h1{font-size:2rem;font-weight:700;letter-spacing:-.5px;position:relative}
.header h1 span{color:var(--primary)}
.header .subtitle{opacity:.8;font-size:.92rem;margin-top:.4rem;position:relative;max-width:800px}
.header .meta{display:flex;gap:.6rem;margin-top:1rem;font-size:.78rem;opacity:.8;
  flex-wrap:wrap;position:relative}
.header .meta span{background:rgba(255,255,255,.1);padding:.25rem .7rem;border-radius:20px}

.container{max-width:1100px;margin:0 auto;padding:2rem 2rem 3rem}
.section{margin-bottom:2.5rem;scroll-margin-top:1rem}
.section h2{font-size:1.2rem;font-weight:700;color:var(--text);margin-bottom:1rem;
  display:flex;align-items:center;gap:.5rem;letter-spacing:-.3px}
.section h2::before{content:'';width:5px;height:1.2em;border-radius:3px;
  background:linear-gradient(180deg,var(--primary),var(--primary-dark));flex-shrink:0}
.section h3{font-size:1rem;font-weight:600;margin:1rem 0 .6rem;color:var(--text)}

.card{background:var(--surface);border:1px solid var(--border);border-radius:16px;
  padding:1.5rem;box-shadow:var(--shadow);transition:transform .3s,box-shadow .3s}
.card + .card{margin-top:1rem}
.card.primary{border-left:5px solid var(--primary)}
.card.secondary{border-left:5px solid var(--secondary)}
.card.accent{border-left:5px solid var(--accent)}
.card h3{margin-top:0;margin-bottom:.8rem;font-weight:600}

.kpi-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:1rem;
  margin-bottom:1.5rem}
.kpi{background:var(--surface);border:1px solid var(--border);border-radius:16px;
  padding:1.2rem 1rem;text-align:center;box-shadow:var(--shadow);
  transition:transform .3s,box-shadow .3s;position:relative;overflow:hidden}
.kpi:hover{transform:translateY(-3px);box-shadow:var(--shadow-lg)}
.kpi::before{content:'';position:absolute;top:0;left:0;right:0;height:3px;
  background:linear-gradient(90deg,var(--primary),var(--secondary))}
.kpi .icon{font-size:1.6rem;margin-bottom:.4rem}
.kpi .value{font-size:1.4rem;font-weight:700;color:var(--primary);
  font-family:'JetBrains Mono',monospace;overflow-wrap:break-word;word-break:break-word}
.kpi .label{font-size:.7rem;color:var(--text3);text-transform:uppercase;
  letter-spacing:.8px;margin-top:.3rem}

.table-wrap{overflow-x:auto;border-radius:12px;border:1px solid var(--border)}
table{width:100%;border-collapse:collapse;font-size:.88rem;background:var(--surface)}
th{background:var(--dark);color:#fff;padding:.65rem .9rem;text-align:left;font-weight:600;
  font-size:.72rem;text-transform:uppercase;letter-spacing:.5px}
[data-theme="dark"] th{background:#252640}
td{padding:.55rem .9rem;border-bottom:1px solid var(--border);vertical-align:top}
tr:last-child td{border-bottom:none}
tr:nth-child(even){background:rgba(0,102,204,.03)}
[data-theme="dark"] tr:nth-child(even){background:rgba(77,158,232,.05)}
td.num{text-align:right;font-family:'JetBrains Mono',monospace;font-size:.82rem}
td code{font-family:'JetBrains Mono',monospace;font-size:.82rem;
  background:rgba(0,102,204,.05);padding:.1rem .3rem;border-radius:4px}

.badge{display:inline-flex;align-items:center;gap:.25rem;padding:.2rem .6rem;
  border-radius:20px;font-size:.72rem;font-weight:600;white-space:nowrap}
.badge-primary{background:rgba(0,102,204,.1);color:var(--primary)}
.badge-secondary{background:rgba(0,168,107,.1);color:var(--secondary)}
.badge-accent{background:rgba(255,107,53,.1);color:var(--accent)}
.badge-none{background:rgba(0,0,0,.04);color:var(--text3)}
[data-theme="dark"] .badge-none{background:rgba(255,255,255,.06)}

.highlight-box{padding:1.2rem 1.5rem;border-radius:0 12px 12px 0;margin:1rem 0}
.highlight-box.info{background:linear-gradient(135deg,#eff6ff,#dbeafe);border-left:5px solid var(--primary)}
.highlight-box.success{background:linear-gradient(135deg,#ecfdf5,#d1fae5);border-left:5px solid var(--secondary)}
.highlight-box.warning{background:linear-gradient(135deg,#fff7ed,#ffedd5);border-left:5px solid var(--accent)}
.highlight-box.error{background:linear-gradient(135deg,#fef2f2,#fee2e2);border-left:5px solid #dc2626}
[data-theme="dark"] .highlight-box.info{background:rgba(0,102,204,.08)}
[data-theme="dark"] .highlight-box.success{background:rgba(0,168,107,.08)}
[data-theme="dark"] .highlight-box.warning{background:rgba(255,107,53,.08)}
[data-theme="dark"] .highlight-box.error{background:rgba(220,38,38,.12)}

.code-block{position:relative;margin:.6rem 0 1rem}
.code-block pre{background:var(--code-bg);color:var(--code-text);padding:1rem 1rem 1rem 1rem;
  border-radius:10px;overflow-x:auto;font-family:'JetBrains Mono',monospace;font-size:.82rem;
  line-height:1.6}
/* Per-line copy buttons — always visible, one per command (no block-level button) */
.code-block pre .code-cmd{display:block;position:relative;border-radius:3px;padding-right:2.2rem}
.code-block pre .code-cmd:hover{background:rgba(255,255,255,.06)}
.code-block pre .code-blank,.code-block pre .code-comment-line{display:block}
.code-block pre .line-copy-btn{position:absolute;top:0;right:.25rem;padding:0 .45rem;font-size:.68rem;
  background:rgba(255,255,255,.1);border:1px solid rgba(255,255,255,.18);border-radius:4px;
  color:#e2e8f0;cursor:pointer;opacity:1;transition:background .15s,border-color .15s;
  font-family:inherit;line-height:1.55}
.code-block pre .line-copy-btn:hover{background:rgba(255,255,255,.22);border-color:rgba(255,255,255,.35)}
.code-block pre .line-copy-btn.copied{background:var(--secondary);border-color:var(--secondary);color:#fff}

/* Visualisation toggle (checkbox for --variable conc_SW) */
.viz-toggle-card{margin-top:.4rem}
.viz-toggle-label{display:inline-flex;align-items:center;gap:.45rem;cursor:pointer;
  font-size:.88rem;color:var(--text);margin:.5rem 0 .3rem}
.viz-toggle-label input[type="checkbox"]{cursor:pointer;accent-color:var(--primary)}
.viz-toggle-label code{font-family:'JetBrains Mono',monospace;font-size:.78rem;
  background:rgba(0,102,204,.08);padding:.05rem .35rem;border-radius:3px;color:var(--primary)}

/* 2D / 3D map view toggle */
.map-view-toggle{display:inline-flex;border:1px solid var(--border);border-radius:8px;
  overflow:hidden;background:var(--surface);flex-shrink:0}
.map-view-toggle .view-btn{border:0;background:transparent;color:var(--text2);
  padding:.4rem .9rem;cursor:pointer;font-family:inherit;font-size:.82rem;
  font-weight:500;transition:background .15s,color .15s}
.map-view-toggle .view-btn:not(:last-child){border-right:1px solid var(--border)}
.map-view-toggle .view-btn:hover{background:rgba(0,102,204,.08);color:var(--primary)}
.map-view-toggle .view-btn.active{background:var(--primary);color:#fff}
[data-theme="dark"] .map-view-toggle .view-btn.active{background:var(--primary);color:#0f1117}

/* Vertical-exaggeration slider (3D view only) */
.ve-slider-wrap{display:flex;align-items:center;gap:.8rem;margin:.2rem 0 .6rem;
  padding:.5rem .8rem;background:var(--surface);border:1px solid var(--border);
  border-radius:8px;flex-wrap:wrap}
.ve-slider-wrap label{font-size:.82rem;color:var(--text2);display:flex;
  align-items:center;gap:.6rem;flex-wrap:wrap}
.ve-slider-wrap .ve-label{font-family:'JetBrains Mono',monospace;font-size:.9rem;
  color:var(--primary);font-weight:700;min-width:3em}
.ve-slider-wrap .ve-note{font-size:.72rem;color:var(--text3);font-weight:400}
.ve-slider-wrap input[type="range"]{flex:1 1 200px;min-width:120px;
  accent-color:var(--primary);cursor:pointer}

.os-badge{display:inline-block;padding:3px 12px;border-radius:14px;font-size:.75rem;
  font-weight:600;color:#fff;margin-right:.5rem}

details.module-details{margin-top:.6rem}
details.module-details>summary{cursor:pointer;padding:.55rem .9rem;border-radius:8px;
  background:var(--surface);border:1px solid var(--border);font-weight:600;font-size:.85rem;
  transition:all .15s;list-style:none;display:flex;align-items:center;gap:.5rem}
details.module-details>summary::before{content:'\25B8';font-size:.85rem;transition:transform .2s}
details[open].module-details>summary::before{transform:rotate(90deg)}
details.module-details>summary:hover{border-color:var(--primary);color:var(--primary)}
.module-content{padding:.8rem 1rem;border:1px solid var(--border);border-top:none;
  border-radius:0 0 10px 10px;background:var(--surface)}

.footer{text-align:center;padding:2rem;color:var(--text3);font-size:.78rem;
  border-top:1px solid var(--border);margin-top:2rem}
.footer .logo{font-weight:700;font-size:.9rem;color:var(--text);margin-bottom:.3rem}
.footer .logo span{color:var(--primary)}

@media(max-width:900px){
  .sidebar{display:none}
  .container{padding:1.2rem}
  .kpi-grid{grid-template-columns:repeat(2,1fr)}
  .header{padding:2rem 1.5rem 1.5rem}
}
"""


# ---------------------------------------------------------------------------
#  JS (theme toggle, copy buttons, scrollspy)
# ---------------------------------------------------------------------------

_JS = r"""
// Theme toggle — dark by default, overridden by localStorage if user changed it
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
    });
  }
})();

// Per-line copy buttons — always visible, one per non-comment, non-blank command.
// Backslash-continuation lines are grouped into a single logical command.
(function(){
  const scratch = document.createElement('span');
  const plainOf = (html) => { scratch.innerHTML = html; return scratch.textContent; };

  function injectLineButtons(pre){
    if (!pre || pre.dataset.perLineProcessed === '1') return;
    pre.dataset.perLineProcessed = '1';
    const rawLines = pre.innerHTML.split('\n');
    const parts = [];
    let i = 0;
    while (i < rawLines.length) {
      const lineHTML = rawLines[i];
      const linePlain = plainOf(lineHTML);
      const trimmed = linePlain.trim();
      if (trimmed === '') {
        parts.push('<span class="code-blank">\u00A0</span>');
        i++;
      } else if (trimmed.startsWith('#')) {
        parts.push('<span class="code-comment-line">' + lineHTML + '</span>');
        i++;
      } else {
        const htmlLines = [lineHTML];
        const plainLines = [linePlain];
        while (plainLines[plainLines.length - 1].trimEnd().endsWith('\\') && i + 1 < rawLines.length) {
          i++;
          htmlLines.push(rawLines[i]);
          plainLines.push(plainOf(rawLines[i]));
        }
        parts.push('<span class="code-cmd">' + htmlLines.join('\n') +
          '<button class="line-copy-btn" type="button" title="Copy this command">\u{1F4CB}</button></span>');
        i++;
      }
    }
    pre.innerHTML = parts.join('');
  }

  // Initial pass over all pre blocks on the page
  document.querySelectorAll('.code-block pre').forEach(injectLineButtons);

  // Expose so dynamic snippets (e.g. visualisation toggle) can re-run
  window.__fluxosInjectLineButtons = injectLineButtons;

  // Per-line button click → copy the command text (strip the button itself)
  document.addEventListener('click', (e) => {
    const btn = e.target.closest('.line-copy-btn');
    if (!btn) return;
    e.stopPropagation();
    e.preventDefault();
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

// Plotly renderer — scans .plot elements carrying a data-fig JSON payload
// and draws them once Plotly.js finishes loading. Theme-aware.
(function(){
  function layoutPatch(){
    const cs = getComputedStyle(document.documentElement);
    const bg = cs.getPropertyValue('--surface').trim() || '#fff';
    const text = cs.getPropertyValue('--text').trim() || '#2d3436';
    const grid = cs.getPropertyValue('--border').trim() || '#e2e8f0';
    return {
      paper_bgcolor: bg,
      plot_bgcolor: bg,
      'font.color': text,
      'xaxis.gridcolor': grid,
      'yaxis.gridcolor': grid,
    };
  }
  function waitForPlotly(cb, tries = 80){
    if (window.Plotly) return cb();
    if (tries <= 0) { console.warn('Plotly CDN failed to load'); return; }
    setTimeout(() => waitForPlotly(cb, tries - 1), 100);
  }
  function renderPlot(el) {
    try {
      const fig = JSON.parse(el.getAttribute('data-fig'));
      // Don't apply the 2D layoutPatch to 3D scene plots — they use
      // layout.scene.{xaxis,yaxis,zaxis} instead of top-level axes.
      const is3d = (fig.data || []).some(t => (t.type || '').includes('3d')
                                            || t.type === 'surface');
      const layoutOpts = Object.assign({}, fig.layout || {},
                                       is3d ? {} : layoutPatch());
      Plotly.newPlot(el, fig.data || [], layoutOpts,
                     { displaylogo: false, responsive: true });
      el._fluxosFig = fig;
      el._rendered = true;
    } catch (e) { console.error('plot render failed', e); }
  }
  waitForPlotly(() => {
    // Skip elements flagged as lazy (rendered on first user interaction)
    document.querySelectorAll('.plot[data-fig]:not([data-lazy])').forEach(renderPlot);
  });
  // Expose so other UI (e.g. 2D/3D toggle) can lazy-render on demand.
  window.__fluxosRenderPlot = renderPlot;
  // Re-theme charts when the user toggles dark / light
  const btn = document.getElementById('theme-btn');
  if (btn) {
    btn.addEventListener('click', () => {
      if (!window.Plotly) return;
      document.querySelectorAll('.plot[data-fig]').forEach(el => {
        if (el._fluxosFig) Plotly.relayout(el, layoutPatch());
      });
    });
  }
})();

// DEM map: 2D / 3D view toggle + vertical-exaggeration slider
(function(){
  const toggle = document.querySelector('.map-view-toggle');
  if (!toggle) return;
  const plot2d = document.getElementById('dem-mesh-plot-2d');
  const plot3d = document.getElementById('dem-mesh-plot-3d');
  const slider = document.getElementById('ve-slider');
  const sliderWrap = document.querySelector('.ve-slider-wrap');
  const veLabel = sliderWrap && sliderWrap.querySelector('.ve-label');
  if (!plot2d || !plot3d) return;

  toggle.addEventListener('click', (e) => {
    const btn = e.target.closest('.view-btn');
    if (!btn) return;
    const want3d = btn.dataset.view === '3d';
    toggle.querySelectorAll('.view-btn').forEach(b => {
      const isActive = b === btn;
      b.classList.toggle('active', isActive);
      b.setAttribute('aria-selected', isActive ? 'true' : 'false');
    });
    plot2d.style.display = want3d ? 'none' : '';
    plot3d.style.display = want3d ? '' : 'none';
    if (sliderWrap) sliderWrap.style.display = want3d ? '' : 'none';
    // Lazy-render 3D on first activation (spares the initial-load WebGL cost)
    if (want3d && !plot3d._rendered && window.__fluxosRenderPlot) {
      plot3d.removeAttribute('data-lazy');
      window.__fluxosRenderPlot(plot3d);
    }
    // Plotly needs a resize notification after a display:none → '' flip
    if (window.Plotly) {
      const target = want3d ? plot3d : plot2d;
      try { Plotly.Plots.resize(target); } catch (_) { /* no-op */ }
    }
  });

  // VE slider — live-update scene.aspectratio.z as the user drags
  if (slider && veLabel) {
    const dz = parseFloat(slider.dataset.dz) || 1.0;
    const horiz = parseFloat(slider.dataset.horiz) || 1.0;
    const apply = () => {
      const ve = parseFloat(slider.value) || 1;
      veLabel.textContent = ve + '×';
      if (plot3d._rendered && window.Plotly) {
        const zAspect = Math.max(0.01, Math.min(3.0, ve * dz / horiz));
        try {
          Plotly.relayout(plot3d, { 'scene.aspectratio.z': zAspect });
        } catch (_) { /* no-op */ }
      }
    };
    slider.addEventListener('input', apply);
  }
})();

// Visualisation toggle: swap the viewer command between base and chemical-transport variants
(function(){
  const cb = document.getElementById('viz-conc-toggle');
  if (!cb) return;
  const baseBlock = document.querySelector('.viz-cmd-base');
  const transportBlock = document.querySelector('.viz-cmd-transport');
  if (!baseBlock || !transportBlock) return;
  const apply = () => {
    baseBlock.style.display = cb.checked ? 'none' : 'block';
    transportBlock.style.display = cb.checked ? 'block' : 'none';
  };
  cb.addEventListener('change', apply);
  apply();
})();

// Scrollspy: highlight the sidebar link for the section currently in view
(function(){
  const links = document.querySelectorAll('.sidebar nav a');
  const map = {};
  links.forEach(a => { const id = a.getAttribute('href').slice(1); map[id] = a; });
  const observer = new IntersectionObserver(entries => {
    entries.forEach(e => {
      if(e.isIntersecting){
        links.forEach(a => a.classList.remove('active'));
        const a = map[e.target.id];
        if(a) a.classList.add('active');
      }
    });
  }, { rootMargin: '-20% 0px -70% 0px' });
  document.querySelectorAll('.section').forEach(s => observer.observe(s));
})();
"""


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

def _esc(s) -> str:
    return html.escape(str(s), quote=True)


def _code_block(cmd: str) -> str:
    """Render a code block. Per-line copy buttons are injected at render time by JS."""
    return (
        '<div class="code-block">'
        f'<pre>{_esc(cmd)}</pre>'
        '</div>'
    )


def _kpi(icon: str, value, label: str) -> str:
    return (
        '<div class="kpi">'
        f'<div class="icon">{icon}</div>'
        f'<div class="value">{_esc(value)}</div>'
        f'<div class="label">{_esc(label)}</div>'
        '</div>'
    )


def _badge(kind: str, text: str) -> str:
    return f'<span class="badge badge-{kind}">{_esc(text)}</span>'


def _detect_os() -> tuple[str, str, str]:
    """Return (key, label, badge_color_hex)."""
    s = platform.system().lower()
    if s.startswith("win"):
        return "windows", "Windows", "#0078d4"
    if s == "darwin":
        return "macos", "macOS", "#111111"
    return "linux", "Linux", "#2d7d46"


def _open_outputs_cmd(os_key: str, results_dir_abs: str) -> str:
    if os_key == "windows":
        return f'explorer "{results_dir_abs}"'
    if os_key == "macos":
        return f'open "{results_dir_abs}"'
    return f'xdg-open "{results_dir_abs}"'


# ---------------------------------------------------------------------------
#  Report sections
# ---------------------------------------------------------------------------

def _header(config: dict, generated_at: str) -> str:
    authors = ", ".join(config.get("authors", [])) or "—"
    return f"""
    <header class="header">
      <h1>FLUX<span>OS</span> — {_esc(config['project_name'])}</h1>
      <p class="subtitle">{_esc(config.get('description') or '')}</p>
      <div class="meta">
        <span>📅 {_esc(config.get('date') or generated_at[:10])}</span>
        <span>👤 {_esc(authors)}</span>
        <span>🧭 UTM zone {_esc(config.get('dem_utm_zone'))}</span>
        <span>🧪 {_esc(config['mesh_type'])} mesh</span>
      </div>
    </header>
    """


def _summary_section(data: dict) -> str:
    config = data["config"]
    dem = data.get("dem") or {}
    mesh = data.get("mesh") or {}
    modules_on = sum([
        bool(config["ade_transport"]["enabled"]),
        bool(config["soil_infiltration"]["enabled"]),
    ])

    kpis = [
        _kpi("🗺️", config["mesh_type"].capitalize(), "Mesh type"),
        _kpi("📐", f'{dem.get("target_resolution_m", "?")} m', "DEM resolution"),
    ]
    if config["mesh_type"] == "triangular":
        kpis.append(_kpi("🔺", f'{mesh.get("cell_count", "—"):,}'
                         if isinstance(mesh.get("cell_count"), int)
                         else "—", "Mesh cells"))
        kpis.append(_kpi("📍", f'{mesh.get("vertex_count", "—"):,}'
                         if isinstance(mesh.get("vertex_count"), int)
                         else "—", "Vertices"))
    else:
        kpis.append(_kpi("🟦", f'{dem.get("output_cols", "—")} × {dem.get("output_rows", "—")}',
                         "Grid size"))
    kpis.append(_kpi("⏱️", f'{config["print_step_s"]} s', "Print step"))
    kpis.append(_kpi("🧩", f'{modules_on}/3', "Modules active"))

    return f"""
    <section class="section" id="summary">
      <h2>Summary</h2>
      <div class="kpi-grid">{''.join(kpis)}</div>
    </section>
    """


def _fig_dem_mesh_map(dem: dict, mesh: dict | None,
                      mesh_type: str) -> dict | None:
    """
    Build two Plotly figure dicts for the preview map:

    - ``fig2d``: top-down ``heatmap`` of the DEM + ``scattergl`` triangle edges
    - ``fig3d``: ``surface`` of the DEM + ``scatter3d`` triangle edges lifted
      onto the vertex elevations

    Returns ``{"fig2d": …, "fig3d": …}`` or ``None`` if the DEM preview is
    missing. The 3D figure is included regardless of mesh type (for a
    regular mesh it is just the DEM surface).
    """
    preview = (dem or {}).get("preview")
    if not preview:
        return None

    # -- 2D heatmap + mesh overlay ----------------------------------------
    data2d = [{
        "type": "heatmap",
        "z": preview["z"],
        "x": preview["x_vals"],
        "y": preview["y_vals"],
        "colorscale": "earth",
        "colorbar": {"title": "Elevation (m)", "thickness": 14, "xpad": 0},
        "hovertemplate": "x: %{x:.0f}<br>y: %{y:.0f}<br>z: %{z:.2f} m<extra></extra>",
        "zsmooth": False,
        "name": "DEM",
    }]

    mesh_preview = (mesh or {}).get("preview") if mesh_type == "triangular" else None
    if mesh_preview and mesh_preview.get("tris"):
        xs = mesh_preview["x"]
        ys = mesh_preview["y"]
        ex: list = []
        ey: list = []
        for a, b, c in mesh_preview["tris"]:
            ex += [xs[a], xs[b], xs[c], xs[a], None]
            ey += [ys[a], ys[b], ys[c], ys[a], None]
        data2d.append({
            "type": "scattergl", "mode": "lines",
            "x": ex, "y": ey,
            "line": {"color": "rgba(255,107,53,0.75)", "width": 0.8},
            "hoverinfo": "skip",
            "name": f"Mesh ({mesh_preview['shown_triangles']:,} of "
                    f"{mesh_preview['total_triangles']:,} triangles)",
        })

    layout2d = {
        "xaxis": {"title": "X (projected CRS, m)",
                  "scaleanchor": "y", "scaleratio": 1},
        "yaxis": {"title": "Y (projected CRS, m)"},
        "height": 540,
        "margin": {"l": 80, "r": 30, "t": 20, "b": 50},
        "showlegend": True,
        "legend": {"orientation": "h", "y": 1.05, "x": 0,
                   "font": {"size": 11}},
    }
    fig2d = {"data": data2d, "layout": layout2d}

    # -- 3D surface + mesh-edge overlay -----------------------------------
    data3d = [{
        "type": "surface",
        "x": preview["x_vals"],
        "y": preview["y_vals"],
        "z": preview["z"],
        "colorscale": "earth",
        "colorbar": {"title": "Elevation (m)", "thickness": 14, "xpad": 0},
        "hovertemplate": "x: %{x:.0f}<br>y: %{y:.0f}<br>z: %{z:.2f} m<extra></extra>",
        "contours": {
            "z": {"show": False},
        },
        "lighting": {"ambient": 0.7, "diffuse": 0.6,
                     "specular": 0.15, "roughness": 0.6},
        "name": "DEM",
    }]

    if mesh_preview and mesh_preview.get("tris") and mesh_preview.get("z"):
        xs = mesh_preview["x"]
        ys = mesh_preview["y"]
        zs = mesh_preview["z"]
        # Bump the edge z a few metres above the surface so lines aren't
        # eaten by the surface polygons due to z-fighting / aliasing.
        elev_range = (max(zs) - min(zs)) if zs else 0.0
        offset = max(0.002 * elev_range, 0.5)
        ex: list = []
        ey: list = []
        ez: list = []
        for a, b, c in mesh_preview["tris"]:
            ex += [xs[a], xs[b], xs[c], xs[a], None]
            ey += [ys[a], ys[b], ys[c], ys[a], None]
            ez += [zs[a] + offset, zs[b] + offset, zs[c] + offset,
                   zs[a] + offset, None]
        data3d.append({
            "type": "scatter3d", "mode": "lines",
            "x": ex, "y": ey, "z": ez,
            "line": {"color": "rgba(255,107,53,0.9)", "width": 1.5},
            "hoverinfo": "skip",
            "name": f"Mesh ({mesh_preview['shown_triangles']:,} of "
                    f"{mesh_preview['total_triangles']:,} triangles)",
        })

    # Pick a z-aspect that makes relief visible even for floodplain DEMs.
    # Without vertical exaggeration a river valley with ~20 m of relief over
    # 2 km horizontal would look completely flat. Strategy:
    #   - compute the "natural" ratio dz / max(dx, dy)
    #   - if that's already >= 0.35 (rough / mountainous terrain), keep it
    #   - otherwise force z to occupy 35% of the horizontal footprint so
    #     the surface has legible relief
    # The annotated "vertical exaggeration: Nx" tells the user how stretched
    # the plot is vs. reality.
    # z-aspect: default to a generous vertical exaggeration so floodplain
    # DEMs (where Δz ≪ Δx) show visible relief. The user can tune this
    # live via the VE slider in the report. Also CLAMP the z-axis range to
    # the actual elevation range — otherwise Plotly extends the axis to
    # include 0 for surface traces, which squashes the terrain to the top.
    zs_all: list = []
    for row in preview["z"]:
        zs_all.extend(v for v in row if v is not None)
    x_vals = preview["x_vals"]; y_vals = preview["y_vals"]
    dx = (max(x_vals) - min(x_vals)) if x_vals else 1.0
    dy = (max(y_vals) - min(y_vals)) if y_vals else 1.0
    horiz = max(dx, dy) or 1.0
    if zs_all:
        z_min = min(zs_all); z_max = max(zs_all)
        dz = max(z_max - z_min, 1.0)
        natural = dz / horiz
        # Default VE: ~20× for flat floodplains, 1× for already-steep terrain.
        default_ve = max(1.0, 20.0 * min(1.0, 0.02 / max(natural, 1e-6)))
        z_aspect = max(natural, default_ve * natural)
        z_aspect = min(z_aspect, 2.0)  # cap
        ve_factor = (z_aspect * horiz) / dz
    else:
        z_min, z_max = 0.0, 1.0
        dz = 1.0
        z_aspect = 0.35
        ve_factor = 1.0

    # Padding on the z-range so the surface doesn't touch the top face
    pad = max(0.05 * dz, 0.5)

    layout3d = {
        "scene": {
            "xaxis": {"title": "X (m)"},
            "yaxis": {"title": "Y (m)"},
            "zaxis": {"title": "Elevation (m)",
                      "range": [z_min - pad, z_max + pad]},
            "aspectmode": "manual",
            "aspectratio": {"x": 1, "y": 1, "z": z_aspect},
            "camera": {"eye": {"x": 1.6, "y": -1.6, "z": 0.9}},
        },
        "height": 620,
        "margin": {"l": 0, "r": 0, "t": 30, "b": 0},
        "showlegend": True,
        "legend": {"orientation": "h", "y": 1.02, "x": 0,
                   "font": {"size": 11}},
    }
    fig3d = {
        "data": data3d,
        "layout": layout3d,
        # Plumbed through to the client so the VE slider can recompute aspect.z
        "meta": {
            "dz": float(dz),
            "horiz": float(horiz),
            "z_aspect_default": float(z_aspect),
            "ve_default": float(ve_factor),
        },
    }

    return {"fig2d": fig2d, "fig3d": fig3d}


def _map_section(data: dict) -> str:
    """HTML section: interactive DEM + mesh preview map (2D / 3D toggle)."""
    dem = data.get("dem") or {}
    mesh = data.get("mesh")
    mesh_type = data["config"]["mesh_type"]

    figs = _fig_dem_mesh_map(dem, mesh, mesh_type)
    if not figs:
        return ""

    preview = dem.get("preview") or {}
    stride = preview.get("stride", 1)
    note_parts = []
    if stride > 1:
        note_parts.append(f"DEM shown at stride {stride} "
                          f"({preview['rows']}×{preview['cols']} cells; "
                          f"keeps the embedded JSON light)")
    if mesh_type == "triangular" and mesh and (mesh.get("preview") or {}).get("stride", 1) > 1:
        mp = mesh["preview"]
        note_parts.append(f"mesh subsampled 1:{mp['stride']} "
                          f"({mp['shown_triangles']:,} of "
                          f"{mp['total_triangles']:,} triangles drawn)")
    note = (f'<p style="color:var(--text2);font-size:.82rem;margin-top:.6rem">'
            + "; ".join(note_parts) + "</p>") if note_parts else ""

    enc_2d = _esc(json.dumps(figs["fig2d"], default=float))
    enc_3d = _esc(json.dumps(figs["fig3d"], default=float))

    body_note = ("Orange lines show the computational mesh."
                 if mesh_type == "triangular"
                 else "Uses a regular Cartesian grid aligned with the DEM cells.")

    meta = figs["fig3d"].get("meta", {})
    ve_default = meta.get("ve_default", 1.0)
    # Expose dz / horiz so the slider JS can convert VE ↔ aspect.z
    dz = meta.get("dz", 1.0)
    horiz = meta.get("horiz", 1.0)
    # Reasonable slider bounds: 1× to 100× (covers flat floodplains up to
    # pretty mountainous terrain without looking silly).
    ve_max = 100
    ve_init = max(1, min(ve_max, round(ve_default)))

    return f"""
    <section class="section" id="dem-mesh-map">
      <h2>DEM &amp; mesh preview</h2>
      <div class="card">
        <div style="display:flex;align-items:center;justify-content:space-between;
                    gap:1rem;margin-bottom:.4rem;flex-wrap:wrap">
          <p style="color:var(--text2);margin:0;flex:1 1 380px">
            Interactive elevation map of the simulation domain.
            {body_note} Zoom / pan with the mouse; hover for the elevation value.
          </p>
          <div class="map-view-toggle" role="tablist" aria-label="Map view">
            <button type="button" class="view-btn active" data-view="2d"
                    aria-selected="true">🗺️ 2D</button>
            <button type="button" class="view-btn" data-view="3d"
                    aria-selected="false">🏔️ 3D</button>
          </div>
        </div>
        <div class="ve-slider-wrap" style="display:none">
          <label for="ve-slider">
            Vertical exaggeration
            <span class="ve-label">{ve_init}×</span>
            <span class="ve-note">
              (Δz = {dz:.1f} m over {horiz / 1000:.2f} km horizontal)
            </span>
          </label>
          <input id="ve-slider" type="range" min="1" max="{ve_max}"
                 step="1" value="{ve_init}"
                 data-dz="{dz:.6f}" data-horiz="{horiz:.6f}">
        </div>
        <div id="dem-mesh-plot-2d" class="plot map-view map-view-2d"
             data-fig="{enc_2d}"
             style="width:100%;min-height:540px"></div>
        <div id="dem-mesh-plot-3d" class="plot map-view map-view-3d"
             data-fig="{enc_3d}" data-lazy="1"
             data-dz="{dz:.6f}" data-horiz="{horiz:.6f}"
             style="width:100%;min-height:600px;display:none"></div>
        {note}
      </div>
    </section>
    """


def _domain_section(data: dict) -> str:
    dem = data.get("dem")
    if not dem:
        return """
        <section class="section" id="domain">
          <h2>Domain</h2>
          <div class="highlight-box error">DEM step failed — see errors below.</div>
        </section>
        """
    bbox = dem["bbox"]
    rows = []
    dl = dem.get("download")
    if dl:
        bbox_wgs = dl.get("bbox_wgs84") or [None] * 4
        vr = dl.get("vertical_rmse_m")
        rows.append(("Source", f'<strong>Auto-downloaded</strong> — {_esc(dl["label"])}'
                               f' ({_esc(dl["provider"])})'))
        if bbox_wgs[0] is not None:
            rows.append(("Requested bbox (WGS84 lon/lat)",
                         f'{bbox_wgs[0]:.4f}, {bbox_wgs[1]:.4f}, '
                         f'{bbox_wgs[2]:.4f}, {bbox_wgs[3]:.4f}'))
        rows.append(("Provider coverage", _esc(dl.get("coverage") or "—")))
        rows.append(("Native resolution", f'{dl["native_res_m"]} m'))
        if vr is not None:
            rows.append(("Reported vertical RMSE", f'{vr} m'))
        rows.append(("Downloaded GeoTIFF", f'<code>{_esc(dl["output_path"])}</code>'))
    else:
        rows.append(("Source GeoTIFF", f'<code>{_esc(dem["source"])}</code>'))
    rows.append(("Output .asc", f'<code>{_esc(dem["output_path"])}</code>'))
    rows.append(("CRS", _esc(dem["crs"])))
    rows.append(("Source resolution", f'{dem["source_resolution_m"]:.3f} m'))
    rows.append(("Target resolution", f'{dem["target_resolution_m"]:.3f} m'))
    rows.append(("Source grid", f'{dem["source_cols"]} × {dem["source_rows"]}'))
    rows.append(("Output grid", f'{dem["output_cols"]} × {dem["output_rows"]}'))
    rows.append(("Bounding box (xmin, ymin, xmax, ymax)",
                 f'{bbox[0]:.2f}, {bbox[1]:.2f}, {bbox[2]:.2f}, {bbox[3]:.2f}'))
    if dem.get("elev_min") is not None:
        rows.append(("Elevation range",
                     f'{dem["elev_min"]:.2f} – {dem["elev_max"]:.2f} m'))
    if dem.get("dummy_asc_for_trimesh"):
        rows.append(("Note",
                     "Triangular mesh selected — DEM .asc is a 10×10 placeholder; "
                     "real elevations live in the .msh vertex z-coordinates."))
    # Warn if target resolution is finer than native (interpolation ≠ real data)
    if dl and dem["target_resolution_m"] < dl["native_res_m"] * 0.9:
        rows.append(("⚠️ Resolution warning",
                     f'Target {dem["target_resolution_m"]:.1f} m is finer than '
                     f'the native {dl["native_res_m"]} m — this is interpolation, '
                     f'not added information. Use LiDAR-derived DEMs for sub-10 m '
                     f'urban flood modelling.'))

    tbl = ''.join(f'<tr><td>{_esc(k)}</td><td>{v}</td></tr>' for k, v in rows)
    return f"""
    <section class="section" id="domain">
      <h2>Domain</h2>
      <div class="card">
        <div class="table-wrap"><table>
          <tbody>{tbl}</tbody>
        </table></div>
      </div>
    </section>
    """


def _mesh_section(data: dict) -> str:
    mesh = data.get("mesh")
    config = data["config"]
    if config["mesh_type"] != "triangular":
        return """
        <section class="section" id="mesh">
          <h2>Mesh</h2>
          <div class="card"><p>Regular Cartesian grid — no Gmsh step.</p></div>
        </section>
        """
    if not mesh:
        return """
        <section class="section" id="mesh">
          <h2>Mesh</h2>
          <div class="highlight-box error">Mesh step failed — see errors below.</div>
        </section>
        """
    kpis = [
        _kpi("🔺", f'{mesh["cell_count"]:,}', "Triangles"),
        _kpi("📍", f'{mesh["vertex_count"]:,}', "Vertices"),
        _kpi("🧱", f'{mesh.get("boundary_edge_count", 0):,}', "Boundary edges"),
        _kpi("⏬", f'{mesh["min_size_m"]:.1f} m', "Min edge"),
        _kpi("⏫", f'{mesh["max_size_m"]:.1f} m', "Max edge"),
        _kpi("📉", f'{mesh.get("max_slope", 0):.3f}', "Max slope (|grad z|)"),
    ]
    return f"""
    <section class="section" id="mesh">
      <h2>Mesh</h2>
      <div class="kpi-grid">{''.join(kpis)}</div>
      <div class="card">
        <p>Adaptive triangular mesh written to
        <code>{_esc(mesh['output_path'])}</code>. Refinement scales with local
        slope (factor = {mesh['slope_factor']:.2f}).</p>
      </div>
    </section>
    """


def _config_section(data: dict) -> str:
    out = data.get("config_out")
    if not out:
        return """
        <section class="section" id="config">
          <h2>Configuration</h2>
          <div class="highlight-box error">modset.json was not generated — see errors below.</div>
        </section>
        """
    modset = out["modset"]
    pretty = json.dumps(modset, indent=2)
    rows = [
        ("Output path", f'<code>{_esc(out["output_path"])}</code>'),
        ("DEM file", f'<code>{_esc(modset["DEM_FILE"])}</code>'),
        ("Mesh type", _esc(modset.get("MESH_TYPE", "regular"))),
    ]
    if modset.get("MESH_FILE"):
        rows.append(("Mesh file", f'<code>{_esc(modset["MESH_FILE"])}</code>'))
    if modset.get("METEO_FILE"):
        rows.append(("Meteo file", f'<code>{_esc(modset["METEO_FILE"])}</code>'))
    if modset.get("INFLOW_FILE"):
        rows.append(("Inflow file",
                     f'<code>{_esc(modset["INFLOW_FILE"]["FILENAME"])}</code>'))
    rows.append(("Simulation start", _esc(modset["SIM_DATETIME_START"])))
    rows.append(("Roughness", f'{modset["ROUGNESS_HEIGHT"]} m'))
    rows.append(("Output folder", _esc(modset["OUTPUT"]["OUTPUT_FOLDER"])))
    rows.append(("Print step", f'{modset["OUTPUT"]["PRINT_STEP"]} s'))
    rows.append(("h_min to print", f'{modset["OUTPUT"]["H_MIN_TO_PRINT"]} m'))

    tbl = ''.join(f'<tr><td>{_esc(k)}</td><td>{v}</td></tr>' for k, v in rows)
    return f"""
    <section class="section" id="config">
      <h2>Configuration</h2>
      <div class="card">
        <div class="table-wrap"><table><tbody>{tbl}</tbody></table></div>
        <details class="module-details" style="margin-top:1rem">
          <summary>View full modset.json</summary>
          <div class="module-content">{_code_block(pretty)}</div>
        </details>
      </div>
    </section>
    """


def _modules_section(data: dict) -> str:
    config = data["config"]

    def row(name: str, enabled: bool, params: dict | None) -> str:
        if enabled:
            pill = _badge("secondary", "Enabled")
            body = ('<div class="module-content">'
                    + ''.join(f'<p><strong>{_esc(k)}</strong>: '
                              f'<code>{_esc(v)}</code></p>'
                              for k, v in (params or {}).items())
                    + '</div>') if params else ''
            summary = f'<summary>{_esc(name)} {pill}</summary>'
            return f'<details class="module-details" open>{summary}{body}</details>'
        pill = _badge("none", "Disabled")
        return (f'<details class="module-details"><summary>{_esc(name)} {pill}'
                f'</summary></details>')

    ade = config["ade_transport"]
    soil = config["soil_infiltration"]
    items = [
        row("ADE Transport (advection–dispersion)", ade["enabled"],
            {"D_COEF": ade.get("d_coef")} if ade["enabled"] else None),
        row("Soil Infiltration", soil["enabled"],
            {"DEFAULT_KS_MM_HR": soil.get("default_ks_mm_hr")}
            if soil["enabled"] else None),
    ]
    return f"""
    <section class="section" id="modules">
      <h2>Optional Modules</h2>
      <div class="card">{''.join(items)}</div>
    </section>
    """


def _errors_section(data: dict) -> str:
    errs = data.get("errors") or []
    if not errs:
        return ""
    items = []
    for step, msg, tb in errs:
        items.append(
            '<div class="highlight-box error">'
            f'<strong>{_esc(step)} step failed:</strong> {_esc(msg)}'
            + (f'<details style="margin-top:.5rem"><summary>Traceback</summary>'
               f'<pre style="margin-top:.5rem;font-size:.75rem;white-space:pre-wrap">'
               f'{_esc(tb)}</pre></details>' if tb else '')
            + '</div>'
        )
    return f"""
    <section class="section" id="errors">
      <h2>Errors</h2>
      {''.join(items)}
    </section>
    """


def _next_steps_section(data: dict) -> str:
    config = data["config"]
    repo_root = data["repo_root"]
    modset_out = data.get("config_out")
    os_key, os_label, os_color = _detect_os()

    if not modset_out:
        return f"""
        <section class="section" id="nextsteps">
          <h2>Next Steps</h2>
          <div class="highlight-box warning">
            modset.json was not produced — fix the errors above and re-run.
          </div>
        </section>
        """

    modset_rel = os.path.relpath(modset_out["output_path"], repo_root).replace(os.sep, "/")
    results_dir = os.path.join(repo_root, "Results")

    cd_repo = f'cd "{repo_root}"'

    # Step 2: build the container image (installs deps, ships source).
    # The Dockerfile no longer compiles FLUXOS — that's step 4, in the shell.
    build_prefix = "USE_MPI=ON " if config.get("use_mpi") else ""
    build_cmd = (f'{build_prefix}'
                 f'docker compose -f containers/docker-compose.yml build')

    # Step 3: open a shell inside the container
    shell_cmd = 'docker compose -f containers/docker-compose.yml run --rm fluxos'

    # Step 4: inside the shell, compile FLUXOS and run the simulation.
    # -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/work/bin drops the compiled binary
    # onto the host's ./bin/ via the docker-compose bind mount.
    cmake_flags = ["-DMODE_release=ON"]
    if config.get("use_trimesh_build"):
        cmake_flags.append("-DUSE_TRIMESH=ON")
    if config.get("use_mpi"):
        cmake_flags.append("-DUSE_MPI=ON")
    cmake_flags.append("-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/work/bin")
    binary_name = "fluxos_mpi" if config.get("use_mpi") else "fluxos"
    run_line = (
        f'mpirun -n {int(config["mpi_np"])} ./bin/{binary_name} {modset_rel}'
        if config.get("use_mpi")
        else f'./bin/{binary_name} {modset_rel}'
    )
    compile_cmd = (
        "# You are now inside the container\n"
        "cd /opt/fluxos\n"
        "mkdir -p build && cd build\n"
        f"cmake {' '.join(cmake_flags)} ..\n"
        "make -j$(nproc)\n"
        "\n"
        "# Run the model — cd to /work so the modset's relative paths\n"
        "# resolve against the bind-mounted Working_example/ directory\n"
        "cd /work\n"
        f"{run_line}"
    )

    check_cmd = f'ls -la "{results_dir}"\n{_open_outputs_cmd(os_key, results_dir)}'

    viewer_script = os.path.join(repo_root, "supporting_scripts",
                                 "2_Read_Outputs", "output_supporting_lib",
                                 "fluxos_viewer.py")
    viz_base = (f'python "{viewer_script}" \\\n'
                f'    --results-dir "{results_dir}" \\\n'
                f'    --dem "{os.path.join(repo_root, modset_out["modset"]["DEM_FILE"])}" \\\n'
                f'    --utm-zone {int(config.get("dem_utm_zone", 10))}')
    viz_with_transport = viz_base + ' \\\n    --variable conc_SW'

    ade_enabled = bool(config["ade_transport"].get("enabled"))
    checked_attr = " checked" if ade_enabled else ""
    base_display = "none" if ade_enabled else "block"
    transport_display = "block" if ade_enabled else "none"

    # Statistics report snippet — runs the 2_Read_Outputs template.
    stats_template = os.path.join(repo_root, "supporting_scripts",
                                  "2_Read_Outputs", "read_output_template.py")
    stats_cmd = (
        f'# Edit results_dir / modset_file at the top of the file, then:\n'
        f'python "{stats_template}"'
    )

    return f"""
    <section class="section" id="nextsteps">
      <h2>Next Steps</h2>
      <div class="card primary">
        <p style="margin-bottom:.4rem">
          <span class="os-badge" style="background:{os_color}">💻 {_esc(os_label)}</span>
          Copy-paste these snippets into your terminal. The container ships the
          build dependencies and the source tree &mdash; FLUXOS is compiled
          inside the container shell in step 4.
        </p>

        <h3>1. Move into the FLUXOS repo</h3>
        {_code_block(cd_repo)}

        <h3>2. Build the container image <span style="font-weight:400;color:var(--text2)">(deps only)</span></h3>
        {_code_block(build_cmd)}

        <h3>3. Open a shell inside the container</h3>
        {_code_block(shell_cmd)}

        <h3>4. Compile FLUXOS and run the simulation <span style="font-weight:400;color:var(--text2)">(inside the shell)</span></h3>
        {_code_block(compile_cmd)}

        <h3>5. Check outputs <span style="font-weight:400;color:var(--text2)">(back on the host)</span></h3>
        {_code_block(check_cmd)}

        <h3>6. Visualise results</h3>
        <p>Two complementary post-processing tools are provided in
        <code>2_Read_Outputs/</code>.</p>

        <h4 style="margin-top:1rem;font-size:.95rem">6a. Interactive WebGL animation (KML / MP4 / WebGL)</h4>
        <div class="viz-toggle-card">
          <label class="viz-toggle-label">
            <input type="checkbox" id="viz-conc-toggle"{checked_attr}>
            Include chemical transport in the animation
            (adds <code>--variable conc_SW</code>)
          </label>
          <div class="viz-cmd-base" style="display:{base_display}">{_code_block(viz_base)}</div>
          <div class="viz-cmd-transport" style="display:{transport_display}">{_code_block(viz_with_transport)}</div>
        </div>

        <h4 style="margin-top:1.5rem;font-size:.95rem">6b. Statistics report
          <span style="font-weight:400;color:var(--text2)">
            — flood volume / flooded-area time series, max-inundation map,
            hazard classification (ARR-2019 H·V), depth histogram, and first-inundation map
          </span>
        </h4>
        {_code_block(stats_cmd)}
      </div>
    </section>
    """


def _footer(generated_at: str) -> str:
    return f"""
    <footer class="footer">
      <div class="logo">FLUX<span>OS</span></div>
      Report generated {_esc(generated_at)} —
      <a href="https://github.com/ue-hydro/FLUXOS_cpp">github.com/ue-hydro/FLUXOS_cpp</a>
    </footer>
    """


# ---------------------------------------------------------------------------
#  Main entry
# ---------------------------------------------------------------------------

def generate_report(report_data: dict, output_path: str) -> str:
    """
    Build a self-contained HTML report at `output_path`.

    Returns the output path (so the caller can log or open it).
    """
    config = report_data["config"]
    generated_at = report_data.get("generated_at") \
        or datetime.now().isoformat(timespec="seconds")

    nav_items = [
        ("summary", "Summary"),
        ("domain", "Domain"),
    ]
    # Only show the DEM & mesh map if we captured preview data for it.
    if (report_data.get("dem") or {}).get("preview"):
        nav_items.append(("dem-mesh-map", "DEM & mesh map"))
    nav_items += [
        ("mesh", "Mesh"),
        ("config", "Configuration"),
        ("modules", "Modules"),
    ]
    if report_data.get("errors"):
        nav_items.append(("errors", "Errors"))
    nav_items.append(("nextsteps", "Next Steps"))

    nav_html = ''.join(f'<a href="#{sid}">{_esc(label)}</a>'
                       for sid, label in nav_items)

    html_doc = f"""<!DOCTYPE html>
<html lang="en" data-theme="dark">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>FLUXOS report — {_esc(config['project_name'])}</title>
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
      {_header(config, generated_at)}
      <div class="container">
        <section class="section" id="project">
          <h2>Project</h2>
          <div class="card primary">
            <h3>{_esc(config['project_name'])}</h3>
            <p>{_esc(config.get('description') or '')}</p>
          </div>
        </section>
        {_summary_section(report_data)}
        {_domain_section(report_data)}
        {_map_section(report_data)}
        {_mesh_section(report_data)}
        {_config_section(report_data)}
        {_modules_section(report_data)}
        {_errors_section(report_data)}
        {_next_steps_section(report_data)}
      </div>
      {_footer(generated_at)}
    </main>
  </div>
  <script>{_JS}</script>
</body>
</html>
"""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_doc)
    return output_path
