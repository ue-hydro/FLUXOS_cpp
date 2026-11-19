// Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
// This file is part of the FLUXOS model.

// This program, FLUXOS, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Time-varying water-level (stage) boundary — see water_level_bc.h.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "water_level_bc.h"

namespace {

std::string lower(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// Read "t,level" pairs; non-numeric lines (header, blanks, comments) skipped.
bool load_series(const std::string& fname, std::vector<double>& t,
                 std::vector<double>& lv, std::string& msg) {
    std::ifstream f(fname);
    if (!f.is_open()) {
        msg = "cannot open WATER_LEVEL file: " + fname;
        return true;
    }
    std::string line;
    while (std::getline(f, line)) {
        for (auto& c : line) if (c == ',' || c == ';' || c == '\t') c = ' ';
        std::istringstream ss(line);
        double a, b;
        if (ss >> a >> b) {
            if (!t.empty() && a <= t.back()) {
                msg = "WATER_LEVEL times must be strictly increasing (t=" +
                      std::to_string(a) + ") in " + fname;
                return true;
            }
            t.push_back(a);
            lv.push_back(b);
        }
    }
    if (t.empty()) {
        msg = "no numeric 't_seconds,level_m' rows in WATER_LEVEL file: " + fname;
        return true;
    }
    return false;
}

struct Segment {
    std::string edge;
    long from = -1, to = -1;   // 1-based, inclusive; -1 = whole edge
};

}  // namespace

bool read_water_level_bc(GlobVar& ds, std::ofstream& logFLUXOSfile)
{
    ds.wl_bc_enabled = false;

    auto bc_it = ds.master_MODSET.find("BOUNDARY_CONDITIONS");
    if (bc_it == ds.master_MODSET.end() || !bc_it->is_object())
        return false;
    auto wl_it = bc_it->find("WATER_LEVEL");
    if (wl_it == bc_it->end())
        return false;

    std::string msg;
    auto fail = [&](const std::string& m) {
        std::string full = "ERROR (BOUNDARY_CONDITIONS.WATER_LEVEL): " + m;
        std::cout << full << std::endl;
        logFLUXOSfile << full << "\n";
        return true;
    };

    const nlohmann::json& wl = *wl_it;
    if (!wl.is_object())
        return fail("must be a JSON object");

    // Reject unknown keys so typos do not silently disable the BC.
    for (auto it = wl.begin(); it != wl.end(); ++it) {
        const std::string& k = it.key();
        if (k != "STATUS" && k != "FILENAME" && k != "EDGE" &&
            k != "VELOCITY" && k != "COMMENT")
            return fail("unknown key '" + k + "' (allowed: STATUS, FILENAME, EDGE, VELOCITY, COMMENT)");
    }

    bool status = false;
    try {
        status = wl.value("STATUS", false);
    } catch (...) {
        return fail("STATUS must be true/false");
    }
    if (!status) {
        logFLUXOSfile << "Water-level boundary: DISABLED\n";
        return false;
    }

    if (ds.mesh_type == "triangular")
        return fail("only implemented for the regular mesh (MESH_TYPE \"regular\")");
#ifdef USE_CUDA
    return fail("not implemented in the CUDA build — rebuild without USE_CUDA "
                "(OpenMP) to use a water-level boundary");
#endif

    // ---- options ---------------------------------------------------------
    try {
        ds.wl_bc_file = wl.at("FILENAME").get<std::string>();
    } catch (...) {
        return fail("FILENAME (string) is required");
    }
    if (wl.contains("VELOCITY")) {
        std::string v;
        try { v = lower(wl["VELOCITY"].get<std::string>()); }
        catch (...) { return fail("VELOCITY must be a string"); }
        if (v == "zero") ds.wl_bc_zero_velocity = true;
        else if (v == "zero_gradient") ds.wl_bc_zero_velocity = false;
        else return fail("VELOCITY must be \"zero_gradient\" or \"zero\"");
    }

    // ---- segments --------------------------------------------------------
    if (!wl.contains("EDGE"))
        return fail("EDGE is required");
    std::vector<Segment> segs;
    const nlohmann::json& je = wl["EDGE"];
    const double dxy = static_cast<double>(ds.dxy);
    auto parse_one = [&](const nlohmann::json& j, Segment& s) -> bool {
        if (j.is_string()) {
            s.edge = lower(j.get<std::string>());
            return false;
        }
        if (!j.is_object() || !j.contains("EDGE") || !j["EDGE"].is_string()) {
            msg = "each EDGE entry must be a string or {\"EDGE\": ..., \"FROM\"/\"TO\" | \"FROM_COORD\"/\"TO_COORD\"}";
            return true;
        }
        s.edge = lower(j["EDGE"].get<std::string>());
        for (auto it = j.begin(); it != j.end(); ++it) {
            const std::string& k = it.key();
            if (k != "EDGE" && k != "FROM" && k != "TO" &&
                k != "FROM_COORD" && k != "TO_COORD") {
                msg = "unknown key '" + k + "' in EDGE segment";
                return true;
            }
        }
        const bool along_rows = (s.edge == "west" || s.edge == "east");
        const double ll = along_rows ? ds.YLLCORNER : ds.XLLCORNER;
        try {
            if (j.contains("FROM_COORD") || j.contains("TO_COORD")) {
                double c0 = j.at("FROM_COORD").get<double>();
                double c1 = j.at("TO_COORD").get<double>();
                if (c1 < c0) std::swap(c0, c1);
                s.from = static_cast<long>(std::floor((c0 - ll) / dxy)) + 1;
                s.to   = static_cast<long>(std::floor((c1 - ll) / dxy)) + 1;
            } else if (j.contains("FROM") || j.contains("TO")) {
                s.from = j.at("FROM").get<long>();
                s.to   = j.at("TO").get<long>();
                if (s.to < s.from) std::swap(s.from, s.to);
            }
        } catch (...) {
            msg = "segment needs both FROM and TO (cell index) or both FROM_COORD and TO_COORD (numbers)";
            return true;
        }
        return false;
    };
    if (je.is_array()) {
        for (const auto& j : je) {
            Segment s;
            if (parse_one(j, s)) return fail(msg);
            segs.push_back(s);
        }
    } else {
        Segment s;
        if (parse_one(je, s)) return fail(msg);
        segs.push_back(s);
    }
    if (segs.empty())
        return fail("EDGE list is empty");

    // ---- build cell list ---------------------------------------------------
    const long NR = static_cast<long>(ds.NROWS), NC = static_cast<long>(ds.NCOLS);
    ds.wl_bc_mask = std::unique_ptr<arma::Mat<float>>(
        new arma::fmat(ds.MROWS, ds.MCOLS, arma::fill::zeros));
    arma::Mat<float>& mask = *ds.wl_bc_mask;
    const arma::Mat<float>& nodata = *ds.innerNeumannBCWeir;
    std::ostringstream summary;
    size_t n_skipped_nodata = 0;

    for (const auto& s : segs) {
        long n_along;
        if (s.edge == "west" || s.edge == "east") n_along = NR;
        else if (s.edge == "south" || s.edge == "north") n_along = NC;
        else return fail("EDGE must be west/east/south/north, got '" + s.edge + "'");
        long a = (s.from < 0) ? 1 : s.from;
        long b = (s.to   < 0) ? n_along : s.to;
        a = std::max(a, 1L);
        b = std::min(b, n_along);
        if (a > b)
            return fail("segment on '" + s.edge + "' lies outside the grid (1.." +
                        std::to_string(n_along) + ")");
        summary << " " << s.edge << "[" << a << ".." << b << "]";
        for (long k = a; k <= b; k++) {
            unsigned int ir, ic, ir_in, ic_in;
            if (s.edge == "west")       { ir = k; ic = 1;  ir_in = k;      ic_in = 2; }
            else if (s.edge == "east")  { ir = k; ic = NC; ir_in = k;      ic_in = NC - 1; }
            else if (s.edge == "south") { ir = 1; ic = k;  ir_in = 2;      ic_in = k; }
            else                        { ir = NR; ic = k; ir_in = NR - 1; ic_in = k; }
            if (nodata(ir, ic) != 0.0f) { n_skipped_nodata++; continue; }
            if (mask(ir, ic) != 0.0f) continue;   // corner shared by two edges
            mask(ir, ic) = 1.0f;
            ds.wl_bc_irow.push_back(ir);
            ds.wl_bc_icol.push_back(ic);
            ds.wl_bc_irow_in.push_back(ir_in);
            ds.wl_bc_icol_in.push_back(ic_in);
        }
    }
    if (ds.wl_bc_irow.empty())
        return fail("no valid (non-NODATA) cells on the requested EDGE");

    // ---- series ------------------------------------------------------------
    if (load_series(ds.wl_bc_file, ds.wl_bc_t, ds.wl_bc_level, msg))
        return fail(msg);

    ds.wl_bc_enabled = true;
    std::ostringstream info;
    info << "Water-level boundary: ENABLED | file = " << ds.wl_bc_file
         << " | " << ds.wl_bc_t.size() << " records, t = " << ds.wl_bc_t.front()
         << ".." << ds.wl_bc_t.back() << " s, level = "
         << *std::min_element(ds.wl_bc_level.begin(), ds.wl_bc_level.end()) << ".."
         << *std::max_element(ds.wl_bc_level.begin(), ds.wl_bc_level.end())
         << " m | cells = " << ds.wl_bc_irow.size() << " on" << summary.str()
         << (n_skipped_nodata ? " (" + std::to_string(n_skipped_nodata) + " NODATA cells skipped)" : "")
         << " | velocity = " << (ds.wl_bc_zero_velocity ? "zero" : "zero_gradient");
    std::cout << info.str() << std::endl;
    logFLUXOSfile << info.str() << "\n";
    return false;
}

double water_level_at(const GlobVar& ds, double t)
{
    const auto& T = ds.wl_bc_t;
    const auto& L = ds.wl_bc_level;
    if (t <= T.front()) return L.front();
    if (t >= T.back())  return L.back();
    size_t k = std::upper_bound(T.begin(), T.end(), t) - T.begin();  // T[k-1] <= t < T[k]
    const double w = (t - T[k - 1]) / (T[k] - T[k - 1]);
    return L[k - 1] + w * (L[k] - L[k - 1]);
}

void apply_water_level_bc(GlobVar& ds, int nchem)
{
    if (!ds.wl_bc_enabled) return;

    const double level = water_level_at(ds, ds.tim);
    ds.wl_bc_level_now = level;

    arma::Mat<double>& z = *ds.z;
    arma::Mat<double>& zb = *ds.zb;
    arma::Mat<double>& h = *ds.h;
    arma::Mat<double>& h0 = *ds.h0;
    arma::Mat<double>& qx = *ds.qx;
    arma::Mat<double>& qy = *ds.qy;
    arma::Mat<float>& ldry = *ds.ldry;
    const double hdry = ds.hdry;
    const double area = static_cast<double>(ds.arbase);

    double vin = 0.0, vout = 0.0;
    const size_t n = ds.wl_bc_irow.size();
    for (size_t k = 0; k < n; k++) {
        const unsigned int ir = ds.wl_bc_irow[k], ic = ds.wl_bc_icol[k];
        const unsigned int jr = ds.wl_bc_irow_in[k], jc = ds.wl_bc_icol_in[k];

        const double zbp = zb(ir, ic);
        const double h_old = std::fmax(0.0, z(ir, ic) - zbp);
        const double h_new = std::fmax(0.0, level - zbp);
        const double dv = (h_new - h_old) * area;
        if (dv > 0.0) vin += dv; else vout -= dv;

        z(ir, ic) = zbp + h_new;
        h(ir, ic) = h_new;
        h0(ir, ic) = h_new;   // forcing-style: ADE sees no artificial dilution

        // Momentum: zero-gradient velocity from the first interior cell.
        double qxn = 0.0, qyn = 0.0;
        if (!ds.wl_bc_zero_velocity && h_new > hdry) {
            const double h_in = std::fmax(0.0, z(jr, jc) - zb(jr, jc));
            if (h_in > hdry) {
                qxn = qx(jr, jc) / h_in * h_new;
                qyn = qy(jr, jc) / h_in * h_new;
            }
        }
        qx(ir, ic) = qxn;
        qy(ir, ic) = qyn;
        ldry(ir, ic) = (h_new > hdry) ? 0.0f : 1.0f;

        // ADE: zero-gradient concentration (no external solute supplied).
        if (ds.ade_solver) {
            for (int ichem = 0; ichem < nchem; ichem++)
                (*ds.conc_SW)[ichem](ir, ic) = (*ds.conc_SW)[ichem](jr, jc);
        }
    }
    ds.wl_bc_vol_in += vin;
    ds.wl_bc_vol_out += vout;
    ds.wl_bc_vol_cum += vin - vout;
}

void write_water_level_bc_balance(const GlobVar& ds, double t_print)
{
    if (!ds.wl_bc_enabled) return;
    static bool first = true;
    const std::string fname = ds.output_folder + "/water_level_bc_balance.csv";
    std::ofstream f(fname, first ? std::ios::trunc : std::ios::app);
    if (!f.is_open()) return;
    if (first) {
        f << "t_s,level_m,vol_in_m3,vol_out_m3,vol_net_m3\n";
        first = false;
    }
    f.precision(10);
    f << t_print << "," << ds.wl_bc_level_now << "," << ds.wl_bc_vol_in << ","
      << ds.wl_bc_vol_out << "," << ds.wl_bc_vol_cum << "\n";
}
