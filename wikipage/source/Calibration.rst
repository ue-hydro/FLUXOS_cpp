Model Calibration
==================================

.. note::

   This section is under development.

Overview
--------

Model calibration is the process of adjusting model parameters to minimize the discrepancy between simulated and observed results. For FLUXOS, the primary calibration parameter is the roughness height (Manning's roughness), which controls the friction forces in the shallow water equations.

Calibration Parameters
-----------------------

The following parameters can be adjusted during calibration:

.. list-table::
   :widths: 30 20 50
   :header-rows: 1

   * - Parameter
     - Config Key
     - Description
   * - Roughness height
     - ``ROUGNESS_HEIGHT``
     - Manning's roughness coefficient controlling bed friction [m]
   * - Minimum water depth
     - ``H_MIN_TO_PRINT``
     - Threshold for wet/dry classification and output filtering [m]

Calibration Workflow
---------------------

A typical calibration workflow for FLUXOS involves:

1. **Define objective function**: Select metrics to compare simulated vs. observed data (e.g., Nash-Sutcliffe Efficiency, RMSE, KGE)
2. **Select observation data**: Identify available field measurements (water levels, discharge at cross-sections, inundation extent)
3. **Parameter sampling**: Define parameter ranges and sampling strategy (manual, grid search, or automated optimization)
4. **Batch simulation**: Run FLUXOS with different parameter sets
5. **Post-processing**: Extract simulated values at observation locations using ``cross_section_extract.py`` or custom scripts
6. **Evaluation**: Compare simulated vs. observed data using the objective function
7. **Iteration**: Refine parameter ranges and repeat

The supporting scripts in ``fluxos_supporting_scripts/`` provide tools for steps 5-6, including cross-section extraction (``cross_section_extract.py``), graphing functions for simulated vs. observed comparison (``graphing_functions.py``), and batch result analysis (``Analyze_Results.m``).
