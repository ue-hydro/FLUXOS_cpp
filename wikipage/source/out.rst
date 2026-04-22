Output
==================================

FLUXOS writes one file per print step into ``Results/`` — ASCII ``.txt`` for
regular meshes, VTK ``.vtu`` for triangular meshes (plus a ``.pvd`` time
series index). For analysis and visualisation of these files, you almost
never need to parse them by hand: the ``read_output_template.py`` template
(see :doc:`SupportingScripts`) streams them into an HTML flood-statistics
report, and ``fluxos_viewer.py`` converts them to animated KML / WebGL / MP4.

The page below documents the on-disk column layout for when you do need to
consume the files directly.

.. toctree::
   :maxdepth: 2

   file

