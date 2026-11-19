# CyFj11 0.99.114

* Support for FlowJo 11.3 workspaces (analysis schema 4.0.0). Gates on
  parameters derived by FlowJo algorithms (event quality filter, FlowSOM,
  t-SNE, UMAP) no longer stop `fj11_to_gatingset()` with
  `colname not found: QC of Experiment Data`. The per-event values stored as
  CSV files inside the `.flowjo` archive are appended to each sample's
  cytoframe before gating. On the FlowJo 11.3 reference workspace all 30
  population counts agree with FlowJo's own to within 44 events.
* The analysis `schemaVersion` is now parsed once (`ws$schema_version`, see
  `fj11_schema_version()`) and passed to every conversion step, so version
  specific behaviour can be added where needed. `read_flowjo11_workspace()`
  warns when a workspace uses a schema version newer than those tested
  (3.0.0 and 4.0.0).
* Fixed ellipse gates for schema 4.0.0 (FlowJo 11.3) and newer: the
  covariance matrix swapped the semi-axes and reversed the rotation, so ellipse
  gates had too few events (113,018 instead of about 140,900 on the 11.3
  reference workspace). Workspaces with an older or unknown schema version keep
  the previous behaviour so that earlier results stay reproducible.
* Known issue: `rotationAngle` of ellipse gates in the older bundled
  `test.data.flowjo` (gateResolution 1024) still gives 4,671 events instead of
  FlowJo's 4,523.
* Known issue: acquisition-defined `$SPILLOVER` matrices are not picked up when
  the workspace stores no compensation ("No compensation found").

# CyFj11 0.1.0

* Initial release of the package
* Added functions to parse FlowJo v11 workspace files
* Implemented gate extraction functionality
* Created population extraction methods
* Added basic validation functions for gate consistency
* Included FCS file discovery and resolution utilities
* Provided export functionality for extracted data
* Fixed issue with handling multiple FCS file matches in fj11_to_gatingset() function
* Improved error handling for stop_on_multiple parameter
* Fixed "Error: 1 not found!" issue in compute_boolean_gates function
* Added support for additional gate types: spiderQuad, range, curlyQuad, quad
* Improved boolean gate conversion robustness with better error handling
* Enhanced validation for population paths in boolean gate computation