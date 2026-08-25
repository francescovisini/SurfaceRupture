# Backbone Surface Rupture (BSR)

**Current release:** MATLAB v1.0.1 · QGIS plugin v1.0.1

**Reference:** Visini (2026), manuscript under review.

The **Backbone Surface Rupture (BSR)** software reconstructs one continuous, ordered representation of a principal surface rupture from mapped **principal-fault traces**.

The method is not restricted to a specific database. **SURE 2.0 rank-1 traces are the reference application used in Visini (2026)**, but the QGIS implementation can operate on any polyline dataset in which the principal rupture traces can be identified.

Two implementations are provided:

- **MATLAB reference implementation (v1.0.1)** — the research implementation used for method development, validation, diagnostics, and manuscript results.
- **QGIS Processing tool (v1.0.1)** — a user-facing implementation that runs directly on a polyline layer.

The BSR is not a replacement for the complete rupture map. It is a **derived geometric object** that follows observed principal-fault traces whenever possible and uses explicit straight artificial connections only where required to preserve continuity.

> **QGIS users:** after installing the plugin, open  
> **Processing Toolbox → Backbone Surface Rupture → Build Backbone Surface Rupture**.  
> The Getting Started dialog shown after installation can also open the Processing Toolbox directly.

---

## 1. What the software does

In simple terms, the software:

1. reads the principal-fault lines belonging to an earthquake;
2. keeps their original mapped geometry unchanged;
3. constructs a graph along the observed lines;
4. allows controlled connections between nearby lines;
5. tests plausible paths between opposite ends of the rupture;
6. selects the path that first represents the global rupture extent and then maximizes support from the observed principal-fault geometry;
7. calculates cumulative distance along the BSR and normalized position `x/L`;
8. exports the BSR, the artificial connections, ordered vertices, and reconstruction metrics.

The final BSR contains two geometric components:

- **TRACE** — a BSR segment that follows an observed principal-fault polyline and preserves its original geometry;
- **GAP** — an explicit artificial straight connection introduced between observed principal-fault traces.

Internal graph junctions are only topological objects. They do not modify the original mapped polylines and are converted to explicit GAP geometry in the final output whenever they span a finite distance.

---

## 2. Which implementation should I use?

### Use MATLAB when:

- you want the complete reference workflow used for the manuscript;
- you need graph, candidate-path, and diagnostic outputs;
- you want to reproduce the SURE 2.0 paper calculations;
- you want to inspect or modify the research code.

### Use QGIS when:

- you want to select a layer and click **Run**;
- your dataset is not necessarily SURE 2.0;
- you do not want to prepare external event lists;
- you want the BSR and `x/L` results directly in a GIS project;
- you want to use the algorithm from the Processing Toolbox, Model Designer, or batch interface.

---

# 3. Input data

## 3.1 Minimum geometric requirements

The input must be a **polyline layer**. Each feature must represent a mapped portion of a surface rupture.

For the QGIS tool, the input requires:

| Information | Example | Required? | Meaning |
|---|---|---:|---|
| Event identifier field | `IdE`, `EventID`, `eq_id` | Yes | Identifies which earthquake each line belongs to. |
| Principal-fault classification field | `Comp_rank`, `Class`, `Type`, `MainFault` | No | Used only when the input contains both principal and secondary/distributed traces. |
| Principal-fault value | `1`, `Principal`, `Main` | Only if a classification field is selected | Value identifying the traces that may form the BSR. |
| Line geometry | — | Yes | The mapped surface-rupture polylines. |

The QGIS field names are not fixed: the user selects the appropriate column from the interface.

The MATLAB reference implementation now supports the same general input logic. The event field, principal-fault classification field, and principal-fault value are defined in the user settings. For SURE 2.0, the reference settings are `IdE` and `Comp_rank = 1`; for other datasets, these can be replaced by the corresponding field and value (for example, `Class = Main`). If the input already contains only principal-fault traces, the classification filter can be disabled.

## 3.2 Selecting the principal-fault traces in QGIS

The QGIS tool supports two input workflows.

### Case A — your layer already contains only principal-fault traces

Leave:

```text
Principal-fault classification field
```

**blank**.

The tool will use every line in the input layer. You do **not** need to create a ranking column.

Example: you have already extracted only the principal ruptures of one earthquake into a separate shapefile or GeoPackage layer.

### Case B — your layer contains principal and secondary/distributed traces

Choose the attribute column that identifies the principal traces and specify the corresponding value.

Examples:

```text
Comp_rank = 1
Class = Principal
Type = Main
```

Numeric and text classifications are both accepted.

### Which option should I use?

| My dataset | What to do |
|---|---|
| It contains only principal/main rupture traces | Leave the classification field blank. |
| It contains principal and secondary/distributed ruptures | Select the classification column and the value identifying the principal traces. |
| It is SURE 2.0 | Select `Comp_rank` and use value `1`. |

The classification filter changes only **which input features are passed to the BSR algorithm**. It does not change the BSR algorithm itself.

## 3.3 Important: treatment of secondary/distributed traces

In the current implementation, secondary or distributed ruptures are **not used to build or weight the BSR**.

They:

- cannot become TRACE portions of the BSR;
- are not used to reduce the cost of an artificial connection;
- are not used as intermediate stepping stones across a gap.

Artificial connections are evaluated using the selected principal-fault network only.

This is intentional in the present BSR definition. A possible future extension is to use secondary/distributed traces as auxiliary spatial evidence when evaluating artificial connections, without allowing them to become part of the final backbone. Such an extension would change the scientific formulation and is not included in the current reference algorithm.

## 3.4 Coordinate reference system

A geographic layer such as EPSG:4326 is accepted. For each event, both implementations convert coordinates to a local metric system centered on the event before calculating distances and orientations.

The reference conversion is:

```text
x = R * delta_longitude * cos(latitude_0)
y = R * delta_latitude
```

where `R = 6371008.8 m` and angular differences are expressed in radians.

The QGIS output geometry is returned in the original input CRS.

---

# 4. QGIS tool — quick start

After installation, the algorithm is found at:

```text
Processing Toolbox
└── Backbone Surface Rupture
    └── Build Backbone Surface Rupture
```

You can also type **Backbone Surface Rupture** in the Processing Toolbox search box.

The plugin welcome/About dialog contains an **Open Processing Toolbox** button.

## 4.1 First run

1. Load a rupture polyline layer in QGIS.
2. Open **Processing Toolbox → Backbone Surface Rupture → Build Backbone Surface Rupture**.
3. Select the input rupture layer.
4. Select the **Event ID field**.
5. Decide whether a classification filter is needed:
   - principal-fault-only layer → leave **Principal-fault classification field** blank;
   - mixed layer → select the classification field and enter the principal-fault value.
6. Enter one **Event ID**, or leave it blank to process all events in the layer.
7. Keep the recommended defaults for the first run.
8. Choose temporary or permanent output destinations.
9. Click **Run**.

The log reports the event, number of tested endpoint pairs, BSR length, global extent, trace fraction, number of gaps, and confidence class.

---

# 5. Outputs

## 5.1 Main QGIS outputs

The QGIS tool creates four outputs.

### A. Backbone Surface Rupture

A line layer representing the reconstructed BSR.

Typical attributes:

| Field | Plain-language meaning |
|---|---|
| `event_id` | Earthquake identifier. |
| `length_m` | Total BSR length, including artificial connections. |
| `network_cov` | Fraction of the total mapped principal-fault network length used by the selected BSR. |
| `trace_frac` | Fraction of the final BSR length that follows observed principal-fault traces. `0.85` means 85%. |
| `n_gaps` | Number of artificial connections. |
| `max_gap_m` | Length of the longest artificial connection. |
| `total_gap_` | Sum of all artificial-connection lengths. |
| `extent_cov` | Fraction of the global rupture extent represented by the BSR. A value near `1` means the path reaches both ends. |
| `confidence` | Qualitative reconstruction-confidence class (A–D). |
| `geom_flag` | Diagnostic geometry flag (for example, `OK` or `LOW_TRACE_COVERAGE`). |

### B. Artificial connections

A line layer containing only GAP portions.

Typical attributes:

| Field | Meaning |
|---|---|
| `event_id` | Earthquake identifier. |
| `gap_id` | Sequential identifier of the artificial connection. |
| `length_m` | Connection length in metres. |

### C. Ordered BSR vertices and x/L

A point layer ordered from the BSR start to end.

Typical attributes:

| Field | Meaning |
|---|---|
| `event_id` | Earthquake identifier. |
| `vertex_id` | Position of the vertex in the ordered sequence. |
| `distance_m` | Cumulative distance from the BSR start, including GAP segments. |
| `x_over_l` | Normalized position: `0` at the start and `1` at the end. |
| `trace_dist_m` | Cumulative distance travelled only on observed TRACE portions; it does not increase across a GAP. |
| `x_l_trace` | Trace-only normalized cumulative distance. |
| `edge_type` | `TRACE` for observed principal-fault geometry or `GAP` for an artificial connection. |

### D. BSR event summary

A non-spatial table containing the main reconstruction metrics for every processed event.

Artificial connections are densified in the QGIS vertex output as well as in the MATLAB reference output. Consequently, `distance_m` and `x_over_l` are sampled progressively along GAP segments rather than being reported only at their endpoints.

## 5.2 MATLAB event outputs

For each earthquake, the reference MATLAB implementation writes an event folder under the selected output directory, for example:

```text
BSR_RESULTS_FINAL/IdE_<EVENT_ID>/
```

Principal files include:

| File | Purpose |
|---|---|
| `BSR_IdE_<ID>.csv` | Ordered BSR vertices, cumulative distance, `x/L`, edge type, and source-feature information. |
| `BSR_edges_IdE_<ID>.csv` | Ordered graph edges selected for the BSR. |
| `BSR_IdE_<ID>.png` | Diagnostic figure. |
| `BSR_IdE_<ID>.fig` | Editable MATLAB figure. |
| `BSR_IdE_<ID>.mat` | Complete MATLAB variables for that event. |
| `BSR_selection_IdE_<ID>.txt` | Explanation of the hierarchical candidate selection. |

Development/diagnostic outputs may also include graph nodes, all graph edges, and all valid candidates.

A global summary is written as:

```text
BSR_summary.csv
```

## 5.3 MATLAB `BSR_IdE_<ID>.csv` fields

The simplified reference output contains:

| Field | Meaning |
|---|---|
| `IdE` | Earthquake identifier. |
| `VertexID` | Ordered vertex number. |
| `X_m`, `Y_m` | Local metric coordinates. |
| `Longitude`, `Latitude` | Geographic coordinates. |
| `AlongBSRDistance_m` | Cumulative BSR distance, including GAP segments. |
| `XoverL` | `AlongBSRDistance_m` divided by total BSR length. |
| `AlongTraceDistance_m` | Cumulative distance travelled only on observed TRACE portions. This value does not increase across a GAP. |
| `XoverL_TraceOnly` | Trace-only cumulative distance normalized by total observed distance used by the BSR. |
| `EdgeID` | Identifier of the graph edge represented by the row. |
| `EdgeType` | `TRACE` or `GAP`. |
| `OriginalPartIndex` | Index of the source observed polyline; empty/NaN for GAP rows. |
| `IdS` | Original SURE feature identifier; empty/NaN for GAP rows. |

Artificial connections are densified in the simplified MATLAB output, so `AlongBSRDistance_m` and `XoverL` are sampled along GAP segments rather than only at their endpoints.

---

# 6. Parameters explained "for dummies"

The recommended defaults are the values used during method development. **For a first run, do not change them.**

## 6.1 Parameters visible in the QGIS tool

### Maximum artificial connection distance (m) — default `15000`

**What it means:** the maximum straight-line distance that an artificial connection is allowed to span.

**Simple version:** two principal traces farther apart than 15 km cannot be connected directly.

**Increase it if:** a genuine rupture contains a larger unmapped interval and no complete BSR can be found.

**Decrease it if:** the algorithm connects obviously unrelated rupture clusters.

**Be careful:** increasing it makes the graph denser, slower, and more permissive.

---

### Off-trace travel-time multiplier — default `5`

**What it means:** travel across an artificial connection is five times more expensive than travel over the same distance along an observed TRACE.

```text
1 km on TRACE  = reference cost 1
1 km on GAP    = reference cost 5
```

**Increase it if:** you want an even stronger preference for observed principal-fault geometry.

**Decrease it if:** the mapped principal traces are very incomplete and continuity across gaps needs more freedom.

**Important:** this is an algorithmic weight, not a physical seismic rupture velocity.

---

### Maximum allowed turn (degrees) — default `140`

**What it means:** the largest directional change accepted between consecutive path elements.

**Simple version:** it prevents strong reversals and U-turn-like solutions.

**Increase it if:** a strongly curved but geologically reasonable rupture is rejected.

**Decrease it if:** the reconstructed path makes implausibly sharp turns.

---

### Tip-to-feature junction tolerance (m) — default `300`

**What it means:** when a principal-trace tip lies within 300 m of the interior of another principal trace, the graph may create an internal topological junction.

**Why it exists:** the path can enter the receiving feature near the geometrically meaningful point instead of travelling to one of its original ends and reversing direction.

**Important:** the original mapped polylines are not moved or redrawn.

---

### Maximum connections per significant exit — default `50`

**What it means:** from each significant exit location, the graph keeps at most the 50 nearest admissible principal-trace tips.

**Why it exists:** without this limit, dense rupture networks can generate thousands of unnecessary possible connections.

**Increase it if:** a necessary connection appears to be missing in a very dense network.

**Decrease it if:** processing is too slow and the geometry is simple.

---

### Minimum global extent for completion — default `0.98`

**What it means:** a candidate is considered complete when its endpoints span at least 98% of the global principal-fault extent.

**Simple version:** the BSR should practically reach both ends of the rupture.

**Recommendation:** normally leave this unchanged.

---

### Adaptive endpoint levels per side — default `10,20,30,50`

**What it means:** the algorithm first considers 10 external tips on each side. If no complete path is found, it expands the search to 20, then 30, then 50 tips per side.

**Why it exists:** testing every possible terminal pair can be unnecessarily slow.

**Increase the final level if:** the input contains many short fragments and the true outer endpoints may not be represented among the first 50 tips.

---

## 6.2 Advanced geometry parameters in the MATLAB reference implementation

### `minimumSplitDistance_m = 25`

An internal projected junction must lie at least 25 m from the receiving feature ends. This prevents nearly redundant graph nodes from being created next to an existing tip.

### `splitMergeTolerance_m = 5`

Split positions less than 5 m apart are treated as the same graph location.

### `nodeSnapTolerance_m = 1`

Graph nodes are clustered only when they are effectively coincident within 1 m. This is **not** a general endpoint-snapping operation and does not move the observed trace geometry.

### `duplicateVertexTolerance_m = 0.05`

Consecutive vertices closer than 5 cm are considered duplicate numerical vertices and removed.

### `minimumCleanFeatureLength_m = 1`

A cleaned feature shorter than 1 m is treated as geometrically degenerate.

## 6.3 Douglas-Peucker significant vertices

Douglas-Peucker simplification is used **only to identify significant graph access/exit nodes**. It does not replace the observed polyline geometry.

For a feature of length `L`, the reference tolerance is:

```text
max(50 m, min(500 m, 0.002 * L))
```

Reference parameters:

- `minimumSimplificationTolerance_m = 50`
- `maximumSimplificationTolerance_m = 500`
- `simplificationLengthFraction = 0.002`

**Smaller tolerance:** more graph nodes and more possible exits, but slower calculations.

**Larger tolerance:** fewer nodes and faster calculations, but potentially fewer useful exit locations.

## 6.4 Final continuity check

### `assemblyContinuityTolerance_m = 1`

After the preferred path is selected, the output geometry is checked for continuity. If consecutive selected geometries remain more than 1 m apart, the finite interval is exported explicitly as a GAP.

This normalization does not rerank the candidates.

---

# 7. How candidate selection works

The BSR does not combine all criteria into a single weighted score. Candidate paths are ranked **hierarchically (lexicographically)**.

In simple terms, the algorithm asks:

1. Does the path represent the full rupture extent?
2. Among equally complete paths, which uses the most observed principal-fault geometry?
3. Among those, which requires fewer artificial connections?
4. Then: which has the shorter maximum gap?
5. Then: which has the shorter total gap length?
6. Then: which has more consistent directional progression?
7. Then: which has less backtracking?
8. Finally: which has the lower travel time?

A lower-priority criterion cannot compensate for a clearly worse higher-priority criterion.

---

# 8. Reconstruction confidence

The confidence class summarizes **observational support for the selected BSR**. It is assigned after the preferred BSR has been selected and therefore does not control path optimization.

It is not a probability and is not a statement that the reconstruction is "right" or "wrong".

General interpretation:

| Class | Interpretation |
|---|---|
| **A** | Reconstruction strongly supported by observed principal-fault geometry, with absent or very limited inferred continuity. |
| **B** | Good observational support with moderate inferred continuity. |
| **C** | Substantial inferred continuity is required because observations are fragmented or contain longer gaps. |
| **D** | The available principal-fault geometry provides comparatively weak or competing support for a unique backbone. |

Always interpret confidence together with the event map, trace coverage, gap statistics, and global extent.

---

# 9. How to interpret the result

A successful BSR should normally:

- cover almost the complete global rupture extent;
- start and end on observed principal-fault traces;
- preserve mapped geometry along TRACE portions;
- avoid strong backtracking;
- use artificial connections only where continuity is required.

Useful warning signs:

| Observation | Possible meaning |
|---|---|
| `extent_cov < 0.98` | The path did not reach the intended global rupture extent. |
| Very low `trace_frac` | Much of the final BSR consists of inferred connections rather than observed principal-fault geometry. |
| Low `network_cov` | The selected BSR uses only a limited fraction of the complete mapped principal-fault network; this may occur in branching or parallel networks. |
| Very high `n_gaps` | The mapped principal rupture is composed of many short fragments. |
| Very long `max_gap_m` | A large unmapped interval is required to preserve continuity. Inspect it visually. |
| Nonzero backtracking | The selected route locally opposes the dominant start-to-end progression. |
| Low confidence | The selected geometry is less strongly supported by the available principal-fault observations; it is not automatically incorrect. |

Always inspect the map together with the numerical metrics.

---


# 10. Citation and license

The software is distributed under the **MIT License**.

The MIT License permits free use, modification, redistribution, and commercial use provided that the copyright and license notice are retained.

Scientific users are requested to cite the publication describing the method. Provisional citation:

```text
Visini (2026), Backbone Surface Rupture: a vector framework for
reconstructing the principal surface rupture from mapped principal-fault traces.
Manuscript under review.
```

The final journal citation and DOI will replace this provisional reference after publication.


---

# 11. AI-assisted development disclosure

Generative AI was used to assist with code implementation,  debugging, documentation, and preparation of the MATLAB and Python/QGIS codes. Methodology, algorithm design, parameter selection, validation, and interpretation are under the responsibility of the author.
