# cellucid 0.9.8 <!-- CELLUCID_VERSION -->

Version 0.9.8 is the CRAN submission release.

It is a patch release on 0.9.7 and stays in step with the Python package of
the same version. Stores written by `cellucid_prepare()` and sessions saved
(format 4) with 0.9.7 are read as before. The viewer protocol is still schema
5: use this release with the web app of this release.

## The R package

- A cells x genes `DelayedArray` (`gene_expression =`, or an assay
  transposed with `t()`) is read under its transposition, each block
  transposed in memory, never through DelayedArray's transposed view: on
  R 4.1 and 4.2 (Bioconductor 3.14 to 3.16) with the matrixStats 1.x CRAN
  serves them, a sparse block read through that view stopped
  `cellucid_prepare()` with "Argument 'useNames' must be either TRUE or
  FALSE".
- The viewer token stays out of what is logged or pasted. A printed viewer
  showed its address with the token (`address:`), and a wait that timed out
  said `No viewer is connected: open <the address with its token>`; a printed
  viewer lands in knitted reports and logs, an error in bug reports. Both now
  say the address without the token and name `viewer$viewer_url`, the
  address that opens it (`.public_viewer_url()`, which
  `cellucid_debug_connection()` uses too). The address with it is said only
  where you are told to open it: `cellucid_serve()`'s message,
  `cellucid_view(pane = "none")` (which printed nothing, though its help
  said it prints the address) or a viewer opened outside an interactive
  session, and `cellucid_mcp(open = )` on stderr.
- A Zarr store whose declared id or name cannot be read (a damaged chunk)
  is refused with that sentence ("… Write the store again, or give
  dataset_id = yourself."), unless the caller names that field:
  `cellucid_serve()` and `cellucid_view()` no longer need from the store a
  field given as `dataset_id =` or `dataset_name =`, as cellucid-python.

- A viewer names what opened it to the web app (`_cellucid/info`'s
  `runtime`, and the agents folder's record): `cellucid_view()` and
  `cellucid_serve()` say `r`, `cellucid_mcp(open = )` says `r-mcp` (it said
  `r`), so the viewer's sentences name "the R session" or "the agent app's
  cellucid_mcp() server" ("… which opened this viewer, stopped (its server
  no longer answers). … Restart the cellucid_mcp() server in the agent
  app."), not "the notebook or R session". The shared viewer protocol file
  lists every runtime with its label and language.
- `cellucid_field_summary()` takes `within`: summarize over the cells of some
  categories of another field, or in a range of a continuous field or gene,
  without filtering the viewer.
- `cellucid_history_summary()` returns the notes on the tracker, each branch,
  group and turn in full (each with its id, text, author and time), not their
  counts; a turn's former `summary` is the first of its `notes`. A received
  session's branch notes were otherwise unreadable to an agent or a script.
- `cellucid_field_restore()` and `cellucid_field_purge()` take `which`: when
  two deleted fields have one name, which of them, in the order the web app's
  **Deleted fields** lists them (named by the name alone, the request was
  refused with no way to name either).
- The identity a Zarr store declares (`uns/cellucid` id and name) is read
  checking every count and length its vlen-utf8 chunk states against the
  bytes it holds: four bytes stating 2^31 - 1 strings allocated a
  character vector that long before anything was read. Such a chunk is
  refused as a damaged store.
- A viewer stopped with `cellucid_stop()` says **Viewer stopped** on its
  frozen last view (the web app of this release), as one whose R session went
  away does; it froze without a word.
- Offline, the cached web app is served as before, now saying so also when
  `quiet`, and the error when nothing can be served names what is missing
  ("no copy of it is cached in …" or "the copy cached in … does not check
  out (…)") and how to get the app, in the words cellucid-python now uses too.
- A new build of the web app is staged from the files the cache already holds
  unchanged, and the others are downloaded six at once (one at a time, the
  1,209 files of a build took minutes); the build the cache holds is used as
  before, as cellucid-python now does too.
- `cellucid_prepare()` refuses a `bit64::integer64` column saying to write
  identifiers `as.character()` (the viewer reads integers past 2^53 as
  text) or counts within 2^53 `as.numeric()`, and refuses an element name
  longer than 255 bytes in UTF-8, the longest file name every platform takes.
- `cellucid_apply_session()` takes a logical column's categories as `false`
  and `true`, as the viewer shows them: a label renamed in the viewer is
  written (the column becomes text) instead of "holds numbers, not
  categories".
- The embeddings a store opens with treat every width alike, as the web app
  and cellucid-python now do (shared conformance fixture, version 9): the
  default rule names the UMAP's width and every other width takes its
  best-ranked array, so `X_umap` beside `X_umap_3d` names both 2D and 3D
  (`.resolve_embedding_slots`, checked against every `embeddingSlots` case).
- The default neighbor graph rule takes a MuData modality's own graph
  (`'<modality>:connectivities'` of the modality a dataset opens with) after
  the dataset's `connectivities`, as the web app and cellucid-python now do
  (shared data-mapping conformance fixture, version 8): a MuData whose graph
  is in a modality, as scanpy and cellucid-python's `prepare()` write it,
  opens with that graph in the viewer.
- `cellucid_mapping_set()` (the viewer's `mapping.set`) moves the genes groups
  it does not name with a named expression matrix or modality, as the
  viewer's Data mapping panel does: a gene-name or highly-variable column the
  new matrix's var axis has is kept, else that axis's default. A session
  saved on `X` with `highly_variable` restores onto a file that opens on
  `raw.X`.
- `cellucid_freeze()` is removed. `cellucid_stop()` tells an open page to
  freeze before it stops serving it; called on its own, `cellucid_freeze()`
  froze a page R still served, and cellucid-python has no such call.
- `cellucid_prepare(total_counts_over =)` (new) says what the column
  `total_counts` names covers: `"all genes"` (the default) or
  `"source panel"`, every gene of the source's gene panel -- a targeted or
  integration panel, more genes than are stored but not the whole
  transcriptome -- recorded as the `over` of `uns['cellucid']['totalCounts']`
  (the specification's third value) and kept when an object carrying the
  record is prepared again. With several modalities, one value or
  `c(RNA = "source panel")`. cellucid-python's `prepare(total_counts_over=)`
  is the same.
- `cellucid_prepare()` writes a column's colors (`obs_colors =`, or those
  an object holds) only when each category's is hex (`#rrggbb`, `#rgb`,
  `#rrggbbaa`), the forms the viewer reads: a palette holding a named color
  (`"red"`) was written and the viewer opened the column in its own colors.
  Given as `obs_colors =`, it is refused; held by the object, left out with a
  message. cellucid-python's `prepare()` checks the same.
- `cellucid_view()` and `cellucid_serve()` refuse a file Cellucid does not
  open in the web app's words, as cellucid-python does: "<path> is not a file
  Cellucid opens (an .h5ad file, …, an .rds file, … or a .zarr.zip archive)."
  It listed one hand-picked ending per format (`.rds.gz`, `.rdata.gz`, and the
  refused `.h5seurat`).
- `cellucid_view()` and `cellucid_serve()` say why a viewer could not bind
  its address in cellucid-python's sentence, with R's spelling: a `host`
  this machine cannot serve from ("Host 'x' is not an address this machine
  can serve from. Use host = \"127.0.0.1\" ...") and a chosen `port` already
  in use ("Port 9000 is already in use on 127.0.0.1, ..."). A bad host was
  said as a port taken ("Could not serve on not a host:8765 (Failed to
  create server). Choose another port, ...").
- `cellucid_prepare()` and `cellucid_view()` write a dataset's name,
  description, citation text and source name with white space at either end
  removed, as the web app now reads them (it refused a whole file whose
  description ended in a space); identifiers keep the strict rule. The same
  in cellucid-python.
- `cellucid_prepare(data = seurat)` names every layer of an assay it does not
  write: `counts` beside the `data` it writes, and `scale.data`, in a v3
  assay as in a v5 one (`assays$RNA layers not read: c("counts",
  "scale.data")`), as cellucid-python names the layers and `raw` of an
  AnnData. A v3 assay's were left out unsaid, and a v5 assay's `counts`.
- US spelling throughout the package's messages and help ("neighbor graph",
  "normalizes", "catalog", "labeled"), as in the web app and the Python
  package. Names are US too: `cellucid_expression_transform_set(normalize =)`,
  `cellucid_vector_field_summary(neighbors =)` and its result's `neighbors`,
  `center` in `cellucid_camera_framing()`'s groups, and the expression values
  `list(normalize, log1p)` a session and `cellucid_apply_session()`'s record
  hold (a step recorded with `normalise` replays as `normalize`).
- A change the viewer applied without a step of its history now says why.
  The viewer's `action_result` carries `unrecorded` (`"unchanged"`: it was so
  already; `"retracted"`: it put back the step just before it, which was
  taken out; `"off"`: the History settings leave this kind of change out;
  `"paused"`; `"elsewhere"`: another tracker was shown meanwhile;
  `"no-history"`: the dataset's history has not started), kept as
  `request$unrecorded` and printed with the request; `cellucid_action_result()`
  returns the outcome and warns it (class `cellucid_action_not_recorded`), as
  in "cellucid_color_by(): it changed nothing (it was so already), so no step
  was recorded." A read records no step by nature and warns nothing. Before,
  `request$step` was `NULL` with no reason, whether the change happened or
  not (cellucid-python's `ActionRequest.unrecorded`; the shared viewer
  protocol fixture says so).
- `cellucid_dataset_summary()` lists every embedding the viewer's data
  mapping offers, on or off, with its obsm `key`, `columns`, `width`,
  `enabled`, `usable` and `reason`: the names `cellucid_view_embedding_set()`
  and `cellucid_mapping_set()` take (it listed the embeddings shown by name
  and dimension only).
- `cellucid_apply_session()` writes every categorical user field with the
  colors it is shown in, where the object keeps what an AnnData keeps in
  `uns`, under Scanpy's name `<key>_colors` (Seurat `misc`, `metadata(x)` --
  which zellkonverter carries to `uns` -- or a data frame's attribute; one
  `#rrggbb` a level, every level): in the view in focus at the applied step,
  the color the person or an agent gave a level (the session's `colors`
  value), else the color it opens with -- its `palette` (read and checked: a
  field made from highlight pages opens in its pages' colors), else
  Cellucid's own category colors by position (the web app's
  `getCategoryColor`, ported in `R/category-colors.R` and pinned by the
  shared fixture `category_default_colors.json` in the three repositories).
  The plan's colors row stays `not-applied` -- every other color is how the
  data was shown -- and says so when a field's colors went with its column.
  A column's colors follow it in every write of an apply and a pull: a
  column without categories (a highlight group, a continuous field) has none
  (stale colors under its name are removed); a renamed column carries its
  colors; merged levels keep their first level's color (colors that did not
  cover the levels are removed); a text column relabeled loses them; a
  column a viewer pull removes takes them along. `cellucid_pull()` of a
  session file writes a field's colors the same way, under the column name
  it writes. Before, no colors were written but a highlight-page field's
  palette (a grouped or merged field's chosen colors never reached the
  object), and a rename or a merge left colors that described another
  column. cellucid-python and the web app's **Save AnnData with this state…**
  write the same.
- A column of the object (`seurat_clusters`, `cell_type`; a factor) recolored
  in the view in focus is written with the colors it is shown in too, in the
  same place: one `column-colors` row of the apply plan per column, the levels
  recolored, the rest in the colors the column opens with (its own
  `<column>_colors` as the viewer draws them, else Cellucid's), under the
  column's name once the apply's renames are made; `unchanged` when the
  object holds them already; `not-applied`, saying why, for a column that is
  not a factor or one with none of the levels recolored. Before, the colors a
  person gave a dataset column never reached the object. And
  `cellucid_pull()` from a live viewer writes a categorical field's colors:
  each categorical result of the viewer protocol carries the colors the view
  in focus shows it in (`colors`, one `#rrggbb` a category), checked when it
  arrives. Before, a viewer pull wrote no colors. cellucid-python and the web
  app's copy do the same.
- A viewer given another `dataset_id =` than the one its object declares
  (`cellucid_view(object, dataset_id = )`, or a Zarr store's own id) takes
  its own session back: `cellucid_apply_session(object, request)` refused
  it, since the session is saved under the id served and the object declares
  its own. A viewer's request now names the id it served the object under
  for an object declaring the id the viewer's data declared (the viewer
  keeps it, `declared_id`); `cellucid_apply_session()` and `cellucid_pull()`
  of a session file take `dataset_id =` among the arguments
  `cellucid_view()` was given, the id the object was shown under, and the
  refusal says to name it. A session of another dataset is refused as
  before, and the cell count, data mapping and cell order are checked as
  always.
- The object's cellucid record is read by its exact name: `misc$cellucid` and
  `metadata(x)$cellucid` were read with `$`, which matches a partial name, so
  an object holding a `cellucid_highlight__..._colors` entry (as Scanpy writes
  after plotting a highlight column, carried over by zellkonverter) and no
  record was read as if that were its record: every apply and pull of it was
  refused ("The object's cellucid record must be a list"), and
  `cellucid_prepare()` refused it too ("misc$cellucid must be a list ...").
- `cellucid_category_delete()` moves the cells into `"Unassigned"` by default
  (the one spelling of cells without a label, as Create categorical's; it was
  `"unassigned"`).
- `cellucid_pull()` of a field's obs-column CSV keeps the field's categories
  and their order: the web app writes a third row, `categories`, holding each
  categorical column's categories in the field's order as a JSON array of
  strings, and the factor's levels are that list. They came back in the order
  the rows first showed them, so Pancreas `clusters` lost its trajectory
  order and a category no cell held was lost. A label the row does not list,
  a categorical column without a list and a CSV without the row (the two-row
  layout of earlier exports) are refused; export the field again.
- `cellucid_view_rename()` (new, generated from the viewer's `view.rename`)
  names a view or, with `NULL`, takes its name away; it shows beside the
  view's number in the viewer and titles its panel in an all-views figure.
  `cellucid_view_keep(name = ...)` keeps and names a view in one step;
  `cellucid_views_state()` and `cellucid_view_picture()` give each view's
  name (`label`).
- `cellucid_vector_field_summary()` (new, generated from the viewer's
  `vector-field.summary`): where the cells of each category are heading
  along a vector field (RNA velocity) in the embedding a view shows: their
  speed against the average cell, how aligned their arrows are, the mean
  direction, and the categories their arrows lead into.
  `cellucid_view_picture()` takes `source = "screen"` (the view as drawn,
  the vector field's flow included) or `source = "export"` (the whole figure
  as the viewer's Export would write it now, with its Figure export
  settings, nothing written) and `cellucid_annotation_suggestions()`
  takes `pull` (it pulls the repository first when the browser never did).
- `cellucid_annotation_publish(viewer, votes = list(...), suggestions = list(...))`
  (new, generated from the viewer's `annotation.publish`) votes and suggests
  community annotation labels together, published as one commit of the
  person's file and asked of them once (cellucid-python's
  `annotation_publish()`); `cellucid_annotation_vote()` and
  `cellucid_annotation_suggest()` are removed with the viewer's actions.
- `cellucid_history_saved_open()` (new, generated from the viewer's
  `history.saved.open`, an agent's opening of a saved state to continue it):
  the viewer refuses it from R, as `cellucid_dataset_open()`, since a
  notebook's viewer shows the dataset its notebook serves.
- `cellucid_layout_panel_dock_all()` (new, generated from the viewer's
  `layout.panel.dock-all`, the sidebar's new **Dock all** button): docks
  every floating sidebar section back into its place, as one step.
- `cellucid_prepare()` records each cell's total count (its library size)
  when `X` holds counts: `obs['total_counts']`, the sum of each row of `X` as
  written, with `uns['cellucid']['totalCounts']` (`over = "stored genes"`),
  so the web app normalizes a gene from one column instead of reading every
  gene once. `total_counts` names the object's own column computed over every
  gene (`"nCount_RNA"`; several modalities: `c(RNA = "nCount_RNA")`),
  recorded `over = "all genes"` and refused when it is below a cell's stored
  genes' sum or `X` holds no counts; each counts modality of a multimodal
  object gets its own `mod/<name>/obs/total_counts`. The same as
  cellucid-python's `prepare(total_counts=)`, byte for byte.
- `cellucid_prepare()`, `cellucid_view()` and `cellucid_serve()` of an object
  keep Scanpy's colors of the columns they write: each categorical `obs`
  column's colors the object holds (`misc[["<column>_colors"]]` of a Seurat
  object, `metadata(x)[["<column>_colors"]]`, where zellkonverter keeps
  `uns`, and where `cellucid_apply_session()` writes them) are written as
  `uns['<column>_colors']`, so the viewer opens the column in them instead of
  Cellucid's colors. `obs_colors` (new) gives or replaces them
  (`list(cell_type = c("#1f77b4", ...))`; `list()` writes none). Colors that
  cannot color the written column (not text, fewer than its categories, a
  column not written as a categorical) are left out with a message, or
  refused when given. Before, every column's colors were dropped. The same
  rule as cellucid-python's `prepare(obs_colors=)`; the Cellucid-ready
  AnnData specification states it.
- `cellucid_apply_session()` records the expression values the genes were
  shown in at the step applied: each record of `cellucid$history_applied`
  (`cellucid_history_applied()`) has `expression`, `list(normalize, log1p)`
  as the step's **Expression values** held it (`NULL` for a session without
  them), so a notebook knows what the person looked at. Recorded, never
  applied: the object's expression values stay as stored. The same field as
  cellucid-python's `apply_cellucid_session_to_anndata()` and the web app's
  **Save AnnData with this state…** (history model, shared golden
  provenance).
- One spelling, as in the web app: the help and messages say color, colored,
  coloring (they said "colour" in places), and so do the generated functions'
  titles and the step names a saved history holds.
- US spelling everywhere else too: DESCRIPTION, the help, the vignettes and
  the messages said `neighbour graph`, `catalogue`, `labelled` and `centre`
  in places; now "neighbor", "catalog", "labeled", "center".
- A refusal the viewer answers when the person stops a dataset's or a saved
  state's loading has the code `canceled` (it was `cancelled`), as the web app
  spells it; `history.branch.cancel` still answers `cancelled`.
- The methods paragraph writes a share of the cells as the web app's
  `formatShare` does: whole from 10% (never 100% short of every cell), two
  significant figures below ("5.5%", where it said "6%"; "99%" for 19,999 of
  20,000 cells, where it said "100%").
- The action catalog shipped with the package (`inst/action_catalog.json`)
  is the web app's of this release: an action may declare what it captures
  (`captures`), so a camera path keyframe is never marked stale when the
  camera moves on.
- `cellucid_figure_export_get()` reads the viewer's Figure export settings,
  by the names `cellucid_figure_export_set()` takes.
- `cellucid_figure_export_set(frame = "cells")` frames every cell the figure
  draws, as the viewer's **Fit the cells** does, so a plot takes the cells'
  shape instead of a narrow pane's; `svg_strategy` is never `NULL` (an SVG is
  full vector until another strategy is chosen), as the viewer's
  `figure-export.set` version 2.
- `cellucid_figure_export_set()` sizes the figure as it is printed, as the
  viewer's `figure-export.set` version 3: `width` is the figure's width in
  millimetres (or `size_preset = "single-column"`, `"one-and-half-column"`,
  `"double-column"`), which the figure fills at its own height, and
  `font_sizes` are points; `height` is gone. "180 mm wide" is `width = 180`.
- `cellucid_highlight_category()` takes `categories`, one or more labels of
  one field (was `category`, one): each becomes a group of the page in one
  step of the viewer's history, one undo for all of them, as the viewer's
  action `highlight.category` version 2 does.
- `cellucid_highlight_cells()` without `label` makes a group named
  "Highlight" (the viewer's list shows its cell count), not "Highlight (N
  cells)"; its help page and the `cellucid_mcp()` tool say so.
- `cellucid_analysis_export()` writes an analysis's result table as the
  viewer's expanded view's CSV button does (marker genes, with `all = TRUE`,
  every gene of every group); `cellucid_category_merge()` takes `label`, the
  merged category's name (giving the `into` category's own name keeps it).
- `cellucid_camera_framing()` says where cells sit in a view and gives
  cameras framing the whole view, each category of a field and the
  categories named together, from an azimuth and elevation;
  `cellucid_camera_path_play()`, `cellucid_camera_path_stop()` and
  `cellucid_camera_path_get()` play, stop and read the viewer's camera path.
- `cellucid_field_group()` makes a new categorical field of groups of a
  field's categories in one step (`groups`: a list of `list(label =,
  categories =)`; `rest`, `name`); `cellucid_dataset_open()` takes `count` or
  `fraction` (and `method`, `seed`, `field`, `allocation`, `embedding`) to
  open a subsample at once, and `cellucid_cells_open()` takes `fraction`
  beside `count`.
- `cellucid_camera_path_keyframe_add()` takes `hold` and
  `cellucid_camera_path_hold_set()` sets a keyframe's hold (or every
  keyframe's): the seconds the camera rests there before moving on.
- `cellucid_camera_path_stop()` stops a path `cellucid_camera_path_play()` is
  still playing, at once (it waited for the path to end and was then
  refused), and returns `stopped`; `FALSE` when nothing played. It and
  `cellucid_analysis_cancel()` are documented as stopping work under way
  (the catalog's new effect `control`).
- A camera path that `cellucid_camera_path_play()` (or an agent) is playing
  when the person stops the turn, or the call is canceled, stays where the
  film stopped, as `cellucid_camera_path_stop()` leaves it, and the play's
  step holds that camera; any other change stopped is still rolled back.
- `cellucid_camera_path_playback_set()` takes `captions`: a playing camera
  path shows each named keyframe's label as a caption over the plot (on by
  default); `cellucid_camera_path_get()` reads it.
- A looping camera path that `cellucid_camera_path_play()` started belongs
  to the turn that started it: the person's Stop of that turn stops it where
  the camera is. A stopped path's playback bar shows where the film stopped.
- A parameter the viewer refuses for its form says the form ("must be a
  name with no space at its start or end", "must be a view name such as
  \"View 2\"", "must be like \"s3\""), in the error a call raises, instead of
  "does not have the expected form".
- `cellucid_view(path, ...)` and `cellucid_serve(path, ...)` take
  `knn_key` (the `obsp` graph drawn) and `expression_source` (`"X"`,
  `"raw"` or `"layer:<name>"`, the matrix genes are read from), the launch
  defaults cellucid-python's `serve_anndata()` states for a file; R refused
  both.
- `cellucid_wait_event(viewer, event, timeout, predicate)` returns the next
  event of a kind, as Python's `viewer.wait_for_event()`.
- `cellucid_on()` takes `"message"`, every event before its own handlers with
  its kind as `event` (Python's `on_message`), and refuses an event the viewer
  never sends with the list of those it does, the names Python's hooks take
  (pinned to the shared protocol fixture); it partially matched (`"sel"`
  was `"selection"`) and refused with R's generic words.
- `cellucid_debug_connection()` reports `viewer_url` without the viewer
  token (the report is what you paste when asking for help; the page's own
  debug snapshot left it out already); `viewer$viewer_url` is the address to
  open. Python's `debug_connection()` does the same.
- `cellucid_mcp()` advertises `tools.listChanged: false` and never sends
  `notifications/tools/list_changed`: the agent's tool list is the same
  whatever the person's settings (the viewer tells the agent those after its
  tool results), so changing a setting or reloading the page no longer
  makes the agent app rewrite its cached prompt and drop the agent's earlier
  reasoning (a permission change made the next turn cost about 4× its cached
  price on Opus 5.5). The viewer channel's `agent_tools_changed` event is
  gone, as in cellucid-python.
- The agents folder's owner is checked by user id, as cellucid-python checks
  it: a user with no name (a container's or a cluster job's uid without a
  passwd entry) had its own folder refused as "another user owns it", so no
  viewer was published and `cellucid_mcp()` found none.
- A session file holding a value nested more deeply than R's stacks follow is
  refused with "This session file holds a value nested more deeply than this
  package reads …" (R 4.4 and later; cellucid-python says the same), not R's
  "evaluation nested too deeply" error.
- The tests compare each file this package shares with the web app (the
  action catalog, the MCP tools, the history model, the input formats, the
  viewer protocol, the conformance fixtures, the category
  colors, the golden history, the round-trip fixtures) with the web
  repository's file, and each file it shares with cellucid-python alone (the
  notebook widget's front end, `local_agents.json`, the transition-drift
  fixture, the MuData cross-writer harness) with cellucid-python's, byte for
  byte, through one helper (`expect_shared_file()`). Where a sibling is not
  checked out (an installed package, a CRAN check) the comparison is skipped
  and says why; CI fetches both siblings at `main` and sets
  `CELLUCID_REQUIRE_SIBLINGS=1`, which makes a missing one a failure. A
  difference names the one command that makes every copy agree
  (`node scripts/sync-shared-files.mjs` in the web repository).

## In the web app

- `cellucid_field_create_from_pages()` takes `rest`, a categorical column
  whose category the cells on no page keep (splitting one label is one
  call), and `cellucid_highlight_page_delete()` takes `pages`, several pages
  deleted in one step. A session file given to `cellucid_serve()` is refused
  with a sentence naming the web app's way in: a dropped session file now
  loads there, its dataset opened with it.

The viewer of this release opens a dataset with no starting view (an object
served from R among them) colored by its best categorical field: a cell
type, another annotation, a lineage, clusters, the tissue, a condition or the
cell-cycle phase, never a sample, batch, id or quality flag; the analyses'
default fields and the agent chat's starters follow the same order. The
Filters list names the smaller side of a legend, a highlight group is named
by what it holds ("NK + ILCP") with its count beside it, a category name cut
short in the legend shows whole on a click, a background you never picked
follows the theme, a drag brought back to where it began leaves no history
step, and Load... names an autosave by its latest step. A dataset opened
whole is shared whole, and a sample opened in volumetric smoke keeps its
starting view. Analyses name
every gene and group of a dot plot, rank tied volcano-plot genes by their
test statistic and list their display options in the expanded view; PNG
figures record their DPI and a Hybrid SVG draws its points at 300 dpi. A
local server (`cellucid_serve()`, `cellucid_view()`) is never paced as a host
that blocks networks. A gene colors a kept view (it did nothing), and a
cell near the camera is drawn up to 160 pixels across (was 128), whole and
opaque however near the camera comes. See the
documentation and the cellucid-python changelog
for the details.

# cellucid 0.9.7

Version 0.9.7 is the CRAN submission release.

It is a patch release on 0.9.6 and stays in step with the Python
package of the same version. Stores written by `cellucid_prepare()` in 0.9.6
are read as before. Session files are format 4: a session saved with 0.9.6 or
earlier (format 3 or before) is refused, saying which format it is. The viewer
protocol is schema 5: use this release with the web app of this release.

## The R package

- **Every action of the viewer, from R.** `cellucid_act(viewer, "<type>",
  ...)` sends any action of the catalogue (`inst/action_catalog.json`, byte
  for byte the web app's, pinned in the three repositories), and each of its
  177 actions has its own generated function with a help page
  (`cellucid_color_by()`, `cellucid_highlight_cells()`,
  `cellucid_category_rename()`, `cellucid_history_get()`,
  `cellucid_view_picture()`, `cellucid_dataset_open()`, ...), among them
  `cellucid_history_anndata_save(viewer, file)`, which has the viewer write
  a new Cellucid-ready AnnData `.zarr.zip` of its dataset with the state
  shown applied (downloaded by the browser), and
  `cellucid_history_saved_delete(viewer, saved)`, which deletes a saved
  state of the viewer's browser by its id. A change is one step
  of the viewer's history, labelled as R's; `cellucid_history_undo()`,
  `cellucid_history_redo()` and `cellucid_history_jump()` move the history
  and record no step. The functions of acts outside the viewer's state and
  history record no step either, and Undo does not take them back (their
  help pages say so): writing a file the browser downloads
  (`cellucid_history_share()`, `cellucid_history_report()`,
  `cellucid_history_anndata_save()`, `cellucid_figure_export()`,
  `cellucid_cells_export()`, `cellucid_field_export()`), publishing
  (`cellucid_annotation_vote()`, `cellucid_annotation_suggest()`,
  `cellucid_dataset_share_link()`; a vote or suggestion publishes only itself
  and is asked of the person first), and what
  the browser keeps (`cellucid_history_save()`,
  `cellucid_history_saved_delete()`, `cellucid_my_datasets_forget()`,
  `cellucid_stored_remove()`). `cellucid_history_tour_start()` and
  `cellucid_history_tour_stop()` play the viewer's guided tour;
  `cellucid_stored_list()` and `cellucid_annotation_suggestions()` read;
  `cellucid_datasets_list()` gives each dataset's size. A request returns at once; `cellucid_action_result()`
  collects the outcome or stops with the viewer's reason. Cells are given as
  R gives them (1-based positions, `TRUE`/`FALSE`, cell names). The
  functions, their help pages and the NAMESPACE exports are written by
  `Rscript tools/generate-actions.R`. A notebook's viewer refuses
  `dataset.open`.
- **Removed, no aliases:** `cellucid_highlight()` (use
  `cellucid_highlight_cells()`), `cellucid_clear_highlights()`
  (`cellucid_highlight_clear_all()`), `cellucid_color_by(viewer, "name")`
  (`cellucid_color_by(viewer, field = list(source = "obs", key = "name"))`),
  `cellucid_set_visibility()` (`cellucid_cells_visibility_set()`),
  `cellucid_reset_camera()` (`cellucid_camera_reset(viewer, view = "live")`),
  `cellucid_set_data_mapping()` (`cellucid_mapping_set()`) and
  `cellucid_push_gene_set()` (`cellucid_gene_set_save()`). `cellucid_push()`
  sends `field.push` and returns its request.
- **Reading web sessions.** `cellucid_read_session()` reads every format-4
  file the web app writes (the history, the current state only, or the steps
  only), checked whole against the history model
  (`inst/history_model.schema.json`, pinned in the three repositories), in
  time near linear in the steps. A damaged file says what is wrong, and a
  file from a newer Cellucid asks you to update the package. The history as
  data frames and records: `cellucid_trackers()`, `cellucid_tracker()`,
  `cellucid_branches()`, `cellucid_branch()`, `cellucid_groups()`,
  `cellucid_group()`, `cellucid_steps()` (with the columns `kind`, `v` and
  `group`), `cellucid_step()`, `cellucid_step_state()`,
  `cellucid_step_diff()`, `cellucid_stale()`, `cellucid_notes()`,
  `cellucid_turns()`, `cellucid_journal()`, `cellucid_conversations()`, and
  `cellucid_methods()` and `cellucid_pseudocode()` in the web app's wording.
  `session$chrome` holds the selection tool, read only. A squashed step's
  record says why it was made (`squashed`), and the trackers table whether
  each records (`recording`, `paused_at`). Removed, no aliases:
  `cellucid_session_trackers()`, `cellucid_session_steps()`,
  `cellucid_session_state()`.
- **Applying a history, a plan first.** `cellucid_apply_session(data,
  session, tracker = NULL, branch = NULL, step = NULL, ..., dry_run = TRUE,
  overwrite = FALSE, keep_session = FALSE)` returns a `cellucid_apply_plan`
  (one row per change, with its status and why), following the same rules as
  cellucid-python row for row. `cellucid_apply_plan(plan)` (or
  `dry_run = FALSE`) writes it onto a Seurat, SingleCellExperiment or
  SummarizedExperiment object and records the apply in
  `cellucid$history_applied` (`cellucid_history_applied()`). Applying again
  changes nothing, and a column you made is never replaced. A steps-only file
  is read, not applied. Each gene set written or replaced is also a
  `TRUE`/`FALSE` column `cellucid_gene_set__<name>` of the feature metadata
  of the gene axis the data mapping reads (a Seurat assay's, `rowData()` of
  the main experiment or of an `altExp()`), `TRUE` on the set's genes, as
  cellucid-python and the web app's "Save AnnData with this state..." write
  it (pinned by the shared golden history's `apply.sequence[].var_columns`);
  missing genes, a mapping reading the raw matrix and a column of that name
  the object owns are warned. A plan that would write, rename or replace the
  cell-metadata column `cellucid_original_position` is refused, dry runs
  included. Removed, no aliases: `expected_dataset_id`,
  `apply_labels`, `add_highlights`, `highlights_prefix`,
  `add_user_defined_fields`, `user_defined_prefix`,
  `include_deleted_user_defined_fields`, `store_metadata`, `return_summary`
  (of `cellucid_apply_session()`). `cellucid_pull()` lost
  `expected_dataset_id`. The packaged demo session is
  `extdata/history_golden.cellucid-session`.
- **The viewer's history, live.** `cellucid_turn(viewer, label, code)` makes
  every action `code` sends one turn of the history and returns the turn's
  id, for `cellucid_history_turn_undo()` or a note on the turn. Follow every
  change with
  `cellucid_last_history_change()` and `cellucid_on(viewer,
  "history_changed", ...)`. `cellucid_load_session()` and
  `cellucid_load_session_result()` push a session file into the viewer.
  New events: `action_result`, `history_changed` and `session_loaded`;
  `command_error` names the request it refused.
- **`cellucid_mcp()`**, run by an agent app as
  `Rscript -e "cellucid::cellucid_mcp()"`: a local MCP server for Claude
  Code, Claude Desktop and Cursor, with the web connector's tools,
  permissions, confirm cards (also asked in the app as an MCP elicitation),
  turns, Stop and Undo turn. It drives R's and Python's viewers alike, and
  Python's `cellucid mcp` drives R's. The protocol is written in the package:
  the handshake of 2024-11-05 to 2025-11-25 and the per-request 2026-07-28
  protocol, with elicitation, progress, cancellation and
  `tools/list_changed`. Every R viewer with a server is published to agent
  apps on this computer through a record readable by your user only (a folder
  others own or may write in is refused), and `CELLUCID_AGENTS=off` publishes
  nothing. A dataset's "Agents see names and summaries only" switch, set in
  the viewer, applies to these agents too. New Imports: promises (already
  httpuv's). Its arguments: `open` serves a file in a viewer of its own
  (`port`, `browser`, `web_source_url`, `web_cache_dir` as for
  `cellucid_view()`), `viewer` names a published viewer by viewer id,
  dataset id or dataset name, and with neither it drives the most recently
  opened viewer whose page is open. A page reloaded or closed mid-call fails
  the call with one sentence, and a page back makes the app list the tools
  again. A notebook widget's viewer is not published. A Windows console is
  refused: an agent app starts the server over a pipe.
- `cellucid_serve()` and `cellucid_view()` refuse a `.cellucid-session` or
  `.cellucid-library` file, saying what it is.
- A session's cell sets are read in the three encodings the web app writes;
  the reader of category-relative encodings, which nothing writes, is
  removed.
- Downloading the web app for a viewer accepts a JavaScript file served as
  `text/javascript` (RFC 9239) where the asset list says
  `application/javascript`; a file of another kind is still refused.
- Releasing: publishing a GitHub Release `vX.Y.Z` is the one way to release
  the package. `.github/workflows/release.yaml`, now "Publish to CRAN", runs
  on a published release only (its manual run is gone) and refuses, before
  building, a draft, a release marked pre-release, a `DESCRIPTION` version
  other than `X.Y.Z`, a tag other than `v` and that version, and a commit not
  on `main`; then it checks, attaches and submits the tarball as before. A
  failed run is re-run with "Re-run jobs" (`publishing.md`).

## In the web app

The viewer of this release adds History v2 (a branch graph with branches,
trackers, recording on and off, groups and squash, notes and bookmarks,
Recently deleted, replay on another dataset), Share... and Load... of
format-4 session files in State > Session, "Save AnnData with this
state..." in the History menu (a new `.zarr.zip` of the dataset with the
state shown applied by the same rules as `cellucid_apply_session()`), My
datasets with export and import, and agents: a chat bar with your own key, a connector for Claude and
ChatGPT, and the per-dataset switch "Agents see names and summaries only".
A replay's report offers other names for fields, categories and genes not
found ("Use instead", "Replay again with these names"), results left out of
a shared file are computed again when it opens, and the history warns when
it holds more memory than the device has room for. Agents can also play the
guided tour, make the report, export figures, cells and fields, manage what
the browser keeps, and vote on and suggest community annotation labels.
See the documentation and the cellucid-python changelog for the details.

# cellucid 0.9.6

Version 0.9.6 is the CRAN submission release.

It is a patch release on 0.9.5 and stays in step with the Python package of
the same version. Stores written by `cellucid_prepare()` in 0.9.5 and sessions
saved with 0.9.5 are read as before. The viewer protocol is schema 2: use this
release with the web app of this release.

## The R package

- **`cellucid_subsample_cells()`** (native code) gives exactly the cells the
  viewer and Python keep.
- **`cellucid_prepare()` stores cells in priority order** (`priority_seed`,
  default 0; `NULL` keeps the order given), with the same record, position
  column and prefix counts as Python, byte for byte; smaller chunks for
  shares; a cells × genes matrix is ordered without holding the whole block
  (36M values: peak 1.31 GB, 1.22 GB without priority order).
- **`cellucid_restore_original_order()`** for data frames, Seurat (every
  assay, reduction, graph and identity reordered) and SingleCellExperiment.
- **The notebook widget from R.** In Jupyter (IRkernel), `cellucid_view()`
  shows the same widget as Python over IRkernel's comm, with no server or
  port (needs IRkernel, IRdisplay and anywidget). Tested end to end in
  JupyterLab; Google Colab's R runtime is supported but not yet tested
  there. `cellucid_serve()` and `pane = "browser"` say in a notebook why they
  cannot answer there.
- **Served in priority order.** `cellucid_view()` / `cellucid_serve()` serve
  an object in priority order without copying its matrix or neighbour graph
  (1M cells: 0.3 GB beyond the object, was 1.4 GB); the same object served by
  R and Python holds identical bytes. A file given by path is served as it
  is, with a strong `ETag`; for a file not in priority order R says, as it
  starts, how to send only a share's cells.
- **BPCells and DelayedArray matrices are never read whole.** A matrix held
  outside R's memory (a Seurat v5 BPCells layer, an HDF5Array or other
  DelayedArray assay, or one passed as `gene_expression`) is read a block at a
  time by `cellucid_prepare()`, `cellucid_view()` and `cellucid_serve()`: one
  pass counts its values, then each gene's values are placed by cell and sorted
  in C, a quarter of the machine's memory at a time, each pass reading only its
  genes' rows. It is no longer joined into one `dgCMatrix` of doubles (12-24
  bytes a value), and the 2^31 - 1 value limit is gone. The bytes written are
  those of the same values in memory. Measured (16 GB Apple M1 Pro): 1M cells
  x 2.2 billion values from BPCells prepared in 413 s at 7.3 GB peak (refused
  before; joined it needed 26-53 GB); 1M cells x 400M values served (laid out)
  in 18 s; 100k cells x 220M values prepared in 64 s at 3.8 GB peak, against
  67 s and 6.1 GB for the same matrix converted to a dgCMatrix first.
  A DelayedArray block's conversion no longer warns
  "omp_pause_resource() failed to relinquish resources" once a block (SparseArray
  on an OpenMP runtime that cannot release its threads, such as macOS's): that
  one warning is muffled there, and every other warning is still raised.
  A DelayedArray block is read with `read_block()` (sparse when its backend
  is) rather than coerced to a `dgCMatrix`: up to Bioconductor 3.18 (R 4.1 to
  4.3) that coercion dropped `NA` and `NaN`, so they were written as zeros
  instead of refused.
- **Genes asked while the viewer lays a matrix out.** An object's expression
  matrices are laid out for the viewer on every core, and a gene asked before
  the layout ends is made from its progress: the values already laid out are
  copied and those of the cells not reached yet found by one scan of the
  matrix's values in its own order, on every core, the layout standing
  meanwhile; a burst of genes waits for the layout once that is sooner. `cellucid_prepare()` makes a genes x cells matrix's gene-major
  arrays on every core as well. In a notebook, the widget has at most two of
  the page's fetch messages with IRkernel at once, the rest gathered into the
  next (IRkernel spends about 1.7 ms on each message it reads).
- **Round trip.** `cellucid_pull(object, viewer or paths)` (Seurat
  `meta.data`, SCE/SE `colData`, data frames; `NA` for cells not open),
  `cellucid_push()`, `cellucid_push_gene_set()`, `cellucid_subset()`,
  `cellucid_gene_sets()`, `cellucid_results()`, `cellucid_write_gmt()` /
  `cellucid_read_gmt()`; `cellucid_highlight()` takes obs names and a label.
- **One session API.** `cellucid_request_session()`, then
  `cellucid_session_result()` or `cellucid_apply_session(object, request)`;
  `cellucid_session()` is removed and `cellucid_apply_session()` takes a file
  or a request, no longer a viewer. Sessions saved on a share apply.
- `cellucid_set_data_mapping()` no longer waits or takes `timeout`; a
  `command_error` no hook handles is a warning.
- **The same subsample size on every platform.** `cellucid_subsample_cells(
  fraction = )` takes the fraction's shortest decimal in exact integer
  arithmetic, a tie going to the even digit as in Python and the web app: on
  Windows the C library wrote 2^-25 as `...313`, and 2^-25 of 50,331,648
  cells came out 2 cells instead of 1.
- **`cellucid_stop()` delivers the close after the freeze.** It sends the open
  page `freeze`, then closes the channel (1001) and waits for the page to
  answer the close (at most a second) before stopping the server: httpuv
  writes the second of two frames sent together a turn of its loop later,
  and stopping at once dropped the close.
- **A refused WebSocket is closed.** httpuv opens a WebSocket whose handshake
  the server refused (a page on another origin, 403; a `Host` that does not
  name the server, 421) all the same; the server now closes it at once (1008)
  instead of leaving it open until it failed to say hello.
- Windows: the package builds again with R 4.1's 32-bit toolchain (Rtools 4.0
  does not declare `GetTickCount64()`; the layouts' clock is the performance
  counter), with no compiler warning.

## In the web app

- **Subsampled reading.** *Cells to open* (Session › Open data) opens a seeded
  share of a dataset's cells by the file's default method or by `uniform`,
  `blocks`, `stratified` or `spatial`; only the parts of the file that hold
  them are read, and the view, filters, analyses, genes, exports and sessions
  work on those cells. A link can name a share (`cells=10%&sample=blocks&seed=0`).
  `cellucid_subsample_cells()` gives exactly the same cells.
- **Too large for the device, progressive opening, Apply in place.** A
  dataset that would not fit opens fitted, with **Open all cells**; slow
  openings show a first view at once and grow in place; a larger or smaller
  share is applied in place, reading only the cells added.
- **R objects and MatrixMarket folders** opened in the browser are read whole
  once, converted in priority order and kept in the browser's storage, so they
  open again in milliseconds; every kept file is checked as it is read.
- **Analyses on a share** match an independent reference on the same cells;
  results go *Stale* when the open cells change and are not saved, exported or
  sent to R until run again.
- **Up to 100 million cells in 3D** on a 16 GB laptop, with levels of detail
  that are uniform samples of the cells.
- **Genes and failures said plainly.** *Gray for zero* in the legend, *Genes
  read for* under Cells to open, and a gene chosen is shown or the app says
  why; a data host's failure names the provider, its host and the status,
  with **Try again** where it may help.

# cellucid 0.9.5

Version 0.9.5 is the CRAN submission release.

It is a patch release on 0.9.4 and stays in step with the Python package of
the same version. It breaks with every earlier release: the prepared export
format is removed, and `cellucid_prepare()` writes a Cellucid-ready AnnData
instead, a standard AnnData Zarr store that Python's anndata reads and the
Cellucid web app opens directly. The contract is the documentation's
"Cellucid-ready AnnData specification" page. An export written by 0.9.4 is not
read by the matching web app; prepare the data again with this release.

## Writing a Cellucid-ready AnnData (breaking)

- `cellucid_prepare(out_path, format, embeddings, dataset_id, dataset_name,
  ...)` takes a required `format` right after `out_path`: `"zarr"` (a `.zarr`
  folder of small files, the way to publish on GitHub: no file over 100 MB, no
  LFS) or `"zarr.zip"` (the same store as one `.zarr.zip` file, for hosts that
  serve single files or a local file). It has no default: a missing `format`
  is reported with the other missing required arguments, naming both choices.
  A value is matched as `match.arg()` matches, anything else is refused with
  the choices, and `out_path` must end in the suffix of `format`. `zarr_format
  = 3` writes Zarr v3 (`zarr.json` documents, `c/...` chunk keys, consolidated
  metadata in the root `zarr.json`); Zarr v2 is the default. There is no
  `"h5ad"`: the package writes no HDF5 and depends on no HDF5 library;
  Python's `cellucid.prepare(format = "h5ad")` writes `.h5ad`.
- The store holds `X` gene-major (CSC from a sparse matrix, stored zeros
  dropped) in float32, or float64 when the genes kept hold integer counts
  above 2^24; `obs` and `var` with anndata's own encodings and pandas'
  category code widths (int8 below 127 categories, int16 below 32,767, int32
  beyond); every embedding whole in `obsm` (float32, unnormalized);
  `obsp/connectivities` as symmetric float32 CSR; and `uns[["cellucid"]]`
  (`schema`, `id`, `name`, `description`, `citation`, `source`, and only the
  launch `defaults` that were set) with `uns[["title"]]`. Chunks hold 2^16
  values for `X` and 2^18 for arrays read whole, deflated at level 6 after a
  byte shuffle on every multi-byte array except the expression values: in Zarr
  v2 numcodecs' `shuffle` filter and zlib; in Zarr v3 only codecs of the v3
  specification, the `blosc` codec (byte shuffle, zlib level 6, the frames
  c-blosc writes) on shuffled arrays and gzip on the others. A `.zarr.zip`
  stores chunks first in grid order and the consolidated metadata last, with
  fixed timestamps and permissions. The store is written beside its target and
  moved into place when complete.
- R and Python write identical Zarr v2 and v3 stores from the same input:
  every metadata document (encodings, dtypes, shapes, chunks, codecs, fill
  values) and every value, checked by `tests/testthat/test-cross-writer.R`.
  Every `NaN` and `NA` in a float array is written as numpy's canonical quiet
  NaN on every platform (R's `NaN` has its sign bit set on x86-64).
- New arguments: `out_path`, `format`, `zarr_format`, `citation_text`,
  `citation_url`, `latent_key` (the `embeddings` entry the outlier filter
  measures in) and `normalize_embeddings`. Choose the columns an embedding
  opens on with `embedding_defaults = list("3d" = list(key = "X_pca", columns
  = 1:3))`.
- Removed with the export format: the `latent_space` (name an embedding with
  `latent_key`), `out_dir`, `obs_categorical_dtype`, `var_quantization`,
  `obs_continuous_quantization`, `compression`, `var_gene_name_columns`,
  `source_citation` (now `citation_text` and `citation_url`), `created_at` and
  `centroid_outlier_quantile` arguments, `cellucid_embedding()` (store the
  whole array), and the cross-process export lock and transaction journal.

## Seurat and Bioconductor objects, several modalities

- `cellucid_prepare(data = )` writes a Seurat (Seurat 3 and later,
  SeuratObject 4 and 5, `Assay` and `Assay5`, split layers joined on a copy),
  SingleCellExperiment (1.8 and later) or SummarizedExperiment object
  directly; every other argument defaults to what the object holds and an
  argument given replaces that component. The mapping is the web app's for the
  same object saved with `saveRDS()`: `meta.data`/`colData` as `obs`,
  reductions/`reducedDims` as `X_<name>`, `data` (else `counts`) or
  `logcounts` (else `counts`) as `X` (one matrix per modality, as the
  Cellucid-ready AnnData contract has it), `var.features` as
  `highly_variable`, the first symmetric SNN graph or `colPair` as
  `obsp/connectivities`, and `misc$cellucid`/`metadata(x)$cellucid` as the
  identity and launch defaults. BPCells and DelayedArray (HDF5Array) matrices
  are read in blocks of whole cells, joined once at the end. Objects too old
  to read are refused with the line that updates them, and one without an
  embedding with the line that computes one. An assay holding the dataset's
  cells in another order than its `meta.data` (a Seurat v5 layer with its own
  cells, `modalities = "ADT"`) is reordered to the dataset's cells, and an
  assay holding only some of the cells is refused with what to do (keep those
  cells, or write it beside an assay of every cell as a MuData). An AnnData's
  `X` is checked to have a row per cell before anything is written. A factor
  with an `NA` level (`addNA()`) is written with those cells missing (code
  -1), as pandas stores them; they get no outlier quantile.
- Several Seurat assays, or a SingleCellExperiment's alternative experiments,
  are written as a Cellucid-ready MuData (Zarr only; R writes no HDF5):
  MuData's own layout, which `mudata.read_zarr()` reads and Python's
  `cellucid.prepare(data = mdata)` writes identically
  (`tests/testthat/test-cross-writer.R`). `modalities` chooses the modalities
  and their order (the active assay or main experiment first by default), and
  the launch default `modality` names the one the viewer opens with; without
  it the viewer applies the shared rule (the first named `rna`, `gex`, `gene
  expression`, `gene_expression` or `sct`, else the first), pinned by the
  data-mapping conformance fixture.
- New arguments `data`, `modalities` and `modality`; `embeddings`,
  `dataset_id` and `dataset_name` are now required only without `data`.
- `gene_identifiers` keeps genes of an object's genes x cells matrix without
  copying its rows: the kept genes are positions the gene-major arrays read
  in place.
- The viewer opens an object on its UMAP under the R world's conventional
  names, as on Scanpy's `X_umap`: after `X_umap_1d`/`X_umap_2d`/`X_umap_3d`,
  the first 1- to 3-column one of `X_umap` (a Seurat `umap` reduction),
  `X_UMAP` (a SingleCellExperiment `UMAP` reducedDim), `X_wnn.umap` and
  `X_ref.umap` is the default UMAP. The rule is the shared data-mapping
  conformance fixture's (version 6).

## The viewer

- `cellucid_view()` opens a Seurat, SingleCellExperiment or
  SummarizedExperiment object, or any file the web app reads (registry:
  `inst/input_formats.json`), in the RStudio or Positron Viewer pane or a
  browser; `cellucid_serve()` serves it for a browser elsewhere (an SSH
  tunnel) and prints the address. No Python is involved. An object is served
  from memory as the Cellucid-ready AnnData (or MuData) Zarr store
  `cellucid_prepare()` writes for it -- the same plan, the same arrays and
  values, every array read back from both and compared (dtype, shape, values)
  by `tests/testthat/test-served-store.R`, with zarr-python where one is found
  -- uncompressed: it declares no compressor and no shuffle, since the browser
  reads it over loopback, where deflating each chunk on R's one main thread
  cost far more than the bytes it saved (files `cellucid_prepare()` writes
  stay compressed). Every chunk is made only when the viewer asks for it, in
  one native pass from its R vector (`src/chunks.c`); only chunks made by a
  transform (text, a dense matrix block) are kept, in a least-recently-used
  cache of `cache_bytes`. An object's sparse expression matrices are laid out
  gene-major once, in the background from the moment the viewer starts (`X`
  first), on a native thread into memory-mapped scratch files removed as
  soon as they are mapped (never a copy in R's memory), as cellucid-python
  does; a gene asked for before is made at once from the object (the block
  of genes around it, one binary search a cell), exactly as it will be laid
  out. A served expression matrix's entries come in chunks of 2^18 values
  (1 MiB: a browser's `fetch` over loopback moves 46-79 MB/s in requests of
  64 KiB, 300-510 MB/s in requests of 1 MiB), and a `Range` of a chunk is
  made alone. On 1,000,000 cells (100 values a cell) in Chromium, the first
  picture takes 1.7 s, 100 genes 0.63 s and a differential expression read
  of 50 genes over 133,334 cells 0.17 s (blocks of genes made within
  `cache_bytes` had taken 2.6 s, 38.8 s and 30.5 s). A 100,000-cell object
  starts serving in 0.7 s and shows its first view (cells, annotations, UMAP)
  in 0.1 s more. A file is
  served from disk with exact `Range` answers. The server (`httpuv`) runs in
  the background, so the console stays usable; it mirrors the routes of
  cellucid-python's server, refuses a `Host` that does not name it (421;
  checked unless the server is bound to a non-loopback IP address with no
  `allowed_hosts`, cellucid-python's rule) and a WebSocket from a foreign
  origin (403), and defaults to port 8765 (the address Google sign-in
  accepts), taking the next free port when that one is taken. Several viewers
  run at once; `cellucid_stop()`, `cellucid_stop_all()` and
  `cellucid_viewers()` manage them.
- The viewer channel: one WebSocket per viewer (`/_cellucid/ws`, the protocol
  of the web app's `data/viewer-protocol.js`, pinned by the shared fixture
  `tests/testthat/fixtures/viewer_protocol.json`). `cellucid_highlight()`,
  `cellucid_clear_highlights()`, `cellucid_color_by()`,
  `cellucid_set_visibility()`, `cellucid_reset_camera()`, `cellucid_freeze()`
  and `cellucid_set_data_mapping()` send commands (queued until the viewer
  connects); a highlight and a visibility change are binary frames, their
  cells little-endian `uint32` checked and written in one native pass (two
  million cells: 7.6 MB in 0.02 s, highlighted in the viewer 0.21 s after the
  call; ten million: 0.07 s and 0.93 s); `cellucid_on()`, `cellucid_off()`, `cellucid_selection()` (cell
  positions and names), `cellucid_wait_ready()` and `cellucid_session()` (the
  viewer's `.cellucid-session` file) take its events. Events (`hover`,
  `click`, `selection`) name cells by 1-based position, as every command
  does, and also by name (`cell_name`) for an object. An event naming a cell the
  dataset does not have closes the channel (1008), and nothing a socket sends
  after such a close is delivered. The channel is pinged every 20 seconds (a
  `ping` command, whose `pong` is not delivered to hooks), as
  cellucid-python's is, so a proxy does not close an idle viewer; a socket
  that does not say hello within 60 seconds is closed.
- The viewer id, token and request ids come from the operating system's
  cryptographic generator (`src/random.c`), not R's: `set.seed()` cannot
  predict the token, and opening a viewer leaves the caller's `.Random.seed`
  as it was.
- The viewer page's address names a host program's viewer with `?viewer=true`,
  and `/_cellucid/info` names the program serving it: `runtime = "r"`, which
  the viewer shows in its Dataset box ("R"). Both are pinned in the shared
  viewer protocol file (`page`, `runtimes`).
- The viewer's server answers the store's bytes with `Content-Encoding:
  identity` (httpuv would gzip a range on R's main thread and drop its
  `Content-Length`), and a `HEAD` carries the `Content-Length` of
  its `GET`, so the web app's Zarr reader can size every key. A range of a
  file is never read whole into memory: the whole file is sent from disk by
  httpuv, a range of up to 8 MiB is read with one seek, and a longer one is
  copied block by block (8 MiB) into a file httpuv sends and removes. A
  request path that is not valid percent-encoded UTF-8 (`%00`, a lone `%ff`)
  is answered 400, and a failed request 500 without the internal error text
  (which goes to the R console).
- A dense genes x cells assay is served through its transpose block by block,
  never transposed whole, and so is a sparse one (its gene-major blocks,
  above); `cellucid_view()`'s help page states what `cache_bytes` holds.
- The web app is downloaded at first use from www.cellucid.com into
  `file.path(tools::R_user_dir("cellucid", "cache"), "web")`, checked file by
  file against the build's SHA-256 inventory with cellucid-python's rules, and
  replaced whole when a new build is published; a verified cache is used
  offline. `cellucid_clear_web_cache()` removes it. Each file is asked for as
  stored (`Accept-Encoding: identity`, as Python's web cache asks), since
  www.cellucid.com (GitHub Pages) gzips what a client accepts gzipped and
  declares the compressed length; the test source gzips as GitHub Pages does.
- New `cellucid_debug_connection(viewer, timeout = 5)`: the report
  cellucid-python's `AnnDataViewer.debug_connection()` gives (the server's
  routes, the store's `HEAD` and first byte, a ping and a debug snapshot over
  the channel), which makes the `pong` and `debug_snapshot` events reachable.
- `quiet = TRUE` also silences the notes reading an object writes; cells given
  as `TRUE`/`FALSE` with an `NA` are refused; a `misc$cellucid$citation` that
  is not a list is named in a note; a failed `VariableFeatures()` is named
  rather than giving no highly variable genes; a Zarr store whose own identity
  cannot be read is said to be damaged; a server that does not stop is
  reported and kept listed, so `cellucid_stop_all()` can try again; unloading
  the package stops every viewer first; web assets are checked against their
  digest once and then served from disk.

## Sessions

- New `cellucid_apply_session(data, session, ...)`: a saved session's
  highlight groups (one `TRUE`/`FALSE` column each) and user-defined fields
  (factors, or copied numbers) added to a Seurat object's `meta.data` or an
  experiment's `colData`, with the session recorded in `misc$cellucid$session`
  / `metadata(x)$cellucid$session`. `session` is a `.cellucid-session` file or
  an open viewer. The rules are cellucid-python's
  `apply_cellucid_session_to_anndata()`, check for check, held to the same
  web-app-written session file (`tests/testthat/fixtures/session_format2/`,
  pinned in both repositories); every check runs before anything is added. The
  session reader is new: the file's framing and manifest, the data mapping
  rules (now run against the shared conformance fixture's `dataMapping` and
  `mappingModality` cases in R too), and the cell-set and user-defined-code
  codecs in C.

## Outlier quantiles, vector fields and graphs

- With `latent_key`, `cellucid_prepare()` stores the latent-space outlier
  quantiles of every factor, character or logical `obs` column in
  `uns[["cellucid"]][["outlierQuantiles"]]` (`latentKey`, `minPoints` =
  `centroid_min_points` or 10, and one float32 array per field, `NaN` where a
  cell has none), so the viewer filters outliers without reading the latent
  space. They are computed in C (`src/latent_outliers.c`), byte for byte what
  Python's `prepare()` stores and what the web app computes, pinned by the
  shared conformance fixture
  (`tests/testthat/fixtures/latent_outlier_conformance.json`).
- Vector fields follow scVelo's naming on any embedding: `vector_fields =
  list(velocity = list(X_umap = v))` writes `obsm[["velocity_umap"]]`, with as
  many columns as the embedding. An embedding whose key does not start with
  `X_`, an array of another width, two fields that would share a key, and a
  key the viewer would read as another field or as an embedding are refused
  before anything is written. `vector_field_default` names the field
  (`"velocity"`) and must be drawn on an embedding the viewer shows when the
  store opens.
- New `cellucid_transition_drift(transition_matrix, embedding,
  normalize_rows)` computes each cell's drift in an embedding's coordinates
  from a cell-cell transition matrix (CellRank's, for example): `T %*% E - E`,
  with the rows of `T` divided by their sums first when `normalize_rows =
  TRUE` (required; a zero row sum is refused). It is Python's
  `compute_transition_drift()` operation for operation, held to the same
  results by a shared conformance fixture
  (`tests/testthat/fixtures/transition_drift_conformance.json`). Store the
  result with `cellucid_prepare(vector_fields = list(T_fwd = list(X_umap =
  drift)))`.
- A sparse neighbour graph's stored zeros are dropped (`Matrix::drop0()`)
  instead of refusing the graph; duplicate coordinates are still refused.

## Windows R 4.1: 32-bit R and sessions that are not UTF-8

- On 32-bit R (the i386 half of Windows R 4.1) every gzip or zlib stream was
  refused with "inflate: a stream of 4 GiB or more is not one session frame",
  so `cellucid_apply_session()` and serving a `.zarr.zip` failed: the check
  cast `UINT32_MAX` to `R_xlen_t`, an `int` there, which made it -1. The
  inflater now hands zlib its input and output in pieces and refuses only an
  output longer than one R raw vector, and every size the C code takes from R
  is checked as a double before it is cast to `R_xlen_t` or multiplied.
- In an R session whose native encoding is not UTF-8 (Windows R before 4.2,
  which runs in CP1252; a C or ISO-8859 locale elsewhere), every message,
  warning and error reaches its handler with the package's exact text.
  `message()`, `warning()` and `stop()` re-encoded it to the code page, so the
  arrow (U+2192) of the note on repeated names ("1 duplicate name was made
  unique, e.g. ...") arrived as `<U+2192>`. They are raised as condition
  objects now, and a test refuses a `stop()`, `warning()` or `message()` given
  strings.
- In such a session, a Seurat or Bioconductor object whose `meta.data`,
  `colData` or `rowData` column is named with a character the code page lacks
  (Greek, Cyrillic, CJK) was refused as "obs column '<U+03B1>' contains a
  character that cannot appear in a file name": the columns were rebuilt with
  `as.data.frame()`, which turns each name into a symbol. They are kept as
  they are.
- Repeated cell and gene names holding such characters are made unique
  without a warning per name, and a name is never taken for its escape
  (`"<U+03B1>"`): the hash set behind the rule, the served store's keys and
  its chunk cache are keyed by each string's UTF-8 bytes.
- Text read from bytes is read as UTF-8 in every session: a Zarr folder's
  metadata (and the consolidated metadata the server makes for a folder that
  has none), the browser's text frames and binary event headers, the web
  build's inventory and the server's own JSON. In a code-page session a
  non-ASCII character was misread as several.
- The viewer's server decodes a request path as UTF-8, so a served key with a
  non-ASCII character (an `obs` column named in Greek) is found, and hands
  httpuv each file it streams by its UTF-8 name, the only one httpuv opens: in
  a CP1252 session a store, the web app's cache or the temporary folder under
  a user folder with an accented name was not served.
- Categories of a character column are sorted in code-point order whatever
  each string's encoding, as Python's `prepare()` sorts them: a radix sort
  compares the bytes as R holds them, so a latin1 e-acute sorted after a UTF-8
  Greek letter.
- A dataset name, description, source or citation held in latin1 or the
  native encoding is reported in the usual words when it carries a defect (a
  trailing space); reading its characters as UTF-8 failed with "missing value
  where TRUE/FALSE needed". Data-mapping keys and session field names are
  compared by code point the same way.

## Package

- The installation guide names the one extra step SeuratObject needs on R
  older than 4.4: install Matrix 1.6-5 first
  (`remotes::install_version("Matrix", version = "1.6-5")`), because
  SeuratObject 5 needs Matrix 1.6-4 or later and CRAN's current Matrix needs
  R 4.4. CI checks those R versions the same way.
- Messages, errors, warnings and `print()` of a viewer make a count and its
  noun agree ("1 row", "1 cell"; they read "1 rows"): the noun comes from one
  internal `.plural()`, and a test refuses `paste0(n, " cells")`.
- `cellucid_serve()` and `cellucid_view()` refuse an h5Seurat file before
  serving it, in the words of the input format registry the web app and
  `cellucid serve` use too: it is not opened; read it with
  `SeuratDisk::LoadH5Seurat()` and pass the object, or save it with
  `saveRDS()`.
- R 4.1 or later (was 4.3); CI checks R 4.1, 4.2, 4.3, oldrel, release and
  devel, with every suggested package (SeuratObject 4 and 5,
  SingleCellExperiment) and with none. Imports `httpuv`, `later` and `curl`
  for the viewer; SeuratObject, SingleCellExperiment, SummarizedExperiment,
  S4Vectors and DelayedArray are suggested, with testthat 3.2.1 or later (the
  tests mock values as well as functions). A test refuses base R functions and
  syntax newer than R 4.1, checked against the list R 4.1 itself writes; the
  fixture SHA-256 pins and the citation checks use the package's own SHA-256,
  so they run on every supported R (they needed `tools::sha256sum()`, which
  is R 4.5, or a system tool).
- `SystemRequirements: zlib`: the C code links zlib (`-lz`) for the session
  codecs and the exact gzip and zlib inflater, so a Linux build needs its
  development headers (`zlib1g-dev`, `zlib-devel`).
- Large matrices are checked in one pass in C (finiteness, float32 range,
  stored zeros, the dtype rule) with no copy, and a SHA-256 in C checks the
  web build on every R version.
- Every JSON text read (Zarr metadata, the viewer's frames, the web build's
  inventory) is parsed with `jsonlite::parse_json()`, which never reads a
  string that names a file or URL from that file or URL, as `fromJSON()` does.
- A gzip or zlib stream is inflated in C with exact checks; base R's
  `memDecompress()` never returns on two concatenated gzip members, which a
  damaged Zarr chunk could hold.
- Tutorials on docs.cellucid.com, each run end to end: "Explore a Seurat
  object from R" (`cellucid_view()`, commands from the console, the cells you
  select back in R), "Share a SingleCellExperiment on GitHub"
  (`cellucid_prepare(data = sce, format = "zarr")`) and "CITE-seq: a Seurat
  object with RNA and ADT" (a Cellucid-ready MuData, the Modality select).

# cellucid 0.9.4

Version 0.9.4 is the CRAN submission release.

It is a patch release on 0.9.3, which was tagged but never published on CRAN,
stays in step with the Python package of the same version, and writes the same
export format 3: an export written by 0.9.3 does not need to be prepared again.

- `cellucid_prepare(obs_categorical_dtype = "uint32")` writes 32-bit
  categorical codes (`obs/{index}.codes.u32`, missing marker 4,294,967,295)
  for fields with more than 65,535 categories, such as barcodes, byte for byte
  as the Python package does. The error for a field that does not fit the
  declared width names the width that does, for example
  "Pass obs_categorical_dtype='uint32'.".
- `cellucid_prepare()` exports a categorical field where no category reaches
  `centroid_min_points` with quantized outlier quantiles, publishing every
  quantile as missing, instead of stopping partway through the export with
  "it contains no finite values".

# cellucid 0.9.3

Version 0.9.3 is the CRAN submission release.

It carries every change of 0.9.2, which was tagged but never published on
CRAN, and stays in step with the Python package of the same version.

`cellucid_prepare()` writes export format v3, byte for byte the format
`cellucid-python` 0.9.3 writes. Opening a v3 export needs the Cellucid web app
release that reads format 3, and an export written by an earlier version is
regenerated with this one.

- Exports any number of embeddings per dimension. `cellucid_prepare()` takes a
  required named list `embeddings` in place of `X_umap_1d`, `X_umap_2d`, and
  `X_umap_3d`, which are removed. A matrix of 1, 2, or 3 columns is plotted
  whole; the new `cellucid_embedding(values, columns, label, source_key)`
  names the 1-based columns of a wider array (published 0-based), a display
  label, and the key the array came from. `embedding_defaults`, for example
  `c("2d" = "X_umap")`, names the embedding each dimension shows first and is
  required for every dimension with more than one.
- Writes each embedding to `embeddings/{index}_{dimension}d.bin`, ordered by
  dimension and then by id in code-point order, and declares it in
  `dataset_identity.json` version 3 under `embeddings.items` with its label,
  dimension, file, and `source` (`key` and 0-based `columns`). The export root
  holds no payload files any more; `embeddings/` is reconciled against its
  declaration as a sixth surface, beside `obs/`, `var/`, `connectivity/`,
  `vectors/`, and the root itself.
- Binds vector fields to embedding ids: `vector_fields = list(velocity =
  list(X_umap = v2, X_umap_3d = v3))`. Each array is scaled by the factor of
  the embedding it is bound to and written to
  `vectors/{fieldIndex}_{embeddingIndex}.bin`; a field id no longer has to
  follow the `<field>_umap_<n>d` grammar.
- Writes `obs_manifest.json` and `var_manifest.json` as `compact_v2`. A
  categorical entry's centroids are keyed by embedding id
  (`centroidsByEmbedding`).
- Names genes from an integer `var` column by its digits (`7157L` is
  `"7157"`), as the Cellucid viewer and `cellucid serve` read the same column;
  a double or logical column is refused with the reason ("floating-point
  values cannot name genes").
- Computes latent-space outlier quantiles in C with the rule `cellucid-python`
  and the web app share, bit for bit: values rounded to float32 and widened to
  double, centroid sums in cell order, squared distances in column order, the
  share as float32. A shared conformance fixture
  (`tests/testthat/fixtures/latent_outlier_conformance.json`) pins the bytes
  in all three. Quantiles can differ from earlier releases only where two
  distances tie within rounding.
- Accepts repeated gene names. They are made unique the way anndata's
  `var_names_make_unique()` does it -- the second `MATR3` becomes `MATR3-1` --
  with the one rule `cellucid-python` and the web app also apply, and an empty
  name (`NA` or `""`) takes the row's name from `rownames(var)`. What changed
  is reported as a message and recorded in `var_manifest.json`.
  `gene_identifiers` matches the names after they were made unique.
- Adds `var_gene_name_columns`: further `var` columns the viewer offers as
  alternative gene names. Each one is stored as metadata in the new
  `gene_name_sets` list of `var_manifest.json`, which replaces the manifest's
  `var_gene_id_column` key (the argument of that name is unchanged and selects
  the default set); no expression payload is duplicated.

- `citation("cellucid")` now lists the Cellucid preprint first, followed by
  the package manual entry. The preprint is a placeholder until it is published
  on bioRxiv and says so; the entry is generated from the `preferred-citation`
  block that every Cellucid repository's `CITATION.cff` shares.

# cellucid 0.9.2

Version 0.9.2 is the CRAN submission release.

It carries the changes prepared as 0.9.1, which was never published; the number
moved to stay in step with the Python package.

- Exposes `cellucid_prepare()` as the single public preparation entry point.
- Publishes complete prepared generations transactionally and rejects mixed or
  stale output directories.
- Serializes independent R and Python exporters with one persistent exact-target
  lock and recovers ownership after process death without leaking native
  handles.
- Requires exact embedding, observation, gene, vector, and categorical-storage
  inputs before writing any scientific artifact.
- Preserves symmetric positive connectivity weights, exact empty graphs, and
  portable edge-index widths without reinterpreting the source graph.
- Rejects non-finite or non-representable float32 values and duplicate
  scientific identifiers before publication, and says how many and where: a
  refusal counts the offending values and names the first few positions.
  `cellucid-python` reports the same counts; each writer prints positions the
  way its own caller indexes them, and names the categories its own language
  distinguishes.
- Applies the float32 range check wherever a value can enter, not only at the
  write. A quantized field below the smallest subnormal was published as zero,
  and an embedding or `latent_space` coordinate was normalized into range before
  the writer saw it; both ends are now one predicate, checked at the source, as
  in `cellucid-python`.
- Names every payload by its integer position on its axis, so no exported
  filename carries dataset content and every dataset has the same directory
  listing: `var/0.values.f32`, `obs/1.codes.u8`, `vectors/0_2d.bin`, and the
  fixed neutral names under `connectivity/`. Each manifest entry declares its
  own index as its first element -- `[index, name]` and
  `[index, name, minValue, maxValue]` for `fields`, `[index, key]` and
  `[index, key, minValue, maxValue]` for `_continuousFields`, and
  `[index, key, categories, dtype, sentinel, centroids]` with the two outlier
  bounds appended for `_categoricalFields` -- and `pathPattern`,
  `codesPathPattern`, and `outlierPathPattern` substitute `{index}`. The writer
  proves that each axis directory uses exactly the indices `0` through `N-1`
  once each, and that the directory holds exactly the payloads its manifest
  declares, before the generation is published. `obs/` is written by both
  manifest arrays, so those two arrays share one index space. The `cellucid`
  Python package writes the identical layout.
- Reconciles the export root against the generation as well, so the
  `points_<n>d.bin` coordinate payloads are covered by the same guarantee as
  every other payload. They are declared in `dataset_identity.json` rather than
  in an axis manifest, and their declared name and their written name are two
  independent spellings of the `compression` setting, so an export that
  published coordinates under a name the viewer is not told to fetch used to
  succeed and then fail in the browser with no points. The root must now hold
  exactly the manifests, the axis directories, and the declared point payloads
  the generation created, and nothing else.
- Writes `_obsSchemas` as a JSON object in every generation, including one that
  carries no observation field at all. `jsonlite` renders an empty *unnamed*
  list as `[]`, so an export prepared from an `obs` with no columns, or with
  `obs_keys = character(0)`, published `"_obsSchemas":[]`; the viewer calls
  `requireRecord()` on that field and refused the dataset with "expected an
  object", while the `cellucid` Python package wrote `{}` from the same input
  and loaded. The manifest validator now asserts the kind rather than only the
  contents, so the two writers cannot part company here again.
- Reads every manifest back out of the staging directory before publishing and
  requires it to parse to exactly the payload that was validated. Everything
  else in the package validates the payload in memory, which proves nothing
  about the file a reader opens: `jsonlite::toJSON()` writes a non-finite double
  as the JSON *string* `"NaN"`, and an empty unnamed list as an array, so a
  manifest could be well-formed, publish cleanly, and still say something other
  than what was checked. A file that is not valid UTF-8, does not parse, or
  disagrees with the payload on any node's kind, keys, length, or value now
  fails the export and names the node, instead of reaching a browser.
- Drops every filename rule from gene names, `obs_keys`, and vector field ids,
  because none of them names a file any more: no portable-ASCII restriction, no
  case-insensitive collision rule, and no Windows device-name rule. `HLA-DRB1/2`
  and `CON` are exported and recorded exactly as supplied. What each identifier
  must still be is what it is for -- a non-empty string, distinct within its
  axis, and text the viewer can draw exactly as stored. `dataset_id` names the
  export directory, so it alone remains a portable filename component.
- Accepts any vector field key matching `<field>_umap_<1|2|3>d`, with no
  character restriction on `<field>`, matching the `cellucid` Python package
  exactly. Vector fields are emitted in code-point order of their ids on both
  writers, so the same input receives the same payload index in either
  language.
- Requires every string the viewer prints verbatim to read on screen as the
  value it stores. A string category label, an exported gene name, an exported
  `obs` key, a vector field id, `dataset_name`, `dataset_description`,
  `source_name`, `source_url`, and `source_citation` are rejected when they
  carry a control character, one of the zero-width
  characters `U+200B`, `U+2060`, or `U+FEFF`, or leading or trailing whitespace
  of any kind including `U+00A0` NO-BREAK SPACE. An empty category label is
  rejected, as are two labels in one field that a whitespace-collapsing
  renderer draws identically. Nothing is trimmed: trimming would rewrite an
  annotation the caller never asked to change, and would merge `"Liver "` into
  a separate `"Liver"` category and move cells between them. The message names
  every offending label in the field at once and gives the one-line repair. The
  `cellucid` Python package enforces the identical rule.
- Extends the `dataset_name` check, which previously missed `U+0080`-`U+009F`
  control characters and non-ASCII whitespace such as `U+00A0`, to the same
  shared rule.
- Derives quantized continuous payloads and bounds from the viewer's exact
  float32 value domain.
- Publishes a constant continuous field instead of refusing the export. A gene
  expressed at one level in every exported cell -- very often zero, once an
  atlas is subset to one lineage -- and an `obs` column a subset flattened are
  both ordinary data, and quantization now has a named case for them:
  `minValue == maxValue` with every code `0`, which the viewer decodes back to
  the exact constant rather than to an approximation of it. `.quantize_continuous()`
  takes an explicit branch for it and never derives a scale, so nothing divides
  by `maxValue - minValue`. Native-double variation finer than float32
  resolution is one float32 value and is published the same way, in place of
  the previous rejection. The `cellucid` Python package implements the
  identical case.
- Documents the two scopes the gene identifier rules have here, matching the
  `cellucid` Python package. Being drawable is a property of a name the viewer
  shows, so it covers the genes `gene_identifiers` selects; a `var` row left out
  reaches no manifest and is not checked. Distinctness spans the whole `var`,
  because `gene_identifiers` addresses `var` rows by identifier and a repeated
  one names no single row.
- Records whatever `var_gene_id_column` selects faithfully. `cellucid_prepare()`
  performs no symbol lookup and ships no mapping, so the caller decides what a
  gene is called.
- Reports an identifier defect in the same words the `cellucid` Python package
  uses. Each axis passes one singular noun -- `Gene`, `Observation field`,
  `Vector field` -- so the sentence composes as
  `Gene identifier at position 1 ...` instead of
  `Gene identifiers identifier at position 1 ...`, and a repeated key is
  reported as `Gene key 'ACTB' is duplicated.` on both writers. The checks that
  speak about a whole axis add their own plural, so no caller has to guess
  which number a message will need.
- Prints a set of values in a message as a list the caller can paste back into
  the failing call, `c(0, 1, 1)` and `c("score", "n_counts")`. Without the
  boundary, `columns not in obs: a, b. Available columns: x` cannot be read
  back as two lists; with Python's `['score', 'n_counts']` the boundary was
  there but the syntax inside it belonged to the other language. Every message
  that shows a set of values now renders it through one function, so
  `gene_identifiers`, `var_gene_id_column`, and the payload-manifest checks no
  longer each print a set their own way. The `cellucid` Python package prints
  the same sets as Python lists, for the same reason.
- Requires `obs_categorical_dtype`, `dataset_name`, and `dataset_id` as
  arguments with no default, and reports one that was never supplied as
  missing. `cellucid_prepare(dataset_name = "Atlas")` answered
  `dataset_id must be exactly one string.`, which describes a value the caller
  never passed and names nothing valid; it now answers
  `cellucid_prepare() is missing 1 required argument:` followed by
  `dataset_id` and what a valid identifier is. Every argument left out is named
  in that one message, in signature order, as the `cellucid` Python package
  names them in one `TypeError`. A value that *was* supplied and is wrong still
  reports the value: `dataset_id = NULL` remains
  `dataset_id must be exactly one string.`
- Defaults `dataset_description` to `NULL` and accepts `NULL` and `""` alike,
  both publishing `""`. It previously defaulted to `""` and rejected `NULL`,
  while every other optional identity argument in the same signature --
  `source_name`, `source_url`, `source_citation`, `created_at` -- took `NULL`
  for "not supplied". The `cellucid` Python package accepts `None` and `""` for
  it and publishes `""` for both.
- Reports an empty or padded `dataset_id` against the rule it actually breaks.
  `dataset_id = ""` answered
  `dataset_id must be one non-empty string without leading or trailing whitespace.`,
  naming whitespace a value with no characters cannot carry; it now answers
  `dataset_id '' is not a portable identifier. Use 1-180 ASCII letters, numbers, '.', '_', or '-', beginning with a letter or number and not ending with '.'.`
  Nothing that was accepted before is rejected now, and nothing that was
  rejected before is accepted.
- Writes one fixed gzip header for every compressed payload, so a `.gz` no
  longer records the machine that produced it. `gzfile()` leaves the header to
  zlib, which stamps its own build platform into the member's `OS` field --
  `0x03` on Unix, `0x0b` on Windows -- so the same input produced different
  export bytes depending on the operating system the export ran on. Each member
  now carries the ten bytes the `cellucid` Python package writes: no filename,
  a Unix-epoch timestamp, the extra-flags value RFC 1952 assigns to the
  requested level, and `OS = 0xff`, the code that names no operating system.
  The compressed payload, the level the caller chose, and the CRC32 and length
  trailer are untouched, and the deflate stream is still whatever the zlib
  underneath R produces at that level -- so a member matches the Python
  package's byte for byte when the two zlib builds agree, and differs only in
  the deflate bytes when they do not.
- Adds `tests/testthat/test-gzip-header.R`, which reads the ten header bytes of
  a written `.gz` directly at every compression level, and proves the member
  still decompresses to exactly the payload that went in.
- Adds `tests/testthat/test-writer-parity.R`, which holds this writer to the
  half of the two-writer contract it owns: a `vectors/` payload is the scaled
  value rounded to float32 once and not the twice-rounded value, a
  `points_<dim>d.bin` payload is the normalized value rounded once, a
  categorical centroid is measured from coordinates that were never rounded,
  and a `dataset_identity.json` vector field entry carries its keys in the
  documented order.
- Refuses to select row numbers as gene identifiers. A `data.frame` always has
  row names, so the old `is.null(rownames(var))` test could not fire, and
  `var <- data.frame(symbol = c("CD8A", "MS4A1"))` with the default
  `var_gene_id_column = NULL` exported two genes named `"1"` and `"2"` -- the
  automatic row-name sequence that `rownames()` materializes. Those strings are
  unique and drawable, so every later check passed and nothing reported the
  substitution; a wet-lab user then searched `CD8A` in the viewer and found
  nothing. `cellucid_prepare()` now distinguishes automatic row names from row
  names a caller set, and answers `var has only automatic row names, so
  rownames(var) would name the genes '1' to '2'. Set rownames(var) to the gene
  identifiers, or pass var_gene_id_column.` Explicit row names are untouched,
  including ones that happen to read as numbers. The `cellucid` Python package
  already refused the same input, because a default pandas `RangeIndex` yields
  integers where its gene identifiers must be strings.
- Refuses a `data.frame` whose columns carry a class wherever it accepts one in
  place of a matrix. `as.matrix()` keeps the numbers and drops the attribute
  that says what they mean, so a frame of one `units` column in metres and one
  in kilometres became a matrix of equal numbers. `latent_space`,
  `gene_expression`, and the embeddings already rejected such a frame;
  `vector_fields` accepted it and published a 45-degree arrow for data pointing
  0.06 degrees off axis, with no error and nothing downstream able to notice.
  All four axes now apply the one rule. `is.numeric()` answers `FALSE` for
  `Date`, `POSIXct`, `difftime`, and `factor`, so those were already refused;
  what this reaches is every classed numeric with no `is.numeric()` method,
  such as `units::units` and `bit64::integer64`.
- Refuses a classed value for every argument that names an identifier or a
  numeric setting, as the string and logical arguments already did. A classed
  character vector passed as `obs_keys` or `gene_identifiers`, or held in the
  column `var_gene_id_column` selects, reached the check that compares the
  written var manifest against the staged gene names -- `identical()` compares
  attributes -- and failed it with an internal message naming a staging path.
  A classed number passed as `compression`, `var_quantization`,
  `obs_continuous_quantization`, or `centroid_outlier_quantile` reached the
  manifest writer and failed with jsonlite's `No method asJSON S3 class:`,
  which names no argument. Both now report the argument and its rule before any
  file is staged, and the identifier message is
  `<argument> must be a native character vector.`

# cellucid 0.9.0

- Initial GitHub package release.
