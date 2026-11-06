# Session files of format 4 (the History v2 shape), made in the tests.
#
# The shared golden fixture (fixtures/history_golden) holds every kind of
# record the history model has; these helpers write files of the same format
# with other content (other datasets, labels, fields, subsamples, branches)
# and take files apart to damage them, so the reader's and the applier's rules
# are tested on every case -- cellucid-python's tests/_session_files.py, helper
# for helper. They follow the format the web app writes (app/session/README.md,
# app/session/history-file.js): steps, parts and values as JSON cut into
# pieces, each chunk gzip'd where that saves a sixteenth of its bytes (else
# stored as it is), every array once by its digest, every sub-value of at least
# 512 characters of skeleton text once as a part, the values the active head
# shows and the parts and arrays they name eager and everything else lazy,
# every record the model's (inst/history_model.schema.json).

session_created_at <- "2026-11-03T12:00:00.000Z"
# The time of a test tracker's opening (2026-01-01T00:00:00Z).
session_start <- 1767225600000

json_object <- function(...) {
  out <- list(...)
  if (length(out) == 0L) names(out) <- character(0)
  out
}

session_array_sizes <- c(Uint8Array = 1L, Uint16Array = 2L, Uint32Array = 4L, Float32Array = 4L, Float64Array = 8L)

# A typed array as a value holds one (Uint32Array, Float32Array...).
typed <- function(type, values) {
  size <- session_array_sizes[[type]]
  values <- unlist(values)
  bytes <- if (type %in% c("Float32Array", "Float64Array")) {
    writeBin(as.double(values), raw(), size = size, endian = "little")
  } else if (identical(type, "Uint32Array")) {
    cellucid:::.encode_uint32(as.double(values))
  } else {
    writeBin(as.integer(values), raw(), size = size, endian = "little")
  }
  structure(list(type = type, bytes = bytes, digest = paste0(type, ":", length(values), ":", cellucid:::.native_sha256(bytes))),
            class = "session_test_array")
}

session_cells <- function(values) typed("Uint32Array", values)

session_mapping <- function(...) {
  out <- list(
    kind = "anndata",
    embeddings = list(enabled = list(list(key = "X_umap", columns = list(0L, 1L))),
                      defaults = list(`1d` = NULL, `2d` = list(key = "X_umap", columns = list(0L, 1L)), `3d` = NULL)),
    knn = NULL, geneNames = "index", expression = "X", log1p = json_object(), highlyVariable = NULL, latent = NULL,
    vectorFields = list()
  )
  changes <- list(...)
  for (key in names(changes)) out[key] <- list(changes[[key]])
  out
}

# A dataset fingerprint; `digest` is the cell-order digest of `embedding`'s
# columns 0 and 1 over the cells (an apply checks it when the object has that
# embedding).
session_fingerprint <- function(cell_count, var_count, dataset_id = "demo", data_mapping = session_mapping(),
                                embedding = "X_umap", digest = "0123456789abcdef") {
  list(datasetId = dataset_id, name = "Demo", cellCount = cell_count, varCount = var_count,
       cellOrder = list(embedding = list(key = embedding, columns = list(0L, 1L)), dimension = 2L, digest = digest),
       dataMapping = data_mapping)
}

session_overlays <- function(renamed_fields = list(), renamed_categories = list(), deleted_fields = list(), user_fields = list()) {
  list(renamedFields = renamed_fields, renamedCategories = renamed_categories, deletedFields = deleted_fields,
       userFields = user_fields)
}

session_user_category <- function(field_id, key, categories, codes, cells_of = NULL, code_type = "Uint8Array", deleted = FALSE) {
  list(id = field_id, source = "obs", kind = "category", key = key, deleted = deleted, purged = FALSE,
       createdAt = session_start, sourceField = NULL, operation = list(type = "annotation"), categories = as.list(categories),
       codes = list(cells = if (is.null(cells_of)) NULL else session_cells(cells_of), values = typed(code_type, codes)),
       sourcePages = list(), overlapStrategy = "first", overlapLabel = NULL, intersectionLabels = NULL, uncoveredLabel = NULL)
}

session_user_values <- function(field_id, key, values, cells_of = NULL, deleted = FALSE) {
  list(id = field_id, source = "obs", kind = "continuous", key = key, deleted = deleted, purged = FALSE,
       createdAt = session_start, sourceField = NULL, operation = list(type = "notebook-push"),
       values = list(cells = if (is.null(cells_of)) NULL else session_cells(cells_of), values = typed("Float32Array", values)))
}

session_user_copy <- function(field_id, key, source, source_key) {
  list(id = field_id, source = source, kind = "continuous", key = key, deleted = FALSE, purged = FALSE,
       createdAt = session_start, sourceField = list(kind = "continuous-copy", sourceKey = source_key),
       operation = list(type = "copy"))
}

session_highlights <- function(..., active = NULL) {
  pages <- list(...)
  list(active = active %||% if (length(pages) > 0L) pages[[1L]]$id, pages = pages)
}

session_page <- function(page_id, name, ..., color = "#2563eb") list(id = page_id, name = name, color = color, groups = list(...))

session_group <- function(group_id, label, members, kind = "lasso", enabled = TRUE, ...) {
  c(list(id = group_id, type = kind, label = label, enabled = enabled, cells = session_cells(members)), list(...))
}

# Every domain an application reads, at its opening value, with `changes`
# (domain = value).
session_state <- function(...) {
  out <- list(
    cells = list(open = NULL, genes = NULL),
    mapping = session_mapping(),
    `field-overlays` = session_overlays(),
    `gene-sets` = list(sets = list()),
    highlights = session_highlights(),
    camera = list(views = list(list(zoom = -0, ratio = NaN)))
  )
  changes <- list(...)
  for (key in names(changes)) out[[key]] <- changes[[key]]
  out
}

# The skeleton text from which an object or a list is stored as a part
# (value-parts.js PART_MIN_CHARS), and the text a JSON chunk is cut after
# (history-file.js CHUNK_TEXT_CHARS).
session_part_min_chars <- 512L
session_chunk_text_chars <- 16 * 1024 * 1024

# A value as the file stores it: each typed array a {"$array": digest}
# placeholder (kept in `arrays`), each number JSON cannot hold a {"$number"},
# and, when `parts` is given, each object or list (not the value itself) whose
# skeleton text is at least session_part_min_chars a {"$part": <sha256 of that
# text>} (kept in `parts`, bottom-up).
session_split <- function(value, arrays, parts = NULL, root = TRUE) {
  if (inherits(value, "session_test_array")) {
    if (!exists(value$digest, envir = arrays, inherits = FALSE)) {
      assign(value$digest, value, envir = arrays)
      arrays$.order <- c(arrays$.order, value$digest)
    }
    return(list(`$array` = value$digest))
  }
  if (is.double(value) && length(value) == 1L && !is.na(value) && value == 0 && 1 / value < 0) return(list(`$number` = "-0"))
  if (is.double(value) && length(value) == 1L && is.nan(value)) return(list(`$number` = "NaN"))
  if (is.double(value) && length(value) == 1L && is.infinite(value)) return(list(`$number` = if (value > 0) "Infinity" else "-Infinity"))
  if (is.list(value)) {
    out <- lapply(value, session_split, arrays = arrays, parts = parts, root = FALSE)
    if (!is.null(names(value))) names(out) <- names(value)
    if (is.null(parts) || root) return(out)
    text <- enc2utf8(cellucid:::.js_json(out))
    if (nchar(text, type = "chars") < session_part_min_chars) return(out)
    name <- cellucid:::.native_sha256(charToRaw(text))
    if (!exists(name, envir = parts, inherits = FALSE)) {
      assign(name, out, envir = parts)
      parts$.order <- c(parts$.order, name)
    }
    return(list(`$part` = name))
  }
  value
}

# The part names a skeleton names directly.
session_part_names <- function(skeleton) {
  if (is.list(skeleton)) {
    if (identical(names(skeleton), "$part")) return(skeleton[["$part"]])
    return(unlist(lapply(skeleton, session_part_names), use.names = FALSE))
  }
  character()
}

# JSON text cut after `limit` characters (history-file.js jsonPieces): a list
# (kind "list") or a map (kind "map", a named list); one item is never split.
session_pieces <- function(kind, items, limit = session_chunk_text_chars) {
  texts <- if (identical(kind, "list")) {
    vapply(items, cellucid:::.js_json, "", USE.NAMES = FALSE)
  } else {
    paste0(cellucid:::.js_strings(names(items)), ":", vapply(items, cellucid:::.js_json, "", USE.NAMES = FALSE))
  }
  brackets <- if (identical(kind, "list")) c("[", "]") else c("{", "}")
  out <- list()
  current <- character()
  size <- 0
  for (text in texts) {
    if (length(current) > 0L && size + nchar(text, type = "chars") > limit) {
      out[[length(out) + 1L]] <- charToRaw(enc2utf8(paste0(brackets[[1L]], paste(current, collapse = ","), brackets[[2L]])))
      current <- character()
      size <- 0
    }
    current <- c(current, text)
    size <- size + nchar(text, type = "chars") + 1
  }
  if (length(current) > 0L) out[[length(out) + 1L]] <- charToRaw(enc2utf8(paste0(brackets[[1L]], paste(current, collapse = ","), brackets[[2L]])))
  out
}

# A chunk's manifest entry and stored bytes (history-file.js storeChunk): with
# codec "auto" gzip'd where that saves at least a sixteenth, else as it is.
session_store <- function(meta, payload, codec = "auto") {
  zipped <- session_gzip(payload)
  keep <- identical(codec, "gzip") || (identical(codec, "auto") && length(zipped) * 16 <= length(payload) * 15)
  body <- if (keep) zipped else payload
  list(meta = c(meta, list(codec = if (keep) "gzip" else "none", storedBytes = length(body), uncompressedBytes = length(payload))),
       bytes = body)
}

# A stored chunk's payload.
session_payload <- function(entry, chunk) {
  if (identical(entry$codec, "gzip")) cellucid:::.native_inflate(chunk, TRUE) else chunk
}

session_digests <- function(skeleton) {
  if (is.list(skeleton)) {
    if (identical(names(skeleton), "$array")) return(skeleton[["$array"]])
    return(unlist(lapply(skeleton, session_digests), use.names = FALSE))
  }
  character()
}

session_json <- function(value) charToRaw(enc2utf8(cellucid:::.js_json(value)))

session_gzip <- function(bytes) cellucid:::.framed(list(zarr_format = 3L), bytes, memCompress(bytes, "gzip"))

# A tracker id (t-<uuid>) for a test's short name ("t1"), the same every time
# (cellucid-python's tracker_id()).
session_tracker_id <- function(label) {
  if (grepl("^t-[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", label)) return(label)
  digest <- cellucid:::.native_sha256(charToRaw(enc2utf8(label)))
  part <- function(from, to) substr(digest, from, to)
  paste0("t-", part(1, 8), "-", part(9, 12), "-4", part(14, 16), "-8", part(18, 20), "-", part(21, 32))
}

session_by <- function(origin = "user", agent = NULL, model = NULL, host = NULL, turn = NULL) {
  list(origin = origin, agent = agent, model = model, host = host, person = NULL, turn = turn)
}

# A session file's manifest and stored chunks, framed.
session_frame <- function(manifest, chunks) {
  text <- session_json(manifest)
  c(cellucid:::.session_magic, cellucid:::.encode_uint32(length(text)), text,
    unlist(lapply(chunks, function(chunk) c(cellucid:::.encode_uint32(length(chunk)), chunk)), use.names = FALSE))
}

# A session file's manifest and stored chunks.
session_unframe <- function(bytes) {
  at <- length(cellucid:::.session_magic)
  length <- readBin(bytes[at + 1:4], "integer", size = 4L, endian = "little")
  manifest <- jsonlite::parse_json(rawToChar(bytes[at + 4L + seq_len(length)]), simplifyVector = FALSE)
  at <- at + 4L + length
  chunks <- list()
  for (entry in manifest$chunks) {
    size <- readBin(bytes[at + 1:4], "integer", size = 4L, endian = "little")
    chunks[[length(chunks) + 1L]] <- if (size == 0L) raw(0) else bytes[at + 4L + seq_len(size)]
    at <- at + 4L + size
  }
  list(manifest = manifest, chunks = chunks)
}

# Write a session file of the History v2 shape (every record the model's).
# Each tracker is list(id, name, root = list(domain = value), steps = list(...),
# head?, branchNames? = list(b1 = "name")); each step list(set = list(domain =
# value), parent?, branch?, type?, label?, status?, turn?, agent?, host?,
# group?, squashed?): its `after` the values it sets that differ from its
# parent's, its `before` its parent's; a step with `squashed` (list(cause,
# count, by)) is a squashed step, with no action. A tracker's `groups` are
# list(list(id, name)), its `recording = FALSE` pauses it (`paused` since its
# opening, by the person); a squashed step's `carried` is empty unless given. content = "steps"
# writes the steps only (no values); "state" marks a state without its history
# (its one tracker the root and one squashed step of cause "share").
session_file <- function(path, fingerprint, trackers, active = NULL, requires = list(obs = list(), var = list()),
                         cells_record = NULL, app_version = "0.9.6", content = "history", times = "exact",
                         piece_chars = session_chunk_text_chars, codec = "auto") {
  values <- !identical(content, "steps")
  domains <- cellucid:::.history_domains
  arrays <- new.env(parent = emptyenv())
  arrays$.order <- character()
  parts <- new.env(parent = emptyenv())
  parts$.order <- character()
  held <- new.env(parent = emptyenv())
  hash_of <- function(value) {
    hash <- cellucid:::.value_hash(session_split(value, new.env(parent = emptyenv())))
    assign(hash, value, envir = held)
    hash
  }
  number <- function(id) as.integer(substring(id, 2L))
  records <- list()
  types <- character()
  for (index in seq_along(trackers)) {
    tracker <- trackers[[index]]
    created <- session_start + 1000000 * (index - 1L)
    root_after <- json_object()
    for (domain in intersect(domains, names(tracker$root))) root_after[[domain]] <- hash_of(tracker$root[[domain]])
    root <- list(id = "s0", parent = NULL, children = list(), branch = "b0", kind = "root", type = NULL, params = NULL, v = NULL,
                 code = NULL, label = "Opened", by = session_by(), at = created, until = created, merged = 1L, ms = NULL, touches = list(),
                 before = json_object(), after = root_after, reads = json_object(), status = "applied", reason = NULL,
                 outcome = NULL, coalesceKey = NULL, replayOf = NULL, squashed = NULL, notes = list(), bookmark = NULL, group = NULL)
    steps <- list(s0 = root)
    hashes_at <- list(s0 = root_after)
    for (position in seq_along(tracker$steps)) {
      step <- tracker$steps[[position]]
      parent <- step$parent %||% steps[[length(steps)]]$id
      status <- step$status %||% "applied"
      after <- json_object()
      before <- json_object()
      for (domain in intersect(domains, names(step$set))) {
        hash <- hash_of(step$set[[domain]])
        if (identical(status, "applied") && identical(hashes_at[[parent]][[domain]], hash)) next
        after[[domain]] <- hash
        before[[domain]] <- hashes_at[[parent]][[domain]]
      }
      if (identical(status, "unapplicable")) for (domain in names(after)) after[[domain]] <- hashes_at[[parent]][[domain]]
      squashed <- step$squashed
      if (!is.null(squashed) && is.null(squashed$carried)) squashed$carried <- list(notes = list(), bookmarks = list())
      action <- if (is.null(squashed)) step$type %||% "test.step" else NULL
      types <- c(types, action)
      origin <- if (!is.null(step$agent)) "agent" else if (!is.null(step$host)) "host" else "user"
      id <- step$id %||% paste0("s", position)
      record <- list(
        id = id, parent = parent, children = list(), branch = step$branch %||% "b0",
        kind = if (is.null(squashed)) "action" else "squashed", type = action,
        params = if (is.null(squashed)) session_split(step$params %||% json_object(), arrays, parts) else NULL,
        v = if (is.null(squashed)) 1L else NULL,
        code = if (is.null(squashed)) paste0(gsub(".", "_", action, fixed = TRUE), "()") else NULL,
        label = step$label %||% paste("Step", position),
        by = session_by(origin, agent = step$agent, host = step$host, turn = step$turn),
        at = created + position * 1000, until = created + position * 1000, merged = 1L, ms = 1.5,
        touches = as.list(intersect(domains, names(after))), before = before, after = after, reads = json_object(),
        status = status, reason = step$reason, outcome = NULL, coalesceKey = NULL, replayOf = NULL, squashed = squashed,
        notes = list(), bookmark = step$bookmark, group = step$group
      )
      steps[[parent]]$children <- c(steps[[parent]]$children, list(id))
      steps[[id]] <- record
      hashes_at[[id]] <- utils::modifyList(hashes_at[[parent]], after)
    }
    ordered <- steps[order(vapply(steps, function(step) number(step$id), 0L))]
    head <- tracker$head %||% steps[[length(steps)]]$id
    branch_ids <- unique(vapply(ordered, function(step) step$branch, ""))
    branches <- lapply(branch_ids, function(branch) {
      own <- Filter(function(step) identical(step$branch, branch), ordered)
      list(id = branch, name = tracker$branchNames[[branch]] %||% if (identical(branch, "b0")) "Main" else paste("Branch", substring(branch, 2L)),
           base = own[[1L]]$parent, tip = own[[length(own)]]$id, createdAt = own[[1L]]$at, by = session_by(),
           lane = number(branch) %% 8L, locked = FALSE, archived = FALSE, notes = list())
    })
    branches <- branches[order(vapply(branches, function(branch) number(branch$id), 0L))]
    turn_ids <- unique(unlist(lapply(ordered, function(step) step$by$turn)))
    turns <- lapply(turn_ids, function(turn) {
      taken <- Filter(function(step) identical(step$by$turn, turn), ordered)
      list(id = turn, by = taken[[1L]]$by, label = "A turn", conversation = NULL, message = NULL, startedAt = taken[[1L]]$at,
           endedAt = taken[[length(taken)]]$at, branch = taken[[1L]]$branch, notes = list())
    })
    records[[index]] <- list(
      manifest = list(
        id = session_tracker_id(tracker$id), name = tracker$name, dataset = list(id = "demo", name = "Demo", cellOrder = NULL),
        createdAt = created, by = session_by(), start = list(kind = "opened", from = NULL, report = NULL, mapping = NULL),
        removed = NULL, rootId = "s0", headId = head, activeBranch = steps[[head]]$branch, mainBranch = "b0",
        paused = if (identical(tracker$recording, FALSE)) list(at = created, by = list(session_by()), draft = NULL) else NULL,
        branches = branches,
        groups = lapply(tracker$groups, function(group) {
          list(id = group$id, name = group$name, createdAt = created, by = session_by(), notes = list())
        }),
        turns = turns, notes = list(),
        nextStep = max(vapply(ordered, function(step) number(step$id), 0L)) + 1L,
        nextBranch = max(vapply(branches, function(branch) number(branch$id), 0L)) + 1L,
        nextGroup = if (length(tracker$groups) == 0L) 0L else max(vapply(tracker$groups, function(group) number(group$id), 0L)) + 1L,
        nextNote = 0L, nextJournal = 0L, steps = list(), journal = NULL, chats = list()
      ),
      steps = unname(ordered),
      head_hashes = hashes_at[[head]]
    )
  }
  ids <- vapply(records, function(record) record$manifest$id, "")
  active_id <- if (is.null(active)) ids[[1L]] else session_tracker_id(active)
  active_record <- records[[match(active_id, ids)]]
  referenced <- unique(unlist(lapply(records, function(record) {
    lapply(record$steps, function(step) c(unlist(step$before), unlist(step$after)))
  })))
  eager_hashes <- unique(unlist(active_record$head_hashes))
  skeleton_of <- function(hash) session_split(get(hash, envir = held), arrays, parts)
  eager <- if (values) stats::setNames(lapply(eager_hashes, skeleton_of), eager_hashes) else json_object()
  lazy_hashes <- sort(setdiff(referenced, names(eager)), method = "radix")
  lazy <- if (values && length(lazy_hashes) > 0L) stats::setNames(lapply(lazy_hashes, skeleton_of), lazy_hashes) else json_object()
  # What the eager values reach, through the parts.
  eager_parts <- character()
  eager_arrays <- character()
  queue <- unname(eager)
  while (length(queue) > 0L) {
    skeleton <- queue[[1L]]
    queue <- queue[-1L]
    eager_arrays <- union(eager_arrays, session_digests(skeleton))
    for (name in setdiff(session_part_names(skeleton), eager_parts)) {
      eager_parts <- c(eager_parts, name)
      queue[[length(queue) + 1L]] <- get(name, envir = parts)
    }
  }
  entries <- list()
  add <- function(meta, payload) entries[[length(entries) + 1L]] <<- list(meta = meta, payload = payload)
  for (index in seq_along(records)) {
    pieces <- session_pieces("list", unname(records[[index]]$steps), piece_chars)
    for (piece in seq_along(pieces)) {
      id <- paste0("steps/", index - 1L, "/", piece - 1L)
      records[[index]]$manifest$steps[[piece]] <- id
      add(list(id = id, kind = "steps", tracker = ids[[index]], priority = "eager"), pieces[[piece]])
    }
  }
  for (priority in c("eager", "lazy")) {
    own_parts <- Filter(function(name) (name %in% eager_parts) == identical(priority, "eager"), parts$.order)
    own_parts <- stats::setNames(lapply(own_parts, get, envir = parts), own_parts)
    if (length(own_parts) > 0L) {
      pieces <- session_pieces("map", own_parts, piece_chars)
      for (piece in seq_along(pieces)) add(list(id = paste0("parts/", priority, "/", piece - 1L), kind = "parts", priority = priority), pieces[[piece]])
    }
    for (digest in Filter(function(digest) (digest %in% eager_arrays) == identical(priority, "eager"), arrays$.order)) {
      add(list(id = paste0("array/", digest), kind = "array", digest = digest, encoding = "raw", priority = priority),
          get(digest, envir = arrays)$bytes)
    }
    own_values <- if (identical(priority, "eager")) eager else lazy
    if (length(own_values) > 0L) {
      pieces <- session_pieces("map", own_values, piece_chars)
      for (piece in seq_along(pieces)) add(list(id = paste0("values/", priority, "/", piece - 1L), kind = "values", priority = priority), pieces[[piece]])
    }
  }
  stored <- lapply(entries, function(entry) session_store(entry$meta, entry$payload, codec))
  chunks <- lapply(stored, function(entry) entry$bytes)
  metas <- lapply(stored, function(entry) entry$meta)
  manifest <- c(
    list(
      createdAt = session_created_at,
      app = list(version = app_version, catalog = strrep("ab", 32L)),
      actions = lapply(sort(unique(types), method = "radix"), function(type) list(type = type, version = 1L, params = strrep("cd", 32L))),
      scope = list(
        trackers = as.list(ids),
        branches = lapply(records, function(record) {
          list(tracker = record$manifest$id, branches = lapply(record$manifest$branches, function(branch) branch$id))
        }),
        content = content, results = if (values) "include" else "recompute", notes = !identical(content, "state"),
        chats = FALSE, journal = FALSE, times = times, source = "full"
      ),
      datasetFingerprint = fingerprint,
      source = NULL,
      requires = requires,
      chrome = list(pressed = json_object())
    ),
    if (!is.null(cells_record)) list(cells = cells_record),
    list(activeTracker = active_id, trackers = lapply(records, function(record) record$manifest), chunks = metas)
  )
  writeBin(session_frame(manifest, chunks), path)
  path
}

# The manifest `cells` of a uniform subsample, and its kept cells (0-based).
session_subsample_record <- function(n_obs, k, seed) {
  kept <- cellucid:::.native_subsample_uniform(n_obs, k, seed)
  list(
    record = list(spec = list(method = "uniform", nObs = n_obs, k = k, seed = seed),
                  kept = jsonlite::base64_enc(cellucid:::.subsample_encode_cell_set(kept)),
                  digest = cellucid:::.subsample_digest(n_obs, kept)),
    kept = kept
  )
}
