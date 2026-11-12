# Actions: everything R changes or reads in an open viewer. The viewer's
# action catalog (inst/action_catalog.json, written by the web app's
# scripts/generate-action-catalog.mjs, byte for byte the web app's
# tests/fixtures/action-catalog.json and cellucid-python's copy, SHA-256 pinned
# in all three) names every action a person can take in the viewer, with its
# parameter schema; cellucid_act() sends any of them, and
# R/actions-generated.R holds one function per action, generated from the
# catalog by tools/generate-actions.R (with its help page and its export),
# so an action the viewer adds reaches R by regenerating, never by hand.
#
# An action is one `action` command (R/viewer-protocol.R): sent at once, or
# kept until the viewer's channel opens, as every command is. The viewer
# applies it as a person's control would -- one step of its history,
# labeled as R's, undone like any other -- and answers with one
# action_result (its outcome, a read's result, and the step it made) or one
# command_error with the reason it refused. Nothing waits for the answer:
# the request comes back at once, and cellucid_action_result() collects the
# answer where R can wait, or in a later notebook cell (a Jupyter notebook
# hears the viewer only between cells).
#
# Parameters are converted to the action's schema here, by annotation and
# never by action: an `x-cells` array takes what a person names cells by in R
# (1-based positions, TRUE/FALSE per cell, the object's cell names) and is
# sent as 0-based original positions in a binary frame; an `x-values` column
# goes in its value type; an array a vector of any length; an object a named
# list. Whether the parameters fit the action is the viewer's to say: its
# schema is the authority, and a viewer newer than this package may know
# more actions or parameters, which reach it unchanged.

.action_catalog_cache <- new.env(parent = emptyenv())

# The actions of the catalog this package carries, by dotted type.
.action_catalog <- function() {
  actions <- .action_catalog_cache$actions
  if (is.null(actions)) {
    path <- system.file("action_catalog.json", package = "cellucid")
    if (!nzchar(path)) .abort("The installed cellucid package has no action_catalog.json: reinstall it.")
    document <- .read_json_file(path)
    actions <- stats::setNames(document$actions, vapply(document$actions, function(entry) entry$type, ""))
    .action_catalog_cache$actions <- actions
  }
  actions
}

# R's reserved words: a parameter of that name is the argument <name>_ in R
# (`in` is `in_`).
.r_reserved_words <- c(
  "if", "else", "repeat", "while", "function", "for", "next", "break", "TRUE", "FALSE", "NULL", "Inf", "NaN",
  "NA", "NA_integer_", "NA_real_", "NA_complex_", "NA_character_", "in"
)

# The parameter an R argument name stands for.
.action_param_name <- function(argument) {
  stripped <- sub("_$", "", argument)
  if (endsWith(argument, "_") && stripped %in% .r_reserved_words) stripped else argument
}

# --- Converting R values to an action's parameters ------------------------------

# The value's shape, as JSON would carry it.
.json_shape <- function(value) {
  if (is.null(value)) return("null")
  if (is.list(value) && !is.data.frame(value)) return(if (is.null(names(value))) "array" else "object")
  if (length(value) == 1L && !inherits(value, "AsIs")) "scalar" else "array"
}

# The form of a oneOf a value takes: the first that admits its shape (an
# object whose constant properties it matches, before one it does not); NULL
# when none does (the value is then sent as it is, and the viewer says why it
# does not fit).
.pick_form <- function(value, forms) {
  shape <- .json_shape(value)
  admits <- function(form) {
    if (isTRUE(form[["x-cells"]]) || isTRUE(form[["x-values"]])) return(shape %in% c("array", "scalar"))
    types <- unlist(form$type)
    if (is.null(types)) return(identical(shape, "scalar"))
    switch(shape,
      null = "null" %in% types,
      object = "object" %in% types,
      array = "array" %in% types,
      scalar = any(c("string", "integer", "number", "boolean") %in% types)
    )
  }
  candidates <- Filter(admits, forms)
  if (length(candidates) == 0L) return(NULL)
  if (identical(shape, "object") && length(candidates) > 1L) {
    matching <- Filter(function(form) {
      for (key in names(form$properties)) {
        constant <- form$properties[[key]]$const
        if (!is.null(constant) && !identical(value[[key]], constant)) return(FALSE)
      }
      TRUE
    }, candidates)
    if (length(matching) > 0L) return(matching[[1L]])
  }
  candidates[[1L]]
}

# One parameter's R value as the JSON value (with typed arrays) its schema
# describes; `schema` NULL: a parameter the catalog does not describe, sent
# as JSON of its shape.
.action_value <- function(viewer, value, schema, label) {
  if (is.null(value)) return(NULL)
  if (inherits(value, "cellucid_typed_array")) return(value)
  if (is.function(value) || is.environment(value) || is.data.frame(value)) {
    .abort(label, " must be text, numbers, TRUE/FALSE, a vector, a list or NULL.")
  }
  if (!is.null(schema)) {
    if (isTRUE(schema[["x-cells"]])) return(.typed_array(.cell_positions(viewer, value, label), "uint32"))
    if (isTRUE(schema[["x-values"]])) return(.column_values(value, label))
    if (!is.null(schema$oneOf)) return(.action_value(viewer, value, .pick_form(value, schema$oneOf), label))
  }
  if (is.factor(value)) value <- stats::setNames(as.character(value), names(value))
  shape <- .json_shape(value)
  types <- if (is.null(schema)) NULL else unlist(schema$type)
  if (identical(shape, "object")) {
    properties <- schema$properties
    out <- structure(list(), names = character())
    for (key in names(value)) {
      if (!nzchar(key)) .abort(label, " names every item of its list.")
      out[key] <- list(.action_value(viewer, value[[key]], properties[[key]], paste0(label, "$", key)))
    }
    return(out)
  }
  if (identical(shape, "array") || (identical(shape, "scalar") && "array" %in% types &&
                                     !any(c("string", "integer", "number", "boolean") %in% types))) {
    items <- if (is.list(value)) unname(value) else as.list(unname(value))
    return(lapply(seq_along(items), function(at) {
      .action_value(viewer, items[[at]], schema$items, paste0(label, "[[", at, "]]"))
    }))
  }
  value <- unname(value)
  if (is.character(value) && is.na(value) || is.logical(value) && is.na(value)) {
    .abort(label, " is NA; give a value, or NULL for none.")
  }
  if (is.numeric(value) && is.na(value)) .abort(label, " is NA; give a value, or NULL for none.")
  if (is.numeric(value) && !is.finite(value)) .abort(label, " must be a finite number.")
  if (!is.character(value) && !is.numeric(value) && !is.logical(value)) {
    .abort(label, " must be text, numbers, TRUE/FALSE, a vector, a list or NULL.")
  }
  value
}

# The parameters of `action`, converted to its schema (`entry` the
# catalog's; NULL for an action it does not describe).
.action_params <- function(viewer, entry, params) {
  if (length(params) > 0L && (is.null(names(params)) || any(!nzchar(names(params))))) {
    .abort("Name every parameter of an action, for example field = list(source = \"obs\", key = \"cell_type\").")
  }
  names(params) <- vapply(names(params), .action_param_name, "", USE.NAMES = FALSE)
  if (anyDuplicated(names(params))) {
    .abort("The parameter ", .format_value(names(params)[anyDuplicated(names(params))]), " is given twice.")
  }
  properties <- entry$params$properties
  out <- structure(list(), names = character())
  for (key in names(params)) {
    out[key] <- list(.action_value(viewer, params[[key]], properties[[key]], key))
  }
  out
}

# Cells (an x-cells parameter) as 0-based original positions, in the order
# given, each once: 1-based positions, TRUE/FALSE for every cell, or the cell
# names of the object the viewer shows. The viewer names a cell by its
# original position, the object's own column (viewer$cell_names).
.cell_positions <- function(viewer, cells, label) {
  n <- .dataset_cell_count(viewer)
  refuse <- function() {
    .abort(label, " must be cell names, 1-based positions", if (!is.null(n)) paste0(" from 1 to ", n),
           ", or TRUE/FALSE per cell.")
  }
  if (is.list(cells)) {
    if (!all(vapply(cells, function(cell) is.atomic(cell) && length(cell) == 1L, TRUE)) ||
        length(unique(vapply(cells, function(cell) if (is.numeric(cell)) "numeric" else class(cell)[[1L]], ""))) > 1L) {
      refuse()
    }
    cells <- unlist(cells, use.names = FALSE)
    if (is.null(cells)) cells <- numeric(0)
  }
  if (is.factor(cells)) cells <- as.character(cells)
  if (is.logical(cells)) {
    if (!is.null(n) && length(cells) != n) .abort(label, " as TRUE/FALSE needs one value per cell (", n, ").")
    if (anyNA(cells)) {
      .abort(label, " as TRUE/FALSE holds NA at cell ", which(is.na(cells))[[1L]], "; say TRUE or FALSE for every cell.")
    }
    cells <- which(cells)
  } else if (is.character(cells)) {
    names <- viewer$cell_names
    if (is.null(names)) .abort(label, " by name needs a viewer of an R object; give positions.")
    positions <- match(cells, names)
    if (anyNA(positions)) {
      .abort(label, " names cells the dataset does not have: ", .format_value_list(utils::head(cells[is.na(positions)], 5L)), ".")
    }
    ambiguous <- cells[cells %in% names[duplicated(names)]]
    if (length(ambiguous) > 0L) {
      .abort(label, " names the cell ", .format_value(ambiguous[[1L]]), ", which the object holds more than once; give positions.")
    }
    cells <- positions
  } else if (!is.numeric(cells)) {
    refuse()
  }
  cells <- unname(cells)
  if (length(cells) > 0L) {
    # One pass, no copy (src/extent.c): no NA, every position a whole number
    # from 1 to the number of cells (at most 2^32, a uint32 each).
    extent <- .native_extent(cells, whole = is.double(cells))
    if (extent[[3L]] > 0 || extent[[1L]] < 1 || (!is.null(n) && extent[[2L]] > n) ||
        extent[[2L]] > 4294967296 || (is.double(cells) && extent[[5L]] != 1)) {
      refuse()
    }
  }
  repeated <- anyDuplicated(cells)
  if (repeated > 0L) {
    .abort(label, " names the cell ", format(cells[[repeated]], scientific = FALSE), " twice; name each cell once.")
  }
  as.double(cells) - 1
}

# A column (an x-values parameter) in its value type: a factor as its codes
# (0-based, NA the code type's largest value: no category), TRUE/FALSE as
# 1/0, whole numbers from 0 in the smallest unsigned type that holds them,
# other numbers as float64 (NA: no value).
.column_values <- function(values, label) {
  if (is.list(values)) values <- unlist(values, use.names = FALSE)
  if (is.null(values)) values <- numeric(0)
  if (is.factor(values)) {
    levels <- nlevels(values)
    type <- if (levels < 255) "uint8" else if (levels < 65535) "uint16" else "uint32"
    codes <- as.double(as.integer(values) - 1L)
    codes[is.na(codes)] <- 2^(8 * .viewer_value_types[[type]]) - 1
    return(.typed_array(codes, type))
  }
  if (is.logical(values)) {
    if (anyNA(values)) {
      .abort(label, " holds NA: give a factor (NA: no category), or numbers (NA: no value).")
    }
    return(.typed_array(as.integer(values), "uint8"))
  }
  if (!is.numeric(values)) .abort(label, " takes numbers, TRUE/FALSE, or a factor (its codes).")
  values <- unname(values)
  if (any(is.infinite(values))) .abort(label, " holds an infinity; give finite numbers (NA: no value).")
  if (is.integer(values)) {
    if (anyNA(values)) {
      .abort(label, " holds NA: give a factor (NA: no category), or numbers (as.numeric(); NA: no value).")
    }
    if (length(values) == 0L || min(values) >= 0L) {
      largest <- if (length(values) == 0L) 0 else max(values)
      type <- if (largest <= 255) "uint8" else if (largest <= 65535) "uint16" else "uint32"
      return(.typed_array(values, type))
    }
  }
  .typed_array(as.double(values), "float64")
}

# --- Requests -----------------------------------------------------------------

# Send `action` with `params` (R values, converted here); `caller` is the R
# call that sent it, for what a refusal says. Returns the request at once,
# invisibly.
.act <- function(viewer, action, params, caller) {
  .require_viewer(viewer)
  if (!.is_exact_text(action) || !grepl(.viewer_action_type_pattern, action, perl = TRUE)) {
    .abort("An action is named by its dotted type, such as \"color.by\"; not ", .format_value(action), ".")
  }
  entry <- .action_catalog()[[action]]
  converted <- .action_params(viewer, entry, params)
  id <- paste0("a", .random_hex(8L))
  # Inside cellucid_turn(), the action is one step of that turn.
  frame <- .encode_viewer_action(id, action, converted, viewer$turn$id, viewer$turn$label)
  if (!isTRUE(viewer$running)) .abort("This viewer is stopped; open a new one with cellucid_view().")
  record <- new.env(parent = emptyenv())
  record$action <- action
  record$caller <- caller
  record$sent <- FALSE
  record$waiting <- FALSE
  record$orphan <- FALSE
  viewer$action_requests[[id]] <- record
  .viewer_send_frame(viewer, frame, record)
  request <- new.env(parent = emptyenv())
  request$viewer <- viewer
  request$id <- id
  request$record <- record
  class(request) <- "cellucid_action_request"
  # The answer is kept while the request is: once nothing holds it, its
  # answer (a read's result can be large) is let go.
  reg.finalizer(request, .let_go_action_request)
  # Invisibly, as a command's call: an action is mostly sent for what it
  # does; print it, or collect it, to see it.
  invisible(request)
}

.let_go_action_request <- function(request) {
  viewer <- get("viewer", envir = request, inherits = FALSE)
  id <- get("id", envir = request, inherits = FALSE)
  record <- viewer$action_requests[[id]]
  if (is.null(record)) return(invisible(NULL))
  if (.action_answered(record)) {
    viewer$action_requests[[id]] <- NULL
  } else {
    # Answered later: the answer is dropped as it arrives.
    record$orphan <- TRUE
  }
  invisible(NULL)
}

.action_answered <- function(record) !is.null(record$state)

# A generated function called without a parameter its action requires.
.missing_action_param <- function(name, argument) {
  .abort(name, "() needs ", argument, ": see ?", name, ".")
}

cellucid_act <- function(.viewer, .action, ...) {
  .act(.viewer, .action, list(...), paste0("cellucid_act(\"", .action, "\")"))
}

cellucid_action_result <- function(request, timeout = 60) {
  if (!inherits(request, "cellucid_action_request")) {
    .abort("request must be what cellucid_act() or a cellucid_<action>() function returned.")
  }
  if (!is.numeric(timeout) || length(timeout) != 1L || is.na(timeout) || timeout < 0) {
    .abort("timeout must be one non-negative number of seconds.")
  }
  viewer <- get("viewer", envir = request, inherits = FALSE)
  record <- get("record", envir = request, inherits = FALSE)
  if (!.action_answered(record)) {
    if (!is.null(viewer$widget)) {
      .abort("The viewer has not answered ", record$caller, " yet: in a Jupyter notebook it answers between cells. ",
             "Run cellucid_action_result() again in a later cell.")
    }
    record$waiting <- TRUE
    on.exit(record$waiting <- FALSE, add = TRUE)
    .wait_for(viewer, function() .action_answered(record), timeout, paste0("the answer to ", record$caller))
  }
  switch(record$state,
    answered = {
      # A change the viewer applied without a step of its history says why
      # (a read records none by nature, and says nothing); the outcome is
      # returned all the same.
      if (!is.null(record$unrecorded)) .warn_not_recorded(record)
      record$outcome
    },
    refused = .abort("The viewer did not take ", record$caller, ": ", record$reason),
    stopped = .abort("The viewer was stopped before it answered ", record$caller, ".")
  )
}

# Fields a request answers by name: its id, its action, whether it was
# answered (done), the history step it made (step: list(id, label), or NULL
# for a read or a change that made none), why a change made none
# (unrecorded: a name of .viewer_action_unrecorded, or NULL) and its state.
`$.cellucid_action_request` <- function(x, name) {
  record <- get("record", envir = x, inherits = FALSE)
  switch(name,
    id = get("id", envir = x, inherits = FALSE),
    action = record$action,
    done = .action_answered(record),
    step = record$step,
    unrecorded = record$unrecorded,
    state = .action_state(x),
    NULL
  )
}

.action_state <- function(request) {
  viewer <- get("viewer", envir = request, inherits = FALSE)
  record <- get("record", envir = request, inherits = FALSE)
  if (!is.null(record$state)) return(record$state)
  if (isTRUE(record$sent)) "sent" else if (isTRUE(viewer$running)) "waiting for the viewer" else "stopped"
}

print.cellucid_action_request <- function(x, ...) {
  record <- get("record", envir = x, inherits = FALSE)
  state <- .action_state(x)
  detail <- if (identical(state, "answered") && !is.null(record$step)) {
    paste0(", step ", record$step$id, ": ", record$step$label)
  } else if (identical(state, "answered") && !is.null(record$unrecorded)) {
    paste0(", no step: ", record$unrecorded)
  } else if (identical(state, "refused")) {
    paste0(": ", record$reason)
  } else {
    ""
  }
  cat("<cellucid action ", record$action, " (", get("id", envir = x, inherits = FALSE), "): ", state, detail, ">\n",
      sep = "")
  invisible(x)
}

# An action's answer, or its refusal, as the viewer sent it (.deliver_event()):
# kept on its request for cellucid_action_result(); an answer nobody holds a
# request for any more is dropped. Returns whether the one asking is waiting
# for it now (a refusal is then said by that wait, not by a warning).
.settle_action <- function(viewer, id, outcome = NULL, step = NULL, unrecorded = NULL, reason = NULL) {
  record <- viewer$action_requests[[id]]
  if (is.null(record)) return(FALSE)
  if (is.null(reason)) {
    record$state <- "answered"
    record$outcome <- outcome
    record$step <- step
    record$unrecorded <- unrecorded
  } else {
    record$state <- "refused"
    record$reason <- reason
  }
  if (isTRUE(record$orphan)) viewer$action_requests[[id]] <- NULL
  isTRUE(record$waiting)
}

# A change answered without a step of the history, said when its outcome is
# collected: a warning of class cellucid_action_not_recorded (its
# `unrecorded` names the reason), muffled like any warning, or by its class.
.warn_not_recorded <- function(record) {
  warning(.classed_warning(
    "cellucid_action_not_recorded", list(unrecorded = record$unrecorded),
    record$caller, ": ", .viewer_action_unrecorded[[record$unrecorded]]
  ))
}

# Every request still open when the viewer stops: none will be answered.
.stop_action_requests <- function(viewer) {
  for (id in names(viewer$action_requests)) {
    record <- viewer$action_requests[[id]]
    if (is.null(record$state)) record$state <- "stopped"
  }
  viewer$action_requests <- list()
  invisible(NULL)
}
