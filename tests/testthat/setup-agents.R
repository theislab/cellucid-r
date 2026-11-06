# Every viewer a test starts is published to agents on this computer
# (R/agent-link.R): into a folder of this test run's own, never the user's.
# Child processes (cellucid_mcp() run over stdio, Python's cellucid mcp)
# inherit it.
Sys.setenv(CELLUCID_AGENTS_DIR = normalizePath(file.path(tempdir(), "cellucid-agents"), winslash = "/", mustWork = FALSE))
