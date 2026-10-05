## ============================================================
## 00_config.R -- paths and environment, sourced by every script
## ============================================================
## Nothing here is specific to the machine the analysis was run on.
## Set the roots as environment variables, or edit the defaults.
##
##   HMS_ANALYSIS   this directory (holds *.R and results/)
##   HMS_OUT        where the pipeline writes its CSV outputs
##   HMS_DATA       directory holding HMS_<wave>_PUBLIC_instchars.dta   [RESTRICTED]
##   CFPS_DATA      directory holding cfps<year>person.dta              [RESTRICTED]
##   R_LIBS_EXTRA   optional extra R library path
##
## The two survey directories are not in this repository. Both cohorts are
## released under a data-use agreement; see README.md "Data availability".
## ============================================================

ev <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (nzchar(v)) v else default
}

WORKSPACE    <- ev("HMS_WORKSPACE", "")
ANALYSIS_DIR <- ev("HMS_ANALYSIS", if (nzchar(WORKSPACE)) file.path(WORKSPACE, "analysis") else ".")
OUT          <- ev("HMS_OUT", file.path(ANALYSIS_DIR, "outputs"))
FIG_DIR      <- ev("HMS_FIG", file.path(ANALYSIS_DIR, "figures"))
HMS_DIR      <- ev("HMS_DATA", "")     # required by 01, 14, 15, 19a-c, 30, 99c, 99e, probe_*
CFPS_DIR     <- ev("CFPS_DATA", "")    # required by 01, 09, 26, 30, _cfps_n_probe
R_LIBS_EXTRA <- ev("R_LIBS_EXTRA", "")
OPENXLSX_TMP <- ev("OPENXLSX_TMP", tempdir())

## Applied here, centrally: only 31 of the 53 scripts used to carry their own
## .libPaths() call, so the other 22 failed the moment the project library was
## not already on the search path. Found by a full rerun from the raw files.
if (nzchar(R_LIBS_EXTRA)) .libPaths(c(R_LIBS_EXTRA, .libPaths()))

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## Called by scripts that need restricted input, so failure is legible
## rather than a "cannot open file" error from read_dta().
need_dir <- function(x, what) {
  if (!nzchar(x) || !dir.exists(x))
    stop(sprintf("%s is unset or missing (%s). Set it in 00_config.R or as an ",
                 what, if (nzchar(x)) x else "<empty>"),
         "environment variable -- see README.md.", call. = FALSE)
  invisible(TRUE)
}
